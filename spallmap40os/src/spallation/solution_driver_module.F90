MODULE SOLUTION_DRIVER_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for temperature, stress, oxided-induced stresses. 
    !  it will determine the temperature profile in the substrate and coating
    !  in the steady state and during ramp up (down)
    !
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: DRIVER_SOLUTION

 CONTAINS

  SUBROUTINE DRIVER_SOLUTION(N, cycle_no, Temp_out, Temp_in, p_out, p_in, &
     & thickness_oxide_layer, time_now, dt, id_high_or_low, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the solution including temperature, thermal expansion,
    !         strain, and stresses profile
    ! 
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module, only: moxide, mave, id_low, id_high
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & total_oxide_thickness_ref
  use induced_stress_module, only: SET_OX_STRAINS, GET_U_JUMP_SURFACE
  use creep_data_module, only  : if_creep_metal, if_creep_ox_scale
  use creep_module, only: OUT_VON_MISES
  use boiler_data_module, only: TEMP_BOILER_data_no
  use solution_data_module, only: rad_int, eps0, cc_st, bc_st
  use oxide_data_module,   only: first_oxide_layer
  use solver_data_module, only: update_sola_type

  ! Argument List
  integer, intent(IN)  :: N, id_high_or_low, cycle_no
  real,    intent(IN)  :: Temp_out, Temp_in, time_now, dt
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale
  real,    intent(IN)  :: p_out, p_in

  ! Local Variables
  integer  :: iter, i
  logical  :: convergent = .false., if_oh_solution

  ! oxide thickness without thermal expasion; just given by the kinetics
  ! total_oxide_thickness_ref = SUM(thickness_oxide_layer(2:N))
  
  ! set the oxidation strains only once; does not work for lateral oxide growth
  call SET_OX_STRAINS(N, cycle_no)

  call GET_U_JUMP_SURFACE(N, thickness_oxide_layer, iter)

  if (if_creep_metal .or. if_creep_ox_scale)  then

    if (update_sola_type == 0)  then

      call SOLUTION_TH_GS_CR(N, cycle_no, Temp_out, Temp_in, p_out, p_in, &
        & thickness_oxide_layer, time_now, dt, id_high_or_low, if_ave_scale)

    else if (update_sola_type == 10)  then
      ! consider only the metal to assess the creep effect only
      call SOLUTION_TH_GS_CR(N-N+1, cycle_no, Temp_out, Temp_in, &
        & p_out, p_in, thickness_oxide_layer, time_now, dt, &
        & id_high_or_low, if_ave_scale)

    else if (update_sola_type >= 1 .and. update_sola_type <= 3)  then

      call SOLUTION_NR_GS_CR(N, cycle_no, Temp_out, Temp_in, p_out, p_in, &
        & thickness_oxide_layer, time_now, dt, id_high_or_low, if_ave_scale)

    endif

  else

    ! no creep here
    GS_IF: if (if_growth_induce_stress)  then
  
      call SOLUTION_TH_GS(N, cycle_no, Temp_out, Temp_in, p_out, p_in, &
       & thickness_oxide_layer, time_now, dt, id_high_or_low, if_ave_scale)

    else if (.not. if_growth_induce_stress)  then

      call SOLUTION_ONLY_TH(N, cycle_no, Temp_out, Temp_in, p_out, p_in, &
       & thickness_oxide_layer, time_now, dt, id_high_or_low, if_ave_scale)

    endif GS_IF

    write(aux_lun, 17) cycle_no, time_now, eps0, &
         & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
         & (bc_st(i), i= first_oxide_layer, N), (cc_st(i), i= first_oxide_layer, N)
 17    format('const_conv_sol_e_b=', i5, 1x, 30(1pe13.6, 1x))

  endif

  ! print out von Mises stress; exclude test problem output
  if (cycle_no < TEMP_BOILER_data_no)  then
    call OUT_VON_MISES(N, id_high_or_low, cycle_no, time_now, dt)
  endif

  return

  END SUBROUTINE DRIVER_SOLUTION

  SUBROUTINE SOLUTION_NR_GS_CR(N, cycle_no, Temp_out, Temp_in, p_out, p_in, &
     & thickness_oxide_layer, time_now, dt, id_high_or_low, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   thermal expansion, growth induced stresses, and creep
    !       
    !   Newton-Rhaphson for the solution including temperature, thermal expansion,
    !         strain, and stresses profile
    ! 
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module, only: moxide, mave, id_low, id_high
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int, no_stress_var
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, if_cylindrical, no_st_rad_points_tube, &
      & no_st_rad_points_oxide, update_sola_type, &
      & if_different_ref_temp, Temp_reference_tube
  use solver_data_module, only  : ratio_time_first_pulse, &
      & no_first_pulse_creep_int, start_creep_cycles, &
      & iter_stress_old, max_iterat_solver
  use property_module, only: OPERAND_ARRAY
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & total_oxide_thickness_ref
  use thermal_cycle_module, only: TEMP_PROFILE, BETA_TAU_TH
  use stress_solution_module, only: STRESS_PROFILE
  use mesh_module, only: GET_RADIUS, GET_MESH
  use induced_stress_module, only: GET_U_JUMP_SURFACE
  use creep_module, only: CREEP_STRAIN, CREEP_CONV, OUT_VON_MISES, &
         & CREEP_INIT_LOCATION
  use oxide_utility_module, only: AVERAGE_PROFILE
  use stress_nr_module, only:  STRESS_PROFILE_NR, SOLVER_DELTA, &
         & SETUP_STRESS_MATRIX, SOLVER_VAR, START_CREEP_STRAIN_MEAN

  ! Argument List
  integer, intent(IN)  :: N, id_high_or_low, cycle_no
  real,    intent(IN)  :: Temp_out, Temp_in, time_now, dt
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale
  real,    intent(IN)  :: p_out, p_in

  ! Local Variables
  integer  :: iter, k, nc_it, ncreep_it
  logical  :: convergent = .false., if_oh_solution = .false.
  real, dimension(200) :: dt_factor
  real     :: dt_cr, time_now_cr
  real, dimension(no_stress_var)      :: var, fval, delta_var 
  real, dimension(no_stress_var, no_stress_var) :: coeff
  
  ! default for large cycle after creep was initiated   
   ncreep_it = 1
   dt_factor = 1.0

   ! cycle_no - counts all the transition points not the number of pulses
   ! if (cycle_no <= start_creep_cycles)  then
   ! start_creep_cycles counts the number of pulses
   ! the end of first pulse has cycle_no = 4 and second cycle_no = 8
   if (cycle_no <= 4 * start_creep_cycles)  then

     if (id_high_or_low == id_high .or. id_high_or_low == id_low)  then
       ncreep_it = no_first_pulse_creep_int
     endif

     if (ABS(ratio_time_first_pulse - 1.0) < 1.0e-8)  then

       dt_factor = 1.0 / no_first_pulse_creep_int

     else

       dt_factor(1) = (ratio_time_first_pulse - 1.0) / &
            & (ratio_time_first_pulse**no_first_pulse_creep_int - 1.0)

       do k = 2, no_first_pulse_creep_int
         dt_factor(k) = dt_factor(k-1) * ratio_time_first_pulse
       end do

     endif

   else if (cycle_no > 4 * start_creep_cycles)  then

     if (iter_stress_old >= 6)  then
       ncreep_it = 1
       dt_factor = 1.0 / ncreep_it
       if (ncreep_it > 1)  then
         write(out_lun, 3) cycle_no, time_now, ncreep_it, iter_stress_old
 3       format(i7, 1x, 1pe13.6, ' interval divided in ', i2, 1x, &
            & 'creep_intervals ', i2)
       endif
       ! write(tty_lun, 3) cycle_no, time_now, ncreep_it, iter_stress_old
        
     endif

   endif

  ! oxide thickness without thermal expasion; just given by the kinetics
  ! total_oxide_thickness_ref = SUM(thickness_oxide_layer(2:N))
    ! updatet the mesh for low and high cycles as well
    ! if (id_high_or_low == id_high)  then

        ! if_th_oxide_growth = .false.  ! no thermal expansion on h^2=Aexp(Q/RT)
        ! if_th_tube = .false.    ! no thermal expansion for the tube 
        ! neglect the effect of displacement of the tube due to 
        ! thermal expansion and applied pressure on the temperature field
        ! get the radius where temperature and stress will be computed for this
        ! cycle
      call GET_RADIUS(N, thickness_oxide_layer, time_now, &
          & if_ave_scale, iter)

      call GET_MESH(N, thickness_oxide_layer, if_ave_scale)

      ! call GET_U_JUMP_SURFACE(N, thickness_oxide_layer, iter)

    ! endif

      ! temperature eqn need current dimensions
      ! temperature is in the loop since the oxide growth temperature 
      ! depends on location
    call TEMP_PROFILE(N, Temp_out, Temp_in, thickness_oxide_layer, &
             & cycle_no, time_now, dt, id_high_or_low, if_ave_scale)

    call BETA_TAU_TH(N, if_ave_scale, id_high_or_low)

    ! consider that creep does NOT affects dimensions;
    ! do not get temperature for the start-up creep time intervals
    ! since the oxide thickness is very small

    ! intialize time - useful for intiating creep
    time_now_cr = time_now - dt

    CREEP_ST1: do nc_it = 1, ncreep_it

      dt_cr = dt * dt_factor(nc_it)
      time_now_cr = time_now_cr + dt_cr

      if (cycle_no > 4)  then  ! initial there will be no previous mesh
        ! initialize the OLD creep strains for nodes in NEW oxide regions
        call CREEP_INIT_LOCATION(N, cycle_no, time_now_cr)
      endif

      ! update stresses and strains
      ! displacements and dimensions are computed in the current and
      ! reference states
      STRESS_LOOP: do iter = 1, max_iterat_solver

        ! initialize von Mises stress since it oscillates between high/low loads
        ! unlike creep strain which accumulates in time and is continous
        call START_CREEP_STRAIN_MEAN(iter, cycle_no, N, id_high_or_low, dt_cr)

        ! get coeff and fval using var at previous iteration
        call SETUP_STRESS_MATRIX(N, cycle_no, iter, dt_cr, p_out, p_in, &
          & id_high_or_low, if_ave_scale)

        if (update_sola_type == 1)  then

          ! solve A*x=rhs, x = (eps0, b, c)
          ! update s_r, s_h, s_z
          ! update vm0 and delta_e_cr
          ! update e_r, e_h, e_z
          call SOLVER_VAR(iter, cycle_no, time_now_cr, dt_cr, convergent)

        else if (update_sola_type == 2)  then

          ! solve A*x=rhs, x = (eps0, b, c, s_r, s_h, s_z, e_r, e_h, e_z)
          ! update vm0 and delta_e_cr
          call SOLVER_VAR(iter, cycle_no, time_now_cr, dt_cr, convergent)

        else if (update_sola_type == 3)  then

          ! solve coeff*delta_var = -fval; obtaine var = var + delta_var
          ! x = (eps0, b, c, s_r, s_h, s_z, e_r, e_h, e_z, vm0, delta_e_cr)
          call SOLVER_DELTA(iter, cycle_no, time_now_cr, dt_cr, convergent)

        endif

        ! get displacements and dimensions; only for OH algorithm
        ! call GET_OG_DIM_DISPL(N, poisson, th_exp, Youngs, &
        ! & poisson_star, th_exp_star, Youngs_star, if_oh_solution)

        if (convergent)  exit STRESS_LOOP

      end do STRESS_LOOP

      call STRESS_PROFILE_NR(N, cycle_no, iter, nc_it, dt, time_now_cr, &
       & p_out, p_in, id_high_or_low, if_ave_scale, if_oh_solution, &
       & convergent)

      ! update creep strain, creep integrals
      ! von Mises stress is not updated it depends on old stresses
      ! call CREEP_STRAIN(iter, cycle_no, N, id_high_or_low, dt_cr)

      ! store the maximum number of stress iterations 
      if (nc_it == 1)  then
        iter_stress_old = iter
      else if (nc_it > 1)  then
        if (iter > iter_stress_old)  iter = iter_stress_old
      endif

      if (ncreep_it > 1)  then !  .and. cycle_no <= 4 * start_creep_cycles)  then
        call AVERAGE_PROFILE(N, time_now_cr, id_high_or_low, if_ave_scale)
        call OUT_VON_MISES(N, id_high_or_low, cycle_no, time_now_cr, dt)
      endif

    end do CREEP_ST1

    ! make sure that you call this one within the stress iteration loop
    ! to study the convergence
    call AVERAGE_PROFILE(N, time_now, id_high_or_low, if_ave_scale)

  return

  END SUBROUTINE SOLUTION_NR_GS_CR

  SUBROUTINE SOLUTION_TH_GS_CR(N, cycle_no, Temp_out, Temp_in, p_out, p_in, &
     & thickness_oxide_layer, time_now, dt, id_high_or_low, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   thermal expansion, growth induced stresses, and creep
    !       
    !   Obtain the solution including temperature, thermal expansion,
    !         strain, and stresses profile
    ! 
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module, only: moxide, mave, id_low, id_high
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer, time_boiler_idle, time_boiler_pulse
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, if_cylindrical, no_st_rad_points_tube, &
      & no_st_rad_points_oxide, &
      & if_different_ref_temp, Temp_reference_tube
  use solver_data_module, only  : ratio_time_first_pulse, &
      & no_first_pulse_creep_int, start_creep_cycles, &
      & iter_stress_old, max_iterat_solver
  use property_module, only: OPERAND_ARRAY
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & total_oxide_thickness_ref
  use thermal_cycle_module, only: TEMP_PROFILE, BETA_TAU_TH
  use stress_solution_module, only: STRESS_PROFILE
  use mesh_module, only: GET_RADIUS, GET_MESH
  use induced_stress_module, only: GET_U_JUMP_SURFACE
  use solver_module, only: NEWTON_RAPHSON
  use creep_module, only: CREEP_STRAIN, CREEP_CONV, OUT_VON_MISES, &
         & CREEP_INIT_LOCATION
  use oxide_utility_module, only: AVERAGE_PROFILE

  ! Argument List
  integer, intent(IN)  :: N, id_high_or_low, cycle_no
  real,    intent(IN)  :: Temp_out, Temp_in, time_now, dt
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale
  real,    intent(IN)  :: p_out, p_in

  ! Local Variables
  integer  :: iter, k, nc_it, ncreep_it
  logical  :: convergent = .false., if_oh_solution = .false.
  real, dimension(200) :: dt_factor
  real     :: dt_cr, time_now_cr, dt_limit_cr, dt_increment_cr
  logical, save :: if_ncreep_it_first = .true. 
  real, save    :: dt_limit_low_load, dt_limit_full_load, &
      & dt_increment_full_load, dt_increment_low_load
 
  if (if_ncreep_it_first)  then

    if (time_boiler_idle(1) > 10.0)  then
      ! weekly schedule
      dt_limit_low_load = 20.0
      dt_increment_low_load = 5.0
    else
      ! daily schedule
      dt_limit_low_load = 6.0
      dt_increment_low_load = 3.0 ! 0.5 ! 1.0 ! 2.0
    endif

    if (time_boiler_pulse(1) > 20.0)  then
      ! weekly schedule
      dt_limit_full_load = 20.0
      dt_increment_full_load = 5.0
    else
      ! daily schedule
      dt_limit_full_load = 6.0
      dt_increment_full_load = 4.0 ! 0.5  ! 1.0 ! 2.0
    endif

    if_ncreep_it_first = .false.

  endif

  if (id_high_or_low == id_high)  then
    dt_limit_cr = dt_limit_full_load 
    dt_increment_cr = dt_increment_full_load
  else if (id_high_or_low == id_low)  then
    dt_limit_cr = dt_limit_low_load 
    dt_increment_cr = dt_increment_low_load
    
    ! the transitions f2l and l2f are considered to be low load
    if (dt < 3.0) then
      dt_limit_cr = 3.0 ! 0.5
      dt_increment_cr = 1.0 ! 0.25
    endif

  endif

  ! default for large cycle after creep was initiated   
   ncreep_it = 1
   dt_factor = 1.0

   ! cycle_no - counts all the transition points not the number of pulses
   ! if (cycle_no <= start_creep_cycles)  then
   ! start_creep_cycles counts the number of pulses
   ! the end of first pulse has cycle_no = 4 and second cycle_no = 8
   if (cycle_no <= 4 * start_creep_cycles)  then

     if (id_high_or_low == id_high .or. id_high_or_low == id_low)  then
       ncreep_it = no_first_pulse_creep_int
     endif

     if (ABS(ratio_time_first_pulse - 1.0) < 1.0e-8)  then

       dt_factor = 1.0 / no_first_pulse_creep_int

     else

       dt_factor(1) = (ratio_time_first_pulse - 1.0) / &
            & (ratio_time_first_pulse**no_first_pulse_creep_int - 1.0)

       do k = 2, no_first_pulse_creep_int
         dt_factor(k) = dt_factor(k-1) * ratio_time_first_pulse
       end do

     endif

   else if (cycle_no > 4 * start_creep_cycles)  then

     if (dt > dt_limit_cr)  then !  (dt > 20.0)

       ncreep_it = AINT(dt / dt_increment_cr) !  AINT(dt / 5.0)
       if (ncreep_it == 0) ncreep_it = 5

       dt_factor(1:ncreep_it-1) = 1.0 / ncreep_it

       dt_factor(ncreep_it)= dt - dt_factor(1) * (ncreep_it - 1)

     else if (dt <= dt_limit_cr)  then

       ! do not use many intervals for solving the creep problem
       if (iter_stress_old >= 6)  then
         ncreep_it = 1
         dt_factor = 1.0 / ncreep_it
         if (ncreep_it > 1)  then
           write(out_lun, 3) cycle_no, time_now, ncreep_it, iter_stress_old
 3         format(i7, 1x, 1pe13.6, ' interval divided in ', i2, 1x, &
            & 'creep_intervals ', i2)
         endif
         ! write(tty_lun, 3) cycle_no, time_now, ncreep_it, iter_stress_old
        
       endif

     endif

   endif

  ! oxide thickness without thermal expasion; just given by the kinetics
  ! total_oxide_thickness_ref = SUM(thickness_oxide_layer(2:N))
    ! updatet the mesh for low and high cycles as well
    ! if (id_high_or_low == id_high)  then

        ! if_th_oxide_growth = .false.  ! no thermal expansion on h^2=Aexp(Q/RT)
        ! if_th_tube = .false.    ! no thermal expansion for the tube 
        ! neglect the effect of displacement of the tube due to 
        ! thermal expansion and applied pressure on the temperature field
        ! get the radius where temperature and stress will be computed for this
        ! cycle
      call GET_RADIUS(N, thickness_oxide_layer, time_now, &
          & if_ave_scale, iter)

      call GET_MESH(N, thickness_oxide_layer, if_ave_scale)

      ! call GET_U_JUMP_SURFACE(N, thickness_oxide_layer, iter)

    ! endif

      ! temperature eqn need current dimensions
      ! temperature is in the loop since the oxide growth temperature 
      ! depends on location
    call TEMP_PROFILE(N, Temp_out, Temp_in, thickness_oxide_layer, &
             & cycle_no, time_now, dt, id_high_or_low, if_ave_scale)

    call BETA_TAU_TH(N, if_ave_scale, id_high_or_low)

    ! consider that creep does NOT affects dimensions;
    ! do not get temperature for the start-up creep time intervals
    ! since the oxide thickness is very small

    ! intialize time - useful for intiating creep
    time_now_cr = time_now - dt

    CREEP_ST1: do nc_it = 1, ncreep_it

      dt_cr = dt * dt_factor(nc_it)
      time_now_cr = time_now_cr + dt_cr

      if (cycle_no > 4)  then  ! initial there will be no previous mesh
        ! initialize the OLD creep strains for nodes in NEW oxide regions
        call CREEP_INIT_LOCATION(N, cycle_no, time_now_cr)
      endif

      STRESS_LOOP: do iter = 1, max_iterat_solver

        ! update creep strain, creep integrals
        ! von Mises stress is not updated it depends on old stresses
        call CREEP_STRAIN(iter, cycle_no, N, id_high_or_low, dt_cr, time_now_cr)

        ! update stresses and strains
        ! displacements and dimensions are computed in the current and
        ! reference states
        call STRESS_PROFILE(N, cycle_no, p_out, p_in, &
            & id_high_or_low, if_ave_scale, if_oh_solution)

        ! update von Mises now
        call CREEP_CONV(iter, nc_it, N, cycle_no, time_now_cr, &
            & dt_cr, id_high_or_low, convergent)

        if (convergent)  exit STRESS_LOOP

      end do STRESS_LOOP

      ! store the maximum number of stress iterations 
      if (nc_it == 1)  then
        iter_stress_old = iter
      else if (nc_it > 1)  then
        if (iter > iter_stress_old)  iter = iter_stress_old
      endif

      if (ncreep_it > 1)  then !  .and. cycle_no <= 4 * start_creep_cycles)  then
        call AVERAGE_PROFILE(N, time_now_cr, id_high_or_low, if_ave_scale)
        call OUT_VON_MISES(N, id_high_or_low, cycle_no, time_now_cr, dt)
      endif

    end do CREEP_ST1

    ! make sure that you call this one within the stress iteration loop
    ! to study the convergence
    call AVERAGE_PROFILE(N, time_now, id_high_or_low, if_ave_scale)

  return

  END SUBROUTINE SOLUTION_TH_GS_CR

  SUBROUTINE SOLUTION_TH_GS(N, cycle_no, Temp_out, Temp_in, p_out, p_in, &
     & thickness_oxide_layer, time_now, dt, id_high_or_low, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   thermal expansion and growth induced stresses
    !       (NO creep)
    !   Obtain the solution including temperature, thermal expansion,
    !         strain, and stresses profile
    ! 
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module, only: moxide, mave, id_low, id_high
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, if_cylindrical, no_st_rad_points_tube, &
      & no_st_rad_points_oxide, &
      & if_different_ref_temp, Temp_reference_tube
  use property_module, only: OPERAND_ARRAY
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & total_oxide_thickness_ref
  use thermal_cycle_module, only: TEMP_PROFILE, BETA_TAU_TH
  use stress_solution_module, only: STRESS_PROFILE
  use mesh_module, only: GET_RADIUS, GET_MESH
  use induced_stress_module, only: GET_U_JUMP_SURFACE
  use solver_module, only: NEWTON_RAPHSON
  use oxide_utility_module, only: AVERAGE_PROFILE

  ! Argument List
  integer, intent(IN)  :: N, id_high_or_low, cycle_no
  real,    intent(IN)  :: Temp_out, Temp_in, time_now, dt
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale
  real,    intent(IN)  :: p_out, p_in

  ! Local Variables
  integer  :: iter
  logical  :: convergent = .false., if_oh_solution = .false.

  ! oxide thickness without thermal expasion; just given by the kinetics
  ! total_oxide_thickness_ref = SUM(thickness_oxide_layer(2:N))
  
    STRESS_LOOP: do iter = 1, 50

      if (id_high_or_low == id_high)  then

        if_oh_solution = .false.

        if (if_oh_solution)  then

          call OH_SOLUTION(N, Temp_out, Temp_in, p_out, p_in, &
               & thickness_oxide_layer, time_now, id_high_or_low, &
               & if_ave_scale)

        endif

        ! if_th_oxide_growth = .false.  ! no thermal expansion on h^2=Aexp(Q/RT)
        ! if_th_tube = .false.    ! no thermal expansion for the tube 
        ! neglect the effect of displacement of the tube due to 
        ! thermal expansion and applied pressure on the temperature field
        ! get the radius where temperature and stress will be computed for this
        ! cycle
        call GET_RADIUS(N, thickness_oxide_layer, time_now, &
          & if_ave_scale, iter)

        call GET_MESH(N, thickness_oxide_layer, if_ave_scale)

        ! call GET_U_JUMP_SURFACE(N, thickness_oxide_layer, iter)

      endif

      ! temperature eqn need current dimensions
      ! temperature is in the loop since the oxide growth temperature 
      ! depends on location
      call TEMP_PROFILE(N, Temp_out, Temp_in, thickness_oxide_layer, &
             & cycle_no, time_now, dt, id_high_or_low, if_ave_scale)

      call BETA_TAU_TH(N, if_ave_scale, id_high_or_low)

      ! also the displacements are computed and dimensions in the current and
      ! reference states
      call STRESS_PROFILE(N, cycle_no, p_out, p_in, &
          & id_high_or_low, if_ave_scale, if_oh_solution)

      call AVERAGE_PROFILE(N, time_now, id_high_or_low, if_ave_scale)
 
      if (if_oh_solution .and. convergent)  then
        exit STRESS_LOOP
      else
        ! do only one loop if steiner
        exit STRESS_LOOP
      endif

    end do STRESS_LOOP

  return

  END SUBROUTINE SOLUTION_TH_GS

  SUBROUTINE SOLUTION_ONLY_TH(N, cycle_no, Temp_out, Temp_in, p_out, p_in, &
     & thickness_oxide_layer, time_now, dt, id_high_or_low, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   only thermal expansion here 
    !       (NO growth induced stresses; NO creep)
    !
    !   Obtain the solution including temperature, thermal expansion,
    !         strain, and stresses profile
    ! 
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module, only: moxide, mave, id_low, id_high
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, if_cylindrical, no_st_rad_points_tube, &
      & no_st_rad_points_oxide, &
      & if_different_ref_temp, Temp_reference_tube
  use property_module, only: OPERAND_ARRAY
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & total_oxide_thickness_ref
  use thermal_cycle_module, only: TEMP_PROFILE, BETA_TAU_TH
  use stress_solution_module, only: STRESS_PROFILE
  use mesh_module, only: GET_RADIUS, GET_MESH
  use induced_stress_module, only: GET_U_JUMP_SURFACE
  use solver_module, only: NEWTON_RAPHSON
  use oxide_utility_module, only: AVERAGE_PROFILE

  ! Argument List
  integer, intent(IN)  :: N, id_high_or_low, cycle_no
  real,    intent(IN)  :: Temp_out, Temp_in, time_now, dt
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale
  real,    intent(IN)  :: p_out, p_in

  ! Local Variables
  integer  :: iter
  logical  :: convergent = .false., if_oh_solution = .false.

  ! oxide thickness without thermal expasion; just given by the kinetics
  ! total_oxide_thickness_ref = SUM(thickness_oxide_layer(2:N))
  
    iter = 1  ! one iteration 

    if (id_high_or_low == id_high)  then

      ! if_th_oxide_growth = .false.  ! no thermal expansion on h^2=Aexp(Q/RT)
      ! if_th_tube = .false.    ! no thermal expansion for the tube 
      ! neglect the effect of displacement of the tube due to 
      ! thermal expansion and applied pressure on the temperature field
        ! get the radius where temperature and stress will be computed for this
        ! cycle
      call GET_RADIUS(N, thickness_oxide_layer, time_now, &
        & if_ave_scale, iter)

      call GET_MESH(N, thickness_oxide_layer, if_ave_scale)

      ! just for consistency; du_jump = 0
      ! no jump conditions when oxide induced stresses are excluded
      ! call GET_U_JUMP_SURFACE(N, thickness_oxide_layer, iter)

    endif

    ! temperature eqn need current dimensions
    call TEMP_PROFILE(N, Temp_out, Temp_in, thickness_oxide_layer, &
             & cycle_no, time_now, dt, id_high_or_low, if_ave_scale)

    call BETA_TAU_TH(N, if_ave_scale, id_high_or_low)

    call STRESS_PROFILE(N, cycle_no, p_out, p_in, &
         & id_high_or_low, if_ave_scale, if_oh_solution)

    call AVERAGE_PROFILE(N, time_now, id_high_or_low, if_ave_scale)

  return

  END SUBROUTINE SOLUTION_ONLY_TH

  SUBROUTINE OH_SOLUTION(N, Temp_out, Temp_in, p_out, p_in, &
     & thickness_oxide_layer, time_now, id_high_or_low, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the solution including temperature, thermal expansion,
    !         strain, and stresses profile
    ! 
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module, only: moxide, mave, id_low, id_high
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, if_cylindrical, no_st_rad_points_tube, &
      & no_st_rad_points_oxide, &
      & if_different_ref_temp, Temp_reference_tube
  use property_module, only: OPERAND_ARRAY
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & u_now_int_jump, total_oxide_thickness_ref
  use thermal_cycle_module, only: TEMP_PROFILE, BETA_TAU_TH
  use stress_solution_module, only: STRESS_PROFILE
  use mesh_module, only: GET_RADIUS_OH, GET_MESH
  use induced_stress_module, only: OXIDE_INDUCED_UPDATE, OXIDE_INDUCED_SETUP
  use solver_module, only: NEWTON_RAPHSON

  ! Argument List
  integer, intent(IN)  :: N, id_high_or_low
  real,    intent(IN)  :: Temp_out, Temp_in, time_now
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale
  real,    intent(IN)  :: p_out, p_in

  ! Local Variables
  integer  :: iter
  logical  :: convergent = .false., if_oh_solution = .true.

  ! oxide thickness without thermal expasion; just given by the kinetics
  ! total_oxide_thickness_ref = SUM(thickness_oxide_layer(2:N))
  
    STRESS_LOOP: do iter = 1, 50

      if (id_high_or_low == id_high)  then

        ! if_th_oxide_growth = .false.  ! no thermal expansion on h^2=Aexp(Q/RT)
        ! if_th_tube = .false.    ! no thermal expansion for the tube 
        ! neglect the effect of displacement of the tube due to 
        ! thermal expansion and applied pressure on the temperature field
        ! get the radius where temperature and stress will be computed for this
        ! cycle
        call GET_RADIUS_OH(N, thickness_oxide_layer, if_ave_scale, iter)

        call GET_MESH(N, thickness_oxide_layer, if_ave_scale)

      endif

      ! temperature eqn need current dimensions
      ! temperature is in the loop since the oxide growth temperature 
      ! depends on location
      ! this is not done since cycle_no is not here so I used N
      !      also I used time_now-time_now insted of dt
      call TEMP_PROFILE(N, Temp_out, Temp_in, thickness_oxide_layer, &
          & N, time_now, time_now-time_now, id_high_or_low, &
          & if_ave_scale)

      call BETA_TAU_TH(N, if_ave_scale, id_high_or_low)

      ! also the displacements are computed and dimensions in the current and
      ! reference states
      call STRESS_PROFILE(N, 1, p_out, p_in, id_high_or_low, &
         & if_ave_scale, if_oh_solution)

      if (convergent)  then
        exit STRESS_LOOP
      endif

      ! setup coefficients for the constraint equations
      call OXIDE_INDUCED_SETUP(iter, id_high_or_low)

      ! solve the equations and increment variables solved for
      call NEWTON_RAPHSON(iter, convergent)

      ! update new dimensions and displacement conditions
      call OXIDE_INDUCED_UPDATE(iter, id_high_or_low)

    end do STRESS_LOOP

  return

  END SUBROUTINE OH_SOLUTION

END MODULE SOLUTION_DRIVER_MODULE
