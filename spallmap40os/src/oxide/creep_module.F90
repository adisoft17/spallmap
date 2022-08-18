MODULE CREEP_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for creep. 
    ! 
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================


  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: CREEP_STRAIN, CREEP_CONV, CREEP_TERM_S_RAD, CREEP_TERM_U, &
          & CREEP_TERM_DUDR, CREEP_TERM_S_HOOP, CREEP_TERM_INT_AXIAL, &
          & INIT_CREEP, OUT_VON_MISES, CREEP_INIT_LOCATION, &
          & DCREEP_DVM, DCREEP, VON_MISES, OUT_CREEP_CONV

 CONTAINS

  SUBROUTINE CREEP_INIT_LOCATION(N, cycle_no, time_now)

    !=======================================================================
    ! Purpose(s):
    !
    ! Reinitialize the creep strain based on the new oxide position
    ! Rule: new oxide just created has ZERO creep strain 
    !=======================================================================

  use parameter_module,   only: moxide, mave
  use boiler_data_module, only: tube_thickness, tube_outer_radius
  use solution_data_module, only: rad_temp, npr_temp, &
      & npr_st, rad_st, rad_st_old, rad_int
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value, first_oxide_layer
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, no_st_rad_points_tube, &
      & no_st_rad_points_oxide 
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use creep_data_module, only: eps_oxide_grid, il_now, ip_now, &
      & creep_strain_old, creep_strain_new, i_cr_start, i_cr_end, &
      & new_oxide_creep_option 

  integer, intent(IN)  :: N, cycle_no
  real, intent(IN)     :: time_now

  ! Local Variables
  integer  :: i, j, io, jo
  integer, save  :: ip = 0
  real     :: rad_diff, displ_max, slope, summ1, rad_now
  logical  :: reset_creep_strain_old, metal_oxide_interface, &
      & oxide_steam_interface

  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  ! check only the displacement between the previous and current position
  ! for the radii only in the oxide layers

  rad_diff = 0.0
  displ_max = ABS(rad_st(first_oxide_layer, 1) - &
        & rad_st(first_oxide_layer, npr_st(first_oxide_layer)))

  do i = first_oxide_layer, N
   
    summ1 = SUM((rad_st(i, 1:npr_st(i)) - rad_st_old(i, 1:npr_st(i)))**2) 
    rad_diff = rad_diff + summ1

  end do

  if (first_oxide_layer < N)  then
    rad_diff = SQRT(rad_diff / SUM(npr_st(first_oxide_layer:N)))
  else
    rad_diff = 0.0
  endif

  if (rad_diff < eps_oxide_grid * displ_max)  then

    reset_creep_strain_old = .false.

  else

    reset_creep_strain_old = .true.

  endif

  if (.not. reset_creep_strain_old)  RETURN

  ip_now = -10
  il_now = -10

  LAY1: do i = MAX(i_cr_start, first_oxide_layer), i_cr_end    ! 2, N

    metal_oxide_interface = .false.
    oxide_steam_interface = .false.

    if (i == MAX(i_cr_start, first_oxide_layer))  &
        & metal_oxide_interface = .true.
    if (i == i_cr_end)  oxide_steam_interface = .true.

    if (.not. (metal_oxide_interface .or. oxide_steam_interface))  then
      ! only the first and last oxide are considered
      ! keep the old creep strains since identical nodes
      creep_strain_new(i, 1:npr_st(i))%rad = &
           & creep_strain_old(i, 1:npr_st(i))%rad
      creep_strain_new(i, 1:npr_st(i))%hoop = &
           & creep_strain_old(i, 1:npr_st(i))%hoop
      creep_strain_new(i, 1:npr_st(i))%ax = &
           & creep_strain_old(i, 1:npr_st(i))%ax
      CYCLE LAY1

    endif

    ! LOC1: do j = 1, npr_st(i)
    ! no need for j loop since we look only for the first and last points    

    ! locations (2, 1) for m-ox and (N, npr_st(N)) will be updated below
    ! this assignment is just for the interal nodes
    creep_strain_new(i, 1:npr_st(i))%rad = &
           & creep_strain_old(i, 1:npr_st(i))%rad
    creep_strain_new(i, 1:npr_st(i))%hoop = &
           & creep_strain_old(i, 1:npr_st(i))%hoop
    creep_strain_new(i, 1:npr_st(i))%ax = &
           & creep_strain_old(i, 1:npr_st(i))%ax

      ! for each new position find the interval indices in the old mesh
      ! rad_now = rad_st(i, j)
      
      ! how about if the oxide changes its type; not a concern for creep
      ! creep in the metal is not transferred to the oxide; creep in the 
      ! metal is low any way (but could be otherwise)
      ! or new oxide created at metal-oxide interface is creep free
      ! 
      ! no need for loop since only m-ox and ox-steam are the problem
      ! LAY2: do io = 2, N

      if (metal_oxide_interface)  then

        io = first_oxide_layer
        j = 1
        rad_now = rad_st(i, j)

      write(aux_lun, 3) i, time_now, rad_st(i, j), rad_st_old(io, j), &
          & rad_now - rad_st_old(io, j)
 3    format('metal_check_reset_cr ', i2, 1x, 20(1pe13.6, 1x))

        if (rad_now > rad_st_old(io, j))  then
          ! in the inner oxide growth mesh, r(1) > r(2) and so on
          ! this is the spine-metal interface
            ip_now(i, j) = 0
            il_now(i, j) = io
            ! use creep_strain_new as dummies for creep_strain_old
          ! npr, ... 3, 2, 1 grid order within each layer
            creep_strain_new(i, j)%rad = &
       &       CREEP_NEW_OXIDE(new_oxide_creep_option, &
       &   creep_strain_old(io, 2)%rad, &
       &   creep_strain_old(io, 1)%rad, rad_st_old(io, 2), &
       &   rad_st_old(io, 1), rad_st(io, 2), &
       &   rad_st(io, 1))
            creep_strain_new(i, j)%hoop = &
       &       CREEP_NEW_OXIDE(new_oxide_creep_option, &
       &   creep_strain_old(io, 2)%hoop, &
       &   creep_strain_old(io, 1)%hoop, rad_st_old(io, 2), &
       &   rad_st_old(io, 1), rad_st(io, 2), &
       &   rad_st(io, 1))
            creep_strain_new(i, j)%ax = &
       &       CREEP_NEW_OXIDE(new_oxide_creep_option, &
       &   creep_strain_old(io, 2)%ax, &
       &   creep_strain_old(io, 1)%ax, rad_st_old(io, 2), &
       &   rad_st_old(io, 1), rad_st(io, 2), &
       &   rad_st(io, 1))

           write(aux_lun, 1) cycle_no, time_now, &
      &    (rad_now - rad_st_old(io, 1)) / &
       &   (rad_st_old(io, 1) - rad_st_old(io, 2)),  &
       &   creep_strain_new(i, j)%rad - creep_strain_old(io, 1)%rad, & 
       &   creep_strain_new(i, j)%hoop - creep_strain_old(io, 1)%hoop, & 
       &   creep_strain_new(i, j)%ax - creep_strain_old(io, 1)%ax 
 1         format('met_ox_str_cr_corr ', i4, 1x, 20(1pe13.6, 1x))

        else if (ABS(rad_now - rad_st_old(io, 1)) < displ_max * 1.0e-6)  then

          ip_now(i, j) = 0
          il_now(i, j) = io
          ! keep the old creep strains since identical nodes
          creep_strain_new(i, j)%rad = creep_strain_old(i, j)%rad
          creep_strain_new(i, j)%hoop = creep_strain_old(i, j)%hoop
          creep_strain_new(i, j)%ax = creep_strain_old(i, j)%ax

        endif

      endif

      if (oxide_steam_interface)  then

        io = i_cr_end  ! N
        j = npr_st(i)
        rad_now = rad_st(i, j)

      write(aux_lun, 4) i, time_now, rad_st(i, j), &
          & rad_st_old(io, npr_st(io)), &
          & rad_now - rad_st_old(io, npr_st(io))
 4    format('steam_check_reset_cr ', i2, 1x, 20(1pe13.6, 1x))

        if (rad_now < rad_st_old(io, npr_st(io)))  then

          ! in the inner oxide growth mesh, r(1) > r(2) and so on
          ! this is the magnetite-steam interface
          ip_now(i, j) = npr_st(i) + 1
          il_now(i, j) = io
          ! npr, npr -1, ... 3, 2, 1 grid order within each layer
            creep_strain_new(i, j)%rad = &
       &       CREEP_NEW_OXIDE(new_oxide_creep_option, &
       &   creep_strain_old(io, j-1)%rad, &
       &   creep_strain_old(io, j)%rad, rad_st_old(io, j-1), &
       &   rad_st_old(io, j), rad_st(io, j-1), &
       &   rad_st(io, j))
            creep_strain_new(i, j)%hoop = &
       &       CREEP_NEW_OXIDE(new_oxide_creep_option, &
       &   creep_strain_old(io, j-1)%hoop, &
       &   creep_strain_old(io, j)%hoop, rad_st_old(io, j-1), &
       &   rad_st_old(io, j), rad_st(io, j-1), &
       &   rad_st(io, j))
            creep_strain_new(i, j)%ax = &
       &       CREEP_NEW_OXIDE(new_oxide_creep_option, &
       &   creep_strain_old(io, j-1)%ax, &
       &   creep_strain_old(io, j)%ax, rad_st_old(io, j-1), &
       &   rad_st_old(io, j), rad_st(io, j-1), &
       &   rad_st(io, j))

           write(aux_lun, 2) cycle_no, time_now, &
       &    (rad_st_old(io, j) - rad_now) / &
       &   (rad_st_old(io, j-1) - rad_st_old(io, j)),  &
       &   creep_strain_new(i, j)%rad - creep_strain_old(io, j)%rad, & 
       &   creep_strain_new(i, j)%hoop - creep_strain_old(io, j)%hoop, & 
       &   creep_strain_new(i, j)%ax - creep_strain_old(io, j)%ax 
 2         format('ox_steam_str_cr_corr ', i4, 20(1pe13.6, 1x))

        else if (ABS(rad_now - rad_st_old(io, npr_st(io)))  < &
             & displ_max * 1.0e-6) then

          ip_now(i, j) = npr_st(i) + 1
          il_now(i, j) = io
          ! keep the old creep strains since identical nodes
          creep_strain_new(i, j)%rad = creep_strain_old(i, j)%rad
          creep_strain_new(i, j)%hoop = creep_strain_old(i, j)%hoop
          creep_strain_new(i, j)%ax = creep_strain_old(i, j)%ax

        endif

      endif

      ! if (il_now(i, j) > 0)  EXIT LAY2

      ! end do LAY2

    ! end do LOC1

  end do LAY1

  ! creep_strain_new was used as dummy here
  creep_strain_old = creep_strain_new

  ! obtain the new creep strain

  ! stop

  return

  END SUBROUTINE CREEP_INIT_LOCATION

  SUBROUTINE CREEP_CONV(iter_st_creep, no_creep_intervals, N, &
      & cycle_no, time_now, dt_now, id_high_or_low, convergent)

    !=======================================================================
    ! Purpose(s):
    !
    ! stress_eq_epsilon = a small number to test the 
    ! convergence of 
    ! SUM |stress_eq - stress_eq_old| < stress_eq_epsilon * stress_eq
    !=======================================================================
  use parameter_module, only: moxide, mave, oxide_thick2si, &
       & id_high, id_low, mgrid
  use creep_data_module, only: creep_strain_old, creep_strain_new, &
       & creep_sum, creep_dif, creep_integral_sum, creep_integral_dif, &
       & if_creep_ox_scale, if_creep_metal, stress_eq, stress_eq_old, &
       & i_cr_start, i_cr_end, stress_eq_epsilon, creep_strain_ave, &
       & stress_eq_ave, creep_strain_eq, creep_strain_incr_eq, &
       & creep_strain_eq_ave, stress_eq_conv, operation_2_initial_id
  use creep_data_module, only: no_intervals_now, operation_id_current, &
       & interval_op_now
  use oxide_data_module,   only: no_oxide, no_layer, &
       & no_oxide_ave, cond_value, first_oxide_layer
  use solution_data_module, only: npr_temp, npr_st, Temp, rad_int, &
       & stress_axial, stress_hoop, stress_rad, eps0, cc_st, bc_st
  use property_module, only: OPERAND_ARRAY
  use solver_data_module, only  : max_iterat_solver
  use output_module,      only: tty_lun, out_lun, aux_lun

  ! Arguments
  integer, intent(IN)  :: iter_st_creep, no_creep_intervals
  integer, intent(IN)  :: id_high_or_low, cycle_no
  integer, intent(IN)  :: N
  real, intent(IN)     :: time_now, dt_now
  logical, intent(OUT)  :: convergent

  ! local variables
  real    :: norm_lhs, norm_rhs, local_norm_lhs, local_norm_rhs, norm_local
  integer  :: i, j, n_norm, n_local_norm, loc1, k, j_temp
  logical  :: debug_creep_local = .true.
  character(LEN = 80), dimension(mave), save    :: operation_type

  if (i_cr_start > i_cr_end)  then
    ! there is no creep
    convergent = .true.
    RETURN
  endif

  if ((.not. if_creep_metal) .and. (.not. if_creep_ox_scale))  then
    ! there is no creep
    convergent = .true.
    RETURN
  endif

 ! update von Mises stress
  do i = i_cr_start, i_cr_end

    do j = 1, npr_st(i)

      ! get the equivalent Von Mises stress
      stress_eq(i, j) = VON_MISES(stress_rad(i, j), stress_hoop(i, j), &
               & stress_axial(i, j))
    end do

  end do

  local_norm_lhs = 0.0
  local_norm_rhs = 0.0
  n_local_norm = 0
  
  ! evaluate the equivalent creep strain
  do i = i_cr_start, i_cr_end

    do j = 1, npr_st(i)

      if (stress_eq(i, j) > 0.0)  then

        n_local_norm = n_local_norm + 1

        local_norm_lhs = local_norm_lhs + (1.0 - stress_eq_old(i, j) / &
                 & stress_eq(i, j))**2

      endif

    end do

  end do

  if (n_local_norm > 0)  then

    ! norm_rhs = stress_eq_epsilon
    local_norm_lhs = SQRT(local_norm_lhs / n_local_norm)

  endif

  norm_lhs = 0.0
  norm_rhs = 0.0
  n_norm = 0
  
  ! evaluate the equivalent creep strain
  do i = i_cr_start, i_cr_end

    do j = 1, npr_st(i)

      if (stress_eq(i, j) > 0.0)  then

        n_norm = n_norm + 1

        norm_lhs = norm_lhs + (1.0 - stress_eq_old(i, j) / &
                 & stress_eq(i, j))**2

      endif

    end do

  end do

  if (n_norm > 0)  then

    ! norm_rhs = stress_eq_epsilon
    norm_lhs = SQRT(norm_lhs / n_norm)

    if (norm_lhs < stress_eq_epsilon)  then  !  norm_rhs)  then

      convergent = .true.

    else

      convergent = .false.

    endif

  else

    ! write(2, *) 'ERROR: no von Mises points for convergence norm'
    ! STOP
    write(2, 6) time_now, iter_st_creep, no_creep_intervals
 6    format('stress_eq_zero ', 1pe13.6, 2(1x, i2), 1x, 50(1pe13.6, 1x))

    if (iter_st_creep == 1)  then
      ! force one more iteration since von Mises stress may be zero initially
      convergent = .false.
    else
      convergent = .true.
    endif

    ! RETURN

  endif

  if (convergent)  then

    if (n_norm > 0) then

      ! convergence was attained; update the old creep
      creep_strain_old%rad = creep_strain_new%rad 
      creep_strain_old%hoop = creep_strain_new%hoop 
      creep_strain_old%ax = creep_strain_new%ax 

      ! this was enough when the only 2-3 and 4-1 was considered
      ! stress_eq_conv(id_high_or_low, :, :) = stress_eq
      if (operation_id_current == 3 .or. &
         & operation_id_current == 7)  then

        ! store the vm stress at full load only at the last interval
        if (interval_op_now == no_intervals_now(operation_id_current))  then
          stress_eq_conv(id_high, :, :) = stress_eq
        endif

      else if (operation_id_current == 1 .or. &
         & operation_id_current == 5)  then

        ! store the vm stress at partial load only at the last interval
        if (interval_op_now == no_intervals_now(operation_id_current))  then
          stress_eq_conv(id_low, :, :) = stress_eq
        endif

      else if (operation_id_current == 4 .or. &  ! f2l transition
         & operation_id_current == 8)  then

        ! store the vm stress at partial load only at the last interval
        ! do k = 1, no_intervals_now(operation_id_current)
        ! if (interval_op_now == no_intervals_now(operation_id_current))  then
          ! do not assigned anything
          !  stress_eq_conv(id_low, :, :) = stress_eq
        ! endif

        if (interval_op_now <= no_intervals_now(operation_id_current)-1)  then
          stress_eq_conv(2+interval_op_now, :, :) = stress_eq
        endif

      ! do not assign the von mises at the ramp up since these will be used
      ! in reverse order

      endif

      ! increment the equivalent creep strain only when convergence was attained
      creep_strain_eq = creep_strain_eq + creep_strain_incr_eq

    else
      ! this was enough when the only 2-3 and 4-1 was considered
      ! stress_eq_conv(id_high_or_low, :, :) = 0.0
      stress_eq_conv = 0.0

    endif

  endif

  call OUT_CREEP_CONV(iter_st_creep, no_creep_intervals, N, &
      & cycle_no, time_now, dt_now, id_high_or_low, &
      & norm_lhs, norm_rhs, n_norm, convergent)

  RETURN

  END SUBROUTINE CREEP_CONV

  SUBROUTINE OUT_CREEP_CONV(iter_st_creep, no_creep_intervals, N, &
      & cycle_no, time_now, dt_now, id_high_or_low, &
      & norm_lhs, norm_rhs, n_norm, convergent)

    !=======================================================================
    ! Purpose(s):
    !
    ! stress_eq_epsilon = a small number to test the 
    ! convergence of 
    ! SUM |stress_eq - stress_eq_old| < stress_eq_epsilon * stress_eq
    !=======================================================================
  use parameter_module, only: moxide, mave, oxide_thick2si, &
       & id_high, id_low, mgrid
  use creep_data_module, only: creep_strain_old, creep_strain_new, &
       & creep_sum, creep_dif, creep_integral_sum, creep_integral_dif, &
       & if_creep_ox_scale, if_creep_metal, stress_eq, stress_eq_old, &
       & i_cr_start, i_cr_end, stress_eq_epsilon, creep_strain_ave, &
       & stress_eq_ave, creep_strain_eq, creep_strain_incr_eq, &
       & creep_strain_eq_ave, stress_eq_conv, operation_2_initial_id
  use creep_data_module, only: no_intervals_now, operation_id_current, &
       & interval_op_now
  use oxide_data_module,   only: no_oxide, no_layer, &
       & no_oxide_ave, cond_value, first_oxide_layer
  use solution_data_module, only: npr_temp, npr_st, Temp, rad_int, &
       & stress_axial, stress_hoop, stress_rad, eps0, cc_st, bc_st, &
       & n_mean_odd_vm, rad_st, n_hf_out
  use property_module, only: OPERAND_ARRAY
  use solver_data_module, only  : max_iterat_solver
  use output_module,      only: tty_lun, out_lun, aux_lun
  use boiler_data_module, only: Time_Boiler_p, TEMP_BOILER_data_no

  ! Arguments
  integer, intent(IN)  :: iter_st_creep, no_creep_intervals
  integer, intent(IN)  :: id_high_or_low, cycle_no
  integer, intent(IN)  :: N, n_norm
  real, intent(IN)     :: time_now, norm_lhs, norm_rhs, dt_now
  logical, intent(IN)  :: convergent

  ! local variables
  real    ::  local_norm_lhs, local_norm_rhs, norm_local
  integer  :: i, j, n_local_norm, loc1, k, j_temp
  logical  :: debug_creep_local = .true.
  character(LEN = 80), dimension(mave), save    :: operation_type
  real, dimension(N)  :: stress_vm_mean, temp_cr_mean

  integer, save :: n_sh_vm_prof = 0
  real, save    :: dtime_sh_vm_prof, time_sh_vm_prof, time_previous = 0.0
  real, save, dimension(20, mgrid)  :: vm_prof, strain_hoop_prof
  real, save, dimension(20)         :: time_sh_vm
  real, dimension(13)  :: time_sh_vm_zarrabi = (/0.0, 5.0, 10.0, 50.0, 100.0, &
      & 200.0, 300.0, 400.0, 500.0, 1000.0, 2000.0, 3000.0, 4000.0/)

  ! initialize the operation type
  operation_type(1) = 'max'
  operation_type(2) = 'min'
  operation_type(3) = 'ave_sum'

  do i = i_cr_start, i_cr_end

      j = npr_st(i)  !  location where stress/strains are maximum in metal

      ! temperature point 
      j_temp = 2 * j - 1

      if (stress_eq(i, j) > 0.0)  then

        norm_local = DABS(1.0 - stress_eq_old(i, j) / stress_eq(i, j))

      else

        norm_local = -1.0

      endif

      write(2, 13) i, iter_st_creep, no_creep_intervals, &
         & time_now, stress_eq(i, j), stress_eq_old(i, j), norm_local, &
         & Temp(i, j_temp), creep_strain_new(i, j)%rad, &
         & creep_strain_new(i, j)%rad - creep_strain_old(i, j)%rad
 13    format('metal_creep_conv ', i1, 1x, i4, 1x, i2, 1x, 50(1pe13.6, 1x))
      
      if (stress_eq(i, j) > 1.0e+8)  then
         write(6, 15) time_now, time_now / 24.0, time_now / &
           & (24.0*365.0), (rad_int(first_oxide_layer) - rad_int(N+1)) / &
           & oxide_thick2si, cycle_no
 15      format('ERROR: divergent; too much creep at t= ', &
             & 1pe13.6, ' h; ', 1pe13.6, ' days;', 1pe13.6, &
             & ' years;', 1pe13.6, ' [micron]', i6, '= cycle')
         stop
      endif

      do k = 1, 3

        ! get maximum, minimum, and average creep in each layer
        call OPERAND_ARRAY(1, mgrid, creep_strain_new(i, :)%rad,  &
              & 1, npr_st(i), operation_type(k), &
              & creep_strain_ave(i, k)%rad, loc1)
        call OPERAND_ARRAY(1, mgrid, creep_strain_new(i, :)%hoop,  &
              & 1, npr_st(i), operation_type(k), &
              & creep_strain_ave(i, k)%hoop, loc1)
        call OPERAND_ARRAY(1, mgrid, creep_strain_new(i, :)%ax,  &
              & 1, npr_st(i), operation_type(k), &
              & creep_strain_ave(i, k)%ax, loc1)

        call OPERAND_ARRAY(1, mgrid, stress_eq(i, :),  &
              & 1, npr_st(i), operation_type(k), &
              & stress_eq_ave(i, k), loc1)

      end do

  end do

  if (convergent)  then

    if (n_norm > 0) then

      write(2, 7) time_now, &
         & (rad_int(first_oxide_layer) - rad_int(N+1)) / oxide_thick2si, &
         & iter_st_creep, no_creep_intervals, norm_lhs, norm_rhs, &
         & (stress_eq_ave(i, 2), stress_eq_ave(i, 3), &
         & stress_eq_ave(i, 1), i = i_cr_start, i_cr_end), &
         & (Temp(i, 1), Temp(i, npr_temp(i)), i = i_cr_start, i_cr_end)
 7    format('creep_conv_norm ', 2(1pe13.6, 1x), 2(i2, 1x), &
         & 70(1pe13.6, 1x))

      do i = 1, N
        stress_vm_mean(i) = VON_MISES(stress_rad(i, n_mean_odd_vm(i)), &
            & stress_hoop(i, n_mean_odd_vm(i)), &
            & stress_axial(i, n_mean_odd_vm(i)))
        temp_cr_mean(i) = Temp(i,  npr_st(i))
      end do

       write(2, 27) time_now, & 
         & (rad_int(first_oxide_layer) - rad_int(N+1)) / oxide_thick2si, &
         & iter_st_creep, no_creep_intervals, norm_lhs, &
         & (stress_vm_mean(i), creep_strain_new(i, n_mean_odd_vm(i))%hoop, &
         & Temp(i,  npr_st(i)), i = i_cr_start, i_cr_end)
 27    format('conv_norm_mean1 ', 2(1pe13.6, 1x), 2(i2, 1x), &
         & 70(1pe13.6, 1x))

      ! write(6, *) 'time shs ', n_sh_vm_prof, time_now -dt_now, &
      ! & time_sh_vm_prof, time_now

      do k = 1, 12

        if (time_now > time_sh_vm_zarrabi(k).and. &
         & time_previous <= time_sh_vm_zarrabi(k))  then

          n_sh_vm_prof = n_sh_vm_prof + 1
          time_sh_vm(n_sh_vm_prof) = time_now

          vm_prof(n_sh_vm_prof, 1:npr_st(1)) = stress_eq(1, 1:npr_st(1))
          strain_hoop_prof(n_sh_vm_prof, 1:npr_st(1)) = &
             & stress_hoop(1, 1:npr_st(1))

          ! write(tty_lun, *) 'zarrabi out ', n_sh_vm_prof, time_now

      ! write everything out at the last output time
      if (n_sh_vm_prof == 11)  then  !  n_hf_out - 1)  then

        write(tty_lun, 51) (time_sh_vm(i), i = 1, 11) ! n_hf_out - 1)
        write(out_lun, 52) (INT(time_sh_vm(i)), INT(time_sh_vm(i)), i = 1, 11)
 52     format('sh_vm_prof radius_cm ', 100(i4, 1x))
        do j = 1, npr_st(1)
            write(out_lun, 51) rad_st(1, j)*100.0, (vm_prof(i, j), &
              & strain_hoop_prof(i, j), i = 1, 11) ! n_hf_out - 1)
 51       format('sh_vm_prof ', 100(1pe13.6, 1x))
        end do

      endif

        endif

      end do

      time_previous = time_now

      do i = i_cr_start, i_cr_end

        do k = 1, 3

          ! get maximum, minimum, and average equivalent creep in each layer
          call OPERAND_ARRAY(1, mgrid, creep_strain_eq(i, :),  &
              & 1, npr_st(i), operation_type(k), &
              & creep_strain_eq_ave(i, k), loc1)

        end do

      end do

      write(aux_lun, 17) cycle_no, time_now, eps0, &
         & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
         & (bc_st(i), i= first_oxide_layer, N), &
         & (cc_st(i), i= first_oxide_layer, N)
 17    format('const_conv_sol_e_b=', i5, 1x, 30(1pe13.6, 1x))

    else

      write(2, 11) time_now, iter_st_creep, time_now - time_now, &
         & time_now - time_now, &
         & (stress_eq_ave(i, 2), stress_eq_ave(i, 3), &
         & stress_eq_ave(i, 1), i = i_cr_start, i_cr_end)
 11    format('creep_conv_no_norm ', 1pe13.6, 1x, i2, 1x, 70(1pe13.6, 1x))

    endif

      write(2, 8) time_now, iter_st_creep, no_creep_intervals, &
        & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
        & (creep_strain_ave(i, 1)%rad*1.0e+3, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 1)%hoop*1.0e+3, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 1)%ax*1.0e+3, i = i_cr_start, i_cr_end)

 8    format('creep_conv_max ', 1pe13.6, 1x, 2(i2, 1x), 50(1pe13.6, 1x))

      write(2, 9) time_now, iter_st_creep, no_creep_intervals, &
        & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
        & (creep_strain_ave(i, 2)%rad*1.0e+3, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 2)%hoop*1.0e+3, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 2)%ax*1.0e+3, i = i_cr_start, i_cr_end)

 9    format('creep_conv_min ', 1pe13.6, 1x, 2(i2, 1x), 50(1pe13.6, 1x))

      write(2, 10) time_now, iter_st_creep, no_creep_intervals, &
        & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
        & (creep_strain_ave(i, 3)%rad*1.0e+3, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 3)%hoop*1.0e+3, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 3)%ax*1.0e+3, i = i_cr_start, i_cr_end)

 10    format('creep_conv_ave ', 1pe13.6, 1x, 2(i2, 1x), 50(1pe13.6, 1x))

  else

    if (debug_creep_local)  then

    if (n_norm > 0) then

      write(2, 2) time_now, iter_st_creep, no_creep_intervals, &
         & norm_lhs, norm_rhs, &
         & (stress_eq_ave(i, 2), stress_eq_ave(i, 3), &
         & stress_eq_ave(i, 1), i = i_cr_start, i_cr_end)
 2    format('creep_iter_norm ', 1pe13.6, 1x, 2(i2, 1x), 50(1pe13.6, 1x))

    else

      write(2, 12) time_now, iter_st_creep, time_now - time_now, &
         & time_now - time_now, &
         & (stress_eq_ave(i, 2), stress_eq_ave(i, 3), &
         & stress_eq_ave(i, 1), i = i_cr_start, i_cr_end)
 12    format('creep_iter_no_norm ', 1pe13.6, 1x, i2, 1x, 50(1pe13.6, 1x))

    endif

      write(2, 3) time_now, iter_st_creep, &
        & (creep_strain_ave(i, 1)%rad, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 1)%hoop, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 1)%ax, i = i_cr_start, i_cr_end)

 3    format('creep_iter_max ', 1pe13.6, 1x, i2, 1x, 50(1pe13.6, 1x))

      write(2, 4) time_now, iter_st_creep, &
        & (creep_strain_ave(i, 2)%rad, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 2)%hoop, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 2)%ax, i = i_cr_start, i_cr_end)

 4    format('creep_iter_min ', 1pe13.6, 1x, i2, 1x, 50(1pe13.6, 1x))

      write(2, 5) time_now, iter_st_creep, &
        & (creep_strain_ave(i, 3)%rad, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 3)%hoop, i = i_cr_start, i_cr_end), &
        & (creep_strain_ave(i, 3)%ax, i = i_cr_start, i_cr_end)

 5    format('creep_iter_ave ', 1pe13.6, 1x, i2, 1x, 50(1pe13.6, 1x))

    endif

  endif
  
  if (.not. convergent .and. iter_st_creep >= max_iterat_solver)  then
    write(out_lun, 16) time_now, time_now / 24.0, time_now / &
        & (24.0*365.0), (rad_int(first_oxide_layer) - rad_int(N+1)) / &
        & oxide_thick2si, cycle_no
 16 format('ERROR: stress iterations > max_iterat_solver at t= ', &
             & 1pe13.6, ' h; ', 1pe13.6, ' days;', 1pe13.6, &
             & ' years;', 1pe13.6, ' [micron]', i6, '= cycle')
    write(tty_lun, 16) time_now, time_now / 24.0, time_now / &
        & (24.0*365.0), (rad_int(first_oxide_layer) - rad_int(N+1)) / &
        & oxide_thick2si, cycle_no
    STOP

  endif

  RETURN

  END SUBROUTINE OUT_CREEP_CONV

  SUBROUTINE INIT_CREEP(N)

    !=======================================================================
    ! Purpose(s):
    !
    !   get the creep strain and integral of creep strain needed 
    !       for u, e_r, s_r, s_hoop, s_z
    !   if average properties are needed for the scale just run 
    !     another case indicating in oxide_layer namelist the scale compounds 
    !   creep routine will break down if if_ave_scale = .true.
    !=======================================================================
 
  use creep_data_module, only: creep_strain_old, creep_strain_new, &
       & creep_sum, creep_dif, creep_integral_sum, creep_integral_dif, &
       & if_creep_ox_scale, if_creep_metal, creep_strain_eq, &
       & i_cr_start, i_cr_end
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value, no_metal_layers, first_oxide_layer
  ! use solver_data_module, only: var, var_conv, update_sola_type
  ! use solution_data_module, only: ib, ic, i_cr_r, i_cr_h, i_cr_z, i_s_r, &
  !    & i_s_h, i_s_z, i_s_vm, i_cr, ieps0, no_stress_var

  ! Arguments
  integer, intent(IN)  :: N
 
  ! local variables
  logical, save :: first_creep_strain = .true.
  real    :: Delta_int, stress_p
  integer  :: i, j, j_temp
  logical    :: debug_local = .false.

  FIRST_IF: if (first_creep_strain)  then

    ! intialize
    creep_strain_old%rad = 0.0
    creep_strain_old%hoop = 0.0
    creep_strain_old%ax = 0.0

    creep_strain_eq = 0.0

    first_creep_strain = .false.

    if (if_creep_metal)  then

      i_cr_start = 1

    else

      creep_strain_new(1:no_metal_layers,:)%rad = 0.0
      creep_strain_new(1:no_metal_layers,:)%hoop = 0.0
      creep_strain_new(1:no_metal_layers,:)%ax = 0.0

      creep_integral_sum(1:no_metal_layers,:) = 0.0
      creep_integral_dif(1:no_metal_layers,:) = 0.0

      i_cr_start = first_oxide_layer

    endif

    if (if_creep_ox_scale)  then

      i_cr_end = N

    else

      creep_strain_new(first_oxide_layer:N, :)%rad = 0.0
      creep_strain_new(first_oxide_layer:N, :)%hoop = 0.0
      creep_strain_new(first_oxide_layer:N, :)%ax = 0.0

      creep_integral_sum(first_oxide_layer:N, :) = 0.0
      creep_integral_dif(first_oxide_layer:N, :) = 0.0

      i_cr_end = 1

    endif 

    ! this initialization is required for the first iteration since there is 
    ! a return of iter = 1
    if (i_cr_start > i_cr_end)  then
      creep_integral_sum(i_cr_start:i_cr_end, :) = 0.0
      creep_integral_dif(i_cr_start:i_cr_end, :) = 0.0
    else
      ! there is no creep in any layer
      creep_integral_sum = 0.0
      creep_integral_dif = 0.0
    endif

   endif FIRST_IF 

  RETURN

  END SUBROUTINE INIT_CREEP

  SUBROUTINE CREEP_STRAIN(iter_st_creep, cycle_no, N, id_high_or_low, &
       & dt, time_now)

    !=======================================================================
    ! Purpose(s):
    !
    !   get the creep strain and integral of creep strain needed 
    !       for u, e_r, s_r, s_hoop, s_z
    !   if average properties are needed for the scale just run 
    !     another case indicating in oxide_layer namelist the scale compounds 
    !   creep routine will break down if if_ave_scale = .true.
    !=======================================================================
 
  use parameter_module, only: id_low, id_high
  use creep_data_module, only: creep_strain_old, creep_strain_new, &
       & creep_sum, creep_dif, creep_integral_sum, creep_integral_dif, &
       & if_creep_ox_scale, if_creep_metal, stress_eq, stress_eq_old, &
       & creep_strain_incr_eq, i_cr_start, i_cr_end, stress_eq_conv, &
       & no_intervals_now, operation_id_current, &
      & interval_op_now, creep_integral_z
  use solution_data_module, only: rad_temp, Temp, npr_temp, &
      & npr_st, rad_st, eps0, rad_int, strain_rad, strain_hoop, &
      & displ_rad, stress_axial, stress_hoop, stress_rad, &
      & strain_th_ave, tau_th_st, &
      & strain_hoop_gen, n_mean_odd_vm
  use oxide_data_module,   only: poisson_ratio, no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use boiler_data_module, only: no_pulse_per_cycle

  ! Arguments
  integer, intent(IN)  :: iter_st_creep
  integer, intent(IN)  :: N
  integer, intent(IN)  :: id_high_or_low, cycle_no
  real, intent(IN)     :: dt, time_now

  ! local variables
  real    :: Delta_int, creep_strain_incr_factor, stress_p, stress_vm_new, &
           & stress_vm_mean, temp_cr_mean
  real    :: vm_relax_factor, zarrabi_creep_factor, exp1_zarrabi
  integer  :: i, j, j_temp, id
  logical    :: debug_local = .false., if_creep_mean_radius

   if (iter_st_creep == 1)  then

     ! for the first stress iteration use the old creep
     creep_strain_new%rad = creep_strain_old%rad 

     creep_strain_new%hoop = creep_strain_old%hoop 

     creep_strain_new%ax = creep_strain_old%ax 

     ! just to force at least two creep iterations 
     ! since criteria of convergence will depend on |stress_eq - stress_eq_old|
     ! also the von Mises for full and partial load are different
     ! stress_eq_old = 0.0
     
     if (no_pulse_per_cycle == 1 .and. cycle_no <= 4)  then 

       if (id_high_or_low == id_low)  then
         ! previous stress level was higher
         stress_eq_old = stress_eq / 3.3
       else if (id_high_or_low == id_high)  then
         ! previous stress level was low
         stress_eq_old = stress_eq / 1.3
       endif

       RETURN

     endif

     if (no_pulse_per_cycle == 2 .and. cycle_no <= 2)  then 

       if (id_high_or_low == id_low)  then
         ! previous stress level was higher
         stress_eq_old = stress_eq / 3.3
       else if (id_high_or_low == id_high)  then
         ! previous stress level was low
         stress_eq_old = stress_eq / 1.3
       endif

       RETURN

     endif

     ! the following is for larger, regular cycles i.e., cycle_no > 2(4)

     if (operation_id_current == 3 .or. &
       & operation_id_current == 7)  then

       stress_eq_old = stress_eq_conv(id_high_or_low, :, :)

     else if (operation_id_current == 1 .or. &
       & operation_id_current == 5)  then

       stress_eq_old = stress_eq_conv(id_high_or_low, :, :)

     else if (operation_id_current == 4 .or. &  ! f2l transition
         & operation_id_current == 8)  then

       if (interval_op_now == &
           & no_intervals_now(operation_id_current))  then
         stress_eq_old = stress_eq_conv(id_low, :, :)
       else
         stress_eq_old = stress_eq_conv(2+interval_op_now, :, :)
       endif

     else if (operation_id_current == 2 .or. &  ! f2l transition
         & operation_id_current == 6)  then

       if (interval_op_now == &
           & no_intervals_now(operation_id_current))  then
         stress_eq_old = stress_eq_conv(id_high, :, :)
       else
         id = 2 + no_intervals_now(operation_id_current) - interval_op_now
         stress_eq_old = stress_eq_conv(id, :, :)
       endif

     endif

     ! creep_integral_sum and creep_integral_dif remain unchanged
     ! make sure they are not reset

     RETURN

   else if (iter_st_creep > 1)  then

     ! prepare for the next iteration
     stress_eq_old = stress_eq

   endif

  vm_relax_factor = 0.0
  if_creep_mean_radius = .false.

  ! evaluate the equivalent creep strain
  do i = i_cr_start, i_cr_end

    stress_vm_mean = VON_MISES(stress_rad(i, n_mean_odd_vm(i)), &
       & stress_hoop(i, n_mean_odd_vm(i)), stress_axial(i, n_mean_odd_vm(i)))
    temp_cr_mean = Temp(i,  npr_st(i))

    do j = 1, npr_st(i)

      ! temperature point 
      j_temp = 2 * j - 1

        ! get the equivalent Von Mises stress
        stress_vm_new = VON_MISES(stress_rad(i, j), stress_hoop(i, j), &
               & stress_axial(i, j))
        stress_eq(i, j) = vm_relax_factor * stress_eq_old(i, j) + &
               & (1.0 - vm_relax_factor) * stress_vm_new
    
        if (stress_eq(i, j) > 0.0)  then

          zarrabi_creep_factor = 1.0
          ! exp1_zarrabi = (367.79-stress_eq(i, j)) / 57.774
          ! zarrabi_creep_factor = 1.0 + exp(time_now - 10.0**exp1_zarrabi)

          ! equivalet creep increment
          creep_strain_incr_eq(i, j) = DCREEP(stress_eq(i, j), dt, &
               & Temp(i, j_temp), i) * zarrabi_creep_factor
          
          if (if_creep_mean_radius)  then
            creep_strain_incr_eq(i, j) = DCREEP(stress_vm_mean, dt, &
               & temp_cr_mean, i)
          endif

          ! add new constants
          creep_strain_incr_factor = 1.5 * creep_strain_incr_eq(i, j) / &
               & stress_eq(i, j)

          ! isostatic stress
          stress_p = (stress_rad(i, j) + stress_hoop(i, j) + &
               & stress_axial(i, j)) / 3.0

          ! new creep strain
          creep_strain_new(i, j)%rad = creep_strain_old(i, j)%rad + &
            & (stress_rad(i, j) - stress_p) * creep_strain_incr_factor

          creep_strain_new(i, j)%hoop = creep_strain_old(i, j)%hoop + &
            & (stress_hoop(i, j) - stress_p) * creep_strain_incr_factor

          creep_strain_new(i, j)%ax = creep_strain_old(i, j)%ax + &
            & (stress_axial(i, j) - stress_p) * creep_strain_incr_factor

        else

          ! Von Mises was set to zero for first iteration
          creep_strain_new(i, j)%rad = creep_strain_old(i, j)%rad
          creep_strain_new(i, j)%hoop = creep_strain_old(i, j)%hoop
          creep_strain_new(i, j)%ax = creep_strain_old(i, j)%ax

        endif

    end do

    ! get creep sum; needed for integration
    creep_sum(i, :) = creep_strain_new(i, :)%rad + &
       & creep_strain_new(i, :)%hoop + &
       & 2.0 * poisson_ratio(i) * creep_strain_new(i, :)%ax

    ! get creep difference
    creep_dif(i, :) = creep_strain_new(i, :)%rad - &
       & creep_strain_new(i, :)%hoop

  end do

  ! I_i = int(low limit r_i+1, upper limit r) (F * x) dx 
  !      since r_i+1 < r
  
  ! indices are backwards, i.e., lower indices point to larger radius

  do i = i_cr_start, i_cr_end

    creep_integral_sum(i, npr_st(i)) = 0.0
    creep_integral_dif(i, npr_st(i)) = 0.0

    do j = npr_st(i) - 1, 1, -1

      ! temperature point 
      j_temp = 2 * j

      ! fractional integral; F * r * dr
      ! the trapezium integral should be (a+b)*h/2
      Delta_int = 0.5 * (creep_sum(i, j) * rad_st(i, j) + &
                & creep_sum(i, j+1) * rad_st(i, j+1)) * &
                & (rad_st(i, j) - rad_st(i, j+1))

      ! integral F*x dx
      creep_integral_sum(i, j) = creep_integral_sum(i, j + 1) + Delta_int

      ! fractional integral; (F / r) * dr
      ! the trapezium integral should be (a+b)*h/2
      Delta_int = 0.5 * (creep_dif(i, j) / rad_st(i, j) + &
                & creep_dif(i, j+1) / rad_st(i, j+1)) * &
                & (rad_st(i, j) - rad_st(i, j+1))

      ! integral (F / r) dr
      creep_integral_dif(i, j) = creep_integral_dif(i, j + 1) + Delta_int

      if (debug_local)  then

        ! write(aux_lun, 33) i, j, j_temp, rad_st(i, j), rad_st(i, j+1), &
        ! & rad_st(i, j) - rad_st(i, j+1), rad_temp(i, j_temp), &
        ! & tau_th(i, j_temp), DI_tau_r_dr
 33     format('dintegral ', i1, 2(1x, i2), 1x, 10(1pe13.6, 1x))

      endif

    end do

  end do

  ! obtain the creep strain integral in the z direction
  do i = i_cr_start, i_cr_end

    creep_integral_z(i) = 0.0

    do j = npr_st(i) - 1, 1, -1

      ! fractional integral; F * r * dr
      ! the trapezium integral should be (a+b)*h/2
      Delta_int = 0.5 * (creep_strain_new(i, j)%ax * rad_st(i, j) + &
           & creep_strain_new(i, j+1)%ax * rad_st(i, j+1)) * &
           & (rad_st(i, j) - rad_st(i, j+1))

      ! integral F*x dx
      creep_integral_z(i) = creep_integral_z(i) + Delta_int

    end do

  end do

  RETURN

  END SUBROUTINE CREEP_STRAIN

  SUBROUTINE CREEP_STRAIN_MEAN(iter_st_creep, cycle_no, N, &
      & id_high_or_low, dt)

    !=======================================================================
    ! Purpose(s):
    !   obtain the mean creep 
    !   get the creep strain and integral of creep strain needed 
    !       for u, e_r, s_r, s_hoop, s_z
    !   if average properties are needed for the scale just run 
    !     another case indicating in oxide_layer namelist the scale compounds 
    !   creep routine will break down if if_ave_scale = .true.
    !=======================================================================
 
  use parameter_module, only: id_low, id_high
  use creep_data_module, only: creep_strain_old, creep_strain_new, &
       & creep_sum, creep_dif, creep_integral_sum, creep_integral_dif, &
       & if_creep_ox_scale, if_creep_metal, stress_eq, stress_eq_old, &
       & creep_strain_incr_eq, i_cr_start, i_cr_end, stress_eq_conv, &
       & no_intervals_now, operation_id_current, &
      & interval_op_now, creep_integral_z
  use solution_data_module, only: rad_temp, Temp, npr_temp, &
      & npr_st, rad_st, eps0, rad_int, strain_rad, strain_hoop, &
      & displ_rad, stress_axial, stress_hoop, stress_rad, &
      & strain_th_ave, tau_th_st, &
      & strain_hoop_gen
  use oxide_data_module,   only: poisson_ratio, no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use boiler_data_module, only: no_pulse_per_cycle

  ! Arguments
  integer, intent(IN)  :: iter_st_creep
  integer, intent(IN)  :: N
  integer, intent(IN)  :: id_high_or_low, cycle_no
  real, intent(IN)     :: dt

  ! local variables
  real    :: Delta_int, creep_strain_incr_factor, stress_p, stress_vm_new
  real    :: vm_relax_factor
  integer  :: i, j, j_temp, id
  logical    :: debug_local = .false.

   if (iter_st_creep == 1)  then

     ! for the first stress iteration use the old creep
     creep_strain_new%rad = creep_strain_old%rad 

     creep_strain_new%hoop = creep_strain_old%hoop 

     creep_strain_new%ax = creep_strain_old%ax 

     ! just to force at least two creep iterations 
     ! since criteria of convergence will depend on |stress_eq - stress_eq_old|
     ! also the von Mises for full and partial load are different
     ! stress_eq_old = 0.0
     
     if (no_pulse_per_cycle == 1 .and. cycle_no <= 4)  then 

       if (id_high_or_low == id_low)  then
         ! previous stress level was higher
         stress_eq_old = stress_eq / 3.3
       else if (id_high_or_low == id_high)  then
         ! previous stress level was low
         stress_eq_old = stress_eq / 1.3
       endif

       RETURN

     endif

     if (no_pulse_per_cycle == 2 .and. cycle_no <= 2)  then 

       if (id_high_or_low == id_low)  then
         ! previous stress level was higher
         stress_eq_old = stress_eq / 3.3
       else if (id_high_or_low == id_high)  then
         ! previous stress level was low
         stress_eq_old = stress_eq / 1.3
       endif

       RETURN

     endif

     ! the following is for larger, regular cycles i.e., cycle_no > 2(4)

     if (operation_id_current == 3 .or. &
       & operation_id_current == 7)  then

       stress_eq_old = stress_eq_conv(id_high_or_low, :, :)

     else if (operation_id_current == 1 .or. &
       & operation_id_current == 5)  then

       stress_eq_old = stress_eq_conv(id_high_or_low, :, :)

     else if (operation_id_current == 4 .or. &  ! f2l transition
         & operation_id_current == 8)  then

       if (interval_op_now == &
           & no_intervals_now(operation_id_current))  then
         stress_eq_old = stress_eq_conv(id_low, :, :)
       else
         stress_eq_old = stress_eq_conv(2+interval_op_now, :, :)
       endif

     else if (operation_id_current == 2 .or. &  ! f2l transition
         & operation_id_current == 6)  then

       if (interval_op_now == &
           & no_intervals_now(operation_id_current))  then
         stress_eq_old = stress_eq_conv(id_high, :, :)
       else
         id = 2 + no_intervals_now(operation_id_current) - interval_op_now
         stress_eq_old = stress_eq_conv(id, :, :)
       endif

     endif

     ! creep_integral_sum and creep_integral_dif remain unchanged
     ! make sure they are not reset

     RETURN

   else if (iter_st_creep > 1)  then

     ! prepare for the next iteration
     stress_eq_old = stress_eq

   endif

  vm_relax_factor = 0.0

  ! evaluate the equivalent creep strain
  do i = i_cr_start, i_cr_end

    do j = 1, npr_st(i)

      ! temperature point 
      j_temp = 2 * j - 1

        ! get the equivalent Von Mises stress
        stress_vm_new = VON_MISES(stress_rad(i, j), stress_hoop(i, j), &
               & stress_axial(i, j))
        stress_eq(i, j) = vm_relax_factor * stress_eq_old(i, j) + &
               & (1.0 - vm_relax_factor) * stress_vm_new
    
        if (stress_eq(i, j) > 0.0)  then

          ! equivalet creep increment
          creep_strain_incr_eq(i, j) = DCREEP(stress_eq(i, j), dt, &
               & Temp(i, j_temp), i)

          ! add new constants
          creep_strain_incr_factor = 1.5 * creep_strain_incr_eq(i, j) / &
               & stress_eq(i, j)

          ! isostatic stress
          stress_p = (stress_rad(i, j) + stress_hoop(i, j) + &
               & stress_axial(i, j)) / 3.0

          ! new creep strain
          creep_strain_new(i, j)%rad = creep_strain_old(i, j)%rad + &
            & (stress_rad(i, j) - stress_p) * creep_strain_incr_factor

          creep_strain_new(i, j)%hoop = creep_strain_old(i, j)%hoop + &
            & (stress_hoop(i, j) - stress_p) * creep_strain_incr_factor

          creep_strain_new(i, j)%ax = creep_strain_old(i, j)%ax + &
            & (stress_axial(i, j) - stress_p) * creep_strain_incr_factor

        else

          ! Von Mises was set to zero for first iteration
          creep_strain_new(i, j)%rad = creep_strain_old(i, j)%rad
          creep_strain_new(i, j)%hoop = creep_strain_old(i, j)%hoop
          creep_strain_new(i, j)%ax = creep_strain_old(i, j)%ax

        endif

    end do

    ! get creep sum; needed for integration
    creep_sum(i, :) = creep_strain_new(i, :)%rad + &
       & creep_strain_new(i, :)%hoop + &
       & 2.0 * poisson_ratio(i) * creep_strain_new(i, :)%ax

    ! get creep difference
    creep_dif(i, :) = creep_strain_new(i, :)%rad - &
       & creep_strain_new(i, :)%hoop

  end do

  ! I_i = int(low limit r_i+1, upper limit r) (F * x) dx 
  !      since r_i+1 < r
  
  ! indices are backwards, i.e., lower indices point to larger radius

  do i = i_cr_start, i_cr_end

    creep_integral_sum(i, npr_st(i)) = 0.0
    creep_integral_dif(i, npr_st(i)) = 0.0

    do j = npr_st(i) - 1, 1, -1

      ! temperature point 
      j_temp = 2 * j

      ! fractional integral; F * r * dr
      ! the trapezium integral should be (a+b)*h/2
      Delta_int = 0.5 * (creep_sum(i, j) * rad_st(i, j) + &
                & creep_sum(i, j+1) * rad_st(i, j+1)) * &
                & (rad_st(i, j) - rad_st(i, j+1))

      ! integral F*x dx
      creep_integral_sum(i, j) = creep_integral_sum(i, j + 1) + Delta_int

      ! fractional integral; (F / r) * dr
      ! the trapezium integral should be (a+b)*h/2
      Delta_int = 0.5 * (creep_dif(i, j) / rad_st(i, j) + &
                & creep_dif(i, j+1) / rad_st(i, j+1)) * &
                & (rad_st(i, j) - rad_st(i, j+1))

      ! integral (F / r) dr
      creep_integral_dif(i, j) = creep_integral_dif(i, j + 1) + Delta_int

      if (debug_local)  then

        ! write(aux_lun, 33) i, j, j_temp, rad_st(i, j), rad_st(i, j+1), &
        ! & rad_st(i, j) - rad_st(i, j+1), rad_temp(i, j_temp), &
        ! & tau_th(i, j_temp), DI_tau_r_dr
 33     format('dintegral ', i1, 2(1x, i2), 1x, 10(1pe13.6, 1x))

      endif

    end do

  end do

  ! obtain the creep strain integral in the z direction
  do i = i_cr_start, i_cr_end

    creep_integral_z(i) = 0.0

    do j = npr_st(i) - 1, 1, -1

      ! fractional integral; F * r * dr
      ! the trapezium integral should be (a+b)*h/2
      Delta_int = 0.5 * (creep_strain_new(i, j)%ax * rad_st(i, j) + &
           & creep_strain_new(i, j+1)%ax * rad_st(i, j+1)) * &
           & (rad_st(i, j) - rad_st(i, j+1))

      ! integral F*x dx
      creep_integral_z(i) = creep_integral_z(i) + Delta_int

    end do

  end do

  RETURN

  END SUBROUTINE CREEP_STRAIN_MEAN

  SUBROUTINE OUT_VON_MISES(N, id_high_or_low, cycle_no, time_now, dt)

    !=======================================================================
    ! Purpose(s):
    !
    ! print out von Mises stress
    !=======================================================================
 
  use parameter_module, only: moxide, mave, id_low, id_high, mgrid
  use creep_data_module, only: stress_eq, stress_eq_ave, &
      & creep_strain_eq_ave, creep_strain_incr_eq, creep_strain_ave
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value, Youngs_modul, nyoungs, youngs_temp, &
      & poisson_ratio, first_oxide_layer, no_metal_layers
  use solution_data_module, only: npr_temp, npr_st, Temp, rad_int, &
      & stress_axial, stress_hoop, stress_rad, eps0, &
      & creep_eq_hf_out, vm_hf_out, creep_eq_hf_out_ave, vm_hf_out_ave, &
      & oxide_thickness_hf_out, time_hf_out, creep_rate_out, &
      & id_full_load_hf, n_hf_out, creep_hf_time, creep_hf_dox, &
      & creep_strain_ave_hoop_out
  use property_module, only: OPERAND_ARRAY, LINEAR_PROPERTY_SIMPLE
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun, &
      & gen_lun, prefix
  use boiler_data_module, only: oxide_thickness

  ! Arguments
  integer, intent(IN)  :: N, id_high_or_low, cycle_no
  real, intent(IN)     :: time_now, dt

  ! local variables
  character(LEN = 80), dimension(mave), save    :: operation_type
  character(LEN = 26)  ::  ch_out, ch_out_ratio, ch_out_eps
  integer  :: i, k, j, loc1, ihf
  real     :: shear_module1, Young_mod, Temperature
  real, dimension(moxide, mave) :: creep_rate_ave
  real, dimension(mave)         :: creep_rate_ave_metal, &
     & stress_eq_ave_metal, creep_strain_eq_ave_metal, &
     & creep_strain_hoop_ave_metal

  ! initialize the operation type
  operation_type(1) = 'max'
  operation_type(2) = 'min'
  operation_type(3) = 'ave_sum'

  do i = 1, N

    do j = 1, npr_st(i)

    ! j = npr_st(i)  !  location where stress/strains are maximum in metal
    ! get the equivalent Von Mises stress
    stress_eq(i, j) = VON_MISES(stress_rad(i, j), stress_hoop(i, j), &
               & stress_axial(i, j))

    end do

  end do

  do i = 1, N

    do k = 1, 3

      ! get maximum, minimum, and average in each layer
      call OPERAND_ARRAY(1, mgrid, stress_eq(i, :),  &
              & 1, npr_st(i), operation_type(k), &
              & stress_eq_ave(i, k), loc1)
      ! creep rate
      call OPERAND_ARRAY(1, mgrid, creep_strain_incr_eq(i, :)/dt,  &
              & 1, npr_st(i), operation_type(k), &
              & creep_rate_ave(i, k), loc1)

    end do

  end do

  do k = 1, 3

    ! get maximum, minimum, and average in each layer
    call OPERAND_ARRAY(1, moxide, stress_eq_ave(:, k),  &
         & 1, no_metal_layers, operation_type(k), &
         & stress_eq_ave_metal(k), loc1)
    ! creep rate
    call OPERAND_ARRAY(1, moxide, creep_rate_ave(:, k),  &
         & 1, no_metal_layers, operation_type(k), &
         & creep_rate_ave_metal(k), loc1)
    ! creep strain in the metal
    call OPERAND_ARRAY(1, moxide, creep_strain_eq_ave(:, k),  &
         & 1, no_metal_layers, operation_type(k), &
         & creep_strain_eq_ave_metal(k), loc1)
    ! hoop creep strain in the metal
    call OPERAND_ARRAY(1, moxide, creep_strain_ave(:, k)%hoop,  &
         & 1, no_metal_layers, operation_type(k), &
         & creep_strain_hoop_ave_metal(k), loc1)

  end do

  if (id_high_or_low == id_high)  then
       
     ch_out = 'von_Mises_min_ave_max_hig '
     ch_out_ratio = 'vm_shearm_min_ave_max_hig '
     ch_out_eps = 'eps0_hig '

  else if (id_high_or_low == id_low)  then

     ch_out = 'von_Mises_min_ave_max_low '
     ch_out_ratio = 'vm_shearm_min_ave_max_low '
     ch_out_eps = 'eps0_low '

  endif
       
  ! rad_int is in meters
  write(2, 7) ch_out, time_now, &
         & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
         & stress_eq_ave_metal(2), stress_eq_ave_metal(3), &
         & stress_eq_ave_metal(1), &
         & (stress_eq_ave(i, 2), stress_eq_ave(i, 3), &
         & stress_eq_ave(i, 1), i = first_oxide_layer, N)
 7    format(a26, 1x, 2(1pe13.6, 1x), &
         & 50(1pe13.6, 1x))

  write(2, 7) ch_out_eps, time_now, &
         & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
         & eps0*1.0e+3

  ! write the ratio between the vm_stress and bulk modulus for creep mechanism
  ! get the Youngs modulus for layer 1
  Temperature = 0.5 * (Temp(1, 1) + &
      & Temp(no_metal_layers, npr_temp(no_metal_layers)))
  call LINEAR_PROPERTY_SIMPLE (Temperature, nyoungs(1) - 1, 0, moxide, &
        & youngs_modul(1, :), youngs_temp(1, :), Young_mod)
  shear_module1 = Young_mod / (2.0 * (1.0 + poisson_ratio(1)))

  write(2, 8) ch_out_ratio, time_now, &
         & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
         & stress_eq_ave_metal(2) / shear_module1, &
         & stress_eq_ave_metal(3) / shear_module1, &
         & stress_eq_ave_metal(1) / shear_module1
 8    format(a26, 1x, 2(1pe13.6, 1x), &
         & 50(1pe13.6, 1x))

  if (dt > 0.0)  then

    ! dt = 0 when the oxide thickness is computed by iteration initially
  do ihf = 1, n_hf_out

    if (ihf == 1 .and. cycle_no == id_full_load_hf(ihf))  then
      ! write the legend out
      write(out_lun, 14) 'creep_evolution time dox hf_out_f temp_out_f hf_in_f temp_in_f ', &
          & 'time dox hf_out_p temp_out_p hf_in_p temp_in_p '
 14    format(2(a))

    endif

    ! end of full load or beginning of partial load
    if (cycle_no == id_full_load_hf(ihf) .or. &
      & cycle_no == id_full_load_hf(ihf) + 1)  then

      creep_hf_time(id_high_or_low, ihf) = time_now
      creep_hf_dox(id_high_or_low, ihf) = &
        & oxide_thickness(cycle_no)

      ! maximum equivalent creep strain
      ! works only for no_metal_layers = 1
    ! creep_eq_hf_out(id_high_or_low, ihf, 1:N) = creep_strain_eq_ave(1:N, 1)
      creep_eq_hf_out(id_high_or_low, ihf, 1:N) = creep_strain_eq_ave(1:N, 1)
      ! maximum Von Mises
      vm_hf_out(id_high_or_low, ihf, 1) = stress_eq_ave_metal(1)
      vm_hf_out(id_high_or_low, ihf, 2:N-first_oxide_layer+2) = &
          & stress_eq_ave(first_oxide_layer:N, 1)
 
      ! average equivalent creep strain
      ! creep_eq_hf_out_ave(id_high_or_low, ihf, 1:N) = &
      !     & creep_strain_eq_ave(1:N, 3)
      creep_eq_hf_out_ave(id_high_or_low, ihf, 1) = &
          & creep_strain_eq_ave_metal(3)
      creep_eq_hf_out_ave(id_high_or_low, ihf, 2:N-first_oxide_layer+2) = &
          & creep_strain_eq_ave(first_oxide_layer:N, 3)
      
      creep_strain_ave_hoop_out(id_high_or_low, ihf, 1) = &
    &     creep_strain_hoop_ave_metal(3)
      creep_strain_ave_hoop_out(id_high_or_low, ihf, &
    & 2:N-first_oxide_layer+2) = creep_strain_ave(first_oxide_layer:N, 3)%hoop

      ! average Von Mises
      ! vm_hf_out_ave(id_high_or_low, ihf, 1:N) = stress_eq_ave(1:N, 3) 
      vm_hf_out_ave(id_high_or_low, ihf, 1) = stress_eq_ave_metal(3)
      vm_hf_out_ave(id_high_or_low, ihf, 2:N-first_oxide_layer+2) = &
          & stress_eq_ave(first_oxide_layer:N, 3)

      ! store the creep rate, creep rate gradient in the metal [1/h]
      creep_rate_out(id_high_or_low, ihf, 1:3) = creep_rate_ave_metal(1:3)
      do i = first_oxide_layer, N
        creep_rate_out(id_high_or_low, ihf, 3*(i-1)+1:3*i) = &
           & creep_rate_ave(i, 1:3)
      end do

    endif

    ! write after both the full and low load data were obtained
    if (id_high_or_low == id_low .and. &
      & cycle_no == id_full_load_hf(ihf) + 1)  then

      write(out_lun, 17)  creep_hf_time(id_high, ihf), &
        & creep_hf_dox(id_high, ihf), &
        & creep_eq_hf_out(id_high, ihf, 1:N-first_oxide_layer+2), &
        & creep_strain_ave_hoop_out(id_high, ihf, 1:N-first_oxide_layer+2), &
        & vm_hf_out(id_high, ihf, 1:N-first_oxide_layer+2)
 17   format('creep_evolution ',50(1pe13.6, 1x))

      if (ihf == n_hf_out)  then

        ! write only beginning and end
        if (N == no_metal_layers + 1)  then   ! one oxide layer
          write(tty_lun, 14) 'problem_creep time dox cr_m cr_ox ', &
          & ' von_Mis_m von_Mis_ox'
        else if (N == no_metal_layers + 2)  then   ! two oxide layers
          write(tty_lun, 14) 'problem_creep time dox cr_m cr_sp cr_mag ', &
          & ' von_Mis_m von_Mis_sp von_Mis_mag'
        else if (N == no_metal_layers + 3)  then   ! three oxide layers
          write(tty_lun, 14) 'problem_creep time dox cr_m cr_sp cr_mag ', &
          & ' von_Mis_m von_Mis_sp von_Mis_mag'
        endif

        do j = 1, n_hf_out
          if (j == 1 .or. j == n_hf_out)  then
            write(tty_lun, 18) TRIM(prefix), creep_hf_time(id_high, j), &
              & creep_hf_dox(id_high, j), &
              & creep_eq_hf_out(id_high, j, 1:N-first_oxide_layer+2), &
              & vm_hf_out(id_high, j, 1:N-first_oxide_layer+2)
          endif
        end do

      endif

      ! write only at the end of run
      if (gen_lun > 0 .and. ihf == n_hf_out)  then

        write(gen_lun, *)

        ! write only beginning and end
        if (N == no_metal_layers + 1)  then   ! one oxide layer
          ! ah indicate average hoop
          write(gen_lun, 14) 'problem_creep time dox cr_m cr_ox ', &
          & ' cr_ah_m cr_ah_ox von_Mis_m von_Mis_ox'
        else if (N == no_metal_layers + 2)  then   ! two oxide layers
          write(gen_lun, 14) 'problem_creep time dox cr_m cr_sp cr_mag ', &
          & ' cr_ah_m cr_ah_sp cr_ah_mag von_Mis_m von_Mis_sp von_Mis_mag'
        else if (N == no_metal_layers + 3)  then   ! three oxide layers
          write(gen_lun, 14) 'problem_creep time dox cr_m cr_sp cr_mag ', &
          & ' cr_ah_m cr_ah_sp cr_ah_mag von_Mis_m von_Mis_sp von_Mis_mag'
        endif

        do j = 1, n_hf_out
          if (j == 1 .or. j == n_hf_out)  then
            write(gen_lun, 18) TRIM(prefix), creep_hf_time(id_high, j), &
              & creep_hf_dox(id_high, j), &
              & creep_eq_hf_out(id_high, j, 1:N-first_oxide_layer+2), &
         & creep_strain_ave_hoop_out(id_high, j, 1:N-first_oxide_layer+2), &
              & vm_hf_out(id_high, j, 1:N-first_oxide_layer+2)
          endif
        end do

        write(gen_lun, *) ' average creep strain and von Mises stress'
        if (N == no_metal_layers + 1)  then   ! one oxide layer
          write(gen_lun, 14) 'ave_problem_creep time dox cr_m cr_ox cr_h_m', &
          & ' cr_h_ox von_Mis_m von_Mis_ox'
        else if (N == no_metal_layers + 2)  then   ! two oxide layers
          write(gen_lun, 14) 'ave_problem_creep time dox cr_m cr_sp cr_mag ', &
          & ' cr_h_m cr_h_sp cr_h_mag von_Mis_m von_Mis_sp von_Mis_mag'
        else if (N == no_metal_layers + 3)  then   ! three oxide layers
          write(gen_lun, 14) 'ave_problem_creep time dox cr_m cr_sp cr_mag ', &
          & ' cr_h_m cr_h_sp cr_h_mag von_Mis_m von_Mis_sp von_Mis_mag'
        endif
        do j = 1, n_hf_out
          if (j == 1 .or. j == n_hf_out)  then
            write(gen_lun, 18) TRIM(prefix), creep_hf_time(id_high, j), &
              & creep_hf_dox(id_high, j), &
         & creep_eq_hf_out_ave(id_high, j, 1:N-first_oxide_layer+2), &
         & creep_strain_ave_hoop_out(id_high, j, 1:N-first_oxide_layer+2), &
         & vm_hf_out_ave(id_high, j, 1:N-first_oxide_layer+2)
          endif
        end do

        write(gen_lun, *)

        ! write all the data
        if (N == no_metal_layers + 1)  then   ! one oxide layer
          write(gen_lun, 14) 'ave_problem_creep time dox cr_m cr_ox cr_h_m', &
          & ' cr_h_ox von_Mis_m von_Mis_ox'
        else if (N == no_metal_layers + 2)  then   ! two oxide layers
          write(gen_lun, 14) 'ave_problem_creep time dox cr_m cr_sp cr_mag ', &
          & ' cr_h_m cr_h_sp cr_h_mag von_Mis_m von_Mis_sp von_Mis_mag'
        else if (N == no_metal_layers + 3)  then   ! three oxide layers
          write(gen_lun, 14) 'ave_problem_creep time dox cr_m cr_sp cr_mag ', &
          & ' cr_h_m cr_h_sp cr_h_mag von_Mis_m von_Mis_sp von_Mis_mag'
        endif

        do j = 1, n_hf_out
          write(gen_lun, 18) TRIM(prefix), creep_hf_time(id_high, j), &
            & creep_hf_dox(id_high, j), &
            & creep_eq_hf_out_ave(id_high, j, 1:N-first_oxide_layer+2), &
          & creep_strain_ave_hoop_out(id_high, j, 1:N-first_oxide_layer+2), &
            & vm_hf_out_ave(id_high, j, 1:N-first_oxide_layer+2)
        end do

        ! write all the data
        if (N == no_metal_layers + 1)  then   ! one oxide layer
          write(gen_lun, 14) 'time dox cr_rt_m_max cr_rt_m_min ', &
             & 'cr_rt_m_ave cr_rt_ox_max cr_rt_ox_min cr_rt_ox_ave'
        else if (N == no_metal_layers + 2)  then   ! one oxide layer
          write(gen_lun, 14) 'time dox cr_rt_m_max cr_rt_m_min ', &
             & 'cr_rt_m_ave cr_rt_sp_max cr_rt_sp_min cr_rt_sp_ave', &
             & 'cr_rt_mag_max cr_rt_mag_min cr_rt_mag_ave'
        else if (N == no_metal_layers + 3)  then   ! one oxide layer
          write(gen_lun, 14) 'time dox cr_rt_m_max cr_rt_m_min ', &
             & 'cr_rt_m_ave cr_rt_sp_max cr_rt_sp_min cr_rt_sp_ave', &
             & 'cr_rt_mag_max cr_rt_mag_min cr_rt_mag_ave'
        endif
        do j = 1, n_hf_out
          write(gen_lun, 13)  creep_hf_time(id_high, j), &
            & creep_hf_dox(id_high, j), &
            & creep_rate_out(id_high_or_low, j, 1:3*(N-no_metal_layers+1))
        end do

 13      format(30(1pe13.6, 1x))
 18      format((a),1x, 30(1pe13.6, 1x))

      endif

    endif

  end do  

  endif

  return

  END SUBROUTINE OUT_VON_MISES

  REAL FUNCTION CREEP_NEW_OXIDE(option, creep_old1, creep_old2, &
        & r1_old, r2_old, r1_new, r2_new)

  ! 2 is at the interface, 1 is inside

  real, Intent(IN)  :: creep_old1, creep_old2, &
        & r1_old, r2_old, r1_new, r2_new
  character(LEN = 80), intent(IN) :: option

  !local vars
  real  :: fnew, new_ox_term, inside_term

  SELECT CASE (TRIM(option))

    CASE ('old')

      CREEP_NEW_OXIDE = creep_old2

    CASE ('zero')

      CREEP_NEW_OXIDE = 0.0

    CASE ('half_old')

      CREEP_NEW_OXIDE = creep_old2 / 2.0

    CASE ('interpolation_zero')

      fnew = 0.0
      new_ox_term = 2.0 * fnew * (r2_new - r2_old)
      ! p1-o2 term; creep_old1 - should be replaced with p1 interpolated value
      inside_term = (creep_old1 + creep_old2) * (r2_old - r1_new)
      CREEP_NEW_OXIDE = -creep_old1 + &
        & (new_ox_term + inside_term) / (r2_new - r1_new)

    CASE ('interpolation_half')

      fnew = creep_old2 / 2.0
      new_ox_term = 2.0 * fnew * (r2_new - r2_old)
      ! p1-o2 term; creep_old1 - should be replaced with p1 interpolated value
      inside_term = (creep_old1 + creep_old2) * (r2_old - r1_new)
      CREEP_NEW_OXIDE = -creep_old1 + &
        & (new_ox_term + inside_term) / (r2_new - r1_new)

    CASE DEFAULT

      CREEP_NEW_OXIDE = creep_old2
      write(6, *) 'new_oxide_creep_option must be: ', &
       & 'old, zero, half_old, interpolation_zero, or ', &
       & 'interpolation_half '
      stop

  END SELECT
    
  RETURN

  END FUNCTION CREEP_NEW_OXIDE

  REAL FUNCTION CREEP_TERM_INT_AXIAL(Poisson, Youngs_module, &
          & Integral_creep_etr_sum, Integral_creep_eor_dif, &
          & Integral_creep_etr_z, r2)

  real, Intent(IN)  :: Poisson, Youngs_module, Integral_creep_etr_z, &
          & Integral_creep_etr_sum, Integral_creep_eor_dif, r2

  !local vars
  real  :: A4, A5

  ! RHS term for CREEP; for the generalized plain strain case
  ! A5= 2*nu*E1*(1+nu*) = 2 * nu * E^* / (1-nu*); we use A5/2
  A5 = Poisson / ((1+Poisson) * (1.0 - 2.0*Poisson)) 

  CREEP_TERM_INT_AXIAL = Youngs_module * A5 * &
    & (Integral_creep_etr_sum - Integral_creep_eor_dif * r2)

  ! we should use A4 and add A6 term 4/22/2009
  ! A4 * E^* / 2 = 0.5 * E * nu / (1 - nu**2)
  A4 = 0.5 * Poisson / (1.0 - Poisson**2)

  CREEP_TERM_INT_AXIAL = Youngs_module * (A4 * &
    & (Integral_creep_etr_sum - Integral_creep_eor_dif * r2) + &
    & Integral_creep_etr_z)

  RETURN

  END FUNCTION CREEP_TERM_INT_AXIAL

  REAL FUNCTION CREEP_TERM_S_RAD(Integral_creep_etr_sum, &
          & Integral_creep_eor_dif, r2)

  real, Intent(IN)  :: Integral_creep_etr_sum, &
          & Integral_creep_eor_dif, r2

  ! RHS term for CREEP on the inner side of the interface (LHS of interface)
  CREEP_TERM_S_RAD = 0.5 * (Integral_creep_etr_sum / r2 - &
          & Integral_creep_eor_dif)

  RETURN

  END FUNCTION CREEP_TERM_S_RAD

  REAL FUNCTION CREEP_TERM_S_HOOP(Integral_creep_etr_sum, &
          & Integral_creep_eor_dif, r2)

  real, Intent(IN)  :: Integral_creep_etr_sum, &
          & Integral_creep_eor_dif, r2

  ! RHS term for CREEP on the inner side of the interface (LHS of interface)
  CREEP_TERM_S_HOOP = 0.5 * (Integral_creep_etr_sum / r2 + &
          & Integral_creep_eor_dif)

  RETURN

  END FUNCTION CREEP_TERM_S_HOOP

  REAL FUNCTION CREEP_TERM_DUDR(Integral_creep_etr_sum, &
          & Integral_creep_eor_dif, poisson_star, r2)

  real, Intent(IN)  :: Integral_creep_etr_sum, &
          & Integral_creep_eor_dif, poisson_star, r2

  ! RHS term for CREEP on the inner side of the interface (LHS of interface)
  CREEP_TERM_DUDR = 0.5 * (-(1.0 + poisson_star) * &
          & Integral_creep_etr_sum / r2 + &
          & (1.0 - poisson_star) * Integral_creep_eor_dif)

  RETURN

  END FUNCTION CREEP_TERM_DUDR

  REAL FUNCTION CREEP_TERM_U(Integral_creep_etr_sum, &
          & Integral_creep_eor_dif, poisson_star, r)

  real, Intent(IN)  :: Integral_creep_etr_sum, &
          & Integral_creep_eor_dif, poisson_star, r

  ! RHS term for CREEP on the inner side of the interface (LHS of interface)
  CREEP_TERM_U = -0.5 * ((1.0 + poisson_star) * &
          & Integral_creep_etr_sum / r + &
          & (1.0 - poisson_star) * Integral_creep_eor_dif * r)

  RETURN

  END FUNCTION CREEP_TERM_U

  REAL FUNCTION DCREEP(sigma_eq, dt, Temp, ilay)
  
    use creep_data_module, only: creep_const, creep_exponent, &
       & creep_activation_energy, creep_activ_energy_stress, &
       & Temp_creep_onset, creep_const_ln
    use creep_data_module, only: if_rupture_time_creep_mat, &
       & rupture_time_c1, rupture_time_t1, rupture_time_exp, &
       & rupture_time_t1s1, rupture_time_t1s2, &
       & if_rupture_time_creep
    ! data for the mixture
    use oxide_data_module, only: no_layer_mat, layer_mat_id, fraction_material
    use creep_data_module, only: creep_const_mat_ln, creep_const_mat, creep_exponent_mat, &
       & creep_activation_energy_mat, creep_activ_energy_stress_mat

  ! equivalent creep strain
  real,   Intent(IN)  :: sigma_eq, dt, Temp
  integer, Intent(IN) :: ilay

  real  :: sigma_expon, exponent, log10_sigma, sigma_eq_log
  integer :: i, k
  logical, save :: if_first_ln1 = .true.

if (no_layer_mat(ilay) == 1)  then

    DCREEP = DCREEP_ONE(sigma_eq, dt, Temp, ilay)

else if (no_layer_mat(ilay) > 1)  then

  if (if_first_ln1)  then

     do i = 1, no_layer_mat(ilay)

       k = layer_mat_id(ilay, i)
       creep_const_mat_ln(k) = LOG(creep_const_mat(k))

        write(2, 2) ilay, i, k, fraction_material(ilay, i), &
             & creep_exponent_mat(k), creep_activation_energy_mat(k), &
             & exp(creep_const_mat_ln(k))
 2  format('dcreep_check3 ', 3(i1, 1x), 20(1pe13.6, 1x))

     end do
  
     if_first_ln1 = .false.

  endif

  if (Temp > Temp_creep_onset(ilay))  then

    ! DCREEP = creep_const(ilay) * (sigma_eq**creep_exponent(ilay)) * &
    !       & EXP(-creep_activation_energy(ilay) / &
    !       & (Temp + 273.0)) * dt

    if (sigma_eq > 0.0)  then

      DCREEP = 0.0
      sigma_eq_log = LOG(sigma_eq)

      do i = 1, no_layer_mat(ilay)

        k = layer_mat_id(ilay, i)

        if (k == 0)  then

          write(6, *) 'ERROR: layer_mat_id not defined in DCREEP'
          STOP

        endif

        if (creep_exponent_mat(k) > 0.0)  then
          sigma_expon = creep_exponent_mat(k) * sigma_eq_log
        else
          sigma_expon = 0.0
        endif

        exponent = creep_const_mat_ln(k) + sigma_expon - &
          & (creep_activation_energy_mat(k) + &
          &  creep_activ_energy_stress_mat(k) * sigma_eq) / (Temp + 273.0)

        DCREEP = DCREEP + fraction_material(ilay, i) * EXP(exponent) * dt

      end do

      if (if_rupture_time_creep(ilay))  then

        write(6, *) 'ERROR mixture CREEP involving rupture_time is only for metal'
        stop
 
      endif

      ! write(2, 1) ilay, sigma_eq, creep_exponent(ilay), &
      ! & creep_const_ln(ilay), creep_activation_energy(ilay), &
      ! & Temp, dt, sigma_expon, exponent, DCREEP
 1  format('dcreep_check2 ', i1, 1x, 10(1pe13.6, 1x))

    else

      ! if creep would be uni-axial, then the sign of the stress would have 
      ! to be considered and we should have
      ! sign(sigma)*sigma**creep_exponent(ilay)
      DCREEP = 0.0

    endif

  else

    DCREEP = 0.0

  endif

endif

  return

  END FUNCTION DCREEP

  REAL FUNCTION DCREEP_ONE(sigma_eq, dt, Temp, ilay)
  
    use creep_data_module, only: creep_const, creep_exponent, &
       & creep_activation_energy, creep_activ_energy_stress, &
       & Temp_creep_onset, creep_const_ln
    use creep_data_module, only: if_rupture_time_creep_mat, &
       & rupture_time_c1, rupture_time_t1, rupture_time_exp, &
       & rupture_time_t1s1, rupture_time_t1s2, &
       & if_rupture_time_creep

  ! equivalent creep strain
  real,   Intent(IN)  :: sigma_eq, dt, Temp
  integer, Intent(IN) :: ilay

  real  :: sigma_expon, exponent, log10_sigma

  if (Temp > Temp_creep_onset(ilay))  then

    ! DCREEP_ONE = creep_const(ilay) * (sigma_eq**creep_exponent(ilay)) * &
    !       & EXP(-creep_activation_energy(ilay) / &
    !       & (Temp + 273.0)) * dt

    if (sigma_eq > 0.0)  then

      if (creep_exponent(ilay) > 0.0)  then
        sigma_expon = creep_exponent(ilay) * LOG(sigma_eq)
      else
        sigma_expon = 0.0
      endif

      exponent = creep_const_ln(ilay) + sigma_expon - &
        & (creep_activation_energy(ilay) + &
        &  creep_activ_energy_stress(ilay) * sigma_eq) / (Temp + 273.0)

      DCREEP_ONE = EXP(exponent) * dt

      if (if_rupture_time_creep(ilay))  then

        log10_sigma = LOG10(sigma_eq)

        ! exponent is used as dummy for log10(rupture time)
        exponent = rupture_time_c1(ilay) + (rupture_time_t1(ilay) + &
          & log10_sigma * (rupture_time_t1s1(ilay) + &
          & rupture_time_t1s2(ilay) * log10_sigma)) / &
          & (Temp + 273.0)
        DCREEP_ONE = DCREEP_ONE * 10.0**(rupture_time_exp(ilay) * exponent)

      endif

      ! write(2, 1) ilay, sigma_eq, creep_exponent(ilay), &
      ! & creep_const_ln(ilay), creep_activation_energy(ilay), &
      ! & Temp, dt, sigma_expon, exponent, DCREEP
 1  format('dcreep_check1 ', i1, 1x, 10(1pe13.6, 1x))

    else

      ! if creep would be uni-axial, then the sign of the stress would have 
      ! to be considered and we should have
      ! sign(sigma)*sigma**creep_exponent(ilay)
      DCREEP_ONE = 0.0

    endif

  else

    DCREEP_ONE = 0.0

  endif

  return

  END FUNCTION DCREEP_ONE

  REAL FUNCTION DCREEP_DVM(sigma_eq, Temp, DCREEP, ilay)
  
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the d(e_cr_increm)/d(sigma_eq)
    ! 
    !=======================================================================
    use creep_data_module, only: creep_const, creep_exponent, &
       & creep_activation_energy, Temp_creep_onset, creep_const_ln
    use creep_data_module, only: if_rupture_time_creep_mat, &
       & rupture_time_c1, rupture_time_t1, rupture_time_exp, &
       & rupture_time_t1s1, rupture_time_t1s2, &
       & if_rupture_time_creep

  ! equivalent creep strain
  real,   Intent(IN)  :: sigma_eq, Temp, DCREEP
  integer, Intent(IN) :: ilay

  real  :: sigma_expon, exponent, log10_sigma

  if (Temp > Temp_creep_onset(ilay))  then

    ! DCREEP = creep_const(ilay) * (sigma_eq**creep_exponent(ilay)) * &
    !       & EXP(-creep_activation_energy(ilay) / &
    !       & (Temp + 273.0)) * dt

    if (sigma_eq > 0.0)  then

      if (ABS(creep_exponent(ilay)) > 1.0e-9)  then 
        DCREEP_DVM = DCREEP * creep_exponent(ilay) / sigma_eq
      endif

      if (ABS(creep_exponent(ilay)) <= 1.0e-9 .and. &
        & (.not. if_rupture_time_creep(ilay)))  then
        DCREEP_DVM = 0.0
      endif

      if (if_rupture_time_creep(ilay))  then

        log10_sigma = LOG10(sigma_eq)

        DCREEP_DVM = DCREEP * rupture_time_exp(ilay) * &
          & (rupture_time_t1s1(ilay) + 2.0 * rupture_time_t1s2(ilay) * &
          & log10_sigma) / (sigma_eq * (Temp + 273.0))

      endif

    else

      DCREEP_DVM = 0.0

    endif

  else

    DCREEP_DVM = 0.0

  endif

  RETURN

  END FUNCTION DCREEP_DVM

  REAL FUNCTION VON_MISES(stress_r, stress_theta, stress_z)
  
  ! Von Mises sequivalent stress
  real,   Intent(IN)  :: stress_r, stress_theta, stress_z

  real :: dummy

  dummy = (stress_r - stress_theta)**2 + &
     & (stress_r - stress_z)**2 + (stress_z - stress_theta)**2

  if (dummy > 0.0)  then

    VON_MISES = SQRT(0.5 * dummy) 

  else

    VON_MISES = 0.0

  endif

  return

  END FUNCTION VON_MISES
  

END MODULE CREEP_MODULE
