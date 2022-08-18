MODULE THERMAL_CYCLE_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for temperature profile. 
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
  public :: TEMP_SS_HOLLOW_CYL, TEMP_PROFILE, BETA_TAU_TH

 CONTAINS

  SUBROUTINE TEMP_SS_HOLLOW_CYL (Temp_subs_out, Temp_coat, &
       & r_subs_coat, r_subs_out, Oxide_Thickness, cond1, cond2, r_now, Temp)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the critical strain for each failure mode
    !  Level II of EPRI diagram using:
    !  eta - material properties and geometry effects
    !  omega - damage kinetics parameter dependent on the operation history
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, pi
    use output_module,      only: tty_lun, out_lun, aux_lun
    use solver_data_module, only: no_thick_interval, no_temp_strain
    use property_module,    only: LINEAR_PROPERTY_SIMPLE

    ! Argument List
    real,    intent(IN)  :: Temp_subs_out, Temp_coat
    real,    intent(IN)  :: r_now, r_subs_coat, r_subs_out, Oxide_Thickness
    real,    intent(IN)  :: cond1, cond2
    real,    intent(OUT) :: Temp

    ! Local Variables

    real   :: den1, nom1, deltaT, r_coat

    ! 
    ! coating is on the inside
    r_coat = r_subs_coat - oxide_thickness

    den1 = cond1 * LOG(r_subs_out / r_subs_coat) + cond2 * LOG(r_subs_coat / r_coat)
    deltaT = Temp_subs_out - Temp_coat

    if (r_now >= r_coat .and. r_now < r_subs_coat)  then
      ! coating
      nom1 = cond2 * LOG(r_now / r_coat)

    else if (r_now >= r_subs_coat .and. r_now <= r_subs_out)  then
      ! substrate
     
      nom1 = cond2 * LOG(r_subs_coat / r_coat) + cond1 * LOG(r_now / r_subs_coat)

    else

      write(out_lun, *) 'ERROR: radius not available'
      stop

    endif

    Temp = Temp_coat + deltaT * nom1 / den1

   return

  END SUBROUTINE TEMP_SS_HOLLOW_CYL

  SUBROUTINE TEMP_PROFILE(N, Temp_out, Temp_in, thickness_oxide_layer, &
     & cycle_no, time_now, dt, id_high_or_low, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the temperature profile
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module, only: moxide, mave, id_low, id_high, mgrid
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer, oxide_thickness
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int, heat_flux_temp_out, oxide_thickness_hf_out, time_hf_out, &
      & id_full_load_hf, n_hf_out, Temp_mean
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value, first_oxide_layer, no_metal_layers
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, if_cylindrical, no_st_rad_points_tube, &
      & no_st_rad_points_oxide, &
      & if_different_ref_temp, Temp_reference_tube
  use property_module, only: OPERAND_ARRAY
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun, &
      & gen_lun, prefix

  integer, intent(IN)  :: N, id_high_or_low, cycle_no
  real,    intent(IN)  :: Temp_out, Temp_in, time_now, dt
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale

  ! Local Variables
  integer  :: i, j, loc1, k, ihf
  real, dimension(N+1) :: r
  real, dimension(N)   :: cond
  real, dimension(N-1) :: htc
  real  :: htc_in, htc_out, heat_flux_out, heat_flux_in, Temp_metal_ave
  real, dimension(N)   :: ac, bc
  character(LEN = 80), dimension(mave), save    :: operation_type
  character(LEN = 13) :: ch_temp1
  logical, save   :: if_first_temp_rad = .true.
  logical    :: debug_local = .false.

  ! N = total number of layers; minimum 2 (substrate + average scale)

  if (id_high_or_low > 3 .or. id_high_or_low < 0)  return

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  ! initialize
  r(1:N+1) = rad_int(1:N+1)

  ! initialize the operation type
  operation_type(1) = 'max'
  operation_type(2) = 'min'
  operation_type(3) = 'ave_sum'

   ! if_first_temp_rad = .false.

  ! endif

  ! write(aux_lun, *) 'ri ', (r(i), i = 1, N+1)

  htc_in = htc_tube_inner
  htc_out = htc_tube_outer
  htc = 0.0

  ! thermal conductivity
  ! ignore the temperature dependence for now
  cond(1) = cond_value(1, 1)
  do i = 2, N

    if (if_ave_scale)  then
      ! in this case N = 2 (tube and one average oxide)
      cond(i) = cond_value(no_oxide_ave, 1)   ! only for average
    else
      ! this will work in general
      cond(i) = cond_value(i, 1)
    endif

  end do

  call TEMP_HTC_CYL(N, Temp_out, Temp_in, r, &
         & htc, htc_in, htc_out, cond, ac, bc)

  ! store the ac and bc temperature constants
  ac_temp(id_high_or_low, 1:N) = ac(1:N)
  bc_temp(id_high_or_low, 1:N) = bc(1:N)

  ! write(aux_lun, *) 'ac ', (ac(i), i= 1, N)
  ! write(aux_lun, *) 'bc ', (bc(i), i= 1, N)

  do i = 1, N

    ! get temperatures in the oxide
    do j = 1, npr_temp(i)
      Temp(i, j) = TEMP_VAL(if_cylindrical, rad_temp(i, j), &
                 & Temp_in, ac(i), bc(i))
    end do

    do k = 1, 3
      ! get maximum temperature
      ! get minimum temperature
      ! get average temperature
      call OPERAND_ARRAY(1, mgrid, Temp(i, :), 1, npr_temp(i), &
       & operation_type(k), Temp_ave(i, k), loc1)
    end do

    Temp_mean(i) = Temp_ave(i, 3)

    if (debug_local)  then
       write(aux_lun, 2) i, Temp_ave(i, 1), Temp_ave(i, 2), Temp_ave(i, 3)
 2   format('temp_max_min_av ', i2, 1x, 10(1pe13.6, 1x))

     endif

  end do

  if (id_high_or_low == id_high)  then
    ch_temp1 = 'temp_ox_hig '
  else if (id_high_or_low == id_low)  then
    ch_temp1 = 'temp_ox_low '       
  endif
  write(out_lun, 9) ch_temp1, time_now, &
    & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
    & (Temp(i, 1), i = first_oxide_layer, N)
 9 format(a13, 1x, 50(1pe13.6, 1x))

  heat_flux_out = htc_out * (Temp_out - Temp(1,1))
  heat_flux_in = htc_in * (Temp(N, npr_temp(N)) - Temp_in)

  if (dt > 0.0)  then
    ! dt = 0 when the oxide thickness is computed by iteration initially
  do ihf = 1, n_hf_out

    if (ihf == 1 .and. cycle_no == id_full_load_hf(ihf))  then
      ! write the legend out
      write(out_lun, 4) 'heat_flux_evolution time dox hf_out_f ', &
          & 'temp_out_f temp_mean_f temp_m_ox_f temp_in_f hf_in_f ', &
          & 'time dox hf_out_p temp_out_p temp_mean_p temp_m_ox_p ', &
          & 'temp_in_p hf_in_p '
 4    format(4(a))

    endif

    ! end of full load or beginning of partial load
    if (cycle_no == id_full_load_hf(ihf) .or. &
      & cycle_no == id_full_load_hf(ihf) + 1)  then

      ! get the average temperature in the metal
      Temp_metal_ave = 0.0
      do i = 1, no_metal_layers
        Temp_metal_ave = Temp_metal_ave + Temp_ave(i, 3)
      end do
      Temp_metal_ave = Temp_metal_ave / no_metal_layers

      heat_flux_temp_out(id_high_or_low, ihf, 1) = time_now
      heat_flux_temp_out(id_high_or_low, ihf, 2) = &
        & oxide_thickness(cycle_no)
      heat_flux_temp_out(id_high_or_low, ihf, 3) = heat_flux_out
      heat_flux_temp_out(id_high_or_low, ihf, 4) = Temp(1, 1)
      heat_flux_temp_out(id_high_or_low, ihf, 5) = Temp_metal_ave
      heat_flux_temp_out(id_high_or_low, ihf, 6) = &
             & Temp(no_metal_layers, npr_temp(no_metal_layers))
      heat_flux_temp_out(id_high_or_low, ihf, 7) = &
        & Temp(N, npr_temp(N))
      heat_flux_temp_out(id_high_or_low, ihf, 8) = heat_flux_in

    endif

    ! write after both the full and low load data were obtained
    if (id_high_or_low == id_low .and. &
      & cycle_no == id_full_load_hf(ihf) + 1)  then

      write(out_lun, 7)  (heat_flux_temp_out(id_high, ihf, k), &
          & k = 1, 8), (heat_flux_temp_out(id_low, ihf, k), &
          & k = 1, 8)
 7      format('heat_flux_evolution ',30(1pe13.6, 1x))

      if (ihf == n_hf_out)  then

        ! write only beginning and end
        write(tty_lun, 4) 'problem time dox hf_out_f temp_out_f ', &
          & 'hf_out_p temp_out_p temp_furnace temp_steam htc_out htc_in'
        do j = 1, n_hf_out
          if (j == 1 .or. j == n_hf_out)  then
            write(tty_lun, 8) TRIM(prefix), (heat_flux_temp_out(id_high, j, k), &
            & k = 1, 4), (heat_flux_temp_out(id_low, j, k), &
            & k = 3, 4), Temp_out, htc_out, Temp_in, htc_in
          endif
        end do

      endif

      ! write only at the end of run
      if (gen_lun > 0 .and. ihf == n_hf_out)  then

        write(gen_lun, *)

        ! write only beginning and end
        write(gen_lun, 4) 'problem time dox hf_out_f temp_out_f ', &
          & 'hf_out_p temp_out_p temp_furnace temp_steam htc_out htc_in'
        do j = 1, n_hf_out
          if (j == 1 .or. j == n_hf_out)  then
            write(gen_lun, 8) TRIM(prefix), (heat_flux_temp_out(id_high, j, k), &
            & k = 1, 4), (heat_flux_temp_out(id_low, j, k), &
            & k = 3, 4), Temp_out, htc_out, Temp_in, htc_in
          endif
        end do

        write(gen_lun, *)

        write(gen_lun, 4) 'heat_flux_evolution time dox hf_out_f ', &
          & 'temp_out_f temp_mean_f temp_m_ox_f temp_in_f hf_in_f ', &
          & 'time dox hf_out_p temp_out_p temp_mean_p temp_m_ox_p ', &
          & 'temp_in_p hf_in_p '
        do j = 1, n_hf_out
          write(gen_lun, 3)  (heat_flux_temp_out(id_high, j, k), &
            & k = 1, 8), (heat_flux_temp_out(id_low, j, k), &
            & k = 1, 8)
        end do

 3      format(30(1pe13.6, 1x))
 8      format((a),1x, 30(1pe13.6, 1x))

      endif

    endif

  end do  

  if (id_high_or_low == id_low)  then

    Temp_low_rad = Temp
    write(out_lun, 5) time_now, heat_flux_out, heat_flux_in, &
      & Temp_out - Temp(1,1), Temp(N, npr_temp(N)) - Temp_in, &
      & Temp_in, Temp(N, npr_temp(N)), Temp(1,1), Temp_out
 5  format('heat_flux_dT_low ', 10(1pe13.6, 1x))

  else if (id_high_or_low == id_high)  then

    Temp_high_rad = Temp
    write(out_lun, 6) time_now, heat_flux_out, heat_flux_in, &
      & Temp_out - Temp(1,1), Temp(N, npr_temp(N)) - Temp_in, &
      & Temp_in, Temp(N, npr_temp(N)), Temp(1,1), Temp_out
 6  format('heat_flux_dT_high ', 10(1pe13.6, 1x))

  endif

  endif

  return

  END SUBROUTINE TEMP_PROFILE

  SUBROUTINE TEMP_HTC_CYL(N, Temp_out, Temp_in, r, &
         & htc, htc_in, htc_out, cond, ac, bc)
    !=======================================================================
    ! Purpose(s):
    !
    !  set up coeff for calculating temperature constants
    !  calculate temperature constants ac and bc 
    ! 
    !=======================================================================

  ! T_i = Tin + A_i + B_i * ln(r)

  ! 
  use solver_data_module, only : if_cylindrical, internal_htc, external_htc
  use solver_module, only      : SOLVER_SIMPLE
  use output_module,      only: tty_lun, out_lun, aux_lun

  ! Argument List
  integer,    intent(IN)        :: N
  real,    intent(IN)  :: Temp_out, Temp_in
  real, intent(IN)  :: htc_in, htc_out
  real, dimension(N+1), intent(IN) :: r
  real, dimension(N), intent(IN)   :: cond
  real, dimension(N-1), intent(INOUT) :: htc
  real, dimension(N), intent(OUT)    :: ac, bc

  ! Local Variables
  integer  :: i, ieq, j
  integer, save  :: ip = 0
  real(kind = 8), dimension(2*N, 2*N) :: coeff
  real(kind = 8), dimension(2*N)      :: rhs
  real, dimension(N)   :: cond_htc  ! conduction for imposing htc internal
  real :: cond_1, cond_N, htc_in_n, htc_out_n
  logical  :: if_cyl
  integer, dimension(N) :: ia, ib

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  if_cyl = if_cylindrical

  ! unknown indices for A_i
  do i = 1, N
    ia(i) = i
  end do

  do i = 1, N
    ib(i) = ia(N) + i
  end do

  ! initialize
  coeff = 0.0
  rhs = 0.0

  ! write(aux_lun, *) 'htc1 ', external_htc, internal_htc, htc_in, htc_out

  if (.not. external_htc)  then

    htc_out_n = 1.0  ! local in this routine
    htc_in_n = 1.0
    cond_1 = 0.0
    cond_N = 0.0

  else
  
    htc_out_n = htc_out  ! local in this routine
    htc_in_n = htc_in
    cond_1 = cond(1)
    cond_N = cond(N)

  endif

  if (.not. internal_htc)  then

      cond_htc = 0.0
      htc = 1.0

  else

      cond_htc = cond

  endif
  
  ! equation 1, heat flux balance at tube outer surface
  ieq = 1
  coeff(ieq, ia(1)) = htc_out_n
  coeff(ieq, ib(1)) = htc_out_n * TP(if_cyl, r(1)) + &
     & cond_1 * DTP_DR(if_cyl, r(1))
  rhs(ieq) = htc_out_n * (Temp_out - Temp_in)

  ! equation 2, heat flux balance at the inner tube surface
  ieq = ieq + 1
  coeff(ieq, ia(N)) = htc_in_n

  coeff(ieq, ib(N)) = htc_in_n * TP(if_cyl, r(N+1)) - &
     & cond_N * DTP_DR(if_cyl, r(N+1))

  ! second type of equations
  ! heat balance at the interface, thermal gradient
  do i = 2, N
    ieq = ieq + 1
    coeff(ieq, ib(i)) = cond(i)
    coeff(ieq, ib(i-1)) = - cond(i-1)
  end do

    if (ip <= 6)  then
      write(aux_lun, 7) htc_in_n, cond_N, r(N+1), TP(if_cyl, r(N+1)), &
     & DTP_DR(if_cyl, r(N+1))
 7    format('second_eq ', 20(1pe13.6, 1x))
    endif

  ! third type of equations
  ! heat balance at the interface, thermal gradient with htc
  do i = 2, N
    ieq = ieq + 1
    coeff(ieq, ia(i)) = htc(i-1)
    coeff(ieq, ia(i-1)) = - htc(i-1)

    coeff(ieq, ib(i)) = htc(i-1) * TP(if_cyl, r(i)) + &
       & cond_htc(i) * DTP_DR(if_cyl, r(i))
    coeff(ieq, ib(i-1)) = - htc(i-1) * TP(if_cyl, r(i))

    if (ip <= 6)  then
      write(aux_lun, 3) ieq, i, r(i), htc(i-1), TP(if_cyl, r(i)), &
         & cond_htc(i), DTP_DR(if_cyl, r(i))
 3    format('bfs_b ', 2(i2, 1x), 20(1pe13.6, 1x))
    endif

  end do
  
  if (ip <= 6)  then

    do i = 1, 2*N
      write(aux_lun, 2)  i, (coeff(i, j), j = 1, 2*N), rhs(i)
 2    format('bfs ', i2, 1x, 20(1pe13.6, 1x))
    end do

  endif

  ! determine the temperature constants
  call SOLVER_SIMPLE(2*N, coeff, rhs)

  ! store the results
  do i = 1, N
    ac(i) = rhs(ia(i))
    bc(i) = rhs(ib(i))
  end do

  if (ip <= 6)  then
    write(aux_lun, 4)  (rhs(i), i = 1, 2*N)
 4    format('bfs_sol ', 20(1pe13.6, 1x))
  endif

  return

  END SUBROUTINE TEMP_HTC_CYL

  SUBROUTINE BETA_TAU_TH(N, if_ave_scale, id_high_or_low)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the increment in integral of tau
    !    as the integral of th_expansion from 
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    ! 
    !=======================================================================

    use parameter_module,   only: moxide, id_low, id_high, mgrid
    use output_module,      only: tty_lun, out_lun, aux_lun
    use oxide_data_module,  only: th_exp_temp, &
          & th_exp_coeff, nth_exp, no_oxide_ave, first_oxide_layer, &
          & no_metal_layers
    use property_module,    only: LINEAR_PROPERTY_SIMPLE, OPERAND_TWO_ARRAY
    use solution_data_module, only: Temp_High_Rad, rad_st, &
          & Temp_Low_Rad, rad_int,  integral_th_mean, tau_th_mean, &
          & integral_th, npr_temp, tau_th, npr_st, rad_temp, tau_th_st
    use solver_data_module, only : if_cylindrical, &
          & if_different_ref_temp, Temp_reference_tube, &
          & Temp_reference_oxide, if_oxide_and_tube_temp

    ! Argument List
    integer, intent(IN)  :: N, id_high_or_low
    logical, intent(IN)  :: if_ave_scale

    ! Local Variables

    real   :: DI_tau_r_dr, ave_tau, var1, temp_now, temp_ref
    integer :: i, j, j_temp, id_beta_ref, id_beta_now, id, &
             & indent_low, append_high, id1, case_tau, n_mean, n_mean1
    real, dimension(0:moxide) :: th_exp_coeff_loc, th_exp_temp_loc
    character(LEN = 80)       :: operation_type = 'Integral_Var'
    character(LEN = 18)       :: ch_temp1
    logical    :: debug_local = .false., switch_beta
    real, dimension(mgrid) :: current_temp, reference_temp

    ! beta is the thermal expansion coefficient at the temperature location
    real, dimension(moxide, mgrid)  :: beta_ref, beta_now
    ! index where thermal expansion at each layer, thicknesss is
    ! located in the initial table id_beta_ref, id_beta_now

    if (id_high_or_low == id_high .and. &
      & (.not. if_different_ref_temp))  then
      ! high temperature case
      ! do not include thermal expansion effects; since the
      ! the coating is at the reference temperature; i.e., 
      ! grows unrestrained of thermal expansion
      tau_th = 0.0
      integral_th = 0.0
      tau_th_st = 0.0

      return

    endif

  do i = 1, N

    if (i <= no_metal_layers)  then
      ! subtrate or tube
      id = i

      if (if_different_ref_temp)  then

        ! reference temperature for the tube is the room temperature
        reference_temp = Temp_reference_tube

        if (id_high_or_low == id_high)  then
          current_temp = Temp_High_Rad(i, :)
        else if (id_high_or_low == id_low)  then
          current_temp = Temp_Low_Rad(i, :)
        endif

      else

        ! this is the classical approach
        ! this case is valid only if id_high_or_low == id_low
        ! reference temperature is the high temperature
        reference_temp = Temp_High_Rad(i, :)
        current_temp = Temp_Low_Rad(i, :)

      endif

    else if (i >= first_oxide_layer)  then

      ! these are the oxide layers

      if (if_ave_scale)  then
        ! use average properties for the oxide scale
        id = no_oxide_ave
      else
        id = i
      endif

      if (if_oxide_and_tube_temp)  then

        reference_temp = Temp_reference_oxide

      else 

        ! Temp_High_Rad and Temp_Low_Rad are not stored in the no_oxide_ave
        ! but the layer counting used by the solver
        ! reference temperature for the oxide is the high temperature
        reference_temp = Temp_High_Rad(i, :)

      endif

      if (id_high_or_low == id_high)  then
        current_temp = Temp_High_Rad(i, :)
      else if (id_high_or_low == id_low)  then
        current_temp = Temp_Low_Rad(i, :)
      endif

    endif

    do j = 1, npr_temp(i)

      ! get beta and its corresponding index; reference temperature
      call LINEAR_PROPERTY_SIMPLE (reference_temp(j), &
           & nth_exp(id), 0, moxide, th_exp_coeff(id, :), &
           & th_exp_temp(id, :), beta_ref(i, j), id_beta_ref)

      ! get beta and its corresponding index; low temperature
      call LINEAR_PROPERTY_SIMPLE (current_temp(j), &
           & nth_exp(id), 0, moxide, th_exp_coeff(id, :), &
           & th_exp_temp(id, :), beta_now(i, j), id_beta_now)

      ! the aim is to compute (Tc-Tref) * alfa; even when Tc < Tref
      ! input data is given for ascending order for T
      ! thus compute all the time the thermal expansion but keep track
      ! if the current temperature was switched with the reference temp

      ! thermal expansion dummy array 
      if (id_beta_ref == id_beta_now)  then

        ! (Tc-Tref) * alfa (for the same segment)
        ! all the time integrate from low to high such that 
        ! we have expansion
        if (reference_temp(j) < current_temp(j))  then

           tau_th(i, j) = 0.5 * (current_temp(j) - &
             & reference_temp(j)) * (beta_ref(i, j) + beta_now(i, j))
           switch_beta = .false.  

        else if (reference_temp(j) >= current_temp(j))  then

           tau_th(i, j) = 0.5 * (reference_temp(j) - &
             & current_temp(j)) * (beta_ref(i, j) + beta_now(i, j))
           switch_beta = .true.  

        endif
 
      else 

        ! switch the data such that all the time "low" < "high"
        ! switch data only if ref > current
        if (abs(id_beta_ref) > abs(id_beta_now))  then

          ! now current state is stored in high data structures

          id1 = id_beta_ref
          id_beta_ref = id_beta_now
          id_beta_now = id1

          var1 = beta_ref(i, j)
          beta_ref(i, j) = beta_now(i, j)
          beta_now(i, j) = var1

          ! temperatures must be switched, too
          temp_ref = current_temp(j)
          temp_now = reference_temp(j) 

          switch_beta = .true.

        else

          temp_ref = reference_temp(j) 
          temp_now = current_temp(j)

          switch_beta = .false.

        endif

      if (debug_local)  then

        if ((i == first_oxide_layer .and. j == 3) .or. &
          & (i ==  first_oxide_layer+ 1 .and. j == 3))   then
          ! done to accomodate i=2 and i=3
          write(aux_lun, 34) i, j, id_beta_ref, id_beta_now, &
            & reference_temp(j), current_temp(j), &
            & temp_ref, temp_now
 34       format('index_after_switch_low ', i1, 1x, i2, 1x, 2(i3, 1x), &
            & 10(1pe13.6, 1x))

        endif

      endif

      ! the low index is all the time _ref
        if (id_beta_ref > 0)  then

          th_exp_coeff_loc = th_exp_coeff(id, :)
          th_exp_temp_loc = th_exp_temp(id, :)

          th_exp_coeff_loc(id_beta_ref) = beta_ref(i, j)
          th_exp_temp_loc(id_beta_ref) = temp_ref

          indent_low = 0  ! no data set was inserted

        else if (id_beta_ref < 0)  then

          ! "stop" is not an option when the metal reference temperature 
          ! is room temperature
          ! write(out_lun, *) 'ERROR: Tau: th expansion at low temp not defined' 
          ! stop
          id_beta_ref = - id_beta_ref  ! id_beta_ref = 1

          if (id_beta_ref > 1)  then

            write(tty_lun, 31) i, j, id_beta_now, id_beta_ref, &
                & current_temp(j), reference_temp(j), &
                & beta_now(i, j), beta_ref(i, j), tau_th(i, j), ave_tau

            write(out_lun, *) 'ERROR: Tau: th exp at low temp not defined' 
            stop

          endif

          th_exp_coeff_loc(2:nth_exp(id)+1) = th_exp_coeff(id, 1:nth_exp(id))
          th_exp_temp_loc(2:nth_exp(id)+1) = th_exp_temp(id, 1:nth_exp(id))

          ! insert the first interval with constant thermal expansion
          th_exp_coeff_loc(1) = beta_ref(i, j)
          th_exp_temp_loc(1) = temp_ref

          indent_low = 1  ! one data set was inserted

        endif

        ! _now is the high value
        if (id_beta_now > 0)  then

          append_high = 0  ! do not append any data set at the end

        else if (id_beta_now < 0)  then
          ! outside the intervals
          ! add another interval with constant thermal expansion
          id_beta_now = - id_beta_now ! id_beta_now = nth_exp(i)+1

          ! n_th_exp_i = nth_exp(i) + 1
          append_high = 0  ! append one data set at the end

        endif

        ! consider the new data sets for the first or last data sets
        id_beta_now = id_beta_now + indent_low + append_high + 1

        th_exp_coeff_loc(id_beta_now) = beta_now(i, j)
        th_exp_temp_loc(id_beta_now) = temp_now

        call OPERAND_TWO_ARRAY(0, moxide, th_exp_coeff_loc, &
               & th_exp_temp_loc, id_beta_ref, &
               & id_beta_now, operation_type, &
               & tau_th(i, j))

      endif

      ! switch the data if there was an inversion
      if (switch_beta)  then
        ! Tc < Tref
        ! (Tref - Tc) * alfa was computed
        ! need to change signs
        case_tau = -1
        tau_th(i, j) = - tau_th(i, j)

      else
        ! Tc > Tref 
        ! (Tc-Tref) * alfa was computed
        ! do not need to change signs here
        case_tau = 1

      endif

      ave_tau = -0.5 * (beta_now(i, j)+ beta_ref(i, j)) * &
              & (temp_ref - temp_now)

      if (debug_local)  then

        if ((i == first_oxide_layer .and. j == 3) .or. &
          & (i ==  first_oxide_layer+ 1 .and. j == 3))   then

          write(aux_lun, 31) i, j, id_beta_now, id_beta_ref, &
            & case_tau, current_temp(j), reference_temp(j), &
            & beta_now(i, j), beta_ref(i, j), &
            & tau_th(i, j), ave_tau
 31       format('index_low ', i1, 1x, i2, 1x, 2(i3, 1x), i2, 1x, &
            & 10(1pe13.6, 1x))

        endif

      endif

    end do

  end do

  ! change sign since T reference > T cold
  ! do not change signs here since for one segment (Tc-Tref) * alfa
  ! was computed, i.e., with the right sign
  ! tau_th = - tau_th

  do i = 1, N

    do j = 1, npr_st(i)

      ! temperature point 
      j_temp = 2 * j - 1

      tau_th_st(i, j) = tau_th(i, j_temp)

    end do

  end do

  ! I_i = int(low limit r_i+1, upper limit r) tau x dx 
  !      since r_i+1 < r
  
  ! indices are backwards, i.e., lower indices point to larger radius

  do i = 1, N

    integral_th(i, npr_st(i)) = 0.0

    do j = npr_st(i) - 1, 1, -1

      ! temperature point 
      j_temp = 2 * j

      ! fractional integral; tau_i_r * r * dr
      DI_tau_r_dr = tau_th(i, j_temp) * rad_temp(i, j_temp) * &
                  & (rad_st(i, j) - rad_st(i, j+1))

      if (debug_local)  then

        write(aux_lun, 33) i, j, j_temp, rad_st(i, j), rad_st(i, j+1), &
        & rad_st(i, j) - rad_st(i, j+1), rad_temp(i, j_temp), &
        & tau_th(i, j_temp), DI_tau_r_dr
 33     format('dintegral ', i1, 2(1x, i2), 1x, 10(1pe13.6, 1x))

      endif

      ! update the radial integral
      integral_th(i, j) = integral_th(i, j + 1) + DI_tau_r_dr

    end do

  end do

  ! find integral_th_mean at rad_mean location; 
  ! careful as the stress points may be even or odd numbers
  do i = 1, N

    n_mean = INT(npr_st(i)/2)
    if (npr_st(i) == 2 * n_mean)  then
      ! npr_st is an even number
      n_mean1 = INT(npr_st(i)/2)
      integral_th_mean(i) = (integral_th(i, n_mean) + &
              & integral_th(i, n_mean+1)) / 2.0

    else

      n_mean1 = INT((npr_st(i)+1)/2)
      if (npr_st(i) + 1 == 2 * n_mean1)  then
        ! npr_st is odd
        integral_th_mean(i) = integral_th(i, n_mean1)
      endif

    endif

    ! get the average thermal expansion
    n_mean = npr_st(i)  ! for the temperatures
    tau_th_mean(i) = tau_th(i, n_mean)

  end do

  if (id_high_or_low == id_high)  then
    ch_temp1 = 'integralth_ox_hig '
  else if (id_high_or_low == id_low)  then
    ch_temp1 = 'integralth_ox_low '
  endif
  write(out_lun, 9) ch_temp1, DI_tau_r_dr-DI_tau_r_dr, &
    & (rad_int(first_oxide_layer) - rad_int(N+1)) * 1.0e+6, &
    & (integral_th(i, 1), i = first_oxide_layer, N)
 9 format(a18, 1x, 50(1pe13.6, 1x))

  if (debug_local)  then
  
  do i = 1, N
    do j = 1, npr_st(i)
      j_temp = 2 * j - 1
      write(aux_lun, 32) i, j, rad_st(i, j), &
        & Temp_Low_Rad(i, j_temp), Temp_High_Rad(i, j_temp), &
        & beta_now(i, j_temp), beta_ref(i, j_temp), &
        & tau_th_st(i, j), integral_th(i, j)
 32   format('integral_low ', i1, 1x, i2, 1x, 10(1pe13.6, 1x))
    end do
  end do

  end if

  return

  END SUBROUTINE BETA_TAU_TH

  REAL FUNCTION TP(cylinder_or_plane, r)
 
  ! function for the temperature profile for a cylinder or plane geometry
  real,   Intent(IN)  :: r
  logical, Intent(IN) :: cylinder_or_plane

  if (cylinder_or_plane)  then
    TP = LOG(ABS(r))
  else
    TP = r
  endif

  return

  END FUNCTION TP
  
  REAL FUNCTION TEMP_VAL(cylinder_or_plane, r, Tin, ac, bc)

  ! temperature profile for a cylinder or plane geometry
  real,   Intent(IN)  :: r, Tin, ac, bc
  logical, Intent(IN) :: cylinder_or_plane
    
    TEMP_VAL = Tin + ac + bc * TP(cylinder_or_plane, r)

  return

  END FUNCTION TEMP_VAL

  REAL FUNCTION DTP_DR(cylinder_or_plane, r)
 
  ! derivative of the TP, dTP/dr
  ! function for the temperature profile for a cylinder or plane geometry
  real,   Intent(IN)  :: r
  logical, Intent(IN) :: cylinder_or_plane

  if (cylinder_or_plane)  then
    DTP_DR = 1.0/r
  else
    DTP_DR = 1.0
  endif

  return

  END FUNCTION DTP_DR

  REAL FUNCTION Int_TP_R_DR(cylinder_or_plane, r)
 
  ! integral of the TP * r * dr, i.e., d(Int_TP_R_DR)/dr = TP * r
  ! function for the temperature profile for a cylinder or plane geometry
  real,   Intent(IN)  :: r
  logical, Intent(IN) :: cylinder_or_plane

  if (cylinder_or_plane)  then
    Int_TP_R_DR = 0.5 * (LOG(ABS(r)) - 0.5) * r**2
  else
    Int_TP_R_DR = (1.0/3.0) * r**3
  endif

  return

  END FUNCTION Int_TP_R_DR

END MODULE THERMAL_CYCLE_MODULE
