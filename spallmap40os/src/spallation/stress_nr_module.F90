MODULE STRESS_NR_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for STRESS profile using Newton-Rhapson solver. 
    !
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: STRESS_PROFILE_NR, SOLVER_DELTA, SETUP_STRESS_MATRIX, &
     & SOLVER_VAR, START_CREEP_STRAIN_MEAN

  CONTAINS

 SUBROUTINE SOLVER_VAR(iter, cycle_now, time_now, dt, convergent)

    !=======================================================================
    ! Purpose(s):
    !
    !  solve linear system A*x=rhs; vm_n = vm(sr, sh, sz), test convergence 
    !  
    !
    !=======================================================================

    use solver_data_module, only       : dTemp_eps, fval, var, &
        & delta_var, rhs, coeff, update_sola_type
    use solver_module,  only           : SOLVER_SIMPLE
    use output_module,        only: aux_lun, blank_line,   &
                                  & out_lun, Output_String, prefix, &
                                  & tty_lun, WRITE_STRING
    use solution_data_module, only: no_stress_var, eq_vm, i_s_vm, &
                                  & vm_stress_mean, i_cr_r, i_cr_h, i_cr_z
    use creep_data_module, only: stress_eq_epsilon, i_cr_start, &
                                  & i_cr_end
    use stress_util_nr_module, only: UPDATE_VAR

    implicit none

    ! arguments
    integer,    intent(IN)  :: iter, cycle_now
    logical,    intent(INOUT) :: convergent
    real,    intent(IN)  :: dt, time_now

    ! local variables
    logical, dimension(no_stress_var)  :: Mask_all, Mask_eq_conv, &
         & Mask_var_conv
    integer            :: neq_solver, i, n_norm
    integer, dimension(1)           :: func_id_max, delta_id_max
    real               :: func_ave, func_max, delta_ave, delta_max, &
         & dvar_over_var, norm_lhs
    ! real(kind = 8), dimension(no_stress_var, no_stress_var)  :: matrix
    ! real(kind = 8), dimension(no_stress_var)  :: solution
    
    convergent = .false.
    neq_solver = no_stress_var

    Mask_all = .false.
    do i = 1, neq_solver
      Mask_all(i) = .true.
    end do

    Mask_var_conv = .false.
    Mask_eq_conv = .false.
    do i = i_cr_start, i_cr_end
      Mask_var_conv(i_s_vm(i)) = .true.
      Mask_eq_conv(eq_vm(i)) = .true.
    end do    

    ! solver returns delta_var; the increment
    call SOLVER_SIMPLE(no_stress_var, coeff, rhs)

    if (i_cr_start > i_cr_end) then
      convergent = .true.
      RETURN
    endif

    call UPDATE_VAR(dt, rhs)

    ! rhs holds the new solution, 
    delta_var = rhs - var

    ! update the solution
    var = rhs

    write(aux_lun, 4) iter, time_now, (var(i_s_vm(i)), i = i_cr_start, i_cr_end), &
         & (delta_var(i_s_vm(i)), i = i_cr_start, i_cr_end)
 4  format('vm_check ', i3, 1x, 50(1pe13.6, 1x))

    write(aux_lun, 14) iter, time_now, (var(i_cr_h(i)), i = i_cr_start, i_cr_end), &
         & (delta_var(i_cr_h(i)), i = i_cr_start, i_cr_end)
 14  format('cr_h_check ', i3, 1x, 50(1pe13.6, 1x))

    if (iter > 1)  then
      dvar_over_var = MAXVAL(abs(delta_var / var), Mask_var_conv .and. &
         & abs(var) > 1.0e-12)
      ! dvar_over_var2 = SUM((1.0-(() / var)**2), Mask_var_conv)

    else
      dvar_over_var = 1.0
    endif

    if (iter > 1)  then
      norm_lhs = 0.0
      n_norm = 0
  
      ! evaluate the equivalent creep strain
      do i = i_cr_start, i_cr_end

        if (ABS(var(i_s_vm(i))) > 1.0e-13 .and. &
           &  ABS( vm_stress_mean(i))  > 1.0e-13)  then

          n_norm = n_norm + 1

          norm_lhs = norm_lhs + (1.0 - vm_stress_mean(i) / &
                 & var(i_s_vm(i)))**2

        endif

      end do

      if (n_norm > 0)  then

        ! norm_rhs = stress_eq_epsilon
        norm_lhs = SQRT(norm_lhs / n_norm)

      else

        norm_lhs = 0.0

      endif

    else

      norm_lhs = 3.0

    endif

    ! check for convergence
    delta_ave = SUM(abs(delta_var), Mask_var_conv)
    delta_max = MAXVAL(abs(delta_var), Mask_var_conv)
    delta_id_max = MAXLOC(abs(delta_var), Mask_var_conv)
    
    ! previous convergence criteria
    ! if (func_max < dTemp_eps .or. dvar_over_var < 1.0e-4)  &

    if (cycle_now <=10 .and. iter == 10)  then

      if (n_norm == 0)  convergent = .true.

    else

      if (norm_lhs < stress_eq_epsilon .and. dvar_over_var < 1.0e-4)  &
         & convergent = .true.

    endif

    if (norm_lhs < stress_eq_epsilon .or. dvar_over_var < 1.0e-4)  &
         & convergent = .true.

    if (convergent)  then

      Output_String = blank_line
      write (Output_String, 1) time_now, iter, delta_id_max, &
       & norm_lhs, stress_eq_epsilon, delta_max, dvar_over_var
 1   format('time= ', 1pe13.6, ' conv_it= ', i3, ' did= ', i3, &
      & ' norm= ', 1pe13.6, ' eps= ', 1pe13.6, &
      & ' dx_max=', 1pe13.6, ' max(dx/x)=', 1pe13.6)
      call WRITE_STRING (Output_String, aux_lun)  ! tty_lun, out_lun)

    else if (.not. convergent)  then

      Output_String = blank_line
      write (Output_String, 2) iter, delta_id_max, &
       & norm_lhs, stress_eq_epsilon, delta_max, dvar_over_var
 2   format('it= ', i3, ' did= ', i3, &
      & ' norm= ', 1pe13.6, ' eps= ', 1pe13.6, &
      & ' dx_max=', 1pe13.6, ' max(dx/x)=', 1pe13.6)
      call WRITE_STRING (Output_String, aux_lun)  ! tty_lun, out_lun)

    endif

    return

 END SUBROUTINE SOLVER_VAR

 SUBROUTINE SOLVER_DELTA(iter, cycle_now, time_now, dt, convergent)
    !=======================================================================
    ! Purpose(s):
    !
    !  solve linear system A*dx=-fval, xn = x + dx; test convergence 
    !  update_sola_type = 3 - full blown Newton Rhapson
    !
    !=======================================================================

    use solver_data_module, only       : dTemp_eps, fval, var, &
        & delta_var, rhs, coeff, update_sola_type
    use solver_module,  only           : SOLVER_SIMPLE
    use output_module,        only: aux_lun, blank_line,   &
                                  & out_lun, Output_String, prefix, &
                                  & tty_lun, WRITE_STRING
    use solution_data_module, only: no_stress_var, eq_vm, i_s_vm, &
                                  & vm_stress_mean, i_s_r, i_s_h, i_cr_h, &
                                  & i_cr_r
    use creep_data_module, only: stress_eq_epsilon, i_cr_start, &
                                  & i_cr_end

    implicit none

    ! arguments
    integer,    intent(IN)  :: iter, cycle_now
    real, intent(IN)        :: time_now, dt
    logical,    intent(INOUT) :: convergent

    ! local variables
    logical, dimension(no_stress_var)  :: Mask_all, Mask_eq_conv, &
         & Mask_var_conv
    integer            :: neq_solver, i, n_norm
    integer, dimension(1)           :: func_id_max, delta_id_max
    real               :: func_ave, func_max, delta_ave, delta_max, &
         & dvar_over_var, norm_lhs, var_max, func_max_all
    ! real(kind = 8), dimension(no_stress_var, no_stress_var)  :: matrix
    ! real(kind = 8), dimension(no_stress_var)  :: solution
    
    convergent = .false.
    neq_solver = no_stress_var

    Mask_all = .false.
    do i = 1, neq_solver
      Mask_all(i) = .true.
    end do

    Mask_var_conv = .false.
    Mask_eq_conv = .false.
    do i = i_cr_start, i_cr_end
      Mask_var_conv(i_s_vm(i)) = .true.
      Mask_eq_conv(eq_vm(i)) = .true.
    end do    

    ! obtain the average and max of the absolute value for the function
    func_ave = SUM(abs(fval), Mask_eq_conv)
    func_max = MAXVAL(abs(fval), Mask_eq_conv)
    func_id_max = MAXLOC(abs(fval), Mask_eq_conv)

    func_max_all = MAXVAL(abs(fval))

    delta_var = -fval

    ! solver returns delta_var; the increment
    call SOLVER_SIMPLE(no_stress_var, coeff, delta_var)

    ! update solution
    var = var + delta_var

    var_max = MAXVAL(abs(var), Mask_var_conv)

    if (iter > 1)  then
      dvar_over_var = MAXVAL(abs(delta_var / var), Mask_var_conv .and. &
         & abs(var) > 1.0e-12)
      dvar_over_var = MAXVAL(abs(delta_var / var),  abs(var) > 1.0e-12)

    else
      dvar_over_var = 1.0
    endif

    if (iter > 1)  then
      norm_lhs = 0.0
      n_norm = 0
  
      ! evaluate the equivalent creep strain
      do i = i_cr_start, i_cr_end

        if (ABS(var(i_s_vm(i))) > 1.0e-13 .and. &
           &  ABS( vm_stress_mean(i))  > 1.0e-13)  then

          n_norm = n_norm + 1

          norm_lhs = norm_lhs + (1.0 - vm_stress_mean(i) / &
                 & var(i_s_vm(i)))**2

        endif

      end do

      if (n_norm > 0)  then

        ! norm_rhs = stress_eq_epsilon
        norm_lhs = SQRT(norm_lhs / n_norm)

      else

        norm_lhs = 1.0

      endif

    else

      norm_lhs = 3.0

    endif

    ! check for convergence
    delta_ave = SUM(abs(delta_var), Mask_var_conv)
    delta_max = MAXVAL(abs(delta_var), Mask_var_conv)
    delta_id_max = MAXLOC(abs(delta_var), Mask_var_conv)
    
    ! previous convergence criteria
    ! if (func_max < dTemp_eps .or. dvar_over_var < 1.0e-4)  &
    ! if (norm_lhs < stress_eq_epsilon .or. dvar_over_var < 1.0e-4)  &
    if (func_max_all < 1.0e-4 .and. dvar_over_var < 1.0e-4)  &
         & convergent = .true.

      Output_String = blank_line
      write (Output_String, 61) iter, func_id_max, delta_id_max, &
       & func_max, delta_max, dvar_over_var
 61   format('it= ', i3, ' fid= ', i3, ' did= ', i3, &
      & ' fmax= ', 1pe13.6, ' dx_max=', 1pe13.6, ' max(dx/x)=', 1pe13.6)
      call WRITE_STRING (Output_String, aux_lun) !  tty_lun, out_lun)

    if (convergent)  then

      Output_String = blank_line
      write (aux_lun, 1) time_now, cycle_now, iter, &
       & var_max, var(i_cr_h(1)), norm_lhs, &
       & delta_id_max, delta_max, dvar_over_var
 1   format(1pe13.6, ' conv_cycle_it= ', i5, 1x, i3,&
      & ' vm_max_ecr_h= ', 1pe13.6, 1x, 1pe13.6, ' norm= ', 1pe13.6, &
      & ' did= ', i3, &
      & ' dx_max=', 1pe13.6, ' max(dx/x)=', 1pe13.6)
   !   call WRITE_STRING (Output_String, aux_lun)  ! tty_lun, out_lun)

    else if (.not. convergent)  then

      Output_String = blank_line
      write (Output_String, 2) iter, delta_id_max, &
       & norm_lhs, stress_eq_epsilon, delta_max, dvar_over_var
 2   format('it= ', i3, ' did= ', i3, &
      & ' norm= ', 1pe13.6, ' eps= ', 1pe13.6, &
      & ' dx_max=', 1pe13.6, ' max(dx/x)=', 1pe13.6)
      call WRITE_STRING (Output_String, aux_lun)  ! tty_lun, out_lun)
 
      if (var_max > 1.0e+9)  then

        Output_String = blank_line
        write(Output_String, 4) iter, time_now, &
         & (var(i_s_vm(i)), i = i_cr_start, i_cr_end), &
         & (var(i_s_r(i)), i = i_cr_start, i_cr_end), &
         & (var(i_s_h(i)), i = i_cr_start, i_cr_end)
 4      format('diverge_vm_r_h ', i3, 1x, 50(1pe13.6, 1x))
        call WRITE_STRING (Output_String, aux_lun, tty_lun, out_lun)

        Output_String = blank_line
        write(Output_String, 14) iter, time_now, &
         & (var(i_cr_h(i)), i = i_cr_start, i_cr_end), &
         & (var(i_cr_r(i)), i = i_cr_start, i_cr_end)
 14  format('diverge_ecr_h_r ', i3, 1x, 50(1pe13.6, 1x))
        call WRITE_STRING (Output_String, aux_lun, tty_lun, out_lun)

        STOP

      endif

    endif

    return

 END SUBROUTINE SOLVER_DELTA


  SUBROUTINE STRESS_PROFILE_NR(N, cycle_no, nr_iter, no_creep_intervals, &
       & dt, time_now, p_out, p_in, id_high_or_low, if_ave_scale, &
       & if_update_displ_dim, convergent)

    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the stress-strain profiles
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================
  use parameter_module, only: id_low, id_high
  use stress_utility_module, only: GET_OG_DIM_DISPL
  use creep_data_module, only: creep_strain_old, creep_strain_new, &
      & creep_sum, creep_dif, creep_integral_sum, creep_integral_dif
  use creep_data_module, only: no_intervals_now, operation_id_current, &
      & interval_op_now
  use solution_data_module, only: ib, ic, i_cr_r, i_cr_h, i_cr_z, i_s_r, &
      & i_s_h, i_s_z, i_s_vm, i_cr, ieps0, no_stress_var, &
      & mean_creep_layer, no_oxide_layers
  use solution_data_module, only: hoop_stress_mean, ax_stress_mean, &
      & rad_stress_mean, vm_stress_mean, creep_eq_incr_mean, &
      & bc_st, cc_st, eps0
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int, integral_th, tau_th, tau_th_st, eps0, bc_st, cc_st, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low
  use solver_data_module, only: var, delta_var, update_sola_type, &
      & var_conv
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use creep_data_module, only: i_cr_start, i_cr_end,&
      & creep_strain_mean, creep_strain_old, creep_strain_new, &
      & stress_eq, stress_eq_old, creep_integral_z
  use creep_module,      only: OUT_CREEP_CONV
  use oxide_utility_module, only : SET_POISSON_YOUNGS
  use stress_solution_module, only: GET_STRESS_STRAIN
  use stress_util_nr_module, only: VAR_CONV_SET

  ! Argument List
  integer,    intent(IN)  :: N, id_high_or_low, cycle_no, nr_iter, &
       & no_creep_intervals
  real,    intent(IN)  :: dt, p_out, p_in, time_now
  logical, intent(IN)  :: if_ave_scale, if_update_displ_dim, &
       & convergent

  ! Local Variables
  logical  :: if_debug_print
  integer  :: i, j, n_norm
  integer, save  :: ip = 0
  real, dimension(N)   :: poisson, Youngs
  real, dimension(N)   :: th_exp, poisson_star, th_exp_star, &
       & Youngs_star
  real, dimension(N+1) :: r, r2, r2dif, r1log
  real, dimension(N)   :: Integral
  real                 :: norm_lhs, norm_rhs

  ! begin setup of star* quantities
  if_debug_print = .false.
  ! check if print out at different oxide thicknesses
  do i = 1, no_full2low_output

    ! print hoop stress strain profile at different oxide thicknesses
    if (cycle_no == no_point_full2low(i) .or. &
        & cycle_no - 1 == no_point_full2low(i) )  then
      if_debug_print = .true.
    endif

  end do

  ! cycle_no - counts all the transition points not the number of pulses
  ! start_creep_cycles counts the number of pulses
  ! the end of first pulse has cycle_no = 4 and second cycle_no = 8
  ! do not consider stresses due to thermal expansion for the first cycle
  ! such that creep can be initiated
  if (cycle_no <= 1)  then  ! this will be in effect only for > 4
    integral_th = 0.0
    tau_th_st = 0.0
    tau_th = 0.0
  endif
 
  ! initialize the radius at the interface between layers
  r(1:N+1) = rad_int(1:N+1)
  r2 = r**2
  r2dif(1:N) = 0.5 * (r2(1:N) - r2(2:N+1))
  r1log = LOG(ABS(r(1:N) / r(2:N+1)))

  ! also good for average layer
  call SET_POISSON_YOUNGS(N, if_ave_scale, poisson, Youngs)

  do i = 1, N
    ! integral_th was computed at the "stress" points
      ! this is to account for the use of alfa* not alfa to 
      ! integrate the thermal expansion, Eq. 6.2 in Nada
    Integral(i) = (1.0 + poisson(i)) * integral_th(i, 1)
    th_exp(i) = 1.0  ! th_exp is factored in the tau_th and Integral_th
  end do
  
  ! define star quantities for plain strain problems
  ! this is good; since the tau did not contain the (1+poisson) factor
  ! this is not needed now; it is accounted in the 
  th_exp_star = th_exp * (1 + poisson)  ! is considered above

  poisson_star = poisson / (1.0 - poisson)
  Youngs_star = Youngs / (1.0 - poisson**2)

  ! END setup of star* quantities

  ! scale down the solution
  ! fval = fval !  / scale_solver

  ! store the results; this may not need it since VAR_CONV_SET is used
  eps0 = var(1)
  do i = 1, N
    bc_st(i) = var(ib(i))
    cc_st(i) = var(ic(i))

    if (i >= i_cr_start .and. i <= i_cr_end)  then
      creep_strain_mean(i)%rad = var(i_cr_r(i))
      creep_strain_mean(i)%hoop = var(i_cr_h(i))
      creep_strain_mean(i)%ax = var(i_cr_z(i))
      vm_stress_mean(i) = var(i_s_vm(i))
      rad_stress_mean(i) = var(i_s_r(i))
      hoop_stress_mean(i) = var(i_s_h(i))
      ax_stress_mean(i) = var(i_s_z(i))
      creep_eq_incr_mean(i) = var(i_cr(i))
    endif

  end do

 ! this was enough when the only 2-3 and 4-1 was considered
  ! stress_eq_conv(id_high_or_low, :, :) = stress_eq
  if (operation_id_current == 3 .or. &
     & operation_id_current == 7)  then

    ! store the vm stress at full load only at the last interval
    if (interval_op_now == no_intervals_now(operation_id_current))  then
      ! stress_eq_conv(id_high, :, :) = stress_eq
      call VAR_CONV_SET(id_high, N)  

    endif

  else if (operation_id_current == 1 .or. &
     & operation_id_current == 5)  then

    ! store the vm stress at partial load only at the last interval
    if (interval_op_now == no_intervals_now(operation_id_current))  then
      ! stress_eq_conv(id_low, :, :) = stress_eq
      call VAR_CONV_SET(id_low, N)  

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
      ! stress_eq_conv(2+interval_op_now, :, :) = stress_eq
      call VAR_CONV_SET(2+interval_op_now, N)  

    endif

  ! do not assign the von mises at the ramp up since these will be used
  ! in reverse order

  endif  

  if (if_debug_print)  then

    write(aux_lun, 5) cycle_no, eps0, (bc_st(i), i= 1, N)
 5 format('const_sol_e_b=', i5, 1x, 30(1pe13.6, 1x))
    write(aux_lun, 6) cycle_no, eps0, (cc_st(i), i= 1, N)
 6 format('const_sol_e_c=', i5, 1x, 30(1pe13.6, 1x))

  endif

  ! update the OLD creep strain when average creep was used
  creep_strain_old%rad = creep_strain_new%rad
  creep_strain_old%hoop = creep_strain_new%hoop
  creep_strain_old%ax = creep_strain_new%ax

  ! updated OLD von Mises stress
  stress_eq_old = stress_eq

  ! update the NEW variables
  do i = 1, N

    if (i >= i_cr_start .and. i <= i_cr_end)  then
      creep_strain_new(i, :)%rad = creep_strain_mean(i)%rad
      creep_strain_new(i, :)%hoop = creep_strain_mean(i)%hoop
      creep_strain_new(i, :)%ax = creep_strain_mean(i)%ax
      stress_eq(i, :) = vm_stress_mean(i) 
    endif

    creep_sum(i, :) = creep_strain_new(i, :)%rad + &
       & creep_strain_new(i, :)%hoop + &
       & 2.0 * poisson(i) * creep_strain_new(i, :)%ax
    creep_dif(i, :) = creep_strain_new(i, :)%rad - &
       & creep_strain_new(i, :)%hoop

    creep_integral_z(i) = creep_strain_new(i, 1)%ax * &
       & (r2(i) - r2(i+1)) / 2.0

    do j = 1, npr_st(i)
      ! fractional integral; F * r * dr
      creep_integral_sum(i, j) = creep_sum(i, 1) * &
       & (rad_st(i, j)**2 - r2(i+1)) / 2.0
      
      ! fractional integral; (F / r) * dr
      creep_integral_dif(i, j) = creep_sum(i, 1) * &
       & LOG(rad_st(i, j) / r(i+1))

    end do

  end do

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

    ! if (norm_lhs < stress_eq_epsilon)  then  !  norm_rhs)  then convergent

  else

    ! write(2, *) 'ERROR: no von Mises points for convergence norm'
    ! STOP
    write(2, 7) time_now, nr_iter, no_creep_intervals
 7    format('stress_eq_zero ', 1pe13.6, 2(1x, i2), 1x, 50(1pe13.6, 1x))

    ! RETURN

  endif

  call OUT_CREEP_CONV(nr_iter, no_creep_intervals, N, &
      & cycle_no, time_now, dt, id_high_or_low, norm_lhs, norm_rhs, n_norm, &
      & convergent)

  call GET_STRESS_STRAIN(N, cycle_no, poisson, th_exp, Youngs, &
       & poisson_star, th_exp_star, Youngs_star)

  ! get displacements and dimensions
  call GET_OG_DIM_DISPL(N, poisson, th_exp, Youngs, &
       & poisson_star, th_exp_star, Youngs_star, if_update_displ_dim)

  return

  END SUBROUTINE STRESS_PROFILE_NR

  SUBROUTINE START_CREEP_STRAIN_MEAN(iter_st_creep, cycle_no, N, &
       & id_high_or_low, dt)

    !=======================================================================
    ! Purpose(s):
    !
    !   assign variables at start of creep-stress iterations
    !=======================================================================
 
  use parameter_module, only: id_low, id_high
  use creep_data_module, only: creep_strain_old, creep_strain_new, &
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
  use stress_util_nr_module, only: VAR_INIT2CONV

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
     ! not ok since old^k = new^k and new^k updated in STRESS_PROFILE_NR was used
     ! creep_strain_new%rad = creep_strain_old%rad 
     ! creep_strain_new%hoop = creep_strain_old%hoop 
     ! creep_strain_new%ax = creep_strain_old%ax 

     ! var(i_cr_r(i)), i_cr_h, i_cr_z do not need to be updated since
     ! the creep strain is accumulating in time and var is all the time stored
 
     ! just to force at least two creep iterations 
     ! since criteria of convergence will depend on |stress_eq - stress_eq_old|
     ! also the von Mises for full and partial load are different
     ! stress_eq_old = 0.0
     
     if (no_pulse_per_cycle == 1 .and. cycle_no <= 4)  then 

       if (id_high_or_low == id_low)  then
         ! previous stress level was higher
         ! stress_eq_old = stress_eq / 3.3
         call VAR_INIT2CONV(id_high_or_low, N)  

       else if (id_high_or_low == id_high)  then
         ! previous stress level was low
         ! stress_eq_old = stress_eq / 1.3
         call VAR_INIT2CONV(id_high_or_low, N)  

       endif

       RETURN

     endif

     if (no_pulse_per_cycle == 2 .and. cycle_no <= 2)  then 

       if (id_high_or_low == id_low)  then
         ! previous stress level was higher
         ! stress_eq_old = stress_eq / 3.3
         call VAR_INIT2CONV(id_high_or_low, N)  

       else if (id_high_or_low == id_high)  then
         ! previous stress level was low
         ! stress_eq_old = stress_eq / 1.3
         call VAR_INIT2CONV(id_high_or_low, N)  

       endif

       RETURN

     endif

     ! the following is for larger, regular cycles i.e., cycle_no > 2(4)

     if (operation_id_current == 3 .or. &
       & operation_id_current == 7)  then

       ! stress_eq_old = stress_eq_conv(id_high_or_low, :, :)
       call VAR_INIT2CONV(id_high_or_low, N)  

     else if (operation_id_current == 1 .or. &
       & operation_id_current == 5)  then

       ! stress_eq_old = stress_eq_conv(id_high_or_low, :, :)
       call VAR_INIT2CONV(id_high_or_low, N)  

     else if (operation_id_current == 4 .or. &  ! f2l transition
         & operation_id_current == 8)  then

       if (interval_op_now == &
           & no_intervals_now(operation_id_current))  then
         ! stress_eq_old = stress_eq_conv(id_low, :, :)
         call VAR_INIT2CONV(id_low, N)  

       else
         ! stress_eq_old = stress_eq_conv(2+interval_op_now, :, :)
         call VAR_INIT2CONV(2+interval_op_now, N)  

       endif

     else if (operation_id_current == 2 .or. &  ! f2l transition
         & operation_id_current == 6)  then

       if (interval_op_now == &
           & no_intervals_now(operation_id_current))  then
         ! stress_eq_old = stress_eq_conv(id_high, :, :)
         call VAR_INIT2CONV(id_high, N)  

       else
         id = 2 + no_intervals_now(operation_id_current) - interval_op_now
         ! stress_eq_old = stress_eq_conv(id, :, :)
         call VAR_INIT2CONV(id, N)  

       endif

     endif

     ! creep_integral_sum and creep_integral_dif remain unchanged
     ! make sure they are not reset

     RETURN

   else if (iter_st_creep > 1)  then

     ! prepare for the next iteration
     ! update_var will take care of it
 
   endif

  RETURN

  END SUBROUTINE START_CREEP_STRAIN_MEAN

  SUBROUTINE SETUP_STRESS_MATRIX(N, cycle_no, nr_iter, dt, p_out, p_in, &
         & id_high_or_low, if_ave_scale)

    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the stress-strain profiles
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================
  ! Integral(i) = integral of tau_i(r) * r dr (limits are r_i to r)
  ! unknowns are in this order
  ! eps0, b1, b2, .. bN, c1, c2, ... cN
  ! setting coeff for d, keeping eps0 as unknown
  ! N = total number of layers; minimum 2 (substrate + average scale)
  ! get stress constants
  ! get stress and strain
 
  use solver_data_module, only  : generalized_plane_strain
  use solver_data_module, only: fval, var, &
        & delta_var, rhs, coeff, update_sola_type
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int, integral_th, tau_th, tau_th_st, eps0, bc_st, cc_st, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low
  use solution_data_module, only: ib, ic, i_cr_r, i_cr_h, i_cr_z, i_s_r, &
      & i_s_h, i_s_z, i_s_vm, i_cr, ieps0, no_stress_var, &
      & mean_creep_layer, no_oxide_layers
  use oxide_data_module,   only: Youngs_modul, Youngs_temp, &
      & th_exp_coeff, th_exp_temp, no_oxide, no_layer, &
      & nth_exp, nyoungs, ncond, ncp_prop, id_mat, &
      & poisson_ratio, no_oxide_ave, pilling_bedworth_ratio, &
      & first_oxide_layer
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & u_now_int_jump, ox_strain_r, ox_strain_hoop, &
      & ox_strain_axial, u_int_jump_now 
  use oxide_stress_data_module,  only: oxidation_strain_mode, &
      & growth_ox_location
  use induced_oxide_stress_module, only: STRAIN_OX_IND_TERM_RADIAL, &
       & U_TERM_STRAIN_OX_IND, STRESS_AXIAL_OX_IND_TERM
  use creep_data_module, only: creep_strain_new, &
       & creep_integral_sum, creep_integral_dif, creep_integral_z, &
       & if_creep_ox_scale, if_creep_metal, i_cr_start, i_cr_end
  use creep_module, only: CREEP_TERM_S_RAD, CREEP_TERM_U, &
       & CREEP_TERM_INT_AXIAL
  use oxide_utility_module, only: SET_POISSON_YOUNGS
  use stress_util_nr_module, only: STRESS_MEAN_JACOBIAN

  ! Argument List
  integer,    intent(IN)  :: N, id_high_or_low, cycle_no, nr_iter
  real,    intent(IN)  :: dt, p_out, p_in
  logical, intent(IN)  :: if_ave_scale
  ! if_update_displ_dim = true only for oh algorithm

  ! Local Variables
  integer  :: i, ieq, j, icr
  integer, save  :: ip = 0
  real     :: Y_alfa_i, Y_alfa_ip1, rf, pf_eps, pf_b, Youngs_ratio, &
        & scale_solver, ox_term_i, ox_term_ip1, r_lnf, creep_term_ip1, &
        & rhs_cr, rhs_th, rhs_ox, rhs_u, a1, b1, rhs_th_total, &
        & rhs_ox_total, rhs_cr_total, coeff_eps0
  ! system to solve coeff * dx = -f ; or coeff * dx = fval
  real, dimension(no_stress_var)      :: creep_term 
  real, dimension(N)   :: poisson, Youngs
  real, dimension(N)   :: th_exp, poisson_star, th_exp_star, &
       & Youngs_star
  real, dimension(N+1) :: r, r2, r2dif, r1log
  real, dimension(N)   :: Integral
  logical  :: if_debug_print

  ! if_simple_update_sola = .true. then solution algorithm will be just
  ! simple updates of the unknowns, and Newton Rhapson will not be used
  ! update_sola_type = 1 - use simple update
  ! update_sola_type = 2 - use more complex update, introducing
  ! creep strains er, eh, ez, (means stresses) sr0, sh0, sz0 as unknown
  ! while keeping de^cr and vm0 at the previous iteration 
  ! update_sola_type = 3 - full blown Newton Rhapson

  if_debug_print = .false.
  ! check if print out at different oxide thicknesses
  do i = 1, no_full2low_output

    ! print hoop stress strain profile at different oxide thicknesses
    if (cycle_no == no_point_full2low(i) .or. &
        & cycle_no - 1 == no_point_full2low(i) )  then
      if_debug_print = .true.
    endif

  end do

  ! cycle_no - counts all the transition points not the number of pulses
  ! start_creep_cycles counts the number of pulses
  ! the end of first pulse has cycle_no = 4 and second cycle_no = 8
  ! do not consider stresses due to thermal expansion for the first cycle
  ! such that creep can be initiated
  if (cycle_no <= 1)  then  ! this will be in effect only for > 4
    integral_th = 0.0
    tau_th_st = 0.0
    tau_th = 0.0
  endif

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1
 
  ! ic - tracks the creep terms
  icr = 0

  ! the constants are too small; use a scaling constant
  ! to multipy rhs and then to divide the results
  scale_solver = 1.0e+4

  ! initialize the radius at the interface between layers
  r(1:N+1) = rad_int(1:N+1)
  r2 = r**2
  r2dif(1:N) = 0.5 * (r2(1:N) - r2(2:N+1))
  r1log = LOG(ABS(r(1:N) / r(2:N+1)))

  ! also good for average layer
  call SET_POISSON_YOUNGS(N, if_ave_scale, poisson, Youngs)

  do i = 1, N

    ! write(out_lun, *) 'r1log ', i, r1log(i), r2dif(i)
    if (i >= i_cr_start .and. i <= i_cr_end)  then

    endif

    ! integral_th was computed at the "stress" points
      ! this is to account for the use of alfa* not alfa to 
      ! integrate the thermal expansion, Eq. 6.2 in Nada
    Integral(i) = (1.0 + poisson(i)) * integral_th(i, 1)
    th_exp(i) = 1.0  ! th_exp is factored in the tau_th and Integral_th
  end do
  
  ! define star quantities for plain strain problems
  ! this is good; since the tau did not contain the (1+poisson) factor
  ! this is not needed now; it is accounted in the 
  th_exp_star = th_exp * (1 + poisson)  ! is considered above

  poisson_star = poisson / (1.0 - poisson)
  Youngs_star = Youngs / (1.0 - poisson**2)

  ! initialize
  coeff = 0.0
  fval = 0.0
  rhs = 0.0

  ! using fractional integrals
  ! u = r * b + c/r + (1+nu*) * J_ri,r/r + u_ox_term + creep_term_u
  ! s_r = E^* (b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2) - 
  !       E^* (Fr + creep_term_s_r)
  ! s_teta = E^* (b/(1-nu*) - tau + c/(r^2 * (1+ nu*)) + J_ri,r/r^2) -
  !          E^* Fhoop + creep_term_s_hoop

  ! tube or first layer of tube
  i = 1

  ! first equation
  ! pressure on the outside surface of the tube
  ieq = 1
  ! eps0 equations
  coeff(ieq, ieps0) = poisson(i) / (1.0 - poisson_star(i))
  !  coeff(ieq, ieps0) = 0.0

  ! b and c are ok
  coeff(ieq, ib(i)) = 1.0 / (1.0 - poisson_star(i))
  coeff(ieq, ic(i)) = - 1.0 / (r2(i) * (1.0 + poisson_star(i)))

  ! no oxide on the outer surface; no oxidation induced stresses
  rhs_th = -p_out / Youngs_star(i) + Integral(i) / r2(i)

  rhs_ox = 0.0  ! no oxide growth induced terms for the metal

  ! for extreme cases, at high temperature there will be creep on the 
  ! outer surface; consider creep effects for a fraction of the metal
  if (i_cr_start == 1)  then

      ! coefficients for 0.5 * J_i(e_cr_dif) - 0.5 * J_i(e_cr_sum)

      ! d_creep_rad
      coeff(ieq, i_cr_r(i)) = 0.5 * (r1log(i)-r2dif(i) / r2(i))
      ! d_creep_hoop
      coeff(ieq, i_cr_h(i)) = 0.5 * (-r1log(i)-r2dif(i) / r2(i))
      ! d_creep_hoop
      coeff(ieq, i_cr_z(i)) = -0.5 * 2.0 * poisson(i) * r2dif(i) / r2(i)

      rhs_cr = var(i_cr_r(i)) * coeff(ieq, i_cr_r(i)) + &
             & var(i_cr_h(i)) * coeff(ieq, i_cr_h(i)) + &
             & var(i_cr_z(i)) * coeff(ieq, i_cr_z(i))
      ! be careful at the sign; to be consitent with rhs_th
      rhs_cr = - rhs_cr

      ! evaluate the NR function f=A*x - rhs at the end change the sign "-f"
      fval(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
           & var(ib(i)) * coeff(ieq, ib(i)) + &
           & var(ic(i)) * coeff(ieq, ic(i))- rhs_th - rhs_cr - rhs_ox

      if (update_sola_type == 1)  then

        rhs(ieq) = rhs_ox + rhs_th + rhs_cr
        coeff(ieq, i_cr_r(i)) = 0.0  ! since they moved to RHS
        coeff(ieq, i_cr_h(i)) = 0.0 
        coeff(ieq, i_cr_z(i)) = 0.0 

      else if (update_sola_type == 2)  then

        rhs(ieq) = rhs_ox + rhs_th

      endif

      ! get the other equations for the first layer
      call STRESS_MEAN_JACOBIAN(i, ieq, N, nr_iter, dt, &
           & r, r2, poisson, Youngs, poisson_star, Youngs_star)

      ! etr_sum == summ of (strain * r)
      ! eor_dif == (e_r - e_hoop) / r
      ! keep it for future reference
      ! creep_term_ip1 = CREEP_TERM_S_RAD(creep_integral_sum(1, 1), &
      !    & creep_integral_dif(1, 1), r2(1))

      icr = icr + 1
      creep_term(icr) = rhs_cr  !  creep_term_ip1

  else if (i_cr_start > 1)  then

    ! evaluate the NR function f=A*x - rhs at the end change the sign "-f"
    fval(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
        & var(ib(i)) * coeff(ieq, ib(i)) + &
        & var(ic(i)) * coeff(ieq, ic(i))- rhs_th - rhs_ox

     rhs(ieq) = rhs_ox + rhs_th


  endif

  ! conditions for interface between materials 'i' and 'i+1'
  do i = 1, N -1

    ! write(out_lun, *) 'no_var ', i, no_stress_var, ieq

    Youngs_ratio = Youngs_star(i) / Youngs_star(i+1)

    ! second type of equations
    ! continuation of radial stresses
    ! s_r = E^* [ b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2 + &
    !             eps0 * nu/(1-nu*) - Fr(strain_ox, r) - &
    !            Jcr_ri,r(creep_sum*r) / (2*r2) + Jcr_ri,r(creep_dif/r) / 2]

    ieq = ieq + 1
    ! coefficents for eps0
    coeff(ieq, ieps0) = poisson(i+1) / (1.0 - poisson_star(i+1)) - &
       & Youngs_ratio * poisson(i) / (1.0 - poisson_star(i))
    ! coeff(ieq, ieps0) = 0.0

    ! coeff for b
    coeff(ieq, ib(i)) = - Youngs_ratio / (1.0 - poisson_star(i))
    coeff(ieq, ib(i+1)) = 1.0 / (1.0 - poisson_star(i+1))

    ! coeff for c
    coeff(ieq, ic(i)) = Youngs_ratio / ((1.0 + poisson_star(i))* r2(i+1))
    coeff(ieq, ic(i+1)) = - 1.0 / &
                        & ((1.0 + poisson_star(i+1))* r2(i+1))

    ! rhs
    rhs_th = Integral(i+1) / r2(i+1)
      
    if (i+1 >= first_oxide_layer)  then
      ! contribution from 'i+1' layer
      ox_term_ip1 = STRAIN_OX_IND_TERM_RADIAL(poisson(i+1),  &
        & poisson_star(i+1), r(i+1), &
        & ox_strain_r(i+1), ox_strain_hoop(i+1), ox_strain_axial(i+1))

    else
      ox_term_ip1 = 0.0
    endif

    if (i >= first_oxide_layer)  then
      ! contribution from 'i' layer: metal does not have any oxidation strains

      ox_term_i = STRAIN_OX_IND_TERM_RADIAL(poisson(i), poisson_star(i), &
        & r(i+1), ox_strain_r(i), ox_strain_hoop(i), ox_strain_axial(i))

    else
      ox_term_i = 0.0
    endif

    rhs_ox = ox_term_ip1 - Youngs_ratio * ox_term_i

    ! add the creep terms only when needed
    if (i + 1 >= i_cr_start .and. i + 1 <= i_cr_end)  then

        ! coefficients for 0.5 * J_(i+1)(e_cr_dif) - 0.5 * J_(i+1)(e_cr_sum)

        ! d_creep_rad
        coeff(ieq, i_cr_r(i+1)) = 0.5 * (r1log(i+1)-r2dif(i+1) / r2(i+1))
        ! d_creep_hoop
        coeff(ieq, i_cr_h(i+1)) = 0.5 * (-r1log(i+1)-r2dif(i+1) / r2(i+1))
        ! d_creep_hoop
        coeff(ieq, i_cr_z(i+1)) = -0.5 * 2.0 * poisson(i+1) * &
           & r2dif(i+1) / r2(i+1)

        rhs_cr = var(i_cr_r(i+1)) * coeff(ieq, i_cr_r(i+1)) + &
             & var(i_cr_h(i+1)) * coeff(ieq, i_cr_h(i+1)) + &
             & var(i_cr_z(i+1)) * coeff(ieq, i_cr_z(i+1))
        ! be careful at the sign; to be consitent with rhs_th
        rhs_cr = - rhs_cr

        ! evaluate the NR function f=A*x - rhs
        fval(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
             & var(ib(i)) * coeff(ieq, ib(i)) + &
             & var(ib(i+1)) * coeff(ieq, ib(i+1)) + &
             & var(ic(i)) * coeff(ieq, ic(i)) + &
             & var(ic(i+1)) * coeff(ieq, ic(i+1)) - &
             & rhs_ox - rhs_th - rhs_cr

        if (update_sola_type == 1)  then

          rhs(ieq) = rhs_ox + rhs_th + rhs_cr
          coeff(ieq, i_cr_r(i+1)) = 0.0  ! since they moved to RHS
          coeff(ieq, i_cr_h(i+1)) = 0.0 
          coeff(ieq, i_cr_z(i+1)) = 0.0 

        else if (update_sola_type == 2)  then

          rhs(ieq) = rhs_ox + rhs_th

        endif

        call STRESS_MEAN_JACOBIAN(i+1, ieq, N, nr_iter, &
           & dt, r, r2, poisson, Youngs, poisson_star, Youngs_star)

        ! etr_sum == summ of (strain * r)
        ! eor_dif == (e_r - e_hoop) / r
        ! keep it for future reference
        ! creep_term_ip1 = CREEP_TERM_S_RAD(creep_integral_sum(i +1, 1), &
        !   & creep_integral_dif(i +1, 1), r2(i+1))

        icr = icr + 1
        creep_term(icr) = creep_term_ip1

    else  if (i + 1 < i_cr_start .or. i + 1 > i_cr_end)  then

      ! evaluate the NR function f=A*x - rhs
      fval(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
           & var(ib(i)) * coeff(ieq, ib(i)) + &
           & var(ib(i+1)) * coeff(ieq, ib(i+1)) + &
           & var(ic(i)) * coeff(ieq, ic(i)) + &
           & var(ic(i+1)) * coeff(ieq, ic(i+1)) - &
           & rhs_ox - rhs_th

        rhs(ieq) = rhs_ox + rhs_th

    endif
    
    ! third type of equations
    ! continuation of radial displacement
    ! u = r * b + c/r + (1+nu*) * J_ri,r/r + ox_ind_term_u +
    !     (1+nu*) * Jcr_ri,r(creep_sum*r) / (2*r) + 
    !     (1-nu*) * r * Jcr_ri,r(creep_dif/r) / 2

    ieq = ieq + 1

    ! no coefficients for eps0

    ! coeff for b
    coeff(ieq, ib(i)) = - r(i+1)
    coeff(ieq, ib(i+1)) = r(i+1)

    ! coeff for c
    coeff(ieq, ic(i)) = - 1.0 / r(i+1)
    coeff(ieq, ic(i+1)) = 1.0 / r(i+1)

    ! rhs
    rhs_th = - Integral(i+1) * (1.0 + poisson_star(i+1)) / r(i+1)
    
    ! previously envisioned for Oh (2006) method
    ! if (i == 2 .and. if_growth_induce_stress)  then
    !  fval(ieq) = fval(ieq) - u_now_int_jump(id_high_or_low)
    ! endif

    if (i >= first_oxide_layer)  then

      ox_term_i = U_TERM_STRAIN_OX_IND(poisson_star(i), &
        & r(i+1), ox_strain_r(i), ox_strain_hoop(i))

    else
      ox_term_i = 0.0
    endif

    if (i + 1 >= first_oxide_layer)  then

      ox_term_ip1 = U_TERM_STRAIN_OX_IND(poisson_star(i+1), r(i+1), &
        & ox_strain_r(i+1), ox_strain_hoop(i+1))

      ! at metal-oxide interface metal does not have any oxidation strains

    else
      ox_term_ip1 = 0.0
    endif

    rhs_ox = ox_term_i - ox_term_ip1

    ! check the sign for u_now_int_jump(id_high_or_low) 
    ! + sign for u_inner - u_outer = U_delta
    ! fval(ieq) = fval(ieq) + ox_term_i + u_now_int_jump(id_high_or_low)
    ! u_int_jump_now(1, i) - interface between i and 'i-1' layers
    rhs_u = u_int_jump_now(1, i+1)

    if (ABS(u_int_jump_now(1, i+1)) > 0.0)  then
 
      write(aux_lun, *) 'u_jump ', u_int_jump_now(1, i+1)

    endif

    ! add the creep terms only when needed
    if (i + 1 >= i_cr_start .and. i + 1 <= i_cr_end)  then

        ! (1+mu*) * J_(i+1)(e_cr_sum) / r + (1-mu*) * J_(i+1)(e_cr_dif) * r
        a1 = (1.0 + poisson_star(i+1)) * r2dif(i+1) / r(i+1)
        b1 = (1.0 - poisson_star(i+1)) * r1log(i+1) * r(i+1)
        coeff(ieq, i_cr_r(i+1)) = (a1 + b1) / 2.0
        coeff(ieq, i_cr_h(i+1)) = (a1 - b1) / 2.0
        coeff(ieq, i_cr_z(i+1)) = poisson(i+1) * a1
        rhs_cr = var(i_cr_r(i+1)) * coeff(ieq, i_cr_r(i+1)) + &
             & var(i_cr_h(i+1)) * coeff(ieq, i_cr_h(i+1)) + &
             & var(i_cr_z(i+1)) * coeff(ieq, i_cr_z(i+1))
        ! be careful at the sign; to be consitent with rhs_th
        rhs_cr = - rhs_cr

        ! evaluate the NR function f=A*x - rhs
        fval(ieq) = var(ib(i)) * coeff(ieq, ib(i)) + &
             & var(ib(i+1)) * coeff(ieq, ib(i+1)) + &
             & var(ic(i)) * coeff(ieq, ic(i)) + &
             & var(ic(i+1)) * coeff(ieq, ic(i+1))- &
             & rhs_ox - rhs_th - rhs_cr - rhs_u

        if (update_sola_type == 1)  then

          rhs(ieq) = rhs_ox + rhs_th + rhs_cr + rhs_u
          coeff(ieq, i_cr_r(i+1)) = 0.0  ! since they moved to RHS
          coeff(ieq, i_cr_h(i+1)) = 0.0 
          coeff(ieq, i_cr_z(i+1)) = 0.0 

        else if (update_sola_type == 2)  then

            rhs(ieq) = rhs_ox + rhs_th + rhs_u

        endif

      ! etr_sum == summ of (strain * r)
      ! eor_dif == (e_r - e_hoop) / r
      ! keep it for future reference
      ! creep_term_ip1 = CREEP_TERM_U(creep_integral_sum(i +1, 1), &
      !    & creep_integral_dif(i +1, 1), poisson_star(i+1), r(i+1))

      icr = icr + 1
      ! creep_term(icr) = creep_term_ip1

    else if (i + 1 < i_cr_start .or. i + 1 > i_cr_end)  then

      ! evaluate the NR function f=A*x - rhs
      fval(ieq) = var(ib(i)) * coeff(ieq, ib(i)) + &
           & var(ib(i+1)) * coeff(ieq, ib(i+1)) + &
           & var(ic(i)) * coeff(ieq, ic(i)) + &
           & var(ic(i+1)) * coeff(ieq, ic(i+1))- &
           & rhs_ox - rhs_th - rhs_u
      rhs(ieq) = rhs_ox + rhs_th + rhs_u

    endif

  end do

  ! second equation
  ! pressure on the inside surface of the last oxide layer (on the steam side)
  ieq = ieq + 1
  coeff(ieq, ieps0) = poisson(N) / (1.0 - poisson_star(N))
  ! coeff(ieq, ieps0) = 0.0

  coeff(ieq, ib(N)) = 1.0 / (1.0 - poisson_star(N))
  coeff(ieq, ic(N)) = - 1.0 / (r2(N+1) * (1.0 + poisson_star(N)))

  rhs_th = -p_in / Youngs_star(N)
  rhs_ox = STRAIN_OX_IND_TERM_RADIAL(poisson(N), poisson_star(N), &
        & r(N+1), ox_strain_r(N), ox_strain_hoop(N), ox_strain_axial(N))
  
  ! lucky on the internal surface, that of the oxide surface
  ! creep integrals are ZERO on this first location (the lower integration limit)
  rhs_cr = 0.0

  ! evaluate the NR function f=A*x - rhs
  fval(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
         & var(ib(N)) * coeff(ieq, ib(N)) + &
         & var(ic(N)) * coeff(ieq, ic(N))- rhs_th - rhs_cr - rhs_ox

  rhs(ieq) = rhs_ox + rhs_th + rhs_cr

  ! the last equation; axial balance or set eps0 to zero
  ieq = ieq + 1

  if (generalized_plane_strain)  then

    ! write(tty_lun, *) 'FIX IT according to new I_i (computed from the inner radius' 
    ! stop

    ! initialize
    rhs_th_total = 0.0
    rhs_ox_total = 0.0
    rhs_cr_total = 0.0

    do i = 1, N

      ! radius factor 
      rf = 0.5 * (r2(i) - r2(i+1))  ! r2dif

      ! poisson factor
      pf_eps = (1.0 - poisson(i)) / ((1.0 + poisson(i)) * &
             & (1.0 - 2.0 * poisson(i)))
      pf_b = 2.0 * poisson(i) / (1.0 - poisson_star(i))

      coeff_eps0 = rf * pf_eps * Youngs(i)

      coeff(ieq, ieps0) = coeff(ieq, ieps0) + coeff_eps0

      coeff(ieq, ib(i)) = rf * pf_b * Youngs_star(i)

      rhs_th = (Youngs(i) / (1.0 + poisson(i)) + &
               & Youngs_star(i) * poisson(i)) * Integral(i)
      rhs_th_total = rhs_th_total + rhs_th

      if (i >= first_oxide_layer)  then
        ! only for the oxide layers

        ! A - term
        rhs_ox = rf * STRESS_AXIAL_OX_IND_TERM( &
            & poisson(i), poisson_star(i), Youngs(i), Youngs_star(i), &
            & ox_strain_r(i), ox_strain_hoop(i), ox_strain_axial(i))

        ! B - term
        r_lnf = 0.5 * (r2(i) * LOG(ABS(r(i))) - r2(i+1) * LOG(ABS(r(i+1))))

        rhs_ox = rhs_ox - Youngs_star(i) * poisson(i) * &
            & (ox_strain_r(i) - ox_strain_hoop(i)) * &
            & (r_lnf - 0.5 * rf)

        rhs_ox_total = rhs_ox_total + rhs_ox

      endif

      ! add the creep terms only when needed
      if (i >= i_cr_start .and. i <= i_cr_end)  then

          ! 0.5*A4*(J_i(e_cr_sum)-r1,i^2*J_i(e_cr_dif)) + A6*Int(r*e_cr_z)
          a1 = poisson(i) * Youngs_star(i) * r2dif(i) / 2.0
          b1 = poisson(i) * Youngs_star(i) * r1log(i) * r2(i) / 2.0

          ! d_creep_rad
          coeff(ieq, i_cr_r(i)) = -(a1 - b1)
          ! d_creep_hoop
          coeff(ieq, i_cr_h(i)) = -(a1 + b1)
          ! d_creep_hoop
          coeff(ieq, i_cr_z(i)) = - 2.0 * poisson(i) * a1 - &
            & Youngs(i) * r2dif(i)
          rhs_cr = var(i_cr_r(i)) * coeff(ieq, i_cr_r(i)) + &
             & var(i_cr_h(i)) * coeff(ieq, i_cr_h(i)) + &
             & var(i_cr_z(i)) * coeff(ieq, i_cr_z(i))
          ! be careful at the sign; to be consitent with rhs_th
          rhs_cr = - rhs_cr
          rhs_cr_total = rhs_cr_total + rhs_cr

          ! evaluate the NR function f=A*x - rhs
          fval(ieq) = fval(ieq) + var(ieps0) * coeff_eps0 + &
             & var(ib(i)) * coeff(ieq, ib(i)) - &
             & rhs_th - rhs_cr - rhs_ox

          ! keep it for future reference
          ! creep_term_ip1 = CREEP_TERM_INT_AXIAL(Youngs(i), poisson(i), &
          !   & creep_integral_sum(i, 1), creep_integral_dif(i, 1), &
          !   & creep_integral_z(i), r2(i))
          !  rhs_cr_total = rhs_cr_total + creep_term_ip1

          ! icr = icr + 1
          ! creep_term(icr) = creep_term_ip1

          if (update_sola_type == 1)  then
            coeff(ieq, i_cr_r(i)) = 0.0  ! since they moved to RHS
            coeff(ieq, i_cr_h(i)) = 0.0 
            coeff(ieq, i_cr_z(i)) = 0.0 
          endif

          if (i == N)  then

            if (update_sola_type == 1)  then
              rhs(ieq) = rhs_th_total + rhs_ox_total + rhs_cr_total
            else if (update_sola_type == 2)  then
              rhs(ieq) = rhs_th_total + rhs_ox_total
            endif

          endif

      else if (i  < i_cr_start .or. i  > i_cr_end)  then

        fval(ieq) = fval(ieq) + var(ieps0) * coeff_eps0 + &
             & var(ib(i)) * coeff(ieq, ib(i)) - &
             & rhs_th - rhs_ox

        rhs(ieq) = rhs(ieq) + rhs_th + rhs_ox

      endif

    end do

  else if (.not. generalized_plane_strain)  then
    
    coeff(ieq, ieps0) = 1.0
    fval(ieq) = 0.0
    rhs(ieq) = 0.0

  endif

  if (if_debug_print)  then

    do i = 1, no_stress_var
      write(aux_lun, 2)  cycle_no, i, (coeff(i, j), &
         & j = 1, no_stress_var), fval(i)
 2  format('sfs ', i5, 1x, i2, 1x, 30(1pe13.6, 1x))
    end do

    if (icr == 0)  then
      write(aux_lun, *) 'no creep terms ', cycle_no
    else

      write(aux_lun, 7)  cycle_no, (creep_term(i), i = 1, icr)
 7  format('creep_term_st ', i5, 1x, 30(1pe13.6, 1x))

    endif

  endif

  return

  END SUBROUTINE SETUP_STRESS_MATRIX

END MODULE STRESS_NR_MODULE
