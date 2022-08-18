MODULE STRESS_UTIL_NR_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for setting up the Jacobian and RHS for Newton-Rhapson
    !     stress solver
    !
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public ::  STRESS_MEAN_JACOBIAN, VAR_INIT2CONV, VAR_CONV_SET, UPDATE_VAR

  CONTAINS

 SUBROUTINE VAR_CONV_SET(store, N)  

  use creep_data_module, only: i_cr_start, i_cr_end
  use solver_data_module, only: var, var_conv
  use solution_data_module, only: ib, ic, i_cr_r, i_cr_h, i_cr_z, i_s_r, &
      & i_s_h, i_s_z, i_s_vm, i_cr, ieps0, no_stress_var

  ! creep strains are not assigned as they accumulate

  ! Argument List
  integer,    intent(IN)  :: N, store

  ! Local Variables
  integer  :: i, j, n_norm

  var_conv(store, 1) = var(1)
  do i = 1, N

    var_conv(store, ib(i)) = var(ib(i))
    var_conv(store, ic(i)) = var(ic(i))

    if (i >= i_cr_start .and. i <= i_cr_end)  then

      var_conv(store, i_s_vm(i)) = var(i_s_vm(i))
      var_conv(store, i_s_r(i)) = var(i_s_r(i))
      var_conv(store, i_s_h(i)) = var(i_s_h(i))
      var_conv(store, i_s_z(i)) = var(i_s_z(i))

    endif

  end do

  return

 END SUBROUTINE VAR_CONV_SET

 SUBROUTINE VAR_INIT2CONV(retrieve, N)  

  use creep_data_module, only: i_cr_start, i_cr_end
  use solver_data_module, only: var, var_conv
  use solution_data_module, only: ib, ic, i_cr_r, i_cr_h, i_cr_z, i_s_r, &
      & i_s_h, i_s_z, i_s_vm, i_cr, ieps0, no_stress_var

  ! creep strains are not assigned as they accumulate

  ! Argument List
  integer,    intent(IN)  :: N, retrieve

  ! Local Variables
  integer  :: i, j, n_norm

  var(1) =  var_conv(retrieve, 1)
  do i = 1, N

    var(ib(i)) = var_conv(retrieve, ib(i))
    var(ic(i)) = var_conv(retrieve, ic(i))

    if (i >= i_cr_start .and. i <= i_cr_end)  then

      var(i_s_vm(i)) = var_conv(retrieve, i_s_vm(i))
      var(i_s_r(i)) = var_conv(retrieve, i_s_r(i))
      var(i_s_h(i)) = var_conv(retrieve, i_s_h(i))
      var(i_s_z(i)) = var_conv(retrieve, i_s_z(i))

    endif

  end do

  return

  END SUBROUTINE VAR_INIT2CONV

 SUBROUTINE UPDATE_VAR(dt, sol)

    !=======================================================================
    ! Purpose(s):
    !
    !  update nonlinear variables after the solving the linear system 
    !  A*x=rhs; vm_n = vm(sr, sh, sz)
    !  
    !
    !=======================================================================

    use solver_data_module, only       : var, update_coeff, update_const, &
        & update_sola_type, jeps0, jb_i, jc_i, j_s_r, j_s_h
    use output_module,        only: aux_lun, blank_line,   &
                                  & out_lun, Output_String, prefix, &
                                  & tty_lun, WRITE_STRING
    use solution_data_module, only: no_stress_var, eq_vm, i_s_vm, &
                                  & vm_stress_mean
    use solution_data_module, only: ib, ic, i_cr_r, i_cr_h, i_cr_z, i_s_r, &
      & i_s_h, i_s_z, i_s_vm, i_cr, ieps0, no_stress_var, &
      & mean_creep_layer, rad_mean, integral_th_mean, tau_th_mean, &
      & Temp_mean, eq_vm
    use creep_data_module, only: i_cr_start, i_cr_end, creep_strain_mean
    use creep_module, only:  VON_MISES, DCREEP

    implicit none

    ! arguments
    real,    intent(IN)  :: dt
    real(kind = 8),  dimension(no_stress_var), intent(INOUT) :: sol

    ! local variables
    integer            :: j, i
    
    ! update von mises and creep strain increment only if update_sola_type =    ! update_sola_type = 1 - use simple update
    ! update_sola_type = 2 - use more complex update, introducing
    ! creep strains er, eh, ez, (means stresses) sr0, sh0, sz0 as unknown
    ! while keeping de^cr and vm0 at the previous iteration 
 
  if (update_sola_type == 1)  then

    do i = i_cr_start, i_cr_end

      j = i_s_r(i)
      sol(j) = sol(ieps0) * update_coeff(j, jeps0) + &
             & sol(ib(i)) * update_coeff(j, jb_i) + &
             & sol(ic(i)) * update_coeff(j, jc_i) + &
             & update_const(j)
      j = i_s_h(i)
      sol(j) = sol(ieps0) * update_coeff(j, jeps0) + &
             & sol(ib(i)) * update_coeff(j, jb_i) + &
             & sol(ic(i)) * update_coeff(j, jc_i) + &
             & update_const(j)
      j = i_s_z(i)
      sol(j) = sol(ieps0) * update_coeff(j, jeps0) + &
             & sol(i_s_r(i)) * update_coeff(i_s_z(i), j_s_r) + &
             & sol(i_s_h(i)) * update_coeff(i_s_z(i), j_s_h) + &
             & update_const(j)

      j = i_s_vm(i)
      sol(j) = VON_MISES(sol(i_s_r(i)), sol(i_s_h(i)), sol(i_s_z(i)))
      j = i_cr(i)
      sol(j) = DCREEP(sol(i_s_vm(i)), dt, Temp_mean(i), i)

    write(aux_lun, 14) i, sol(i_cr(i)), sol(i_s_vm(i)), Temp_mean(i)
 14  format('cr_vm_temp_check2 ', i2, 1x, 50(1pe13.6, 1x))

      j = i_cr_r(i)
      if (ABS(sol(i_s_vm(i))) > 1.0e-10)  then
        sol(j) = creep_strain_mean(i)%rad + sol(i_cr(i)) * &
           & 0.5 * (2.0 * sol(i_s_r(i)) - &
           & sol(i_s_h(i)) - sol(i_s_z(i))) / sol(i_s_vm(i))
      else
        sol(j) = creep_strain_mean(i)%rad
      endif

      j = i_cr_h(i)
      if (ABS(sol(i_s_vm(i))) > 1.0e-10)  then
        sol(j) = creep_strain_mean(i)%hoop + sol(i_cr(i)) * &
           & 0.5 * (2.0 * sol(i_s_h(i)) - &
           & sol(i_s_r(i)) - sol(i_s_z(i))) / sol(i_s_vm(i))
      else
        sol(j) = creep_strain_mean(i)%hoop

      endif

      j = i_cr_z(i)
      if (ABS(sol(i_s_vm(i))) > 1.0e-10)  then
        sol(j) = creep_strain_mean(i)%ax + sol(i_cr(i)) * &
           & 0.5 * (2.0 * sol(i_s_z(i)) - &
           & sol(i_s_r(i)) - sol(i_s_h(i))) / sol(i_s_vm(i))
      else
        sol(j) = creep_strain_mean(i)%ax
      endif

    end do

  else if (update_sola_type == 2)  then

    do i = i_cr_start, i_cr_end

      j = i_s_vm(i)
      sol(j) = VON_MISES(sol(i_s_r(i)), sol(i_s_h(i)), sol(i_s_z(i)))
      j = i_cr(i)
      sol(j) = DCREEP(sol(i_s_vm(i)), dt, Temp_mean(i), i)

    end do    
    
  endif

  return

 END SUBROUTINE UPDATE_VAR

  SUBROUTINE STRESS_MEAN_JACOBIAN(i, ieq, N, nr_iter, &
    & dt, r, r2, poisson, Youngs, poisson_star, Youngs_star)

  use solver_data_module, only: generalized_plane_strain
  use solver_data_module, only: fval, var, &
      & delta_var, rhs, coeff, update_sola_type, update_coeff, &
      & update_const, jeps0, jb_i, jc_i, j_s_r, j_s_h
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int, integral_th, tau_th, tau_th_st, eps0, bc_st, cc_st, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low
  use solution_data_module, only: ib, ic, i_cr_r, i_cr_h, i_cr_z, i_s_r, &
      & i_s_h, i_s_z, i_s_vm, i_cr, ieps0, no_stress_var, &
      & mean_creep_layer, rad_mean, integral_th_mean, tau_th_mean, &
      & Temp_mean, eq_vm
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
       & if_creep_ox_scale, if_creep_metal, i_cr_start, i_cr_end,&
       & creep_strain_mean

  use creep_module, only:  VON_MISES, DCREEP, DCREEP_DVM

  ! Argument List
  integer,    intent(IN)     :: i, N, nr_iter
  integer,    intent(INOUT)  :: ieq
  real,     intent(IN)       :: dt
  real, intent(IN), dimension(N)   :: poisson, Youngs, poisson_star, &
       & Youngs_star
  real, intent(IN), dimension(N+1) :: r2, r

  ! Local Variables
  real     :: r2mean, r2dif_m, r1log_m, ox_term_i, rhs_ox, rhs_cr, &
    & rhs_th_hoop, ox_term_i_hoop, rhs_ox_hoop, creep_increment_eq, &
    & rhs_th, rhs_cr_hoop
  integer  :: j, icr, ieq_s_r, ieq_s_r0

  r2mean = rad_mean(i)**2
  r1log_m = LOG(ABS(rad_mean(i) / r(i+1)))
  r2dif_m = 0.5 * (r2mean - r2(i+1))

  ieq_s_r = ieq

  ! sigma_r0 = sigma_r(r0) fourth type of equations
  ! f = sigma_r(r0) - sigma_r0
  ieq = ieq + 1
  ieq_s_r0 = ieq

  coeff(ieq, ib(i)) = Youngs_star(i) / (1.0 - poisson_star(i))
  ! coeff(ieq, ib(i)) = Youngs_star(i) * coeff(ieq_s_r, ib(i)) cannot be
  ! used since I have various solver options
  ! it cannot be Youngs_star(i) * coeff(ieq_s_r, ieps0) since at interface
  ! it will not work well
  ! coeff(ieq, ieps0) = Youngs_star(i) * poisson(i) * coeff(ieq_s_r, ib(i))
  coeff(ieq, ieps0) = Youngs_star(i) * poisson(i) / &
        & (1.0 - poisson_star(i))
  coeff(ieq, ic(i)) = -Youngs_star(i) / (r2mean * (1.0 + poisson_star(i)))

  rhs_th = Youngs_star(i) * (1.0 + poisson(i)) * &
       & integral_th_mean(i)  / r2mean

  if (i >= first_oxide_layer)  then
    ! contribution from 'i' layer: metal does not have any oxidation strains

    ox_term_i = STRAIN_OX_IND_TERM_RADIAL(poisson(i), poisson_star(i), &
      & rad_mean(i), ox_strain_r(i), ox_strain_hoop(i), ox_strain_axial(i))
    rhs_ox = Youngs_star(i) * ox_term_i

  else

    ox_term_i = 0.0
    rhs_ox = 0.0

  endif

  ! mean_creep_layer(i) is assumed to be .true.
  ! coefficients for 0.5 * J_i(e_cr_dif) - 0.5 * J_i(e_cr_sum)
  ! d_creep_rad
  coeff(ieq, i_cr_r(i)) = 0.5 * Youngs_star(i) * &
     & (r1log_m - r2dif_m / r2mean)
  ! d_creep_hoop
  coeff(ieq, i_cr_h(i)) = 0.5 * Youngs_star(i) * &
     & (-r1log_m - r2dif_m / r2mean)
  ! d_creep_hoop
  coeff(ieq, i_cr_z(i)) = -0.5 * Youngs_star(i) * &
     & 2.0 * poisson(i) * r2dif_m / r2mean

  rhs_cr = var(i_cr_r(i)) * coeff(ieq, i_cr_r(i)) + &
         & var(i_cr_h(i)) * coeff(ieq, i_cr_h(i)) + &
         & var(i_cr_z(i)) * coeff(ieq, i_cr_z(i))
  ! be careful at the sign; to be consitent with rhs_th
  rhs_cr = - rhs_cr
  
  coeff(ieq, i_s_r(i)) = -1.0

  ! evaluate the NR function f=A*x - rhs
  fval(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
         & var(ib(i)) * coeff(ieq, ib(i)) + &
         & var(ic(i)) * coeff(ieq, ic(i)) + &
         & var(i_s_r(i)) * coeff(ieq, i_s_r(i)) - &
         & rhs_th - rhs_cr - rhs_ox

  if (update_sola_type == 1)  then

    coeff(ieq, i_s_r(i)) = 1.0
    rhs(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
         & var(ib(i)) * coeff(ieq, ib(i)) + &
         & var(ic(i)) * coeff(ieq, ic(i)) - (rhs_ox + rhs_th + rhs_cr)

    write(aux_lun, 16) i, integral_th_mean(i), &
       & var(i_s_vm(i)), var(i_s_r(i)), var(i_s_h(i)), var(i_s_z(i))
 16  format('th_vm_sr_sh_sz_check1 ', i2, 1x, 50(1pe13.6, 1x))

    ! store the coefficients for the update step of the solver
    update_coeff(i_s_r(i), jeps0) = coeff(ieq, ieps0)
    update_coeff(i_s_r(i), jb_i) = coeff(ieq, ib(i))
    update_coeff(i_s_r(i), jc_i) = coeff(ieq, ic(i))
    update_const(i_s_r(i)) = - (rhs_ox + rhs_th + rhs_cr)

    coeff(ieq, ieps0) = 0.0
    coeff(ieq, ib(i)) = 0.0
    coeff(ieq, ic(i)) = 0.0
    coeff(ieq, i_cr_r(i)) = 0.0  ! since they moved to RHS
    coeff(ieq, i_cr_h(i)) = 0.0 
    coeff(ieq, i_cr_z(i)) = 0.0 

  else if (update_sola_type == 2)  then

    rhs(ieq) = rhs_ox + rhs_th

  endif

  ! sigma_hoop0 = sigma_hoop(r0) fifth type of equations
  ! f = sigma_hoop(r0) - sigma_hoop0

  ieq = ieq + 1
  coeff(ieq, ib(i)) = Youngs_star(i) / (1.0 - poisson_star(i))
 ! coeff(ieq_s_r0, ib(i))
  coeff(ieq, ieps0) = Youngs_star(i) * poisson(i) / &
        & (1.0 - poisson_star(i))  ! coeff(ieq_s_r0, ieps0)
  coeff(ieq, ic(i)) = Youngs_star(i) / (r2mean * (1.0 + poisson_star(i)))
        ! - coeff(ieq_s_r0, ic(i))

  rhs_th_hoop = Youngs_star(i) * (1.0 + poisson(i)) * &
     & tau_th_mean(i) - rhs_th

  if (i >= first_oxide_layer)  then
    ! contribution from 'i' layer: metal does not have any oxidation strains

    ox_term_i_hoop = -(ox_strain_r(i) - ox_strain_hoop(i)) / 2.0
    rhs_ox_hoop = Youngs_star(i) * (ox_term_i - ox_term_i_hoop)

  else

    rhs_ox_hoop = 0.0

  endif

  ! coefficients for 
  ! e_cr_hoop + mu*e_cr_z - 0.5 * J_i(e_cr_dif) - 0.5 * J_i(e_cr_sum)
  ! d_creep_rad
  coeff(ieq, i_cr_r(i)) = 0.5 * Youngs_star(i) * &
     & (r1log_m + r2dif_m / r2mean)
  ! d_creep_hoop
  coeff(ieq, i_cr_h(i)) = Youngs_star(i) * (-1.0 + 0.5 * &
     & (-r1log_m + r2dif_m / r2mean))
  ! d_creep_hoop
  coeff(ieq, i_cr_z(i)) = Youngs_star(i) * poisson(i) * &
     & (0.5 * 2.0 *  r2dif_m / r2mean - 1.0)

  rhs_cr_hoop = var(i_cr_r(i)) * coeff(ieq, i_cr_r(i)) + &
         & var(i_cr_h(i)) * coeff(ieq, i_cr_h(i)) + &
       & var(i_cr_z(i)) * coeff(ieq, i_cr_z(i))
  ! be careful at the sign; to be consitent with rhs_th
  rhs_cr_hoop = - rhs_cr_hoop
    
  coeff(ieq, i_s_h(i)) = -1.0

  ! evaluate the NR function f=A*x - rhs
  fval(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
         & var(ib(i)) * coeff(ieq, ib(i)) + &
         & var(ic(i)) * coeff(ieq, ic(i)) + &
         & var(i_s_h(i)) * coeff(ieq, i_s_h(i)) - &
         & rhs_th_hoop - rhs_cr_hoop - rhs_ox_hoop

   if (update_sola_type == 1)  then

     coeff(ieq, i_s_h(i)) = 1.0
     rhs(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
         & var(ib(i)) * coeff(ieq, ib(i)) + &
         & var(ic(i)) * coeff(ieq, ic(i)) - &
         & (rhs_ox_hoop + rhs_th_hoop + rhs_cr_hoop)

     ! store the coefficients for the update step of the solver
     update_coeff(i_s_h(i), jeps0) = coeff(ieq, ieps0)
     update_coeff(i_s_h(i), jb_i) = coeff(ieq, ib(i))
     update_coeff(i_s_h(i), jc_i) = coeff(ieq, ic(i))
     update_const(i_s_h(i)) = - (rhs_ox_hoop + rhs_th_hoop + &
         & rhs_cr_hoop)

     coeff(ieq, ieps0) = 0.0
     coeff(ieq, ib(i)) = 0.0
     coeff(ieq, ic(i)) = 0.0
     coeff(ieq, i_cr_r(i)) = 0.0  ! since they moved to RHS
     coeff(ieq, i_cr_h(i)) = 0.0 
     coeff(ieq, i_cr_z(i)) = 0.0 

   else if (update_sola_type == 2)  then

      rhs(ieq) = rhs_ox_hoop + rhs_th_hoop

   endif

  ! sigma_z0 = sigma_z(r0) six-th type of equation
  ! f = sigma_z(r0) - sigma_z0
  ! f = mu*(sigma_r0+sigma_hoop0) + E*(eps0-eth_mean-e_ox_z-e_cr_z) - sigma_z0
  ieq = ieq + 1
  coeff(ieq, i_s_r(i)) = poisson(i)
  coeff(ieq, i_s_h(i)) = poisson(i)
  coeff(ieq, ieps0) = Youngs(i)
  coeff(ieq, i_cr_z(i)) = - Youngs(i)
  coeff(ieq, i_s_z(i)) = -1.0

  rhs_th = Youngs(i) * tau_th_mean(i) 

  if (i >= first_oxide_layer)  then
    ! contribution from 'i' layer: metal does not have any oxidation strains

    rhs_ox = Youngs(i) * ox_strain_axial(i)

  else

    rhs_ox = 0.0

  endif

  ! evaluate the NR function f=A*x - rhs
  fval(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
       & var(i_s_r(i)) * coeff(ieq, i_s_r(i)) + &
         & var(i_s_h(i)) * coeff(ieq, i_s_h(i)) + &
         & var(i_s_z(i)) * coeff(ieq, i_s_z(i)) + &
         & var(i_cr_z(i)) * coeff(ieq, i_cr_z(i)) - &
         & rhs_th - rhs_ox

  if (update_sola_type == 1)  then

    coeff(ieq, i_s_z(i)) = 1.0
    rhs(ieq) = var(ieps0) * coeff(ieq, ieps0) + &
         & var(i_s_r(i)) * coeff(ieq, i_s_r(i)) + &
         & var(i_s_h(i)) * coeff(ieq, i_s_h(i)) + &
         & var(i_cr_z(i)) * coeff(ieq, i_cr_z(i)) - rhs_th - rhs_ox

    ! store the coefficients for the update step of the solver
    update_coeff(i_s_z(i), jeps0) = coeff(ieq, ieps0)
    update_coeff(i_s_z(i), j_s_r) = coeff(ieq, i_s_r(i))
    update_coeff(i_s_z(i), j_s_h) = coeff(ieq, i_s_h(i))
    update_const(i_s_z(i)) = var(i_cr_z(i)) * coeff(ieq, i_cr_z(i)) - &
        & (rhs_ox + rhs_th)

    coeff(ieq, ieps0) = 0.0
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_cr_z(i)) = 0.0

  else if (update_sola_type == 2)  then

    rhs(ieq) = rhs_ox + rhs_th

  endif

  ! sigma_vm0 = sigma_vm(sigma_r0, sigma_h0, sigma_z0)
  !  seven-th type of equation - the first nonlinear equation
  ! f = sigma_vm(sigma_r0, sigma_h0, sigma_z0) - sigma_vm0
  ieq = ieq + 1
  eq_vm(i) = ieq

  coeff(ieq, i_s_vm(i)) = -1.0

  if (ABS(var(i_s_vm(i))) > 1.0e-10)  then
    coeff(ieq, i_s_r(i)) = 0.5 * (2.0 * var(i_s_r(i)) - &
         & var(i_s_h(i)) - var(i_s_z(i))) / var(i_s_vm(i))
    coeff(ieq, i_s_h(i)) = 0.5 * (2.0 * var(i_s_h(i)) - &
         & var(i_s_r(i)) - var(i_s_z(i))) / var(i_s_vm(i))
    coeff(ieq, i_s_z(i)) = 0.5 * (2.0 * var(i_s_z(i)) - &
         & var(i_s_h(i)) - var(i_s_r(i))) / var(i_s_vm(i))
  else
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_s_z(i)) = 0.0
  endif

  ! evaluate the NR function f= F() - Sigma_vm
  fval(ieq) = VON_MISES(var(i_s_r(i)), var(i_s_h(i)), var(i_s_z(i))) - &
         & var(i_s_vm(i)) 

  if (update_sola_type == 1)  then

    coeff(ieq, i_s_vm(i)) = 1.0
    rhs(ieq) = VON_MISES(var(i_s_r(i)), var(i_s_h(i)), var(i_s_z(i)))
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_s_z(i)) = 0.0

  else if (update_sola_type == 2)  then

    ! use sigma_eq from previous step
    coeff(ieq, i_s_vm(i)) = 1.0
    rhs(ieq) = VON_MISES(var(i_s_r(i)), var(i_s_h(i)), var(i_s_z(i)))
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_s_z(i)) = 0.0

  endif

  ieq = ieq + 1
  ! creep strain increment equation
  ! f = e_cr0 - dt*F(To, VM0)
  coeff(ieq, i_cr(i)) = 1.0
  creep_increment_eq = DCREEP(var(i_s_vm(i)), dt, Temp_mean(i), i)
  coeff(ieq, i_s_vm(i)) = - DCREEP_DVM(var(i_s_vm(i)), Temp_mean(i), &
         & creep_increment_eq, i)
  fval(ieq) = var(i_cr(i)) * coeff(ieq, i_cr(i)) - creep_increment_eq
  
  if (update_sola_type == 1)  then

    rhs(ieq) = creep_increment_eq
    coeff(ieq, i_s_vm(i)) = 0.0

    write(aux_lun, 14) i, var(i_cr(i)), creep_increment_eq, &
       & var(i_s_vm(i)), var(i_s_r(i)), var(i_s_h(i)), var(i_s_z(i))
 14  format('cr_cr_vm_sr_sh_sz_check1 ', i2, 1x, 50(1pe13.6, 1x))
 
  else if (update_sola_type == 2)  then

    rhs(ieq) = creep_increment_eq
    coeff(ieq, i_s_vm(i)) = 0.0

  endif

  ieq = ieq + 1
  ! creep strain increment in the radial direction
  ! f = e_cr - e_cr0 - 1.5*(sigma_r0-p0) * de_cr / sigma0
  coeff(ieq, i_cr_r(i)) = 1.0

  if (ABS(var(i_s_vm(i))) > 1.0e-10)  then
    coeff(ieq, i_s_r(i)) = -var(i_cr(i)) / var(i_s_vm(i))
    coeff(ieq, i_cr(i)) = -0.5 * (2.0 * var(i_s_r(i)) - &
       & var(i_s_h(i)) - var(i_s_z(i))) / var(i_s_vm(i))
  else
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_cr(i)) = 0.0
  endif

  coeff(ieq, i_s_h(i)) = -0.5 * coeff(ieq, i_s_r(i))
  coeff(ieq, i_s_z(i)) = -0.5 * coeff(ieq, i_s_r(i))
  coeff(ieq, i_s_vm(i)) = coeff(ieq, i_cr(i)) * coeff(ieq, i_s_r(i))
  ! make sure that creep_strain_mean(i)%rad is the value at 
  ! the previous time step; not the value at previous iteration
  fval(ieq) = coeff(ieq, i_cr_r(i)) * var(i_cr_r(i)) + &
     & coeff(ieq, i_cr(i)) * var(i_cr(i)) - creep_strain_mean(i)%rad

  if (update_sola_type == 1)  then

    rhs(ieq) = creep_strain_mean(i)%rad - coeff(ieq, i_cr(i)) * var(i_cr(i))
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_cr(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_s_z(i)) = 0.0
    coeff(ieq, i_s_vm(i)) = 0.0

  else if (update_sola_type == 2)  then

    rhs(ieq) = creep_strain_mean(i)%rad - coeff(ieq, i_cr(i)) * var(i_cr(i))
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_cr(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_s_z(i)) = 0.0
    coeff(ieq, i_s_vm(i)) = 0.0

  endif

  ieq = ieq + 1
  ! creep strain increment in the hoop direction
  ! f = e_cr - e_cr0 - 1.5*(sigma_h0-p0) * de_cr / sigma0
  coeff(ieq, i_cr_h(i)) = 1.0

  if (ABS(var(i_s_vm(i))) > 1.0e-10)  then
    coeff(ieq, i_s_h(i)) = -var(i_cr(i)) / var(i_s_vm(i))
    coeff(ieq, i_cr(i)) = -0.5 * (2.0 * var(i_s_h(i)) - &
       & var(i_s_r(i)) - var(i_s_z(i))) / var(i_s_vm(i))
  else
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_cr(i)) = 0.0
  endif

  coeff(ieq, i_s_r(i)) = -0.5 * coeff(ieq, i_s_h(i))
  coeff(ieq, i_s_z(i)) = -0.5 * coeff(ieq, i_s_h(i))
  coeff(ieq, i_s_vm(i)) = coeff(ieq, i_cr(i)) * coeff(ieq, i_s_h(i))
  ! make sure that creep_strain_mean(i)%rad is the value at 
  ! the previous time step; not the value at previous iteration
  fval(ieq) = coeff(ieq, i_cr_h(i)) * var(i_cr_h(i)) + &
     & coeff(ieq, i_cr(i)) * var(i_cr(i)) - creep_strain_mean(i)%hoop

  if (update_sola_type == 1)  then

    rhs(ieq) = creep_strain_mean(i)%hoop - coeff(ieq, i_cr(i)) * &
             & var(i_cr(i))
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_cr(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_s_z(i)) = 0.0
    coeff(ieq, i_s_vm(i)) = 0.0

  else if (update_sola_type == 2)  then

    rhs(ieq) = creep_strain_mean(i)%hoop - coeff(ieq, i_cr(i)) * &
             & var(i_cr(i))
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_cr(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_s_z(i)) = 0.0
    coeff(ieq, i_s_vm(i)) = 0.0

  endif

  ieq = ieq + 1
  ! creep strain increment in the axial direction
  ! f = e_cr - e_cr0 - 1.5*(sigma_z0-p0) * de_cr / sigma0
  coeff(ieq, i_cr_z(i)) = 1.0

  if (ABS(var(i_s_vm(i))) > 1.0e-10)  then
    coeff(ieq, i_s_z(i)) = -var(i_cr(i)) / var(i_s_vm(i))
    coeff(ieq, i_cr(i)) = -0.5 * (2.0 * var(i_s_z(i)) - &
       & var(i_s_h(i)) - var(i_s_r(i))) / var(i_s_vm(i))
  else
    coeff(ieq, i_s_z(i)) = 0.0
    coeff(ieq, i_cr(i)) = 0.0
  endif

  coeff(ieq, i_s_h(i)) = -0.5 * coeff(ieq, i_s_z(i))
  coeff(ieq, i_s_r(i)) = -0.5 * coeff(ieq, i_s_z(i))
  coeff(ieq, i_s_vm(i)) = coeff(ieq, i_cr(i)) * coeff(ieq, i_s_z(i))
  ! make sure that creep_strain_mean(i)%rad is the value at 
  ! the previous time step; not the value at previous iteration
  fval(ieq) = coeff(ieq, i_cr_z(i)) * var(i_cr_z(i)) + &
     & coeff(ieq, i_cr(i)) * var(i_cr(i)) - creep_strain_mean(i)%ax

  if (update_sola_type == 1)  then

    rhs(ieq) = creep_strain_mean(i)%ax - coeff(ieq, i_cr(i)) * var(i_cr(i))
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_cr(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_s_z(i)) = 0.0
    coeff(ieq, i_s_vm(i)) = 0.0

  else if (update_sola_type == 2)  then

    rhs(ieq) = creep_strain_mean(i)%ax - coeff(ieq, i_cr(i)) * var(i_cr(i))
    coeff(ieq, i_s_r(i)) = 0.0
    coeff(ieq, i_cr(i)) = 0.0
    coeff(ieq, i_s_h(i)) = 0.0
    coeff(ieq, i_s_z(i)) = 0.0
    coeff(ieq, i_s_vm(i)) = 0.0

  endif

  RETURN

  END SUBROUTINE STRESS_MEAN_JACOBIAN

END MODULE STRESS_UTIL_NR_MODULE
