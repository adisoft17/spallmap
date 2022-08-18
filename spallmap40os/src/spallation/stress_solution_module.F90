MODULE STRESS_SOLUTION_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for STRESS profile. 
    !
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: STRESS_PROFILE, GET_STRESS_STRAIN

  CONTAINS

  SUBROUTINE STRESS_PROFILE(N, cycle_no, p_out, p_in, &
         & id_high_or_low, if_ave_scale, if_update_displ_dim)

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
  use solver_module, only       : SOLVER_SIMPLE
  use solution_data_module, only: rad_temp, Temp, Temp_ave, npr_temp, &
      & ac_temp, bc_temp, npr_st, rad_st, Temp_low_rad, Temp_high_rad, &
      & rad_int, integral_th, tau_th, tau_th_st, eps0, bc_st, cc_st, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low
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
  use stress_utility_module, only: GET_OG_DIM_DISPL

  ! Argument List
  integer,    intent(IN)  :: N, id_high_or_low, cycle_no
  real,    intent(IN)  :: p_out, p_in
  logical, intent(IN)  :: if_ave_scale, if_update_displ_dim
  ! if_update_displ_dim = true only for oh algorithm

  ! Local Variables
  real, dimension(N)   :: poisson, Youngs
  integer  :: i, ieq, ieps0, j, icr
  integer, dimension(N) :: ib, ic
  integer, save  :: ip = 0
  real     :: Y_alfa_i, Y_alfa_ip1, rf, pf_eps, pf_b, Youngs_ratio, &
        & scale_solver, ox_term_i, ox_term_ip1, r_lnf, creep_term_ip1
  real, dimension(2*N+1, 2*N+1) :: coeff
  real, dimension(2*N+1)      :: rhs, creep_term  
  real, dimension(N)   :: th_exp, &
       & poisson_star, th_exp_star, Youngs_star
  real, dimension(N+1) :: r, r2
  real, dimension(N)   :: Integral
  logical  :: if_debug_print

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

  ! also good for average layer
  call SET_POISSON_YOUNGS(N, if_ave_scale, poisson, Youngs)

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

  ieps0 = 1

  do i = 1, N
    ib(i) = ieps0 + i
  end do
 
  do i = 1, N
    ic(i) = ib(N) + i
  end do

  ! initialize
  coeff = 0.0
  rhs = 0.0

  ! using fractional integrals
  ! u = r * b + c/r + (1+nu*) * J_ri,r/r + u_ox_term + creep_term_u
  ! s_r = E^* (b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2) - 
  !       E^* (Fr + creep_term_s_r)
  ! s_teta = E^* (b/(1-nu*) - tau + c/(r^2 * (1+ nu*)) + J_ri,r/r^2) -
  !          E^* Fhoop + creep_term_s_hoop

  ! first equation
  ! pressure on the outside surface of the tube
  ieq = 1
  ! eps0 equations
  coeff(ieq, ieps0) = poisson(1) / (1.0 - poisson_star(1))
  !  coeff(ieq, ieps0) = 0.0

  ! b and c are ok
  coeff(ieq, ib(1)) = 1.0 / (1.0 - poisson_star(1))
  coeff(ieq, ic(1)) = - 1.0 / (r2(1) * (1.0 + poisson_star(1)))

  ! no oxide on the outer surface; no oxidation induced stresses
  rhs(ieq) = -p_out / Youngs_star(1) + Integral(1) / r2(1)

  ! for extreme cases, at high temperature there will be creep on the 
  ! outer surface; consider creep effects for a fraction of the metal
  if (i_cr_start == 1)  then

    ! etr_sum == summ of (strain * r)
    ! eor_dif == (e_r - e_hoop) / r
    creep_term_ip1 = CREEP_TERM_S_RAD(creep_integral_sum(1, 1), &
          & creep_integral_dif(1, 1), r2(1))

    rhs(ieq) = rhs(ieq) + creep_term_ip1

    icr = icr + 1
    creep_term(icr) = creep_term_ip1

  endif

  ! second equation
  ! pressure on the inside surface of the last oxide layer (on the steam side)
  ieq = ieq + 1
  coeff(ieq, ieps0) = poisson(N) / (1.0 - poisson_star(N))
  ! coeff(ieq, ieps0) = 0.0

  coeff(ieq, ib(N)) = 1.0 / (1.0 - poisson_star(N))
  coeff(ieq, ic(N)) = - 1.0 / (r2(N+1) * (1.0 + poisson_star(N)))

  rhs(ieq) = -p_in / Youngs_star(N) + STRAIN_OX_IND_TERM_RADIAL( &
        & poisson(N), poisson_star(N), r(N+1), &
        & ox_strain_r(N), ox_strain_hoop(N), ox_strain_axial(N))
  
  ! lucky on the internal surface, that of the oxide surface
  ! creep integrals are ZERO on this first location (the lower integration limit)

  ! second type of equations
  ! continuation of radial stresses
  ! s_r = E^* [ b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2 + &
  !             eps0 * nu/(1-nu*) - Fr(strain_ox, r) - &
  !             Jcr_ri,r(creep_sum*r) / (2*r2) + Jcr_ri,r(creep_dif/r) / 2]

  do i = 1, N -1

    ! condition for interface between materials 'i' and 'i+1'

    ieq = ieq + 1

    Youngs_ratio = Youngs_star(i) / Youngs_star(i+1)

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
    rhs(ieq) = Integral(i+1) / r2(i+1)
      
    if (i+1 >= first_oxide_layer)  then
      ox_term_ip1 = STRAIN_OX_IND_TERM_RADIAL(poisson(i+1),  &
        & poisson_star(i+1), r(i+1), &
        & ox_strain_r(i+1), ox_strain_hoop(i+1), ox_strain_axial(i+1))

    else
      ! contribution from 'i+1' layer
      ox_term_ip1 = 0.0
    endif

    if (i >= first_oxide_layer)  then
      ! contribution from 'i' layer: metal does not have any oxidation strains

      ox_term_i = STRAIN_OX_IND_TERM_RADIAL(poisson(i), poisson_star(i), &
        & r(i+1), ox_strain_r(i), ox_strain_hoop(i), ox_strain_axial(i))
    else
      ox_term_i = 0.0
    endif

    rhs(ieq) = rhs(ieq) + ox_term_ip1 - Youngs_ratio * ox_term_i

    ! add the creep terms only when needed
    if (i + 1 >= i_cr_start .and. i + 1 <= i_cr_end)  then

      ! etr_sum == summ of (strain * r)
      ! eor_dif == (e_r - e_hoop) / r
      creep_term_ip1 = CREEP_TERM_S_RAD(creep_integral_sum(i +1, 1), &
          & creep_integral_dif(i +1, 1), r2(i+1))

      rhs(ieq) = rhs(ieq) + creep_term_ip1

      icr = icr + 1
      creep_term(icr) = creep_term_ip1

    endif

  end do

  ! third type of equations
  ! continuation of radial displacement
  ! u = r * b + c/r + (1+nu*) * J_ri,r/r + ox_ind_term_u +
  !     (1+nu*) * Jcr_ri,r(creep_sum*r) / (2*r) + 
  !     (1-nu*) * r * Jcr_ri,r(creep_dif/r) / 2

  do i = 1, N -1

    ieq = ieq + 1

    ! no coefficients for eps0

    ! coeff for b
    coeff(ieq, ib(i)) = - r(i+1)
    coeff(ieq, ib(i+1)) = r(i+1)

    ! coeff for c
    coeff(ieq, ic(i)) = - 1.0 / r(i+1)
    coeff(ieq, ic(i+1)) = 1.0 / r(i+1)

    ! rhs
    rhs(ieq) = - Integral(i+1) * (1.0 + poisson_star(i+1)) / r(i+1)
    
    ! previously envisioned for Oh (2006) method
    ! if (i == 2 .and. if_growth_induce_stress)  then
    !  rhs(ieq) = rhs(ieq) - u_now_int_jump(id_high_or_low)
    ! endif

    if (i + 1 >= first_oxide_layer)  then

      ox_term_ip1 = U_TERM_STRAIN_OX_IND(poisson_star(i+1), r(i+1), &
        & ox_strain_r(i+1), ox_strain_hoop(i+1))

      ! at metal-oxide interface metal does not have any oxidation strains
    else
      ox_term_ip1 = 0.0

    endif

    if (i >= first_oxide_layer)  then

      ox_term_i = U_TERM_STRAIN_OX_IND(poisson_star(i), &
        & r(i+1), ox_strain_r(i), ox_strain_hoop(i))

    else
      ox_term_i = 0.0

    endif

    rhs(ieq) = rhs(ieq) + ox_term_i - ox_term_ip1

      ! check the sign for u_now_int_jump(id_high_or_low) 
      ! + sign for u_inner - u_outer = U_delta
      ! rhs(ieq) = rhs(ieq) + ox_term_i + u_now_int_jump(id_high_or_low)
      ! u_int_jump_now(1, i) - interface between i and 'i-1' layers
      rhs(ieq) = rhs(ieq) + u_int_jump_now(1, i+1)

      if (ABS(u_int_jump_now(1, i+1)) > 0.0)  then
 
        write(aux_lun, *) 'u_jump ', u_int_jump_now(1, i+1)

      endif

    ! add the creep terms only when needed
    if (i + 1 >= i_cr_start .and. i + 1 <= i_cr_end)  then

      ! etr_sum == summ of (strain * r)
      ! eor_dif == (e_r - e_hoop) / r
      creep_term_ip1 = CREEP_TERM_U(creep_integral_sum(i +1, 1), &
          & creep_integral_dif(i +1, 1), poisson_star(i+1), r(i+1))

      rhs(ieq) = rhs(ieq) + creep_term_ip1

      icr = icr + 1
      creep_term(icr) = creep_term_ip1

    endif

  end do

  ! the last equation; axial balance or set eps0 to zero
  ieq = ieq + 1

  if (generalized_plane_strain)  then

    ! write(tty_lun, *) 'FIX IT according to new I_i (computed from the inner radius' 
    ! stop

    ! initialize
    coeff(ieq, ieps0) = 0.0
    
    ! rhs(ieq) = 66.26*(r2(1)-r2(N+1))/2.0 for validation_zarrabi1.txt
    rhs(ieq) = 0.0

    do i = 1, N

      ! radius factor 
      rf = 0.5 * (r2(i) - r2(i+1))

      ! poisson factor
      pf_eps = (1.0 - poisson(i)) / ((1.0 + poisson(i)) * &
             & (1.0 - 2.0 * poisson(i)))
      pf_b = 2.0 * poisson(i) / (1.0 - poisson_star(i))

      coeff(ieq, ieps0) = coeff(ieq, ieps0) + rf * pf_eps * Youngs(i)

      coeff(ieq, ib(i)) = rf * pf_b * Youngs_star(i)

      rhs(ieq) = rhs(ieq) + (Youngs(i) / (1.0 + poisson(i)) + &
               & Youngs_star(i) * poisson(i)) * &
               & Integral(i)

      if (i >= first_oxide_layer)  then
        ! only for the oxide layers

        ! A - term
        rhs(ieq) = rhs(ieq) + rf * STRESS_AXIAL_OX_IND_TERM( &
            & poisson(i), poisson_star(i), Youngs(i), Youngs_star(i), &
            & ox_strain_r(i), ox_strain_hoop(i), ox_strain_axial(i))

        ! B - term
        r_lnf = 0.5 * (r2(i) * LOG(ABS(r(i))) - r2(i+1) * LOG(ABS(r(i+1))))

        rhs(ieq) = rhs(ieq) - Youngs_star(i) * poisson(i) * &
            & (ox_strain_r(i) - ox_strain_hoop(i)) * &
            & (r_lnf - 0.5 * rf)

      endif

      ! add the creep terms only when needed
      if (i >= i_cr_start .and. i <= i_cr_end)  then

        creep_term_ip1 = CREEP_TERM_INT_AXIAL(Youngs(i), poisson(i), &
          & creep_integral_sum(i, 1), creep_integral_dif(i, 1), &
          & creep_integral_z(i), r2(i))

        rhs(ieq) = rhs(ieq) + creep_term_ip1

        icr = icr + 1
        creep_term(icr) = creep_term_ip1

      endif

    end do

  else if (.not. generalized_plane_strain)  then
    
    coeff(ieq, ieps0) = 1.0
    rhs(ieq) = 0.0

  endif

  ! scale rhs
  rhs = scale_solver * rhs

  if (if_debug_print)  then

    do i = 1, 2*N+1
      write(aux_lun, 2)  cycle_no, i, (coeff(i, j), j = 1, 2*N+1), rhs(i)
 2  format('sfs ', i5, 1x, i2, 1x, 30(1pe13.6, 1x))
    end do

    if (icr == 0)  then
      write(aux_lun, *) 'no creep terms ', cycle_no
    else

      write(aux_lun, 7)  cycle_no, (creep_term(i), i = 1, icr)
 7  format('creep_term_st ', i5, 1x, 30(1pe13.6, 1x))

    endif

  endif

  ! determine the stress constants

  call SOLVER_SIMPLE(2*N + 1, coeff, rhs)

  ! scale down the solution
  rhs = rhs / scale_solver

  ! store the results
  eps0 = rhs(1)
  do i = 1, N
    bc_st(i) = rhs(ib(i))
    cc_st(i) = rhs(ic(i))
  end do

  if (if_debug_print)  then

    write(aux_lun, 5) cycle_no, eps0, (bc_st(i), i= 1, N)
 5 format('const_sol_e_b=', i5, 1x, 30(1pe13.6, 1x))
    write(aux_lun, 6) cycle_no, eps0, (cc_st(i), i= 1, N)
 6 format('const_sol_e_c=', i5, 1x, 30(1pe13.6, 1x))

  endif

  call GET_STRESS_STRAIN(N, cycle_no, poisson, th_exp, Youngs, &
       & poisson_star, th_exp_star, Youngs_star)

  ! get displacements and dimensions
  call GET_OG_DIM_DISPL(N, poisson, th_exp, Youngs, &
       & poisson_star, th_exp_star, Youngs_star, if_update_displ_dim)

  return

  END SUBROUTINE STRESS_PROFILE

  SUBROUTINE GET_STRESS_STRAIN(N, cycle_no, poisson, th_exp, Youngs, &
       & poisson_star, th_exp_star, Youngs_star)

  use solution_data_module, only:  npr_st, rad_st,  &
      & rad_int, integral_th, tau_th, tau_th_st, eps0, &
      & strain_rad, strain_hoop, displ_rad, &
      & stress_axial, stress_hoop, stress_rad, strain_hoop_gen, &
      & strain_hoop_el, strain_rad_el, strain_axial_el, bc_st, cc_st, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low 
  use oxide_data_module,   only: Youngs_modul, Youngs_temp, &
      & th_exp_coeff, first_oxide_layer
  use oxide_stress_data_module,   only: trace_total, trace_without_th, &
      & ox_strain_r, ox_strain_hoop, &
      & ox_strain_axial
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use induced_oxide_stress_module, only: STRAIN_OX_IND_TERM_RADIAL, &
       & STRAIN_OX_IND_TERM_HOOP, &
       & U_TERM_STRAIN_OX_IND, DUDR_TERM_STRAIN_OX_IND, &
       & STRESS_AXIAL_OX_IND_TERM
  use creep_data_module, only: creep_strain_new, &
       & creep_integral_sum, creep_integral_dif, &
       & if_creep_ox_scale, if_creep_metal, i_cr_start, i_cr_end
  use creep_module, only: CREEP_TERM_S_RAD, CREEP_TERM_U, &
       & CREEP_TERM_INT_AXIAL, CREEP_TERM_S_HOOP, CREEP_TERM_DUDR
  use stress_utility_module, only: RADIAL_STRESS, HOOP_STRESS, &
       & AXIAL_STRESS, RADIAL_DISPL, HOOP_STRAIN, RADIAL_STRAIN

  ! Argument List
  integer,    intent(IN)          :: N, cycle_no
  real, dimension(N), intent(IN)   :: poisson, th_exp, Youngs, &
       & poisson_star, th_exp_star, Youngs_star

  ! Local Variables
  integer  :: i, j
  real     :: Y_alfa, pf_c, pf_b, df, r, r2, strain_hoop_eq, &
            & Integral, tau, ox_growth_induce_term, creep_term
  ! current radius
  integer, save  :: ip = 0
  logical  :: if_debug_print

  ! using fractional integrals
  ! u = r * b + c/r + (1+nu*) * J_ri,r/r
  ! s_r = E^* (b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2)
  ! with eps0 becomes:
  ! s_r = E^* [ b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2 + eps0 * nu/(1-nu*)]
  ! s_teta = E^* (b/(1-nu*) - tau + c/(r^2 * (1+ nu*)) + J_ri,r/r^2)
  ! with eps0 becomes:
  ! s_teta = E^* [b/(1-nu*) - tau + c/(r^2 * (1+ nu*)) + J_ri,r/r^2 + eps0 * nu/(1-nu*)]

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  if_debug_print = .false.
  ! check if print out at different oxide thicknesses
  do i = 1, no_full2low_output

    ! print hoop stress strain profile at different oxide thicknesses
    if (cycle_no == no_point_full2low(i) .or. &
        & cycle_no - 1 == no_point_full2low(i) )  then
      if_debug_print = .true.
    endif

  end do

  do i = 1, N

    do j = 1, npr_st(i)

      r = rad_st(i, j) 
      r2 = r**2

      ! huge simplification brought by including the 
      ! thermal expansion into the expression for tau and integral
      ! this is to account for the use of alfa* not alfa to 
      ! integrate the thermal expansion, Eq. 6.2 in Nada
      Integral = (1.0 + poisson(i)) * Integral_th(i, j)
      tau = (1.0 + poisson(i)) * tau_th_st(i, j)

      ! get the stresses
      ! second part is the due to oxide-growth induced strains
      ! s_r = E^* [ b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2 + &
      !         eps0 * nu/(1-nu*) - Fr(strain_ox, r)- &
      !         Jcr_ri,r(creep_sum*r) / (2*r2) + Jcr_ri,r(creep_dif/r) / 2]

      if (i >= first_oxide_layer)  then

        ox_growth_induce_term = STRAIN_OX_IND_TERM_RADIAL(poisson(i), &
         & poisson_star(i), r, ox_strain_r(i), &
         & ox_strain_hoop(i), ox_strain_axial(i))
      else

        ox_growth_induce_term = 0.0

      endif

      ! add the creep terms only when needed
      if (i >= i_cr_start .and. i <= i_cr_end)  then
 
        creep_term = CREEP_TERM_S_RAD(creep_integral_sum(i, j), &
          & creep_integral_dif(i, j), r2)

      else
        ! creep term should be zero if no creep is present
        creep_term = 0.0
      endif

      stress_rad(i, j) = RADIAL_STRESS(Youngs_star(i), Integral, &
         & r2, bc_st(i) / (1.0 - poisson_star(i)), &
         & cc_st(i) / (1.0 + poisson_star(i)), &
         & eps0 * poisson(i) / (1.0 - poisson_star(i))) - &
         & Youngs_star(i) * (ox_growth_induce_term + creep_term)

      ! s_teta / E^* = b/(1-nu*) - tau + c/(r^2 * (1+ nu*)) + J_ri,r/r^2 + 
      !                eps0 * nu/(1-nu*) - Fhoop

      if (i >= first_oxide_layer)  then

        ox_growth_induce_term = STRAIN_OX_IND_TERM_HOOP(poisson(i), &
         & poisson_star(i), r, ox_strain_r(i), ox_strain_hoop(i), &
         & ox_strain_axial(i))

      else

        ox_growth_induce_term = 0.0

      endif

      ! add the creep terms only when needed
      if (i >= i_cr_start .and. i <= i_cr_end)  then
        creep_term = CREEP_TERM_S_HOOP(creep_integral_sum(i, j), &
          & creep_integral_dif(i, j), r2) - &
          & (creep_strain_new(i, j)%hoop + &
          & poisson(i) * creep_strain_new(i, j)%ax)

      else
        ! creep term should be zero if no creep is present
        creep_term = 0.0
      endif

      stress_hoop(i, j) = HOOP_STRESS(Youngs_star(i), Integral, &
         & r2, bc_st(i) / (1.0 - poisson_star(i)), &
         & cc_st(i) / (1.0 + poisson_star(i)), &
         & eps0 * poisson(i) / (1.0 - poisson_star(i)), tau) - &
         & Youngs_star(i) * (ox_growth_induce_term - creep_term)

      ! hoop_strain from stress-strain equations; 
      ! last terms due to oxide growth
      strain_hoop_eq = tau - eps0 * poisson(i) + &
          & (stress_hoop(i, j) - poisson_star(i) * stress_rad(i, j)) / &
          & Youngs_star(i) + &
          & ox_strain_hoop(i) + poisson(i) * ox_strain_axial(i) + &
          & creep_strain_new(i, j)%hoop + &
          & poisson(i) * creep_strain_new(i, j)%ax

      ! alpha * tau = tau_th in the code
      ! used to be Youngs(i) * (eps0 - tau) it was wrong
      ! last terms due to oxide growth
      stress_axial(i, j) = Youngs(i) * (eps0 - tau_th_st(i, j)) + &
         & poisson(i) * (stress_rad(i, j) + stress_hoop(i, j)) - &
         & Youngs(i) * (ox_strain_axial(i) + creep_strain_new(i, j)%ax)

      ! get the strains
      if (i >= first_oxide_layer)  then

        ! last terms due to oxide growth, 
        ! e_hoop = u/r = b + c/r**2 + (1+nu*) * J_ri,r / r**2 + u_ox_term/r
        ox_growth_induce_term = U_TERM_STRAIN_OX_IND(poisson_star(i), r, &
          & ox_strain_r(i), ox_strain_hoop(i))

      else

        ox_growth_induce_term = 0.0

      endif

      ! add the creep terms only when needed
      if (i >= i_cr_start .and. i <= i_cr_end)  then
 
        creep_term = -CREEP_TERM_U(creep_integral_sum(i, j), &
          & creep_integral_dif(i, j), poisson_star(i), r) 

      else
        ! creep term should be zero if no creep is present
        creep_term = 0.0
      endif

      strain_hoop(i, j) = HOOP_STRAIN(Integral, r, &
          & bc_st(i), cc_st(i), 1.0 + poisson_star(i)) + &
          & (ox_growth_induce_term + creep_term) / r

      ! elastic strain to be used for obtaining elastic storage energy
      strain_hoop_el(i, j) = strain_hoop(i, j) - tau_th_st(i, j) - &
        & ox_strain_hoop(i) - creep_strain_new(i, j)%hoop

      if (if_debug_print)  then
      ! if (ip <= 6)  then

        write(aux_lun, 2) i, j, strain_hoop_eq, strain_hoop(i, j), &
          & abs(strain_hoop_eq - strain_hoop(i, j))
 2      format('eq_strain_hoop ', i1, 1x, i2, 1x, 20(1pe13.6, 1x))

      endif

      if (i >= first_oxide_layer)  then

        ! last terms due to oxide growth, 
        ! du/dr = b + (1+nu*)*tau - c/r2 - (1+nu*) * J_ri,r/r2 + 
        !         d(u_ox_term)/dr
        ox_growth_induce_term = DUDR_TERM_STRAIN_OX_IND(poisson_star(i), r, &
          & ox_strain_r(i), ox_strain_hoop(i))

      else

        ox_growth_induce_term = 0.0

      endif

      ! add the creep terms only when needed
      if (i >= i_cr_start .and. i <= i_cr_end)  then
 
        creep_term = CREEP_TERM_DUDR(creep_integral_sum(i, j), &
          & creep_integral_dif(i, j), poisson_star(i), r) + &
          & creep_strain_new(i, j)%rad + &
          & poisson_star(i) * creep_strain_new(i, j)%hoop + &
          & poisson(i) * (1.0 + poisson_star(i)) * creep_strain_new(i, j)%ax

      else
        ! creep term should be zero if no creep is present
        creep_term = 0.0
      endif

      strain_rad(i, j) = RADIAL_STRAIN(Integral, &
          & tau, r2, bc_st(i), cc_st(i), 1.0 + poisson_star(i)) + &
          & ox_growth_induce_term + creep_term

      ! elastic strain to be used for obtaining elastic storage energy
      strain_rad_el(i, j) = strain_rad(i, j) - tau_th_st(i, j) - &
        & ox_strain_r(i) - creep_strain_new(i, j)%rad

      ! elastic strain to be used for obtaining elastic storage energy
      strain_axial_el(i, j) = eps0 - tau_th_st(i, j) - &
        & ox_strain_axial(i) - creep_strain_new(i, j)%ax

      ! total trace = e_r + e_hoop + e_z
      trace_total(i, j) = strain_hoop(i, j) + strain_rad(i, j) + eps0
      
      ! exclude stress-free thermal strains
      ! for hoop and radial strain: e_th = tau_th_st(i, :) * (1.0 + poisson(i)) - poisson(i) * eps0
      ! for axial: e_th = should be 0 or eps0?? do not know this one now
      ! 
      trace_without_th(i, j) = trace_total(i, j) - &
          & 2.0 * tau_th_st(i, j) * (1.0 + poisson(i)) + &
          & 2.0 * poisson(i) * eps0

    end do

    ! strain_hoop_gen(i, :) = strain_hoop(i, :) - &
    !     & tau_th_st(i, :) * (1.0 + poisson(i)) + poisson(i) * eps0
    strain_hoop_gen(i, :) = (stress_hoop(i, :) - &
         & poisson_star(i) * stress_rad(i, :)) / Youngs_star(i)

    ! write(aux_lun, 1) i, Y_alfa, Integral_th(i, j), &
    !     & bc_st(i), cc_st(i), df, (stress_rad(i, j), j = 1, npr_st(i))
 1   format('stress_rad1 ', i2, 1x, 20(1pe13.6, 1x))

  end do

  if (if_debug_print)  then
  ! if (ip <= 6)  then

    write(aux_lun, 3) (stress_hoop(i, npr_st(i)), stress_hoop(i, 1), i=2, N), &
        & stress_hoop(1, npr_st(1)), stress_hoop(1, 1)
 3   format('stress_hoop oxide, ox-metal ', 20(1pe13.6, 1x))

  endif

  return

  END SUBROUTINE GET_STRESS_STRAIN

END MODULE STRESS_SOLUTION_MODULE
