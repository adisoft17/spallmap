MODULE INDUCED_STRESS_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for temperature, stress, oxided-induced stresses. 
    !  it will determine the temperature profile in the substrate and coating
    !  in the steady state and during ramp up (down)
    !
    ! sabaua@ornl.gov induced_stress_module.F90
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================


  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: OXIDE_INDUCED_UPDATE, OXIDE_INDUCED_SETUP, SET_OX_STRAINS, &
      & GET_U_JUMP_SURFACE

 CONTAINS

 !  SUBROUTINE STRESS_CYL_OX_INDUCED(N, p_out, p_in, poisson, Youngs, &
 !      & id_high_or_low, if_ave_scale)
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the coefficients for the stress-strain and dimensions 
    !     using the oxide-induced growth stresses
    !     valid only for three layers (metal, spinel, magnetite)
    !   based on Oh et al. (2006)
    !=======================================================================


  ! return

  ! END SUBROUTINE STRESS_CYL_OX_INDUCED

  SUBROUTINE OXIDE_INDUCED_SETUP(iter_og, id_high_or_low)

    !=======================================================================
    ! Purpose(s):
    !
    !   setup coefficiencs for oxide induced growth constraints
    !   based on Oh et al. (2006)
    !          
    !  
    !=======================================================================

  use boiler_data_module, only: tube_thickness, tube_outer_radius
  use solution_data_module, only:  npr_st, rad_st,  &
      & rad_int, integral_th, tau_th, tau_th_st, eps0, &
      & strain_rad, strain_hoop, displ_rad, &
      & stress_axial, stress_hoop, stress_rad, bc_st, cc_st
  use oxide_data_module,   only: Youngs_modul, Youngs_temp, &
      & th_exp_coeff, pilling_bedworth_ratio, poisson_ratio
  use output_module,             only: tty_lun, out_lun, inp_lun, aux_lun
  use solver_data_module, only       : neq_solver, fval, jacf, var, &
      & delta_var
  use oxide_stress_data_module,  only: imetal, ispinel, imagnetite, &
      & iout, imet_spin, ispin_mag, i_in
  use oxide_stress_data_module,  only: rc_og, u_og, rr, trace_strain, &
      & h_ox_sp, h_ox_mag, h_ox_total, ir_c_i, ir_c_m_out, ir_r_i_s, &
      & ir_r_o_ma, pbr_sp, pbr_mag, &
      & rr_metal_in_initial, area_fr_metal_2_spinel

  ! Argument List
  integer, intent(IN)  :: iter_og
  integer, intent(IN)  :: id_high_or_low

  ! Local Variables
  integer  :: i, j, neq, nvar
  integer  :: ir_r_i_ma, ir_c_ms, ir_r_i_m, ir_r_o_s, ir_c_sm, &
       & ir_r_m_out
  real     :: Y_alfa, pf_c, pf_b, df, r, r2, strain_hoop_eq, &
            & Integral, tau ! current radius
  integer, save  :: ip = 0

  if (iter_og == 1)  then

    ! initialize the initial outer radius of the tube; reference radius
    rr_metal_in_initial = (tube_outer_radius - tube_thickness) * &
      & (1.0 + tau_th_st(imetal, 1)*(1.0 + poisson_ratio(imetal)))

  endif

    ! setup jacobian for imposing constraints due to oxide growth
    neq = 0
    nvar = 0
    fval = 0.0
    jacf = 0.0

    ! ir_c - stores id for radius in the current system (deformed)
    ! ir_r - stores id for radius in the reference system (stress free)

    !  magnetite surface displacement
    neq = neq + 1
    nvar = nvar + 1
    ir_c_i = nvar  ! index for r_i
    nvar = nvar + 1
    ir_r_i_ma = nvar
    ! r_i  - first variable
    ! R_i,ma = r_i - u_i, ma

    if (iter_og == 1) then
      var(ir_c_i) = rc_og(i_in)
      var(ir_r_i_ma) = var(ir_c_i) - u_og(i_in, 1)
    endif

    fval(neq) = var(ir_r_i_ma) - var(ir_c_i) + u_og(i_in, 1)  
    jacf(neq, ir_r_i_ma) = 1.0
    jacf(neq, ir_c_i) = -1.0
    
    ! current radius at metal-spinel interface
    neq = neq + 1
    nvar = nvar + 1
    ir_c_ms = nvar

    if (iter_og == 1) then
      var(ir_c_ms) = var(ir_c_i) + h_ox_total
    endif

    fval(neq) = var(ir_c_ms) - var(ir_c_i) - h_ox_total   ! r_ms = r_i + h_s + h_ma
    jacf(neq, ir_c_ms) = 1.0
    jacf(neq, ir_c_i) = -1.0

    ! reference radius of the metal
    neq = neq + 1
    nvar = nvar + 1
    ir_r_i_m = nvar
    ir_r_o_s = ir_r_i_m         ! R_o,s = R_i,m
    ! R_i,m = r_i + h_s + h_ma - u_ms

    if (iter_og == 1) then
      var(ir_r_i_m) = var(ir_c_i) + h_ox_total - u_og(imet_spin, 1) 
    endif

    fval(neq) = var(ir_r_i_m) - var(ir_c_i) - h_ox_total + u_og(imet_spin, 1) 
    jacf(neq, ir_r_i_m) = 1.0
    jacf(neq, ir_c_i) = -1.0

    ! current radius at spinel-magnetite interface
    neq = neq + 1
    nvar = nvar + 1
    ir_c_sm = nvar
    ! r_sm = r_i + h_ma

    if (iter_og == 1) then
      var(ir_c_sm) = var(ir_c_i) + h_ox_mag
    endif

    fval(neq) = var(ir_c_sm) - var(ir_c_i) - h_ox_mag   
    jacf(neq, ir_c_sm) = 1.0
    jacf(neq, ir_c_i) = -1.0

    ! current radius at the outer surface of metal
    neq = neq + 1
    nvar = nvar + 1
    ir_c_m_out = nvar  ! index for r_o
    nvar = nvar + 1
    ir_r_m_out = nvar

    if (iter_og == 1) then
      var(ir_c_m_out) = rc_og(iout)
      var(ir_r_m_out) = var(ir_c_m_out) - u_og(iout, 1)  
    endif
    ! 
    fval(neq) = var(ir_r_m_out) - var(ir_c_m_out) + u_og(iout, 1)  
    jacf(neq, ir_r_m_out) = 1.0
    jacf(neq, ir_c_m_out) = -1.0

    ! fm equation density ratio - dilation relationship for metal
    neq = neq + 1
    ! no new variables
    fval(neq) = var(ir_r_m_out)**2 - var(ir_r_i_m)**2 - &
      & (1.0 - trace_strain(imetal)) * (var(ir_c_m_out)**2 - var(ir_c_ms)**2)
    jacf(neq, ir_r_m_out) = 2.0 * var(ir_r_m_out)
    jacf(neq, ir_r_i_m) = - 2.0 * var(ir_r_i_m)
    jacf(neq, ir_c_m_out) = - 2.0 * (1.0 - trace_strain(imetal)) * var(ir_c_m_out)
    jacf(neq, ir_c_ms) = 2.0 * (1.0 - trace_strain(imetal)) * var(ir_c_ms)

    ! fs equation; density ratio - dilation relationship for spinel
    neq = neq + 1
    nvar = nvar + 1
    ir_r_i_s = nvar  ! index for R_i,s

    if (iter_og == 1) then
      var(ir_r_i_s) = rr(ispin_mag, 1)
    endif

    fval(neq) = var(ir_r_o_s)**2 - var(ir_r_i_s)**2 - &
        & (1.0 - trace_strain(ispinel)) * (var(ir_c_ms)**2 - var(ir_c_sm)**2)
    jacf(neq, ir_r_o_s) = 2.0 * var(ir_r_o_s)
    jacf(neq, ir_r_i_s) = - 2.0 * var(ir_r_i_s)
    jacf(neq, ir_c_ms) = - 2.0 * (1.0 - trace_strain(ispinel)) * var(ir_c_ms)
    jacf(neq, ir_c_sm) = 2.0 * (1.0 - trace_strain(ispinel)) * var(ir_c_sm)

    ! f_mag equation; density ratio - dilation relationship for magnetite
    neq = neq + 1
    nvar = nvar + 1
    ir_r_o_ma = nvar  ! index for R_o,ma

    if (iter_og == 1) then
      var(ir_r_o_ma) = rr(ispin_mag, 2)
    endif

    fval(neq) = var(ir_r_o_ma)**2 - var(ir_r_i_ma)**2 - &
        & (1.0 - trace_strain(imagnetite)) * (var(ir_c_sm)**2 - var(ir_c_i)**2)
    jacf(neq, ir_r_o_ma) = 2.0 * var(ir_r_o_ma)
    jacf(neq, ir_r_i_ma) = - 2.0 * var(ir_r_i_ma)
    jacf(neq, ir_c_sm) = - 2.0 * (1.0 - trace_strain(imagnetite)) * var(ir_c_sm)
    jacf(neq, ir_c_i) = 2.0 * (1.0 - trace_strain(imagnetite)) * var(ir_c_i)

    ! pbr is the same for spinel and magnetite
    pbr_sp = pilling_bedworth_ratio
    pbr_mag = pilling_bedworth_ratio

    ! pbr equation
    neq = neq + 1
    ! no new variables
    fval(neq) = pbr_mag * (var(ir_r_o_s)**2 - var(ir_r_i_s)**2) + &
           & pbr_sp * (var(ir_r_o_ma)**2 - var(ir_r_i_ma)**2) - &
           & pbr_mag * pbr_sp * (var(ir_r_o_s)**2 - rr_metal_in_initial**2)
    jacf(neq, ir_r_o_s) = 2.0 * var(ir_r_o_s) * pbr_mag * (1.0 - pbr_sp)
    jacf(neq, ir_r_i_s) = - 2.0 * pbr_mag * var(ir_r_i_s)
    jacf(neq, ir_r_o_ma) = 2.0 * pbr_sp * var(ir_r_o_ma)
    jacf(neq, ir_r_i_ma) = - 2.0 * pbr_sp * var(ir_r_i_ma)

    neq_solver = neq

    if (.not. neq_solver == nvar)  then

      write(tty_lun, *) 'ERROR: oxide-induced-stress nvar ne neq'
      stop

    endif

    do i = 1, neq_solver
      write(aux_lun, 2) i, (jacf(i, j), j= 1, neq_solver), fval(i)
 2    format('oxis_jac ', i2, 1x, 20(1pe13.6, ' '))
    end do

    ! store the fractional area of the metal that converts to spinel
    ! f_sg
    area_fr_metal_2_spinel = (var(ir_r_o_s)**2 - var(ir_r_i_s)**2) / &
        & (pbr_sp * (var(ir_r_o_s)**2 - rr_metal_in_initial**2))

    write(aux_lun, 3) iter_og, area_fr_metal_2_spinel
 3    format('a_metal2spinel ', i2, 1x, 20(1pe13.6, ' '))

  return

  END SUBROUTINE OXIDE_INDUCED_SETUP

  SUBROUTINE SET_OX_STRAINS(N, cycle_no)

    !=======================================================================
    ! Purpose(s):
    !
    !   set oxidation strains in the oxide
    !  
    !=======================================================================

  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & ox_strain_r, ox_strain_hoop, ox_strain_hoop_old, &
      & ox_strain_axial, total_oxide_thickness_ref, rr_metal_in_initial, &
      & d_oxide_thickness_4_strain_hoop
  use oxide_stress_data_module,  only: oxidation_strain_mode, &
      & growth_ox_location
  use boiler_data_module, only: tube_outer_radius, tube_thickness, &
      & Time_Boiler_p, oxide_thickness
  use oxide_data_module,   only: pilling_bedworth_ratio, &
      & pbr_bernstein_factor, thickness_fr, lateral_ox_strain_length, &
      & ox_hoop_strain_factor, ox_axial_strain_factor, &
      & first_oxide_layer
  use parameter_module,   only: moxide

  ! Argument List
  integer,    intent(IN)  :: N, cycle_no
  ! cycle_no == i from i= 2, TEMP_BOILER_data_no ! no of data intervals

  ! Local Variables
  integer  :: i, k, j
  real, dimension(N+1)      :: psi
  real     :: b_term, c_term, r_metal_oxide, u_jump_fr, &
           & total_oxide_thickness_ref_si, delta_h, dstrain_ox
  logical, save  :: first_set_ox_strains = .true.

  ox_strain_r = 0.0
  ox_strain_hoop = 0.0
  ox_strain_axial = 0.0

  if (.not. if_growth_induce_stress)  RETURN

  ! delta_h = oxide_thickness(cycle_no) - oxide_thickness(cycle_no-1)
  delta_h = d_oxide_thickness_4_strain_hoop

  write(out_lun, 10) cycle_no, oxide_thickness(cycle_no), &
      & oxide_thickness(cycle_no) - oxide_thickness(cycle_no-1), &
      & d_oxide_thickness_4_strain_hoop, first_set_ox_strains
  10 format('check_set_ox_str ', i6, 1x, 3(1pe13.6, 1x), 1(L2, 1x))

  if (first_set_ox_strains)  then

    ox_strain_hoop_old = 0.0
    first_set_ox_strains = .false.

  endif 

  ! set the oxidation strains
  do i = first_oxide_layer, N

    if (TRIM(oxidation_strain_mode(i)) == 'radial')  then

      ox_strain_r(i) = pbr_bernstein_factor * &
         & (pilling_bedworth_ratio - 1.0)

    else if (TRIM(oxidation_strain_mode(i)) == 'dilational')  then

      ox_strain_r(i) = pbr_bernstein_factor * &
         & (pilling_bedworth_ratio**(1.0/3.0) - 1.0)
      ox_strain_hoop(i) = ox_strain_r(i)
      ox_strain_axial(i) = ox_strain_r(i)

    else if (TRIM(oxidation_strain_mode(i)) == 'radial_dominant')  then

      ox_strain_r(i) = pbr_bernstein_factor * &
         & (pilling_bedworth_ratio - 1.0)
      ox_strain_hoop(i) = ox_hoop_strain_factor * ox_strain_r(i)
      ox_strain_axial(i) = ox_axial_strain_factor * ox_strain_r(i)

    else if (TRIM(oxidation_strain_mode(i)) == 'radial_lateral_hoop')  then

      ! delta_h = oxide_thickness(cycle_no) - oxide_thickness(cycle_no-1)
      ! bug; for transition cycle_no is not the absolute cycle number
      dstrain_ox = delta_h / lateral_ox_strain_length

      ! lateral_ox_strain_length [microns]
      ox_strain_hoop(i) = ox_strain_hoop_old(i) + dstrain_ox
      ! store previous value for the next cycle
      ox_strain_hoop_old(i) = ox_strain_hoop(i)

      ox_strain_axial(i) = 0.0

      ox_strain_r(i) = pbr_bernstein_factor * &
         & (pilling_bedworth_ratio - 1.0) - ox_strain_hoop(i) - &
         & ox_strain_axial(i)

      if (ox_strain_r(i) < 1.0e-9)  then
        write(6, *) 'ERROR: oxidation_strain_mode(i) == radial_lateral_hoop'
        write(6, *) 'ox_strain_r < 1.0e-9 '
        STOP
      endif

    else if (TRIM(oxidation_strain_mode(i)) == 'radial_lateral_hoop_axial')  then

      ! delta_h = oxide_thickness(cycle_no) - oxide_thickness(cycle_no-1)
      ! bug; for transition cycle_no is not the absolute cycle number
      dstrain_ox = delta_h / lateral_ox_strain_length

      ! lateral_ox_strain_length [microns]
      ox_strain_hoop(i) = ox_strain_hoop_old(i) + dstrain_ox
      ! store previous value for the next cycle
      ox_strain_hoop_old(i) = ox_strain_hoop(i)

      ox_strain_axial(i) = ox_strain_hoop(i)

      ox_strain_r(i) = pbr_bernstein_factor * &
         & (pilling_bedworth_ratio - 1.0) - ox_strain_hoop(i) - &
         & ox_strain_axial(i)

      if (ox_strain_r(i) < 1.0e-9)  then
        write(6, *) 'ERROR: oxidation_strain_mode(i) == radial_lateral'
        write(6, *) 'ox_strain_r < 1.0e-9 '
        STOP
      endif

    else if (TRIM(oxidation_strain_mode(i)) == 'lateral_hoop')  then

      ! delta_t = Time_Boiler_p(cycle_no) - Time_Boiler_p(cycle_no-1)
      ! h_ox = oxide_thickness(cycle_no) * oxide_thick2si
      ! delta_h2 = oxide_thickness(cycle_no)**2 - oxide_thickness(cycle_no-1)**2
      ! dh_ox = delta_h2 / (2.0 * oxide_thickness(cycle_no))
      !
      ! delta_h = oxide_thickness(cycle_no) - oxide_thickness(cycle_no-1)
      ! bug; for transition cycle_no is not the absolute cycle number
      dstrain_ox = delta_h / lateral_ox_strain_length

      ! lateral_ox_strain_length [microns]
      ox_strain_hoop(i) = ox_strain_hoop_old(i) + dstrain_ox
      ox_strain_axial(i) = 0.0
      ox_strain_r(i) = 0.0

      ! store previous value for the next cycle
      ox_strain_hoop_old(i) = ox_strain_hoop(i)

    else if (TRIM(oxidation_strain_mode(i)) == 'lateral_hoop_axial')  then

      ! delta_t = Time_Boiler_p(cycle_no) - Time_Boiler_p(cycle_no-1)
      ! h_ox = oxide_thickness(cycle_no) * oxide_thick2si
      ! delta_h2 = oxide_thickness(cycle_no)**2 - oxide_thickness(cycle_no-1)**2
      ! dh_ox = delta_h2 / (2.0 * oxide_thickness(cycle_no))
      !
      ! delta_h = oxide_thickness(cycle_no) - oxide_thickness(cycle_no-1)
      ! bug; for transition cycle_no is not the absolute cycle number
      dstrain_ox = delta_h / lateral_ox_strain_length

      ! lateral_ox_strain_length [microns]
      ox_strain_hoop(i) = ox_strain_hoop_old(i) + dstrain_ox

      ox_strain_axial(i) = ox_strain_hoop(i)

      ox_strain_r(i) = 0.0

      ! store previous value for the next cycle
      ox_strain_hoop_old(i) = ox_strain_hoop(i)

    else 

      write(out_lun, *) 'ERROR: oxidation mode ', i, first_oxide_layer
      stop

    endif

      if (i == 2)  then
        write(out_lun, 2) Time_Boiler_p(cycle_no), &
          & oxide_thickness(cycle_no), ox_strain_hoop(i), &
          & ox_strain_axial(i), ox_strain_r(i)
 2      format('ox_ind_growth_lateral ', 10(1pe13.6, 1x))
      endif

  end do 

  return

  END SUBROUTINE SET_OX_STRAINS

  SUBROUTINE GET_U_JUMP_SURFACE(N, thickness_oxide_layer, iter)

    !=======================================================================
    ! Purpose(s):
    !
    !   get the jump in displacements if there surface oxidation
    !  
    !=======================================================================

  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & u_now_int_jump, ox_strain_r, ox_strain_hoop, &
      & ox_strain_axial, total_oxide_thickness_ref, rr_metal_in_initial
  use oxide_stress_data_module,  only: oxidation_strain_mode, &
      & growth_ox_location, u_int_jump_now, displacement_jump_type, &
      & displacement_jump_fr
  use oxide_data_module,   only: pilling_bedworth_ratio, thickness_fr, &
      & if_temp_steam_or_oxide, first_oxide_layer
  use boiler_data_module, only: tube_outer_radius, tube_thickness, &
      & oxide_location
  use parameter_module,   only: moxide

  ! Argument List
  integer,    intent(IN)  :: N, iter
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer

  ! Local Variables
  integer  :: i, k, j
  real     :: u_jump_fr, total_oxide_thickness_ref_si, sign_u

  ! set the displacement for surface growth
  u_int_jump_now = 0.0

  if (TRIM(oxide_location) == 'inner_metal_surface')  then
 
    sign_u = 1.0

  else if (TRIM(oxide_location) == 'outer_metal_surface')  then

    sign_u = -1.0

  else

    write(6, *) 'ERROR: oxide_location U_JUMP'
    STOP

  endif

  if (.not. if_growth_induce_stress)  RETURN

  ! total oxide thickness in SI units
  total_oxide_thickness_ref_si = &
     & SUM(thickness_oxide_layer(first_oxide_layer:N))

  do k = first_oxide_layer, N  ! interface between k and k - 1 layers

    if (TRIM(growth_ox_location(k)) == 'surface')  then
      ! surface growth on this interface

      if (TRIM(displacement_jump_type(k)) == 'steiner06')  then

        if (N == 2)  then  ! one oxide layer

          u_jump_fr = ox_strain_r(k) / (1.0 + ox_strain_r(k))
          u_int_jump_now(1, k) = sign_u * u_jump_fr * thickness_oxide_layer(k)

          write(aux_lun, *) 'u_jump_ass1 ', ox_strain_r(k), sign_u, &
           & u_jump_fr, thickness_oxide_layer(k), u_int_jump_now(1, k)

        else if (N == 3)  then

          ! most problable this is the second oxide layer, k = N = 3
          u_jump_fr = (ox_strain_r(k) * thickness_oxide_layer(k) - &
               & ox_strain_r(k-1) * thickness_oxide_layer(k-1)) / &
               & (1.0 + ox_strain_r(k))
          u_int_jump_now(1, k) = sign_u * u_jump_fr 

        else if (N > 3)  then
          write(out_lun, *) 'ERROR: steiner06 for three-oxide is not setup '
          STOP

        endif

      else if (TRIM(displacement_jump_type(k)) == 'hsueh83')  then

        u_jump_fr = 1.0
        u_int_jump_now(1, k) = sign_u * u_jump_fr * thickness_oxide_layer(k)

      else if (TRIM(displacement_jump_type(k)) == 'input' .or. &
        & TRIM(displacement_jump_type(k)) == 'value' .or. &
        & TRIM(displacement_jump_type(k)) == 'input_value')  then

        u_jump_fr = displacement_jump_fr(k)
        u_int_jump_now(1, k) = sign_u * u_jump_fr * thickness_oxide_layer(k)

      else 
        write(out_lun, *) 'ERROR: displacement_jump_type must be ', &
           & ' steiner06, hsueh83, or input'
        STOP
      endif

    else if (TRIM(growth_ox_location(k)) == 'interface')  then
      !  growth right on this interface
      u_int_jump_now(1, k) = 0.0
    endif

    ! write(aux_lun, *) 'u_jump_ass2 ', ox_strain_r(k), sign_u, &
    !       & u_jump_fr, thickness_oxide_layer(k), u_int_jump_now(1, k)

  end do

  return

  END SUBROUTINE GET_U_JUMP_SURFACE

  SUBROUTINE OXIDE_INDUCED_SETUP_OH(iter_og, id_high_or_low)

    !=======================================================================
    ! Purpose(s):
    !
    !   setup coefficiencs for oxide induced growth constraints
    !   based on Oh et al. (2006)
    !          
    !  
    !=======================================================================

  use boiler_data_module, only: tube_thickness, tube_outer_radius
  use solution_data_module, only:  npr_st, rad_st,  &
      & rad_int, integral_th, tau_th, tau_th_st, eps0, &
      & strain_rad, strain_hoop, displ_rad, &
      & stress_axial, stress_hoop, stress_rad, bc_st, cc_st
  use oxide_data_module,   only: Youngs_modul, Youngs_temp, &
      & th_exp_coeff, pilling_bedworth_ratio, poisson_ratio
  use output_module,             only: tty_lun, out_lun, inp_lun, aux_lun
  use solver_data_module, only       : neq_solver, fval, jacf, var, &
      & delta_var
  use oxide_stress_data_module,  only: imetal, ispinel, imagnetite, &
      & iout, imet_spin, ispin_mag, i_in
  use oxide_stress_data_module,  only: rc_og, u_og, rr, trace_strain, &
      & h_ox_sp, h_ox_mag, h_ox_total, ir_c_i, ir_c_m_out, ir_r_i_s, &
      & ir_r_o_ma, pbr_sp, pbr_mag, &
      & rr_metal_in_initial, area_fr_metal_2_spinel

  ! Argument List
  integer, intent(IN)  :: iter_og
  integer, intent(IN)  :: id_high_or_low

  ! Local Variables
  integer  :: i, j, neq, nvar
  integer  :: ir_r_i_ma, ir_c_ms, ir_r_i_m, ir_r_o_s, ir_c_sm, &
       & ir_r_m_out
  real     :: Y_alfa, pf_c, pf_b, df, r, r2, strain_hoop_eq, &
            & Integral, tau ! current radius
  integer, save  :: ip = 0

  if (iter_og == 1)  then

    ! initialize the initial outer radius of the tube; reference radius
    rr_metal_in_initial = (tube_outer_radius - tube_thickness) * &
      & (1.0 + tau_th_st(imetal, 1)*(1.0 + poisson_ratio(imetal)))

  endif

    ! setup jacobian for imposing constraints due to oxide growth
    neq = 0
    nvar = 0
    fval = 0.0
    jacf = 0.0

    ! ir_c - stores id for radius in the current system (deformed)
    ! ir_r - stores id for radius in the reference system (stress free)

    !  magnetite surface displacement
    neq = neq + 1
    nvar = nvar + 1
    ir_c_i = nvar  ! index for r_i
    nvar = nvar + 1
    ir_r_i_ma = nvar
    ! r_i  - first variable
    ! R_i,ma = r_i - u_i, ma

    if (iter_og == 1) then
      var(ir_c_i) = rc_og(i_in)
      var(ir_r_i_ma) = var(ir_c_i) - u_og(i_in, 1)
    endif

    fval(neq) = var(ir_r_i_ma) - var(ir_c_i) + u_og(i_in, 1)  
    jacf(neq, ir_r_i_ma) = 1.0
    jacf(neq, ir_c_i) = -1.0
    
    ! current radius at metal-spinel interface
    neq = neq + 1
    nvar = nvar + 1
    ir_c_ms = nvar

    if (iter_og == 1) then
      var(ir_c_ms) = var(ir_c_i) + h_ox_total
    endif

    fval(neq) = var(ir_c_ms) - var(ir_c_i) - h_ox_total   ! r_ms = r_i + h_s + h_ma
    jacf(neq, ir_c_ms) = 1.0
    jacf(neq, ir_c_i) = -1.0

    ! reference radius of the metal
    neq = neq + 1
    nvar = nvar + 1
    ir_r_i_m = nvar
    ir_r_o_s = ir_r_i_m         ! R_o,s = R_i,m
    ! R_i,m = r_i + h_s + h_ma - u_ms

    if (iter_og == 1) then
      var(ir_r_i_m) = var(ir_c_i) + h_ox_total - u_og(imet_spin, 1) 
    endif

    fval(neq) = var(ir_r_i_m) - var(ir_c_i) - h_ox_total + u_og(imet_spin, 1) 
    jacf(neq, ir_r_i_m) = 1.0
    jacf(neq, ir_c_i) = -1.0

    ! current radius at spinel-magnetite interface
    neq = neq + 1
    nvar = nvar + 1
    ir_c_sm = nvar
    ! r_sm = r_i + h_ma

    if (iter_og == 1) then
      var(ir_c_sm) = var(ir_c_i) + h_ox_mag
    endif

    fval(neq) = var(ir_c_sm) - var(ir_c_i) - h_ox_mag   
    jacf(neq, ir_c_sm) = 1.0
    jacf(neq, ir_c_i) = -1.0

    ! current radius at the outer surface of metal
    neq = neq + 1
    nvar = nvar + 1
    ir_c_m_out = nvar  ! index for r_o
    nvar = nvar + 1
    ir_r_m_out = nvar

    if (iter_og == 1) then
      var(ir_c_m_out) = rc_og(iout)
      var(ir_r_m_out) = var(ir_c_m_out) - u_og(iout, 1)  
    endif
    ! 
    fval(neq) = var(ir_r_m_out) - var(ir_c_m_out) + u_og(iout, 1)  
    jacf(neq, ir_r_m_out) = 1.0
    jacf(neq, ir_c_m_out) = -1.0

    ! fm equation density ratio - dilation relationship for metal
    neq = neq + 1
    ! no new variables
    fval(neq) = var(ir_r_m_out)**2 - var(ir_r_i_m)**2 - &
      & (1.0 - trace_strain(imetal)) * (var(ir_c_m_out)**2 - var(ir_c_ms)**2)
    jacf(neq, ir_r_m_out) = 2.0 * var(ir_r_m_out)
    jacf(neq, ir_r_i_m) = - 2.0 * var(ir_r_i_m)
    jacf(neq, ir_c_m_out) = - 2.0 * (1.0 - trace_strain(imetal)) * var(ir_c_m_out)
    jacf(neq, ir_c_ms) = 2.0 * (1.0 - trace_strain(imetal)) * var(ir_c_ms)

    ! fs equation; density ratio - dilation relationship for spinel
    neq = neq + 1
    nvar = nvar + 1
    ir_r_i_s = nvar  ! index for R_i,s

    if (iter_og == 1) then
      var(ir_r_i_s) = rr(ispin_mag, 1)
    endif

    fval(neq) = var(ir_r_o_s)**2 - var(ir_r_i_s)**2 - &
        & (1.0 - trace_strain(ispinel)) * (var(ir_c_ms)**2 - var(ir_c_sm)**2)
    jacf(neq, ir_r_o_s) = 2.0 * var(ir_r_o_s)
    jacf(neq, ir_r_i_s) = - 2.0 * var(ir_r_i_s)
    jacf(neq, ir_c_ms) = - 2.0 * (1.0 - trace_strain(ispinel)) * var(ir_c_ms)
    jacf(neq, ir_c_sm) = 2.0 * (1.0 - trace_strain(ispinel)) * var(ir_c_sm)

    ! f_mag equation; density ratio - dilation relationship for magnetite
    neq = neq + 1
    nvar = nvar + 1
    ir_r_o_ma = nvar  ! index for R_o,ma

    if (iter_og == 1) then
      var(ir_r_o_ma) = rr(ispin_mag, 2)
    endif

    fval(neq) = var(ir_r_o_ma)**2 - var(ir_r_i_ma)**2 - &
        & (1.0 - trace_strain(imagnetite)) * (var(ir_c_sm)**2 - var(ir_c_i)**2)
    jacf(neq, ir_r_o_ma) = 2.0 * var(ir_r_o_ma)
    jacf(neq, ir_r_i_ma) = - 2.0 * var(ir_r_i_ma)
    jacf(neq, ir_c_sm) = - 2.0 * (1.0 - trace_strain(imagnetite)) * var(ir_c_sm)
    jacf(neq, ir_c_i) = 2.0 * (1.0 - trace_strain(imagnetite)) * var(ir_c_i)

    ! pbr is the same for spinel and magnetite
    pbr_sp = pilling_bedworth_ratio
    pbr_mag = pilling_bedworth_ratio

    ! pbr equation
    neq = neq + 1
    ! no new variables
    fval(neq) = pbr_mag * (var(ir_r_o_s)**2 - var(ir_r_i_s)**2) + &
           & pbr_sp * (var(ir_r_o_ma)**2 - var(ir_r_i_ma)**2) - &
           & pbr_mag * pbr_sp * (var(ir_r_o_s)**2 - rr_metal_in_initial**2)
    jacf(neq, ir_r_o_s) = 2.0 * var(ir_r_o_s) * pbr_mag * (1.0 - pbr_sp)
    jacf(neq, ir_r_i_s) = - 2.0 * pbr_mag * var(ir_r_i_s)
    jacf(neq, ir_r_o_ma) = 2.0 * pbr_sp * var(ir_r_o_ma)
    jacf(neq, ir_r_i_ma) = - 2.0 * pbr_sp * var(ir_r_i_ma)

    neq_solver = neq

    if (.not. neq_solver == nvar)  then

      write(tty_lun, *) 'ERROR: oxide-induced-stress nvar ne neq'
      stop

    endif

    do i = 1, neq_solver
      write(aux_lun, 2) i, (jacf(i, j), j= 1, neq_solver), fval(i)
 2    format('oxis_jac ', i2, 1x, 20(1pe13.6, ' '))
    end do

    ! store the fractional area of the metal that converts to spinel
    ! f_sg
    area_fr_metal_2_spinel = (var(ir_r_o_s)**2 - var(ir_r_i_s)**2) / &
        & (pbr_sp * (var(ir_r_o_s)**2 - rr_metal_in_initial**2))

    write(aux_lun, 3) iter_og, area_fr_metal_2_spinel
 3    format('a_metal2spinel ', i2, 1x, 20(1pe13.6, ' '))

  return

  END SUBROUTINE OXIDE_INDUCED_SETUP_OH

  SUBROUTINE OXIDE_INDUCED_UPDATE(iter_og, id_high_or_low)

    !=======================================================================
    ! Purpose(s):
    !
    !   update geometry and displacement condition for a new iteration on
    !     induced growth constraints
    !          
    !  
    !=======================================================================

  use solution_data_module, only:  npr_st, rad_st,  &
      & rad_int, integral_th, tau_th, tau_th_st, eps0, &
      & strain_rad, strain_hoop, displ_rad, &
      & stress_axial, stress_hoop, stress_rad, bc_st, cc_st
  use oxide_data_module,   only: Youngs_modul, Youngs_temp, &
      & th_exp_coeff, pilling_bedworth_ratio
  use output_module,             only: tty_lun, out_lun, inp_lun, aux_lun
  use solver_data_module, only       : neq_solver, fval, jacf, var, &
      & delta_var
  use oxide_stress_data_module,  only: imetal, ispinel, imagnetite, &
      & iout, imet_spin, ispin_mag, i_in
  use oxide_stress_data_module,  only: rc_og, u_og, rr, trace_strain, &
      & h_ox_sp, h_ox_mag, h_ox_total, ir_c_i, ir_c_m_out, ir_r_i_s, &
      & ir_r_o_ma, pbr_sp, pbr_mag, if_growth_induce_stress, &
      & u_now_int_jump, rr_metal_in_initial

  ! Argument List
  integer, intent(IN)  :: iter_og
  integer, intent(IN)  :: id_high_or_low

  ! Local Variables
  integer  :: i, j, neq, nvar
  integer  :: ir_r_i_ma, ir_c_ms, ir_r_i_m, ir_r_o_s, ir_c_sm, &
       & ir_r_m_out
  real     :: Y_alfa, pf_c, pf_b, df, r, r2, strain_hoop_eq, &
            & Integral, tau ! current radius
  integer, save  :: ip = 0

    ! ir_c - stores id for radius in the current system (deformed)
    ! ir_r - stores id for radius in the reference system (stress free)

    rc_og(i_in) = var(ir_c_i)                       ! r_i
    rc_og(iout) = var(ir_c_m_out)                   ! r_o
    rc_og(imet_spin) = rc_og(i_in) + h_ox_total     ! metal-spinel interface
    rc_og(ispin_mag) = rc_og(i_in) + h_ox_mag       ! spinel-magnetite interface

    ! u_int_jump = u_spinel - u_magnetite at spinel-magnetite interface
    !            = u(2) - u(3) at the interface between layers 2 and 3
    u_now_int_jump(id_high_or_low) = var(ir_r_o_ma) - var(ir_r_i_s)
    
    ! the fractional area of the metal that converts to spinel
    ! was determined in the setup part since not all indices are stored f_sg
    ! area_fr_metal_2_spinel

  return

  END SUBROUTINE OXIDE_INDUCED_UPDATE


END MODULE INDUCED_STRESS_MODULE
