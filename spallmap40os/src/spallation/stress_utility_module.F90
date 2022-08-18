MODULE STRESS_UTILITY_MODULE
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
  public :: RADIAL_STRESS, HOOP_STRESS, AXIAL_STRESS, RADIAL_DISPL, &
          & HOOP_STRAIN, RADIAL_STRAIN, GET_OG_DIM_DISPL 

  CONTAINS

  REAL FUNCTION RADIAL_STRESS(Y_alfa, Integral_i_r, r2, bc, cc, df)

  ! plain strain relationships 
  ! s_r = E^* (b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2)
  ! with eps0 becomes:
  ! s_r = E^* [ b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2 + eps0 * nu/(1-nu*)]

  real, Intent(IN)  :: Y_alfa, Integral_i_r, r2, bc, cc, df

  RADIAL_STRESS = Y_alfa * (bc - cc / r2 - Integral_i_r / r2 + df)

  return

  END FUNCTION RADIAL_STRESS

  REAL FUNCTION HOOP_STRESS(Y_alfa, Integral_i_r, r2, bc, cc, df, tau_i)

  real, Intent(IN)  :: Y_alfa, Integral_i_r, r2, bc, cc, df, tau_i
  ! s_teta = E^* (b/(1-nu*) - tau + c/(r^2 * (1+ nu*)) + J_ri,r/r^2)
  ! with eps0 becomes:
  ! s_teta = E^* [b/(1-nu*) - tau + c/(r^2 * (1+ nu*)) + J_ri,r/r^2 + eps0 * nu/(1-nu*)]

  HOOP_STRESS = Y_alfa * (bc - tau_i + cc / r2 + Integral_i_r / r2  + df)

  return

  END FUNCTION HOOP_STRESS

  REAL FUNCTION AXIAL_STRESS(Ymod, Total_Axial_Strain, &
       & pf_eps, pf_b, pf_tau, bc, tau_i)

      ! poisson factor
      ! pf_eps = (1.0 - poisson(i)) / ((1.0 - poisson(i)) * &
      !       & (1.0 - 2.0 * poisson(i)))
      ! 
      ! pf_b = 2.0 * poisson(i) * th_exp(i) / (1.0 - poisson(i))
      !
      ! pf_tau = th_exp(i) / (1.0 - poisson(i))
      ! 
      ! Total_Axial_Strain = eps0

  real, Intent(IN)  :: Ymod, Total_Axial_Strain, &
       & pf_eps, pf_b, pf_tau, bc, tau_i

  AXIAL_STRESS = Ymod * (Total_Axial_Strain * pf_eps + &
       & pf_b * bc - pf_tau * tau_i)

  return

  END FUNCTION AXIAL_STRESS

  REAL FUNCTION RADIAL_DISPL(Integral_i_r, r, bc, cc, pf_int)

  ! plain strain relationships 
  ! u = r * b + c/r + (1+nu*) * J_ri,r/r

  real, Intent(IN)  :: Integral_i_r, r, bc, cc, pf_int

  RADIAL_DISPL = bc * r + (pf_int * Integral_i_r + cc) / r 

  return

  END FUNCTION RADIAL_DISPL

  REAL FUNCTION HOOP_STRAIN(Integral_i_r, r, bc, cc, pf_int)

  ! plain strain relationships 
  ! Y_alfa = Youngs * th_exp (these are both star values)

  real, Intent(IN)  :: Integral_i_r, r, bc, cc, pf_int

  HOOP_STRAIN = RADIAL_DISPL(Integral_i_r, r, bc, cc, pf_int) / r

  return

  END FUNCTION HOOP_STRAIN

  REAL FUNCTION RADIAL_STRAIN(Integral_i_r, tau_i_r, r2, bc, cc, pf_tau)

  ! plain strain relationships 
  ! du/dr = b + (1+nu*)*tau - c/r2 - (1+nu*) * J_ri,r/r2
 
  real, Intent(IN)  :: Integral_i_r, tau_i_r, r2, bc, cc, pf_tau

  RADIAL_STRAIN = pf_tau * tau_i_r + bc - (cc + pf_tau * Integral_i_r) / r2

  return

  END FUNCTION RADIAL_STRAIN

  SUBROUTINE GET_OG_DIM_DISPL(N, poisson, th_exp, Youngs, &
       & poisson_star, th_exp_star, Youngs_star, if_update_displ_dim)
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the dimensions and displacements to be used for oxide induced 
    !          growth
    !  
    !=======================================================================
  use parameter_module,   only: moxide, oxide_thick2si, id_low, id_high, &
      & mgrid
  use solution_data_module, only:  npr_st, rad_st,  &
      & rad_int, integral_th, tau_th, tau_th_st, eps0, &
      & strain_rad, strain_hoop, displ_rad, &
      & stress_axial, stress_hoop, stress_rad, bc_st, cc_st
  use oxide_data_module,   only: Youngs_modul, Youngs_temp, &
      & th_exp_coeff, thickness_fr
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,   only: rc_og, u_og, rr, trace_strain, &
      & h_ox_sp, h_ox_mag, h_ox_total, rr_metal_in_initial, &
      & total_oxide_thickness_ref, trace_total, trace_without_th
  use oxide_stress_data_module,  only: imetal, ispinel, imagnetite, &
      & iout, imet_spin, ispin_mag, i_in, if_th_oxide_growth, &
      & if_exclude_th_strain_trace, stress_max_oxide_manning
  use property_module, only: OPERAND_ARRAY

  ! Argument List
  integer,    intent(IN)          :: N
  real, dimension(N), intent(IN)  :: poisson, th_exp, Youngs, &
       & poisson_star, th_exp_star, Youngs_star
  logical, intent(IN)             :: if_update_displ_dim

  ! Local Variables
  integer  :: i, j, id_layer, loc1
  real     :: Y_alfa, pf_c, pf_b, df, r, r2, strain_hoop_eq, &
            & Integral, tau   ! current radius
  integer, save  :: ip = 0
  character(LEN = 80)   :: operation_type = 'ave_sum'

  if (.not. if_update_displ_dim)  RETURN  ! update displacements only for oh
  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  ! initialize
  ! if (if_first_temp_rad)   then cannot do it since oxide thickness
  ! varies

  if (N > 3)  then
    write(tty_lun, *) 'ERROR: growth induces stresses for one or two oxide '
    stop
  else if (N < 3)  then
    write(tty_lun, *) 'ERROR: growth induces stresses for one or two oxide '
    stop
  endif

  ! set radius for the constraints due to oxide-induces stress
  ! one metal, two oxide layers
    ! rc_og - radius in current state for oxide growth
    ! radius is not double valued in the current system of reference

    rc_og(iout) = rad_int(1)      ! outer radius
    rc_og(imet_spin) = rad_int(2) ! metal-spinel interface
    rc_og(ispin_mag) = rad_int(3) ! spinel-magnetite interface
    rc_og(i_in) = rad_int(4)       ! inner radius 

    ! write(aux_lun, 1) rc_og(i_in), rc_og(ispin_mag), rc_og(imet_spin), rc_og(iout)

    ! set displacements
    ! 1 - outer side of the surface
    ! 2 - inner side of the surface
    id_layer = 1   ! metal now
    Integral = (1.0 + poisson(id_layer)) * Integral_th(id_layer, 1)
    u_og(iout, 1) = RADIAL_DISPL(Integral, rc_og(iout), &
          & bc_st(id_layer), cc_st(id_layer), 1.0 + poisson_star(id_layer))

    ! on the metal side; same id_layer as above
    Integral = (1.0 + poisson(id_layer)) * Integral_th(id_layer, npr_st(id_layer))
    u_og(imet_spin, 1) = RADIAL_DISPL(Integral, rc_og(imet_spin), &
          & bc_st(id_layer), cc_st(id_layer), 1.0 + poisson_star(id_layer))

    ! on the spinel side (metal-spinel)
    id_layer = 2   ! spinel now
    Integral = (1.0 + poisson(id_layer)) * Integral_th(id_layer, 1)
    u_og(imet_spin, 2) = RADIAL_DISPL(Integral, rc_og(imet_spin), &
          & bc_st(id_layer), cc_st(id_layer), 1.0 + poisson_star(id_layer))

    ! on the spinel side (spinel-magnetite)
    Integral = (1.0 + poisson(id_layer)) * Integral_th(id_layer, npr_st(id_layer))
    u_og(ispin_mag, 1) = RADIAL_DISPL(Integral, rc_og(ispin_mag), &
          & bc_st(id_layer), cc_st(id_layer), 1.0 + poisson_star(id_layer))

    ! magnetite side (spinel-magnetite)
    id_layer = 3   ! magnetite now
    Integral = (1.0 + poisson(id_layer)) * Integral_th(id_layer, 1)
    u_og(ispin_mag, 2) = RADIAL_DISPL(Integral, rc_og(ispin_mag), &
          & bc_st(id_layer), cc_st(id_layer), 1.0 + poisson_star(id_layer))

    ! magnetite surface
    Integral = (1.0 + poisson(id_layer)) * Integral_th(id_layer, npr_st(id_layer))
    u_og(i_in, 1) = RADIAL_DISPL(Integral, rc_og(i_in), &
          & bc_st(id_layer), cc_st(id_layer), 1.0 + poisson_star(id_layer))

    ! dimensions of the reference frame
    rr(iout, 1) = rc_og(iout) - u_og(iout, 1)                  ! R_m^o
    rr(imet_spin, 1) = rc_og(imet_spin) - u_og(imet_spin, 1)   ! R_m^i
    rr(imet_spin, 2) = rc_og(imet_spin) - u_og(imet_spin, 2)   ! R_sp^o
    rr(ispin_mag, 1) = rc_og(ispin_mag) - u_og(ispin_mag, 1)   ! R_sp^i
    rr(ispin_mag, 2) = rc_og(ispin_mag) - u_og(ispin_mag, 2)   ! R_mag^o
    rr(i_in, 1) = rc_og(i_in) - u_og(i_in, 1)                  ! R_mag^o

    write(aux_lun, 1) -(rr(i_in, 1) - rr(ispin_mag, 2)), &
         & -(rr(ispin_mag, 1) - rr(imet_spin, 2)), &
         & -(rr(imet_spin, 1) - rr(iout, 1))
 1  format('mag_thick= ', 1pe13.6, ' sp_thick= ', 1pe13.6, ' metal_thick=', &
         & 1pe13.6)

    ! get the dilatational components; they are constant within each layer
    ! e_r + e_theta + e_z - 3 * alfa * delta_T; how about if you have creep; subtract that, too
    ! how about the displacements; 

    ! get average trace per each layer
    ! if you just use average for each strain may not be ok; use the average of the sum
    ! get the trace at each stress location, then average it
  
    do i = 1, 3

      if (if_exclude_th_strain_trace)  then

        call OPERAND_ARRAY(1, mgrid, trace_without_th(i, :), 1, npr_st(i), &
           & operation_type, trace_strain(i), loc1)

      else

        call OPERAND_ARRAY(1, mgrid, trace_total(i, :), 1, npr_st(i), &
           & operation_type, trace_strain(i), loc1)

      endif

    end do

    !  get the oxide thickness at this NEW current temperature
    ! spinel is the second layer
    ! total_oxide_thickness_ref = h_oxide_total(RT)

    !       tau = (1.0 + poisson(i)) * tau_th_st(i, j)

    if (if_th_oxide_growth)  then

      h_ox_sp = total_oxide_thickness_ref * thickness_fr(ispinel) * oxide_thick2si * &
            & (1.0 + tau_th_st(ispinel, 1)*(1.0 + poisson(ispinel)))
      h_ox_mag = total_oxide_thickness_ref * thickness_fr(imagnetite) * oxide_thick2si * &
            & (1.0 + tau_th_st(imagnetite, npr_st(imagnetite))*(1.0 + poisson(imagnetite)))

    else

      h_ox_sp = total_oxide_thickness_ref * thickness_fr(ispinel) * oxide_thick2si
      h_ox_mag = total_oxide_thickness_ref * thickness_fr(imagnetite) * oxide_thick2si

    endif

    h_ox_total = h_ox_sp + h_ox_mag

  return

  END SUBROUTINE GET_OG_DIM_DISPL

END MODULE STRESS_UTILITY_MODULE
