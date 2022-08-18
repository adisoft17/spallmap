MODULE INDUCED_OXIDE_STRESS_MODULE

    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for oxided-induced growth stresses. 
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================


  implicit none

  ! Private Module
  private

  ! Public Procedures
  public ::   STRAIN_OX_IND_TERM_RADIAL, STRAIN_OX_IND_TERM_HOOP, &
       & U_TERM_STRAIN_OX_IND, DUDR_TERM_STRAIN_OX_IND, &
       & STRESS_AXIAL_OX_IND_TERM

 CONTAINS

  REAL FUNCTION STRAIN_OX_IND_TERM_RADIAL(poisson, poisson_star, r, &
     & ox_strain_r, ox_strain_hoop, ox_strain_axial)

  ! contribution to the oxidation induced strains to the radial stress
  ! s_r / E^* = b/(1-nu*) - c/(r^2 * (1+ nu*)) - J_ri,r/r^2 + 
  !       eps0 * nu/(1-nu*) - Fr

  real, Intent(IN)  :: poisson, poisson_star, r, &
     & ox_strain_r, ox_strain_hoop, ox_strain_axial

  real :: ln1

  ln1 = (1.0 - poisson_star) * LOG(ABS(r))

  STRAIN_OX_IND_TERM_RADIAL = (poisson * ox_strain_axial + &
     & 0.5 * (1.0 - ln1) * ox_strain_r + &
     & 0.5 * (1.0 + ln1) * ox_strain_hoop) / (1.0 - poisson_star)

  return

  END FUNCTION STRAIN_OX_IND_TERM_RADIAL

  REAL FUNCTION STRAIN_OX_IND_TERM_HOOP(poisson, poisson_star, r, &
     & ox_strain_r, ox_strain_hoop, ox_strain_axial)

  ! contribution to the oxidation induced strains to the radial stress
  ! s_teta / E^* = b/(1-nu*) - tau + c/(r^2 * (1+ nu*)) + J_ri,r/r^2 + 
  !                eps0 * nu/(1-nu*) - Fhoop

  real, Intent(IN)  :: poisson, poisson_star, r, &
     & ox_strain_r, ox_strain_hoop, ox_strain_axial

  real :: ln1

  ln1 = (1.0 - poisson_star) * LOG(ABS(r))

  STRAIN_OX_IND_TERM_HOOP = (poisson * ox_strain_axial + 0.5 * &
     & (poisson_star - ln1) * ox_strain_r + 0.5 * &
     & (2.0 - poisson_star + ln1) * ox_strain_hoop) / (1.0 - poisson_star)

  return

  END FUNCTION STRAIN_OX_IND_TERM_HOOP

  REAL FUNCTION U_TERM_STRAIN_OX_IND(poisson_star, r, &
     & ox_strain_r, ox_strain_hoop)

  ! plain strain relationships 
  ! u = r * b + c/r + (1+nu*) * J_ri,r/r + u_ox_term

  real, Intent(IN)  :: poisson_star, r, &
     & ox_strain_r, ox_strain_hoop

  U_TERM_STRAIN_OX_IND = (1.0 - poisson_star) * &
     & (ox_strain_r - ox_strain_hoop) * 0.5 * r * LOG(ABS(r))

  return

  END FUNCTION U_TERM_STRAIN_OX_IND

  REAL FUNCTION DUDR_TERM_STRAIN_OX_IND(poisson_star, r, &
     & ox_strain_r, ox_strain_hoop)

  ! du/dr term due to oxide growth 
  ! u = r * b + c/r + (1+nu*) * J_ri,r/r + u_ox_term
  ! du/dr = b + (1+nu*)*tau - c/r2 - (1+nu*) * J_ri,r/r2 + 
  !         d(u_ox_term)/dr

  real, Intent(IN)  :: poisson_star, r, &
     & ox_strain_r, ox_strain_hoop

  DUDR_TERM_STRAIN_OX_IND = (1.0 - poisson_star) * &
     & (ox_strain_r - ox_strain_hoop) * 0.5 * (LOG(ABS(r)) + 1.0)

  return

  END FUNCTION DUDR_TERM_STRAIN_OX_IND

  REAL FUNCTION STRESS_AXIAL_OX_IND_TERM(poisson, poisson_star, &
     & Youngs, Youngs_star, ox_strain_r, ox_strain_hoop, &
     & ox_strain_axial)

  ! contribution to the oxidation induced strains to the radial stress
  ! s_teta / E^* = b/(1-nu*) - tau + c/(r^2 * (1+ nu*)) + J_ri,r/r^2 + 
  !                eps0 * nu/(1-nu*) - Fhoop
  ! sigma_z_ox = - A + B * ln(r); this routine is for A

  real, Intent(IN)  :: poisson, poisson_star, &
     & Youngs, Youngs_star, ox_strain_r, ox_strain_hoop, ox_strain_axial

  real :: ax1, rad1

  ax1 = Youngs + 2.0 * Youngs_star * poisson**2 / (1.0 - poisson_star)
  rad1 = 0.5 * Youngs_star * poisson / (1.0 - poisson_star)

  STRESS_AXIAL_OX_IND_TERM = ax1 * ox_strain_axial + &
     & rad1 * ((poisson_star + 1) * ox_strain_r + &
     & (3.0 - poisson_star) * ox_strain_hoop)

  return

  END FUNCTION STRESS_AXIAL_OX_IND_TERM

  END MODULE INDUCED_OXIDE_STRESS_MODULE
