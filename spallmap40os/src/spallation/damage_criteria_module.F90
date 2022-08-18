MODULE DAMAGE_CRITERIA_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for spallation. 
    !
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================


  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: DAMAGE_SCALE_OUTAGE

 CONTAINS

  SUBROUTINE DAMAGE_SCALE_OUTAGE()
    !=======================================================================
    ! Purpose(s):
    !
    !  find the intersection between damage map criteria Ak/sqrt(d), and
    !    linear curves reprezenting the min or maximum for the strain 
    !    evolution during an outage
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, mdamage, pi
    use output_module,      only: tty_lun, out_lun, aux_lun
    use curve_fit_module,   only: LINEAR_FIT
    use oxide_data_module,  only: const_damage_map


    ! call LINEAR_FIT(X,Y,NDATA,SIG,MWT,A,B,SIGA,SIGB,CHI2,Q)

  RETURN

  END SUBROUTINE DAMAGE_SCALE_OUTAGE

END MODULE DAMAGE_CRITERIA_MODULE
