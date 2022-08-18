MODULE PROPERTY_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all combustion data structures integer and scalar parameters.
  !
  !  Author: Adrian S. Sabau, sabaua@ornl.gov
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================
  use parameter_module, only: mfuel, mcomb, mcp, mtref, mtint, ms, mh, &
                      & mspecie, mtrans, mgas, mstates
  implicit none

  save
  ! Public Module
  public

  ! the name of the specie from the CEA-NASA transport database
  character(LEN = 80), dimension(0:mspecie) ::  specie_trans_name

  ! the name of the interaction specie from the CEA-NASA transport database
  character(LEN = 80), dimension(0:mspecie) ::  specie_trans_inter_name

  ! actual number of species introduced
  integer ::   no_trans_species

  ! air id for transport properties
  integer ::   no_trans_air

  ! actual number of temperature intervals for each property
  integer, dimension(0:mspecie) ::  no_visc_temp_domain, no_conduct_temp_domain

  ! coefficients for viscosity and thermal conductivity from 
  ! CEA-NASA transport database
  real, dimension(0:mspecie, mtref, mtrans)    ::  viscosity_specie_coeff, &
      & thermal_conduc_specie_coeff

  ! parameters for CEA-NASA model
  ! temperature points for which each property is given
  real, dimension(0:mspecie, mtref+1) :: Tref_specie_cond, Tref_specie_visc

END MODULE PROPERTY_DATA_MODULE
