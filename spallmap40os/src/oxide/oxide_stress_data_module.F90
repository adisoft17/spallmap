MODULE OXIDE_STRESS_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all oxide data structures integer and scalar parameters.
  !
  !
  !
  !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================
  use parameter_module, only: mfuel, mcomb, mcp, mtref, mtint, ms, mh, &
                            & mspecie, string_dim, string_len, mdiluent, &
                            & moxidant
  use parameter_module, only: moxide_layer, moxide, mtime, mrate, mgrid

  implicit none
  save
  ! Public Module
  public

  ! rc_og - radius in current state for oxide growth
  ! radius is not double valued in the current system of reference
  real, dimension(1:moxide_layer) :: rc_og, trace_strain

  ! u_og - displacemens at interface to be used for oxide growth
  ! u_og(k, 1) - displacement at k interface on the outer side of interface
  ! u_og(k, 2) - displacement at k interface on the inner side of interface
  ! u_og(k, 2) | u_og(k, 1)
  real, dimension(1:moxide_layer, 1:2) :: u_og, rr

  real, dimension(moxide, mgrid)  :: trace_total, trace_without_th

  ! oxide thickness in spinel, magnetite, and total
  real       :: h_ox_sp, h_ox_mag, h_ox_total
  
  ! pilling_bedworth_ratio per each layer
  real       :: pbr_sp, pbr_mag

  ! indices for different layers
  integer    :: imetal = 1, ispinel = 2, imagnetite = 3

  ! indices for radius in current, deformed, state
  integer    :: iout = 1, imet_spin = 2, ispin_mag = 3, i_in = 4

  ! if_exclude_th_strain_trace = .true. then e_th is subtracted from total_trace
  ! if_exclude_th_strain_trace = .false. then total trace is used in 
  ! induced_oxide stresses
  logical    :: if_exclude_th_strain_trace

  ! if_th_oxide_growth = .true., h_oxide(T) = h_oxide(RT) * (1 + CTE*(T-RT))
  ! if_th_oxide_growth = .false., h_oxide(T) = h_oxide(RT)
  ! if_th_tube = .false., rad_int(1) = tube_outer_radius, rad_int(2) = rad_int(1) - tube_thickness
  ! if_th_tube = .true., tube_outer_radius and tube_thickness need to account for 
  ! thermal expansion
  logical    :: if_th_oxide_growth, if_th_tube

  ! if_growth_induce_stress = .true. then consider the oxide growth induced stresses
  ! if_growth_induce_stress = .false. no oxide growth induced stresses
  logical    :: if_growth_induce_stress

  integer    :: ir_c_i     ! index for r_i
  integer    :: ir_c_m_out ! index for r_o
  integer    :: ir_r_i_s   ! index for R_i,s
  integer    :: ir_r_o_ma  ! index for R_o,ma

  ! rr_metal_in_initial = radius of the outer surface of the metal in the
  !  reference state 
  ! at high temperature R_initial(T) = R_initial(RT) * (1+CTE*(T-RT))
  ! at low temperature what R_initial do we use
  ! R_initial(T_high) or R_initial(T_low)

  real       :: rr_metal_in_initial 

  ! total_oxide_thickness_ref = h_oxide_total(RT)
  real       :: total_oxide_thickness_ref 

  ! u_int_jump = u_spinel - u_magnetite at spinel-magnetite interface
  ! u_int_jump(k) - k is either id_low or id_high to indicate full or partial  load
  real, dimension(1:2)       :: u_now_int_jump, u_old_int_jump
  real, dimension(1:2, 1:moxide_layer)   :: u_int_jump_now, u_int_jump_old

  ! area_fr_metal_2_spinel = fractional area of the metal that converts to spinel
  real       :: area_fr_metal_2_spinel

  ! oxidation strains; eigen strains
  real, dimension(moxide_layer)   :: ox_strain_r, ox_strain_hoop, &
     & ox_strain_hoop_old, ox_strain_axial

  ! oxidation_strain_mode = 'radial' or 'dilational' to set ox_strain_?
  ! growth_ox_location = 'surface', 'outer_interface', 'inner_interface', 'interface'
  ! top oxide can grow at the surface; 
  ! itermediate oxides grow at the outer_interface or inner interface
  character(LEN = 80), dimension(0:moxide) ::  oxidation_strain_mode, &
      & growth_ox_location

  ! type of jump in displacements at interface due to 
  ! growth of the oxide at the surface 
  ! displacement_jump_type must be steiner06, hsueh83, or 
  ! input (value, input_value)
  character(LEN = 80), dimension(0:moxide) ::  displacement_jump_type

  ! delta_u = displacement_jump_fr * oxide_thickness 
  ! when displacement_jump_type = input (value, input_value)
  real, dimension(0:moxide_layer)   :: displacement_jump_fr

  ! maximum stress in the oxide due to oxide-growth-induced stresses
  ! Manning (1981); as given in formulae 50 by Steiner
  real :: stress_max_oxide_manning

  real :: d_oxide_thickness_4_strain_hoop

 END MODULE OXIDE_STRESS_DATA_MODULE
