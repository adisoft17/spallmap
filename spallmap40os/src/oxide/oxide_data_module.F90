MODULE OXIDE_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all oxide data structures integer and scalar parameters.
  !
  !! Author: Adrian S. Sabau, sabaua@ornl.gov oxide_data_module.F90
  !
  !=======================================================================
  use parameter_module, only: mfuel, mcomb, mcp, mtref, mtint, ms, mh, &
                            & mspecie, string_dim, string_len, mdiluent, &
                            & moxidant
  use parameter_module, only: moxide_layer, moxide, mtime, mrate, mgrid, &
       & mdamage

  implicit none
  save
  ! Public Module
  public

  ! h is for enthalpy
  ! s is for entropy

  ! if_temp_steam_or_oxide = .false., use surface oxide temperature to 
  ! grow the oxide
  ! if_temp_steam_or_oxide = .true., use steam temperature to 
  ! grow the oxide
  logical  :: if_temp_steam_or_oxide

  ! number of oxide layers, number of oxide species; the first species is 
  ! the metal by default; thus 
  ! isubstrate stores the id of the substrate or tube material
  integer :: no_layer, no_oxide, isubstrate

  ! number of metal layers; no_metal_layers = 1 by default
  ! first_oxide_layer = 2 by default when metal is only one layer
  integer :: no_metal_layers, first_oxide_layer, no_layers_oxide

  ! radius_layer_inside, radius_layer_outside used to define the radii 
  ! metal layer k is between radius_layer_inside(k) < radius_layer_outside(k)
  real, dimension(0:moxide)    :: radius_layer_inside, radius_layer_outside

  ! material index for average scale properties
  integer :: no_oxide_ave

  ! the name of the material; this helps with the properties
  ! oxide_material is to keep track of layer compositions
  character(LEN = 80), dimension(0:moxide) ::  material_name, layer_name

  ! lists materials (compounds) in the current layer
  character(LEN = 80), dimension(0:moxide) ::  layer_materials_name

  ! layer_mat_id - stores material (compound) id for averaging purposes
  integer, dimension(0:moxide, 0:moxide)    :: layer_mat_id

  ! no_layer_mat - numbers of compounds per each layer
  integer, dimension(0:moxide)              :: no_layer_mat

  ! what profile will be used for averaging; linear or log
  character(LEN = 80), dimension(0:moxide) :: profile_function_ave

  character(LEN = 80)   :: oxide_thick_units, oxide_rate_units, const_energy_units

  ! temp_measure_ox_thickness - used to specified the temperature at which
  ! oxide thickness was measured from experiments
  ! temp_measure_ox_thickness = 'room_temperature', 
  ! 'corrected_for_high_temperature', 'average_temperature'
  ! temp_measure_ox_thickness - not implemented
  character(LEN = 80)   :: temp_measure_ox_thickness

  ! growth_ox_temperature = inner_surface, outer_surface, average_layer, 
  ! inner_surface = use temperature at the inner surface (closed to metal)
  ! inner_surface = could be used for spinel
  ! outer_surface = use temperature at the outer surface (closed to steam)
  ! outer_surface = for magnetite
  ! average_layer = use the average temperature in the layer
  ! average_all_layers = use the average temperature in all the layers
  ! steam = use steam temperature; equivalent to if_temp_steam_or_oxide = .true.
  ! oxide_surface = use oxide surface temperature; equivalent to if_temp_steam_or_oxide = .false.,
  character(LEN = 80), dimension(0:moxide)    :: growth_ox_temperature

  ! thickness_growth_fr - fraction by which this oxide grows
  ! since const_energy_oxide is for the entire thickness
  ! for this layer, A = const_energy_oxide * thickness_growth_fr**2
  real, dimension(0:moxide)    :: thickness_growth_fr

  ! fraction_mat_dmin(i) - fraction of material (compound) i for current layer
  ! at the location closest to the metal 
  ! _dmax is at the location furthest to the metal
  !         used to obtain average material properties
  ! spatial dependence cannot be considered for analytical solution
  ! use fraction_material
  real, dimension(0:moxide, 0:moxide)  :: fraction_mat_dmin, &
      & fraction_mat_dmax, fraction_material
  
  ! unit factor for oxide thickness [conversion to microns]
  real  :: oxide_thick_uf, const_energy_uf

  ! kelvin degrees 
  real  :: temp_kelvin_unit

  ! real, dimension(1:moxide_layer, 1:moxide) :: oxide_fraction

  ! moxide_layer should be the same as moxide 
  ! since there is one oxide per layer
  ! thickness_fr_test = thickness_fr at input to be used for testing only
  ! thickness_fr = thickness fraction of each layer of oxide
  real, dimension(0:moxide)     :: thickness_fr, thickness_fr_test, &
       & oxide_thickness
  ! e_th_mat is the thermal strain due to the layer material, 
  ! not combined with the substrate
  ! e_th is the thermal strain in the layer, include the substrate effect
  real, dimension(1:moxide)     :: e_th_mat, sigma_th, e_th

  ! points to materials in each layer, tube is first
  integer, dimension(1:moxide)     :: id_mat

  ! oxide growth rate based on activation energy law
  ! kp = const_energy_oxide * exp(- activ_energy_oxide / (R * T))
  real :: activ_energy_oxide, const_energy_oxide

  ! exponential form of kp
  ! activ_energy_oxide_r = activ_energy_oxide / R
  ! log_const_energy_oxide = 
  ! kp = exp(log_const_energy_oxide - activ_energy_oxide_r / T)
  real :: log_const_energy_oxide, activ_energy_oxide_r

  ! if_function_or_table = .false. then use oxide_rate_log at discrete 
  ! temperatures, i.e., table look up
  ! if_function_or_table = .true. then use the function to evalulate kp
  logical :: if_function_or_table_oxide

  ! pilling_bedworth_ratio = h_oxide/h_metal consumed
  real  :: pilling_bedworth_ratio

  ! Bernstein (1987) factor for reducing the volumetric strains due to pbr
  real  :: pbr_bernstein_factor

  ! lateral_ox_strain_length [microns] length scale to be used with
  ! (oxidation_strain_mode(i) == 'lateral') 
  ! lateral_ox_strain_length = 50-1,000 microns Panicaud et al. (2006)
  ! lateral_ox_strain_length = 1,000 micron recommended
  real  :: lateral_ox_strain_length

  ! ox_hoop_strain_factor, ox_axial_strain_factor valid for 
  ! oxidation_strain_mode(i) == 'radial_dominant'
  ! ox_strain_hoop(i) = ox_hoop_strain_factor * ox_strain_r(i)
  real  :: ox_hoop_strain_factor, ox_axial_strain_factor

  ! damage, flaws, porosity, defects data
  real  :: amplitude_roughness, &
      & wavelength_roughness, delamination_size_existent, &
      & separated_area, porosity_fraction, &
      & number_cracks_volumetric

  ! if this is > 0, then 
  ! delamination_size_existent = thickness * delamination_size_factor
  real  :: delamination_size_factor

  ! oxide growth rate constant as input for cases in which the 
  ! rate law does not apply
  real, dimension(mrate)  :: oxide_rate_value, oxide_rate_temp, &
       & oxide_rate_log, oxide_rate_1temp

  ! number of data points in the oxide_rate_value
  integer :: nrate

  integer  :: id_fe2o3, id_fe3o4

  ! data used in the code for each oxide layer; density at room temperature
  real, dimension(0:moxide) :: rho, poisson_ratio, &
      & surf_fracture_energy, fracture_toughness

  ! data used in the code for each oxide layer; 
  ! property data for tube and other oxide
  ! used to be moxide, mcp ; choose the largest
  real, dimension(0:moxide, 0:moxide) :: cp_value, cp_temp, &
      & cond_value, cond_temp, rho_value, rho_temp, &
      & Youngs_modul, Youngs_temp, &
      & th_exp_coeff, th_exp_temp

  ! data used in the code for each oxide layer; 
  ! number of intervals for the thermal expansion, youngs, cond, cp
  ! number of data points = number of intervals + 1
  integer, dimension(0:moxide) :: nth_exp, nyoungs, ncond, ncp_prop

  ! data read at input for each material or compound; 
  ! data used to obtain properties for each layer based on averaging
  ! density at room temperature
  real, dimension(0:moxide) :: rho_mat, poisson_ratio_mat, &
      & surf_fracture_energy_mat, fracture_toughness_mat

  ! data read at input for each material or compound; 
  ! data used to obtain properties for each layer based on averaging
  ! property data for tube and other oxide
  ! used to be moxide, mcp ; choose the largest
  real, dimension(0:moxide, 0:moxide) :: cp_value_mat, cp_temp_mat, &
      & cond_value_mat, cond_temp_mat, rho_value_mat, rho_temp_mat, &
      & Youngs_modul_mat, Youngs_temp_mat, &
      & th_exp_coeff_mat, th_exp_temp_mat

  ! data read at input for each material or compound; 
  ! data used to obtain properties for each layer based on averaging
  ! number of intervals for the thermal expansion, youngs, cond, cp
  ! number of data points = number of intervals + 1
  integer, dimension(0:moxide) :: nth_exp_mat, nyoungs_mat, &
      & ncond_mat, ncp_prop_mat

  ! when averaging layer properties based on individual compounds
  ! the data needs to have the same temperature intervals

  !    flaw_oxide_size_factor * scale thickness = gives the flaw size of within oxide
  ! flaw_defect_factor = factor dependent on shape, size, and position of the void
  real :: flaw_defect_factor, flaw_oxide_size_factor

  ! initial size of flaws
  real, dimension(10)  :: flaw_size_init

  ! current size of flaws
  real, dimension(10)  :: flaw_size_now
  
  ! type of flaws (porosity or crack); 
  character(LEN = 80), dimension(10):: flaw_name

  ! number of flaws
  integer  :: nflaws

  ! id for flaw type
  integer  :: flaw_id_burried_crack, flaw_id_pore, flaw_id_interf_crack, &
            & flaw_id_surf_crack, flaw_id_delamin

  ! actual number of damage types
  integer :: ndamage  

  ! eta_crit = const_damage_map / sqrt(d_ox)
  real, dimension(mdamage, 0:moxide)    :: const_damage_map

  ! ==============

  ! real, dimension(moxide_layer, mtime) :: 

  ! the name of the specie from the CEA-NASA thermodynamic database
  character(LEN = 80), dimension(0:mspecie) ::  specie_thermo_name

  ! the air_fr_st and gas_comb_fr are given for O2 comubstion; but then 
  ! easly corrected for dry air

  ! actual number of fuel compounds, actual number of combustion compounds
  integer :: nfuel, ncomb, ndiluent

  ! actual number of species introduced
  integer ::   no_thermo_species

  ! data for the oxidant; only one oxidant; only one name list
  ! the oxidant is fed through the compressor and 
  ! the oxidant species are given at input

  character(LEN = 80)  :: oxidant_name
  ! when oxidant_thermo_database_id > 0; we have standard air
  ! the transport and thermodynamic properties are not recomputed
  ! the oxidant species id is used to get properties only for combustion
  ! 
  integer              :: oxidant_thermo_database_id, &
          & oxidant_transport_database_id

  ! actual number of species in the oxidant; must be small than moxidant
  integer              :: noxidant

  character(LEN = 80), dimension(1:moxidant)   :: oxidant_specie_name
  real,                dimension(1:moxidant)   :: oxidant_mole_fraction
  integer,            dimension(1:moxidant)    :: oxidant_sp_th_data_id
  ! actual number of moles per each species of oxidant
  real,                dimension(1:moxidant)   :: oxidant_nmoles

  real                 :: mw_oxidant
  real                 :: hf0_oxidant
  real                 :: dh_oxidant_t0_zero
  real, dimension(mtint, mcp) ::   cp_oxidant_const
  real, dimension(mtint) ::   h1_oxidant_const, s1_oxidant_const
  real, dimension(mtref) ::   Tref_oxidant
  real, dimension(mcp)   ::   cp_oxidant_exp
  ! these are internal variables for the oxidant; most of the time it is air
  real                        :: hf0_oxidant_adj
  real, dimension(mtint, mcp) :: h_oxidant_const
  real, dimension(mcp)        :: h_oxidant_exp
  real, dimension(mtint, mcp) :: hp_oxidant_const
  real, dimension(mcp)        :: hp_oxidant_exp
  real, dimension(ms)         :: s_oxidant_exp
  real, dimension(mtint, ms)  :: s_oxidant_const

  ! dry_air_fuel_id = id for dry air as reactant (in fuel category)
  ! dry_air_comb_id = id for excess dry air in combustion products
  ! n2_comb_id = id of N2 in combustion products
  integer :: dry_air_fuel_id, dry_air_comb_id, n2_comb_id

  ! units of fuel_fr, could be 'mass', 'volumetric', or 'molar', 
  ! default is 'mass'
  ! they have to be consistent; or convert everything to use the 
  ! molar fractions
  character(LEN = 80), dimension(0:mfuel)  ::  fuel_fr_units
  character(LEN = 80), dimension(0:mdiluent)  ::  diluent_fr_units

  ! input file prepared for combustion using NASA-CEA and CANTERA
  character(LEN = string_len) :: cea_combustion_input_file, &
        & cantera_combustion_input_file

  ! relationship between the index and species it points to is in the 
  !  cea_database_id.txt input file
  !  DO NOT CHANGE THIS FILE UNLESS NEW MATERIALS ARE ADDED TO IT
  ! 
  !  fuel id in the thermo database provided by CEA-NASA software
  integer, dimension(0:mfuel)         :: fuel_thermo_database_id

  !  fuel id in the transport database provided by CEA-NASA software
  integer, dimension(0:mfuel)         :: fuel_transport_database_id

  !  diluent id in the thermo database provided by CEA-NASA software
  integer, dimension(0:mdiluent)         :: diluent_thermo_database_id

  !  diluent id in the transport database provided by CEA-NASA software
  integer, dimension(0:mdiluent)         :: diluent_transport_database_id

  !  comb id in the thermo database provided by CEA-NASA software
  integer, dimension(0:mcomb)         :: comb_thermo_database_id

  !  comb id in the transport database provided by CEA-NASA software
  integer, dimension(0:mcomb)         :: comb_transport_database_id

  ! fuel parameters
       ! fuel fraction - fuel_fr
       ! fuel molecular weight - mw_fuel
       ! nmoles_fuel - total number of moles of each fuel
       ! Temp_fuel_inlet - temperature of the fuel as it enters the 
                         ! combustion chamber
  real, dimension(0:mfuel) :: fuel_fr, mw_fuel, Temp_fuel_inlet, nmoles_fuel
  
  ! the MW of fuel (as averaged); equivalent numbers of moles 
  ! in the entire fuel
  real                    :: mw_fuel_ave, nmoles_fuel_ave, nmoles_comb_ave

  ! the MW of combustion products
  real                    :: mw_comb_ave

  ! diluent parameters
       ! diluent fraction - diluent_fr
       ! diluent molecular weight - mw_diluent
       ! Temp_diluent_inlet - temperature of the diluent as it enters the 
                         ! combustion chamber
  real, dimension(0:mdiluent) :: diluent_fr, mw_diluent, Temp_diluent_inlet

  real, dimension(0:mspecie) :: mw_specie

  ! air_to_fuel ratio for total combustion in air for each fuel
  real, dimension(0:mfuel) :: air_fr_st

  ! number of moles of oxygen, O2, required in the reactants 
  ! for total combustion in pure oxygen for each fuel
  real, dimension(0:mfuel) :: oxygen_mole_st

  ! total number of moles of oxygen, O2, for stoichiometric combustion
  real                     :: oxygen_mole_st_ave

  ! total number of moles of air, for stoichiometric combustion
  real                     :: air_mole_st_ave
  
  ! hf0_fuel - is the enthalpy of formation for each fuel
  real, dimension(0:mfuel) :: hf0_fuel
  real, dimension(0:mdiluent) :: hf0_diluent
  real, dimension(0:mspecie) :: hf0_specie

  ! hf0_fuel_adj - is the "hf0_fuel - enthalpy_at_298.15K"  used to save time
  real, dimension(mfuel) :: hf0_fuel_adj
  real, dimension(mdiluent) :: hf0_diluent_adj

       ! dh_fuel_t0_zero - "enthalpy_at_RT - enthalpy_at_0K" 
  real, dimension(0:mfuel) :: dh_fuel_t0_zero
  real, dimension(0:mdiluent) :: dh_diluent_t0_zero
  real, dimension(0:mspecie) :: dh_specie_t0_zero

  ! parameters for CEA-NASA model
  ! temperature points for which each property is given
  real, dimension(0:mfuel, mtref) :: Tref_fuel
  real, dimension(0:mdiluent, mtref) :: Tref_diluent

  real, dimension(0:mspecie, mtref) :: Tref_specie_thermo

  ! additional constans in the CEA-NASA model for enthalpy and entropy
  real, dimension(0:mfuel, mtint) :: h1_fuel_const, s1_fuel_const
  real, dimension(0:mdiluent, mtint) :: h1_diluent_const, s1_diluent_const
  real, dimension(0:mspecie, mtint) :: h1_specie_const, s1_specie_const

  ! exponents for the CEA-NASA model for the specific heat
  real, dimension(0:mfuel, mcp) :: cp_fuel_exp
  real, dimension(0:mdiluent, mcp) :: cp_diluent_exp
  real, dimension(0:mspecie, mcp) :: cp_specie_exp

  ! exponents for the CEA-NASA model for the Dh/DT to be used in jacobian
  ! calculations; not used now since they are identical to that for specific heat
  real, dimension(0:mfuel, mcp) :: hp_fuel_exp
  real, dimension(0:mdiluent, mcp) :: hp_diluent_exp

  ! exponents for the CEA-NASA model for the enthalpy
  real, dimension(0:mfuel, mh) :: h_fuel_exp
  real, dimension(0:mdiluent, mh) :: h_diluent_exp

  ! exponents for the CEA-NASA model for the entropy
  real, dimension(0:mfuel, ms) :: s_fuel_exp
  real, dimension(0:mdiluent, ms) :: s_diluent_exp

  ! constants for the CEA-NASA model for the specific heat
  real, dimension(0:mfuel, mtint, mcp) :: cp_fuel_const
  real, dimension(0:mdiluent, mtint, mcp) :: cp_diluent_const
  real, dimension(0:mspecie, mtint, mcp) :: cp_specie_const

  ! constants for the CEA-NASA model for the Dh/DT to be used in jacobian
  ! calculations; not used now since they are identical to that for specific heat
  real, dimension(0:mfuel, mtint, mcp) :: hp_fuel_const
  real, dimension(0:mdiluent, mtint, mcp) :: hp_diluent_const

  ! constants for the CEA-NASA model for the enthalpy
  real, dimension(0:mfuel, mtint, mh) :: h_fuel_const
  real, dimension(0:mdiluent, mtint, mh) :: h_diluent_const

  ! constants for the CEA-NASA model for the entropy
  real, dimension(0:mfuel, mtint, ms) :: s_fuel_const
  real, dimension(0:mdiluent, mtint, ms) :: s_diluent_const

  ! number of moles for each combustion products that resulted for the
  ! stoichiometric combustion in pure oxygen
  ! per each reactant, i.e., fuel, 
  real, dimension(0:mfuel, mcomb) :: mole_comb_prod_o2_st
  
  ! name of the combustion product for which the number of moles was given
  !   per each fuel
  character(LEN = 80), dimension(0:mfuel, mcomb) ::  name_comb_prod_o2_st

  ! name of combustion products such that they can be identified
  character(LEN = 80), dimension(mcomb) ::  combustion_products_name

  ! per each reactant, i.e., fuel, 
  ! the fraction of combustion products of each species for the combustion 
  ! in AIR
  real, dimension(0:mfuel, mcomb) :: gas_comb_fr_air
  ! the fraction of combustion products of each species for the combustion 
  ! in O2
  real, dimension(0:mfuel, mcomb) :: gas_comb_fr_oxygen

  ! the same variables defined now for the combustion side
  real, dimension(0:mcomb)  :: mw_comb, hf0_comb, dh_comb_t0_zero, comb_nmoles
  real, dimension(0:mcomb)  :: hf0_comb_adj

  real, dimension(0:mcomb, mtref) :: Tref_comb
  real, dimension(0:mcomb, mtint) :: h1_comb_const, s1_comb_const

  real, dimension(0:mcomb, mcp) :: cp_comb_exp, hp_comb_exp
  real, dimension(0:mcomb, mh) :: h_comb_exp
  real, dimension(0:mcomb, ms) :: s_comb_exp
  ! mtint gives the temperature range domain
  real, dimension(0:mcomb, mtint, mcp) :: cp_comb_const, hp_comb_const
  ! h is for enthalpy
  real, dimension(0:mcomb, mtint, mh) :: h_comb_const
  ! s is for entropy
  real, dimension(0:mcomb, mtint, ms) :: s_comb_const

  ! fuel_air_ratio = fuel-air ratio
  real  :: fuel_air_ratio

  ! units of fuel_air_ratio, could be 'mass', 'volumetric', or 'molar', default is 'mass'
  character(LEN = 80)  :: fuel_air_ratio_units
  
  ! far_st =  means stoichiometric fuel-air ratio
  ! air_flow_rate = the flow rate through compressor
  ! gas_flow_rate = fuel_flow_rate + air_flow_rate, hot gas flow rate through turbine
  ! Temp3ad_init = initial guess for the adiabatic flame, or gas temperature 
  ! dTemp_eps = acceptable tolerance for the adiabatic flame temperature computations
  ! Temp3ad  = adiabatic flame temperature
  !
  real  :: far_st, air_st, fuel_flow_rate, diluent_flow_rate, &
           & gas_flow_rate, dTemp_eps, Temp_comb_init, Temp3ad
  real  :: pressure_fuel_inlet

  real  :: combustor_efficiency, pressure_drop_fraction

  ! amount of air required from compressor is given fuel-air ratio
  ! fuel_fr will contain the fraction of the each fuel, i.e. n_fuel/total_n_fuel

  ! one calculation; impose the adiabatic temperature and get the fuel to air ratio
  ! or, impose the fuel-air-ratio and get the combustion temperature

  ! the function id which has the natural log multiplying the power law fuction
  ! max_iterat_solver = the maximum number of iterations
  integer :: cp_log_fuel, cp_log_comb, h_log_fuel, h_log_comb, &
           & hp_log_fuel, hp_log_comb, s_log_fuel, s_log_comb, &
           & cp_log_diluent, h_log_diluent, &
           & hp_log_diluent, s_log_diluent, &
           & cp_log_oxidant, h_log_oxidant, &
           & hp_log_oxidant, s_log_oxidant

  ! fuel_mask is true if not air
  logical, dimension(0:mfuel)  :: Fuel_Mask

  ! if dilutants are considered
  logical                      :: if_diluent

  ! if coolant bleeding from compressor is considered
  logical                      :: if_coolant_compressor

  ! if flow rate is a variable in the model are considered
  logical                      :: if_fuel_flow_rate

END MODULE OXIDE_DATA_MODULE
