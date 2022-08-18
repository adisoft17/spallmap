MODULE THERMO_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all combustion data structures integer and scalar parameters.
  !
  ! Author: Adrian S. Sabau, sabaua@ornl.gov
  !
  !=======================================================================
  use parameter_module, only: mfuel, mcomb, mcp, mtref, mtint, ms, mh, &
                      & mspecie, mtrans, mgas, mstates, mcool, &
                      & mchord_cool, mspan_cool, mggt
  implicit none

  save
  ! Public Module
  public

  ! true if static calculations
  ! false if stagnation properties are used for the thermodynamic model
  logical              :: if_static_or_stagnation

   ! air_id; the id of air mixture in the thermodynamic cycle
   integer, parameter   :: air_id = 1, gas_id = 2, coolant_id = 3

  ! each gas mixture (air, hot gas, and after cooling air+hot gas)
  ! has a state and some properties; the fraction of air/hot gas 
  ! will change as the coolant mixes with hot gas at different blade stages
  
  type ID_STATE
    ! number of moles; use this as moles per unit second
    integer   :: moles  
    ! molecular weight
    integer   :: mw
    ! mass 
    integer   :: mass
    ! static
    integer   :: Temp
    ! pressure
    integer   :: press
    ! total enthalpy
    integer     :: h
    ! molar enthalpy
    integer     :: h_m
    ! total entropy
    integer     :: s
    ! molar entropy
    integer     :: s_m
    ! specific heat; molar basis; for the air + hot gasses; too complex
    ! integer     :: cp_m
    ! stagnation temperature
    integer     :: Temp_st
    ! stagnation pressure
    integer     :: press_st
    ! stagnation enthalpy
    integer     :: h_st
    ! specific molar stagnation enthalpy
    integer     :: h_st_m
    ! stagnation entropy (only the thermal component)
    ! there is one entropy made up of 
    ! thermal component, pressure component, etc
    integer    :: s_st
    ! specific molar stagnation entropy
    integer    :: s_st_m
  end type ID_STATE

  ! inlet fuel
  type(ID_STATE), dimension (1)                       :: id_fuel
  ! inlet dilutant
  type(ID_STATE), dimension (1)                       :: id_dil
  ! inlet compressor is 1, intermidiate and exit
  type(ID_STATE), dimension (10)                      :: id_air_compr
  ! inlet into combustor
  type(ID_STATE), dimension (1)                       :: id_air_comb

  ! coolant extraction from compressor
  ! this would be divided later into the stator only, rotor only
  ! or stator and rotor
  type(ID_STATE), dimension (mcool)                   :: id_air_bleed
  
  ! the air_bleed could split into two; to the rotor and stator
  type(ID_STATE), dimension (2 * mcool)                   :: id_cool

  ! coolant from compressor to the rotor turbine blades
  type(ID_STATE), dimension (mcool)                   :: id_air_rot
  ! coolant from compressor to the stator turbine vanes
  type(ID_STATE), dimension (mcool)                   :: id_air_van

  ! everywhere air 
  type(ID_STATE), dimension (mstates)                 :: id_air

  type ID_CP1
    ! specific heat on a molar basis
    integer :: m_air
    integer :: m_gas

    ! specific heat on a mass basis
    integer :: air
    integer :: gas

  end type ID_CP1

  type(ID_CP1), dimension (mstates)              :: id_cp

  ! everywhere hot gas
  type(ID_STATE), dimension (mstates)                 :: id_gas
  ! all the time store the temperature and pressure into the id_air for 
  ! the id_gas in order not to introduce more unknownws
  ! mixture of everywhere upstream combustion
  ! type(ID_STATE), dimension (mstates)                 :: id_gas

  type VEL1
    ! absolute velocity; at inlet and outlet of the blade or vane
    real   :: abs
    ! relative velocity
    real   :: rel
    ! peripheral velocity component due solely to blade rotation
    ! peripheral velocity; blade velocity at the reference radius;
    ! omega * radius; use terminology of rotational velocity
    real   :: rot
    ! axial velocity
    real   :: abs_axial
    real   :: rel_axial
    ! circumferential or tangential velocity
    real   :: abs_circ
    real   :: rel_circ
  end type VEL1

  type(VEL1), dimension(10)  :: compr_vel, comb_vel

  ! actual number of points on the blade in the chord direction
  ! from the leading edge towards the trailing edge
  ! to be used for expansion and film cooling
  integer   :: nchord
  ! actual number of points on the blade in the span direction
  ! from the root towards the blade tip to be used for blade cooling analysis
  integer   :: nspan

  ! actual rows cooled internally, using TBC, and Film cooling techniques
  integer   :: nint_cool_row, ntbc_cool_row, nfc_cool_row

  ! actual number would be ggt_no_stage
  ! type(ID_VEL2), dimension(2, mggt)   :: id_gas_vel

  type(VEL1), dimension(2, mggt, mchord_cool, mspan_cool)  :: vel_gas1

  ! velocity of the gas in the turbine
  ! inlet vanes (1, stage); outlet vanes (2, stage)=inlet blades (3, nstage)
  ! ; and outlet blades(3, stage)
  type(VEL1), dimension(3, mggt)  :: gas_vel

  real, dimension(mggt)          :: degree_reaction_stage, loading_factor_stage

  ! use this for different expansion steps
  ! one expansion steps is equivalent to average values
  type ID_STAGN_STATE_TURB
    ! mean radius across the span of the blade
    ! 1 - at the root of the blade
    ! nspan_cool - at the blade tip
    integer, dimension(mspan_cool)   :: radius
    ! velocity
    integer, dimension(mchord_cool, mspan_cool)   :: abs_vel
    ! swirl velocity
    integer, dimension(mchord_cool, mspan_cool)   :: swirl_vel
    ! peripheral velocity; omega * radius; use terminology of rotational velocity
    integer, dimension(mspan_cool)   :: rot_vel
    ! relative velocity
    integer, dimension(mchord_cool, mspan_cool)   :: rel
    ! axial velocity
    integer, dimension(mchord_cool, mspan_cool)   :: axial
    ! chord position
    integer, dimension(mchord_cool)    :: axial_chord
    ! stagnation temperature
    integer , dimension(mchord_cool, mspan_cool)  :: Temp_st
    ! stagnation pressure
    integer, dimension(mchord_cool, mspan_cool)   :: press_st
    ! stagnation enthalpy
    integer, dimension(mchord_cool, mspan_cool)     :: h_st
    ! stagnation entropy
    integer, dimension(mchord_cool, mspan_cool)    :: s_st
  end type ID_STAGN_STATE_TURB


  type ID_PROP
    !  density
    integer     :: rho
    !  specific heat
    integer     :: cp
    !  viscosity
    integer     :: Visc
    !  thermal conductivity
    integer     :: Cond
    !  Prandtl Number
    integer     :: Prandtl
    !  Nusselt Number
    integer     :: Nusselt
    !  Peclet number
    integer     :: Peclet 
    ! Reynolds number
    integer     :: Reynolds
    ! heat transfer coefficient
    integer     :: htc
  end type ID_PROP

  ! coolant from compressor to the rotor turbine blades
  type(ID_PROP), dimension (mcool)                   :: cool_prop_rot
  ! coolant from compressor to the stator turbine vanes
  type(ID_PROP), dimension (mcool)                   :: cool_prop_van
  ! hot gases
  type(ID_PROP), dimension (mcool)                   :: gas_prop  

  type ID_TEMP_PROFILE
     ! static mainstream gas recovery temperature
     integer   :: GR
     ! stagnation mainstream gas recovery temperature
     integer   :: GR_st
     ! stagnation mainstream gas recovery temperature for rotors (based on relative temperature)
     integer   :: GR_rel
     ! mean adiabadic wall temperature; this is due to film cooling
     integer   :: aw
     ! temperature of the TBC surface at the inlet
     integer   :: TBC_in
     ! temperature of the TBC surface at the exit of coolant (maximum temp)
     integer   :: TBC_out
     ! temperature of the OUTER surface of the blade at the coolant inlet, downstream 
     integer   :: Metal_TBC_in
     ! temperature of the OUTER surface of the blade at the coolant OUTlet, upstream 
     integer   :: Metal_TBC_out
     ! temperature of the INNER surface of the blade at the coolant inlet, downstream 
     integer   :: Metal_Cool_in
     ! temperature of the INNER surface of the blade at the coolant OUTlet, upstream 
     integer   :: Metal_Cool_out
     ! temperature of the coolant at the inlet
     integer   :: Cool_in
     ! temperature of the coolant at the outlet
     integer   :: Cool_out
     ! enthalpy of the coolant at inlet
     integer   :: h_cool_in
     ! enthalpy of the coolant at inlet
     integer   :: h_cool_out
     ! heat flux through the blade; YW02
     integer   :: flux
  end type ID_TEMP_PROFILE

  ! temperature profile across the thickness of the rotor
  type(ID_TEMP_PROFILE), dimension (mstates)                   :: id_temp

  ! temperature profile across the thickness of the stator
  ! type(ID_TEMP_PROFILE), dimension (mstates)                   :: id_temp_van

  type ID_TH_CYCLE
    ! 
    integer   :: work  
    !
    integer   :: heat
    ! 
    integer   :: exergy
    ! 
    integer   :: irrever
    ! 
  end type ID_TH_CYCLE

  type(ID_TH_CYCLE), dimension (mstates)                 :: id_eff

  ! number of actual components in the turbine for which the work or heat is defined
  integer     ::  no_work_components

  ! index for the total compression work; and combustion component  
  integer     ::  no_work_compressor, no_combustion_heat, no_work_pt

  ! index for cycle efficiency of the entire thermodynamic cycle
  integer     :: cycle_effiency_id

  ! index for the total turbine work; this is used when the
  ! power turbine is not considered
  ! eta_turbine = no_work_useful / no_cumbustion_heat
  ! no_work_last_stage = work in the last stage
  ! it is the power turbine, if there is one,
  ! or it is the last stage of the turbine
  integer     :: no_work_turbine_total, no_work_ggt, no_work_useful, no_work_last_stage

  ! index for the work performed per each stage
  integer, dimension(mggt)     :: no_work_stage

  ! define property structure
  type MASS_STATE 
    ! molecular weight
    real     :: mw
    ! mass; should be mass flow rate if the time frame is 1sec; uniform flow
    real     :: mass
    ! fraction (molar = volumetric)
    real     :: fraction
    !  velocity
    real     :: Velocity
  end type MASS_STATE

  type TH_STATIC_STATE 
    !  static temperature
    real     :: Temp
    ! static pressure
    real     :: pressure
  end type TH_STATIC_STATE

  type TH_STAGNATION_STATE 
    !  total (stagnation) temperature
    real     :: Temp
    ! total (stagnation) pressure
    real     :: pressure
  end type TH_STAGNATION_STATE

  type TH_STATIC_PROP
    ! enthalpy
    real     :: enthalpy
    ! entropy
    real     :: entropy
  end type TH_STATIC_PROP

  type TH_STAGNATION_PROP
    ! total (stagnation) enthalpy
    real     :: enthalpy
    ! total (stagnation) entropy
    real     :: entropy
  end type TH_STAGNATION_PROP

  type TR_STATIC_PROP
    !  viscosity
    real     :: Visc
    !  thermal conductivity
    real     :: Cond
  end type TR_STATIC_PROP

  type TR_STAGNATION_PROP
    !  total (stagnation) viscosity
    real     :: Visc
    !  total (stagnation) thermal conductivity
    real     :: Cond
  end type TR_STAGNATION_PROP

  type REL_STATIC_STATE
    !  Peclet number
    real     :: Peclet 
    ! Reynolds number
    real     :: Reynolds
    ! heat transfer coefficient
    real     :: htc
  end type REL_STATIC_STATE

  type REL_STAGNATION_STATE
    !  total (stagnation) Peclet number
    real     :: Peclet
    ! total (stagnation) Reynolds number
    real     :: Reynolds
    ! heat transfer coefficient
    real     :: htc
  end type REL_STAGNATION_STATE

   ! this should be valid after combustion or 
   ! for average properties of mixtures

   ! type(MASS_STATE), dimension (mgas, mstates)           :: flow
   ! type(TH_STATIC_STATE), dimension (mgas, mstates)      :: state
   ! type(TH_STAGNATION_STATE), dimension (mgas, mstates)  :: state_tot
   ! type(TH_STATIC_PROP),  dimension (mgas, mstates)      :: prop_th
   ! type(TH_STAGNATION_PROP),  dimension (mgas, mstates)  :: prop_th_tot
   ! type(TR_STATIC_PROP), dimension (mgas, mstates)       :: prop_tr
   ! type(TR_STAGNATION_PROP), dimension (mgas, mstates)   :: prop_tr_tot
   ! type(REL_STATIC_STATE), dimension (mstates)           :: relation
   ! type(REL_STAGNATION_STATE), dimension (mstates)       :: relation_tot

   ! 4 variables for flow, 2 variables for state, prop_th, and prop_tr
   integer, dimension (4), parameter             :: id_stat = (/4, 2, 2, 2/)

   ! beginning index for the static case
   integer, dimension (4), parameter       :: id_stat_start = (/0, 4, 6, 8/)

   ! ending index for the static case
   integer, dimension (4), parameter        :: id_stat_end = (/4, 6, 8, 10/)

   ! these indices would work if the variable arrangement would be 
   ! air variables then gas variables
   ! indices for the MASS_STATE
   integer     :: mw_id = id_stat_start(1) + 2,  &
                & mass_id = id_stat_start(1) + 1, &
                & fr_id = id_stat_start(1) + 3, vel_id = id_stat_start(1) + 4

   ! indices for the TH_STATE (static or stagnation)
   integer     :: temp_id = id_stat_start(2) + 1, &
                & press_id = id_stat_start(2) + 2

   ! indices for the TH_PROP (static or stagnation)
   integer     :: enth_id = id_stat_start(3) + 1, &
                & entro_id = id_stat_start(3) + 2

   ! indices for the TR_PROP (static or stagnation)
   integer     :: visc_id = id_stat_start(4) + 1, &
                & cond_id = id_stat_start(4) + 2

   ! these indices would work if the variable arrangement would be 
   ! flow variables (air, gas, coolant), state variables (air, gas, coolant)
   ! prop_th variables (air, gas, coolant), 
   ! prop_tr variables (air, gas, coolant)
   ! this arrangement would yield better matrix

   ! the order for unknown numbering would be:
   ! fuel (flow, state, th), dilutant (flow, state, th, tr), 
   ! air(flow, state, th, tr), gas (flow, state, th, tr), 
   ! coolant (flow, state, th, tr) - coolant is considered only when there is 
   ! a coolant flow out or coolant flow in

   ! states definition
   integer     ::  compressor_inlet = 1
   integer     ::     compressor_isentropic = 2

   integer     ::     fuel_solver_state = 1
   integer     ::     dilutant_solver_state = 2
   integer     ::     air_solver_inlet = 3
  
   
END MODULE THERMO_DATA_MODULE
