MODULE COOLING_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all combustion data structures integer and scalar parameters.
  !
  ! Author: Adrian S. Sabau, sabaua@ornl.gov
  !
  !=======================================================================
  use parameter_module, only: mfuel, mcomb, mcp, mtref, mtint, ms, mh, &
          & mcool, mgas, mstates, mcool, mggt, mrow_ggt, mcool_problem
  implicit none

  save
  ! Public Module
  public

  ! fluctuation of the combustion temperature
  real, dimension(0:2, 0:mggt)   ::   delta_temp_combust
  
  ! the total gas temperature fluctuates as delta_temp_combust * factor_temp_combust
  real, dimension(0:2, 0:mggt)   ::   factor_temp_combust

   ! transversal area of the internal cooling passages
   real, dimension(0:2, 0:mggt)  :: Area_internal_cool_passage, area_total_internal_cool

   ! perimeter for the internal cooling passage
   real, dimension(0:2, 0:mggt)  :: perimeter_internal_cool_passage 

   ! number of cooling channels
   integer, dimension(0:2, 0:mggt) :: no_cool_channels 

   ! number of passes or turns per each channel
   integer, dimension(0:2, 0:mggt) :: no_passes_per_channel 

   ! height of each pass; this should account for the bent; U-turn like
  ! length or height of coolant pass; each channel could have more than one 
  ! pass or loop that goes up to the blade tip and down back to the blade root
   real, dimension(0:2, 0:mggt)     :: height_cool_internal_pass 
   
  ! factor that allows for heat transfer between channels
  ! appendix A1, TH05
  real, dimension(0:2, 0:mggt)      :: ht_channel_to_channel_factor

   ! Enhancemen factor for ribbs in the internal channels
   real, dimension(0:2, 0:mggt)     :: ribs_enhance_factor_ref, &
         & ribs_enhance_factor 

   ! diameter at which the coolant is injected 
   real, dimension(0:2, 0:mggt)     ::   diameter_coolant_injection

   ! data for the film cooling now
   logical  :: if_film_cooling

   ! slot width for injection of coolant out of the blade
   real, dimension(0:2, 0:mggt)     :: slot_fc_width
   ! angle at which the coolant is injected out of the blade; film cooling
   real, dimension(0:2, 0:mggt)     ::    angle_fc_inject
   ! characteristic length in the chordwise direction for film cooling correlations
   real, dimension(0:2, 0:mggt)     ::    length_chord_fc

  ! thermal conductivity of the blade at high temperatures
  real, dimension(0:2, 0:mggt)          :: metal_conductivity, metal_cp
  ! thermal conductivity of the tbc at high temperatures
  real, dimension(0:2, 0:mggt)          :: tbc_conductivity, tbc_cp
   
  ! maximum temperature allowed in the blade
  real          :: Temp_metal_max

  ! maximum coolant faction extracted for each point in the compressor
  ! should have the dimension of the maximum coolant extraction points
  ! or maximum stages in the compressor
  real, dimension(mcool)    :: coolant_fr_compr_extr_max, coolant_fr_compr_extr_init
 
  ! maximum rotor coolant faction extracted from the each stage coolant
  real, dimension(mcool)    :: coolant_fr_rotor_extr_max, coolant_fr_rotor_extr_init

  ! coolant_compr_extr_stage = i specifies the compressor stage after which the
  ! coolant is extracted; coolant is extracted between stages i and i+1
  integer, dimension(mcool) ::   coolant_compr_extr_stage

  ! coolant_destinationspecifies the destination of the coolant extracted from compressor
  ! 's&r' - stator and rotor
  ! 's'   - only stator 
  ! 'r'   - only rotor
  character(LEN = 80), dimension(mcool) ::  coolant_destination

  ! specifies the turbine stage where coolant will be connected
  integer, dimension(mcool) ::  coolant_dest_turbine_stage

  ! maximum film cooling fraction (film_cooling / total_internal_coolant)
  real, dimension(mrow_ggt)          :: coolant_fr_fc_max, coolant_fr_fc_init

  ! store coolant id
  integer, dimension(mrow_ggt)          :: no_cool_row

  ! store the bleed id from the compressor
  integer, dimension(mrow_ggt)          :: bleed_id_row

  ! fraction of coolant treated implicitly
  ! 1 = all is implicit; 0 = all is explicity
  real             ::  implicit_coolant_fraction 
  
  ! explicit = 1 - implicit
  real             ::  explicit_coolant_fraction 

  ! used to store the temperature difference in order to neglect the
  ! cooling
  real, dimension(mrow_ggt, 200)  :: coolant_to_decrease

  real, dimension(mrow_ggt, 5, 200)  :: dtemp_to_decrease

  ! pressure in the stator before mixing of coolant
  ! used to estimate the desired compressor ratios for 
  ! extracting the coolant
  real, dimension(mggt)  :: press_stator_stage, press_compr_stage_need

  ! constraint_Temp_metal = .true. if maximum metal temperature is fixed
  !       here the coolant fractions are unknowns
  ! constraint_cool_fr = .true. if the coolant fractions are fixed
  !       here the metal temperature is unknown

  ! match_press_bleed_compr_turb = true when the
  ! pressure on the intermediate compressor stages are matched to those of
  ! the turbine points where coolant is injected
  logical    :: constraint_Temp_metal, constraint_cool_fr, &
         & match_press_bleed_compr_turb

  logical, dimension(mrow_ggt)   :: constraint_zero_coolant, &
       & constraint_zero_heat, eliminate_film_cool, &
       & eliminate_temp_constraint

  ! if small_coolant_fr_run = .true., then run first case with very small
  !  coolant fraction; this is known to converge
  logical    ::   small_coolant_fr_run

  ! maximum coolant faction extracted for each point in the compressor
  ! data is stored for different problems
  ! the last problem is with coolant_fr_compr_extr_init
  real, dimension(mcool_problem, mcool)    :: coolant_fr_compr_extr_prob

  ! maximum rotor coolant faction extracted from the each stage coolant
  ! data is stored for different problems
  ! the last problem is with _init variable
  real, dimension(mcool_problem, mcool)    :: coolant_fr_rotor_extr_prob

  ! number of problems to be run 
  ! the second problem will continue the first problem
  integer  :: no_problem

  ! if this is true, then determine the turbine constant, CT 
  ! per Jordal rel. 3.11 pp. 31
  logical  :: if_reference_turbine, if_compr_turbine_match 

  !  reference turbine constant
  real  :: turbine_constant_ref

  ! reference molecular weight at the first vane
  real  :: mw_first_vane_ref

END MODULE COOLING_DATA_MODULE
