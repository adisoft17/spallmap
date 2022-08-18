MODULE BOILER_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all boiler data structures integer and scalar parameters.
  !
  !  Author: Adrian S. Sabau, sabaua@ornl.gov
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================
  use parameter_module, only: npulse
 
  ! mfuel, mcomb, mcp, mtref, mtint, &
  !            & mcompr, mggt, mcool, max_mix, mrotor_expansion

  implicit none
  save
  ! Public Module
  public

  ! pr in handbook; pr = exp(s/R); relative pressure

  ! tube geometry
  real :: tube_outer_radius, tube_inner_radius, tube_curvature, tube_thickness

  ! htc at full load
  real  :: htc_tube_inner, htc_tube_outer

  ! htc at partial load
  real  :: htc_tube_inner_partial, htc_tube_outer_partial

  ! temperature schedule; use this one for simple cases 
  ! by setting no_total_pulse = 0
  real, dimension(100) :: Time_boiler, Temp_steam

  ! deal with pressure and flow rate later
  real  :: gas_flow_rate

  ! data at each point interval; internal to the code
  ! time units are by default hours
  real, dimension(:), pointer   :: Time_Boiler_p, TEMP_STEAM_p, &
       & TEMP_STEAM_Slope_p, Time_Boiler_sec_p, oxide_thickness, &
       & Temp_gas_p, Temp_gas_slope_p
  real, dimension(:, :), pointer   :: thickness_fr_var

  ! thickness_fr_var - when oxide layers grow at different temperature
  ! thickness_fr need to be updated at each time period

  ! id-s for the segments begin (ramp-up, pulse, ramp-down, idle)
  !                        end  (idle, ramp-up, pulse, ramp-down)
  integer, dimension(:), pointer   :: id_temp_op, id_press_op, id_flow_op

  ! pressure in the tubes
  ! supercritical - new boilers both the presssure and temperature are varied
  ! subcritical - older boilers only the temperature is cycled
  real, dimension(:), pointer   :: press_inner_p, press_outer_p, &
       & press_inner_slope_p, press_outer_slope_p

  ! use this data for multiple pulse,
  integer :: no_total_pulse, TEMP_BOILER_data_no

  real, dimension(npulse)  :: temp_steam_pulse, press_inner_pulse, &
          & press_outer_pulse, temp_gas_pulse
  real, dimension(2)  :: time_boiler_idle, time_boiler_pulse, &
          & temp_steam_idle, temp_gas_idle
  real  :: time_pulse_idle, &
          & temp_steam_first_idle, temp_gas_first_idle, &
          & time_boiler_first_idle, &
          & temp_steam_last_idle, temp_gas_last_idle, &
          & time_boiler_last_idle, time_idle_pulse

  ! data for events; an event is bringing the system to a 
  ! a lower temperature than the idle
  real  ::    event_time_interval, time_boiler_event, &
          & temp_steam_event, temp_gas_event, press_inner_event, &
          & press_outer_event, time_pulse_event, time_event_pulse

  ! data on pressure
  real  :: press_inner_first_idle, press_outer_first_idle, &
         & press_inner_last_idle, press_outer_last_idle
  real, dimension(2)  :: press_inner_idle, press_outer_idle

  ! used for outer oxide scales
  character(LEN = 80)  :: oxide_location

  ! no_pulse_per_cycle = 1 for daily cycle; 2 for weekly cycle
  integer :: no_pulse_per_cycle

  ! fraction of heat flux at idle
  real, dimension(2)   ::  heat_flux_fraction_idle

  ! sliding variables; 
  logical  :: if_sliding_load

  ! transition idle - pulse ; pulse idle using load 
  real, dimension(10)  :: time_pulse_idle_load, load_fraction_pi_time

  ! transition idle - pulse ; pulse idle using load 
  real, dimension(10)  :: time_idle_pulse_load, load_fraction_ip_time

  ! load fraction at idle; this variable may be seldome used since
  ! the load may not vary linear from 1 to idle
  real, dimension(2)   ::  load_fraction_idle

  ! steam load variation
  real, dimension(10)  :: load_temp_steam, temp_steam_vs_load, &
      & press_steam_vs_load

  ! flue gas load variation
  real, dimension(10)  ::   load_temp_gas, temp_gas_vs_load, &
      & press_gas_vs_load, heat_flux_fraction_vs_load

  integer              :: high_load_id_history, low_load_id_history

  ! end_state_boiler = 'partial_load' or end_state = 'outage' to prescribe the 
  ! last end point for the boiler operation
  character(LEN = 80)  :: end_state_boiler

END MODULE BOILER_DATA_MODULE
