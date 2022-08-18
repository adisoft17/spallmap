MODULE SHUTDOWN_DATA_MODULE

  !=======================================================================
  ! Purpose(s):
  !
  !   Define all oxide data structures integer and scalar parameters.
  !
  !   Author: Adrian S. Sabau, sabaua@ornl.gov oxide_data_module.F90
  !   Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================
  implicit none
  save
  ! Public Module
  public

  ! number_outage_periods = number of outage regimes
  integer :: number_outage_periods

  ! time_outage_period = intervals for outage periods; must be number_outage_periods + 1
  ! duration_outage_period = DELTA(time_outage_period)
  real, dimension(5)   :: time_outage_period, duration_outage_period

  ! intervals between the output per each period
  real, dimension(5)   :: delta_time_outage_period

  ! no_out_f2l_period * delta_time_outage_period = duration_outage_period
  integer, dimension(5)    :: no_out_f2l_period

  ! no_out_f2l_total = SUMM(no_out_f2l_period)
  integer   :: no_out_f2l_total  ! number of total outputs for large shutdown

  ! steam temperature at time_outage_period
  real, dimension(5)   :: temp_steam_outage_period

  ! just for output
  ! output times and main variables
  real, dimension(500)  :: time_now_f2l, temp_steam_now_f2l, &
        & temp_gas_now_f2l, press_steam_now_f2l, press_gas_now_f2l

END MODULE SHUTDOWN_DATA_MODULE
