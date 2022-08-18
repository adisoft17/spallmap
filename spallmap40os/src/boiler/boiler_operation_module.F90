MODULE BOILER_OPERATION_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Read  namelists. 
    !   sabaua@ornl.gov read_input_module.F90
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: BOILER_SCHEDULE_TP_CYCLES

 CONTAINS

  SUBROUTINE BOILER_SCHEDULE_TP_CYCLES(if_cycle)
  
    use parameter_module,   only: mrow_ggt, pi, zero, mramp, mramp_out
    use input_utilities_module,   only: SEARCH_NML
    use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
    use boiler_data_module, only: tube_thickness, tube_outer_radius, &
          & tube_curvature, &
          & tube_inner_radius, htc_tube_inner, htc_tube_outer, &
          & Time_boiler, Temp_steam, Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, & 
          & no_total_pulse, time_boiler_idle, time_boiler_pulse, &
          & time_pulse_idle, temp_steam_pulse, temp_steam_idle, &
          & temp_steam_first_idle, time_boiler_first_idle, &
          & temp_steam_last_idle, time_boiler_last_idle, time_idle_pulse
    use boiler_data_module, only: temp_gas_pulse, temp_gas_idle, &
          & temp_gas_first_idle, &
          & temp_gas_last_idle, temp_gas_slope_p, temp_gas_p, &
          & gas_flow_rate, oxide_location, no_pulse_per_cycle, &
          & press_inner_p, press_outer_p, press_inner_pulse, &
          & press_outer_pulse, press_inner_first_idle, &
          & press_outer_first_idle, &
          & press_inner_idle, press_outer_idle, &
          & press_inner_last_idle, press_outer_last_idle, &
          & press_inner_slope_p, press_outer_slope_p
    use boiler_data_module, only:    event_time_interval, time_boiler_event, &
          & temp_steam_event, temp_gas_event, press_inner_event, &
          & press_outer_event, time_pulse_event, time_event_pulse
    use boiler_data_module, only:   time_pulse_idle_load, &
      & load_fraction_pi_time, time_idle_pulse_load, &
      & load_fraction_ip_time, load_fraction_idle, &
      & load_temp_steam, temp_steam_vs_load, press_steam_vs_load, &
      & load_temp_gas, temp_gas_vs_load, &
      & press_gas_vs_load, heat_flux_fraction_vs_load
    use boiler_data_module, only: id_temp_op, id_press_op, id_flow_op
    use solution_data_module, only: no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low, time_full2low, &
      & time_full2low, temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low, &
      & no_f2l_event_out, no_event_point_full2low, &
      & no_event_full2low_output

    ! Argument List
    logical, intent(IN) :: if_cycle

    ! Local Variables
    integer :: ioerror, i, k, j, j1, &
        & no_pct_first_pulse = 2, no_pct_pulse = 3, &
        & no_pct_idle = 3, no_pct_joint = 1, no_pct_last_pulse = 3, Status
    integer :: event_local_point
    integer, save :: no_event
    integer, save, dimension(200) :: event_start

    logical :: no_compr_namelist, compr_namelist
    real    :: Temp_ave, total_time, Temp_ave_local, Time_event_now
    
        if (no_pulse_per_cycle == 1)  then
          no_pct_first_pulse = 2
          no_pct_pulse = 3
          no_pct_idle = 3
          no_pct_joint = 1 
          no_pct_last_pulse = 3
        else if (no_pulse_per_cycle == 2)  then 
          no_pct_first_pulse = 2
          no_pct_pulse = 3
          no_pct_idle = 3
          no_pct_joint = 1 
          no_pct_last_pulse = 3
        else 
          write(6, *) 'STOP: no_pulse_per_cycle > 2 or < 1'
          STOP
        endif

        ! roughly define no_event_full2low_output rather than assigning it
        ! first find total time
        if (no_event_full2low_output(1) < 0)  then
          ! this is the default for simplified input
          total_time = time_boiler_first_idle + time_boiler_last_idle + &
            & no_total_pulse * (time_boiler_idle(1) + time_idle_pulse + &
            & time_boiler_pulse(1) + time_pulse_idle)

          ! no_event_full2low_output = 1, 2, 3, 4, 5, 6, 4*0, for
          ! daily cycle and event_time_interval = 7920.0,
          do i = 1, 10
            if (event_time_interval * i < total_time)  then
              no_event_full2low_output(i) = i
            else
              no_event_full2low_output(i) = 0
            endif
          enddo

        endif

        write(out_lun, 1)  no_event_full2low_output
 1      format('no_event_full2low_output = ', 10(1x, i2))

        ! no_pulse_per_cycle is considered
        TEMP_BOILER_data_no = no_pulse_per_cycle * &
          & (no_total_pulse * no_pct_pulse + &
          & (no_total_pulse -1) * no_pct_idle - &
          & (no_total_pulse -1) * 2 * no_pct_joint + 1) + &
          & no_pct_first_pulse - no_pct_joint + &
          & no_pct_last_pulse - no_pct_joint -1

        write(out_lun,*) 'allocate boiler cycle ', TEMP_BOILER_data_no
        ALLOCATE(TEMP_STEAM_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(TEMP_GAS_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(Time_Boiler_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(TEMP_STEAM_Slope_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(TEMP_GAS_Slope_p(1:TEMP_BOILER_data_no), STAT = Status)

        ALLOCATE(id_temp_op(1:TEMP_BOILER_data_no), STAT = Status)

        ! pressure cycle
        ALLOCATE(press_inner_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(press_outer_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(press_inner_slope_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(press_outer_slope_p(1:TEMP_BOILER_data_no), STAT = Status)

        no_full2low_output = 0
        
        do i = mramp_out, 1, -1
          if (no_full2low_output == 0 .and. &
            & no_pulse_full2low_output(i) /= 0) &
            & no_full2low_output = i
        end do

        ! get the event output
        no_f2l_event_out = 0
        do i = 10, 1, -1
          if (no_f2l_event_out == 0 .and. &
            & no_event_full2low_output(i) /= 0) &
            & no_f2l_event_out = i
        end do

        write(out_lun, *) ' no_full2low_output = ', no_full2low_output
        write(out_lun, *) ' no_f2l_event_out = ', no_f2l_event_out
        ! do i = no_full2low_output +1, no_full2low_output + no_f2l_event_out

        ! end do

        ! set the values and times based on number of pulses
        Time_Boiler_p(1) = 0.0
        TEMP_STEAM_p(1) = temp_steam_first_idle
        TEMP_GAS_p(1) = temp_gas_first_idle
        press_inner_p(1) = press_inner_first_idle
        press_outer_p(1) = press_outer_first_idle

        k = 1
        Temp_ave = 0.0
        total_time = 0.0

        ! initialize the time for event now
        Time_event_now = event_time_interval
        event_local_point = 0
        no_event = 0

        do i=1, no_total_pulse

          k = k + 1

          if (k == 2)  then
            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_first_idle
            TEMP_STEAM_p(k) = temp_steam_first_idle
            TEMP_GAS_p(k) = temp_gas_first_idle
            press_inner_p(k) = press_inner_first_idle
            press_outer_p(k) = press_outer_first_idle

            ! id for end of idle and beginning of a new cycle ramp up
            id_temp_op(k) = 1

          else 

            if (event_local_point == 0)  then
              ! no event now
              if (no_pulse_per_cycle == 1)  then
                Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_idle(1)
                TEMP_STEAM_p(k) = temp_steam_idle(1)
                TEMP_GAS_p(k) = temp_gas_idle(1)
                press_inner_p(k) = press_inner_idle(1)
                press_outer_p(k) = press_outer_idle(1)
              else if (no_pulse_per_cycle == 2)  then
                ! the second pulse ends here
                Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_idle(2)
                TEMP_STEAM_p(k) = temp_steam_idle(2)
                TEMP_GAS_p(k) = temp_gas_idle(2)
                press_inner_p(k) = press_inner_idle(2)
                press_outer_p(k) = press_outer_idle(2)
              endif

            else if (event_local_point == 2)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_event
              TEMP_STEAM_p(k) = temp_steam_event
              TEMP_GAS_p(k) = temp_gas_event
              press_inner_p(k) = press_inner_event
              press_outer_p(k) = press_outer_event

              ! increment local event point
              event_local_point = event_local_point + 1

            endif

            ! id for end of idle and beginning of a new cycle ramp up
            id_temp_op(k) = 1

          endif

          k = k + 1

          ! n_pulse = 1 + (k-3)/4, this should be just the pulse number; i

          ! last point in the local event
          if (event_local_point == 0)  then
            ! no event

            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_idle_pulse
            if (no_pulse_per_cycle == 1)  then
              TEMP_STEAM_p(k) = temp_steam_pulse(i)
              TEMP_GAS_p(k) = temp_gas_pulse(i)
              press_inner_p(k) = press_inner_pulse(i)
              press_outer_p(k) = press_outer_pulse(i)
            else if (no_pulse_per_cycle == 2)  then
              ! first pulse
              TEMP_STEAM_p(k) = temp_steam_pulse(1)
              TEMP_GAS_p(k) = temp_gas_pulse(1)
              press_inner_p(k) = press_inner_pulse(1)
              press_outer_p(k) = press_outer_pulse(1)
            endif

            Temp_ave = Temp_ave + 0.5 * time_idle_pulse * &
                       & (TEMP_STEAM_p(k-2) + TEMP_STEAM_p(k-1))

          else if (event_local_point == 3)  then

            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_event_pulse

            if (no_pulse_per_cycle == 1)  then
              TEMP_STEAM_p(k) = temp_steam_pulse(i)
              TEMP_GAS_p(k) = temp_gas_pulse(i)
              press_inner_p(k) = press_inner_pulse(i)
              press_outer_p(k) = press_outer_pulse(i)
            else if (no_pulse_per_cycle == 2)  then
              ! use the first pulse
              TEMP_STEAM_p(k) = temp_steam_pulse(1)
              TEMP_GAS_p(k) = temp_gas_pulse(1)
              press_inner_p(k) = press_inner_pulse(1)
              press_outer_p(k) = press_outer_pulse(1)
            endif

            Temp_ave = Temp_ave + 0.5 * time_event_pulse * &
                       & (TEMP_STEAM_p(k-2) + TEMP_STEAM_p(k-1))

            event_local_point = 0  ! event complete

          endif

          ! id for end of ramp-up and beginning of a new pulse segment
          id_temp_op(k) = 2

          k = k + 1
          Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_pulse(1)

          if (no_pulse_per_cycle == 1)  then
            TEMP_STEAM_p(k) = temp_steam_pulse(i)
            TEMP_GAS_p(k) = temp_gas_pulse(i)
            press_inner_p(k) = press_inner_pulse(i)
            press_outer_p(k) = press_outer_pulse(i)
          else if (no_pulse_per_cycle == 2)  then
            ! end of first pulse
            TEMP_STEAM_p(k) = temp_steam_pulse(1)
            TEMP_GAS_p(k) = temp_gas_pulse(1)
            press_inner_p(k) = press_inner_pulse(1)
            press_outer_p(k) = press_outer_pulse(1)
 
          endif

          ! id for end of pulse and beginning of a new ramp-down segment
          id_temp_op(k) = 3

          ! check for event
          if (k > 4*no_pulse_per_cycle .and. i < no_total_pulse)  then

            if (Time_Boiler_p(k-4*no_pulse_per_cycle) < Time_event_now .and. &
              & Time_Boiler_p(k) >= Time_event_now)  then

              event_local_point = 1
              no_event = no_event + 1
              event_start(no_event) = k
              write(out_lun, *) 'event should occur ', no_event, &
                 & Time_Boiler_p(k-4*no_pulse_per_cycle), &
                 & Time_event_now, Time_Boiler_p(k)

              ! next time event
              Time_event_now = Time_event_now + event_time_interval


            endif

          endif

          if (no_pulse_per_cycle == 1)  then
            Temp_ave = Temp_ave + time_boiler_pulse(1) * temp_steam_pulse(i)
          else if (no_pulse_per_cycle == 2)  then
            ! the first pulse
            Temp_ave = Temp_ave + time_boiler_pulse(1) * temp_steam_pulse(1)
          endif

          k = k + 1

          if (i == no_total_pulse)  then

            if (no_pulse_per_cycle == 1)  then
              ! preparing for the last idle period
              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              TEMP_STEAM_p(k) = temp_steam_last_idle
              TEMP_GAS_p(k) = temp_gas_last_idle
              total_time = total_time + time_idle_pulse + &
                       &  time_boiler_pulse(1) + time_pulse_idle

              press_inner_p(k) = press_inner_last_idle
              press_outer_p(k) = press_outer_last_idle

            else if (no_pulse_per_cycle == 2)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              ! pulse 1 ends; but the cycle is not over yet
              TEMP_STEAM_p(k) = temp_steam_idle(1) 
              TEMP_GAS_p(k) = temp_gas_idle(1) 
              total_time = total_time + time_idle_pulse + &
                    &  time_boiler_pulse(1) + time_pulse_idle + &
                    &  time_boiler_idle(1)
              press_inner_p(k) = press_inner_idle(1)
              press_outer_p(k) = press_outer_idle(1)

              Temp_ave = Temp_ave + time_boiler_idle(1) * temp_steam_idle(1)
            endif
            ! no_pulse_per_cycle = 2

          else 

            if (event_local_point == 0)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              ! irrespective of no_pulse_per_cycle; pulse 1 ends 
              TEMP_STEAM_p(k) = temp_steam_idle(1) 
              TEMP_GAS_p(k) = temp_gas_idle(1) 
              total_time = total_time + time_idle_pulse + &
                    &  time_boiler_pulse(1) + time_pulse_idle + &
                    &  time_boiler_idle(1)
              press_inner_p(k) = press_inner_idle(1)
              press_outer_p(k) = press_outer_idle(1)

              Temp_ave = Temp_ave + time_boiler_idle(1) * temp_steam_idle(1)

            else if (event_local_point == 1)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_event
              TEMP_STEAM_p(k) = temp_steam_event 
              TEMP_GAS_p(k) = temp_gas_event 
              total_time = total_time + time_event_pulse + &
                &  time_boiler_pulse(1) + time_pulse_event + time_boiler_event
              Temp_ave = Temp_ave + time_boiler_event * temp_steam_event

              press_inner_p(k) = press_inner_event
              press_outer_p(k) = press_outer_event

              ! increment local event point
              event_local_point = event_local_point + 1

            endif

          endif 

          ! check if printout for ramp-down segment
          F2L_LOOP1: do j = 1, no_full2low_output

            if (no_pulse_full2low_output(j) == i)  then

              no_point_full2low(j) = k-1  ! start of ramp down
              time_full2low(j, 1) = Time_Boiler_p(k-1) 
              time_full2low(j, mramp) = Time_Boiler_p(k) 
              temp_gas_full2low(j, 1) = TEMP_GAS_p(k-1)
              temp_gas_full2low(j, mramp) = TEMP_GAS_p(k)
              temp_steam_full2low(j, 1) = TEMP_STEAM_p(k-1)
              temp_steam_full2low(j, mramp) = TEMP_STEAM_p(k)
              press_gas_full2low(j, 1) = press_outer_p(k-1)
              press_gas_full2low(j, mramp) = press_outer_p(k)
              press_steam_full2low(j, 1) = press_inner_p(k-1)
              press_steam_full2low(j, mramp) = press_inner_p(k)

              write(out_lun, *) 'no_f2l_cycle ', &
               & (no_point_full2low(j1), j1 = 1, j)  ! no_full2low_output)
              write(out_lun, *) 'no_f2l_time ', &
               & (time_full2low(j1, 1), j1 = 1, j)

              exit F2L_LOOP1

            endif

          end do F2L_LOOP1

          ! check if printout for the EVENT ramp-down segment
          do j = no_full2low_output + 1, &
               & no_full2low_output + no_f2l_event_out

            if (no_event_full2low_output(j - no_full2low_output) == &
              & no_event .and. event_local_point > 0)  then

              no_point_full2low(j) = k-1  ! start of ramp down
              time_full2low(j, 1) = Time_Boiler_p(k-1) 
              time_full2low(j, mramp) = Time_Boiler_p(k) 
              temp_gas_full2low(j, 1) = TEMP_GAS_p(k-1)
              temp_gas_full2low(j, mramp) = TEMP_GAS_p(k)
              temp_steam_full2low(j, 1) = TEMP_STEAM_p(k-1)
              temp_steam_full2low(j, mramp) = TEMP_STEAM_p(k)
              press_gas_full2low(j, 1) = press_outer_p(k-1)
              press_gas_full2low(j, mramp) = press_outer_p(k)
              press_steam_full2low(j, 1) = press_inner_p(k-1)
              press_steam_full2low(j, mramp) = press_inner_p(k)

              write(out_lun, *) 'no_f2l_event ', &
               & (no_point_full2low(j1), j1 = &
               & 1, no_full2low_output + no_f2l_event_out)
               ! & no_full2low_output+1, no_full2low_output + no_f2l_event_out)

            endif

          end do

          ! id for end of ramp-down and beginning of a new idle segment
          id_temp_op(k) = 4

          if (event_local_point == 0)  then
            Temp_ave = Temp_ave + 0.5 * time_pulse_idle * &
                       & (TEMP_STEAM_p(k-1) + TEMP_STEAM_p(k))
          else if (event_local_point == 2)  then
            Temp_ave = Temp_ave + 0.5 * time_pulse_event * &
                       & (TEMP_STEAM_p(k-1) + TEMP_STEAM_p(k))
          endif

          SECOND_PULSE: if (no_pulse_per_cycle == 2)  then

            ! do not consider any events on the second pulse of the cycle
            k = k + 1
         
            ! the first idle starts now (between pulses)
            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_idle(1)
            TEMP_STEAM_p(k) = temp_steam_idle(1)
            TEMP_GAS_p(k) = temp_gas_idle(1)
            press_inner_p(k) = press_inner_idle(1)
            press_outer_p(k) = press_outer_idle(1)

            ! id for end of idle and beginning of a new cycle ramp up
            id_temp_op(k) = 5

            k = k + 1

            ! n_pulse = 1 + (k-3)/4, this should be just the pulse number; i

            ! second pulse
            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_idle_pulse
            TEMP_STEAM_p(k) = temp_steam_pulse(2)
            TEMP_GAS_p(k) = temp_gas_pulse(2)
            press_inner_p(k) = press_inner_pulse(2)
            press_outer_p(k) = press_outer_pulse(2)

            Temp_ave = Temp_ave + 0.5 * time_idle_pulse * &
                       & (TEMP_STEAM_p(k-2) + TEMP_STEAM_p(k-1))

            ! id for end of ramp-up and beginning of a new pulse segment
            id_temp_op(k) = 6

            k = k + 1
            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_pulse(2)

            ! end of second pulse
            TEMP_STEAM_p(k) = temp_steam_pulse(2)
            TEMP_GAS_p(k) = temp_gas_pulse(2)
            press_inner_p(k) = press_inner_pulse(2)
            press_outer_p(k) = press_outer_pulse(2)
 
            ! id for end of pulse and beginning of a new ramp-down segment
            id_temp_op(k) = 7

            ! second pulse
            Temp_ave = Temp_ave + time_boiler_pulse(2) * temp_steam_pulse(2)

            k = k + 1

            if (i == no_total_pulse)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              TEMP_STEAM_p(k) = temp_steam_last_idle
              TEMP_GAS_p(k) = temp_gas_last_idle
              total_time = total_time + time_idle_pulse + &
                       &  time_boiler_pulse(2) + time_pulse_idle

              press_inner_p(k) = press_inner_last_idle
              press_outer_p(k) = press_outer_last_idle

            else 

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              ! pulse 2 ends 
              TEMP_STEAM_p(k) = temp_steam_idle(2) 
              TEMP_GAS_p(k) = temp_gas_idle(2) 
              total_time = total_time + time_idle_pulse + &
                    &  time_boiler_pulse(2) + time_pulse_idle + &
                    &  time_boiler_idle(2)
              press_inner_p(k) = press_inner_idle(2)
              press_outer_p(k) = press_outer_idle(2)

              Temp_ave = Temp_ave + time_boiler_idle(2) * temp_steam_idle(2)

            endif 

            ! id for end of ramp-down of 2 pulse and beginning of a new idle segment
            id_temp_op(k) = 8

            Temp_ave = Temp_ave + 0.5 * time_pulse_idle * &
                       & (TEMP_STEAM_p(k-1) + TEMP_STEAM_p(k))

          endif SECOND_PULSE

        end do

        k = k + 1
        Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_last_idle
        TEMP_STEAM_p(k) = temp_steam_last_idle
        TEMP_GAS_p(k) = temp_gas_last_idle

        press_inner_p(k) = press_inner_last_idle
        press_outer_p(k) = press_outer_last_idle

        Temp_ave = Temp_ave / total_time
        Temp_ave_local =(temp_steam_pulse(1) * time_boiler_pulse(1) + &
                     & temp_steam_idle(1) * time_boiler_idle(1)) / &
                     & (time_boiler_pulse(1) + time_boiler_idle(1))
        if (no_pulse_per_cycle == 2)  then
          Temp_ave_local = Temp_ave_local + (temp_steam_pulse(2) * &
                     & time_boiler_pulse(2) + &
                     & temp_steam_idle(2) * time_boiler_idle(2)) / &
                     & (time_boiler_pulse(2) + time_boiler_idle(2))
        endif

        write(out_lun, *) 'total_time, Temp_ave, &
             & Temp_ave/temp_steam_pulse(1), &
             & Temp_ave_local/temp_steam_pulse(1)'
        write(out_lun, 207) total_time, Temp_ave, &
             & Temp_ave / temp_steam_pulse(1), &
             & Temp_ave_local / temp_steam_pulse(1)
 207    format('Temp_ave ', 10(1pe13.6, 1x))
 
        if (no_event > 0)  then

          write(out_lun, *) 'Events no= ', no_event
          write(out_lun, 210) (Time_Boiler_p(event_start(i)), i = 1, no_event)
 210      format('Start event times ', 12(1pe13.6, 1x))

          no_full2low_output = no_full2low_output + no_f2l_event_out
          write(out_lun, *) ' total no_full2low_output (including events)= ', &
            & no_full2low_output

          write(out_lun, *) 'no_f2l (including event ', &
               & (no_point_full2low(j1), j1 = 1, no_full2low_output)

        else
          write(out_lun, *) 'NO events'
        endif
        
      ! endif

    RETURN

  END SUBROUTINE BOILER_SCHEDULE_TP_CYCLES

  SUBROUTINE BOILER_SCHEDULE_TP_SLIDE(if_cycle)
  
    use parameter_module,   only: mrow_ggt, pi, zero, mramp, mramp_out
    use input_utilities_module,   only: SEARCH_NML
    use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
    use boiler_data_module, only: tube_thickness, tube_outer_radius, &
          & tube_curvature, &
          & tube_inner_radius, htc_tube_inner, htc_tube_outer, &
          & Time_boiler, Temp_steam, Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, & 
          & no_total_pulse, time_boiler_idle, time_boiler_pulse, &
          & time_pulse_idle, temp_steam_pulse, temp_steam_idle, &
          & temp_steam_first_idle, time_boiler_first_idle, &
          & temp_steam_last_idle, time_boiler_last_idle, time_idle_pulse
    use boiler_data_module, only: temp_gas_pulse, temp_gas_idle, &
          & temp_gas_first_idle, &
          & temp_gas_last_idle, temp_gas_slope_p, temp_gas_p, &
          & gas_flow_rate, oxide_location, no_pulse_per_cycle, &
          & press_inner_p, press_outer_p, press_inner_pulse, &
          & press_outer_pulse, press_inner_first_idle, &
          & press_outer_first_idle, &
          & press_inner_idle, press_outer_idle, &
          & press_inner_last_idle, press_outer_last_idle, &
          & press_inner_slope_p, press_outer_slope_p
    use boiler_data_module, only:    event_time_interval, time_boiler_event, &
          & temp_steam_event, temp_gas_event, press_inner_event, &
          & press_outer_event, time_pulse_event, time_event_pulse
    use boiler_data_module, only:   time_pulse_idle_load, &
      & load_fraction_pi_time, time_idle_pulse_load, &
      & load_fraction_ip_time, load_fraction_idle, &
      & load_temp_steam, temp_steam_vs_load, press_steam_vs_load, &
      & load_temp_gas, temp_gas_vs_load, &
      & press_gas_vs_load, heat_flux_fraction_vs_load
    use boiler_data_module, only: id_temp_op, id_press_op, id_flow_op
    use solution_data_module, only: no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low, time_full2low, &
      & time_full2low, temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low, &
      & no_f2l_event_out, no_event_point_full2low, &
      & no_event_full2low_output

    ! Argument List
    logical, intent(IN) :: if_cycle

    ! Local Variables
    integer :: ioerror, i, k, j, j1, &
        & no_pct_first_pulse = 2, no_pct_pulse = 3, &
        & no_pct_idle = 3, no_pct_joint = 1, no_pct_last_pulse = 3, Status
    integer :: event_local_point
    integer, save :: no_event
    integer, save, dimension(200) :: event_start

    logical :: no_compr_namelist, compr_namelist

    real    :: Temp_ave, total_time, Temp_ave_local, Time_event_now

    real    :: speed_axial_over_rotational
    
        if (no_pulse_per_cycle == 1)  then
          no_pct_first_pulse = 2
          no_pct_pulse = 3
          no_pct_idle = 3
          no_pct_joint = 1 
          no_pct_last_pulse = 3
        else if (no_pulse_per_cycle == 2)  then 
          no_pct_first_pulse = 2
          no_pct_pulse = 3
          no_pct_idle = 3
          no_pct_joint = 1 
          no_pct_last_pulse = 3
        else 
          write(6, *) 'STOP: no_pulse_per_cycle > 2 or < 1'
          STOP
        endif

        TEMP_BOILER_data_no = no_pulse_per_cycle * &
          & (no_total_pulse * no_pct_pulse + &
          & (no_total_pulse -1) * no_pct_idle - &
          & (no_total_pulse -1) * 2 * no_pct_joint + 1) + &
          & no_pct_first_pulse - no_pct_joint + &
          & no_pct_last_pulse - no_pct_joint -1

        write(out_lun,*) 'allocate boiler cycle ', TEMP_BOILER_data_no
        ALLOCATE(TEMP_STEAM_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(TEMP_GAS_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(Time_Boiler_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(TEMP_STEAM_Slope_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(TEMP_GAS_Slope_p(1:TEMP_BOILER_data_no), STAT = Status)

        ALLOCATE(id_temp_op(1:TEMP_BOILER_data_no), STAT = Status)

        ! pressure cycle
        ALLOCATE(press_inner_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(press_outer_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(press_inner_slope_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(press_outer_slope_p(1:TEMP_BOILER_data_no), STAT = Status)

        no_full2low_output = 0
        
        do i = mramp_out, 1, -1
          if (no_full2low_output == 0 .and. &
            & no_pulse_full2low_output(i) /= 0) &
            & no_full2low_output = i
        end do

        ! get the event output
        no_f2l_event_out = 0
        do i = 5, 1, -1
          if (no_f2l_event_out == 0 .and. &
            & no_event_full2low_output(i) /= 0) &
            & no_f2l_event_out = i
        end do

        write(out_lun, *) ' no_full2low_output = ', no_full2low_output
        write(out_lun, *) ' no_f2l_event_out = ', no_f2l_event_out
        ! do i = no_full2low_output +1, no_full2low_output + no_f2l_event_out

        ! end do

        ! set the values and times based on number of pulses
        Time_Boiler_p(1) = 0.0
        TEMP_STEAM_p(1) = temp_steam_first_idle
        TEMP_GAS_p(1) = temp_gas_first_idle
        press_inner_p(1) = press_inner_first_idle
        press_outer_p(1) = press_outer_first_idle

        k = 1
        Temp_ave = 0.0
        total_time = 0.0

        ! initialize the time for event now
        Time_event_now = event_time_interval
        event_local_point = 0
        no_event = 0

        do i=1, no_total_pulse

          k = k + 1
          if (k == 2)  then
            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_first_idle
            TEMP_STEAM_p(k) = temp_steam_first_idle
            TEMP_GAS_p(k) = temp_gas_first_idle
            press_inner_p(k) = press_inner_first_idle
            press_outer_p(k) = press_outer_first_idle
          else 

            if (event_local_point == 0)  then
              ! no event now
              if (no_pulse_per_cycle == 1)  then
                Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_idle(1)
                TEMP_STEAM_p(k) = temp_steam_idle(1)
                TEMP_GAS_p(k) = temp_gas_idle(1)
                press_inner_p(k) = press_inner_idle(1)
                press_outer_p(k) = press_outer_idle(1)
              else if (no_pulse_per_cycle == 2)  then
                ! the second pulse ends here
                Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_idle(2)
                TEMP_STEAM_p(k) = temp_steam_idle(2)
                TEMP_GAS_p(k) = temp_gas_idle(2)
                press_inner_p(k) = press_inner_idle(2)
                press_outer_p(k) = press_outer_idle(2)
              endif

            else if (event_local_point == 2)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_event
              TEMP_STEAM_p(k) = temp_steam_event
              TEMP_GAS_p(k) = temp_gas_event
              press_inner_p(k) = press_inner_event
              press_outer_p(k) = press_outer_event

              ! increment local event point
              event_local_point = event_local_point + 1

            endif

            ! id for end of idle and beginning of a new cycle ramp up
            id_temp_op(k) = 1

          endif

          k = k + 1

          ! n_pulse = 1 + (k-3)/4, this should be just the pulse number; i

          ! last point in the local event
          if (event_local_point == 0)  then
            ! no event

            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_idle_pulse
            if (no_pulse_per_cycle == 1)  then
              TEMP_STEAM_p(k) = temp_steam_pulse(i)
              TEMP_GAS_p(k) = temp_gas_pulse(i)
              press_inner_p(k) = press_inner_pulse(i)
              press_outer_p(k) = press_outer_pulse(i)
            else if (no_pulse_per_cycle == 2)  then
              ! first pulse
              TEMP_STEAM_p(k) = temp_steam_pulse(1)
              TEMP_GAS_p(k) = temp_gas_pulse(1)
              press_inner_p(k) = press_inner_pulse(1)
              press_outer_p(k) = press_outer_pulse(1)
            endif

            Temp_ave = Temp_ave + 0.5 * time_idle_pulse * &
                       & (TEMP_STEAM_p(k-2) + TEMP_STEAM_p(k-1))

          else if (event_local_point == 3)  then

            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_event_pulse

            if (no_pulse_per_cycle == 1)  then
              TEMP_STEAM_p(k) = temp_steam_pulse(i)
              TEMP_GAS_p(k) = temp_gas_pulse(i)
              press_inner_p(k) = press_inner_pulse(i)
              press_outer_p(k) = press_outer_pulse(i)
            else if (no_pulse_per_cycle == 2)  then
              ! use the first pulse
              TEMP_STEAM_p(k) = temp_steam_pulse(1)
              TEMP_GAS_p(k) = temp_gas_pulse(1)
              press_inner_p(k) = press_inner_pulse(1)
              press_outer_p(k) = press_outer_pulse(1)
            endif

            Temp_ave = Temp_ave + 0.5 * time_event_pulse * &
                       & (TEMP_STEAM_p(k-2) + TEMP_STEAM_p(k-1))

            event_local_point = 0  ! event complete

          endif

          ! id for end of ramp-up and beginning of a new pulse segment
          id_temp_op(k) = 2

          k = k + 1
          Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_pulse(1)

          if (no_pulse_per_cycle == 1)  then
            TEMP_STEAM_p(k) = temp_steam_pulse(i)
            TEMP_GAS_p(k) = temp_gas_pulse(i)
            press_inner_p(k) = press_inner_pulse(i)
            press_outer_p(k) = press_outer_pulse(i)
          else if (no_pulse_per_cycle == 2)  then
            ! end of first pulse
            TEMP_STEAM_p(k) = temp_steam_pulse(1)
            TEMP_GAS_p(k) = temp_gas_pulse(1)
            press_inner_p(k) = press_inner_pulse(1)
            press_outer_p(k) = press_outer_pulse(1)
 
          endif

          ! id for end of pulse and beginning of a new ramp-down segment
          id_temp_op(k) = 3

          ! check for event
          if (k > 4*no_pulse_per_cycle .and. i < no_total_pulse)  then

            if (Time_Boiler_p(k-4*no_pulse_per_cycle) < Time_event_now .and. &
              & Time_Boiler_p(k) >= Time_event_now)  then

              event_local_point = 1
              no_event = no_event + 1
              event_start(no_event) = k
              Time_event_now = Time_event_now + event_time_interval

            endif

          endif

          if (no_pulse_per_cycle == 1)  then
            Temp_ave = Temp_ave + time_boiler_pulse(1) * temp_steam_pulse(i)
          else if (no_pulse_per_cycle == 2)  then
            ! the first pulse
            Temp_ave = Temp_ave + time_boiler_pulse(1) * temp_steam_pulse(1)
          endif

          k = k + 1

          if (i == no_total_pulse)  then

            if (no_pulse_per_cycle == 1)  then
              ! preparing for the last idle period
              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              TEMP_STEAM_p(k) = temp_steam_last_idle
              TEMP_GAS_p(k) = temp_gas_last_idle
              total_time = total_time + time_idle_pulse + &
                       &  time_boiler_pulse(1) + time_pulse_idle

              press_inner_p(k) = press_inner_last_idle
              press_outer_p(k) = press_outer_last_idle

            else if (no_pulse_per_cycle == 2)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              ! pulse 1 ends; but the cycle is not over yet
              TEMP_STEAM_p(k) = temp_steam_idle(1) 
              TEMP_GAS_p(k) = temp_gas_idle(1) 
              total_time = total_time + time_idle_pulse + &
                    &  time_boiler_pulse(1) + time_pulse_idle + &
                    &  time_boiler_idle(1)
              press_inner_p(k) = press_inner_idle(1)
              press_outer_p(k) = press_outer_idle(1)

              Temp_ave = Temp_ave + time_boiler_idle(1) * temp_steam_idle(1)
            endif
            ! no_pulse_per_cycle = 2

          else 

            if (event_local_point == 0)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              ! irrespective of no_pulse_per_cycle; pulse 1 ends 
              TEMP_STEAM_p(k) = temp_steam_idle(1) 
              TEMP_GAS_p(k) = temp_gas_idle(1) 
              total_time = total_time + time_idle_pulse + &
                    &  time_boiler_pulse(1) + time_pulse_idle + &
                    &  time_boiler_idle(1)
              press_inner_p(k) = press_inner_idle(1)
              press_outer_p(k) = press_outer_idle(1)

              Temp_ave = Temp_ave + time_boiler_idle(1) * temp_steam_idle(1)

            else if (event_local_point == 1)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_event
              TEMP_STEAM_p(k) = temp_steam_event 
              TEMP_GAS_p(k) = temp_gas_event 
              total_time = total_time + time_event_pulse + &
                &  time_boiler_pulse(1) + time_pulse_event + time_boiler_event
              Temp_ave = Temp_ave + time_boiler_event * temp_steam_event

              press_inner_p(k) = press_inner_event
              press_outer_p(k) = press_outer_event

              ! increment local event point
              event_local_point = event_local_point + 1

            endif

          endif 

          ! check if printout for ramp-down segment
          do j = 1, no_full2low_output

            if (no_pulse_full2low_output(j) == i)  then

              no_point_full2low(j) = k-1  ! start of ramp down
              time_full2low(j, 1) = Time_Boiler_p(k-1) 
              time_full2low(j, mramp) = Time_Boiler_p(k) 
              temp_gas_full2low(j, 1) = TEMP_GAS_p(k-1)
              temp_gas_full2low(j, mramp) = TEMP_GAS_p(k)
              temp_steam_full2low(j, 1) = TEMP_STEAM_p(k-1)
              temp_steam_full2low(j, mramp) = TEMP_STEAM_p(k)
              press_gas_full2low(j, 1) = press_outer_p(k-1)
              press_gas_full2low(j, mramp) = press_outer_p(k)
              press_steam_full2low(j, 1) = press_inner_p(k-1)
              press_steam_full2low(j, mramp) = press_inner_p(k)

              write(out_lun, *) 'no_f2l ', &
               & (no_point_full2low(j1), j1 = 1, j)  ! no_full2low_output)

            endif

          end do

          ! check if printout for the EVENT ramp-down segment
          do j = no_full2low_output + 1, &
               & no_full2low_output + no_f2l_event_out

            if (no_event_full2low_output(j - no_full2low_output) == &
              & no_event .and. event_local_point > 0)  then

              no_point_full2low(j) = k-1  ! start of ramp down
              time_full2low(j, 1) = Time_Boiler_p(k-1) 
              time_full2low(j, mramp) = Time_Boiler_p(k) 
              temp_gas_full2low(j, 1) = TEMP_GAS_p(k-1)
              temp_gas_full2low(j, mramp) = TEMP_GAS_p(k)
              temp_steam_full2low(j, 1) = TEMP_STEAM_p(k-1)
              temp_steam_full2low(j, mramp) = TEMP_STEAM_p(k)
              press_gas_full2low(j, 1) = press_outer_p(k-1)
              press_gas_full2low(j, mramp) = press_outer_p(k)
              press_steam_full2low(j, 1) = press_inner_p(k-1)
              press_steam_full2low(j, mramp) = press_inner_p(k)

              write(out_lun, *) 'no_f2l_event ', &
               & (no_point_full2low(j1), j1 = &
               & 1, no_full2low_output + no_f2l_event_out)
               ! & no_full2low_output+1, no_full2low_output + no_f2l_event_out)

            endif

          end do

          ! id for end of ramp-down and beginning of a new idle segment
          id_temp_op(k) = 4

          if (event_local_point == 0)  then
            Temp_ave = Temp_ave + 0.5 * time_pulse_idle * &
                       & (TEMP_STEAM_p(k-1) + TEMP_STEAM_p(k))
          else if (event_local_point == 2)  then
            Temp_ave = Temp_ave + 0.5 * time_pulse_event * &
                       & (TEMP_STEAM_p(k-1) + TEMP_STEAM_p(k))
          endif

          SECOND_PULSE: if (no_pulse_per_cycle == 2)  then

            ! do not consider any events on the second pulse of the cycle
            k = k + 1
         
            ! the first idle starts now (between pulses)
            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_idle(1)
            TEMP_STEAM_p(k) = temp_steam_idle(1)
            TEMP_GAS_p(k) = temp_gas_idle(1)
            press_inner_p(k) = press_inner_idle(1)
            press_outer_p(k) = press_outer_idle(1)

            ! id for end of idle and beginning of a new cycle ramp up
            id_temp_op(k) = 5

            k = k + 1

            ! n_pulse = 1 + (k-3)/4, this should be just the pulse number; i

            ! second pulse
            TEMP_STEAM_p(k) = temp_steam_pulse(2)
            TEMP_GAS_p(k) = temp_gas_pulse(2)
            press_inner_p(k) = press_inner_pulse(2)
            press_outer_p(k) = press_outer_pulse(2)

            Temp_ave = Temp_ave + 0.5 * time_idle_pulse * &
                       & (TEMP_STEAM_p(k-2) + TEMP_STEAM_p(k-1))

            ! id for end of ramp-up and beginning of a new pulse segment
            id_temp_op(k) = 6

            k = k + 1
            Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_pulse(2)

            ! end of second pulse
            TEMP_STEAM_p(k) = temp_steam_pulse(2)
            TEMP_GAS_p(k) = temp_gas_pulse(2)
            press_inner_p(k) = press_inner_pulse(2)
            press_outer_p(k) = press_outer_pulse(2)
 
            ! id for end of pulse and beginning of a new ramp-down segment
            id_temp_op(k) = 7

            ! second pulse
            Temp_ave = Temp_ave + time_boiler_pulse(2) * temp_steam_pulse(2)

            k = k + 1

            if (i == no_total_pulse)  then

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              TEMP_STEAM_p(k) = temp_steam_last_idle
              TEMP_GAS_p(k) = temp_gas_last_idle
              total_time = total_time + time_idle_pulse + &
                       &  time_boiler_pulse(2) + time_pulse_idle

              press_inner_p(k) = press_inner_last_idle
              press_outer_p(k) = press_outer_last_idle

            else 

              Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_pulse_idle
              ! pulse 2 ends 
              TEMP_STEAM_p(k) = temp_steam_idle(2) 
              TEMP_GAS_p(k) = temp_gas_idle(2) 
              total_time = total_time + time_idle_pulse + &
                    &  time_boiler_pulse(2) + time_pulse_idle + &
                    &  time_boiler_idle(2)
              press_inner_p(k) = press_inner_idle(2)
              press_outer_p(k) = press_outer_idle(2)

              Temp_ave = Temp_ave + time_boiler_idle(2) * temp_steam_idle(2)

            endif 

            ! id for end of ramp-down of 2 pulse and beginning of a new idle segment
            id_temp_op(k) = 8

            Temp_ave = Temp_ave + 0.5 * time_pulse_idle * &
                       & (TEMP_STEAM_p(k-1) + TEMP_STEAM_p(k))

          endif SECOND_PULSE

        end do

        k = k + 1
        Time_Boiler_p(k) = Time_Boiler_p(k-1) + time_boiler_last_idle
        TEMP_STEAM_p(k) = temp_steam_last_idle
        TEMP_GAS_p(k) = temp_gas_last_idle

        press_inner_p(k) = press_inner_last_idle
        press_outer_p(k) = press_outer_last_idle

        Temp_ave = Temp_ave / total_time
        Temp_ave_local =(temp_steam_pulse(1) * time_boiler_pulse(1) + &
                     & temp_steam_idle(1) * time_boiler_idle(1)) / &
                     & (time_boiler_pulse(1) + time_boiler_idle(1))
        if (no_pulse_per_cycle == 2)  then
          Temp_ave_local = Temp_ave_local + (temp_steam_pulse(2) * &
                     & time_boiler_pulse(2) + &
                     & temp_steam_idle(2) * time_boiler_idle(2)) / &
                     & (time_boiler_pulse(2) + time_boiler_idle(2))
        endif

        write(out_lun, *) 'total_time, Temp_ave, &
             & Temp_ave/temp_steam_pulse(1), &
             & Temp_ave_local/temp_steam_pulse(1)'
        write(out_lun, 207) total_time, Temp_ave, &
             & Temp_ave / temp_steam_pulse(1), &
             & Temp_ave_local / temp_steam_pulse(1)
 207    format('Temp_ave ', 10(1pe13.6, 1x))
 
        if (no_event > 0)  then

          write(out_lun, *) 'Events no= ', no_event
          write(out_lun, 210) (Time_Boiler_p(event_start(i)), i = 1, no_event)
 210      format('Start event times ', 12(1pe13.6, 1x))

          no_full2low_output = no_full2low_output + no_f2l_event_out
          write(out_lun, *) ' total no_full2low_output (including events)= ', &
            & no_full2low_output

          write(out_lun, *) 'no_f2l (including event ', &
               & (no_point_full2low(j1), j1 = 1, no_full2low_output)

        else
          write(out_lun, *) 'NO events'
        endif
        
      ! endif

    RETURN

  END SUBROUTINE BOILER_SCHEDULE_TP_SLIDE

  END MODULE BOILER_OPERATION_MODULE
