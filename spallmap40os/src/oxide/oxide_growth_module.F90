MODULE OXIDE_GROWTH_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   subroutines for oxide growth. 
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: OXIDE_GROWTH, OXIDE_GROWTH_LAYER

 CONTAINS

  SUBROUTINE OXIDE_GROWTH()

    !=======================================================================
    ! Purpose(s):
    !
   ! Compute the oxide thickness based on kinetics and temperature schedule 
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, oxide_thick2si, id_high, id_low
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, Temp_gas_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, &
          & oxide_thickness, id_temp_op, thickness_fr_var
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, &
          & if_temp_steam_or_oxide, thickness_fr, no_oxide_ave, &
          & growth_ox_temperature
    use solver_data_module, only: no_thick_interval, if_average_oxide
    use thermal_cycle_module, only: TEMP_PROFILE
    use mesh_module, only: GET_RADIUS, GET_MESH
    use solution_data_module, only: Temp, npr_temp, no_oxide_layers

    implicit none

    ! Local Variables
    integer :: i, j, k, Status, n_it_max, n1, n2, n21, ik, it_dummy = 0
    real    :: kp, delta_t, kp_log, temp2, temp1, temp1w, temp2w, &
             & one_over_temp, temp_ave, time2, error_temp, &
             & error_temp_max, dtemp
    real    :: time_dummy = 0.0
    real, dimension(:), pointer      :: dz2, Temp_dummy
    real, dimension(moxide) :: thickness_layer
    logical  :: one_check

    ! units are in microns
    ALLOCATE(oxide_thickness(1:TEMP_BOILER_data_no + 5), STAT = Status)

    ! for metal only do not obtain the oxide thickness
    if (no_oxide_layers == 1)  RETURN

    ALLOCATE(dz2(1:TEMP_BOILER_data_no), STAT = Status)
    ALLOCATE(Temp_dummy(1:TEMP_BOILER_data_no), STAT = Status)

    error_temp_max = 0.1

    n_it_max = 20
    do i = 2, no_oxide_layers
      if (TRIM(growth_ox_temperature(i)) == 'steam')  then
        ! oxide grows at steam temperature
        n_it_max = 1
      else
        n_it_max = 20
      endif
    end do

    Temp_dummy = TEMP_STEAM_p

    ! initiate the error
    error_temp = -10.0
    n1 = TEMP_BOILER_data_no / 7
    n2 = TEMP_BOILER_data_no / 2
    n21 = n2 - n1 + 1

    TEMP_LOOP:  do ik = 1, n_it_max

    dz2 = 0.0

    do i= 2, TEMP_BOILER_data_no ! no of data intervals

      ! get the times and temperatures
      delta_t = (Time_Boiler_p(i) - Time_Boiler_p(i-1)) / no_thick_interval
      temp1 = TEMP_dummy(i-1)

      ! initialize
      dz2(i) = dz2(i-1) 
      time2 = Time_Boiler_p(i-1)

      SMALL_LOOP: do k = 1, no_thick_interval

        ! assume the same slope for the steam as for the Temperature 
        ! of the oxide surface

        temp2 = temp1 + delta_t * TEMP_STEAM_Slope_p(i-1)
        time2 = time2 + delta_t
        ! get the oxidation rate constant for this time interval
        ! linear interpolation based on 1/T 

        ! oxide growth based on T_dummy
        call KP_OXIDE(temp1, temp2, kp)

        dz2(i) = dz2(i) + 2.0 * delta_t * kp

        if (i == 2)  then
          write(out_lun, 31) i, time2, kp, 2.0*sqrt(kp), dz2(i),sqrt(dz2(i)) 
 31  format('dz_check ', i4, 1x, 10(1pe13.6, 1x))
        endif

        temp1 = temp2

      end do SMALL_LOOP

      ! oxide_thickness(i) = oxide_thickness(i-1) + sqrt(dz2)

    end do

    oxide_thickness = sqrt(dz2)

    ! there should be no need for this iteration since for steam temp only one iter will be peformed
    ! if (if_temp_steam_or_oxide) exit TEMP_LOOP
    if (n_it_max == 1)  exit TEMP_LOOP

    ! use dz2 as a dummy variable
    dz2 = -10.0

    do j = n1, n2

      if (id_temp_op(j) == 3)  then

        thickness_layer = thickness_fr(1:moxide) * &
                   & oxide_thickness(j) * oxide_thick2si

        call GET_RADIUS(no_oxide_layers, thickness_layer, &
          & time_dummy, if_average_oxide, it_dummy)

        call GET_MESH(no_oxide_layers, thickness_layer, if_average_oxide)

        ! temp profile for two layers (metal + one oxide layer)
        call TEMP_PROFILE(no_oxide_layers, Temp_gas_p(j), Temp_steam_p(j), &
           & thickness_layer, j, Time_Boiler_p(j), &
           & Time_Boiler_p(j) - Time_Boiler_p(j), id_high, if_average_oxide)  ! high temperature
        ! two layers; asumming growth at the oxide steam interface
        dz2(j) = abs(Temp_dummy(j) - Temp(no_oxide_layers, npr_temp(no_oxide_layers)))

      else if (id_temp_op(j) == 4)  then

         call GET_RADIUS(no_oxide_layers, thickness_layer, time_dummy, &
            & if_average_oxide, it_dummy)

         call GET_MESH(no_oxide_layers, thickness_layer, if_average_oxide)

         call TEMP_PROFILE(no_oxide_layers, Temp_gas_p(j), Temp_steam_p(j), &
           & thickness_layer, j, Time_Boiler_p(j), &
           & Time_Boiler_p(j) - Time_Boiler_p(j), id_low, if_average_oxide)  ! low temperature
        ! two layers; asumming growth at the oxide steam interface
        dz2(j) = abs(Temp_dummy(j) - Temp(no_oxide_layers, npr_temp(no_oxide_layers)))

      endif

    end do

    error_temp = MAXVAL(dz2, dz2 > 0.0)

        write(tty_lun, 42) ik, error_temp / error_temp_max
 42     format(' eps loop', i4, 1x, 10(1pe13.6, 1x))

    if (error_temp < error_temp_max)  then 

      exit TEMP_LOOP

    else

      Temp_dummy(1:16) = TEMP_STEAM_p(1:16)

      do j = 17, TEMP_BOILER_data_no

        thickness_layer = thickness_fr(1:moxide) * &
                   & oxide_thickness(j) * oxide_thick2si

        call GET_RADIUS(no_oxide_layers, thickness_layer, time_dummy, &
             & if_average_oxide, it_dummy)

        call GET_MESH(no_oxide_layers, thickness_layer, if_average_oxide)

        ! temp profile for two layers (metal + one oxide layer)
        call TEMP_PROFILE(no_oxide_layers, Temp_gas_p(j), Temp_steam_p(j), &
           & thickness_layer, j, Time_Boiler_p(j), &
           & Time_Boiler_p(j) - Time_Boiler_p(j), 1, if_average_oxide)  ! do not assign high or low temperature
        Temp_dummy(j) = Temp(no_oxide_layers, npr_temp(no_oxide_layers))

      end do
      
    endif

    end do TEMP_LOOP

    DEALLOCATE(Temp_dummy)
    DEALLOCATE(dz2)

    do i = 2, TEMP_BOILER_data_no
       write(out_lun,209) i, id_temp_op(i), Time_Boiler_p(i), TEMP_STEAM_p(i), &
       & oxide_thickness(i), kp, kp_log, oxide_thickness(i)**2/(2.0*Time_Boiler_p(i))
 208    format('oxide_val ', 10(1pe13.6, 1x))
 209    format('oxide_val ', i4, 1x, i2, 1x, 10(1pe13.6, 1x))
      
    end do

    return

  END SUBROUTINE OXIDE_GROWTH

  SUBROUTINE OXIDE_GROWTH_LAYER()

    !=======================================================================
    ! Purpose(s):
    !
   ! Compute the oxide thickness based on kinetics and temperature schedule 
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, oxide_thick2si, id_high, id_low
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, Temp_gas_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, &
          & oxide_thickness, id_temp_op, thickness_fr_var
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, &
          & if_temp_steam_or_oxide, thickness_fr, no_oxide_ave, &
          & growth_ox_temperature, thickness_growth_fr, first_oxide_layer
    use solver_data_module, only: no_thick_interval, if_average_oxide
    use thermal_cycle_module, only: TEMP_PROFILE
    use mesh_module, only: GET_RADIUS, GET_MESH
    use solution_data_module, only: Temp, npr_temp, no_oxide_layers, &
          & oxide_thickness_hf_out, time_hf_out, &
          & id_full_load_hf, n_hf_out
    use solution_data_module, only: no_full2low_output, no_point_full2low

    implicit none

    ! Local Variables
    integer :: i, j, k, Status, n_it_max, n1, n2, n21, ik, ilay, &
             & ihf, it_dummy = 1
    real    :: delta_t, temp1w, temp2w, error_temp_max, &
             & one_over_temp, temp_ave, time2
    real    :: time_dummy = 0.0
    real, dimension(2:no_oxide_layers)    :: kp, kp_log, temp2, temp1,  &
             & temp_slope, error_temp, dtemp, growth_temp_this_layer
    real, dimension(:, :), pointer      :: dz2, Temp_growth_now, &
             & Temp_growth_update
    real, dimension(moxide) :: thickness_layer
    logical  :: one_check

    ! units are in microns
    ALLOCATE(oxide_thickness(1:TEMP_BOILER_data_no + 5), STAT = Status)

    ! for metal only do not obtain the oxide thickness
    if (no_oxide_layers == 1)  then

      oxide_thickness = 0.0
      ALLOCATE(thickness_fr_var(1:TEMP_BOILER_data_no, 2:2), &
             & STAT = Status)
      thickness_fr_var(:, 2) = 1.0

      RETURN
    endif

    if (SUM(thickness_growth_fr) < 1.0e-6)  then

      write(out_lun, *) 'ERROR: old run needs thickness_growth_fr', &
         'default growth_ox_temperature = oxide_surface'
      STOP

    endif

    ALLOCATE(thickness_fr_var(1:TEMP_BOILER_data_no, 2:no_oxide_layers), &
             & STAT = Status)
    ALLOCATE(dz2(1:TEMP_BOILER_data_no, 2:no_oxide_layers), STAT = Status)
    ALLOCATE(Temp_growth_now(1:TEMP_BOILER_data_no, 2:no_oxide_layers), &
             & STAT = Status)

    error_temp_max = 0.1
    thickness_fr(1:first_oxide_layer) = 0.0

    n_it_max = 20
    do i = first_oxide_layer, no_oxide_layers

      ! initialize the growth temperature per each layer
      Temp_growth_now(1:TEMP_BOILER_data_no, i) = &
        & TEMP_STEAM_p(1:TEMP_BOILER_data_no)

      if (TRIM(growth_ox_temperature(i)) == 'steam')  then
        ! oxide grows at steam temperature
        n_it_max = 1
      else
        n_it_max = 20
      endif

    end do

    if (n_it_max > 1) then
      ALLOCATE(Temp_growth_update(1:TEMP_BOILER_data_no, 2:no_oxide_layers), &
             & STAT = Status)
      Temp_growth_update = Temp_growth_now
    endif

    ! initiate the error
    error_temp = -10.0

    ! this loop is needed since the growth temperature is known apriori
    LARGE_LOOP:  do ik = 1, n_it_max 

      dz2 = 0.0

      ! step 1 - get oxide thickness at all times
      THICKNESS_LOOP: do i= 2, TEMP_BOILER_data_no ! no of data intervals

        ! get the time at the start of the i-th interval and its increment
        delta_t = (Time_Boiler_p(i) - Time_Boiler_p(i-1)) / &
                & no_thick_interval
        time2 = Time_Boiler_p(i-1)

        do ilay = first_oxide_layer, no_oxide_layers 

          ! temperature at interval start
          temp1(ilay) = Temp_growth_now(i-1, ilay)

          temp_slope(ilay) = (Temp_growth_now(i, ilay) - &
              & Temp_growth_now(i-1, ilay)) / &
              & (Time_Boiler_p(i) - Time_Boiler_p(i-1))

          ! initialize to previous
          dz2(i, ilay) = dz2(i-1, ilay) 

        end do

        ! loop over the time increments of this interval
        TIME_INCREMENT: do k = 1, no_thick_interval

          ! get the time
          time2 = time2 + delta_t

          ! assume the same slope for the steam as for the Temperature 
          ! of the oxide surface

          LAYER_LOOP1: do ilay = first_oxide_layer, no_oxide_layers 

            ! get the steam temperature now
            ! temp2(ilay) = temp1(ilay) + delta_t * TEMP_STEAM_Slope_p(i-1) not ok
            temp2(ilay) = temp1(ilay) + delta_t * temp_slope(ilay)

            ! get the oxidation rate constant for this time interval
            ! linear interpolation based on 1/T 

            ! oxide growth based on T_dummy
            call KP_OXIDE(temp1(ilay), temp2(ilay), kp(ilay))

            dz2(i, ilay) = dz2(i, ilay) + 2.0 * delta_t * kp(ilay) * &
                  & thickness_growth_fr(ilay)**2

            temp1(ilay) = temp2(ilay)

          end do LAYER_LOOP1

          ! get the total oxide thickness of all the layers
          oxide_thickness(i) = SUM(SQRT(dz2(i, &
             & first_oxide_layer:no_oxide_layers)))

          if (i == 2)  then
            write(out_lun, 31) i, time2, kp(first_oxide_layer), &
               & 2.0*sqrt(kp(first_oxide_layer)), &
               & dz2(i, first_oxide_layer:no_oxide_layers), &
               & sqrt(dz2(i, first_oxide_layer:no_oxide_layers)), &
               & oxide_thickness(i) 
 31    format('dz_check ', i4, 1x, 30(1pe13.6, 1x))
          endif

          if (oxide_thickness(i) > 1.0e+6 .or. &
             & ANY(ABS(Temp_growth_now(i, :)) > 2.0e+3))  then
            write(tty_lun, *) 'ERROR: oxide_thickness update error; iter= ', ik, i
            write(tty_lun,52) ik, i, oxide_thickness(i), (Temp_growth_now(i, j), &
           & Temp_growth_update(i, j),  j = first_oxide_layer, no_oxide_layers)
            stop
          endif

        end do TIME_INCREMENT

      end do THICKNESS_LOOP

      ! store the oxide thickness in the fraction place for the time being
      ! thickness_fr_var = sqrt(dz2)

      ! store the thickness fraction of each layer
      do i = first_oxide_layer, no_oxide_layers

        ! since dz2(1) = 0.0 and oxide_thickness(1)
        thickness_fr_var(1, i) = SQRT(dz2(2, i)) / &
             & oxide_thickness(2)

        thickness_fr_var(2:TEMP_BOILER_data_no, i) = &
             & SQRT(dz2(2:TEMP_BOILER_data_no, i)) / &
             & oxide_thickness(2:TEMP_BOILER_data_no)
      end do

      ! check for convergence in all the layers
      ! use dz2 as a dummy variable
      ! oxide is too small initially, thus use steam temperature
      dz2(1:4, :) = 0.0
      do i = first_oxide_layer, no_oxide_layers
        Temp_growth_update(1:4, i) = Temp_growth_now(1:4, i)
      end do

      TEMP_GROWTH_LOOP: do j = 5, TEMP_BOILER_data_no

        ! update the thickness_fr value 
        thickness_fr(first_oxide_layer:no_oxide_layers) = &
            & thickness_fr_var(j, first_oxide_layer:no_oxide_layers) 

        ! current oxide thickness
        thickness_layer = thickness_fr(1:moxide) * &
                 & oxide_thickness(j) * oxide_thick2si

        call GET_RADIUS(no_oxide_layers, thickness_layer, &
                   & time_dummy, if_average_oxide, it_dummy)

        call GET_MESH(no_oxide_layers, thickness_layer, if_average_oxide)

        ! temp profile for two layers (metal + one oxide layer)
        call TEMP_PROFILE(no_oxide_layers, Temp_gas_p(j), &
          & Temp_steam_p(j), thickness_layer, j, Time_Boiler_p(j), &
          & Time_Boiler_p(j) - Time_Boiler_p(j), 1, if_average_oxide)  

        ! get the growth temperature
        call GET_GROWTH_TEMP(no_oxide_layers, Temp_steam_p(j), &
           & growth_temp_this_layer)

        ! dz2 stores the maximum temperature deviation
        dz2(j, first_oxide_layer:no_oxide_layers) = &
           & ABS(Temp_growth_now(j, first_oxide_layer:no_oxide_layers) - &
           & growth_temp_this_layer(first_oxide_layer:no_oxide_layers))

        Temp_growth_update(j, first_oxide_layer:no_oxide_layers) = &
           & growth_temp_this_layer(first_oxide_layer:no_oxide_layers)

        write(out_lun,52) ik, j, (thickness_layer(i), Temp_growth_now(j, i), &
           & Temp_growth_update(j, i),  i = first_oxide_layer, no_oxide_layers)
 52     format(' temp_conv1 ', i2, 1x, i4, 10(1pe13.6, 1x))

        if (ANY(thickness_layer > 1.0) .or. &
           & ANY(ABS(Temp_growth_now(j, :)) > 2.0e+3) .or. &
           & ANY(ABS(Temp_growth_update(j, :)) > 2.0e+3))  then
          write(tty_lun, *) 'ERROR: Temp_growth update error; iter= ', ik, j
          write(tty_lun,52) ik, j, (thickness_layer(i), Temp_growth_now(j, i), &
           & Temp_growth_update(j, i),  i = first_oxide_layer, no_oxide_layers)
          stop
        endif

      end do TEMP_GROWTH_LOOP
      
      error_temp(first_oxide_layer:no_oxide_layers) = &
        & MAXVAL(dz2(:, first_oxide_layer:no_oxide_layers), &
        & dz2(:, first_oxide_layer:no_oxide_layers) > 0.0)

      write(tty_lun, 42) ik, &
        & error_temp(first_oxide_layer:no_oxide_layers) / error_temp_max
 42     format('eps_loop ', i4, 1x, 10(1pe13.6, 1x))

      if (MAXVAL(error_temp) < error_temp_max)  then 

        exit LARGE_LOOP

      else

        Temp_growth_now = Temp_growth_update
        ! Temp_growth_now = 0.5 * Temp_growth_update + &
        !   & (1.0 - 0.5) * Temp_growth_now

      endif

    end do LARGE_LOOP

    DEALLOCATE(Temp_growth_now)
    DEALLOCATE(dz2)
    if (n_it_max > 1) then
      DEALLOCATE(Temp_growth_update)
    endif

    do i = 2, TEMP_BOILER_data_no
       write(out_lun,209) i, id_temp_op(i), Time_Boiler_p(i),  &
       & TEMP_STEAM_p(i), oxide_thickness(i), &
       & oxide_thickness(i)**2/(2.0*Time_Boiler_p(i)), &
       & thickness_fr_var(i, first_oxide_layer:no_oxide_layers)
 208    format('oxide_val ', 15(1pe13.6, 1x))
 209    format('oxide_val ', i4, 1x, i2, 1x, 15(1pe13.6, 1x))
      
    end do

       ! get the heat flux output data
  ! store the heat flux information at different oxide thicknesses
  ! by default choose 10 oxide thickness distributed uniformly
  n_hf_out = 10
  do ihf = 1, n_hf_out
    ! current oxide thickness considered at this output
    oxide_thickness_hf_out(ihf) = oxide_thickness(TEMP_BOILER_data_no-5) * &
         & ihf / n_hf_out

    write(out_lun, *) 'start check_hf_out ', ihf, &
          & oxide_thickness_hf_out(ihf)

    do i = 3, TEMP_BOILER_data_no

      if (oxide_thickness(i) <= oxide_thickness_hf_out(ihf) .and. &
        & oxide_thickness(i+4) > oxide_thickness_hf_out(ihf) .and. &
        & id_temp_op(i) == 3)  then

        time_hf_out(ihf) = Time_Boiler_p(i)
        id_full_load_hf(ihf) = i

        do j = 1, no_full2low_output

          if (i == no_point_full2low(j))  then
            ! outage; need to move this further or earlier
           
            time_hf_out(ihf) = Time_Boiler_p(i+4)
            id_full_load_hf(ihf) = i+4
            oxide_thickness_hf_out(ihf) = oxide_thickness(i+4)

          endif

       end do
          
        write(out_lun, 210) ihf, i, &
          & oxide_thickness_hf_out(ihf), time_hf_out(ihf)
 210    format('check_hf_out ', i2, 1x, i6, 1x, 2(1pe13.6, 1x))
      endif

    end do
  end do

    return

  END SUBROUTINE OXIDE_GROWTH_LAYER

  SUBROUTINE GET_GROWTH_TEMP(N, steam_temperature, growth_temp_this_layer)

  ! obtain the growth temperature per each layer 
  use solution_data_module, only: Temp, npr_temp
  use oxide_data_module,    only: growth_ox_temperature, first_oxide_layer

  ! Argument List
  integer,               intent(IN)  :: N
  real,                  intent(IN)  :: steam_temperature
  real, dimension(2:N),  intent(OUT) :: growth_temp_this_layer

  real, dimension(2:N) :: dummy

  integer :: j, i

    do j = first_oxide_layer, N

      SELECT CASE (TRIM(growth_ox_temperature(j)))

        CASE ('inner_surface')

          growth_temp_this_layer(j) = Temp(j, 1)

        CASE ('outer_surface')

          growth_temp_this_layer(j) = Temp(j, npr_temp(j))

        CASE ('average_layer')

          growth_temp_this_layer(j) = SUM(Temp(j, 1:npr_temp(j))) / &
            & npr_temp(j)

        CASE ('average_all_layers')

          do i = first_oxide_layer, N
            dummy(i) = SUM(Temp(i, 1:npr_temp(i))) / npr_temp(i)
          end do
          growth_temp_this_layer(j) = SUM(dummy(first_oxide_layer:N)) / &
             & (N-first_oxide_layer + 1)

        CASE ('steam')

          growth_temp_this_layer(j) = steam_temperature ! TEMP_STEAM_p(cycle_no)

        CASE ('oxide_surface')

          growth_temp_this_layer(j) = Temp(N, npr_temp(N))

        CASE DEFAULT

          write(6, *) 'ERROR: growth_ox_temperature not ok'
          STOP

      END SELECT

    end do

    RETURN

  END SUBROUTINE GET_GROWTH_TEMP

  SUBROUTINE KP_OXIDE(temp1, temp2, KP)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the oxide thickness based on kinetics and temperature schedule 
    ! 
    !=======================================================================
    use parameter_module,   only: mfuel
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, log_const_energy_oxide, &
          & activ_energy_oxide_r, if_function_or_table_oxide
    use solver_data_module, only: no_thick_interval

    ! Arguments
    real, intent(OUT)    :: kp
    real, intent(INOUT)     :: temp1, temp2

    ! Local Variables
    integer :: i, j, k, Status
    real    :: delta_t, kp_log, temp2w, temp1w, &
             & one_over_temp, temp_ave, time2
    real, dimension(:), pointer      :: dz2
    logical   :: one_check

        one_check = .false.

        if (temp2 < temp1)  then
          temp1w = temp2
          temp2w = temp1
        else if (temp1 < temp2)  then
          temp1w = temp1
          temp2w = temp2
        else 
          one_check = .true.  ! if temperatures are identical
          temp1w = temp1
          temp2w = temp1w         
        endif

        LOOP1: do j = 1, nrate - 1

          CHECK1: if (one_check)  then

            ! constant temperature
            if (temp1w >= oxide_rate_temp(j) .and. &
              & temp1w <= oxide_rate_temp(j+1))  then
              temp_ave = temp1w
              one_over_temp = 1.0 / (temp_kelvin_unit + temp_ave)
              kp_log = (oxide_rate_log(j) * &
                     & (oxide_rate_1temp(j+1) - one_over_temp) + &
                     & oxide_rate_log(j+1) * &
                     & (one_over_temp - oxide_rate_1temp(j))) / &
                     & (oxide_rate_1temp(j+1) - oxide_rate_1temp(j))
      !    write(out_lun, 11) i, k, j, oxide_rate_temp(j), temp1w, temp2w, &
      !   & oxide_rate_temp(j+1), Time_Boiler_p(i), temp_ave, kp_log
      !   11 format('check1 ', 3(i4, 1x), 10(1pe13.6, 1x))

              exit LOOP1

            endif

          else

            ! temp variation
            TEMP1_CHECK: if (temp1w >= oxide_rate_temp(j) .and. &
              & temp1w < oxide_rate_temp(j+1))  then

            if (temp2w >= oxide_rate_temp(j) .and. &
              & temp2w <= oxide_rate_temp(j+1))  then

              temp_ave = 0.5 * (temp2 + temp1)
              one_over_temp = 1.0 / (temp_kelvin_unit + temp_ave)
              kp_log = (oxide_rate_log(j) * &
                     & (oxide_rate_1temp(j+1) - one_over_temp) + &
                     & oxide_rate_log(j+1) * &
                     & (one_over_temp - oxide_rate_1temp(j))) / &
                     & (oxide_rate_1temp(j+1) - oxide_rate_1temp(j))

     !     write(out_lun, 12) i, k, j, oxide_rate_temp(j), temp1w, temp2w, &
     !    & oxide_rate_temp(j+1), Time_Boiler_p(i), temp_ave, kp_log
     !    12 format('check2 ', 3(i4, 1x), 10(1pe13.6, 1x))

              exit LOOP1

            else if (temp2w > oxide_rate_temp(j+1))  then

              if (j == nrate -1)  then

                write(out_lun, *) 'ERROR temp too high'
                stop

              else

                ! the local temperature points are in successive intervals for the
                ! oxide growth data
                temp_ave = 0.5 * (temp1w + oxide_rate_temp(j+1))
              one_over_temp = 1.0 / (temp_kelvin_unit + temp_ave)

             kp_log = (oxide_rate_log(j) * &
                    & (oxide_rate_1temp(j+1) - one_over_temp) + &
                    & oxide_rate_log(j+1) * &
                    & (one_over_temp - oxide_rate_1temp(j))) / &
                    & (oxide_rate_1temp(j+1) - oxide_rate_1temp(j))

                temp_ave = 0.5 * (temp2w + oxide_rate_temp(j+1))
              one_over_temp = 1.0 / (temp_kelvin_unit + temp_ave)

             kp_log = kp_log + &
                 & (oxide_rate_log(j+1) * &
                 & (oxide_rate_1temp(j+2) - one_over_temp) + &
                 & oxide_rate_log(j+2) * &
                 & (one_over_temp - oxide_rate_1temp(j+1))) / &
                 & (oxide_rate_1temp(j+2) - oxide_rate_1temp(j+1))

             kp_log = 0.5 * kp_log
      !    write(out_lun, 13) i, k, j, oxide_rate_temp(j), temp1w, temp2w, &
      !   & oxide_rate_temp(j+1), Time_Boiler_p(i), temp_ave, kp_log
      !   13 format('check3 ', 3(i4, 1x), 10(1pe13.6, 1x))

              exit LOOP1

              endif

            endif

            exit LOOP1

          endif TEMP1_CHECK
        
          endif CHECK1

        end do LOOP1

        kp = EXP(kp_log)

    return

  END SUBROUTINE KP_OXIDE

END MODULE OXIDE_GROWTH_MODULE
