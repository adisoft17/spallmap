MODULE OXIDE_TEST_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   just for debugging against exact solutions. 
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: CHECK_EXACT_STRAIN_CYL

 CONTAINS

  SUBROUTINE CHECK_EXACT_STRAIN_CYL(no_total_layers)
    !=======================================================================
    ! Purpose(s):
    !
    !  check the strains at different temperatures and %Fe2O3 
    !  for the tube, i.e., hollow cylinder case
    ! 
    !=======================================================================
    use parameter_module,   only: mstep, moxide, oxide_thick2si, &
          & id_low, id_high, mave, ten_p3, mramp, zero, mramp_out
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & temp_steam_pulse, temp_steam_idle, no_total_pulse, &
          & Temp_gas_p, press_inner_slope_p, press_outer_slope_p, &
          & press_inner_p, press_outer_p
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, e_th, e_th_mat, &
          & id_mat, no_layer, material_name, poisson_ratio
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
           & ntemp_check_strain, &
           & small_oxide_thickness
    use solver_data_module, only  : generalized_plane_strain
    use oxide_data_module,   only: id_fe3o4, no_oxide_ave, &
           & Youngs_modul, thickness_fr, thickness_fr_test
    use boiler_data_module, only: id_temp_op
    use solution_driver_module, only: DRIVER_SOLUTION
    use stress_utility_module, only: HOOP_STRESS
    use solution_data_module, only: Temp, Temp_ave, tau_th_st, &
      & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
      & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
      & strain_th_ave, stress_hoop_th_ave, npr_temp, stress_rad, &
      & strain_th_plate_ave, strain_hoop, stress_hoop, npr_st, &
      & no_pulse_full2low_output, eps0, &
      & no_full2low_output, no_point_full2low, time_full2low, &
      & time_full2low, temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low, &
      & integral_th, rad_int, bc_st, cc_st, strain_hoop_gen_ave
    use solution_data_module, only: no_oxide_layers
    use solver_data_module, only: if_average_oxide
    use oxide_stress_data_module,   only: total_oxide_thickness_ref

    ! Argument List
    integer, intent(IN)  ::  no_total_layers 

    ! Local Variables
    integer  :: i, j, k, step_cycle, nstart, no_out_f2l
    real     :: Temp_High, Temp_Low, Temp_gas_now, Temp_steam_now, &
              & Press_gas_now, Press_steam_now, slope_time, &
              & slope_temp_gas, slope_temp_steam, slope_press_gas, &
              & slope_press_steam, time_now, &
              & press_outer_test, press_inner_test, &
              & Temp_gas_test, Temp_steam_test, &
              & a_check, b_check, poisson_star, eps0_ex, &
              & Youngs_star, eps0_s_factor, time_test, dtime_test
    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'
    real, dimension(moxide) :: e_th_scale, thickness_layer
    real, dimension(moxide) :: dtemp, dstrain_th, &
          & dstrain_hoop, dstress_hoop
    real, dimension(mramp_out)  :: press_full_int
    real, dimension(no_total_layers)     :: bc_ex, cc_ex, J_int, tau_star
    real, dimension(2*no_total_layers)     :: e_hoop, s_hoop, s1_hoop

    if (.not. if_average_oxide)  then
      write(6, *) 'TEST starts for multiple layer case'
    endif

    write(6, *) 'TEST begin'
    write(aux_lun, *) 'test begin'
    write(out_lun, *) 'test begin'

      do i = 1, 1

        oxide_thickness(i) = 200.0 ! microns

        ! uniform high temperature
        Temp_gas_test = 590.0
        Temp_steam_test = 590.0

        ! non-uniform high temperature
        Temp_gas_test = 1800.0
        Temp_steam_test = 560.0

        press_outer_test = 0.0
        press_inner_test = 0.0

        thickness_layer = oxide_thickness(i) * oxide_thick2si

        ! oxide thickness [micron]
        ! without thermal expasion; just given by the kinetics
        total_oxide_thickness_ref = oxide_thickness(i)

        time_test = 0.0
        dtime_test = 8.0  ! for creep

        ! use the test thickness fraction
        thickness_fr = thickness_fr_test

        ! obtain temperature, stress, strain, dimensions -> high temp
        call DRIVER_SOLUTION(no_oxide_layers, TEMP_BOILER_data_no+10, &
              & Temp_gas_test, Temp_steam_test, &
              & press_outer_test, press_inner_test, thickness_layer, &
              & time_test, dtime_test, id_high, if_average_oxide) 

        ! write average temperatures out
  ! mave = 1 - max, 2 - min, 3 - sum average (layer, mave)
  ! strain_th_ave - is the average quantity for the tau 

            ! write the output for the paper at the full load
            ! write temperatures
            write(aux_lun, 19) i, oxide_thickness(i), &
              & Temp(1, 1), Temp(2, npr_temp(2)), Temp(2, 1)
 19     format('test_ull_load_temp ', i4, 1x, 1(1pe13.6, 1x), &
      & 4(1pe10.3, 1x))

            ! write hoop strains and stresses for full load
            write(aux_lun, 20) i, oxide_thickness(i), &
              & strain_hoop_ave(2, 1)*ten_p3, strain_hoop_ave(2, 3)*ten_p3, &
              & strain_hoop_ave(2, 2)*ten_p3, stress_hoop_ave(2, 1), &
              & stress_hoop_ave(2, 3), stress_hoop_ave(2, 2)
 20     format('test_ull_load_hoop ', i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))

        ! uniform low temperature
        Temp_gas_test = 460.0
        Temp_steam_test = 460.0

        ! non-uniform low temperature
        Temp_gas_test = 1600.0
        Temp_steam_test = 400.0
        dtime_test = 8.0  ! for creep

        ! use the test thickness fraction
        thickness_fr = thickness_fr_test

        ! obtain temperature, stress, strain, dimensions -> low temp
        call DRIVER_SOLUTION(no_oxide_layers, TEMP_BOILER_data_no+10, &
              & Temp_gas_test, Temp_steam_test, &
              & press_outer_test, press_inner_test, thickness_layer, &
              & time_test, dtime_test, id_low, if_average_oxide) 

            ! write the output for the paper;
            ! the last call was made for the low load

            ! write temperatures
            write(aux_lun, 21) i, oxide_thickness(i), &
              & Temp(1, 1), Temp(2, npr_temp(2)), Temp(2, 1)
 21     format('test_ow_load_temp ', i4, 1x, 1(1pe13.6, 1x), &
      & 4(1pe10.3, 1x))

            ! write pressures, full load, low load at oxide metal interface
            ! write(aux_lun, 23) i, oxide_thickness(i), &
            !  & -press_full_int(j), -stress_rad(2, 1)
 23     format('test_ull_low_load_press ', i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe10.3, 1x))

            ! write hoop strains and stresses for low load
            write(aux_lun, 22) i, oxide_thickness(i), &
              & strain_hoop_ave(2, 1)*ten_p3, strain_hoop_ave(2, 3)*ten_p3, &
              & strain_hoop_ave(2, 2)*ten_p3, stress_hoop_ave(2, 1), &
              & stress_hoop_ave(2, 3), stress_hoop_ave(2, 2)
 22     format('test_ow_load_hoop ', i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))

   ! strain_hoop_gen_ave  is the effective hoop strain that generates stress
   ! strain_hoop_gen_ave = total_strain - thermal_free_stress_strain
            write(aux_lun, 27) i, oxide_thickness(i), &
       & strain_hoop_gen_ave(2, 1)*ten_p3, strain_hoop_gen_ave(2, 3)*ten_p3, &
       & strain_hoop_gen_ave(2, 2)*ten_p3, stress_hoop_ave(2, 1), &
              & stress_hoop_ave(2, 3), stress_hoop_ave(2, 2)
 27     format('test_ow_load_hoop_effective_strain_generating_stress ', &
                 & i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))

            write(aux_lun, 26) i, oxide_thickness(i), &
              & tau_th_st(1, 1)*ten_p3, &
              & strain_hoop(1, 1)*ten_p3, &  ! metal outer surface 
              & tau_th_st(1, npr_st(1))*ten_p3, &
              & strain_hoop(1, npr_st(1))*ten_p3, &  ! metal-oxide 
              & tau_th_st(2, 1)*ten_p3, &
              & strain_hoop(2, 1)*ten_p3             ! oxide-metal
 26     format('test_ow_load_th_hoop_interface ', i4, 1x, 1(1pe13.6, 1x), &
                 & 6(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))

            write(aux_lun, 25) i, oxide_thickness(i), &
              & strain_rad_ave(2, 1)*ten_p3, strain_rad_ave(2, 3)*ten_p3, &
              & strain_rad_ave(2, 2)*ten_p3, stress_rad_ave(2, 1), &
              & stress_rad_ave(2, 3), stress_rad_ave(2, 2)
 25     format('test_ow_load_rad ', i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))

            write(aux_lun, 24) i, oxide_thickness(i), &
              & strain_th_plate_ave(2, 1)*ten_p3, &
              & strain_th_plate_ave(2, 2)*ten_p3, &
              & strain_th_plate_ave(2, 3)*ten_p3,  &
              & stress_hoop_th_ave(2, 1), stress_hoop_th_ave(2, 3), &
              & stress_hoop_th_ave(2, 2)
 24     format('test_ow_load_hoop_plate ', i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))

        ! exact hoop strain solution

        do j = 1, no_total_layers
          J_int(j) = integral_th(j, 1) * (1.0 + poisson_ratio(1))
          tau_star(j) = tau_th_st(j, j) * (1.0 + poisson_ratio(1))
        end do

        poisson_star = poisson_ratio(1) / (1.0 - poisson_ratio(1))
        Youngs_star = Youngs_modul(i, 1) / (1.0 - poisson_ratio(1)**2)

        a_check = 1.0 / (1.0 - poisson_star)
        b_check = 1.0 / (1.0 + poisson_star)
        eps0_s_factor = poisson_ratio(1) * a_check

        ! bc_ex(1) = J_int(2) / &
        !         & (a_check*(rad_int(1)**2 - rad_int(3)**2) - J_int(1))
        ! bc_ex(2) = bc_ex(1)

        if (no_total_layers == 2)  then

          if (generalized_plane_strain)  then

            eps0_ex = 2.0 * (J_int(1)+ J_int(2)) / &
             & ((1.0 + poisson_ratio(1)) * (rad_int(1)**2 - rad_int(3)**2))
            bc_ex(1) = -eps0_ex * poisson_ratio(1) + &
                 & (J_int(1)+ J_int(2)) / &
                 & (a_check * (rad_int(1)**2 - rad_int(3)**2))
            bc_ex(2) = bc_ex(1)

            cc_ex(2) = (bc_ex(2) + eps0_ex * poisson_ratio(1)) * &
                   & a_check * rad_int(3)**2 / b_check
            cc_ex(1) = (bc_ex(1) + eps0_ex * poisson_ratio(1)) * &
                   & a_check * rad_int(1)**2 / b_check - &
                   & J_int(1) / b_check

          else if (.not. generalized_plane_strain)  then

            eps0_ex = 0.0
            bc_ex(1) = (J_int(1)+ J_int(2)) / &
                 & (a_check*(rad_int(1)**2 - rad_int(3)**2))
            bc_ex(2) = bc_ex(1)

            cc_ex(2) = bc_ex(2) * a_check * rad_int(3)**2 / b_check
            cc_ex(1) = bc_ex(1) * a_check * rad_int(1)**2 / b_check - &
                 & J_int(1) / b_check

          endif

        else if (no_total_layers == 3)  then

          if (generalized_plane_strain)  then

            eps0_ex = 2.0 * (J_int(1)+ J_int(2)+ J_int(3)) / &
             & ((1.0 + poisson_ratio(1)) * (rad_int(1)**2 - rad_int(4)**2))
            bc_ex(1) = -eps0_ex * poisson_ratio(1) + &
                 & (J_int(1)+ J_int(2)+ J_int(3)) / &
                 & (a_check * (rad_int(1)**2 - rad_int(4)**2))
            bc_ex(2) = bc_ex(1)
            bc_ex(3) = bc_ex(1)

            cc_ex(3) = (bc_ex(3) + eps0_ex * poisson_ratio(1)) * &
                   & a_check * rad_int(4)**2 / b_check
            cc_ex(2) = cc_ex(3) + J_int(3) / b_check
            cc_ex(1) = cc_ex(2) + J_int(2) / b_check

          else if (.not. generalized_plane_strain)  then

            eps0_ex = 0.0
            bc_ex(1) = (J_int(1)+ J_int(2)+J_int(3)) / &
                 & (a_check*(rad_int(1)**2 - rad_int(4)**2))
            bc_ex(2) = bc_ex(1)
            bc_ex(3) = bc_ex(1)

            cc_ex(3) = bc_ex(3) * a_check * rad_int(4)**2 / b_check
            cc_ex(2) = cc_ex(3) + J_int(3) / b_check
            cc_ex(1) = cc_ex(2) + J_int(2) / b_check

          endif

        else if (no_total_layers > 3)  then

          write(aux_lun, *) 'END test; no test for stress; too many layers'
          return

        endif

        write(aux_lun, *) 'test_ow_load_eps0 = ', eps0, eps0_ex
        write(aux_lun, *) 'test q ', Youngs_star, eps0_s_factor
        write(aux_lun, *) 'test q1 ', bc_ex(2) * a_check, &
                  & -tau_star(2), b_check * cc_ex(2) / rad_int(3)**2, &
                  & eps0_s_factor * eps0
        ! stop
 
        write(aux_lun, *) 'test J1_check ', (J_int(k), k = 1, no_total_layers)
        write(aux_lun, *) 'test J2_check ', (-J_int(k) / &
          & (130.0 * 0.5 * (1.0 + poisson_ratio(1)) * &
          & (rad_int(k)**2 - rad_int(k+1)**2)), k = 1, no_total_layers)

        if (ABS(bc_ex(1)) > 0.0)  then
          write(aux_lun, *) 'test bc_check ', bc_ex(1), bc_st(1), &
            & 100.0 * abs((bc_st(1) - bc_ex(1))/bc_ex(1))
        endif

        do k = 1, no_total_layers
          if (ABS(cc_ex(1)) > 0.0)  then 
            write(aux_lun, *) 'test cc1_check ', k, cc_ex(k), cc_st(k), &
              & 100.0 * abs((cc_st(k) - cc_ex(k))/cc_ex(k))
          endif
        end do

        do k = 1, no_total_layers + 1

          if (k >= 2)  then

            ! exact solution for hoop strain at oxide-oxide interface
            ! e_hoop(1) = bc_ex(1) * ( 1.0 + a_check / b_check)
            e_hoop(k-1) = bc_ex(k-1) + cc_ex(k-1) / rad_int(k)**2
            s_hoop(k-1) = Youngs_star * &
                  & (bc_ex(k-1) * a_check - tau_star(k-1) + &
                  & b_check * cc_ex(k-1) / rad_int(k)**2 + &
                  & eps0_s_factor * eps0)
            s1_hoop(k-1) = HOOP_STRESS(Youngs_star, eps0 - eps0, &
               & rad_int(k)**2, bc_ex(k-1) * a_check, &
               & b_check * cc_ex(k-1), eps0_s_factor * eps0, &
               & tau_star(k-1))
            if (ABS(s1_hoop(k-1)) > 0.0)  then
              write(aux_lun, *) 's_hoop_test ', k-1, s_hoop(k-1), &
               & s1_hoop(k-1), &
               & 100.0 * abs((1.0-s_hoop(k-1)/s1_hoop(k-1)))
            endif
            if (ABS(e_hoop(k-1)) > 0.0)  then
              write(aux_lun, *) 'e1_test ', e_hoop(k-1), &
              & strain_hoop(k-1, npr_st(k-1)), &
              & 100.0 * abs((1.0-strain_hoop(k-1, npr_st(k-1))/e_hoop(k-1))) 
            endif
          endif

          if (k <= no_total_layers)  then

            ! exact solution for hoop strain at metal-oxide interface
            if (ABS(b_check) > 0.0)  then
              e_hoop(k) = bc_ex(k) + cc_ex(k) / rad_int(k)**2 + &
                   & J_int(k) / (b_check * rad_int(k)**2)
              s_hoop(k) = Youngs_star * &
                  & (bc_ex(k) * a_check - tau_star(k) + &
                  & b_check * cc_ex(k) / rad_int(k)**2 + &
                  & J_int(k) / rad_int(k)**2 + &
                  & eps0_s_factor * eps0)
              s1_hoop(k) = HOOP_STRESS(Youngs_star, J_int(k), &
               & rad_int(k)**2, bc_ex(k) * a_check, &
               & b_check * cc_ex(k), eps0_s_factor * eps0, &
               & tau_star(k))
              if (ABS(s1_hoop(k)) > 0.0)  then
                write(aux_lun, *) 's_hoop_test ', k, s_hoop(k), &
                 & s1_hoop(k), &
                 & 100.0 * abs((1.0-s_hoop(k)/s1_hoop(k)))
              endif
              if (ABS(e_hoop(k)) > 0.0)  then
                write(aux_lun, *) 'e1_test ', e_hoop(k), strain_hoop(k, 1), &
                 & 100.0 * abs((1.0-strain_hoop(k, 1)/e_hoop(k))) 
              endif
            endif

          endif

        end do

   end do

   write(6, *) 'TEST END'
    write(aux_lun, *) 'test end'
    write(out_lun, *) 'test end'
   ! stop

   return

  END SUBROUTINE CHECK_EXACT_STRAIN_CYL

END MODULE OXIDE_TEST_MODULE
