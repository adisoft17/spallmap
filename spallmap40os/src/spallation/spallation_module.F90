MODULE SPALLATION_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for spallation. 
    !
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================


  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: FRACTURE_TENSION, CRITICAL_ETA_OMEGA, &
       & CHECK_MAP_EPRI, CHECK_STRAIN

 CONTAINS

  SUBROUTINE CHECK_MAP_EPRI(average_or_layer_sum)
    !=======================================================================
    ! Purpose(s):
    !
    !  check the EPRI map at different temperatures and %Fe2O3 
    ! 
    !=======================================================================
    use parameter_module,   only: mstep, mdamage, moxide
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & temp_steam_pulse, temp_steam_idle, no_total_pulse
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, e_th, e_th_mat, &
          & id_mat, no_layer, material_name, thickness_fr, ndamage, &
          & no_oxide_ave, first_oxide_layer
    use solver_data_module, only: no_thick_spall_map, map_level, &
           & interval_type_thick_spall_map, thickness_spallation_map
    use solution_data_module, only: no_oxide_layers

    logical, intent(IN) :: average_or_layer_sum

    ! Local Variables
    integer  :: i, j, k, id_fe2o3, id_fe3o4
    real     :: Temp_High, Temp_Low
    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'
    character(LEN = 80), dimension(mdamage)  :: damage_type
    real, dimension(moxide, mdamage, mstep)  :: eta, omega
    real, dimension(mdamage, mstep)  :: eta_ave, omega_ave

    if (no_thick_spall_map == 0) return

    damage_type(1) = 'through_scale_cracking'
    damage_type(2) = 'deflection_delamination'
    damage_type(3) = 'intefacial_crack_growth'
    damage_type(4) = 'buckling'
    damage_type(5) = 'crack_deflection'
    damage_type(6) = 'spalling'

    Temp_Low = 20.0

    ! initialize type of damage

    if (average_or_layer_sum)  then

      ! determine the critical eta_omega for the entire scale
      ! using the average properties for a scale with different layers
      ! compute it based on average property for the entire scale
      do j = 1, no_thick_spall_map

        i = no_oxide_ave
        ! scale has the entire oxide thickness
        ! oxide_thickness(id_mat(i)) = thickness_spallation_map(j) 

        write(aux_lun, *) 'cs deb ', i, id_mat(i), j, &
        & thickness_spallation_map(j), oxide_thickness(id_mat(i))

        do k = 1, ndamage

          call CRITICAL_ETA_OMEGA(id_mat(i), damage_type(k), &
             &  eta_ave(k, j), omega_ave(k, j), &
             &  thickness_spallation_map(j), Temp_Low)

        end do
    
        ! write(out_lun, *) k, ' ',damage_type(k)

       ! could not print the damage out
       ! write(out_lun, 21) thickness_spallation_map(j), &
       ! & damage_type(k), omega_ave(k, j), eta_ave(k, j)
 ! 21 format('damage1 ', 1pe13.6, 1x, A23, 1x, 5(1pe13.6, 1x))

       if (map_level == 1)  then

         if (j == 1)  then
           write(out_lun, 34) ((k, i, i = 1, 1), k=1, ndamage)
 34 format('damage1_ave t_ox ',  50('eta_d', i1, '_L', i1,' '))
         endif

         ! use this one for FIRST generation map
         write(out_lun, 21) thickness_spallation_map(j), &
          &  (1.0e+3 * eta_ave(k, j), k=1, ndamage)
 21 format('damage1_ave ', 1pe13.6, 1x, 20(1pe13.6, 1x))

       else if (map_level == 2)  then

         ! use this one for SECOND generation map
         write(out_lun, 22) thickness_spallation_map(j), &
          &  (omega_ave(k, j), eta_ave(k, j), k=1, ndamage)
 22 format('damage2_ave ', 1pe13.6, 1x, 20(1pe13.6, 1x))

       endif

       ! write(tty_lun, *) 'SP MAP AVE; DONE thickness ', j, &
       !    & thickness_spallation_map(j)

      end do


    else if (.not. average_or_layer_sum)  then

      do j = 1, no_thick_spall_map

        ! for the entire oxide layer
        i = no_oxide_ave
        oxide_thickness(i) = thickness_spallation_map(j) 

        do k = 1, ndamage

          call CRITICAL_ETA_OMEGA(id_mat(i), damage_type(k), &
             &  eta(id_mat(i), k, j), omega(id_mat(i), k, j), &
             &  oxide_thickness(id_mat(i)), Temp_Low)

        end do

        do i = first_oxide_layer, no_layer

          oxide_thickness(id_mat(i)) = thickness_fr(id_mat(i)) * &
           & thickness_spallation_map(j) 

          do k = 1, ndamage

            call CRITICAL_ETA_OMEGA(id_mat(i), damage_type(k), &
             &  eta(id_mat(i), k, j), omega(id_mat(i), k, j), &
             &  oxide_thickness(id_mat(i)), Temp_Low)

          end do

        end do

        if (map_level == 1)  then

         if (j == 1)  then

           if (no_layer == first_oxide_layer)  then  ! one oxide layer

             write(out_lun, 33) ((k, i, i = 2, no_layer), k=1, ndamage)
 33 format('damage1_lay t_ox_total ',  6('eta_d', i1, '_L', i1,' '))

             write(out_lun, 35) (i, (k, i, k=1, ndamage), i = 2, no_layer)
 35 format('damage1_lay t_ox_total ', 't_ox',i1,' ', &
            & 6('eta_d', i1, '_L', i1,' '), 't_ox',i1,' ', &
            & 6('eta_d', i1, '_L', i1,' '))

           else if (no_layer == 3)  then  ! two oxide layers

             write(out_lun, 36) no_oxide_ave, &
                   & (k, no_oxide_ave, k=1, ndamage), &
               & (i, (k, i, k=1, ndamage), i = 2, no_layer)
 36 format('damage1_lay t_ox_total ', 't_ox',i1,' ', &
            & 6('eta_d', i1, '_L', i1,' '), 't_ox',i1,' ', &
            & 6('eta_d', i1, '_L', i1,' '), 't_ox',i1,' ', &
            & 6('eta_d', i1, '_L', i1,' '))

           else if (no_layer == 4)  then  ! three oxide layers

             write(out_lun, 36) (i, (k, i, k=1, ndamage), i = 2, no_layer)

           else

             write(6, *) 'WARNING: more than three oxide layers for spallation'
             ! stop

           endif

         endif

         ! no_oxide_layers == no_layer 
         ! write(6, *) 'check thickness ', no_oxide_layers, no_layer, &
         !  & thickness_spallation_map(j), &
         !  & (id_mat(i), oxide_thickness(id_mat(i)), i= 2, no_layer)
   
         ! use this one for second generation map
         write(out_lun, 31) thickness_spallation_map(j), &
          &  oxide_thickness(id_mat(no_oxide_ave)), &
          &  (1.0e+3 * eta(id_mat(no_oxide_ave), k, j), k=1, ndamage), &
          &  (oxide_thickness(id_mat(i)), (1.0e+3 * eta(id_mat(i), k, j), &
          &  k=1, ndamage), i = first_oxide_layer, no_layer)
 31 format('damage1_lay ', 1pe13.6, 1x, 20(1pe13.6, 1x))

       else

         ! use this one for second generation map
         write(out_lun, 32) thickness_spallation_map(j), &
          &  (omega_ave(k, j), eta_ave(k, j), k=1, ndamage)
 32 format('damage2_lay ', 1pe13.6, 1x, 20(1pe13.6, 1x))

       endif

        ! write(tty_lun, *) 'SP MAP LAY; DONE thickness ', j, &
        ! thickness_spallation_map(j)
    
      end do

    endif

      ! write(out_lun, 11) Temp_High, &
      !  & - e_ave(1:nFe2O3_pct_check_strain) * 1.0e+3, energy(1)+energy(2)
 11   format('epri1_check1 ', 10(1pe13.6, 1x)) 

   return

  END SUBROUTINE CHECK_MAP_EPRI

  SUBROUTINE CHECK_STRAIN()
    !=======================================================================
    ! Purpose(s):
    !
    !  check the strains at different temperatures and %Fe2O3 
    ! 
    !=======================================================================
    use parameter_module,   only: mstep
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & temp_steam_pulse, temp_steam_idle, no_total_pulse
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, e_th, e_th_mat, &
          & id_mat, no_layer, material_name
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
           & Fe2O3_pct_check_strain, ntemp_check_strain, nFe2O3_pct_check_strain
    use oxide_utility_module, only: COOLING_STRAIN_PLANE

    ! Local Variables
    integer  :: i, j, k, id_fe2o3, id_fe3o4
    real     :: Temp_High, Temp_Low
    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'

    if (ntemp_check_strain == 0 .or. nFe2O3_pct_check_strain == 0) return
    Temp_Low = 20.0

    ! identify Fe2O3
    id_fe2o3 = 0
    do k = 1, no_layer
      if (TRIM(material_name(id_mat(k))) == 'Fe2O3')  id_fe2o3 = k
    end do

    ! identify Fe3o4
    id_fe3o4 = 0
    do k = 1, no_layer
      if (TRIM(material_name(id_mat(k))) == 'Fe3O4')  id_fe3o4 = k
    end do

    write(tty_lun, *) 'Fe2O3',  id_fe2o3, 'Fe3O4',  id_fe3o4

    do j = 1, ntemp_check_strain

      Temp_High = temp_check_strain(j)

      ! get the cooling strain for the layers; substrate is layer 1
    
      do i = 1, no_layer
        e_th_mat(i)= COOLING_STRAIN_PLANE(id_mat(i), Temp_High, Temp_Low)
      end do 

      ! include the substrate effect
      e_th = e_th_mat - e_th_mat(1)

      ! careful when fe2O3 is not defined
      if (id_fe2o3 == 0 .or. id_fe3o4 == 0)  then
      
        if (id_fe2o3 == 0)  write(out_lun, *) 'fe2o3 not defined'
        if (id_fe3o4 == 0)  write(out_lun, *) 'fe3o4 not defined'

        return

      end if 

      ! store results for the two layers
      e_check(1) = e_th(id_fe2o3)
      e_check(2) = e_th(id_fe3o4)

      do k = 1, nFe2O3_pct_check_strain
    
        ! fr1 stores he fraction of Fe2O3
        fr(1)= Fe2O3_pct_check_strain(k) / 100.0
        fr(2) = 1.0 - fr(1)
        ! average the strain in the layers
        ! this one neglects the actual stacking of the layers
        e_ave(k) = e_check(1) * fr(1) + e_check(2) * fr(2) 

        energy(1) = STORED_ENERGY(id_fe2o3, fr(1), -e_check(1), mode_stress)
        energy(2) = STORED_ENERGY(id_fe3o4, fr(2), -e_check(2), mode_stress)
  
        energy_ave(k) = energy(1) + energy(2) 

      end do

      write(out_lun, 11) Temp_High, &
       & -e_ave(1:nFe2O3_pct_check_strain) * 1.0e+3, energy(1)+energy(2)
 11   format('echeck1 ', 10(1pe13.6, 1x)) 

   end do

   return

  END SUBROUTINE CHECK_STRAIN

  REAL FUNCTION FRACTURE_TENSION(ilayer, Temp_High, Temp_Low)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the cooling strain based on temperature drop and thermal
    !  expansion for planar biaxial stress distribution
    ! 
    !=======================================================================
    use parameter_module,   only: moxide
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, th_exp_temp, &
          & th_exp_coeff, nth_exp
    use solver_data_module, only: no_thick_interval, no_temp_strain

    ! Argument List
    integer, intent(IN)  :: ilayer
    real,    intent(IN)  :: Temp_High, Temp_Low

    ! Local Variables

    real   :: Temp, delta_temp
    real, dimension(moxide)  :: var, var_ave
    integer :: k, i

    FRACTURE_TENSION = 0.0

    return

  END FUNCTION FRACTURE_TENSION

  REAL FUNCTION STRESS_STRAIN(ilayer, mode_stress, Temperature, strain)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the cooling strain based on temperature drop and thermal
    !  expansion for planar biaxial stress distribution
    ! 
    !=======================================================================
    use parameter_module,   only: moxide
    use output_module,      only: tty_lun, out_lun, aux_lun
    use oxide_data_module,  only: Youngs_modul, Youngs_temp, nyoungs, &
          & poisson_ratio
    use solver_data_module, only: no_thick_interval, no_temp_strain
    use property_module,    only: LINEAR_PROPERTY_SIMPLE

    ! Argument List
    integer, intent(IN)  :: ilayer
    real,    intent(IN)  :: Temperature, strain
    character(LEN = 80),    intent(IN)  :: mode_stress

    ! Local Variables

    real   :: Young_mod
    real, dimension(moxide)  :: var, var_ave
    integer :: k, i

    ! get the Youngs modulus
    call LINEAR_PROPERTY_SIMPLE (Temperature, nyoungs(ilayer) - 1, 0, moxide, &
        & youngs_modul(ilayer, :), youngs_temp(ilayer, :), Young_mod)

    if (TRIM(mode_stress) == 'biaxial_plane')  then

      STRESS_STRAIN = strain * Young_mod / (1.0 - Poisson_ratio(ilayer))

    else if (TRIM(mode_stress) == 'biaxial_axisym') then

      write(out_lun, *) 'ERROR: biaxial_axisym not available'
      stop
 
    else

      write(out_lun, *) 'ERROR: stress_mode not available'
      stop

    endif

   return

  END FUNCTION STRESS_STRAIN

  REAL FUNCTION CRITICAL_STRAIN(ilayer, mode_failure, Temperature)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the critical strain for each failure mode
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, pi
    use output_module,      only: tty_lun, out_lun, aux_lun
    use oxide_data_module,  only: Youngs_modul, Youngs_temp, nyoungs, &
          & poisson_ratio, fracture_toughness, flaw_defect_factor, &
          & flaw_oxide_size_factor, oxide_thickness, surf_fracture_energy
    use solver_data_module, only: no_thick_interval, no_temp_strain
    use property_module,    only: LINEAR_PROPERTY_SIMPLE

    ! Argument List
    integer, intent(IN)  :: ilayer
    real,    intent(IN)  :: Temperature
    character(LEN = 80),    intent(IN)  :: mode_failure

    ! Local Variables

    real   :: Young_mod
    real, dimension(moxide)  :: var, var_ave
    integer :: k, i

    ! get the Youngs modulus
    call LINEAR_PROPERTY_SIMPLE (Temperature, nyoungs(ilayer) - 1, 0, moxide, &
        & youngs_modul(ilayer, :), youngs_temp(ilayer, :), Young_mod)

    if (TRIM(mode_failure) == 'tensile_toughness')  then

      ! rel 11, Nagl end Evans (1993)
      ! this one did not work too well reported by Armitt et al (1978)
      CRITICAL_STRAIN = fracture_toughness(ilayer) / (Young_mod * &
          & flaw_defect_factor * &
          & sqrt(pi * flaw_oxide_size_factor * oxide_thickness(ilayer)))

    else if (TRIM(mode_failure) == 'tensile_surf_energy')  then

      ! rel. 12, Nagl end Evans (1993)
      ! proposed by Armitt et al (1978)
      CRITICAL_STRAIN = sqrt(2.0 * surf_fracture_energy(ilayer) / &
          & (Young_mod * &
          & pi * flaw_oxide_size_factor * oxide_thickness(ilayer))) / &
          & flaw_defect_factor

    else if (TRIM(mode_failure) == 'compress_wedge') then

      write(out_lun, *) 'ERROR: biaxial_axisym not available'
      stop
 
    else if (TRIM(mode_failure) == 'compress_buckle') then

      write(out_lun, *) 'ERROR: biaxial_axisym not available'
      stop
 
    else

      write(out_lun, *) 'ERROR: stress_mode not available'
      stop

    endif

   return

  END FUNCTION CRITICAL_STRAIN

  REAL FUNCTION STORED_ENERGY(ilayer, thickness, strain, mode_stress)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the cooling strain based on temperature drop and thermal
    !  expansion for planar biaxial stress distribution
    ! 
    !=======================================================================
    use parameter_module,   only: moxide
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, th_exp_temp, &
          & th_exp_coeff, nth_exp
    use solver_data_module, only: no_thick_interval, no_temp_strain

    ! Argument List
    integer, intent(IN)  :: ilayer
    real,    intent(IN)  :: thickness, strain
    character(LEN = 80),    intent(IN)  :: mode_stress

    ! Local Variables

    real   :: stress
    real, dimension(moxide)  :: var, var_ave
    integer :: k, i
   
    ! temperature to evaluate the youngs modulus;
    ! use the room temperature for the time being
    stress = STRESS_STRAIN(ilayer, mode_stress, 20.0, strain)
    
    ! this is only in one direction
    ! STORED_ENERGY = 0.5 * thickness * stress * ABS(strain)

    ! for biaxial, Nagl end Evans (1993) ref 10; Steiner and Konys (2006)
    ! there is no factor 0.5 
    STORED_ENERGY = thickness * stress * ABS(strain)

    return

  END FUNCTION STORED_ENERGY

  SUBROUTINE CRITICAL_ETA_OMEGA1(ilayer, mode_failure, &
       & critical_strain, omega, Thickness, Temperature)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the critical strain for each failure mode
    !  Level II of EPRI diagram using:
    !  eta - material properties and geometry effects
    !  omega - damage kinetics parameter dependent on the operation history
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, pi
    use output_module,      only: tty_lun, out_lun, aux_lun
    use oxide_data_module,  only: Youngs_modul, Youngs_temp, nyoungs, &
          & poisson_ratio, fracture_toughness, flaw_defect_factor, &
          & flaw_oxide_size_factor, oxide_thickness, surf_fracture_energy, &
          & amplitude_roughness, &
          & wavelength_roughness, delamination_size_existent, &
          & separated_area, porosity_fraction, &
          & number_cracks_volumetric, delamination_size_factor
    use solver_data_module, only: no_thick_interval, no_temp_strain
    use property_module,    only: LINEAR_PROPERTY_SIMPLE
    ! to deal with damage map at shut down
    use oxide_data_module,  only: const_damage_map

    ! Argument List
    integer, intent(IN)  :: ilayer
    real,    intent(IN)  :: Temperature, Thickness
    character(LEN = 80),    intent(IN)  :: mode_failure
    real,    intent(OUT) :: CRITICAL_STRAIN, omega

    ! Local Variables

    real   :: Young_mod, eta, delamination_size_now
    real, dimension(moxide)  :: var, var_ave
    integer :: k, i, ierr

    ! get the Youngs modulus
    call LINEAR_PROPERTY_SIMPLE (Temperature, nyoungs(ilayer) - 1, 0, moxide, &
        & youngs_modul(ilayer, :), youngs_temp(ilayer, :), Young_mod)

      write(aux_lun, *) 'cs ', ilayer, fracture_toughness(ilayer), Young_mod, &
          & flaw_defect_factor, &
          & pi, flaw_oxide_size_factor, Thickness

    if (delamination_size_factor > 0.0)  then
      delamination_size_now = Thickness * delamination_size_factor
    else
      delamination_size_now = delamination_size_existent
    endif

    ierr = 0

    if (TRIM(mode_failure) == 'through_scale_cracking')  then

      ! rel 11, Nagl end Evans (1993)
      ! this one did not work too well reported by Armitt et al (1978)
  
      CRITICAL_STRAIN = fracture_toughness(ilayer) / (Young_mod * &
          & flaw_defect_factor * &
          & sqrt(pi * flaw_oxide_size_factor * Thickness))
 
      omega = Thickness
      ierr = 1

      const_damage_map(ierr, ilayer) = fracture_toughness(ilayer) / &
          & (Young_mod * flaw_defect_factor * &
          & sqrt(pi * flaw_oxide_size_factor))

    else if (TRIM(mode_failure) == 'deflection_delamination')  then

      CRITICAL_STRAIN = 2.0 * fracture_toughness(ilayer) / (Young_mod * &
          & flaw_defect_factor * &
          & sqrt(pi * flaw_oxide_size_factor * Thickness))
 
      omega = Thickness
      ierr = 2

      const_damage_map(ierr, ilayer) = 2.0 * fracture_toughness(ilayer) / &
          & (Young_mod * flaw_defect_factor * &
          & sqrt(pi * flaw_oxide_size_factor))

    else if (TRIM(mode_failure) == 'intefacial_crack_growth')  then

      CRITICAL_STRAIN = - 0.5 * fracture_toughness(ilayer) * &
         & (1.0 + poisson_ratio(ilayer)) / (Young_mod * &
         & flaw_defect_factor * &
         & sqrt(pi * flaw_oxide_size_factor * Thickness))
      omega = Thickness

      ierr = 3

      ! rel. 12, Nagl end Evans (1993)
      ! proposed by Armitt et al (1978)

      const_damage_map(ierr, ilayer) = -0.5 * fracture_toughness(ilayer) * &
         & (1.0 + poisson_ratio(ilayer)) / (Young_mod * &
         & flaw_defect_factor * &
         & sqrt(pi * flaw_oxide_size_factor))

    else if (TRIM(mode_failure) == 'buckling')  then

      CRITICAL_STRAIN = - 1.22 * (Thickness / &
         & delamination_size_now)**2 / (1.0 - poisson_ratio(ilayer)**2)
      omega = Thickness
      ierr = 4

      ! rel. 12, Nagl end Evans (1993)
      ! proposed by Armitt et al (1978)

      const_damage_map(ierr, ilayer) = -1.22 * (1.0 / &
         & delamination_size_now)**2 / (1.0 - poisson_ratio(ilayer)**2)

    else if (TRIM(mode_failure) == 'crack_deflection')  then

      CRITICAL_STRAIN = - 3.6 * (Thickness / &
         & delamination_size_now)**2
      omega = Thickness
      ierr = 5
      ! rel. 12, Nagl end Evans (1993)
      ! proposed by Armitt et al (1978)

      const_damage_map(ierr, ilayer) = - 3.6 * (1.0 / &
         & delamination_size_now)**2

    else if (TRIM(mode_failure) == 'spalling')  then

      CRITICAL_STRAIN = - sqrt(2.0 * surf_fracture_energy(ilayer) / &
        & (Thickness * Young_mod * &
        & (1.0 - poisson_ratio(ilayer))))

      ! check units for the 0.1 factor
      ! CRITICAL_STRAIN = - sqrt(2.0 * surf_fracture_energy(ilayer) * &
      !  & (1.0 + 0.1 * poisson_ratio(ilayer) * Young_mod / &
      !  & (2.0 * surf_fracture_energy(ilayer)**2))  / &
      !  & (Thickness * Young_mod * &
      !  & (1.0 - poisson_ratio(ilayer))))
      omega = Thickness
      ierr = 6

      const_damage_map(ierr, ilayer) = - sqrt(2.0 * surf_fracture_energy(ilayer) / &
        & (Young_mod * (1.0 - poisson_ratio(ilayer))))

    else

      write(out_lun, *) 'ERROR: stress_mode not available'
      stop

    endif

    eta = CRITICAL_STRAIN

   return

  END SUBROUTINE CRITICAL_ETA_OMEGA1

  SUBROUTINE CRITICAL_ETA_OMEGA2(ilayer, mode_failure, &
       & critical_strain, omega, Thickness, Temperature)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the critical strain for each failure mode
    !  Level II of EPRI diagram using:
    !  eta - material properties and geometry effects
    !  omega - damage kinetics parameter dependent on the operation history
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, pi
    use output_module,      only: tty_lun, out_lun, aux_lun
    use oxide_data_module,  only: Youngs_modul, Youngs_temp, nyoungs, &
          & poisson_ratio, fracture_toughness, flaw_defect_factor, &
          & flaw_oxide_size_factor, oxide_thickness, surf_fracture_energy, &
          & amplitude_roughness, &
          & wavelength_roughness, delamination_size_existent, &
          & separated_area, porosity_fraction, &
          & number_cracks_volumetric, delamination_size_factor
  use oxide_data_module,       only: flaw_size_init, flaw_name, nflaws, &
      & flaw_id_burried_crack, flaw_id_pore, flaw_id_interf_crack, &
      & flaw_id_surf_crack, flaw_id_delamin
    use solver_data_module, only: no_thick_interval, no_temp_strain
    use property_module,    only: LINEAR_PROPERTY_SIMPLE

    ! Argument List
    integer, intent(IN)  :: ilayer
    real,    intent(IN)  :: Temperature, Thickness
    character(LEN = 80),    intent(IN)  :: mode_failure
    real,    intent(OUT) :: critical_strain, omega

    ! Local Variables

    real   :: Young_mod, eta, delamination_size_now, flaw_size
    real, dimension(moxide)  :: var, var_ave
    integer :: k, i, ierr

    ! get the Youngs modulus
    call LINEAR_PROPERTY_SIMPLE (Temperature, nyoungs(ilayer) - 1, 0, moxide, &
        & youngs_modul(ilayer, :), youngs_temp(ilayer, :), Young_mod)

      write(aux_lun, *) 'cs ', ilayer, fracture_toughness(ilayer), Young_mod, &
          & flaw_defect_factor, &
          & pi, flaw_oxide_size_factor, Thickness

    if (delamination_size_factor > 0.0)  then
      delamination_size_now = Thickness * delamination_size_factor
    else
      delamination_size_now = delamination_size_existent
    endif

    ierr = 0

    if (TRIM(mode_failure) == 'through_scale_cracking')  then

      ! rel 11, Nagl end Evans (1993)
      ! this one did not work too well reported by Armitt et al (1978)
  
      ! update flaw_size
      call GET_FLAW_SIZE(flaw_id_burried_crack, &
           & Thickness, Temperature, flaw_size)

      ! flaw_defect_factor = 1 since this is a burried defect
      ! no delamination below the crack
      eta = fracture_toughness(ilayer) / (Young_mod * &
          & flaw_defect_factor * &
          & sqrt(pi * flaw_size_init(flaw_id_burried_crack)))
 
      ! omega = flaw_size_now(flaw_id) / flaw_size_init(flaw_id)
      omega = SQRT(flaw_size / flaw_size_init(flaw_id_burried_crack))
      CRITICAL_STRAIN = eta / omega
      ierr = 1

    else if (TRIM(mode_failure) == 'deflection_delamination')  then

      ! update flaw_size
      call GET_FLAW_SIZE(flaw_id_burried_crack, &
           & Thickness, Temperature, flaw_size)
      ! call GET_FLAW_SIZE('burried_crack', Thickness, Temperature, flaw_size)

      ! flaw_defect_factor = 1 since this is a burried defect
      ! delamination below the burried or internal crack
      eta = 2.0 * fracture_toughness(ilayer) / (Young_mod * &
          & flaw_defect_factor * &
          & sqrt(pi * flaw_size_init(flaw_id_burried_crack)))
 
      omega = SQRT(flaw_size / flaw_size_init(flaw_id_burried_crack))
      CRITICAL_STRAIN = eta / omega
      ierr = 2

    else if (TRIM(mode_failure) == 'intefacial_crack_growth')  then

      ! update flaw_size
      call GET_FLAW_SIZE(flaw_id_interf_crack, &
           & Thickness, Temperature, flaw_size)
      ! call GET_FLAW_SIZE('interface_crack', Thickness, Temperature, flaw_size)

      ! amplitude_roughness = 20 microns for defect size 20 microns and 
      ! scale thickness of 450 microns

      eta = - 0.5 * fracture_toughness(ilayer) * &
         & (1.0 + poisson_ratio(ilayer)) / (Young_mod * &
         & flaw_defect_factor * &
         & sqrt(pi * flaw_size_init(flaw_id_interf_crack)))
      omega = (1.0 + amplitude_roughness / Thickness) * &
            & SQRT(flaw_size / flaw_size_init(flaw_id_interf_crack))
      CRITICAL_STRAIN = eta / omega
      ierr = 3

      ! rel. 12, Nagl end Evans (1993)
      ! proposed by Armitt et al (1978)

    else if (TRIM(mode_failure) == 'buckling')  then

      eta = - 1.22 / (1.0 - poisson_ratio(ilayer)**2)
      omega = (delamination_size_now / Thickness)**2 
      CRITICAL_STRAIN = eta / omega
      ierr = 4

      ! rel. 12, Nagl end Evans (1993)
      ! proposed by Armitt et al (1978)

    else if (TRIM(mode_failure) == 'crack_deflection')  then

      ! before 2008 Michael Schutze's report
      ! CRITICAL_STRAIN = - 3.6 * (Thickness / &
      !  & delamination_size_now)**2
      ! the equations has changed a bit
      eta = - 3.75 / (1.0 - poisson_ratio(ilayer)**2)
      omega = (delamination_size_now / Thickness)**2 
      CRITICAL_STRAIN = eta / omega
      ierr = 5

    else if (TRIM(mode_failure) == 'spalling')  then

      ! before 2008 Michael Schutze's report
      ! CRITICAL_STRAIN = - sqrt(2.0 * surf_fracture_energy(ilayer) / &
      !  & (Thickness * Young_mod * &
      !  & (1.0 - poisson_ratio(ilayer))))

      ! update flaw_size
      call GET_FLAW_SIZE(flaw_id_burried_crack, &
           & Thickness, Temperature, flaw_size)

      ! flaw_defect_factor = 1 since this is a burried defect
      ! no delamination below the crack
      eta = fracture_toughness(ilayer) * SQRT(3.0) / (Young_mod * &
          & flaw_defect_factor * &
          & sqrt(pi * flaw_size_init(flaw_id_burried_crack)))
 
      ! omega = flaw_size_now(flaw_id) / flaw_size_init(flaw_id)
      omega = SQRT(flaw_size / flaw_size_init(flaw_id_burried_crack))
      CRITICAL_STRAIN = eta / omega
      ierr = 6

    else

      write(out_lun, *) 'ERROR: stress_mode not available'
      stop

    endif

   return

  END SUBROUTINE CRITICAL_ETA_OMEGA2

  SUBROUTINE GET_FLAW_SIZE(flaw_id, Thickness, Temperature, flaw_size)

    use oxide_data_module,  only: Youngs_modul, Youngs_temp, nyoungs, &
          & poisson_ratio, fracture_toughness, flaw_defect_factor, &
          & flaw_oxide_size_factor, oxide_thickness, surf_fracture_energy, &
          & amplitude_roughness, &
          & wavelength_roughness, delamination_size_existent, &
          & separated_area, porosity_fraction, &
          & number_cracks_volumetric, delamination_size_factor, &
          & flaw_size_init, flaw_size_now, flaw_name
  use oxide_data_module,       only: flaw_size_init, flaw_name, nflaws, &
      & flaw_id_burried_crack, flaw_id_pore, flaw_id_interf_crack, &
      & flaw_id_surf_crack, flaw_id_delamin

    ! Argument List
    integer, intent(IN)  :: flaw_id
    real, intent(OUT)    :: flaw_size
    real,    intent(IN)  :: Temperature, Thickness

    ! notes on defect sizes M. Schutze 2008 Report
    ! Fig. 5 c_normal = 0 to 20 micron at 0 to 100 h for Ni/NiO 
    ! Fig. 5 c_parallel = 5 micron at 0 to 100 h
    ! Fig. 6a Pore size = 3 to 8 microns at 20 and 100 h for  Ni/NiO at 900 C
    ! Fig. 6b Pore size = 3 to 8 microns at 20 and 100 h for oxide on TiAl at 900 C
    ! Fig. 7 c= 6.20 LOG(t) - 7.13 for Ti50Al (0 to 12.5 microns at 0 and 1,000 h)
    ! Fig. 7 c= 4.5e-3*t +0.68 for Ti50Al2Nb (0 to 5 microns at 0 and 1,000 h)

    if (flaw_id == flaw_id_burried_crack)  then   ! 'burried_crack')  then
    ! if (TRIM(flaw_name(flaw_id)) == 'crack')  then

      flaw_size = flaw_oxide_size_factor * Thickness

    else if (flaw_id == flaw_id_surf_crack)  then  ! 'surface_crack')  then

      flaw_size = flaw_oxide_size_factor * Thickness

    else if (flaw_id == flaw_id_pore)  then    ! 'pore')  then

      flaw_size = flaw_oxide_size_factor * Thickness

    else if (flaw_id == flaw_id_interf_crack) then   ! 'interface_crack')  then

      ! delamination radius to be used with 
      ! mode_failure == 'buckling', 'crack_deflection'
      flaw_size = 10.0 * flaw_oxide_size_factor * Thickness

    else if (flaw_id == flaw_id_delamin)  then   ! 'plane_reference')  then
      ! total size of all physical defects (pores and flaws )
      ! this yields unrealistically low critical strains and is not fully developed
      ! here flaw_size_init(flaw_id) = 1.0 micron = oxide thickness
      ! reference length 500 microns; dotal length of defects 150 microns
      flaw_size = (500.0 / (500.0 - 150.0)) * &
           & Thickness / flaw_size_init(1)  !   flaw_id)

    else

      write(6, *) 'flaw_name not valid; use one of '
      write(6, *) 'crack, pore, delamination, or plane_reference'
      STOP

    endif

    RETURN

  END SUBROUTINE GET_FLAW_SIZE

  SUBROUTINE CRITICAL_ETA_OMEGA(ilayer, mode_failure,  &
       & critical_strain, omega, Thickness, Temperature)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the critical strain for each failure mode
    !  Level II of EPRI diagram using:
    !  eta - material properties and geometry effects
    !  omega - damage kinetics parameter dependent on the operation history
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, pi
    use output_module,      only: tty_lun, out_lun, aux_lun
    use oxide_data_module,  only: Youngs_modul, Youngs_temp, nyoungs, &
          & poisson_ratio, fracture_toughness, flaw_defect_factor, &
          & flaw_oxide_size_factor, oxide_thickness, surf_fracture_energy, &
          & amplitude_roughness, &
          & wavelength_roughness, delamination_size_existent, &
          & separated_area, porosity_fraction, &
          & number_cracks_volumetric, delamination_size_factor
    use solver_data_module, only: no_thick_interval, no_temp_strain
    use property_module,    only: LINEAR_PROPERTY_SIMPLE
    use solver_data_module, only: map_level

    ! Argument List
    integer, intent(IN)  :: ilayer
    real,    intent(IN)  :: Temperature, Thickness
    character(LEN = 80),    intent(IN)  :: mode_failure
    real,    intent(OUT) :: critical_strain, omega

    if (map_level == 1)  then

      call CRITICAL_ETA_OMEGA1(ilayer, mode_failure, &
          & critical_strain, omega, Thickness, Temperature)

    else if (map_level == 2)  then

      call CRITICAL_ETA_OMEGA2(ilayer, mode_failure, &
          & critical_strain, omega, Thickness, Temperature)

    else 

      write(6, *) 'map_level should be 1 or 2 '
      STOP

    endif

    RETURN

  END SUBROUTINE CRITICAL_ETA_OMEGA

END MODULE SPALLATION_MODULE
