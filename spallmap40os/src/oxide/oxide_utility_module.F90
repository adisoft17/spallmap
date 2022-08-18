MODULE OXIDE_UTILITY_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   average properties. 
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: COOLING_STRAIN_PLANE, &  ! STORED_EL_ENERGY_OX,
        & OUTPUT_STRAIN_STRESS_ONE, OUTPUT_STRAIN_STRESS_INT, &
        & OUTPUT_STRAIN_STRESS_F2L
  public :: AVERAGE_PROFILE, SET_POISSON_YOUNGS

 CONTAINS

 SUBROUTINE STORED_EL_ENERGY_OX(ilayer, elastic_energy, &
        & distortion_el_energy)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the elastic strain per unit volume in the oxide layer
    ! 
    !=======================================================================
    use parameter_module,   only: moxide
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, th_exp_temp, &
          & th_exp_coeff, nth_exp
    use solution_data_module, only: npr_st, rad_st, stress_axial, &
          & stress_hoop, stress_rad,  &
          & strain_hoop_el, strain_rad_el, strain_axial_el

    implicit none

    ! Argument List
    integer, intent(IN)  :: ilayer
    real, intent(OUT)    :: elastic_energy, distortion_el_energy

    ! Local Variables

    real   :: r, rp1, delta_energy_hoop, delta_energy_rad, &
     & delta_energy_axial, area, transveral_length, delta_energy_dilat_j, &
     & delta_energy_dilat_jp1, delta_energy_dilat, dilatat_el_energy, &
     & delta_energy
    integer :: k, j
   
    ! Nagl end Evans (1993) ref 10; Steiner and Konys (2006); Hsueh and Evans (1983)

    elastic_energy = 0.0
    dilatat_el_energy = 0.0

    do j = npr_st(ilayer) - 1, 1, -1

      r = rad_st(ilayer, j) 
      rp1 = rad_st(ilayer, j+1)

      delta_energy_hoop = stress_hoop(ilayer, j) * &
          & strain_hoop_el(ilayer, j) * r + &
          & stress_hoop(ilayer, j+1) * strain_hoop_el(ilayer, j+1) * &
          & rp1
      delta_energy_rad = stress_rad(ilayer, j) * &
          & strain_rad_el(ilayer, j) * r + &
          & stress_rad(ilayer, j+1) * strain_rad_el(ilayer, j+1) * &
          & rp1
      delta_energy_axial = stress_axial(ilayer, j) * &
          & strain_axial_el(ilayer, j) * r + &
          & stress_axial(ilayer, j+1) * strain_axial_el(ilayer, j+1) * &
          & rp1

      delta_energy = delta_energy_hoop + &
          & delta_energy_rad + delta_energy_axial

      elastic_energy = elastic_energy + delta_energy * (r - rp1)

      ! vm_stress = VON_MISES(stress_rad(ilayer, j), stress_hoop(ilayer, j), &
      !         & stress_axial(ilayer, j))
      delta_energy_dilat_j = (strain_hoop_el(ilayer, j) + &
          & strain_rad_el(ilayer, j) + strain_hoop_el(ilayer, j)) * &
          & (stress_hoop(ilayer, j) + stress_rad(ilayer, j) + &
          & stress_axial(ilayer, j)) / 3.0
      delta_energy_dilat_jp1 = (strain_hoop_el(ilayer, j+1) + &
          & strain_rad_el(ilayer, j+1) + strain_hoop_el(ilayer, j+1)) * &
          & (stress_hoop(ilayer, j+1) + stress_rad(ilayer, j+1) + &
          & stress_axial(ilayer, j+1)) / 3.0
      delta_energy_dilat = delta_energy_dilat_j * r + &
          & delta_energy_dilat_jp1 * rp1

      dilatat_el_energy = dilatat_el_energy + delta_energy_dilat * (r - rp1)

    end do

    area = rad_st(ilayer, 1)**2 - rad_st(ilayer, npr_st(ilayer))**2
    ! in order to be compatible with Armitt map using strain energy; get
    ! the transversal area; dz is not considered
    transveral_length = 0.5 * (rad_st(ilayer, 1) + &
          & rad_st(ilayer, npr_st(ilayer)))

    elastic_energy = 0.5 * elastic_energy / transveral_length

    dilatat_el_energy = 0.5 * dilatat_el_energy / transveral_length

    distortion_el_energy = elastic_energy - dilatat_el_energy

    return

  END SUBROUTINE STORED_EL_ENERGY_OX
 
  SUBROUTINE STORED_EL_ENERGY_OX_PLATE(ilayer, elastic_energy, &
        & distortion_el_energy)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the elastic strain per unit volume in the oxide layer
    ! 
    !=======================================================================
    use parameter_module,   only: moxide
    use output_module,      only: tty_lun, out_lun, aux_lun
    use oxide_data_module,  only: poisson_ratio, Youngs_modul
    use solution_data_module, only:  npr_st, rad_st, strain_hoop_el_ave

    implicit none

    ! Argument List
    integer, intent(IN)  :: ilayer
    real, intent(OUT)    :: elastic_energy, distortion_el_energy

    ! Local Variables
    real :: thickness, dilatat_el_energy
   
    ! ASD-TDR-63-647 technical report for biaxial state of stress

    ! tau_th_st(1, npr_st(1)) - tau_th_st(i, :)
    ! now following Evans and Schutze E * (tau_metal - tau_oxide)/ (1 - Poiss)
    ! (tau_th_st(1, npr_st(1))-tau_th_st(i, :)) * Youngs(i) / &
    !      & (1.0 - poisson(i))

    elastic_energy = Youngs_modul(ilayer, 1) * (strain_hoop_el_ave(ilayer, 3)**2) / &
                   & (1.0 - poisson_ratio(ilayer))

    dilatat_el_energy = 2.0 * (1.0 - 2.0 * poisson_ratio(ilayer)) / &
                   & (3.0 * (1.0 - poisson_ratio(ilayer)))

    thickness = rad_st(ilayer, 1) - rad_st(ilayer, npr_st(ilayer))

    elastic_energy = elastic_energy * thickness

    dilatat_el_energy = dilatat_el_energy  * thickness

    distortion_el_energy = elastic_energy - dilatat_el_energy

    return

  END SUBROUTINE STORED_EL_ENERGY_OX_PLATE

  SUBROUTINE SET_POISSON_YOUNGS(N, if_ave_scale, poisson, Youngs)
  
    !=======================================================================
    ! Purpose(s):
    !
    !  set average Youngs and POisson ratio
    !=======================================================================

  use parameter_module, only: moxide, mave, id_low, id_high
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, poisson_ratio, Youngs_modul
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

  integer, intent(IN)  :: N
  real, dimension(N),    intent(OUT)  :: poisson, Youngs
  logical, intent(IN)  :: if_ave_scale

  ! Local Variables
  integer  :: i, j, loc1, k
 
  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! Poisson ratio is not temperature dependent
  poisson = poisson_ratio

  ! ignore the variation with temperature for the Youngs modulus
  Youngs(1) = Youngs_modul(1, 1)
  poisson(1) = poisson_ratio(1)

  do i = 2, N
    if (if_ave_scale)  then
      ! in this case N = 2 (tube and one average oxide)
      Youngs(i) = Youngs_modul(no_oxide_ave, 1)   ! only for average
      poisson(i) = poisson_ratio(no_oxide_ave)
    else
      ! this will work in general
      Youngs(i) = Youngs_modul(i, 1) 
      poisson(i) = poisson_ratio(i)
    endif
  end do

  RETURN

  END SUBROUTINE SET_POISSON_YOUNGS 

  SUBROUTINE AVERAGE_PROFILE(N, time_now, id_high_or_low, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the average stress-strain profiles
    !  
    !=======================================================================

  use parameter_module, only: moxide, mave, id_low, id_high, mgrid
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer
  use solution_data_module, only: rad_temp, Temp, npr_temp, &
      & npr_st, rad_st, eps0, rad_int, stress_hoop_th_ave, &
      & strain_th_plate_ave, &
      & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
      & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
      & strain_rad, strain_hoop, displ_rad, strain_hoop_gen_ave, &
      & stress_axial, stress_hoop, stress_rad, strain_th_ave, tau_th_st, &
      & strain_hoop_gen, strain_hoop_el, strain_hoop_el_ave
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, if_cylindrical, no_st_rad_points_tube, &
      & no_st_rad_points_oxide
  use property_module, only: OPERAND_ARRAY, OPERAND_2LAYER_ARRAY
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

  ! arguments
  integer, intent(IN)  :: N, id_high_or_low
  logical, intent(IN)  :: if_ave_scale
  real, intent(IN)     :: time_now

  ! Local Variables
  integer  :: i, j, loc1, k, i23
  real, dimension(N+1) :: r
  real, dimension(N)   :: poisson, Youngs
  real  :: htc_in, htc_out, poisson_star, Youngs_star, poisson_star1, Youngs_star1
  real, dimension(N)   :: bc, cc
  character(LEN = 80), dimension(mave), save    :: operation_type

  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! initialize the operation type
  operation_type(1) = 'max'
  operation_type(2) = 'min'
  operation_type(3) = 'ave_sum'

  ! also good for average layer
  call SET_POISSON_YOUNGS(N, if_ave_scale, poisson, Youngs)

  do i = 1, N

    ! introduce star quantities to deal with hoop strain that generates stress
    ! very efficiently as (S_hoop - nu^* S_r) / E^* rather than invoving
    ! the growth and creep terms
    poisson_star = poisson(i) / (1.0 - poisson(i))
    Youngs_star = Youngs(i) / (1.0 - poisson(i)**2)

    do k = 1, 3

      ! get maximum, minimum, and average of thermal strain
      !  k =   1,       2,      3
      ! in each layer
      call OPERAND_ARRAY(1, mgrid, tau_th_st(i, :), 1, npr_st(i), &
       & operation_type(k), strain_th_ave(i, k), loc1)
   
      ! get maximum, minimum, and thermal hoop stress
      ! this would be valid for flat plate geometry
      ! following Townsend et al (1987); JAP 
      ! call OPERAND_ARRAY(1, mgrid, -tau_th_st(i, :) * Youngs(i) / &
      !    & (1.0 - poisson(i)), 1, npr_st(i), &
      !    & operation_type(k), stress_hoop_th_ave(i, k), loc1)

      ! get maximum, minimum, and average stress
      call OPERAND_ARRAY(1, mgrid, stress_rad(i, :), 1, npr_st(i), &
       & operation_type(k), stress_rad_ave(i, k), loc1)

      call OPERAND_ARRAY(1, mgrid, stress_hoop(i, :), 1, npr_st(i), &
       & operation_type(k), stress_hoop_ave(i, k), loc1)

      ! get maximum, minimum, and average strains
      call OPERAND_ARRAY(1, mgrid, strain_rad(i, :), 1, npr_st(i), &
       & operation_type(k), strain_rad_ave(i, k), loc1)

      call OPERAND_ARRAY(1, mgrid, strain_hoop(i, :), 1, npr_st(i), &
       & operation_type(k), strain_hoop_ave(i, k), loc1)

      ! get the equivalent armitt strains such that the damage map may be 
      ! directly compare the critical and strain evolution

   ! strain_hoop_gen_ave = total_strain - thermal_free_stress_strain
   ! the above is for general stress strain equations
   ! for the plain strain we have to use 
   ! strain_hoop_gen_ave = total_strain - thermal_free_stress_strain * (1+poisson)
   ! without ox induced strain and creep
   ! strain_hoop_gen = strain_hoop(i, :) - &
   !      & tau_th_st(i, :) * (1.0 + poisson(i)) + poisson(i) * eps0
    ! introduce star quantities to deal with hoop strain that generates stress
    ! very efficiently as (S_hoop - nu^* S_r) / E^* rather than invoving
    ! the growth and creep terms
      call OPERAND_ARRAY(1, mgrid, (stress_hoop(i, :) - &
         & poisson_star * stress_rad(i, :)) / Youngs_star, &
         & 1, npr_st(i), operation_type(k), &
         & strain_hoop_gen_ave(i, k), loc1)

      call OPERAND_ARRAY(1, mgrid, strain_hoop_el(i, :), 1, npr_st(i), &
       & operation_type(k), strain_hoop_el_ave(i, k), loc1)

      ! get maximum, minimum, and thermal strain
      ! this would be valid for flat plate geometry
      call OPERAND_ARRAY(1, mgrid, &
       & tau_th_st(1, npr_st(1)) - tau_th_st(i, :), 1, npr_st(i), &
       & operation_type(k), strain_th_plate_ave(i, k), loc1)

      ! get maximum, minimum, and thermal hoop stress
      ! this would be valid for flat plate geometry
      ! following Townsend et al (1987); JAP 
      ! use -E (tau_metal - tau_oxide)/ (1 - Poiss)
      ! now following Evans and Schutze E * (tau_metal - tau_oxide)/ (1 - Poiss)
      call OPERAND_ARRAY(1, mgrid, &
          & (tau_th_st(1, npr_st(1))-tau_th_st(i, :)) * Youngs(i) / &
          & (1.0 - poisson(i)), 1, npr_st(i), &
          & operation_type(k), stress_hoop_th_ave(i, k), loc1)

    end do

  end do

  if (N == 4)  then   ! three oxide layers

    i = N-1  ! = 3
    ! introduce star quantities to deal with hoop strain that generates stress
    ! very efficiently as (S_hoop - nu^* S_r) / E^* rather than invoving
    ! the growth and creep terms
    poisson_star = poisson(i) / (1.0 - poisson(i))
    Youngs_star = Youngs(i) / (1.0 - poisson(i)**2)

    i = N  ! = 4
    poisson_star1 = poisson(i) / (1.0 - poisson(i))
    Youngs_star1 = Youngs(i) / (1.0 - poisson(i)**2)

    i23 = N+1
    if (i23 > moxide)  then 
      write(6, *) 'ERROR: i23 > moxide'
      STOP
    endif

    do k = 1, 3

      ! get maximum, minimum, and average of thermal strain
      !  k =   1,       2,      3
      call OPERAND_2LAYER_ARRAY(1, mgrid, tau_th_st(i-1, :), tau_th_st(i, :), &
       & 1, npr_st(i-1), 1, npr_st(i), operation_type(k), &
       & strain_th_ave(i23, k), loc1)
   
      ! get maximum, minimum, and average stress
      call OPERAND_2LAYER_ARRAY(1, mgrid, stress_rad(i-1, :), stress_rad(i, :), &
       & 1, npr_st(i-1), 1, npr_st(i), operation_type(k), &
       & stress_rad_ave(i23, k), loc1)

      call OPERAND_2LAYER_ARRAY(1, mgrid, stress_hoop(i-1, :), stress_hoop(i, :), &
       & 1, npr_st(i-1), 1, npr_st(i), operation_type(k), &
       & stress_hoop_ave(i23, k), loc1)

      ! get maximum, minimum, and average strains
      call OPERAND_2LAYER_ARRAY(1, mgrid, strain_rad(i-1, :), strain_rad(i, :), &
       & 1, npr_st(i-1), 1, npr_st(i), operation_type(k), &
       & strain_rad_ave(i23, k), loc1)

      call OPERAND_2LAYER_ARRAY(1, mgrid, strain_hoop(i-1, :), strain_hoop(i, :), &
       & 1, npr_st(i-1), 1, npr_st(i), operation_type(k), &
       & strain_hoop_ave(i23, k), loc1)

      ! get the equivalent armitt strains such that the damage map may be 
      ! directly compare the critical and strain evolution

      call OPERAND_2LAYER_ARRAY(1, mgrid, (stress_hoop(i-1, :) - &
         & poisson_star * stress_rad(i-1, :)) / Youngs_star, (stress_hoop(i, :) - &
         & poisson_star1 * stress_rad(i, :)) / Youngs_star1, &
         & 1, npr_st(i-1), 1, npr_st(i), operation_type(k), &
         & strain_hoop_gen_ave(i23, k), loc1)

      call OPERAND_2LAYER_ARRAY(1, mgrid, strain_hoop_el(i-1, :), strain_hoop_el(i, :), &
       & 1, npr_st(i-1), 1, npr_st(i), operation_type(k), &
       & strain_hoop_el_ave(i23, k), loc1)

      ! get maximum, minimum, and thermal strain
      ! this would be valid for flat plate geometry
      call OPERAND_2LAYER_ARRAY(1, mgrid, &
       & tau_th_st(1, npr_st(1)) - tau_th_st(i-1, :), &
       & tau_th_st(1, npr_st(1)) - tau_th_st(i, :), &
       & 1, npr_st(i-1), 1, npr_st(i), operation_type(k), &
       & strain_th_plate_ave(i23, k), loc1)

      call OPERAND_2LAYER_ARRAY(1, mgrid, &
          & (tau_th_st(1, npr_st(1))-tau_th_st(i-1, :)) * Youngs(i-1) / &
          & (1.0 - poisson(i-1)), &
          & (tau_th_st(1, npr_st(1))-tau_th_st(i, :)) * Youngs(i) / &
          & (1.0 - poisson(i)), &
          & 1, npr_st(i-1), 1, npr_st(i), operation_type(k), &
       & stress_hoop_th_ave(i23, k), loc1)

    end do

  endif

  if (if_ave_scale)  then
    ! when average scale properties are used only

     write(aux_lun, 2) (stress_rad_ave(i, 1), i=1, N), &
        & (stress_hoop_ave(i, 1), i=1, N), &
        & (strain_rad_ave(i, 1), i=1, N), &
        & (strain_hoop_ave(i, 1), i=1, N)
 2   format('stress_rad_max ', 30(1pe13.6, 1x))

     write(aux_lun, 1) (stress_rad_ave(i, 2), i=1, N), &
        & (stress_hoop_ave(i, 2), i=1, N), &
        & (strain_rad_ave(i, 2), i=1, N), &
        & (strain_hoop_ave(i, 2), i=1, N)
 1   format('stress_rad_min ', 30(1pe13.6, 1x))

     write(aux_lun, 3) (stress_rad_ave(i, 3), i=1, N), &
        & (stress_hoop_ave(i, 3), i=1, N), &
        & (strain_rad_ave(i, 3), i=1, N), &
        & (strain_hoop_ave(i, 3), i=1, N)
 3   format('stress_rad_ave ', 30(1pe13.6, 1x))

    endif

  return

  END SUBROUTINE AVERAGE_PROFILE

  SUBROUTINE PLANAR_GEOM_STRAIN()
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the oxide thickness based on kinetics and temperature schedule 
    ! 
    !=======================================================================
    use parameter_module,   only: mfuel
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & temp_steam_pulse, temp_steam_idle, no_total_pulse, &
          & no_pulse_per_cycle
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, e_th, e_th_mat, &
          & id_mat, no_layer
    use solver_data_module, only: no_thick_interval

    ! Local Variables
    integer  :: i
    real     :: Temp_High, Temp_Low

    if (no_pulse_per_cycle == 1)  then

      if (no_total_pulse > 0)  then
        Temp_High = temp_steam_pulse(1)
        Temp_Low = temp_steam_idle(1)
      else
        Temp_High = MAXVAL(TEMP_STEAM_p, TEMP_STEAM_p > 0.0)
        Temp_Low = MINVAL(TEMP_STEAM_p, TEMP_STEAM_p > 0.0)
      endif

    else if (no_pulse_per_cycle == 2)  then

      ! may need fixing here
      if (no_total_pulse > 0)  then
        Temp_High = temp_steam_pulse(1)
        Temp_Low = temp_steam_idle(1)
      else
        Temp_High = MAXVAL(TEMP_STEAM_p, TEMP_STEAM_p > 0.0)
        Temp_Low = MINVAL(TEMP_STEAM_p, TEMP_STEAM_p > 0.0)
      endif

    endif

    ! get the cooling strain for the layers; substrate is layer 1
    
    do i = 1, no_layer
      e_th_mat(i)= COOLING_STRAIN_PLANE(id_mat(i), Temp_High, Temp_Low)
    end do 

    ! include the substrate effect
    e_th = e_th_mat - e_th_mat(1)

   return

  END SUBROUTINE PLANAR_GEOM_STRAIN

  REAL FUNCTION COOLING_STRAIN_PLANE(ilayer, Temp_High, Temp_Low)
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
    use property_module,    only: LINEAR_PROPERTY_SIMPLE

    ! Argument List
    integer, intent(IN)  :: ilayer
    real,    intent(IN)  :: Temp_High, Temp_Low

    ! Local Variables

    real   :: Temp, delta_temp
    real, dimension(moxide)  :: var, var_ave
    integer :: k, i

    ! eq. 30c, TBB 1987
    ! call One_DIM_INTEGRAL()
    delta_temp = (Temp_High - Temp_Low) / no_temp_strain

    do i = 1, no_temp_strain + 1

      Temp = Temp_Low + delta_temp * (i - 1)

      call LINEAR_PROPERTY_SIMPLE (Temp, nth_exp(ilayer) - 1, 0, moxide, &
        & th_exp_coeff(ilayer, :), th_exp_temp(ilayer, :), Var(i))

    end do

    ! write(out_lun, *) 'check ', (var(i)-varc(i), i = 1, no_temp_strain + 1)

    var_ave(1: no_temp_strain) = 0.5* (var(1:no_temp_strain) + &
         & var(2:no_temp_strain + 1))
    COOLING_STRAIN_PLANE = delta_temp * SUM(var_ave)

    return

  END FUNCTION COOLING_STRAIN_PLANE

  REAL FUNCTION HF_STRAIN_PLANE(ilayer, heat_flux, Temp, oxide_size, tube_size)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the Heat flux strain based on temperature drop and thermal
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
    real,    intent(IN)  :: heat_flux, oxide_size, tube_size, Temp

    ! Local Variables

    real   :: delta_temp
    real, dimension(moxide)  :: var, var_ave
    integer :: k, i

    HF_STRAIN_PLANE = 0.0

    return

  END FUNCTION HF_STRAIN_PLANE
    
  SUBROUTINE OUTPUT_STRAIN_STRESS_ONE(i, thickness_oxide_layer, &
     & id_high_or_low, N, if_ave_scale)
    !=======================================================================
    ! Purpose(s):
    !
    !  PRINT out the strain, stress, temp for one oxide layer
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
          & id_mat, no_layer, material_name, first_oxide_layer
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
           & ntemp_check_strain, &
           & small_oxide_thickness
    use oxide_data_module,   only: id_fe3o4, no_oxide_ave
    use boiler_data_module, only: id_temp_op
    use solution_data_module, only: Temp, Temp_ave, &
      & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
      & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
      & strain_th_ave, stress_hoop_th_ave, npr_temp, stress_rad, &
      & strain_th_plate_ave, eps0, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low, time_full2low, &
      & time_full2low, temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low, &
      & strain_hoop_gen_ave
    use solution_data_module, only: no_oxide_layers

  ! arguments
  integer, intent(IN)  :: N, i, id_high_or_low
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp

    ! Local Variables
    integer  :: j, k, step_cycle, nstart, no_out_f2l
    real     :: Temp_High, Temp_Low, Temp_gas_now, Temp_steam_now, &
              & Press_gas_now, Press_steam_now, slope_time, &
              & slope_temp_gas, slope_temp_steam, slope_press_gas, &
              & slope_press_steam, time_now
    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'
    real, dimension(moxide) :: e_th_scale, thickness_layer
    real, dimension(moxide) :: dtemp, dstrain_th, &
          & dstrain_hoop, dstress_hoop
    real, dimension(mramp_out)  :: press_full_int
    character(LEN = 20)  ::  ch_ave_var, ch_min_max, ch_max_ave_min, ch_map
    character(LEN = 40)  ::  ch_max_ave_min_strain_stres
    logical, save  :: if_map_first = .true.
 
        ! write average temperatures out
  ! mave = 1 - max, 2 - min, 3 - sum average (layer, mave)
  ! strain_th_ave - is the average quantity for the tau 

        dtemp = 0.5 * (Temp_ave(:, 1) - Temp_ave(:, 2))  ! +- the average
        dstrain_th = 0.5 * (strain_th_ave(:, 1) - strain_th_ave(:, 2))
        dstrain_hoop = 0.5 * (strain_hoop_ave(:, 1) - strain_hoop_ave(:, 2))
        dstress_hoop = 0.5 * (stress_hoop_ave(:, 1) - stress_hoop_ave(:, 2))

     if (id_high_or_low == id_high)  then
       
       ch_ave_var = 'temp_ave_delta_hig '
       ch_min_max = 'temp_min_max_hig '
       ch_max_ave_min = 'max_ave_min_hig '
       ch_max_ave_min_strain_stres = 'max_ave_min_hig_eff_strain_gen_stress '
       ch_map = 'strain_ave_map_hig '

     else if (id_high_or_low == id_low)  then

       ch_ave_var = 'temp_ave_delta_low '
       ch_min_max = 'temp_min_max_low '
       ch_max_ave_min = 'max_ave_min_low '
       ch_max_ave_min_strain_stres = 'max_ave_min_low_eff_strain_gen_stress '
       ch_map = 'strain_ave_map_low '

     endif

     if (if_map_first)  then
       ! write the header
       write(aux_lun, 58) ch_map(1:19)
 58     format(a, 'Time t_ox e_h_max')
       write(aux_lun, 59) ch_map(1:15)
 59     format(a, 'low Time t_ox e_h_max')
       if_map_first = .false.
     endif

     write(aux_lun, 57) ch_map, Time_Boiler_p(i), oxide_thickness(i), &
        & strain_hoop_gen_ave(first_oxide_layer, 1)*ten_p3
 57     format(a20, 1x, &
        & 2(1pe13.6, 1x), &
        & 6(1pe11.4, 1x), 1x, 10(1pe10.3, 1x))

     write(aux_lun, 13) ch_ave_var, i, Time_Boiler_p(i), oxide_thickness(i), &
            & Temp_gas_p(i), Temp_steam_p(i), &
            & (Temp_ave(j, 3), dtemp(j), j=1, 2), &
            & (strain_th_ave(j, 3), dstrain_th(j), j=1, 2), &
            & (strain_hoop_ave(j, 3), dstrain_hoop(j), j=1, 2), &
            & (stress_hoop_ave(j, 3), dstress_hoop(j), j=1, 2), &
            & (stress_hoop_th_ave(j, 3), j= 1, 2)
 13     format(a20, i4, 1x, 30(1pe13.6, 1x))
 ! 13     format('temp_ave_hig ', i4, 1x, 30(1pe13.6, 1x))

        write(aux_lun, 15) ch_min_max, i, Time_Boiler_p(i), oxide_thickness(i), &
            & Temp_gas_p(i), Temp_steam_p(i), &
            & (Temp_ave(j, 1), j=1, 2), &
            & (Temp_ave(j, 2), j=1, 2), &
            & (strain_th_ave(j, 1), j=1, 2), &
            & (strain_th_ave(j, 2), j=1, 2), &
            & (strain_hoop_ave(j, 1), j=1, 2), &
            & (strain_hoop_ave(j, 2), j=1, 2), &
            & (stress_hoop_ave(j, 1), j=1, 2), &
            & (stress_hoop_ave(j, 2), j=1, 2), &
            & (stress_rad_ave(j, 1), j=1, 2), &
            & (stress_rad_ave(j, 2), j=1, 2), &
            & (stress_hoop_th_ave(j, 1), j=1, 2), &
            & (stress_hoop_th_ave(j, 2), j=1, 2)

 15     format(a20, i4, 1x, 30(1pe13.6, 1x))
 ! 15     format('temp_min_max_hig ', i4, 1x, 30(1pe13.6, 1x))
 
        write(aux_lun, 16) ch_max_ave_min, i, Time_Boiler_p(i), oxide_thickness(i), &
      & strain_hoop_ave(2, 1)*ten_p3, strain_hoop_ave(2, 3)*ten_p3, &
      & strain_hoop_ave(2, 2)*ten_p3,  &
      & strain_th_ave(2, 1)*ten_p3, strain_th_ave(2, 2)*ten_p3, &
      & strain_th_ave(2, 3)*ten_p3,  &
      & stress_hoop_ave(2, 1), stress_hoop_ave(2, 3), stress_hoop_ave(2, 2), &
      & stress_hoop_th_ave(2, 1), stress_hoop_th_ave(2, 3), &
      & stress_hoop_th_ave(2, 2)

 16     format(a20, i4, 1x, 2(1pe13.6, 1x), &
      & 6(1pe11.4, 1x), 1x, 10(1pe10.3, 1x))
 ! 16     format('max_ave_min_hig ', i4, 1x, 2(1pe13.6, 1x), &
 !     & 6(1pe11.4, 1x), 1x, 10(1pe10.3, 1x))
! 16     format('max_ave_min_hig ', i4, 1x, 2(1pe13.6, 1x), &
!      & 6(F8.4, 1x), 1x, 10(F6.1, 1x))

        write(aux_lun, 56) ch_max_ave_min_strain_stres, i, &
      & Time_Boiler_p(i), oxide_thickness(i), &
      & strain_hoop_gen_ave(2, 1)*ten_p3, strain_hoop_gen_ave(2, 3)*ten_p3, &
      & strain_hoop_gen_ave(2, 2)*ten_p3,  &
      & strain_th_ave(2, 1)*ten_p3, strain_th_ave(2, 2)*ten_p3, &
      & strain_th_ave(2, 3)*ten_p3,  &
      & stress_hoop_ave(2, 1), stress_hoop_ave(2, 3), stress_hoop_ave(2, 2), &
      & stress_hoop_th_ave(2, 1), stress_hoop_th_ave(2, 3), &
      & stress_hoop_th_ave(2, 2)

 56     format(a40, i4, 1x, &
      & 2(1pe13.6, 1x), &
      & 6(1pe11.4, 1x), 1x, 10(1pe10.3, 1x))
 ! 56     format('max_ave_min_hig_eff_strain_gen_stress ', i4, 1x, &
 !     & 2(1pe13.6, 1x), &
 !     & 6(1pe11.4, 1x), 1x, 10(1pe10.3, 1x))

      if (id_high_or_low == id_low .and. if_ave_scale)  then

         write(aux_lun, 64) j, oxide_thickness(i), &
              & strain_th_plate_ave(2, 1)*ten_p3, &
              & strain_th_plate_ave(2, 2)*ten_p3, &
              & strain_th_plate_ave(2, 3)*ten_p3,  &
              & stress_hoop_th_ave(2, 1), stress_hoop_th_ave(2, 3), &
              & stress_hoop_th_ave(2, 2)
 64     format('low_load_hoop_plate_max_min_ave ', i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))

      endif

      if (id_high_or_low == id_high .and. if_ave_scale)  then

   ! strain_hoop_gen_ave = total_strain - thermal_free_stress_strain
        write(aux_lun, 55) i, Time_Boiler_p(i), oxide_thickness(i), &
            & Temp_gas_p(i), Temp_steam_p(i), &
            & (Temp_ave(j, 1), j=1, 2), &
            & (Temp_ave(j, 2), j=1, 2), &
            & (strain_th_ave(j, 1), j=1, 2), &
            & (strain_th_ave(j, 2), j=1, 2), &
            & (strain_hoop_gen_ave(j, 1), j=1, 2), &
            & (strain_hoop_gen_ave(j, 2), j=1, 2), &
            & (stress_hoop_ave(j, 1), j=1, 2), &
            & (stress_hoop_ave(j, 2), j=1, 2), &
            & (stress_rad_ave(j, 1), j=1, 2), &
            & (stress_rad_ave(j, 2), j=1, 2), &
            & (stress_hoop_th_ave(j, 1), j=1, 2), &
            & (stress_hoop_th_ave(j, 2), j=1, 2)

 55     format('temp_min_max_hig_eff_strain_gen_stress ', i4, 1x, &
        & 30(1pe13.6, 1x))

     endif
       ! print out the stress strain state during ramp down
        ! call OUTPUT_FULL2LOW()
    return

  END SUBROUTINE OUTPUT_STRAIN_STRESS_ONE

  SUBROUTINE OUTPUT_STRAIN_STRESS_INT(i, thickness_oxide_layer, &
     & id_high_or_low, N, if_ave_scale, if_outage_ramp)
    !=======================================================================
    ! Purpose(s):
    !
    !  PRINT out the strain, stress, temp at the interfaces metal-flue gas, 
    !    oxide-metal, oxide-oxide, oxide-steam
    !  Temperature, displacement, and strain are continuous across interfaces
    !    and stress (is there a variable that is discontinuous??)
    ! 
    !=======================================================================
    use parameter_module,   only: mstep, moxide, oxide_thick2si, &
          & id_low, id_high, mave, ten_p3, mramp, zero, mramp_out
    use output_module,      only: tty_lun, out_lun, aux_lun, &
          & strn_lun, enrg_lun, strs_lun, gen_lun, enrg_dist_lun, &
          & enrg_dil_lun, strn_gen_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & no_total_pulse, &
          & Temp_gas_p, press_inner_slope_p, press_outer_slope_p, &
          & press_inner_p, press_outer_p
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, e_th, e_th_mat, &
          & id_mat, no_layer, material_name, thickness_fr, &
          & first_oxide_layer, no_layers_oxide
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
           & ntemp_check_strain, &
           & small_oxide_thickness
    use oxide_data_module,   only: id_fe3o4, no_oxide_ave
    use boiler_data_module, only: id_temp_op
    use solution_data_module, only: Temp, Temp_ave, &
      & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
      & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
      & strain_th_ave, stress_hoop_th_ave, npr_temp, stress_rad, &
      & strain_th_plate_ave, eps0, rad_int, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low, time_full2low, &
      & time_full2low, temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low, &
      & strain_hoop_gen_ave
    use solution_data_module, only:  strain_rad, strain_hoop, displ_rad, &
      & strain_hoop_gen_ave, stress_axial, stress_hoop, stress_rad, &
      & tau_th_st, npr_st, strain_hoop_gen, strain_hoop_el_ave, &
      & strain_hoop_el
    use solution_data_module, only: no_oxide_layers, no_f2l_event_out, &
      & strain_elastic_out, energy_elastic_out, temperature_out, &
      & strain_elastic_event_out, temperature_event_out, &
      & energy_elastic_event_out, max_energy_elastic_event, &
      & min_strain_elast_event, max_strain_elast_event, &
      & strain_elastic_out_gen, min_strain_elast_event_gen, &
      & max_strain_elast_event_gen
    use solution_data_module, only: dist_energy_elastic_out, &
      & max_dist_energy_el_event, dist_energy_elastic_out_plate, &
      & energy_elastic_out_plate, strain_elastic_out_plate
    use property_module, only: OPERAND_2LAYER_ARRAY

    implicit none

  ! arguments
  integer, intent(IN)  :: N, i, id_high_or_low
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale, if_outage_ramp
    ! if_outage_ramp = .true. when elastic energy during outage will be printed
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp

    ! Local Variables
    integer, dimension(moxide)  :: max_min
    ! stored_dist_el_energy = distortional energy
    ! 
    real, dimension(N+1)  :: stored_el_energy, stored_dist_el_energy
    real, dimension(N+1)  :: stored_el_energy_plate, stored_dist_el_energy_plate
    integer  :: j, k, step_cycle, nstart, no_out_f2l
    real     :: Temp_High, Temp_Low, Temp_gas_now, Temp_steam_now, &
              & Press_gas_now, Press_steam_now, slope_time, &
              & slope_temp_gas, slope_temp_steam, slope_press_gas, &
              & slope_press_steam, time_now, stored_el_energy_total, &
              & dist_el_energy_tot, stored_el_energy_total_pl, &
              & dist_el_energy_tot_pl
    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'
    character(LEN = 180)    :: ch_map_lay1
    real, dimension(moxide) :: e_th_scale, thickness_layer
    real, dimension(moxide) :: dtemp, dstrain_th, &
          & dstrain_hoop, dstress_hoop
    real, dimension(mramp_out)  :: press_full_int
    character(LEN = 20)  ::  ch_map, ch_map_energy, ch_map_el_strain
    character(LEN = 25)  ::  ch_ave_var, ch_min_max, ch_max_ave_min, &
              & ch_int1, ch_int2, ch_out_thick1, ch_out_thick2, &
          & ch_int2e
    character(LEN = 30)  ::  ch_int2_dw, ch_int2e_dw
    character(LEN = 35)  ::  ch_int3, ch_int4
    character(LEN = 45)  ::  ch_max_ave_min_strain_stres, ch_out_thick3
    character(LEN = 45)  ::  ch_max_ave_min_strn_strs_dv ! double valued
    character(LEN = 80), dimension(20) :: ch_repeat, ch_repeat1
    character(LEN = 1), dimension(9)   :: ch_index = (/'1','2','3','4','5', &
                  & '6','7','8','9'/)
    logical, save  :: if_map_first = .true., if_first_min_max_f_p = .true.
    character(LEN = 80), dimension(mave), save    :: operation_type

        ! write average temperatures out
  ! mave = 1 - max, 2 - min, 3 - sum average (layer, mave)
  ! initialize the operation type
  operation_type(1) = 'max'
  operation_type(2) = 'min'
  operation_type(3) = 'ave_sum'

  ! strain_th_ave - is the average quantity for the tau 

        dtemp = 0.5 * (Temp_ave(:, 1) - Temp_ave(:, 2))  ! +- the average
        dstrain_th = 0.5 * (strain_th_ave(:, 1) - strain_th_ave(:, 2))
        dstrain_hoop = 0.5 * (strain_hoop_ave(:, 1) - strain_hoop_ave(:, 2))
        dstress_hoop = 0.5 * (stress_hoop_ave(:, 1) - stress_hoop_ave(:, 2))

     if (id_high_or_low == id_high)  then
       
       ch_ave_var = 'int_temp_ave_delta_hig '
       ch_min_max = 'int_temp_min_max_hig '
       ch_max_ave_min = 'int_max_ave_min_hig '
       ch_max_ave_min_strain_stres = 'int_max_ave_min_hig_eff_strain_gen_stress '
       ch_max_ave_min_strn_strs_dv = 'int_max_ave_min_hig_eff_strn_gen_strss_dv '

       ch_int1 = 'int_jump_hig '
       ch_int2 = 'int_jump_profile_hig '
       ch_int2e = 'int_jump_es_profile_hig '
       ch_int2_dw = 'int_jump_dw_profile_hig '
       ch_int2e_dw = 'int_jump_dw_es_profile_hig '
       ch_int3 = 'table_layer_prop_hoop_gen_hig '
       ch_int4 = 'table_prop_layer_hoop_gen_hig '

       ch_out_thick1 = 'full_load_temp '
       ch_out_thick2 = 'full_load_hoop '
       ch_out_thick3 = 'full_load_hoop_eff_strain_generate_stress '

       ! lay for layer
       ch_map = 'strain_lay_map_hig '
       ch_map_energy = 'energy_lay_map_hig '
       ch_map_el_strain = 'elstra_lay_map_hig '

     else if (id_high_or_low == id_low)  then

       ch_ave_var = 'int_temp_ave_delta_low '
       ch_min_max = 'int_temp_min_max_low '
       ch_max_ave_min = 'int_max_ave_min_low '
       ch_max_ave_min_strain_stres = 'int_max_ave_min_low_eff_strain_gen_stress '
       ch_max_ave_min_strn_strs_dv = 'int_max_ave_min_low_eff_strn_gen_strss_dv '

       ch_int1 = 'int_jump_low '
       ch_int2 = 'int_jump_profile_low '
       ch_int2e = 'int_jump_es_profile_low '
       ch_int2_dw = 'int_jump_dw_profile_low '
       ch_int2e_dw = 'int_jump_dw_es_profile_low '
       ch_int3 = 'table_layer_prop_hoop_gen_low '  ! p1, p2, p3 per each layer
       ch_int4 = 'table_prop_layer_hoop_gen_low '  ! p1, p1, p1, for all layers

       ch_out_thick1 = 'low_load_temp '
       ch_out_thick2 = 'low_load_hoop '
       ch_out_thick3 = 'low_load_hoop_eff_strain_generate_stress '

       ch_map = 'strain_lay_map_low '
       ch_map_energy = 'energy_lay_map_low '
       ch_map_el_strain = 'elstra_lay_map_low '

     endif

     if (if_map_first)  then
       ! write the header
       ! ch_map_lay1 = ('e_h_max ', j=2, N) 
       ! write(aux_lun, 58) ch_map(1:19), TRIM(ch_map_lay1)

       if (N == 2) then   ! one oxide layer

         write(aux_lun, 58) ch_map(1:19), &
           & ('e_h_max e_h_ave e_h_min ', j=2, N)
 58      format(a, 'Time t_ox ', 10(a))
         write(aux_lun, 59) ch_map(1:15), &
           & ('e_h_max e_h_ave e_h_min ', j=2, N)
 59      format(a, 'low Time t_ox ', 10(a))

         write(aux_lun, 70) ch_map_energy(1:19), &
           & ('e ', j=2, N)
 70      format(a, 'Time t_ox ', 10(a))
         write(aux_lun, 71) ch_map_energy(1:15), &
           & ('e ', j=2, N)
 71      format(a, 'low Time t_ox ', 10(a))

         ! simplified output
         if (strn_lun > 0)  then
           write(strn_lun, 85) 
           write(strn_gen_lun, 85) 
 85        format(' d_f t_f d_p t_p eh1_ave_f eh1_ave_p ', &
            &   'eh1_max_f eh1_max_p eh1_min_f eh1_min_p ', &
            &   'eh1_plat_ave_p eh1_plat_max_p eh1_plat_min_p')

 ! 85        format('d_f d_p t_f t_p eh_max_f eh_ave_f eh_min_f ', &
 !           &   'eh_max_p eh_ave_p eh_min_p')
         endif

         if (enrg_lun > 0)  then
           write(enrg_lun, 86) 
 86        format(' d_f t_f d_p t_p uel_tot_f ', &
            &   'uel_tot_p')
           write(enrg_dist_lun, 101) 
 101        format(' d_f t_f d_p t_p uel_dist_tot_f ', &
            &   'uel_dist_tot_p')
           write(enrg_dil_lun, 102) 
 102        format(' d_f t_f d_p t_p uel_dil_tot_f ', &
            &   'uel_dil_tot_p')
         endif

       else if (N == 3)  then  ! two oxide layers

         write(aux_lun, 60) ch_map(1:19), (j, &
           & ' e_h_max e_h_ave e_h_min ', j=2, N)
 60      format(a, 'Time t_ox_total ', 2('t_ox',i1, 1(a)))
         write(aux_lun, 61) ch_map(1:15), (j, &
           & ' e_h_max e_h_ave e_h_min ', j=2, N)
 61      format(a, 'low Time t_ox_total ', 2('t_ox',i1, 1(a)))

         write(aux_lun, 72) ch_map_energy(1:19), (j, &
           & ' e ', j=2, N)
 72      format(a, 'Time t_ox_total e_tot ', 2('t_ox',i1, 1(a)))
         write(aux_lun, 73) ch_map_energy(1:15), (j, &
           & ' e ', j=2, N)
 73      format(a, 'low Time t_ox_total e_tot ', 2('t_ox',i1, 1(a)))

         write(aux_lun, 80) ch_map_el_strain(1:19), (j, &
           & ' e_h_max e_h_ave e_h_min ', j=2, N)
 80      format(a, 'Time t_ox_total e_h_met_max e_h_met_min ', &
           & 2('t_ox',i1, 1(a)))
         write(aux_lun, 81) ch_map_el_strain(1:15), (j, &
           & ' e_h_max e_h_ave e_h_min ', j=2, N)
 81      format(a, 'low Time t_ox_total e_h_met_max e_h_met_min ', &
           & 2('t_ox',i1, 1(a)))

         ! simplified output
         if (strn_lun > 0)  then
           write(strn_lun, 90) ' d_f t_f d_p t_p eh1_ave_f eh2_ave_f ', &
            &   'eh1_ave_p eh2_ave_p ', &
            &   'eh1_max_f eh2_max_f eh1_max_p eh2_max_p ', &
            &   'eh1_min_f eh2_min_f eh1_min_p eh2_min_p ', &
            &   'eh1_plat_ave_p eh2_plat_ave_p eh1_plat_max_p eh2_plat_max_p ', &
            &   'eh1_plat_min_p eh2_plat_min_p'
           write(strn_gen_lun, 90) ' d_f t_f d_p t_p eh1_ave_f eh2_ave_f ', &
            &   'eh1_ave_p eh2_ave_p ', &
            &   'eh1_max_f eh2_max_f eh1_max_p eh2_max_p ', &
            &   'eh1_min_f eh2_min_f eh1_min_p eh2_min_p'
 90        format(6(a))
         endif

         if (enrg_lun > 0)  then
           write(enrg_lun, 89) 
 89        format(' d_f t_f d_p t_p uel_sp_f uel_mag_f uel_tot_f ', &
            &   'uel_sp_p uel_mag_p uel_tot_p')
           write(enrg_dist_lun, 103) 
 103        format(' d_f t_f d_p t_p uel_dist_sp_f uel_dist_mag_f ', &
            &   'uel_dist_tot_f uel_dist_sp_p uel_dist_mag_p uel_dist_tot_p')
           write(enrg_dil_lun, 104) 
 104        format(' d_f t_f d_p t_p uel_dil_sp_f uel_dil_mag_f ', &
            &   'uel_dil_tot_f uel_dil_sp_p uel_dil_mag_p uel_dil_tot_p')
         endif

       else if (N == 4)  then  ! three oxide layers

         write(aux_lun, 62) ch_map(1:19), (j, &
           & ' e_h_max e_h_ave e_h_min ', j=2, N)
 62      format(a, 'Time t_ox_total ', 3('t_ox',i1, 1(a)))
         write(aux_lun, 63) ch_map(1:15), (j, &
           & ' e_h_max e_h_ave e_h_min ', j=2, N)
 63      format(a, 'low Time t_ox_total ', 3('t_ox',i1, 1(a)))

         write(aux_lun, 74) ch_map_energy(1:19), (j, &
           & ' e ', j=2, N)
 74      format(a, 'Time t_ox_total e_tot', 3('t_ox',i1, 1(a)))
         write(aux_lun, 75) ch_map_energy(1:15), (j, &
           & ' e ', j=2, N)
 75      format(a, 'low Time t_ox_total e_tot', 3('t_ox',i1, 1(a)))

         write(aux_lun, 82) ch_map_el_strain(1:19), (j, &
           & ' e_h_max e_h_ave e_h_min ', j=2, N)
 82      format(a, 'Time t_ox_total e_h_met_max e_h_met_min ', &
           & 3('t_ox',i1, 1(a)))
         write(aux_lun, 83) ch_map_el_strain(1:15), (j, &
           & ' e_h_max e_h_ave e_h_min ', j=2, N)
 83      format(a, 'low Time t_ox_total e_h_met_max e_h_met_min ', &
           & 3('t_ox',i1, 1(a)))

         ! simplified output
         if (strn_lun > 0)  then
           write(strn_lun, 105) &
            &   ' d_f t_f d_p t_p eh1_ave_f eh2_ave_f eh3_ave_f ', &
            &   'eh1_ave_p eh2_ave_p eh3_ave_p eh1_max_f ', &
            &   'eh2_max_f eh3_max_f eh1_max_p eh2_max_p eh3_max_p ', &
            &   'eh1_min_f eh2_min_f eh3_min_f eh1_min_p eh2_min_p eh3_min_p ', &
            &   'eh23_ave_p eh23_max_p eh23_min_p ', & 
            &   'eh1_plat_ave_p eh23_plat_ave_p eh1_plat_max_p eh23_plat_max_p ', &
            &   'eh1_plat_min_p eh23_plat_min_p'
           write(strn_gen_lun, 105) &
            &   ' d_f t_f d_p t_p eh1_ave_f eh2_ave_f eh3_ave_f ', &
            &   'eh1_ave_p eh2_ave_p eh3_ave_p eh1_max_f ', &
            &   'eh2_max_f eh3_max_f eh1_max_p eh2_max_p eh3_max_p ', &
            &   'eh1_min_f eh2_min_f eh3_min_f eh1_min_p eh2_min_p eh3_min_p ', &
            &   'eh23_ave_p eh23_max_p eh23_min_p '
 105        format(40(a))
         endif

         if (enrg_lun > 0)  then
           write(enrg_lun, 106) 
 106        format(' d_f t_f d_p t_p uel_sp_f uel_mag_f uel_maghe_f uel_mag_maghe_f ', &
            &   'uel_tot_f uel_sp_p uel_mag_p uel_maghe_p uel_mag_maghe_p uel_tot_p')
           write(enrg_dist_lun, 107) 
 107        format(' d_f t_f d_p t_p uel_dist_sp_f uel_dist_mag_f ', &
            &   'uel_dist_maghe_f uel_dist_mag_maghe_f uel_dist_tot_f uel_dist_sp_p ', &
            &   'uel_dist_mag_p uel_dist_maghe_p uel_dist_mag_maghe_p uel_dist_tot_p')
           write(enrg_dil_lun, 108) 
 108        format(' d_f t_f d_p t_p uel_dil_sp_f uel_dil_mag_f ', &
            &   'uel_dil_maghe_f uel_dil_mag_maghe_f uel_dil_tot_f uel_dil_sp_p ', &
            &   'uel_dil_mag_p uel_dil_maghe_p uel_dil_mag_maghe_p uel_dil_tot_p')
         endif

       else 

         write(6, *) 'WARNING: more than three oxide layers for spallation'
         ! stop

       endif
       if_map_first = .false.
     endif

     ! checking that thickness_oxide_layer(j) / oxide_thick2si = &
     !    oxide_thickness(i) * thickness_fr(j)
     ! write(6, *) 'check thick ', (thickness_fr(j), j=1, N), &
     !    & (thickness_oxide_layer(j) / oxide_thick2si, j=1, N), &
     !    & (oxide_thickness(i) * thickness_fr(j), j=1, N)

     max_min = 1  ! default show maximum

     ! show maximum in absolute value (keep the sign)
     do j = 2, N

       if (strain_hoop_gen_ave(j, 1) <= 0.0 .or. &
         & (strain_hoop_gen_ave(j, 1) >= 0.0 .and. &
         & strain_hoop_gen_ave(j, 2) < 0.0 .and. &
         & abs(strain_hoop_gen_ave(j, 2)) > strain_hoop_gen_ave(j, 1)))  then
         ! if the maximum is negativ, min < max < 0 and 
         !  |min| > |max|
         ! this concept does not work 7/18/08
         !                            since there will be jumps in the 
         ! curve when the transition occur; better show all max, ave, min
         max_min(j) = 2  
       endif

     end do

     ! write(aux_lun, 57) ch_map, Time_Boiler_p(i), oxide_thickness(i), &
     !   & (oxide_thickness(i) * thickness_fr(j), &
     !   & strain_hoop_gen_ave(j, max_min(j))*ten_p3, j = 2, N)
     write(aux_lun, 57) ch_map, Time_Boiler_p(i), oxide_thickness(i), &
        & (oxide_thickness(i) * thickness_fr(j), &
        & strain_hoop_gen_ave(j, 1)*ten_p3, &
        & strain_hoop_gen_ave(j, 3)*ten_p3, &
        & strain_hoop_gen_ave(j, 2)*ten_p3, j = 2, N)
 57     format(a20, 1x, 2(1pe13.6, 1x), &
        & 6(1pe12.5, 1x), 10(1pe12.5, 1x))

     ! print the elastic strain; not the strain generating stress
     if (id_high_or_low == id_high .or. &
      & (id_high_or_low == id_low .and. (.not. if_outage_ramp)))  then
       write(aux_lun, 84) ch_map_el_strain, Time_Boiler_p(i), &
        & oxide_thickness(i), &
        & strain_hoop_el_ave(1, 1)*ten_p3, strain_hoop_el_ave(1, 2)*ten_p3, &
        & (oxide_thickness(i) * thickness_fr(j), &
        & strain_hoop_el_ave(j, 1)*ten_p3, &
        & strain_hoop_el_ave(j, 3)*ten_p3, &
        & strain_hoop_el_ave(j, 2)*ten_p3, j = 2, N)
 84     format(a20, 1x, 4(1pe13.6, 1x), &
        & 6(1pe12.5, 1x), 10(1pe12.5, 1x))

     else if (id_high_or_low == id_low .and. if_outage_ramp)  then

       ! when outage events took place at id_low; print the highest strain 
       write(aux_lun, 84) ch_map_el_strain, Time_Boiler_p(i), &
        & oxide_thickness(i), &
        & strain_hoop_el_ave(1, 1)*ten_p3, strain_hoop_el_ave(1, 2)*ten_p3, &
        & (oxide_thickness(i) * thickness_fr(j), &
        & strain_hoop_el_ave(j, 1)*ten_p3, &
        & max_strain_elast_event(j), &  ! valid only for first_oxide_layer:N
        & strain_hoop_el_ave(j, 2)*ten_p3, j = 2, N)

       ! when outage events took place at id_low; print the min strain 
       write(aux_lun, 84) ch_map_el_strain, Time_Boiler_p(i), &
        & oxide_thickness(i), &
        & strain_hoop_el_ave(1, 1)*ten_p3, strain_hoop_el_ave(1, 2)*ten_p3, &
        & (oxide_thickness(i) * thickness_fr(j), &
        & strain_hoop_el_ave(j, 1)*ten_p3, &
        & min_strain_elast_event(j), &  ! valid only for first_oxide_layer:N
        & strain_hoop_el_ave(j, 2)*ten_p3, j = 2, N)

     endif

     if (gen_lun > 0)  then

       if (if_first_min_max_f_p)  then

         if (N == 3)  then  ! two oxide layers

           write(gen_lun, 94) 'el_mmfp_strain_outage d_f t_f d_p t_p ', &
            &   'eh1_ave_f eh1_min eh1_max eh1_ave_p ', &
            &   'eh2_ave_f eh2_min eh2_max eh2_ave_p'

           write(gen_lun, 94) 'max_gen_mmfp_strain_outage d_f t_f d_p t_p ', &
            &   'eh1_ave_f eh1_min eh1_max eh1_ave_p ', &
            &   'eh2_ave_f eh2_min eh2_max eh2_ave_p'

           write(gen_lun, 94) 'min_gen_mmfp_strain_outage d_f t_f d_p t_p ', &
            &   'eh1_ave_f eh1_min eh1_max eh1_ave_p ', &
            &   'eh2_ave_f eh2_min eh2_max eh2_ave_p'

           write(gen_lun, 94) 'ave_gen_mmfp_strain_outage d_f t_f d_p t_p ', &
            &   'eh1_ave_f eh1_min eh1_max eh1_ave_p ', &
            &   'eh2_ave_f eh2_min eh2_max eh2_ave_p'

 94        format(3(a))

         else if (N == 4)  then  ! three oxide layers

           write(gen_lun, 105) 'el_mmfp_strain_outage d_f t_f d_p t_p ', &
            &   'eh1_ave_f eh1_min eh1_max eh1_ave_p ', &
            &   'eh2_ave_f eh2_min eh2_max eh2_ave_p', &
            &   'eh3_ave_f eh3_min eh3_max eh3_ave_p'
 
           write(gen_lun, 105) 'max_gen_mmfp_strain_outage d_f t_f d_p t_p ', &
            &   'eh1_ave_f eh1_min eh1_max eh1_ave_p ', &
            &   'eh2_ave_f eh2_min eh2_max eh2_ave_p', &
            &   'eh3_ave_f eh3_min eh3_max eh3_ave_p'

           write(gen_lun, 105) 'min_gen_mmfp_strain_outage d_f t_f d_p t_p ', &
            &   'eh1_ave_f eh1_min eh1_max eh1_ave_p ', &
            &   'eh2_ave_f eh2_min eh2_max eh2_ave_p', &
            &   'eh3_ave_f eh3_min eh3_max eh3_ave_p'

           write(gen_lun, 105) 'ave_gen_mmfp_strain_outage d_f t_f d_p t_p ', &
            &   'eh1_ave_f eh1_min eh1_max eh1_ave_p ', &
            &   'eh2_ave_f eh2_min eh2_max eh2_ave_p', &
            &   'eh3_ave_f eh3_min eh3_max eh3_ave_p'

         endif

         if_first_min_max_f_p = .false.

       endif

       ! write the minimum, maximum, partial and full load so we can construct
       ! the damage map
       strain_elastic_out(id_high_or_low, 1:no_layers_oxide) = &
         & strain_hoop_el_ave(first_oxide_layer:N, 3)*ten_p3 ! ave

       ! strain generating stress
       do j = first_oxide_layer, N
         strain_elastic_out_gen(1:3, id_high_or_low, j-first_oxide_layer+1) = &
         & strain_hoop_gen_ave(j, 1:3)*ten_p3 ! ave
       end do

       if (id_high_or_low == id_low .and. if_outage_ramp)  then
           ! when outage events took place at id_low; print the highest strain 

           ! & (strain_elastic_event_out(id_outage, 1, j), &
         write(gen_lun, 93) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & (strain_elastic_out(id_high, j-first_oxide_layer+1), &
           & min_strain_elast_event(j), max_strain_elast_event(j), &
           & strain_elastic_out(id_low, j-first_oxide_layer+1), &
           & j=first_oxide_layer,N)

 93      format('el_mmfp_strain_outage ', 50(1pe12.5, 1x))

         ! do k = 1, 3  ! max, min, ave
         k = 1
         write(gen_lun, 95) oxide_thickness(i-1), Time_Boiler_p(i-1), &
             & oxide_thickness(i), Time_Boiler_p(i), &
       & (strain_elastic_out_gen(k, id_high, j-first_oxide_layer+1), &
       & min_strain_elast_event_gen(k, j), &
       & max_strain_elast_event_gen(k, j), &
       & strain_elastic_out_gen(k, id_low, j-first_oxide_layer+1), &
             & j=first_oxide_layer,N)
         k = 2
         write(gen_lun, 96) oxide_thickness(i-1), Time_Boiler_p(i-1), &
             & oxide_thickness(i), Time_Boiler_p(i), &
       & (strain_elastic_out_gen(k, id_high, j-first_oxide_layer+1), &
       & min_strain_elast_event_gen(k, j), &
       & max_strain_elast_event_gen(k, j), &
       & strain_elastic_out_gen(k, id_low, j-first_oxide_layer+1), &
             & j=first_oxide_layer,N)
         k = 3
         write(gen_lun, 97) oxide_thickness(i-1), Time_Boiler_p(i-1), &
             & oxide_thickness(i), Time_Boiler_p(i), &
       & (strain_elastic_out_gen(k, id_high, j-first_oxide_layer+1), &
       & min_strain_elast_event_gen(k, j), &
       & max_strain_elast_event_gen(k, j), &
       & strain_elastic_out_gen(k, id_low, j-first_oxide_layer+1), &
             & j=first_oxide_layer,N)
 95      format('max_gen_mmfp_strain_outage ', 50(1pe12.5, 1x))
 96      format('min_gen_mmfp_strain_outage ', 50(1pe12.5, 1x))
 97      format('ave_gen_mmfp_strain_outage ', 50(1pe12.5, 1x))

         ! end do

       endif

     endif

     if (strn_lun > 0)  then


       ! write only after the low load was completed
       ! write the data with/without outage events

       strain_elastic_out(id_high_or_low, 1:no_layers_oxide) = &
         & strain_hoop_el_ave(first_oxide_layer:N, 3)*ten_p3 ! ave
       ! strain_elastic_out(1, 2*no_layers_oxide+1:3*no_layers_oxide) = &
       strain_elastic_out_plate(id_high_or_low, 1:no_layers_oxide) = &
         & strain_th_plate_ave(first_oxide_layer:N, 3)*ten_p3 ! ave

       if (N >= 4)  then
         ! N+1 oxide layer store the average of the layers [2&3 or 2+3+4 .. etc]
         strain_elastic_out(id_high_or_low, no_layers_oxide+1) = &
           & strain_hoop_el_ave(N+1, 3)*ten_p3 ! ave
         strain_elastic_out_plate(id_high_or_low, no_layers_oxide+1) = &
           & strain_th_plate_ave(N+1, 3)*ten_p3
       endif

       strain_elastic_out(id_high_or_low, &
         & 1*(no_layers_oxide+1)+1:2*(no_layers_oxide+1)-1) = &
         & strain_hoop_el_ave(first_oxide_layer:N, 1)*ten_p3 ! max
           ! strain_elastic_out(1, 4*no_layers_oxide+1:5*no_layers_oxide) = &
       strain_elastic_out_plate(id_high_or_low, &
         & 1*(no_layers_oxide+1)+1:2*(no_layers_oxide+1)-1) = &
         & strain_th_plate_ave(first_oxide_layer:N, 1)*ten_p3 ! max

       if (N >= 4)  then
         ! N+1 oxide layer store the average of the layers [2&3 or 2+3+4 .. etc]
         strain_elastic_out(id_high_or_low, 2*(no_layers_oxide+1)) = &
           & strain_hoop_el_ave(N+1, 1)*ten_p3  ! max
         strain_elastic_out_plate(id_high_or_low, 2*(no_layers_oxide+1)) = &
           & strain_th_plate_ave(N+1, 1)*ten_p3
       endif

       strain_elastic_out(id_high_or_low, &
         & 2*(no_layers_oxide+1)+1:3*(no_layers_oxide+1)-1) = &
         & strain_hoop_el_ave(first_oxide_layer:N, 2)*ten_p3 ! min
       strain_elastic_out_plate(id_high_or_low, &
         & 2*(no_layers_oxide+1)+1:3*(no_layers_oxide+1)-1) = &
         & strain_th_plate_ave(first_oxide_layer:N, 2)*ten_p3 ! min

       if (N >= 4)  then
         ! N+1 oxide layer store the average of the layers [2&3 or 2+3+4 .. etc]
         strain_elastic_out(id_high_or_low, 3*(no_layers_oxide+1)) = &
           & strain_hoop_el_ave(N+1, 2)*ten_p3 ! min
         strain_elastic_out_plate(id_high_or_low, 3*(no_layers_oxide+1)) = &
           & strain_th_plate_ave(N+1, 2)*ten_p3
       endif

       ! write the minimum, maximum, partial and full load so we can construct
       ! the damage map
       ! strain generating stress
       do j = first_oxide_layer, N
         strain_elastic_out_gen(1:3, id_high_or_low, j-first_oxide_layer+1) = &
         & strain_hoop_gen_ave(j, 1:3)*ten_p3  ! max, min, ave
       end do

       if (id_high_or_low == id_low)  then

         if (if_outage_ramp)  then
           ! when outage events took place at id_low; print the highest strain 

       if (N < 4)  then

         write(strn_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & strain_elastic_out(id_high, 1:no_layers_oxide), &  ! ave_f
           & max_strain_elast_event(first_oxide_layer:N), &     ! ave_p
     & strain_elastic_out(id_high, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &  ! max_f
     & strain_elastic_out(id_low, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &   ! max_p
     & strain_elastic_out(id_high, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &  ! min_f
     & strain_elastic_out(id_low, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &   ! min_p
     & strain_elastic_out_plate(id_low, 1:no_layers_oxide), &    ! ave_p plate
     & strain_elastic_out_plate(id_low, 1*(no_layers_oxide+1)+1:2*(no_layers_oxide+1)-1), &    ! max_p plate
     & strain_elastic_out_plate(id_low, 2*(no_layers_oxide+1)+1:3*(no_layers_oxide+1)-1)  ! min_p plate

       else if (N >= 4)  then

         write(strn_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & strain_elastic_out(id_high, 1:no_layers_oxide), &  ! ave_f
           & max_strain_elast_event(first_oxide_layer:N), &     ! ave_p
     & strain_elastic_out(id_high, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &  ! max_f
     & strain_elastic_out(id_low, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &   ! max_p
     & strain_elastic_out(id_high, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &  ! min_f
     & strain_elastic_out(id_low, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &   ! min_p
     & (strain_elastic_out(id_low, k*(no_layers_oxide+1)), k = 1, 3), &  ! 23 ave_max_min_p
     & strain_elastic_out_plate(id_low, 1), &                         ! ave_p plate sp
     & strain_elastic_out_plate(id_low, no_layers_oxide+1), &         ! ave_p plate 23
     & strain_elastic_out_plate(id_low, no_layers_oxide+2), &         ! max_p plate sp
     & strain_elastic_out_plate(id_low, 2*(no_layers_oxide+1)), &     ! max_p plate 23
     & strain_elastic_out_plate(id_low, 2*no_layers_oxide+3), &       ! min_p plate sp
     & strain_elastic_out_plate(id_low, 3*(no_layers_oxide+1))        ! min_p plate 23

        endif

         ! write maximum strain that "generates stress"
         write(strn_gen_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & strain_elastic_out_gen(3, id_high, 1:no_layers_oxide), &
           & max_strain_elast_event_gen(3, first_oxide_layer:N), &
     & strain_elastic_out_gen(1, id_high, 1:no_layers_oxide), &
           & max_strain_elast_event_gen(1, first_oxide_layer:N), &
     & strain_elastic_out_gen(2, id_high, 1:no_layers_oxide), &
           & max_strain_elast_event_gen(2, first_oxide_layer:N)

           ! when outage events took place at id_low; print the min strain 

       if (N < 4)  then

         write(strn_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & strain_elastic_out(id_high, 1:no_layers_oxide), &  ! ave_f
           & min_strain_elast_event(first_oxide_layer:N), &     ! ave_p
     & strain_elastic_out(id_high, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &  ! max_f
     & strain_elastic_out(id_low, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &   ! max_p
     & strain_elastic_out(id_high, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &  ! min_f
     & strain_elastic_out(id_low, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &   ! min_p
     & strain_elastic_out_plate(id_low, 1:no_layers_oxide), &    ! ave_p plate
     & strain_elastic_out_plate(id_low, 1*(no_layers_oxide+1)+1:2*(no_layers_oxide+1)-1), &    ! max_p plate
     & strain_elastic_out_plate(id_low, 2*(no_layers_oxide+1)+1:3*(no_layers_oxide+1)-1)  ! min_p plate

       else if (N >= 4)  then

         write(strn_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & strain_elastic_out(id_high, 1:no_layers_oxide), &  ! ave_f
           & min_strain_elast_event(first_oxide_layer:N), &     ! ave_p
     & strain_elastic_out(id_high, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &  ! max_f
     & strain_elastic_out(id_low, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &   ! max_p
     & strain_elastic_out(id_high, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &  ! min_f
     & strain_elastic_out(id_low, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &   ! min_p
     & (strain_elastic_out(id_low, k*(no_layers_oxide+1)), k = 1, 3), &  ! 23 ave_max_min_p
     & strain_elastic_out_plate(id_low, 1), &                         ! ave_p plate sp
     & strain_elastic_out_plate(id_low, no_layers_oxide+1), &         ! ave_p plate 23
     & strain_elastic_out_plate(id_low, no_layers_oxide+2), &         ! max_p plate sp
     & strain_elastic_out_plate(id_low, 2*(no_layers_oxide+1)), &     ! max_p plate 23
     & strain_elastic_out_plate(id_low, 2*no_layers_oxide+3), &       ! min_p plate sp
     & strain_elastic_out_plate(id_low, 3*(no_layers_oxide+1))        ! min_p plate 23

        endif

         ! write minimum strain that "generates stress"
         write(strn_gen_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & strain_elastic_out_gen(3, id_high, 1:no_layers_oxide), &
           & min_strain_elast_event_gen(3, first_oxide_layer:N), &
     & strain_elastic_out_gen(1, id_high, 1:no_layers_oxide), &
           & min_strain_elast_event_gen(1, first_oxide_layer:N), &
     & strain_elastic_out_gen(2, id_high, 1:no_layers_oxide), &
           & min_strain_elast_event_gen(2, first_oxide_layer:N)

         else if (.not. if_outage_ramp)  then

           ! write the energy at the begining of partial load (end of ramp down)
       if (N < 4)  then

         write(strn_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & strain_elastic_out(id_high, 1:no_layers_oxide), &  ! ave_f
           & strain_elastic_out(id_low, 1:no_layers_oxide), &   ! ave_p
     & strain_elastic_out(id_high, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &  ! max_f
     & strain_elastic_out(id_low, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &   ! max_p
     & strain_elastic_out(id_high, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &  ! min_f
     & strain_elastic_out(id_low, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &   ! min_p
     & strain_elastic_out_plate(id_low, 1:no_layers_oxide), &    ! ave_p plate
     & strain_elastic_out_plate(id_low, 1*(no_layers_oxide+1)+1:2*(no_layers_oxide+1)-1), &    ! max_p plate
     & strain_elastic_out_plate(id_low, 2*(no_layers_oxide+1)+1:3*(no_layers_oxide+1)-1)  ! min_p plate

       else if (N >= 4)  then

         write(strn_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & strain_elastic_out(id_high, 1:no_layers_oxide), &  ! ave_f
           & strain_elastic_out(id_low, 1:no_layers_oxide), &   ! ave_p
     & strain_elastic_out(id_high, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &  ! max_f
     & strain_elastic_out(id_low, no_layers_oxide+2:2*(no_layers_oxide+1)-1), &   ! max_p
     & strain_elastic_out(id_high, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &  ! min_f
     & strain_elastic_out(id_low, 2*no_layers_oxide+3:3*(no_layers_oxide+1)-1), &   ! min_p
     & (strain_elastic_out(id_low, k*(no_layers_oxide+1)), k = 1, 3), &  ! 23 ave_max_min_p
     & strain_elastic_out_plate(id_low, 1), &                         ! ave_p plate sp
     & strain_elastic_out_plate(id_low, no_layers_oxide+1), &         ! ave_p plate 23
     & strain_elastic_out_plate(id_low, no_layers_oxide+2), &         ! max_p plate sp
     & strain_elastic_out_plate(id_low, 2*(no_layers_oxide+1)), &     ! max_p plate 23
     & strain_elastic_out_plate(id_low, 2*no_layers_oxide+3), &       ! min_p plate sp
     & strain_elastic_out_plate(id_low, 3*(no_layers_oxide+1))        ! min_p plate 23

        endif

         ! write the strain that "generates stress"
         write(strn_gen_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
           & oxide_thickness(i), Time_Boiler_p(i), &
           & strain_elastic_out_gen(3, id_high, 1:no_layers_oxide), &
           & strain_elastic_out_gen(3, id_low, 1:no_layers_oxide), &
     & strain_elastic_out_gen(1, id_high, 1:no_layers_oxide), &
     & strain_elastic_out_gen(1, id_low , 1:no_layers_oxide), &
     & strain_elastic_out_gen(2, id_high, 1:no_layers_oxide), &
     & strain_elastic_out_gen(2, id_low,  1:no_layers_oxide)
  92     format(70(1pe12.5, 1x))

         endif

       endif

     endif

     ! obtain the elastic storage energy
     stored_el_energy_total = 0.0
     dist_el_energy_tot = 0.0
     do j = first_oxide_layer, N
       call STORED_EL_ENERGY_OX(j, stored_el_energy(j), &
          & stored_dist_el_energy(j))
       stored_el_energy_total = stored_el_energy_total + stored_el_energy(j)
       dist_el_energy_tot = dist_el_energy_tot + stored_dist_el_energy(j)
     end do

     ! obtain the elastic storage energy for the flat plate assumption
     stored_el_energy_total_pl = 0.0
     dist_el_energy_tot_pl = 0.0
     do j = first_oxide_layer, N
       call STORED_EL_ENERGY_OX_PLATE(j, stored_el_energy_plate(j), &
          & stored_dist_el_energy_plate(j))
       stored_el_energy_total_pl = stored_el_energy_total_pl + stored_el_energy_plate(j)
       dist_el_energy_tot_pl = dist_el_energy_tot_pl + stored_dist_el_energy_plate(j)
     end do

     ! print the stored energy in J/m
     write(aux_lun, 76) ch_map_energy, Time_Boiler_p(i), oxide_thickness(i), &
        & 1.0e+6*stored_el_energy_total, &
        & (oxide_thickness(i) * thickness_fr(j), &
        & 1.0e+6*STORED_EL_ENERGY(j), j = first_oxide_layer, N)
 76     format(a20, 1x, 2(1pe13.6, 1x), &
        & 6(1pe13.6, 1x), 10(1pe13.6, 1x))

     ! print the distortion energy in J/m
     write(aux_lun, 77) ch_map_energy, Time_Boiler_p(i), oxide_thickness(i), &
        & 1.0e+6*dist_el_energy_tot, &
        & (oxide_thickness(i) * thickness_fr(j), &
        & 1.0e+6*stored_dist_el_energy(j), j = first_oxide_layer, N)
 77     format(a20, 1x, 2(1pe13.6, 1x), &
        & 6(1pe13.6, 1x), 10(1pe13.6, 1x))

     if (enrg_lun > 0)  then

       ! write(enrg_lun, 89) 
       ! write only after the low load was completed
       ! write the data with/without outage events

       energy_elastic_out(id_high_or_low, 1:no_layers_oxide) = &
         & 1.0e+6 * STORED_EL_ENERGY(first_oxide_layer:N)
       energy_elastic_out(id_high_or_low, no_layers_oxide + 1) = &
         & 1.0e+6 * stored_el_energy_total
       energy_elastic_out_plate(id_high_or_low, 1:no_layers_oxide) = &
         & 1.0e+6 * STORED_EL_ENERGY_PLATE(first_oxide_layer:N)
       energy_elastic_out_plate(id_high_or_low, no_layers_oxide + 1) = &
         & 1.0e+6 * stored_el_energy_total_pl

       dist_energy_elastic_out(id_high_or_low, 1:no_layers_oxide) = &
         & 1.0e+6 * stored_dist_el_energy(first_oxide_layer:N)
       dist_energy_elastic_out(id_high_or_low, no_layers_oxide + 1) = &
         & 1.0e+6 * dist_el_energy_tot
       dist_energy_elastic_out_plate(id_high_or_low, 1:no_layers_oxide) = &
         & 1.0e+6 * stored_dist_el_energy_plate(first_oxide_layer:N)
       dist_energy_elastic_out_plate(id_high_or_low, no_layers_oxide + 1) = &
         & 1.0e+6 * dist_el_energy_tot_pl

       if (id_high_or_low == id_low)  then

         if (if_outage_ramp)  then
           ! when outage events took place at id_low; print the highest energy  
           write(enrg_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
             & oxide_thickness(i), Time_Boiler_p(i), &
             & energy_elastic_out(id_high, 1:no_layers_oxide + 1), &
             & max_energy_elastic_event(first_oxide_layer:N), &
             & max_energy_elastic_event(1) 

           write(enrg_dist_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
             & oxide_thickness(i), Time_Boiler_p(i), &
             & dist_energy_elastic_out(id_high, 1:no_layers_oxide + 1), &
             & max_dist_energy_el_event(first_oxide_layer:N), &
             & max_dist_energy_el_event(1) 

           write(enrg_dil_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
             & oxide_thickness(i), Time_Boiler_p(i), &
             & energy_elastic_out(id_high, 1:no_layers_oxide + 1) - &
             & dist_energy_elastic_out(id_high, 1:no_layers_oxide + 1), &
             & max_energy_elastic_event(first_oxide_layer:N) - &
             & max_dist_energy_el_event(first_oxide_layer:N), &
             & max_energy_elastic_event(1) - max_dist_energy_el_event(1) 

         else if (.not. if_outage_ramp)  then
           ! write the energy at the begining of partial load (end of ramp down)
           write(enrg_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
             & oxide_thickness(i), Time_Boiler_p(i), &
             & energy_elastic_out(id_high, 1:no_layers_oxide + 1), &
             & energy_elastic_out(id_low, 1:no_layers_oxide + 1)

           write(enrg_dist_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
             & oxide_thickness(i), Time_Boiler_p(i), &
             & dist_energy_elastic_out(id_high, 1:no_layers_oxide + 1), &
             & dist_energy_elastic_out(id_low, 1:no_layers_oxide + 1)

           write(enrg_dil_lun, 92) oxide_thickness(i-1), Time_Boiler_p(i-1), &
             & oxide_thickness(i), Time_Boiler_p(i), &
             & energy_elastic_out(id_high, 1:no_layers_oxide + 1) - &
             & dist_energy_elastic_out(id_high, 1:no_layers_oxide + 1), &
             & energy_elastic_out(id_low, 1:no_layers_oxide + 1) - &
             & dist_energy_elastic_out(id_low, 1:no_layers_oxide + 1)

         endif

       endif

     endif

     write(aux_lun, 13) ch_ave_var, i, Time_Boiler_p(i), oxide_thickness(i), &
            & Temp_steam_p(i), &
            & (Temp(j, 1), j=first_oxide_layer, N), Temp(N, npr_temp(N)), &
            &  Temp_gas_p(i), &
            & (tau_th_st(j, 1), j=first_oxide_layer, N), &
            &  tau_th_st(N, npr_st(N)), & ! strain_th
            & (strain_hoop(j, 1), j=first_oxide_layer, N), &
            &  strain_hoop(N, npr_st(N)), &
            & (stress_hoop(j, 1), j=first_oxide_layer, N), &
            &  stress_hoop(N, npr_st(N))
 13     format(a25, i4, 1x, 30(1pe13.6, 1x))
 ! 13     format('temp_ave_hig ', i4, 1x, 30(1pe13.6, 1x))

      ! column 

        write(aux_lun, 56) ch_max_ave_min_strain_stres, i, &
      & Time_Boiler_p(i), oxide_thickness(i), &
      & (strain_hoop_gen(j, 1)*ten_p3, j=first_oxide_layer, N), &
      & strain_hoop_gen(N, npr_st(N))*ten_p3, &
      & (tau_th_st(j, 1)*ten_p3, j=1, N), tau_th_st(N, npr_st(N))*ten_p3, & ! strain_th
      & (stress_hoop(j, 1), j=first_oxide_layer, N), stress_hoop(N, npr_st(N))

        write(aux_lun, 66) ch_max_ave_min_strn_strs_dv, i, &
      & Time_Boiler_p(i), oxide_thickness(i), &
      & (strain_hoop_gen(j, 1)*ten_p3, &
      & strain_hoop_gen(j, npr_st(j))*ten_p3, j=first_oxide_layer, N), &
      & strain_hoop_gen(N, npr_st(N))*ten_p3, &
      & (tau_th_st(j, 1)*ten_p3, j=first_oxide_layer, N), tau_th_st(N, npr_st(N))*ten_p3, & ! strain_th
      & (stress_hoop(j, 1), j=first_oxide_layer, N), stress_hoop(N, npr_st(N))

      ! reporting stress_hoop_th require more work as the formula for 
      ! multiple layers is different than that for one layer

 56     format(a45, i4, 1x, &
      & 2(1pe13.6, 1x), &
      & 6(1pe11.4, 1x), 1x, 40(1pe11.4, 1x))
 66    format(a45, i4, 1x, &
      & 2(1pe13.6, 1x), &
      & 6(1pe15.8, 1x), 1x, 40(1pe15.8, 1x))
 ! 56     format('max_ave_min_hig_eff_strain_gen_stress ', i4, 1x, &
 !     & 2(1pe13.6, 1x), &
 !     & 6(1pe11.4, 1x), 1x, 10(1pe10.3, 1x))

      ! check jump in strain/stress at interfaces
        write(aux_lun, 14) ch_int1, i, &
      & Time_Boiler_p(i), oxide_thickness(i), &
      & (strain_hoop_gen(j-1, npr_st(j-1))*ten_p3, &
      &  strain_hoop_gen(j, 1)*ten_p3, j=first_oxide_layer, N), &
      & (tau_th_st(j-1, npr_st(j-1))*ten_p3, &
      &  tau_th_st(j, 1)*ten_p3, j=first_oxide_layer, N), & 
      & (stress_hoop(j-1, npr_st(j-1)), &
      &  stress_hoop(j, 1), j=first_oxide_layer, N)

 14     format(a45, i4, 1x, 2(1pe13.6, 1x), 30(1pe11.4, 1x))

   do k = 1, no_full2low_output

     ! print hoop stress strain profile at different oxide thicknesses
     if (i == no_point_full2low(k) .or. i - 1 == no_point_full2low(k) )  then

        if (k==1)  then
          write(aux_lun, *) ch_int4, 'id_time Time t_ox T_gas T_steam ', &
            & ('T_max T_min ', j=2, N), &
            & ('e_hg_max e_hg_min ', j=2, N), &
            & ('s_h_max s_h_min ', j=2, N), &
            & ('e_h_max e_h_min ', j=2, N), &
            & ('e_th_max e_th_min ', j=2, N), &
            & ('s_th_max s_th_min ', j=2, N)
        endif

        write(aux_lun, 16) ch_int4, i, Time_Boiler_p(i), oxide_thickness(i), &
            & Temp_gas_p(i), Temp_steam_p(i), &
            & (Temp_ave(j, 1), Temp_ave(j, 2), j=first_oxide_layer, N), &
            & (strain_hoop_gen_ave(j, 1)*ten_p3, &
            & strain_hoop_gen_ave(j, 2)*ten_p3, j=first_oxide_layer, N), &
            & (stress_hoop_ave(j, 1), stress_hoop_ave(j, 2), j=first_oxide_layer, N), &
            & (strain_hoop_ave(j, 1)*ten_p3, strain_hoop_ave(j, 2)*ten_p3, j=first_oxide_layer, N), &
            & (strain_th_ave(j, 1)*ten_p3, strain_th_ave(j, 2)*ten_p3, j=first_oxide_layer, N), &
            & (stress_hoop_th_ave(j, 1), stress_hoop_th_ave(j, 2), j=first_oxide_layer, N)

 16     format(a35, i4, 1x, 2(1pe13.6, 1x), 30(1pe11.4, 1x))

        write(aux_lun, 16) ch_int3, i, Time_Boiler_p(i), oxide_thickness(i), &
            & Temp_gas_p(i), Temp_steam_p(i), &
            & (Temp_ave(j, 1), Temp_ave(j, 2), &
            & strain_hoop_gen_ave(j, 1)*ten_p3, strain_hoop_gen_ave(j, 2)*ten_p3, &
            & stress_hoop_ave(j, 1), stress_hoop_ave(j, 2), &
            & strain_hoop_ave(j, 1)*ten_p3, strain_hoop_ave(j, 2)*ten_p3, &
            & strain_th_ave(j, 1)*ten_p3, strain_th_ave(j, 2)*ten_p3, &
            & stress_hoop_th_ave(j, 1), stress_hoop_th_ave(j, 2), j=first_oxide_layer, N)

        ! get the plot at each given thickness (or given time)
        if (k==1)  then
          ! this is added with awk at the end
          write(aux_lun, 41) ch_int2(1:24)
          write(aux_lun, 41) ch_int2e(1:25)
 41       format(a, 'id_time Time t_ox r_cm e_hg s_h s_r s_z')
        endif

        ! radius is in cm        

        ! elastic hoop strain that generates stress
        write(aux_lun, 15) ch_int2, i, &
      & Time_Boiler_p(i), oxide_thickness(i), &
      & rad_int(1)*100.0, strain_hoop_gen(1, 1)*ten_p3, &
      & stress_hoop(1, 1), stress_rad(1, 1), stress_axial(1, 1)

        ! elastic hoop strain
        write(aux_lun, 15) ch_int2e, i, &
      & Time_Boiler_p(i), oxide_thickness(i), &
      & rad_int(1)*100.0, strain_hoop_el(1, 1)*ten_p3, &
      & stress_hoop(1, 1), stress_rad(1, 1), stress_axial(1, 1)

        do j = 1, N-1
          ! elastic hoop strain that generates stress
          write(aux_lun, 15) ch_int2, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_gen(j, npr_st(j))*ten_p3, &
      &     stress_hoop(j, npr_st(j)), stress_rad(j, npr_st(j)), &
      &     stress_axial(j, npr_st(j))
          write(aux_lun, 15) ch_int2, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_gen(j+1, 1)*ten_p3, &
      &     stress_hoop(j+1, 1), stress_rad(j+1, 1), &
      &     stress_axial(j+1, 1)

          ! elastic hoop strain
          write(aux_lun, 15) ch_int2e, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_el(j, npr_st(j))*ten_p3, &
      &     stress_hoop(j, npr_st(j)), stress_rad(j, npr_st(j)), &
      &     stress_axial(j, npr_st(j))
          write(aux_lun, 15) ch_int2e, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_el(j+1, 1)*ten_p3, &
      &     stress_hoop(j+1, 1), stress_rad(j+1, 1), &
      &     stress_axial(j+1, 1)
        end do

        j = N
        ! elastic hoop strain that generates stress
        write(aux_lun, 15) ch_int2, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_gen(j, npr_st(j))*ten_p3, &
      &     stress_hoop(j, npr_st(j)), stress_rad(j, npr_st(j)), &
      &     stress_axial(j, npr_st(j))
        ! elastic hoop strain
        write(aux_lun, 15) ch_int2e, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_el(j, npr_st(j))*ten_p3, &
      &     stress_hoop(j, npr_st(j)), stress_rad(j, npr_st(j)), &
      &     stress_axial(j, npr_st(j))

 15     format(a45, i4, 1x, 2(1pe13.6, 1x), 30(1pe11.4, 1x))

       ! this output is before and after the main f2l tansition loop
       ! this was the old output
            ! write the output for the paper at the full load
            ! write temperatures
            write(aux_lun, 19) ch_out_thick1, j, oxide_thickness(i), &
              & Temp(1, 1), Temp(2, npr_temp(2)), Temp(2, 1), eps0*ten_p3
 19     format(a25, i4, 1x, 1(1pe13.6, 1x), &
      & 6(1pe10.3, 1x))

            ! write hoop strains and stresses for full load
            write(aux_lun, 20) ch_out_thick2, j, oxide_thickness(i), &
              & strain_hoop_ave(2, 1)*ten_p3, strain_hoop_ave(2, 3)*ten_p3, &
              & strain_hoop_ave(2, 2)*ten_p3, stress_hoop_ave(2, 1), &
              & stress_hoop_ave(2, 3), stress_hoop_ave(2, 2)
 20     format(a25, i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))

   ! strain_hoop_gen_ave  is the effective hoop strain that generates stress
   ! strain_hoop_gen_ave = total_strain - thermal_free_stress_strain
            write(aux_lun, 50) ch_out_thick3, j, oxide_thickness(i), &
        & strain_hoop_gen_ave(2, 1)*ten_p3, strain_hoop_gen_ave(2, 3)*ten_p3, &
        & strain_hoop_gen_ave(2, 2)*ten_p3, stress_hoop_ave(2, 1), &
        & stress_hoop_ave(2, 3), stress_hoop_ave(2, 2)
 50     format(a45, &
                 & i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))


     endif

   enddo

   ! write the profile at daily or weekly cycle, right before the outage
   ! do it only for the last ones
   do k = no_full2low_output - 3, no_full2low_output

     ! print hoop stress strain profile at different oxide thicknesses
     if (i == no_point_full2low(k) -4 .or. &
       & i - 1 == no_point_full2low(k) -4)  then

        ! get the plot at each given thickness (or given time)
        if (k==no_full2low_output - 3)  then
          ! this is added with awk at the end
          write(aux_lun, 41) ch_int2_dw(1:29)
          write(aux_lun, 41) ch_int2e_dw(1:29)

        endif

        ! radius is in cm        

        ! elastic hoop strain that generates stress
        write(aux_lun, 15) ch_int2_dw, i, &
      & Time_Boiler_p(i), oxide_thickness(i), &
      & rad_int(1)*100.0, strain_hoop_gen(1, 1)*ten_p3, &
      & stress_hoop(1, 1), stress_rad(1, 1), stress_axial(1, 1)
        ! elastic hoop strain
        write(aux_lun, 15) ch_int2e_dw, i, &
      & Time_Boiler_p(i), oxide_thickness(i), &
      & rad_int(1)*100.0, strain_hoop_el(1, 1)*ten_p3, &
      & stress_hoop(1, 1), stress_rad(1, 1), stress_axial(1, 1)

        do j = 1, N-1
          ! elastic hoop strain that generates stress
          write(aux_lun, 15) ch_int2_dw, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_gen(j, npr_st(j))*ten_p3, &
      &     stress_hoop(j, npr_st(j)), stress_rad(j, npr_st(j)), &
      &     stress_axial(j, npr_st(j))
          write(aux_lun, 15) ch_int2_dw, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_gen(j+1, 1)*ten_p3, &
      &     stress_hoop(j+1, 1), stress_rad(j+1, 1), &
      &     stress_axial(j+1, 1)

          ! elastic hoop strain
          write(aux_lun, 15) ch_int2e_dw, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_el(j, npr_st(j))*ten_p3, &
      &     stress_hoop(j, npr_st(j)), stress_rad(j, npr_st(j)), &
      &     stress_axial(j, npr_st(j))
          write(aux_lun, 15) ch_int2e_dw, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_el(j+1, 1)*ten_p3, &
      &     stress_hoop(j+1, 1), stress_rad(j+1, 1), &
      &     stress_axial(j+1, 1)
        end do

        j = N
        ! elastic hoop strain that generates stress
        write(aux_lun, 15) ch_int2_dw, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_gen(j, npr_st(j))*ten_p3, &
      &     stress_hoop(j, npr_st(j)), stress_rad(j, npr_st(j)), &
      &     stress_axial(j, npr_st(j))
        ! elastic hoop strain
        write(aux_lun, 15) ch_int2e_dw, i, &
      &     Time_Boiler_p(i), oxide_thickness(i), &
      &     rad_int(j+1)*100.0, strain_hoop_el(j, npr_st(j))*ten_p3, &
      &     stress_hoop(j, npr_st(j)), stress_rad(j, npr_st(j)), &
      &     stress_axial(j, npr_st(j))

     endif

   enddo

      if (id_high_or_low == id_low .and. if_ave_scale)  then

         ! report only the min and max

      endif

      if (id_high_or_low == id_high .and. if_ave_scale)  then


     endif

    return

  END SUBROUTINE OUTPUT_STRAIN_STRESS_INT

  SUBROUTINE OUTPUT_STRAIN_STRESS_F2L(i_oxide_index, f2l_index, &
     & thickness_oxide_layer, time_now, &
     & temp_gas_now, temp_steam_now, N, if_ave_scale, if_new_output)
    !=======================================================================
    ! Purpose(s):
    !
    !  PRINT out the strain, stress, temp at the interfaces metal-flue gas, 
    !    oxide-metal, oxide-oxide, oxide-steam
    !  Temperature, displacement, and strain are continuous across interfaces
    !    and stress (is there a variable that is discontinuous??)
    ! 
    !=======================================================================
    use parameter_module,   only: mstep, moxide, oxide_thick2si, &
          & id_low, id_high, mave, ten_p3, mramp, zero, mramp_out
    use output_module,      only: tty_lun, out_lun, aux_lun, f2l_lun, &
          & gen_lun, f2l_pl_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & no_total_pulse, &
          & Temp_gas_p, press_inner_slope_p, press_outer_slope_p, &
          & press_inner_p, press_outer_p
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, e_th, e_th_mat, &
          & id_mat, no_layer, material_name, first_oxide_layer
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
      & ntemp_check_strain, no_out_f2l, &
      & small_oxide_thickness
    use oxide_data_module,   only: id_fe3o4, no_oxide_ave, &
           & first_oxide_layer, no_metal_layers
    use boiler_data_module, only: id_temp_op
    use solution_data_module, only: Temp, Temp_ave, &
      & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
      & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
      & strain_th_ave, stress_hoop_th_ave, npr_temp, stress_rad, &
      & strain_th_plate_ave, eps0, rad_int, &
      & strain_hoop_gen_ave, no_full2low_output, time_full2low
    use solution_data_module, only:  strain_rad, strain_hoop, displ_rad, &
      & strain_hoop_gen_ave, stress_axial, stress_hoop, stress_rad, &
      & tau_th_st, npr_st, strain_hoop_gen, strain_hoop_el_ave
    use solution_data_module, only: no_oxide_layers, no_f2l_event_out, &
      & strain_elastic_event_out, energy_elastic_event_out, &
      & temperature_event_out, no_pulse_full2low_output, &
      & max_energy_elastic_event, max_energy_elastic_event_out, &
      & min_strain_elast_event, min_strain_elast_event_out, &
      & max_strain_elast_event, max_strain_elast_event_out, &
      & strain_elastic_event_out_g, min_strain_elast_event_gen, &
      & max_strain_elast_event_gen
    use solution_data_module, only: max_dist_energy_el_event, &
      & energy_dist_el_event_out, max_dist_energy_el_event_out
    use shutdown_data_module, only:  time_outage_period, &
         & duration_outage_period, delta_time_outage_period, &
         & temp_steam_outage_period, number_outage_periods, &
         & no_out_f2l_period, no_out_f2l_total, time_now_f2l, &
         & temp_steam_now_f2l, temp_gas_now_f2l, &
         & press_steam_now_f2l, press_gas_now_f2l, no_out_f2l_total
    ! for flat plate assumption
    use solution_data_module, only: dist_energy_elastic_out_plate, &
      & energy_elastic_out_plate, strain_elastic_out_plate, &
      & strain_elastic_event_out_plate, energy_elastic_event_out_plate, &
      & energy_dist_el_event_out_plate, max_energy_elastic_event_out_pl, &
      & max_energy_elastic_event_pl, max_dist_energy_el_event_out_pl, &
      & max_dist_energy_el_event_pl, min_strain_elast_event_pl, &
      & max_strain_elast_event_pl, max_strain_elast_event_out_pl, &
      & min_strain_elast_event_out_pl

  ! arguments
  integer, intent(IN)  :: N, i_oxide_index, f2l_index  ! i - time index; k - f2l transition index
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  real, intent(IN)  :: time_now, temp_gas_now, temp_steam_now
  logical, intent(IN)  :: if_ave_scale, if_new_output

    ! Local Variables
    integer  :: j, step_cycle, nstart, m, i, k
    real     :: Temp_High, Temp_Low, dtime_f2l_now, &
              & Press_gas_now, Press_steam_now, slope_time, &
              & slope_temp_gas, slope_temp_steam, slope_press_gas, &
              & slope_press_steam, stored_el_energy_total, &
              & dist_el_energy_tot, stored_el_energy_total_pl, &
              & dist_el_energy_tot_pl
    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'
    real, dimension(moxide) :: e_th_scale, thickness_layer
    real, dimension(moxide) :: dtemp, dstrain_th, &
          & dstrain_hoop, dstress_hoop
    real, dimension(mramp_out)  :: press_full_int
    character(LEN = 15)  ::  ch_out_f2l
    character(LEN = 16)  ::  ch_out_f2l_el
    character(LEN = 17)  ::  ch_out_f2l_el_en
    character(LEN = 80), dimension(20) :: ch_repeat, ch_repeat1
    character(LEN = 1), dimension(9)   :: ch_index = (/'1','2','3','4','5', &
                  & '6','7','8','9'/)
    logical, save :: legend_first = .true.
    integer, save :: id_now, id_outage = 0

    real, dimension(N)  :: stored_el_energy, stored_dist_el_energy, &
       & stored_el_energy_plate, stored_dist_el_energy_plate

        ! write average temperatures out
  ! mave = 1 - max, 2 - min, 3 - sum average (layer, mave)
  ! strain_th_ave - is the average quantity for the tau 

     if (if_ave_scale)  then

       ! old output
       if (.not. if_new_output)  then

        write(6, *) f2l_index, 60.0 * time_now, oxide_thickness(i_oxide_index), &
      & strain_hoop_ave(2, 1)*ten_p3, strain_hoop_ave(2, 3)*ten_p3, &
      & strain_hoop_ave(2, 2)*ten_p3
  
        write(aux_lun, 18) f2l_index, 60.0 * time_now, &
      & oxide_thickness(i_oxide_index), &
      & strain_hoop_ave(2, 1)*ten_p3, strain_hoop_ave(2, 3)*ten_p3, &
      & strain_hoop_ave(2, 2)*ten_p3,  &
      & strain_th_ave(2, 1)*ten_p3, strain_th_ave(2, 2)*ten_p3, &
      & strain_th_ave(2, 3)*ten_p3,  &
      & stress_hoop_ave(2, 1), stress_hoop_ave(2, 3), stress_hoop_ave(2, 2), &
      & stress_hoop_th_ave(2, 1), stress_hoop_th_ave(2, 3), &
      & stress_hoop_th_ave(2, 2), &
      & Temp(2, npr_temp(2)), Temp(2, 1), temp_gas_now, temp_steam_now, &
      & stress_rad(2, 1)
      ! temp at the oxide surface, oxide-metal interface, pressure at interface
      ! & temp_gas_now, temp_steam_now, press_gas_now, press_steam_now

 18     format('f2l_max_ave_min ', i4, 1x, 2(1pe13.6, 1x), &
      & 6(1pe11.4, 1x), 1x, 16(1pe10.3, 1x))

        write(aux_lun, 58) f2l_index, 60.0 * time_now, &
      & oxide_thickness(i_oxide_index), &
      & strain_hoop_gen_ave(2, 1)*ten_p3, strain_hoop_gen_ave(2, 3)*ten_p3, &
      & strain_hoop_gen_ave(2, 2)*ten_p3,  &
      & strain_th_ave(2, 1)*ten_p3, strain_th_ave(2, 2)*ten_p3, &
      & strain_th_ave(2, 3)*ten_p3,  &
      & stress_hoop_ave(2, 1), stress_hoop_ave(2, 3), stress_hoop_ave(2, 2), &
      & stress_hoop_th_ave(2, 1), stress_hoop_th_ave(2, 3), &
      & stress_hoop_th_ave(2, 2), &
      & Temp(2, npr_temp(2)), Temp(2, 1), stress_rad(2, 1)

 58     format('f2l_max_ave_min_eff_strain_gen_stress ', i4, 1x, &
        & 2(1pe13.6, 1x), &
      & 6(1pe11.4, 1x), 1x, 14(1pe10.3, 1x))

      endif

     else if (.not. if_ave_scale)  then


     endif

     if (if_new_output)  then

        if (legend_first)  then

          ch_out_f2l = 'f2l_fig_data '

          write(aux_lun, *) ch_out_f2l, 'id_f2l time_f2l t_ox ', &
            & ('e_hg_max ', m=first_oxide_layer, N), &
            & ('e_hg_ave ', m=first_oxide_layer, N), &
            & ('e_hg_min ', m=first_oxide_layer, N), &
            & ('s_h_max ', m=first_oxide_layer, N), &
            & ('s_h_ave ', m=first_oxide_layer, N), &
            & ('s_h_min ', m=first_oxide_layer, N)

          ch_out_f2l_el = 'f2l_el_fig_data '

          write(aux_lun, *) ch_out_f2l_el, 'id_f2l time_f2l t_ox ', &
            & 'Tg Tm Tm Tox Tst ', ('e_he_max ', m=first_oxide_layer, N), &
            & ('e_he_ave ', m=first_oxide_layer, N), &
            & ('e_he_min ', m=first_oxide_layer, N), &
            & ('s_h_max ', m=first_oxide_layer, N), &
            & ('s_h_ave ', m=first_oxide_layer, N), &
            & ('s_h_min ', m=first_oxide_layer, N)

          ch_out_f2l_el_en = 'f2l_el_energ_fig '

          write(aux_lun, *) ch_out_f2l_el_en, 'id_f2l time_f2l t_ox ', &
            & 'Tg Tm Tm Tox Tst ', 'enel_tot ', ('enel ', m=first_oxide_layer, N), &
            & ('e_he_ave ', m=first_oxide_layer, N), &
            & ('s_h_ave ', m=first_oxide_layer, N)

          if (N == 2)  then   ! one oxide layer

            ! must be repeaded no_full2low_output times
            if (f2l_lun > 0)  then
              ch_repeat(1) = ' uel_tot eh_ave'
              ch_repeat(2:10) = ch_repeat(1)
              write(f2l_lun, 87) (ch_repeat(j)(1:15), j = 1, no_full2low_output)
 ! 87         format('temp', (' uel_tot eh_ave'))
 87           format('T_gas_metal T_metal_ox T_ox_steam', 10(a))

            endif

          else if (N == 3)  then  ! two oxide layers

            if (f2l_lun > 0)  then
              ! moved down in order to get the time, cycle, oxide thickness info

            endif

          endif

          legend_first = .false.

        endif

       ch_out_f2l = 'f2l_fig_data '

       write(aux_lun, 50) ch_out_f2l, f2l_index, 60.0 * time_now, &
        & oxide_thickness(i_oxide_index), &
        & (strain_hoop_gen_ave(m, 1)*ten_p3, m=first_oxide_layer, N), &
        & (strain_hoop_gen_ave(m, 3)*ten_p3, m=first_oxide_layer, N), &
        & (strain_hoop_gen_ave(m, 2)*ten_p3, m=first_oxide_layer, N), &
        & (stress_hoop_ave(m, 1), m=first_oxide_layer, N), &
        & (stress_hoop_ave(m, 3), m=first_oxide_layer, N), &
        & (stress_hoop_ave(m, 2), m=first_oxide_layer, N)
 50     format(a15, i4, 1x, 2(1pe13.6, 1x), &
                 & 30(1pe11.4, 1x))  ! , 1x, 3(1pe10.3, 1x))

       ch_out_f2l_el = 'f2l_el_fig_data '

     ! print the elastic strain; not the strain generating stress
     write(aux_lun, 84) ch_out_f2l_el, f2l_index, 60.0 * time_now, &
        & oxide_thickness(i_oxide_index), temp_gas_now, Temp(1, 1), &
        & Temp(no_metal_layers, npr_temp(no_metal_layers)), &
        & Temp(N, npr_temp(N)), temp_steam_now, &
        & (strain_hoop_el_ave(m, 1)*ten_p3, m=first_oxide_layer, N), &
        & (strain_hoop_el_ave(m, 3)*ten_p3, m=first_oxide_layer, N), &
        & (strain_hoop_el_ave(m, 2)*ten_p3, m=first_oxide_layer, N), &
        & (stress_hoop_ave(m, 1), m=first_oxide_layer, N), &
        & (stress_hoop_ave(m, 3), m=first_oxide_layer, N), &
        & (stress_hoop_ave(m, 2), m=first_oxide_layer, N)
 84     format(a16, i4, 1x, 2(1pe13.6, 1x), &
                 & 30(1pe11.4, 1x))

     ! obtain the elastic storage energy
     stored_el_energy_total = 0.0
     dist_el_energy_tot = 0.0
     do j = first_oxide_layer, N
       call STORED_EL_ENERGY_OX(j, stored_el_energy(j), &
          & stored_dist_el_energy(j))
       stored_el_energy_total = stored_el_energy_total + stored_el_energy(j)
       dist_el_energy_tot = dist_el_energy_tot + stored_dist_el_energy(j)
     end do

     ! obtain the elastic storage energy for the flat plate assumption
     stored_el_energy_total_pl = 0.0
     dist_el_energy_tot_pl = 0.0
     do j = first_oxide_layer, N
       call STORED_EL_ENERGY_OX_PLATE(j, stored_el_energy_plate(j), &
          & stored_dist_el_energy_plate(j))
       stored_el_energy_total_pl = stored_el_energy_total_pl + stored_el_energy_plate(j)
       dist_el_energy_tot_pl = dist_el_energy_tot_pl + stored_dist_el_energy_plate(j)
     end do

     ch_out_f2l_el_en = 'f2l_el_energ_fig '

     ! print the stored energy in J/m
     write(aux_lun, 76) ch_out_f2l_el_en, f2l_index, 60.0 * time_now, &
        & oxide_thickness(i_oxide_index), temp_gas_now, Temp(1, 1), &
        & Temp(no_metal_layers, npr_temp(no_metal_layers)), &
        & Temp(N, npr_temp(N)), temp_steam_now, &
        & 1.0e+6*stored_el_energy_total, &
        & (1.0e+6*STORED_EL_ENERGY(j), j = first_oxide_layer, N), &
        & (strain_hoop_el_ave(m, 3)*ten_p3, m=first_oxide_layer, N), &
        & (stress_hoop_ave(m, 3), m=first_oxide_layer, N)
 76     format(a17, i4, 1x, 2(1pe13.6, 1x), &
                 & 30(1pe11.4, 1x))

      if (f2l_lun > 0 .and. no_pulse_full2low_output(f2l_index) == 0)  then
        ! store the energy for the output; print at the last time when this is called
        ! identify outage events by the fact that 
        ! no_pulse_full2low_output(f2l_index) = 0 while 
        ! i_oxide_index > 0 - index in the oxide growth

        if (time_now <= 1.0e-4)  then

          ! id_outage - stores the current number of outage
          id_outage = id_outage + 1

          ! id_now - stores the number of step in the current outage
          id_now = 1

        else

         id_now = id_now + 1

        endif

        ! ids are working fine
        ! id_outage == no_f2l_event_out .and. id_now == no_out_f2l
        ! write(6, 3) f2l_index, &
        !  & no_pulse_full2low_output(f2l_index), &
        !  & id_outage, no_f2l_event_out, id_now, no_out_f2l
 3      format('check index ', 10(i3, 1x))

        energy_elastic_event_out(id_outage, id_now, first_oxide_layer:N) = &
          & 1.0e+6*STORED_EL_ENERGY(first_oxide_layer:N)
        energy_elastic_event_out(id_outage, id_now, 1) = &
          & 1.0e+6*stored_el_energy_total

        ! for the flat plate assumption
        energy_elastic_event_out_plate(id_outage, id_now, first_oxide_layer:N) = &
          & 1.0e+6*STORED_EL_ENERGY_PLATE(first_oxide_layer:N)
        energy_elastic_event_out_plate(id_outage, id_now, 1) = &
          & 1.0e+6*stored_el_energy_total_pl

        strain_elastic_event_out(id_outage, id_now, first_oxide_layer:N) = &
          &  strain_hoop_el_ave(first_oxide_layer:N, 3)*ten_p3

        ! for the flat plate assumption
        strain_elastic_event_out_plate(id_outage, id_now, first_oxide_layer:N) = &
          &  strain_th_plate_ave(first_oxide_layer:N, 3)*ten_p3

        ! strain generating stress
        do j = first_oxide_layer, N
          strain_elastic_event_out_g(:, id_now, j) = &
            &  strain_hoop_gen_ave(j, :) * ten_p3
        end do

        temperature_event_out(id_outage, id_now) = &
          &  Temp(1, npr_temp(1))

        if (id_now == no_out_f2l + 1)  then
          ! end of this current outage was now reached

          do k = 1, N

            max_energy_elastic_event(k) = &
         & MAXVAL(energy_elastic_event_out(id_outage, 1:no_out_f2l + 1, k))
              ! & MAXVAL(energy_elastic_event_out(id_outage, :, k))

            max_energy_elastic_event_pl(k) = &
         & MAXVAL(energy_elastic_event_out_plate(id_outage, 1:no_out_f2l + 1, k))

            max_dist_energy_el_event(k) = &
         & MAXVAL(energy_dist_el_event_out(id_outage, 1:no_out_f2l + 1, k))

            max_dist_energy_el_event_pl(k) = &
         & MAXVAL(energy_dist_el_event_out_plate(id_outage, 1:no_out_f2l + 1, k))

            max_strain_elast_event(k) = &  ! first_oxide_layer:N valid only
         & MAXVAL(strain_elastic_event_out(id_outage, 1:no_out_f2l + 1, k))
            min_strain_elast_event(k) = &  
         & MINVAL(strain_elastic_event_out(id_outage, 1:no_out_f2l + 1, k))

            max_strain_elast_event_pl(k) = &  ! first_oxide_layer:N valid only
         & MAXVAL(strain_elastic_event_out_plate(id_outage, 1:no_out_f2l + 1, k))
            min_strain_elast_event_pl(k) = &  
         & MINVAL(strain_elastic_event_out_plate(id_outage, 1:no_out_f2l + 1, k))

            do j = 1, 3
              max_strain_elast_event_gen(j, k) = & 
         & MAXVAL(strain_elastic_event_out_g(j, 1:no_out_f2l + 1, k))
              min_strain_elast_event_gen(j, k) = &  
         & MINVAL(strain_elastic_event_out_g(j, 1:no_out_f2l + 1, k))
            end do

          end do
          ! to be sure
          ! max_energy_elastic_event(2:N) = &
          !  & MAXVAL(energy_elastic_event_out(id_outage, :, 2:N))

          max_energy_elastic_event_out(id_outage, :) = &
            & max_energy_elastic_event(:)
          max_dist_energy_el_event_out(id_outage, :) = &
            & max_dist_energy_el_event(:)
          
          max_energy_elastic_event_out_pl(id_outage, :) = &
            & max_energy_elastic_event_pl(:)
          max_dist_energy_el_event_out_pl(id_outage, :) = &
            & max_dist_energy_el_event_pl(:)

          max_strain_elast_event_out(id_outage, :) = &
            & max_strain_elast_event(:)
          min_strain_elast_event_out(id_outage, :) = &
            & min_strain_elast_event(:)

          max_strain_elast_event_out_pl(id_outage, :) = &
            & max_strain_elast_event_pl(:)
          min_strain_elast_event_out_pl(id_outage, :) = &
            & min_strain_elast_event_pl(:)

          dtime_f2l_now = (time_full2low(f2l_index, mramp) - &
               & time_full2low(f2l_index, 1)) / no_out_f2l
          ! write(tty_lun, 4) f2l_index, id_outage, id_now, &
          write(tty_lun, 4) id_outage, &
               & oxide_thickness(i_oxide_index), &
               & time_full2low(f2l_index, mramp), &
               & dtime_f2l_now * no_out_f2l, dtime_f2l_now
          write(out_lun, 4) id_outage, &
               & oxide_thickness(i_oxide_index), &
               & time_full2low(f2l_index, mramp), &
               & dtime_f2l_now * no_out_f2l, dtime_f2l_now
          ! write(out_lun, 4) f2l_index, id_outage, id_now, &
          !     & dtime_f2l_now, dtime_f2l_now * no_out_f2l
          if (gen_lun > 0)  then
            write(gen_lun, 4) id_outage, &
               & oxide_thickness(i_oxide_index), &
               & time_full2low(f2l_index, mramp), &
               & dtime_f2l_now * no_out_f2l, dtime_f2l_now
          endif

 4        format('outage ', i2, ' at oxide_thickness= ', &
               & 1pe13.6, ' time= ', 1pe13.6, &
               & ' transition duration= ', 1pe13.6, &
               & ' dt= ', 1pe13.6)  ! , 1(i4, 1x))
        endif

        ! make sure that all the f2l transitions are complete before printing, 
        ! including the last transition
        if (id_outage == no_f2l_event_out .and. &
          & ((id_now == no_out_f2l_total .and. no_out_f2l_total > 0) .or. &
          & (id_now == no_out_f2l + 1 .and. no_out_f2l_total == 0))) then

          ! moved down in order to get the time, cycle, oxide thickness info
          ! 'T_gas_metal T_metal_ox T_ox_steam', 10(a))
          ! be careful what temperatures are plotted as they vary with oxide thickness
         ! delete ave since for 6 outages, two lines of characters are written
         !  write(f2l_lun, 91) (' T'//ch_index(j)//'_ox_met', &
         !   & ' uel'//ch_index(j)//'_tot'// &
         !   & ' uel'//ch_index(j)//'_sp'//&
         !   & ' eh'//ch_index(j)//'_sp_ave', &
         !   & ' uel'//ch_index(j)//'_mag',&
         !   & ' eh'//ch_index(j)//'_mag_ave', j = 1, no_f2l_event_out) 
         if (N == 2) then  ! one-oxide layer

           write(f2l_lun, 91) (' T'//ch_index(j)//'_ox_met', &
            & ' uel'//ch_index(j)//'_tot'// &
            & ' eh'//ch_index(j)//'_tot', j = 1, no_f2l_event_out)

           write(f2l_pl_lun, 91) (' T'//ch_index(j)//'_ox_met', &
            & ' uel'//ch_index(j)//'_tot'// &
            & ' eh'//ch_index(j)//'_tot', j = 1, no_f2l_event_out)

         else if (N == 3) then  ! two-oxide layers 

           write(f2l_lun, 91) (' T'//ch_index(j)//'_ox_met', &
            & ' uel'//ch_index(j)//'_tot'// &
            & ' uel'//ch_index(j)//'_sp'//&
            & ' eh'//ch_index(j)//'_sp_ave', &
            & ' uel'//ch_index(j)//'_mag',&
            & ' eh'//ch_index(j)//'_mag_ave', j = 1, no_f2l_event_out) 

           write(f2l_pl_lun, 91) (' T'//ch_index(j)//'_ox_met', &
            & ' uel'//ch_index(j)//'_tot'// &
            & ' uel'//ch_index(j)//'_sp'//&
            & ' eh'//ch_index(j)//'_sp_ave', &
            & ' uel'//ch_index(j)//'_mag',&
            & ' eh'//ch_index(j)//'_mag_ave', j = 1, no_f2l_event_out) 

         else if (N == 4)  then  ! three oxide layers

           write(f2l_lun, 91) (' T'//ch_index(j)//'_ox_met', &
            & ' uel'//ch_index(j)//'_tot'// &
            & ' uel'//ch_index(j)//'_sp'//&
            & ' eh'//ch_index(j)//'_sp_ave', &
            & ' uel'//ch_index(j)//'_mag',&
            & ' eh'//ch_index(j)//'_mag_ave', & 
            & ' uel'//ch_index(j)//'_maghe',&
            & ' eh'//ch_index(j)//'_maghe_ave', &
            & ' uel'//ch_index(j)//'_mag23he',&
            & ' eh'//ch_index(j)//'_mag23he_ave', j = 1, no_f2l_event_out) 

           write(f2l_pl_lun, 91) (' T'//ch_index(j)//'_ox_met', &
            & ' uel'//ch_index(j)//'_tot'// &
            & ' uel'//ch_index(j)//'_sp'//&
            & ' eh'//ch_index(j)//'_sp_ave', &
            & ' uel'//ch_index(j)//'_mag',&
            & ' eh'//ch_index(j)//'_mag_ave', & 
            & ' uel'//ch_index(j)//'_maghe',&
            & ' eh'//ch_index(j)//'_maghe_ave', &
            & ' uel'//ch_index(j)//'_mag23he',&
            & ' eh'//ch_index(j)//'_mag23he_ave', j = 1, no_f2l_event_out) 

         endif
 91      format(' time', 150(a))
         
          if (no_out_f2l_total == 0)  then
            ! this is the old output sometimes the no_out_f2l_total may not be defined
            ! assume that all the events have the same duration
            ! otherwise need to print the time, too
            ! this will work only for simple shutdown, where p(t) and T(t) are the same
            !      and for one time interval
            dtime_f2l_now = (time_full2low(f2l_index, mramp) - &
                    & time_full2low(f2l_index, 1)) / no_out_f2l
          
            do j = 1, no_out_f2l + 1

              call WRITE_EVENT(f2l_lun, N, j, dtime_f2l_now, temperature_event_out, &
                 & energy_elastic_event_out, strain_elastic_event_out)

              ! flat plate assumption or biaxial state of stress
              call WRITE_EVENT(f2l_pl_lun, N, j, dtime_f2l_now, temperature_event_out, &
                 & energy_elastic_event_out_plate, strain_elastic_event_out_plate)

            end do

          else if (no_out_f2l_total > 0)  then

            ! different p(t) and T(t) over several periods can be handled at shut down
            do j = 1, no_out_f2l_total

              call WRITE_EVENT(f2l_lun, N, j, dtime_f2l_now, temperature_event_out, &
                 & energy_elastic_event_out, strain_elastic_event_out)

              ! flat plate assumption or biaxial state of stress
              call WRITE_EVENT(f2l_pl_lun, N, j, dtime_f2l_now, temperature_event_out, &
                 & energy_elastic_event_out_plate, strain_elastic_event_out_plate)

            end do

          endif

        endif

       endif

     endif

    return

  END SUBROUTINE OUTPUT_STRAIN_STRESS_F2L

  SUBROUTINE WRITE_EVENT(unit, N, j, dtime_f2l_now, temperature, &
    & energy_elastic, strain_elastic)

  use parameter_module, only: moxide_layer
  use solution_data_module, only: mout_in_f2l, no_f2l_event_out
  use oxide_data_module,   only: first_oxide_layer, thickness_fr

  integer, intent(IN)  :: unit, N, j
  real   , intent(IN)  :: dtime_f2l_now
  real, intent(IN), dimension(10, mout_in_f2l)   :: temperature
  real, intent(IN), dimension(10, mout_in_f2l, moxide_layer)   :: &
       & energy_elastic, strain_elastic

  integer :: i, k
  real, dimension(10)   :: strain23, energy23  ! for the 2+3 layers combined

  if (N <= 3)  then

     write(unit, 92) dtime_f2l_now * (j - 1), &
        & (temperature(i, j), energy_elastic(i, j, 1), &               ! u_tot
        &  (energy_elastic(i, j, k), &                                 ! u_sp(mag)
        &  strain_elastic(i, j, k), k = first_oxide_layer, N), & ! e_sp(mag)
        &  i = 1, no_f2l_event_out)

 92           format(1x, 1pe12.5, 70(1pe12.5, 1x))

  else if (N == 4)  then

     do i = 1, no_f2l_event_out

       strain23(i) = SUM(strain_elastic(i, j, 3:N) * thickness_fr(3:N)) / &
        & SUM(thickness_fr(3:N))
       energy23(i) = SUM(energy_elastic(i, j, 3:N) * thickness_fr(3:N)) / &
        & SUM(thickness_fr(3:N))

     end do

     write(unit, 92) dtime_f2l_now * (j - 1), &
        & (temperature(i, j), energy_elastic(i, j, 1), &               ! u_tot
        &  (energy_elastic(i, j, k), &                                 ! u_sp(mag)
        &  strain_elastic(i, j, k), k = first_oxide_layer, N), & ! e_sp(mag)
        &  energy23(i), strain23(i), i = 1, no_f2l_event_out)

  endif

  return

  END SUBROUTINE WRITE_EVENT

END MODULE OXIDE_UTILITY_MODULE
