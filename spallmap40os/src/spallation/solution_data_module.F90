MODULE SOLUTION_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all solution data structures integer and scalar parameters.
  !  !
  !  Author: Adrian S. Sabau, sabaua@ornl.gov
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================
  use parameter_module, only: mtime, mstep, moxide, mave, mramp, &
      & mramp_out, moxide_layer, mgrid

  implicit none
  save
  ! Public Module
  public

  ! number of cycles at which the transition from full load to low load
  ! is printed out
  ! no_point_full2low - contains the number of points in the TEMP_BOILER_data_no
  integer        :: no_full2low_output, no_f2l_event_out
  integer, dimension(mramp_out)    :: no_pulse_full2low_output, &
      & no_point_full2low
  integer, dimension(10)    :: no_event_full2low_output, &
      & no_event_point_full2low

  ! maximum number of outputs during an outage event
  integer, parameter     :: mout_in_f2l = 300

  ! time and oxide thickness at which the ramp down is followed
  real, dimension(mramp_out, mramp)    :: time_full2low, oxide_full2low
  real, dimension(mramp_out, mramp)    :: temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low

   ! max_iterat_solver = the maximum number of iterations
   integer           :: max_iterat_solver

   ! rad = radius of each point
   ! temp = temperature at each point
   ! strain_rad - radial strain
   ! strain_hoop - hoop strain
   ! strain_hoop_gen - effective hoop strain that generates stress
   !                 strain_hoop_gen = total_strain - thermal_free_stress_strain
   ! displ_rad - radial displacement
   ! strain_th = integral of alfa(T)/alfa(To) dT from T_i_hot to T_i_cold
   ! integral_th = integral of strain_th * r dr from r_i to current radius
   !
   ! strain_hoop_el - elastic strain (not equal to strain_hoop_gen)
   ! strain_rad_el - elastic strain
   ! strain_axial_el - elastic strain

   ! rad_st - new radius position for the stress calculations
   ! rad_st_old - old radius position needed for the creep 
   ! real, dimension(0:moxide, 0:mgrid) :: rad_st_old
   real, dimension(moxide, mgrid) :: rad_temp, temp, &
     & rad_st, rad_st_old, Temp_high_rad, Temp_Low_Rad, &
     & strain_rad, strain_hoop, displ_rad, &
     & stress_axial, stress_hoop, stress_rad, &
     & strain_th, integral_th, strain_hoop_gen
  real, dimension(moxide, mgrid) :: strain_hoop_el, strain_rad_el, &
     & strain_axial_el

  ! tau is the equivalent of beta * DT (or thermal strain)
  ! tau_th_st is the tau_th at the stress points
  real, dimension(moxide, mgrid)  :: tau_th, tau_th_st

  ! npr - number of points in the radial direction
  ! npr_temp points for the temp
  ! npr_st points for the stress
  integer, dimension(moxide)  :: npr_temp, npr_st

  ! mave: maximum number of averaging quantities
  ! mave = 1 - max
  ! mave = 2 - min
  ! mave = 3 - sum average
   real, dimension(moxide, mave) :: rad_ave, temp_ave, &
     & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
     & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
     & strain_th_ave, integral_th_ave, &
     & strain_th_plate_ave, stress_hoop_th_ave, &   ! for the plane geometry
     & strain_hoop_gen_ave  ! effective hoop strain that generates stress
   ! strain_hoop_gen_ave = total_strain - thermal_free_stress_strain - growth_strain
   !    - creep_strain (from the plain stress Hook's law)

   ! strain_hoop_el_ave = total strain - thermal_free_stress_strain - growth_strain
   !    - creep_strain (from the Hook's law)
   real, dimension(moxide, mave) :: strain_hoop_el_ave

   ! define armit strains as sigma_hoop/Youngs in order to enable dirrect comparison; use the strains that generate strain - a better quantity to use
   ! real, dimension(moxide, mave) :: strain_hoop_el_av_armit

   ! T_i = Tin + A_i + B_i * ln(r) ; ac_temp and bc_temp are temperature const.
   ! ac_temp(1, i) is the a constant for layer i, after the end of heating
   ! i.e., A - constant at high temperature; same bc
   real, dimension(2, moxide)   :: ac_temp, bc_temp

   ! stress-strain constants
   real, dimension(moxide)   :: bc_st, cc_st

   ! delta_ac_temp = ac_temp(1, :) - ac_temp(2, :)
   real, dimension(moxide)   :: delta_ac_temp, delta_bc_temp

   real, dimension(moxide)   :: rad_int  ! radius at the interface

   real, dimension(moxide)   :: rad_mean  ! average radius per each layer

   ! quantities defined at the mean temperature location
   real, dimension(moxide)   :: Temp_mean, hoop_stress_mean, ax_stress_mean, &
     & rad_stress_mean, vm_stress_mean, creep_eq_incr_mean, &
     & integral_th_mean, tau_th_mean, vm_stress_mean_old

   ! real, dimension(10, moxide)   :: vm_stress_mean_conv  ! similar to stress_eq_conv
   ! real, dimension(10, moxide)   :: bc_st_conv, cc_st_conv, 
   ! real, dimension(10)   :: eps0_conv

   ! use vm_stress_mean_old and vm_stress_mean to initialize the vm stresses

   ! no_stress_var - number of variables in the stress solver
   integer   :: no_stress_var
   integer  :: ieps0

   !  i_cr_r(i) = ic(i) + 1        ! radial creep increment
   !  i_cr_h(i) = i_cr_r(i) + 1    ! hoop   creep increment
   !  i_cr_z(i) = i_cr_h(i) + 1    ! axial  creep increment

   !  i_s_r(i) = i_cr_z(i) + 1     ! radial mean stress
   !  i_s_h(i) = i_s_r(i) + 1      ! hoop   mean stress
   !  i_s_z(i) = i_s_h(i) + 1      ! axial  mean stress
   !  i_s_vm(i) = i_s_z(i) + 1     ! von Mises mean stress

   !  i_cr(i) = i_s_vm(i) + 1      ! creep strain increment

   integer, dimension(moxide) :: ib, ic, i_cr_r, i_cr_h, i_cr_z, i_s_r, i_s_h, &
      & i_s_z, i_s_vm, i_cr, n_mean_odd_vm, n_mean_even_vm

   ! n_mean_odd_vm - indicate mid point the stress mesh to evaluate VM
   ! n_mean_odd_vm - for odd npr_st and n_mean_even_vm for even npr_st

   ! identify equation type to get the convergence criterial
   integer, dimension(moxide) :: eq_vm, eq_cr

   ! mean_creep_layer = .true. - the crep rate and creep strain is considered constant in this layer
   logical, dimension(moxide)    :: mean_creep_layer

   ! delta_Temp_in = Temp_in_hot - Temp_in_cold
   real   :: delta_Temp_in

   ! eps0; the strain for generalized plain strain formulations
   real   :: eps0

   ! no_oxide_layers = 2 if if_average_oxide = .true.
   ! no_oxide_layers = no_oxide if if_average_oxide = .true.
   integer :: no_oxide_layers

   ! output data for the entire duration of oxide growth
   ! 1 - full load, 2 - partial load
   ! strain_elastic_out = strain_elastic_sp, _mag, _tot
   real, dimension(2, 3*moxide_layer)   :: strain_elastic_out, &
    & energy_elastic_out, temperature_out, strain_elastic_out_plate, &
    & energy_elastic_out_plate

   ! dist_energy_elastic_out - distortion energy
   real, dimension(2, 3*moxide_layer)   :: dist_energy_elastic_out, &
    & dist_energy_elastic_out_plate

   ! strain_elastic_out_gen(min_ave_max, high_low, position)
   real, dimension(mave, 2, 3*moxide_layer)   :: strain_elastic_out_gen
   
   ! output data for the outage events
   ! three variables for outage, and up to 10 outages
   ! strain_elastic_event_out = strain_elastic_event_sp, _mag, _tot
   ! energy_elastic_event_out = energy_elastic_event_sp, _mag, _tot
   real, dimension(10, mout_in_f2l)   :: temperature_event_out
   real, dimension(10, mout_in_f2l, moxide_layer)   :: &
       & strain_elastic_event_out, strain_elastic_event_out_plate
   real, dimension(10, mout_in_f2l, moxide_layer)   :: &
       & energy_elastic_event_out, energy_dist_el_event_out, &
       & energy_elastic_event_out_plate, energy_dist_el_event_out_plate

   real, dimension(mave, mout_in_f2l, moxide_layer)   :: &
       & strain_elastic_event_out_g

   ! max_energy_elastic_event_out = MAX(energy_elastic_event_out(event, :, layer)
   real, dimension(10, moxide_layer)   :: max_energy_elastic_event_out, &
      & max_dist_energy_el_event_out
   ! _pl is for flat plate assumption; biaxial stress
   real, dimension(10, moxide_layer)   :: max_energy_elastic_event_out_pl, &
      & max_dist_energy_el_event_out_pl
   ! max_energy_elastic_event - to avoid the event index storage
   ! (1) will store the total elastic energy, then 2,..N will store for each oxide layer
   real, dimension(moxide_layer)       :: max_energy_elastic_event, &
      & max_dist_energy_el_event
   real, dimension(moxide_layer)       :: max_energy_elastic_event_pl, &
      & max_dist_energy_el_event_pl

   ! max_strain_elast_event_out = MAX(strain_elastic_event_out(event, :, layer)
   real, dimension(10, moxide_layer)   :: min_strain_elast_event_out, &
       & max_strain_elast_event_out
   real, dimension(10, moxide_layer)   :: min_strain_elast_event_out_pl, &
       & max_strain_elast_event_out_pl
   ! min_strain_elast_event - to avoid the event index storage
   ! (1) will store the total elastic energy, then 2,..N will store for each oxide layer
   real, dimension(moxide_layer)       :: min_strain_elast_event, &
       & max_strain_elast_event
   real, dimension(moxide_layer)       :: min_strain_elast_event_pl, &
       & max_strain_elast_event_pl

   ! max_strain_elast_event_gen(min_ave_max, layer)
   real, dimension(mave, moxide_layer)       :: min_strain_elast_event_gen, &
       & max_strain_elast_event_gen

   ! heat flux and temperature evolution during operation
   real, dimension(2, 30, 20)   :: heat_flux_temp_out   

   ! number of points for which heat flux information is written out
   integer  :: n_hf_out
   ! oxide thickness, time, and point in boiler operation at which hf is out
   real, dimension(30)  :: oxide_thickness_hf_out, time_hf_out
   ! maximum creep strain and von Mises stress at different thicknesses 
   ! at the same as time as those when heat flux info was written out
   ! oxide thickness, time, and point in boiler operation at which hf is out
   real, dimension(2, 30)  :: creep_hf_time, creep_hf_dox
   real, dimension(2, 30, moxide_layer)  :: creep_eq_hf_out, vm_hf_out
   ! creep_strain_ave_hoop_out - average creep strain in the hoop direction
   real, dimension(2, 30, moxide_layer)  :: creep_strain_ave_hoop_out

   ! average creep strain and von Mises stress at different thicknesses 
   ! creep_rate_out - creep rate (max, min, ave) in metal and each oxide
   real, dimension(2, 30, moxide_layer)  :: creep_eq_hf_out_ave, &
      & vm_hf_out_ave, creep_rate_out

   ! cycle at which hf info is written out
   integer, dimension(30)  :: id_full_load_hf

   ! to deal with the damage map strain4damage(outage, temperature, average_type)
   ! this strain would be only for magnetite as it exfoliates first
   !  or the the second and third layer (magnetite and mag+haematite) if N=4
   ! strain4damage(outage, 1, ) - at full load
   ! strain4damage(outage, 2, ) - at RT
   real, dimension(10, 10, 3)   :: strain_el_4damage, dox_4damage, &
      & strain_el_plate_4damage, strain_el_gen_4damage, temp_4damage

END MODULE SOLUTION_DATA_MODULE
