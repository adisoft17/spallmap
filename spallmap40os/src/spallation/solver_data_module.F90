MODULE SOLVER_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all combustion data structures integer and scalar parameters.
  !
  ! Author: Adrian S. Sabau, sabaua@ornl.gov
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================
  use parameter_module, only: mfuel, mcomb, mcp, mtref, mtint, ms, mh, &
                            & mspecie, string_dim, string_len, mdiluent, &
                            & moxidant, mgas, mstates, nmax
  use parameter_module, only: mtime, mstep, mgrid

  implicit none
  save
  ! Public Module
  public

   ! max_iterat_solver = the maximum number of iterations
   integer           :: max_iterat_solver

   ! no_trans_intervals = number of intervals used for f2l and l2f transitions
   ! no_trans_intervals = needs to be used when creep is considered
   integer           :: no_trans_intervals, no_intervals_f2l, &
      & no_intervals_l2f

   ! time step for growing the oxide scale (at given times)
   real, dimension(mstep) :: dt, time

   real, dimension(mstep)  :: temp_check_strain, Fe2O3_pct_check_strain

   ! actual number of data points for checking the strain
   ! no_thick_spall_map number of thickness points for the spallation map
   integer   :: ntemp_check_strain, nFe2O3_pct_check_strain, &
       & no_thick_spall_map

   ! interval type for thikness in the spallation map
   character(LEN = 80)  :: interval_type_thick_spall_map

   ! at input stores the minimum and maximum
   ! based on no_thick_spall_map and interval_type_thick_spall_map,
   ! the thickness for the spallation map will be stored here
   real, dimension(mstep) :: thickness_spallation_map

   ! level of EPRI map; level 1 include e_c vs thickness
   integer :: map_level

   ! the oxide thickness beyond which strains, stresses, damage is 
   ! computed
   real    :: small_oxide_thickness

   ! reference temperature for the tube when it is installed
   ! considered only when it is positive, i.e., 
   ! if_different_ref_temp = .true. if Temp_reference_tube > 0
   real  :: Temp_reference_tube, Temp_reference_oxide

   ! if_oxide_and_tube_temp = .true. when
   ! both Temp_reference_tube > 0 and Temp_reference_oxide > 0
   logical :: if_different_ref_temp, if_oxide_and_tube_temp

   !  if_average_oxide = .true. if all oxide layers are averaged
   !  if_average_oxide = .false. individual layers are considered
   logical :: if_average_oxide 

   ! actual time
   real              :: atime
   
   ! number at which the data is saved
   integer           :: ntime_save

   ! interval of time [hours] at which data is saved
   real              :: dtime_save
   
   ! actual times at which data is saved
   real, dimension(mtime)   :: time_save

   ! scale_var = scale unit for variable to be used with convergence criteria  
   real, dimension(nmax)   :: const, scale_var
   real(kind = 8), dimension(:), pointer   :: fval, var, delta_var, &
         & rhs, update_const
   real(kind = 8), dimension(:, :), pointer   :: coeff, update_coeff
   ! var_conv to update variables for low and full loads, similar to stress_eq_conv
   real(kind = 8), dimension(:, :), pointer   :: var_conv

   integer  :: update_sola_type

   integer, parameter  :: jeps0 = 1, jb_i = 2, jc_i = 3, j_s_r = 4, j_s_h = 5

   ! jd stores all the variables id's such that you can refer to variables 
   ! from previous iterations
   integer, dimension(nmax)   :: jd

   ! real, dimension(nmax, nmax)   :: jacf
   real, dimension(2, 2)   :: jacf

   ! just a parameter
   real :: dTemp_eps

   ! no_thick_interval - number of temp intervals for oxide scale growth
   ! no_temp_strain - number of temperature intervals for strain calculations
   integer           :: neq_solver, no_thick_interval, no_temp_strain
   ! no_out_f2l number of temp intervals for full-to-partial load transition
   integer           :: no_out_f2l
   integer           :: nrhs = 1

   ! numbers of problems to solve
   integer  :: no_problem

   ! true if generalized plain strain; otherwise just plain strain
   logical  :: generalized_plane_strain

   ! if_cylindrical = .true. then cylindrical plain stress will be used
   ! internal_htc = .true., then after buldging and wedging effects,
   ! the htc should be lower
   ! external_htc - .true. then htc will be used rather than temperatures
   logical  :: if_cylindrical, internal_htc, external_htc

   ! number of points in the radial direction where temperature
   !  will be computed
   integer  :: no_temp_rad_points_tube, no_temp_rad_points_oxide

   ! number of points in the radial direction where stress, strain, displ
   ! will be computed
   integer  :: no_st_rad_points_tube, no_st_rad_points_oxide

   ! metal_recession_type = constant; radius of metal-oxide interface
   ! does not change
   ! metal_recession_type = variable; radius of metal-oxide interface
   ! changes from a volumetric condition taken into account the 
   ! oxide thickness availability as a function of time
   character(LEN = 80) :: metal_recession_type
   
   ! no_first_pulse_creep_int = number of creep time intervals in the first 
   ! full cycle
   integer           :: no_first_pulse_creep_int

   ! ratio between adjacent creep time intervals for first pulse
   real :: ratio_time_first_pulse

   ! start_creep_cycles = how many full/partial load cycles are 
   ! divided in small increments for the onset of creep
   integer           :: start_creep_cycles

   ! store the number of stress iterations 
   integer  :: iter_stress_old

END MODULE SOLVER_DATA_MODULE
