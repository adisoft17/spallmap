MODULE CREEP_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all CREEP data structures integer and scalar parameters.
  !
  !! Author: Adrian S. Sabau, sabaua@ornl.gov creep_data_module.F90
  !
  !=======================================================================
  use parameter_module, only: moxide_layer, moxide, mtime, mrate, mave, &
        & mgrid

  implicit none
  save
  ! Public Module
  public

  logical  :: if_creep_ox_scale, if_creep_metal

  ! strain_rate_creep = creep_const * sigma_eq**creep_exponent * 
  !      EXP(-creep_activation_energy / Temperature)
  ! Temp_creep_onset = threshold temperature at which creep 
  !      starts to be considered
  ! without _mat are used for each layer that is made up of different materials
  real, dimension(0:moxide)     :: creep_const, creep_exponent, &
       & creep_activation_energy, creep_activ_energy_stress, &
       & Temp_creep_onset
  ! use the rupture time data as in Takahashi and EPRI
  real, dimension(0:moxide)     :: rupture_time_c1, rupture_time_t1, &
       & rupture_time_t1s1, rupture_time_t1s2, rupture_time_exp

  ! creep_const_ln = natural log of creep constant to better evaluate 
  ! the creep rate
  real, dimension(0:moxide)     :: creep_const_ln, creep_const_mat_ln

  ! _mat is used at input for each material
  real, dimension(0:moxide)     :: creep_const_mat, creep_exponent_mat, &
       & creep_activation_energy_mat, creep_activ_energy_stress_mat, &
       & Temp_creep_onset_mat

  ! use the rupture time data as in Takahashi and EPRI
  real, dimension(0:moxide)     :: rupture_time_c1_mat, rupture_time_t1_mat, &
       & rupture_time_t1s1_mat, rupture_time_t1s2_mat, rupture_time_exp_mat
  ! if_rupture_time_creep = T, then tr is used to get the creep rate
  logical, dimension(0:moxide)  :: if_rupture_time_creep_mat, &
       & if_rupture_time_creep

  ! sigma_eq = equivalent stress given by von Mises expression
  ! press_eq = average of stresses in the three directions
  real, dimension(moxide, mgrid) :: stress_eq, stress_eq_old

  ! store the stress level after convergence for the next cycle
  ! stress_eq_conv(1, :, :) = high load
  ! stress_eq_conv(2, :, :) = low load
  ! stress_eq_conv(k, :, :) = for the k-2 interval during transitions
  real, dimension(10, moxide, mgrid) :: stress_eq_conv

  ! operation_2_initial_id maps the id_temp_op with the initial point where
  ! the stress_eq_conv need to be stored
  integer, dimension(8), parameter  :: operation_2_initial_id = &
         &  (/1, 2, 2, 1, 1, 2, 2, 1/)

  ! no_intervals_now = total number of intervals for the current segment
  integer, dimension(8) :: no_intervals_now

  ! operation_id_current = stores id_temp_op(i) current
  ! interval_op_now = current interval now considered for the current segment
  integer  :: operation_id_current, interval_op_now

  ! creep_strain_eq = equivalent uniaxial creep strain (valid for von Mises equivalent stress expression)
  ! creep_strain_incr_eq - must be stored such that the creep_strain_eq
  !  can be computed outside the creep_strain routine, when the convergence was established
  real, dimension(moxide, mgrid) :: creep_strain_eq, creep_strain_incr_eq
  

  real, dimension(moxide, mave) :: stress_eq_ave, creep_strain_eq_ave

  ! stress_eq_epsilon = a small number to test the 
  ! convergence of SUM |stress_eq - stress_eq_old| < stress_eq_epsilon * stress_eq
  real  :: stress_eq_epsilon

  ! eps_oxide_grid maximum difference in the grid positions 
  ! between the current and previous
  real  :: eps_oxide_grid

  ! new_oxide_creep_option = 'old' - do not initialize
  ! new_oxide_creep_option = 'zero' - start from zero
  ! new_oxide_creep_option = 'half_old' - start from half of the previous creep strain
  character(LEN = 80)   :: new_oxide_creep_option

  ! il_now, ip_now id for layer and index in the previous layer for the new grid 
  ! with respect to the old grid
  integer, dimension(moxide, mave) :: il_now, ip_now

  ! dcreep = current creep increment = dt * creep_rate 

  type CREEP_TYPE1

    real  ::  rad
    real  ::  hoop
    real  ::  ax

  end type CREEP_TYPE1

  type(CREEP_TYPE1), dimension(moxide, mgrid) :: creep_strain_old, &
       & creep_strain_new
  type(CREEP_TYPE1), dimension(moxide, mave) :: creep_strain_ave

  ! creep_strain_mean - used where average creep strain is considered
  ! creep_strain_mean_old - old time level
  type(CREEP_TYPE1), dimension(moxide) :: creep_strain_mean, &
       & creep_strain_mean_old

  ! creep_sum = e_CR,r + e_CR,hoop + 2*poisson*e_CR,z
  real, dimension(moxide, mgrid) :: creep_sum

  ! creep_integral_sum = integral of r * creep_sum
  real, dimension(moxide, mgrid) :: creep_integral_sum

  ! creep_integral_z = Integral (r*e_z^CR)
  real, dimension(moxide) :: creep_integral_z

  ! creep_dif = e_CR,r - e_CR,hoop
  real, dimension(moxide, mgrid) :: creep_dif

  ! creep_integral_dif = integral of creep_dif / r
  real, dimension(moxide, mgrid) :: creep_integral_dif

  ! i_cr_start, i_cr_end == start and end indices for creep layers
  integer :: i_cr_start, i_cr_end

  ! real, dimension(moxide, mave) :: 

END MODULE CREEP_DATA_MODULE
