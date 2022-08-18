MODULE MESH_MODULE

    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for setting dimensions and meshes. 
    !
    ! sabaua@ornl.gov mesh_module.F90
    !
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================


  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: GET_RADIUS, GET_RADIUS_OH, GET_MESH, SOLVER_STRESS_VARIABLES

 CONTAINS

  SUBROUTINE GET_RADIUS_OH(N, thickness_oxide_layer, if_ave_scale, iter)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the temperature profile
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module,   only: moxide, mave
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & oxide_location
  use solution_data_module, only: rad_temp, npr_temp, &
      & npr_st, rad_st, rad_int
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, no_st_rad_points_tube, &
      & no_st_rad_points_oxide 
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: imetal, ispinel, imagnetite, &
      & iout, imet_spin, ispin_mag, i_in, if_growth_induce_stress, &
      & h_ox_sp, h_ox_mag, rc_og

  integer, intent(IN)  :: N, iter
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  logical, intent(IN)  :: if_ave_scale
  character(LEN = 80)  :: char1  

  ! Local Variables
  integer  :: i, j, loc1, k
  integer, save  :: ip = 0
  real     :: tube_outer_radius_now, tube_thickness_now
  real, dimension(moxide)  :: thickness_oxide_layer_now

  ! N = total number of layers; minimum 2 (substrate + average scale)

  if (TRIM(oxide_location) == 'outer_metal_surface')  then
    write(6, *) 'outer oxide not fixed for oh growth'
    STOP
  endif

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  ! initialize
  ! if (if_first_temp_rad)   then cannot do it since oxide thickness
  ! varies

  if (N + 1 > moxide)  then
    write(tty_lun, *) 'ERROR: radius cannot be computed '
    stop
  endif

  ! tube_outer_radius - is in fact tube radius at RT
  ! tube_thickness - is in fact tube thickness at RT
  ! temperature field is non-uniform and the tube expansions is also non-uniform
  ! more iterations are require to get the expanded radius since temperature 
  ! profile is needed  

  if (iter == 1)  then

    tube_outer_radius_now = tube_outer_radius
    tube_thickness_now = tube_thickness
    thickness_oxide_layer_now = thickness_oxide_layer

  else

    ! iter > 1 for the oxide growth problems
    if (if_growth_induce_stress)  then

      ! r_i, r_o, var(ir_r_o_ma) - var(ir_r_i_s) are known
      tube_outer_radius_now = rc_og(iout)   ! r_o
      tube_thickness_now = rc_og(iout) - rc_og(imet_spin)
      thickness_oxide_layer_now(ispinel) = h_ox_sp
      thickness_oxide_layer_now(imagnetite) = h_ox_mag

    else

      tube_outer_radius_now = tube_outer_radius
      tube_thickness_now = tube_thickness
      thickness_oxide_layer_now = thickness_oxide_layer

    endif

  endif

  rad_int(1) = tube_outer_radius_now
  rad_int(2) = rad_int(1) - tube_thickness_now

  ! deal with the scale that got average properties
  do i = 3, N+1

    if (if_ave_scale)  then

      if (if_growth_induce_stress)  then
        write(out_lun, *) 'ERROR: average scale and oxide-growth-induced stress not yet'
        stop
      endif

      ! in this case N = 2 (tube and one average oxide)
      rad_int(i) = rad_int(i-1) - thickness_oxide_layer(no_oxide_ave)  ! only for average

      if (ip <= 6)  then
        write(aux_lun, 7) i, thickness_oxide_layer(no_oxide_ave), &
           & rad_int(i-1), rad_int(i), &
           & thickness_oxide_layer(no_oxide_ave)
      endif

    else

      ! this will work in general; thickness_oxide_layer(1) is zero by default (metal)
      rad_int(i) = rad_int(i-1) - thickness_oxide_layer_now(i-1)  ! oxide layer i-2

      if (ip <= 6)  then

        write(aux_lun, 7) i, thickness_oxide_layer_now(i-1), rad_int(i-1), rad_int(i), &
           & thickness_oxide_layer(no_oxide_ave)
      ! write(aux_lun, 7) i,  (thickness_oxide_layer(j), j=1, moxide)
 7      format('mat_interf ', i2, 1x, 20(1pe13.6, 1x))

      endif

    endif

  end do

  return

  END SUBROUTINE GET_RADIUS_OH

  SUBROUTINE GET_RADIUS(N, thickness_oxide_layer, time_now, &
       & if_ave_scale, iter)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the temperature profile
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module,   only: moxide, mave
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & oxide_location
  use solution_data_module, only: rad_temp, npr_temp, &
      & npr_st, rad_st, rad_int, rad_mean
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value, no_metal_layers
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, no_st_rad_points_tube, &
      & no_st_rad_points_oxide, update_sola_type
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & rr_metal_in_initial, u_int_jump_now
  use creep_data_module, only: i_cr_start, i_cr_end

  integer, intent(IN)  :: N, iter
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  real, intent(IN)      :: time_now
  logical, intent(IN)  :: if_ave_scale
  character(LEN = 80)  :: char1  

  ! Local Variables
  integer  :: i, j, loc1, k
  integer, save  :: ip = 0
  real     :: tube_outer_radius_now, tube_thickness_now, &
         & first_spacing, ratio_metal_layer
  real, dimension(moxide)  :: thickness_oxide_layer_now

  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  ! initialize
  ! if (if_first_temp_rad)   then cannot do it since oxide thickness
  ! varies

  if (N + 1 > moxide)  then
    write(tty_lun, *) 'ERROR: radius cannot be computed '
    stop
  endif

  ! initialize the initial inner radius of the tube; reference radius
  rr_metal_in_initial = tube_outer_radius - tube_thickness ! * &
      ! & (1.0 + tau_th_st(imetal, 1)*(1.0 + poisson_ratio(imetal)))

  ! get position of metal-oxide interface
  call GET_METAL_DIM(N, tube_thickness_now, &
           & thickness_oxide_layer, time_now, iter)

  ! tube_outer_radius - is in fact tube radius at RT
  ! tube_thickness - is in fact tube thickness at RT
  ! temperature field is non-uniform and the tube expansions is also non-uniform
  ! more iterations are require to get the expanded radius since temperature 
  ! profile is needed  

  ! if (iter == 1)  then

    ! initialize
    ! tube_outer_radius_now = tube_outer_radius
    ! tube_thickness_now = tube_thickness
    ! thickness_oxide_layer_now = thickness_oxide_layer

  ! else

    ! iter > 1 for the oxide growth problems
    if (if_growth_induce_stress)  then

      ! this will be more complex if outer radius of the tube varies
      ! or if the reference dimensions due to oxide growth are considered

      ! r_i, r_o, var(ir_r_o_ma) - var(ir_r_i_s) are known
      tube_outer_radius_now = tube_outer_radius   ! r_o
      thickness_oxide_layer_now = thickness_oxide_layer

    else

      tube_outer_radius_now = tube_outer_radius
      thickness_oxide_layer_now = thickness_oxide_layer

    endif

  ! endif

    if (TRIM(oxide_location) == 'inner_metal_surface')  then
    
      rad_int(1) = tube_outer_radius_now

    else if (TRIM(oxide_location) == 'outer_metal_surface')  then

      rad_int(1) = - tube_outer_radius_now

    else

      write(6, *) 'oxide_location must be inner_metal_surface or ', &
        & ' outer_metal_surface'
      write(6, *) 'oxide_location =', TRIM(oxide_location)
      STOP

    endif
  
    ! valid only for no_metal_layers = 1
    ! rad_int(2) = rad_int(1) - tube_thickness_now

    ratio_metal_layer = 0.8

    if (no_metal_layers == 1)  then

      do i = 1, no_metal_layers
        rad_int(i+1) = rad_int(i) - tube_thickness_now / no_metal_layers
      end do

    else if (no_metal_layers > 1)  then

      if (ABS(ratio_metal_layer-1.0) <= 1.0e-4)  then

        first_spacing = tube_thickness_now / no_metal_layers
        do i = 1, no_metal_layers
          rad_int(i+1) = rad_int(i) - tube_thickness_now / no_metal_layers
        end do

      else

        first_spacing = tube_thickness_now *(1.0 /ratio_metal_layer-1.0) / &
           & (1.0/ratio_metal_layer**no_metal_layers-1.0)
        do i = 1, no_metal_layers
          rad_int(i+1) = rad_int(i) -first_spacing/ratio_metal_layer**(i-1)
        end do

      endif

    endif

    do i = 1, no_metal_layers
 
      rad_mean(i) = 0.5 * (rad_int(i+1) + rad_int(i))

    end do

    ! deal with the scale that got average properties
    do i = no_metal_layers + 2, N+1

      if (if_ave_scale)  then

        if (if_growth_induce_stress)  then
          write(out_lun, *) 'ERROR: average scale and', &
            & ' oxide-growth-induced stress not yet'
          stop
        endif

        ! in this case N = 2 (tube and one average oxide)
        rad_int(i) = rad_int(i-1) - thickness_oxide_layer(no_oxide_ave)  ! only for average

        if (ip <= 6)  then
          write(aux_lun, 7) i, thickness_oxide_layer(no_oxide_ave), &
             & rad_int(i-1), rad_int(i), &
             & thickness_oxide_layer(no_oxide_ave)
        endif

      else

        ! this will work in general; thickness_oxide_layer(1) is zero by default (metal)
        rad_int(i) = rad_int(i-1) - thickness_oxide_layer_now(i-1)  ! oxide layer i-2

        if (ip <= 6)  then

          write(aux_lun, 7) i, thickness_oxide_layer_now(i-1), &
             & rad_int(i-1), rad_int(i), &
             & thickness_oxide_layer(no_oxide_ave)
        ! write(aux_lun, 7) i,  (thickness_oxide_layer(j), j=1, moxide)
 7        format('mat_interf ', i2, 1x, 20(1pe13.6, 1x))

        endif

      endif

    end do

    do i = no_metal_layers + 1, N
 
      rad_mean(i) = 0.5 * (rad_int(i+1) + rad_int(i))

    end do

    if (iter > 0 .and. time_now > 1.0e-3)  then

    write(aux_lun, 8) time_now, rad_int(1) - rad_int(N+1), &
         & (rad_int(i), i=1, N+1)
 8  format('rad_interf ', 20(1pe13.6, 1x))

    endif

  return

  END SUBROUTINE GET_RADIUS

  SUBROUTINE SOLVER_STRESS_VARIABLES(N)

    !=======================================================================
    ! Purpose(s):
    !
    !   set up variables id
    !  
    !=======================================================================

  use solution_data_module, only: ib, ic, i_cr_r, i_cr_h, i_cr_z, i_s_r, &
      & i_s_h, i_s_z, i_s_vm, i_cr, ieps0, no_stress_var, mean_creep_layer
  use solver_data_module, only       : fval, var, update_sola_type, &
        & delta_var, rhs, coeff, update_coeff, update_const, var_conv
  use creep_data_module, only: i_cr_start, i_cr_end, creep_strain_mean
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

  ! Argument List
  integer,    intent(IN)  :: N

  ! Local Variables
  integer  :: i, ivar, Status

  ieps0 = 1

  ivar = 1

  if (update_sola_type == 0)  then
    do i = 1, N
      ivar = ivar + 1
      ib(i) = ivar
      ic(i) = ib(i) + 1
      ivar = ivar + 1
    end do

  else if (update_sola_type > 0)  then

  do i = 1, N
    ivar = ivar + 1
    ib(i) = ivar
    ic(i) = ib(i) + 1

    ! add the creep terms only when needed
    if (i >= i_cr_start .and. i <= i_cr_end)  then

      i_cr_r(i) = ic(i) + 1        ! radial creep increment
      i_cr_h(i) = i_cr_r(i) + 1    ! hoop   creep increment
      i_cr_z(i) = i_cr_h(i) + 1    ! axial  creep increment

      i_s_r(i) = i_cr_z(i) + 1     ! radial mean stress
      i_s_h(i) = i_s_r(i) + 1      ! hoop   mean stress
      i_s_z(i) = i_s_h(i) + 1      ! axial  mean stress

      i_s_vm(i) = i_s_z(i) + 1     ! von Mises mean stress

      i_cr(i) = i_s_vm(i) + 1      ! creep strain increment

      ivar = ivar + 9

      ! Mask_var_conv(i_s_vm(i)) = .true.

    else

      ivar = ivar + 1

    endif

  end do  

  endif

  no_stress_var = ivar
  
  write(out_lun, *) 'no_variables ', no_stress_var

  ALLOCATE(var(1:no_stress_var), STAT = Status)
  ALLOCATE(fval(1:no_stress_var), STAT = Status)
  ALLOCATE(rhs(1:no_stress_var), STAT = Status)
  ALLOCATE(delta_var(1:no_stress_var), STAT = Status)
  ALLOCATE(coeff(1:no_stress_var, 1:no_stress_var), STAT = Status)

  ALLOCATE(update_coeff(1:no_stress_var, 1:5), STAT = Status)
  ALLOCATE(update_const(1:no_stress_var), STAT = Status)
 
  ALLOCATE(var_conv(1:10, 1:no_stress_var), STAT = Status)

    if (update_sola_type == 0)  then

      RETURN

    else if (update_sola_type == 1 .or. &
      & update_sola_type == 2 .or. &
      & update_sola_type == 3)  then

      do i = i_cr_start, i_cr_end

        ! zero creep strains
        var(i_cr_r(i)) = 0.0
        var(i_cr_h(i)) = 0.0
        var(i_cr_z(i)) = 0.0

        creep_strain_mean(i)%rad = 0.0
        creep_strain_mean(i)%hoop = 0.0
        creep_strain_mean(i)%ax = 0.0

        ! zero mean stresses
        var(i_s_r(i)) = 0.0
        var(i_s_h(i)) = 0.0
        var(i_s_z(i)) = 0.0

        var(i_s_vm(i)) = 0.0

        var(i_cr(i)) = 0.0

      end do

      do i = 1, 10
        var_conv(i, :) = var(:)
      end do

    ! else if (update_sola_type == 2)  then

    ! else if (update_sola_type == 3)  then

    endif

  return

  END SUBROUTINE SOLVER_STRESS_VARIABLES

  SUBROUTINE GET_METAL_DIM(N, tube_thickness_now, thickness_oxide_layer, &
      & time_now, iter)

    !=======================================================================
    ! Purpose(s):
    !
    !   get metal dimensions
    !  
    !=======================================================================

  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_stress_data_module,  only: if_growth_induce_stress, &
      & total_oxide_thickness_ref, rr_metal_in_initial
  use oxide_data_module,   only: pilling_bedworth_ratio, thickness_fr, &
      & first_oxide_layer, no_metal_layers
  use boiler_data_module, only: tube_outer_radius, tube_thickness
  use parameter_module,   only: moxide
  use solver_data_module, only  : metal_recession_type

  ! Argument List
  integer,    intent(IN)  :: N, iter
  real,    intent(OUT)  :: tube_thickness_now
  real, intent(IN)      :: time_now
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer

  ! Local Variables
  integer  :: i, k, j, mn, kp
  real, dimension(N+1)      :: psi
  real     :: b_term, c_term, r_metal_oxide, u_jump_fr, &
           & total_oxide_thickness_ref_si, r_approx

  if (TRIM(metal_recession_type) == 'constant')  then
    ! unchanged geometry
    tube_thickness_now = tube_thickness
    RETURN
  endif

  ! total oxide thickness in SI units
  total_oxide_thickness_ref_si = &
       & SUM(thickness_oxide_layer(first_oxide_layer:N))

  if (no_metal_layers > 1)  then
    ! consider that metal is one layer since no details is required here
    mn = N - no_metal_layers + 1
  else if (no_metal_layers == 1)  then
    mn = N
  endif

  ! get psi for each layer
  psi(1:mn+1) = 0.0
  psi(2) = thickness_fr(2) / 2.0
  do k = 3, mn   ! N layers including the metal
    kp = k -2 + first_oxide_layer
    psi(k) = psi(k-1) + 0.5 * (thickness_fr(kp-1) + thickness_fr(kp))
  end do

  ! get ri term and free term
  b_term = 0.0
  c_term = 0.0
  do k = 2, mn
    kp = k -2 + first_oxide_layer
    b_term = b_term + thickness_fr(kp) / pilling_bedworth_ratio ! (k)
    c_term = c_term + psi(k) * thickness_fr(kp) / pilling_bedworth_ratio ! (k)
  end do

  b_term = b_term * total_oxide_thickness_ref_si
  c_term = c_term * 2.0 * total_oxide_thickness_ref_si**2 - &
         & rr_metal_in_initial**2

  ! solve the quadratic equation for the inner metal radius (metal_oxide interf)
  ! r^2 - 2 * b * r + c = 0
  r_metal_oxide = b_term + sqrt(b_term**2 - c_term)  ! this should be positive
  ! r_sol(2) = b_term - sqrt(b_term**2 - c_term), this is negative
  
  ! small oxide approximation
  r_approx = 0.5 * (rr_metal_in_initial**2 - c_term) / &
           & (rr_metal_in_initial - b_term)

  if (iter < 5)  then

    write(out_lun, 2) iter, b_term, c_term, &
        & r_approx, b_term - sqrt(b_term**2 - c_term), &
        & b_term + sqrt(b_term**2 - c_term)
2 format('check_met_recess ', i4, 1x, 20(1pe13.6, 1x))

  endif

  tube_thickness_now = tube_outer_radius - r_metal_oxide

  if (iter > 0)  then
    ! displacement of metal_oxide interface
    write(out_lun, 1) iter, time_now, r_metal_oxide - &
      & (tube_outer_radius - tube_thickness), &
      & tube_thickness_now, tube_thickness, &
      & r_metal_oxide, r_approx, tube_outer_radius, &
      & total_oxide_thickness_ref_si, rr_metal_in_initial
1   format('met_ox_displ ', i4, 1x, 20(1pe13.6, 1x))
  endif

  return

  END SUBROUTINE GET_METAL_DIM

  SUBROUTINE GET_MESH(N, thickness_oxide_layer, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the temperature profile
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module,   only: moxide, mave
  use boiler_data_module, only: tube_thickness, tube_outer_radius
  use solution_data_module, only: rad_temp, npr_temp, &
      & npr_st, rad_st, rad_st_old, rad_int, n_mean_even_vm, n_mean_odd_vm
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, no_st_rad_points_tube, &
      & no_st_rad_points_oxide 
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

  integer, intent(IN)  :: N
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  ! real, intent(IN)      :: time_now
  logical, intent(IN)  :: if_ave_scale
  character(LEN = 80)  :: char1  

  ! Local Variables
  integer  :: i, j, loc1, k
  integer, save  :: ip = 0

  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  do i = 1, N

    ! get number of points in the radial direction
    if (i == 1)  then
      npr_temp(i) = no_temp_rad_points_tube
      npr_st(i) = no_st_rad_points_tube
    else 
      npr_temp(i) = no_temp_rad_points_oxide
      npr_st(i) = no_st_rad_points_oxide
    endif
    
    n_mean_even_vm(i) = INT(npr_st(i)/2)
    n_mean_odd_vm(i) = INT((npr_st(i)+1)/2)

    if (npr_st(i) == 2 * n_mean_even_vm(i))  then
      ! npr_st is an even number
      ! n_mean1 = INT(npr_st(i)/2)
      ! integral_th_mean(i) = (integral_th(i, n_mean) + &
      !        & integral_th(i, n_mean+1)) / 2.0
    else
      ! n_mean1 = INT((npr_st(i)+1)/2)
      if (npr_st(i) + 1 == 2 * n_mean_odd_vm(i))  then
        ! npr_st is odd
        ! integral_th_mean(i) = integral_th(i, n_mean1)
      endif
    endif

    ! get the radius first
    do j = 1, npr_temp(i)
      rad_temp(i, j) = rad_int(i) - (rad_int(i) - rad_int(i+1)) * &
         & (j - 1) / (npr_temp(i) - 1)
    end do

    ! store the old radii for the creep module
    rad_st_old(i, 1:moxide) = rad_st(i, 1:moxide)

    do j = 1, npr_st(i)
      rad_st(i, j) = rad_int(i) - (rad_int(i) - rad_int(i+1)) * &
         & (j - 1) / (npr_st(i) - 1)
    end do

  end do

  return

  END SUBROUTINE GET_MESH

END MODULE MESH_MODULE
