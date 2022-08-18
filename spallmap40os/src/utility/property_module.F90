MODULE PROPERTY_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Estimate the property value and Enthalpy values based on 
    !     functionals given by CEA NASA data structure for Cp, H, and S. 
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: LINEAR_PROPERTY_SIMPLE, ALLOCATE_AVERAGE_PROPERTY, &
         & PROPERTY_SCALE_AVE, OPERAND_ARRAY, OPERAND_2LAYER_ARRAY, &
         & OPERAND_TWO_ARRAY

 CONTAINS

  SUBROUTINE LINEAR_PROPERTY_SIMPLE (Variable, &
     & data_no, data_size1, data_size2, Property_Value, &
     & Property_Variable, Value, Index)

    !======================================================================
    ! Purpose(s):
    !
    !   Evaluate the a property using linear interpolation
    !
    ! Input: property_type   - Thermophysical property to be updated.
    !        material_number - Material number of the material whose
    !                          propery is updated.
    ! Output: Value - Value of the thermophysical property (ncells array).
    !======================================================================
    ! use kind_module,          only: real_kind, int_kind, log_kind
    use parameter_module,     only: zero
    use output_module,      only: tty_lun, out_lun, aux_lun

    implicit none

    ! Argument List
    ! Property_data_no
    integer,         intent(IN) :: data_no, data_size1, data_size2 
    real,           intent(IN)  :: Variable
    real, dimension(data_size1:data_size2), intent(IN) :: Property_Value
    real, dimension(data_size1:data_size2), intent(IN) :: Property_Variable
    real, intent(OUT)           :: Value
    integer, intent(OUT), optional :: Index

    ! Local Variables
    integer  :: n, i
    real :: slope

    if (Variable <= Property_Variable(1))  then

      Value = Property_Value(1)
      if (PRESENT(INDEX))  Index = -1

    else if (Variable >= Property_Variable(data_no+1))  then

      Value = Property_Value(data_no+1)
      if (PRESENT(INDEX))  Index = - (data_no + 1)

    else

    PLIN_LOOP: do i=1, data_no

      if (Variable >= Property_Variable(i) .and. &
        &  Variable < Property_Variable(i+1))  then

        slope = (Property_Value(i+1) - Property_Value(i)) / &
              & (Property_Variable(i+1) - Property_Variable(i))

        Value = Property_Value(i)+ &
           & (Variable - Property_Variable(i)) * Slope

        if (PRESENT(INDEX))  Index = i

        exit PLIN_LOOP

      endif

    enddo PLIN_LOOP

    endif

    return

    END SUBROUTINE LINEAR_PROPERTY_SIMPLE

  SUBROUTINE ALLOCATE_AVERAGE_PROPERTY(id1, id2, id3, weight_ave_fr, &
        & average_type, nprop, property, variable)

    !=======================================================================
    ! Purpose(s):
    !
    !   Calculate average property for the entire scale, layers id1:id2
    !   based on a certain weighting criteria and average type
    !   allocate new data in id3 index
    !   property_i = property(variable_i), i= 1, nprop
    !  
    !=======================================================================

    use parameter_module,   only: moxide, mcp, zero
    use input_utilities_module,   only: SEARCH_NML
    use oxide_data_module,   only: no_oxide, no_layer, &
      & nth_exp, nyoungs, ncond, ncp_prop, id_mat, &
      & poisson_ratio, no_oxide_ave

    use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

    ! Argument List
    integer, intent(IN)     :: id1, id2, id3
    real, dimension(0:moxide), intent(IN) :: weight_ave_fr
    character(LEN = 80), intent(IN)     :: average_type
    integer, dimension(0:moxide), intent(INOUT) :: nprop
    real, dimension(0:moxide, 0:moxide), intent(INOUT) :: &
      & property, variable

    ! local variables
    integer, dimension(1)           :: id_max
    integer, dimension(1:id2-id1+1) :: n_intervals
    real :: prop1, prop2, var1
    integer :: i, j

    ! get the maximum number of the property
    n_intervals(1:id2-id1+1) = nprop(id1:id2)
    ! nprop(id3) = MAXVAL(nprop(id1:id2))
    ! id_max = MAXLOC(nprop(id1:id2))
    nprop(id3) = MAXVAL(n_intervals)
    id_max = MAXLOC(n_intervals) + id1 - 1
    write(aux_lun, *) 'material id with the largest #of intervals ', id_max(1)

    ! store the variable
    variable(id3, :) = variable(id_max(1), :)

    ! initialize
    if (TRIM(average_type) == 'add')  then

      property(id3, :) = 0.0

    else if (TRIM(average_type) == 'reciprocal_add')  then

      property(id3, :) = 0.0

    else if (TRIM(average_type) == 'sqrt')  then

      property(id3, :) = 1.0

    else 
 
      write(6, *) 'ERROR: average_type not ok'
      stop

    endif

    ! nprop holds the number of intevals not the number of data points
    ! number of data points = intervals + 1
    do j = 1, nprop(id3) + 1

      var1 =  variable(id3, j)  
      prop1 = property(id3, j)

      do i = id1, id2

        ! find the value of the property for layer i

        if (i == id_max(1))  then

          prop2 = property(i, j)

        else
          
          ! nprop(i) is the number of intervals
          call LINEAR_PROPERTY_SIMPLE (Var1, nprop(i), 0, &
              & moxide, Property(i, :), Variable(i, :), prop2)

        endif

        ! add the contribution from each layer

        if (TRIM(average_type) == 'add')  then

          prop1 = prop1 + weight_ave_fr(i) * prop2

        else if (TRIM(average_type) == 'reciprocal_add')  then

          prop1 = prop1 + weight_ave_fr(i) / prop2

        else if (TRIM(average_type) == 'geom')  then

          prop1 = prop1 * prop2**weight_ave_fr(i)

        endif

      end do

      ! finalyze
      if (TRIM(average_type) == 'add')  then

        property(id3, j) = prop1 / SUM(weight_ave_fr(id1:id2))

      else if (TRIM(average_type) == 'reciprocal_add')  then

        property(id3, j) = SUM(weight_ave_fr(id1:id2)) / prop1

      else if (TRIM(average_type) == 'geom')  then

        property(id3, j) = prop1**(1.0 / SUM(weight_ave_fr(id1:id2)))

      endif

    end do

  END SUBROUTINE ALLOCATE_AVERAGE_PROPERTY

  SUBROUTINE PROPERTY_SCALE_AVE ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Calculate average property for the entire scale
    !   based on a certain fraction distribution 
    !  average properties are hard to obtain if the temperature array 
    !  is not the same for all the layers 
    !=======================================================================
 
    use parameter_module,   only: moxide, mcp, zero
    use input_utilities_module,   only: SEARCH_NML
    use oxide_data_module,   only: material_name, cp_value, cp_temp, &
      & cond_value, cond_temp, rho, rho_value, rho_temp, &
      & Youngs_modul, Youngs_temp, thickness_fr, isubstrate, &
      & th_exp_coeff, th_exp_temp, no_oxide, no_layer, &
      & nth_exp, nyoungs, ncond, ncp_prop, id_mat, &
      & poisson_ratio, no_oxide_ave, &
      & surf_fracture_energy, fracture_toughness
    use oxide_data_module,   only: id_fe2o3, id_fe3o4
    use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
    use solver_data_module, only:  Fe2O3_pct_check_strain, &
           & ntemp_check_strain, nFe2O3_pct_check_strain

    ! Argument List

    ! Local Variables
    integer :: ioerror, i, id1, id2, id3, k, j
    logical :: no_oxide_namelist, fatal
    logical, save :: tube_found = .false.
    real, dimension(0:moxide) :: dummy1
    real, dimension(0:moxide)   :: weight_ave_fr
    real   :: denominator_add_ave
    character(LEN = 80), dimension(3) :: average_type

  ! START

  ! note; 
  ! fracture_toughness = sqrt(2.0 * surf_fracture_energy * Youngs_modul)

  ! redefine a ne material with average properties
  ! where is the defect located, within which layer of the scale, or 
  ! it spans on the interface between layers, or across multiple layers?

  ! index for the new material
  no_oxide_ave = no_layer + 1  ! changed from no_oxide + 1 

  ! exclude the tube; store beginning and ending indices
  id1 = 2
  id2 = no_layer
  id3 = no_oxide_ave

  if (id2 < id1)  then
    ! only the metal without oxide
    thickness_fr(no_oxide_ave) = 1.0
    material_name(no_oxide_ave) = 'average_scale'
    id_mat(id3) = id3  ! this is the last material
    rho(no_oxide_ave) = rho(1)
    poisson_ratio(no_oxide_ave) = poisson_ratio(1)
    surf_fracture_energy(no_oxide_ave) = 1.0
    fracture_toughness(no_oxide_ave) = 1.0
    cond_value(no_oxide_ave, :) = cond_value(1, :)
    cond_temp(no_oxide_ave, :) = cond_temp(1, :)
    th_exp_coeff(no_oxide_ave, :) = th_exp_coeff(1, :)
    th_exp_temp(no_oxide_ave, :) = th_exp_temp(1, :)
    Youngs_modul(no_oxide_ave, :) = Youngs_modul(1, :)
    Youngs_temp(no_oxide_ave, :) = Youngs_temp(1, :)
    
    RETURN

  endif


  ! when using average properties, the oxide scale needs fr=1
  thickness_fr(no_oxide_ave) = 1.0

  material_name(no_oxide_ave) = 'average_scale'
  id_mat(id3) = id3  ! this is the last material

  ! choose the largest number of properties
  weight_ave_fr(1:moxide) = thickness_fr(1:moxide)  ! still valid with new layer namelist
  ! weight_ave_fr(1:moxide) = 1.0  ! could correct by using "radius" as in weight for cylindrical geometry
  denominator_add_ave = SUM(weight_ave_fr(id1:id2))

  write(aux_lun, *) 'summ1 ', id1, id2, denominator_add_ave, id3, &
    & weight_ave_fr(id1:id2), thickness_fr(id1: id2)

  ! average type
  average_type(1) = 'add'
  average_type(2) = 'reciprocal_add'
  average_type(3) = 'geom'

  ! rho
  dummy1 = rho * weight_ave_fr
  rho(no_oxide_ave) = SUM(dummy1(id1:id2)) / denominator_add_ave

  ! poisson
  dummy1 = poisson_ratio * weight_ave_fr
  poisson_ratio(no_oxide_ave) = SUM(dummy1(id1:id2)) / denominator_add_ave

  ! surf_fracture
  dummy1 = surf_fracture_energy * weight_ave_fr
  surf_fracture_energy(no_oxide_ave) = SUM(dummy1(id1:id2)) / &
         & denominator_add_ave
  
  ! fracture toughness
  dummy1 = fracture_toughness * weight_ave_fr
  fracture_toughness(no_oxide_ave) = SUM(dummy1(id1:id2)) / &
         & denominator_add_ave

  call allocate_average_property(id1, id2, id3, weight_ave_fr, &
        & average_type(1), ncp_prop, cp_value, cp_temp)
  call allocate_average_property(id1, id2, id3, weight_ave_fr, &
        & average_type(1), ncond, cond_value, cond_temp)
  call allocate_average_property(id1, id2, id3, weight_ave_fr, &
        & average_type(1), nth_exp, th_exp_coeff, th_exp_temp)

  write(aux_lun, *) 'Young ', id1-1, nyoungs(1), Youngs_modul(1, 1:nyoungs(1)+1)
  write(aux_lun, *) 'Young ', id1-1, nyoungs(1), Youngs_temp(1, 1:nyoungs(1)+1)
  write(aux_lun, *) 'Young ', id1, nyoungs(id1), Youngs_modul(id1, 1:nyoungs(id1)+1)
  write(aux_lun, *) 'Young ', id1, nyoungs(id1), Youngs_temp(id1, 1:nyoungs(id1)+1)
  write(aux_lun, *) 'Young ', id2, nyoungs(id2), Youngs_modul(id2, 1:nyoungs(id2)+1)
  write(aux_lun, *) 'Young ', id2, nyoungs(id2), Youngs_temp(id2, 1:nyoungs(id2)+1)

  call allocate_average_property(id1, id2, id3, weight_ave_fr, &
        & average_type(1), nyoungs, Youngs_modul, Youngs_temp)

  write(aux_lun, *) 'cp ', id3, ncp_prop(id3), cp_value(id3, 1:ncp_prop(id3)+1)
  write(aux_lun, *) 'cp ', id3, ncp_prop(id3), cp_temp(id3, 1:ncp_prop(id3)+1)
  write(aux_lun, *) 'cond ', id3, ncond(id3), cond_value(id3, 1:ncond(id3)+1)
  write(aux_lun, *) 'cond ', id3, ncond(id3), cond_temp(id3, 1:ncond(id3)+1)
  write(aux_lun, *) 'id thermal_exp_oxide_ave ', id3, nth_exp(id3)
  
  do j = 1, nth_exp(id3) + 1
    write(aux_lun, 31) th_exp_temp(id3, j), th_exp_coeff(id3, j)
  end do
 31 format('th_exp_oxide_ave', 2(1x, 1pe13.6))

  ! write(aux_lun, *) 'th_exp ', id3, nth_exp(id3), th_exp_temp(id3, 1:nth_exp(id3)+1)
  write(aux_lun, *) 'Young ', id3, nyoungs(id3), Youngs_modul(id3, 1:nyoungs(id3)+1)
  write(aux_lun, *) 'Young ', id3, nyoungs(id3), Youngs_temp(id3, 1:nyoungs(id3)+1)

  ! allocate the other average properties for the two layer system

    ! identify Fe2O3
    id_fe2o3 = 0
    do k = 1, no_layer
      if (TRIM(material_name(id_mat(k))) == 'Fe2O3')  id_fe2o3 = k
    end do

    id_fe3o4 = 0
    ! identify Fe3o4
    do k = 1, no_layer
      if (TRIM(material_name(id_mat(k))) == 'Fe3O4')  id_fe3o4 = k
    end do

  ! exclude the tube; store beginning and ending indices
  ! these are unchanged; if id2 - id1 > 1 then do not do this
  id1 = 2
  id2 = no_layer

  if (id2 - id1 > 1)  return

  ! careful when fe2O3 is not defined
  if (id_fe2o3 == 0 .or. id_fe3o4 == 0)  then

    if (id_fe2o3 == 0)  write(out_lun, *) 'fe2o3 not defined'
    if (id_fe3o4 == 0)  write(out_lun, *) 'fe3o4 not defined'

    return

  end if

  do i = 1, nFe2O3_pct_check_strain

    ! index for the new material
    id3 = id3 + 1

    material_name(id3) = 'average_scale'// 'i'
    id_mat(id3) = id3  ! this is the last material

    ! choose the largest number of properties
    weight_ave_fr = 0.0
    weight_ave_fr(id_mat(id_fe2o3)) = Fe2O3_pct_check_strain(i) / 100.0
    weight_ave_fr(id_mat(id_fe3o4)) = 1.0 - weight_ave_fr(id_mat(id_fe2o3))
    denominator_add_ave = SUM(weight_ave_fr(id1:id2))
  
    write(out_lun, *) 'no_oxide_ave ', id3, no_oxide_ave

    ! rho
    dummy1 = rho * weight_ave_fr
    rho(id3) = SUM(dummy1(id1:id2)) / denominator_add_ave

    ! poisson
    dummy1 = poisson_ratio * weight_ave_fr
    poisson_ratio(id3) = SUM(dummy1(id1:id2)) / denominator_add_ave

    ! surf_fracture
    dummy1 = surf_fracture_energy * weight_ave_fr
    surf_fracture_energy(id3) = SUM(dummy1(id1:id2)) / &
         & denominator_add_ave
  
    ! fracture toughness
    dummy1 = fracture_toughness * weight_ave_fr
    fracture_toughness(id3) = SUM(dummy1(id1:id2)) / &
         & denominator_add_ave

    call allocate_average_property(id1, id2, id3, weight_ave_fr, &
        & average_type(1), ncp_prop, cp_value, cp_temp)
    call allocate_average_property(id1, id2, id3, weight_ave_fr, &
        & average_type(1), ncond, cond_value, cond_temp)
    call allocate_average_property(id1, id2, id3, weight_ave_fr, &
        & average_type(1), nth_exp, th_exp_coeff, th_exp_temp)
    call allocate_average_property(id1, id2, id3, weight_ave_fr, &
        & average_type(1), nyoungs, Youngs_modul, Youngs_temp)

  end do

  return

  END SUBROUTINE PROPERTY_SCALE_AVE 

  SUBROUTINE OPERAND_ARRAY(dim1, dim2, property, id1, id2, &
       & operation_type, value_out, location)

    !=======================================================================
    ! Purpose(s):
    !
    !   Calculate maximum of array property(id1:id2) and its location
    !  
    !=======================================================================

    use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

    ! Argument List 
    character(LEN = 80), intent(IN)     :: operation_type
    integer, intent(IN)     :: id1, id2, dim1, dim2
    real, dimension(dim1:dim2), intent(IN) :: property
    real, intent(OUT)    :: value_out
    integer, intent(OUT) :: location

    ! local variables
    integer, dimension(1)        :: id_max
    real, dimension(1:id2-id1+1) :: new_prop
    real :: prop1, prop2, var1
    integer :: i, j

    ! get the maximum number of the property
    new_prop(1:id2-id1+1) = property(id1:id2)
    ! nprop(id3) = MAXVAL(nprop(id1:id2))
    ! id_max = MAXLOC(nprop(id1:id2))

    if (TRIM(operation_type)== 'MAX' .or. &
      & TRIM(operation_type)== 'max' .or. &
      & TRIM(operation_type)== 'Max')  then

      value_out = MAXVAL(new_prop)
      id_max = MAXLOC(new_prop) + id1 - 1

    else if (TRIM(operation_type)== 'MIN' .or. &
      & TRIM(operation_type)== 'min' .or. &
      & TRIM(operation_type)== 'Min')  then

      value_out = MINVAL(new_prop)
      id_max = MINLOC(new_prop) + id1 - 1

    else if (TRIM(operation_type)== 'AVE_SUM' .or. &
      & TRIM(operation_type)== 'ave_sum' .or. &
      & TRIM(operation_type)== 'Ave_Sum')  then

      value_out = SUM(new_prop)/(id2-id1+1)
      id_max = 0

    else if (TRIM(operation_type)== 'AVE_REC_SUM' .or. &
      & TRIM(operation_type)== 'ave_rec_sum' .or. &
      & TRIM(operation_type)== 'Ave_Rec_Sum')  then

      value_out = (id2-id1+1) / SUM(1.0 / new_prop)
      id_max = 0

    else if (TRIM(operation_type)== 'AVE_PROD' .or. &
      & TRIM(operation_type)== 'ave_prod' .or. &
      & TRIM(operation_type)== 'Ave_Prod')  then

      value_out = PRODUCT(new_prop) ** (id2-id1+1)
      id_max = 0

    endif

    location = id_max(1)

    return

  END SUBROUTINE OPERAND_ARRAY

  SUBROUTINE OPERAND_2LAYER_ARRAY(dim1, dim2, property_one, property_two, &
     & id1_one, id2_one, id1_two, id2_two, operation_type, value_out, &
     & location)

    !=======================================================================
    ! Purpose(s):
    !
    !   Calculate maximum of array property(id1:id2) and its location
    !  
    !=======================================================================

    use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

    ! Argument List 
    character(LEN = 80), intent(IN)     :: operation_type
    integer, intent(IN)     :: id1_one, id2_one, id1_two, id2_two, dim1, dim2
    real, dimension(dim1:dim2), intent(IN) :: property_one, property_two
    real, intent(OUT)    :: value_out
    integer, intent(OUT) :: location

    ! local variables
    integer, dimension(1)        :: id_max
    real, dimension(1:id2_one-id1_one+2+id2_two-id1_two) :: new_prop
    real :: prop1, prop2, var1
    integer :: i, j, id_new_max

    ! get the maximum number of the property
    id_new_max = id2_one-id1_one+2+id2_two-id1_two

    ! store new array
    new_prop(1:id2_one-id1_one+1) = property_one(id1_one:id2_one)
    new_prop(id2_one-id1_one+2:id_new_max) = property_two(id1_two:id2_two)

    ! nprop(id3) = MAXVAL(nprop(id1:id2))
    ! id_max = MAXLOC(nprop(id1:id2))

    if (TRIM(operation_type)== 'MAX' .or. &
      & TRIM(operation_type)== 'max' .or. &
      & TRIM(operation_type)== 'Max')  then

      value_out = MAXVAL(new_prop)
      id_max = MAXLOC(new_prop) + id1_one - 1

    else if (TRIM(operation_type)== 'MIN' .or. &
      & TRIM(operation_type)== 'min' .or. &
      & TRIM(operation_type)== 'Min')  then

      value_out = MINVAL(new_prop)
      id_max = MINLOC(new_prop) + id1_one - 1

    else if (TRIM(operation_type)== 'AVE_SUM' .or. &
      & TRIM(operation_type)== 'ave_sum' .or. &
      & TRIM(operation_type)== 'Ave_Sum')  then

      value_out = SUM(new_prop) / id_new_max
      id_max = 0

    else if (TRIM(operation_type)== 'AVE_REC_SUM' .or. &
      & TRIM(operation_type)== 'ave_rec_sum' .or. &
      & TRIM(operation_type)== 'Ave_Rec_Sum')  then

      value_out = 1.0 * id_new_max / SUM(1.0 / new_prop)
      id_max = 0

    else if (TRIM(operation_type)== 'AVE_PROD' .or. &
      & TRIM(operation_type)== 'ave_prod' .or. &
      & TRIM(operation_type)== 'Ave_Prod')  then

      value_out = PRODUCT(new_prop)**id_new_max
      id_max = 0

    endif

    location = id_max(1)

    return

  END SUBROUTINE OPERAND_2LAYER_ARRAY

  SUBROUTINE INTEGRAND_VAR_LIMITS()

    !=======================================================================
    ! Purpose(s):
    !
    !   set up indices (id1:id2) to integrate a property
    ! or integrals of property=property(variable)
    !  
    !=======================================================================
  
  return

  END SUBROUTINE INTEGRAND_VAR_LIMITS


  SUBROUTINE OPERAND_TWO_ARRAY(dim1, dim2, property, variable, &
       & id1, id2, operation_type, value_out)

    !=======================================================================
    ! Purpose(s):
    !
    !   Calculate weighted averages property(id1:id2)
    ! or integrals of property=property(variable)
    !  
    !=======================================================================
   
    use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

    ! Argument List
    character(LEN = 80), intent(IN)     :: operation_type
    integer, intent(IN)     :: id1, id2, dim1, dim2
    real, dimension(dim1:dim2), intent(IN) :: property, variable
    real, intent(OUT)    :: value_out

    ! local variables
    real, dimension(1:id2-id1+1) :: new_prop, new_var
    real, dimension(1:id2-id1) :: new_int
    real :: prop1, prop2, var1
    integer :: i, j, nlocal

    ! number of values 
    nlocal = id2-id1+1

    ! get the maximum number of the property
    new_prop(1:id2-id1+1) = property(id1:id2)
    new_var(1:id2-id1+1) = variable(id1:id2)

    if  (TRIM(operation_type)== 'AVE_SUM' .or. &
      & TRIM(operation_type)== 'ave_sum' .or. &
      & TRIM(operation_type)== 'Ave_Sum')  then

      value_out = SUM(new_prop) / SUM(new_var)

    else if (TRIM(operation_type)== 'AVE_REC_SUM' .or. &
      & TRIM(operation_type)== 'ave_rec_sum' .or. &
      & TRIM(operation_type)== 'Ave_Rec_Sum')  then

      value_out = SUM(new_var) / SUM(new_var / new_prop)

    else if (TRIM(operation_type)== 'AVE_PROD' .or. &
      & TRIM(operation_type)== 'ave_prod' .or. &
      & TRIM(operation_type)== 'Ave_Prod')  then

      value_out = PRODUCT(new_prop) ** (id2-id1+1)

    else if (TRIM(operation_type)== 'INTEGRAL_VAR' .or. &
      & TRIM(operation_type)== 'integral_var' .or. &
      & TRIM(operation_type)== 'Integral_Var')  then

      new_int(1:nlocal-1) = 0.5 * (new_var(2:nlocal) -  &
         & new_var(1:nlocal-1)) * (new_prop(2:nlocal) + &
         & new_prop(1:nlocal-1))
      value_out = SUM(new_int)

    else if (TRIM(operation_type)== 'INTEGRAL_CONST' .or. &
      & TRIM(operation_type)== 'integral_const' .or. &
      & TRIM(operation_type)== 'Integral_Const')  then

      ! multiply later on with delta_x
      new_int(1:nlocal-1) = 0.5 * &
         & (new_prop(2:nlocal) + new_prop(1:nlocal-1))
      value_out = SUM(new_int)

    endif

    return

  END SUBROUTINE OPERAND_TWO_ARRAY

  
END MODULE PROPERTY_MODULE
