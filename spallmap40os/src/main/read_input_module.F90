MODULE READ_INPUT_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Read  namelists. 
    ! sabaua@ornl.gov read_input_module.F90  
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: READ_INPUT

 CONTAINS

  SUBROUTINE READ_PREINPUT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Read the simplified input file and construct the actual input file. 
    !
    !=======================================================================

  use output_module,        only: aux_lun, blank_line, &
                               &   inp_lun, input_file,   &
                               &   out_lun, Output_String, prefix, &
                               &   title, adv_lun, tty_lun,  &
                               &   WRITE_STRING
  use parameter_module, only: string_len

  implicit none

  ! Local Variables
  integer :: input_type

  rewind inp_lun

  ! read the input type
  read (inp_lun, '(a)') title

  Output_String = blank_line
  write (Output_String, 1) TRIM(title)
1 format(a)

  input_type = 0

  if (TRIM(title) == 'advanced' .or. &
     & TRIM(title) == 'ADVANCED')  then
    
    input_type = 0
    write(out_lun, *) 'input_file_type = ', 'advanced'

  else if (TRIM(title) == 'simple' .or. &
     & TRIM(title) == 'SIMPLE')  then

    input_type = 1
    write(out_lun, *) 'input_file_type = ', 'simple'

  else

    ! this could be old file
    input_type = -1
    write(out_lun, *) 'input_file_type = ', 'advanced - no keyword'

  endif

  ! advanced input file
  if (input_type == 0 .or. input_type == -1)  then
    write(out_lun, *) 'input file type: ADVANCED'
    ! if error on reading the input file when input_type == -1, then advise to correct it
    rewind inp_lun
    RETURN
  else if (input_type == 1)  then

    call WRITE_ADVANCED_INPUT()

  endif

  RETURN

  END SUBROUTINE READ_PREINPUT

  SUBROUTINE WRITE_ADVANCED_INPUT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Read the simplified input file and construct the actual input file. 
    !
    !=======================================================================

  use output_module,        only: aux_lun, blank_line, tmp_lun, &
                               &   inp_lun, input_file,   &
                               &   out_lun, Output_String, prefix, &
                               &   title, adv_lun, tty_lun,  &
                               &   tmp_lun, WRITE_STRING
  use output_module,            only: metal_property_file,  &
     & magnetite_property_file, &
     & hematite_property_file, spinel_property_file, &
     & oxide_growth_file, damage_map_file, oxide_layers_file, &
     & solver_options_file
  use parameter_module, only: string_len

  implicit none

  ! Local Variables
  integer :: l, n_lines, i, iend
  
  ! number of lines from the preinput file for each namelist
  integer :: oxide_growth_pre_no, boiler_pre_no 

  ! number of lines from the "input.txt" files
  integer :: oxide_layer_txt_no, oxide_growth_txt_no, damage_txt_no, &
           & property_txt_no, boiler_txt_no, solver_txt_no
  logical  :: fatal

  ! before calling here, the simple keyword was read

  ! Open the long input file
  input_file = Trim(prefix) // '_adv.inp'
  l = LEN_TRIM(input_file)
  OPEN (UNIT = adv_lun, FILE = input_file(:l), STATUS = 'unknown', &
             POSITION = 'rewind')

    ! read the names of the input files
    rewind inp_lun
    fatal = .false.
    ! Read oxide data
    call FILENAME_INPUT(fatal)
    if (fatal)  then
       write (out_lun, 21) fatal
       write (tty_lun, 21) fatal
21     format (/,9x,'FATAL: Error in FILENAME namelist fatal flag = ', L4/)
       stop
    endif

  rewind inp_lun !  does not position to the first line; close and reopen
 
  ! this should be the simple keyword
  read (inp_lun, '(a)') title

  ! read(inp_lun)  ! should be blank
  read (inp_lun, '(a)') title

  ! read the problem name
  read (inp_lun, '(a)') title

  Output_String = blank_line
  write (Output_String, 1) TRIM(title)
1 format(a)
  call WRITE_STRING (Output_String, adv_lun)

  ! read(inp_lun)  ! should be blank
  read (inp_lun, '(a)') title

  ! read alloy name
  read (inp_lun, '(a)') title

  if (TRIM(title) == 'T22' .or. TRIM(title) == 't22' .or. &
    & TRIM(title) == 'T91' .or. TRIM(title) == 't91' .or. &
    & TRIM(title) == 'TP304H' .or. TRIM(title) == 'tp304h' .or. &
    & TRIM(title) == 'TP304h'.or. &
    & TRIM(title) == '347H' .or. TRIM(title) == '347h' .or. &
    & TRIM(title) == '347HFG' .or. TRIM(title) == '347hfg')  then

    ! write empty line
    write(adv_lun, *) ''

    ! boiler_input.txt  fe2o3_input.txt  fecr3o4_input.txt
    ! oxide_layer_input.txt  t22_input.txt
    ! damage_input.txt  fe3o4_input.txt  oxide_growth_input.txt  solver_input.txt

    ! Open the oxide layer namelist files
    input_file = oxide_layers_file ! 'oxide_layer_input.txt'
    l = LEN_TRIM(input_file)
    OPEN (UNIT = tmp_lun, FILE = input_file(:l), STATUS = 'unknown', &
          & POSITION = 'rewind')
    do  ! n_lines = 28

      read (tmp_lun, '(a)', iostat=iend) title

      if (iend < 0) exit   ! end of file was reached

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

    close(tmp_lun)
  
    ! read empty line after alloy name
    read (inp_lun, '(a)') title
    write(adv_lun, *) ''

    ! read now the oxide growth namelist - the part from the simple file
    ! this is an incomplete namelist; does not have the "/" end of namelist 
    ! intentionally left out to indicate that additional file must be used
    n_lines = 6
    do i = 1, n_lines

      read (inp_lun, '(a)') title
      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do
    
    ! Open the oxide growth namelist files
    input_file = oxide_growth_file !  'oxide_growth_input.txt'
    l = LEN_TRIM(input_file)
    OPEN (UNIT = tmp_lun, FILE = input_file(:l), STATUS = 'unknown', &
          & POSITION = 'rewind')
    do ! n_lines = 14

      read (tmp_lun, '(a)', iostat=iend) title

      if (iend < 0) exit   ! end of file was reached

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

    close(tmp_lun)

    write(adv_lun, *) ''  ! write empty line

    ! Open the damage namelist file
    input_file = damage_map_file ! 'damage_input.txt'
    l = LEN_TRIM(input_file)
    OPEN (UNIT = tmp_lun, FILE = input_file(:l), STATUS = 'unknown', &
          & POSITION = 'rewind')
    do !  n_lines = 18

      read (tmp_lun, '(a)', iostat=iend) title

      if (iend < 0) exit   ! end of file was reached

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

    close(tmp_lun)

    write(adv_lun, *) ''  ! write empty line

    ! Open the property namelist files
    input_file = metal_property_file ! 't22_input.txt'
    ! n_lines = 21  without the rupture time data
    ! n_lines = 27  ! including rupture time data
    l = LEN_TRIM(input_file)
    OPEN (UNIT = tmp_lun, FILE = input_file(:l), STATUS = 'unknown', &
          & POSITION = 'rewind')
    do ! i = 1, n_lines

      read (tmp_lun, '(a)', iostat=iend) title

      if (iend < 0) exit   ! end of file was reached

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

    close(tmp_lun)

    write(adv_lun, *) ''  ! write empty line

   ! Open the property namelist files
    input_file = spinel_property_file !  'fecr3o4_input.txt'
    l = LEN_TRIM(input_file)
    OPEN (UNIT = tmp_lun, FILE = input_file(:l), STATUS = 'unknown', &
          & POSITION = 'rewind')
    do ! n_lines = 24

      read (tmp_lun, '(a)', iostat=iend) title

      if (iend < 0) exit   ! end of file was reached

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

    close(tmp_lun)

    write(adv_lun, *) ''  ! write empty line

   ! Open the property namelist files
    input_file = magnetite_property_file  !  'fe3o4_input.txt'
    l = LEN_TRIM(input_file)
    OPEN (UNIT = tmp_lun, FILE = input_file(:l), STATUS = 'unknown', &
          & POSITION = 'rewind')
    do ! n_lines = 24

      read (tmp_lun, '(a)', iostat=iend) title

      if (iend < 0) exit   ! end of file was reached

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

    close(tmp_lun)

    write(adv_lun, *) ''  ! write empty line

   ! Open the property namelist files
    input_file = hematite_property_file  ! 'fe2o3_input.txt'
    l = LEN_TRIM(input_file)
    OPEN (UNIT = tmp_lun, FILE = input_file(:l), STATUS = 'unknown', &
          & POSITION = 'rewind')
    do ! n_lines = 24

      read (tmp_lun, '(a)', iostat=iend) title

      if (iend < 0) exit   ! end of file was reached

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

    close(tmp_lun)

    write(adv_lun, *) ''  ! write empty line

    ! boiler namelist - the part from the simple file
    ! this is an incomplete namelist; does not have the "/" end of namelist 
    ! intentionally left out to indicate that additional file must be used
    n_lines = 27  ! without events
    n_lines = 36  ! with outage events
    n_lines = 38  ! with end of namelist
    n_lines = 40  ! end_state_boiler added
    do i = 1, n_lines

      read (inp_lun, '(a)') title

      Output_String = blank_line
      if (i >= 23 .and. i <=26)  then
        ! these should be the full load lines
        
        ! write (adv_lun, *) TRIM(title), ' 1999*0.0,'
        ! use this one when using old naames
        ! write (Output_String, 2) TRIM(title), '1999*0.0,'
        ! try this one for the weekly schedule and the daily schedule
        write (Output_String, 2) TRIM(title), '1998*0.0,'
        ! for new names we need to use this one
        ! write (Output_String, 1) TRIM(title)

2       format(a,1x, a)

      else 

        ! write (adv_lun, 1) TRIM(title)
        write (Output_String, 1) TRIM(title)
        ! call WRITE_STRING (Output_String, adv_lun)

      endif

      call WRITE_STRING (Output_String, adv_lun)

    end do

    ! solver namelist - the part from the simple file
    ! this is an incomplete namelist; does not have the "/" end of namelist 
    ! intentionally left out to indicate that additional file must be used
    n_lines = 5
    do i = 1, n_lines

      read (inp_lun, '(a)') title

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

   ! Open the solver namelist files
    input_file = solver_options_file  !  'solver_input.txt'
    n_lines = 35
    n_lines = 36  ! one line for update_sola_type
    l = LEN_TRIM(input_file)
    OPEN (UNIT = tmp_lun, FILE = input_file(:l), STATUS = 'unknown', &
          & POSITION = 'rewind')
    do

      read (tmp_lun, '(a)', iostat=iend) title

      if (iend < 0) exit   ! end of file was reached

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

    close(tmp_lun)

    ! output namelist - the part from the simple file
    ! this is an incomplete namelist; does not have the "/" end of namelist 
    ! intentionally left out to indicate that additional file must be used
    n_lines = 9
    do i = 1, n_lines

      read (inp_lun, '(a)') title

      Output_String = blank_line
      write (Output_String, 1) TRIM(title)
      call WRITE_STRING (Output_String, adv_lun)

    end do

  else

    write(out_lun, *) 'ERROR: alloy type not specified'
    write(out_lun, *) 'use T22; otherwise must generate the input files'
    stop

  endif

  ! at the end, switch adv_lun <--> inp_lun such that the new input file to be 
  ! actually read would be adv_lun
  tmp_lun = adv_lun
  adv_lun = inp_lun
  inp_lun = tmp_lun

  rewind inp_lun

  RETURN

  END SUBROUTINE WRITE_ADVANCED_INPUT

  SUBROUTINE READ_INPUT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Read all the namelists. 
    !
    !=======================================================================
    
    use parameter_module,         only: moxide
    use output_module,            only: tty_lun, out_lun, inp_lun
    use oxide_data_module,        only: no_oxide, no_layer, material_name, &
        & layer_name, no_metal_layers, first_oxide_layer, no_layers_oxide
    ! local variables
    logical  :: oxide_namelist, fatal, comb_namelist, diluent_namelist, &
                 & blade_namelist
    integer  :: ib, i

    ! generate the input file from the simplified input file, if needed
    call READ_PREINPUT ()

    fatal = .false.
    ! Read solver data
    call SOLVER_INPUT (fatal)
    if (fatal)  then
       write (out_lun, 29) fatal
       write (tty_lun, 29) fatal
29     format (/,9x,'FATAL: Error in SOLVER namelist fatal flag = ', L4/)
       stop
    endif

    rewind inp_lun
    no_oxide = 0
    fatal = .false.
    ! Read OXIDE data
    do ib = 1, moxide
       oxide_namelist = .false.
       call PROPERTY_INPUT (oxide_namelist)
       if (.not. oxide_namelist) exit
       ! fatal = .true.
       ! write (out_lun, 15) ib, fatal
       ! write (tty_lun, 15) ib, fatal
15     format (/,9x,'FATAL: Error in property namelist ',i2,' fatal flag = ', L4/)
    end do
    ! Make sure one property namelist has been read
    fatal = (no_oxide == 0)

    if (fatal)  then
       write (out_lun, 16) fatal
       write (tty_lun, 16) fatal
16     format (/,9x,'FATAL: Error in property namelist fatal flag = ', L4/)
       stop
    endif

       write (tty_lun, 18) no_oxide
       write (out_lun, 18) no_oxide
       write (tty_lun, *) 'materials ', (TRIM(material_name(i)),' ',  i = 1, no_oxide)
       write (out_lun, *) 'materials ', (TRIM(material_name(i)),' ', i = 1, no_oxide)
18     format (/' DONE Reading ',i2,' Property Namelists ...')

    rewind inp_lun
    no_layer = 0
    fatal = .false.
    ! Read OXIDE data
    do ib = 1, moxide
       oxide_namelist = .false.
       call LAYER_INPUT (oxide_namelist)
       if (.not. oxide_namelist) exit
       ! fatal = .true.
       ! write (out_lun, 15) ib, fatal
       ! write (tty_lun, 15) ib, fatal
31     format (/,9x,'FATAL: Error in layer namelist ',i2,' fatal flag = ', L4/)
    end do
    ! Make sure one property namelist has been read
    fatal = (no_layer == 0)

    if (fatal)  then
       write (out_lun, 32) fatal
       write (tty_lun, 32) fatal
32     format (/,9x,'FATAL: Error in layer namelist fatal flag = ', L4/)
       stop
    endif

       write (tty_lun, 33) no_layer
       write (out_lun, 33) no_layer
       write (tty_lun, *) 'layers ', (TRIM(layer_name(i)),' ',  i = 1, no_layer)
       write (out_lun, *) 'layers ', (TRIM(layer_name(i)),' ',  i = 1, no_layer)
33     format (/' DONE Reading ',i2,' Layer Namelists ...')

    if (first_oxide_layer == 0)  then
      ! only to deal with the case when no oxide is considered, i.e., no_metal_layers = 1
      first_oxide_layer = no_metal_layers + 1
    endif
    
    write(out_lun, *) 'end_metal_layers ', no_metal_layers, first_oxide_layer

    no_layers_oxide = no_layer - no_metal_layers

    rewind inp_lun
    fatal = .false.
    ! Read oxide data
    call OXIDE_INPUT (fatal)
    if (fatal)  then
       write (out_lun, 21) fatal
       write (tty_lun, 21) fatal
21     format (/,9x,'FATAL: Error in OXIDE_GROWTH namelist fatal flag = ', L4/)
       stop
    endif

    rewind inp_lun
    fatal = .false.
    ! Read damage data
    call DAMAGE_INPUT (fatal)
    if (fatal)  then
       write (out_lun, 22) fatal
       write (tty_lun, 22) fatal
22     format (/,9x,'FATAL: Error in DAMAGE namelist fatal flag = ', L4/)
       stop
    endif

    rewind inp_lun
    fatal = .false.
    ! Read oxide data
    call OUTPUT_INPUT (fatal)
    if (fatal)  then
       write (out_lun, 23) fatal
       write (tty_lun, 23) fatal
23     format (/,9x,'FATAL: Error in OUTPUT_INPUT namelist fatal flag = ', L4/)
       stop
    endif

    rewind inp_lun
    fatal = .false.
    ! Read compressor data
    call BOILER_INPUT (fatal)
    if (fatal)  then
       write (out_lun, 20) fatal
       write (tty_lun, 20) fatal
20     format (/,9x,'FATAL: Error in BOILER namelist fatal flag = ', L4/)
       stop
    endif

    ! moved before setting TEMP_GAS_p
    ! call CHECK_PARAM()

    return

  END SUBROUTINE READ_INPUT

  SUBROUTINE CHECK_PARAM()

  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer, temp_gas_pulse, &
      & temp_steam_pulse, temp_gas_idle, temp_steam_idle, &
      & heat_flux_fraction_idle, time_boiler_pulse, htc_tube_inner_partial
  use boiler_data_module, only: temp_gas_slope_p, temp_gas_p, Time_Boiler_p, &
      & TEMP_BOILER_data_no
  use oxide_data_module,  only: cond_value
  use output_module,            only: tty_lun, out_lun, inp_lun, &
      & aux_lun, gen_lun

  real :: rad_out, rad_in, dout, din, temp_max_full, temp_max_idle, a1, &
       & temp_ave_full, temp_ave_idle, temp_min_full, temp_min_idle, &
       & delta_temp_full, a_factor, b_factor, old_heat_flux_fr, &
       & temp_old_gas_idle
  real, dimension(2) :: heat_flux_full, heat_flux_idle
  integer :: i

  rad_out = tube_outer_radius
  rad_in = rad_out - tube_thickness

  a1 = LOG(rad_out / rad_in)
  dout = cond_value(1, 1) / (htc_tube_outer * rad_out)
  din = cond_value(1, 1) / (htc_tube_inner * rad_in)

  do i = 1, 2

    if (heat_flux_fraction_idle(i) > 0.0)  then
          
      if (temp_gas_idle(i) < 625.0)  then

        ! here the htc_inner and htc_outer are kept constant
        htc_tube_inner_partial = htc_tube_inner

            ! gas temperature adjusted to have a fraction of heat flux
            temp_old_gas_idle = temp_gas_idle(i)
            temp_gas_idle(i) = heat_flux_fraction_idle(i) * &
               & (temp_gas_pulse(i) - temp_steam_pulse(i)) + &
               & temp_steam_idle(i)

          write(out_lun, 53) temp_gas_idle(i), heat_flux_fraction_idle(i), &
               & i, temp_old_gas_idle, htc_tube_inner

          if (gen_lun > 0)  then

            write(gen_lun, 53) temp_gas_idle(i), heat_flux_fraction_idle(i), &
              & i, temp_old_gas_idle, htc_tube_inner
 53           format('new furnace_temperature_low_load = ', 1pe13.6, &
              & ' calculated for heat_flux_fraction ', &
              & 1pe13.6, 'at partial load type ', i2, &
              & ' old T=', 1pe13.6 ' same htc_steam=', 1pe13.6)

          endif

        else if (temp_gas_idle(i) >= 625.0)  then 

          ! htc_inner will be changed; temp_gas_idle will be held constant
          ! this is used instead of prescribing  htc_tube_inner_partial
          ! which will break the simplified input       

          a_factor = rad_in * (1.0 / (htc_tube_outer * rad_out) + &
            & a1 / cond_value(1, 1))
          b_factor = (temp_gas_idle(i) - temp_steam_idle(i)) / &
            & ((temp_gas_pulse(i) - temp_steam_pulse(i)) * &
            & heat_flux_fraction_idle(i))
          htc_tube_inner_partial = 1.0 / (b_factor / htc_tube_inner + &
            & (b_factor - 1.0) / a_factor)

            write(out_lun, 54) htc_tube_inner_partial, &
               & heat_flux_fraction_idle(i), i, htc_tube_inner, &
               & temp_gas_idle(i)

            if (gen_lun > 0)  then
              write(gen_lun, 54) htc_tube_inner_partial, &
               & heat_flux_fraction_idle(i), i, htc_tube_inner, &
               & temp_gas_idle(i)
 54           format('new htc_steam= ', 1pe13.6, &
              & ' calculated for heat_flux_fraction ', &
              & 1pe13.6, 'at partial load type ', i2, &
              & ' old htc_steam= ', 1pe13.6, &
              & ' same furnace_temperature_low_load= ', 1pe13.6)

            endif

            write(tty_lun, *) 'ERROR: if furnace_temperature_low_load > 625 C, ', &
            & ' set htc_tube_inner = htc_tube_inner_partial and maybe update oter variables '
            STOP

        endif

      else

        htc_tube_inner_partial = htc_tube_inner

            ! print the heat flux fraction
            if (time_boiler_pulse(i) > 0.0)  then

              old_heat_flux_fr = heat_flux_fraction_idle(i)
              heat_flux_fraction_idle(i) = (temp_gas_idle(i) - &
               & temp_steam_idle(i)) / &
               & (temp_gas_pulse(i) - temp_steam_pulse(i))

              write(out_lun, 4) heat_flux_fraction_idle(i), &
                 & temp_gas_idle(i), i, old_heat_flux_fr
              if (gen_lun > 0)  then
                write(gen_lun, 54) heat_flux_fraction_idle(i), &
                 & temp_gas_idle(i), i, old_heat_flux_fr
 4             format('new heat_flux_fraction = ', 1pe13.6, &
                 & ' calculated for furnace_temperature ', &
                 & 1pe13.6, 'at partial load type ', i2, &
                & ' old heat_flux_fraction=', 1pe13.6)

              endif

            endif

          endif

        end do

  do i = 1, 1

    din = cond_value(1, 1) / (htc_tube_inner * rad_in)
    heat_flux_full(i) = (temp_gas_pulse(i) - temp_steam_pulse(i)) * dout * &
      & htc_tube_outer / (a1 + dout + din)

    temp_max_full = temp_gas_pulse(i) - heat_flux_full(i) / htc_tube_outer
    temp_min_full = temp_gas_pulse(i) - (temp_gas_pulse(i) - &
             & temp_steam_pulse(i)) * (dout + a1) / &
             & (a1 + dout + din)
    temp_ave_full = 0.5 * (temp_max_full + temp_min_full)
    delta_temp_full = temp_max_full - temp_min_full

    din = cond_value(1, 1) / (htc_tube_inner_partial * rad_in)
    heat_flux_idle(i) = (temp_gas_idle(i) - temp_steam_idle(i)) * dout * &
      & htc_tube_outer / (a1 + dout + din)
    temp_max_idle = temp_gas_idle(i) - heat_flux_idle(i) / htc_tube_outer

    write(2, 1)  htc_tube_outer, htc_tube_inner, htc_tube_inner_partial, &
       & temp_gas_pulse(i), temp_steam_pulse(i), heat_flux_full(i), &
       & temp_max_full, temp_ave_full, temp_min_full, delta_temp_full, &
       & temp_gas_idle(i), temp_steam_idle(i), heat_flux_idle(i), &
       & temp_max_idle
 1  format('htc_htc_htc_TgFull_TstFull_hf_TgP_TstP_hf ', 20(1pe13.6, 1x))
  end do

  return

  END SUBROUTINE CHECK_PARAM

  SUBROUTINE BOILER_INPUT(fatal)
  
    use parameter_module,   only: pi, zero, mramp, mramp_out, npulse
    use input_utilities_module,   only: SEARCH_NML
    use output_module,            only: tty_lun, out_lun, inp_lun, &
      & aux_lun, gen_lun
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
          & temp_gas_first_idle, heat_flux_fraction_idle, &
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
    use boiler_data_module, only:   if_sliding_load, time_pulse_idle_load, &
      & load_fraction_pi_time, time_idle_pulse_load, &
      & load_fraction_ip_time, load_fraction_idle, &
      & load_temp_steam, temp_steam_vs_load, press_steam_vs_load, &
      & load_temp_gas, temp_gas_vs_load, end_state_boiler, &
      & press_gas_vs_load, heat_flux_fraction_vs_load, &
      & high_load_id_history, low_load_id_history, end_state_boiler
    use boiler_data_module, only: id_temp_op, id_press_op, id_flow_op
    use solution_data_module, only: no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low, time_full2low, &
      & time_full2low, temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low, &
      & no_f2l_event_out, no_event_point_full2low, &
      & no_event_full2low_output
    use boiler_operation_module, only: BOILER_SCHEDULE_TP_CYCLES

    ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    integer :: ioerror, i, k, j, j1, &
        & no_pct_first_pulse = 2, no_pct_pulse = 3, &
        & no_pct_idle = 3, no_pct_joint = 1, no_pct_last_pulse = 3, Status
    integer :: event_local_point
    integer, save :: no_event
    integer, save, dimension(200) :: event_start

    logical :: no_compr_namelist, compr_namelist

    real    :: Temp_ave, total_time, Temp_ave_local, Time_event_now
    ! variables for the simplified input
    real    :: time_full2low_load, time_low2full_load
    integer :: number_load_cycles
    integer, dimension(mramp_out) :: cycle_number_output_full2low

    real, dimension(2)    :: time_low_load, time_full_load, &
        & steam_temperature_low_load, &
        & heat_flux_fraction_low_load, furnace_temperature_low_load, &
        & pressure_steam_low_load, pressure_furnace_low_load
    real, dimension(npulse)    :: steam_temperature_full_load, &
        & furnace_temperature_full_load, &
        & pressure_steam_full_load, pressure_furnace_full_load

    ! input data for simplified input
  real  ::   outage_time_interval, time_outage_duration, &
        & steam_temperature_outage, furnace_temperature_outage, &
        & pressure_steam_outage, pressure_furnace_outage, &
        & time_full2outage, time_outage2full
  ! integer, dimension(10)    :: 

    namelist /boiler/  tube_thickness, tube_outer_radius, tube_curvature, &
      & htc_tube_inner, htc_tube_outer, temp_steam, time_boiler, &
      & no_total_pulse, time_boiler_idle, time_boiler_pulse, &
      & time_pulse_idle, temp_steam_pulse, temp_steam_idle, &
      & temp_steam_first_idle, time_boiler_first_idle, &
      & temp_steam_last_idle, time_boiler_last_idle, time_idle_pulse, &
      & temp_gas_pulse, temp_gas_idle, temp_gas_first_idle, &
      & temp_gas_last_idle, gas_flow_rate, press_inner_pulse, &
      & press_outer_pulse, press_inner_first_idle, &
      & press_outer_first_idle, press_inner_idle, press_outer_idle, &
      & press_inner_last_idle, press_outer_last_idle, &
      & no_pulse_full2low_output, no_event_full2low_output, &
      & event_time_interval, time_boiler_event, &
      & temp_steam_event, temp_gas_event, press_inner_event, &
      & press_outer_event, time_pulse_event, time_event_pulse, &
      & oxide_location, no_pulse_per_cycle, if_sliding_load, &
      & time_pulse_idle_load, heat_flux_fraction_idle, &
      & load_fraction_pi_time, time_idle_pulse_load, &
      & load_fraction_ip_time, load_fraction_idle, &
      & load_temp_steam, temp_steam_vs_load, press_steam_vs_load, &
      & load_temp_gas, temp_gas_vs_load, end_state_boiler, &
      & press_gas_vs_load, heat_flux_fraction_vs_load, &
      & time_low_load, time_full2low_load, time_full_load, &
      & time_low2full_load, steam_temperature_low_load, &
      & heat_flux_fraction_low_load, furnace_temperature_low_load, &
      & pressure_steam_low_load, pressure_furnace_low_load, &
      & number_load_cycles, steam_temperature_full_load, &
      & furnace_temperature_full_load, pressure_steam_full_load, &
      & pressure_furnace_full_load, cycle_number_output_full2low, &
      & outage_time_interval, time_outage_duration, &
      & steam_temperature_outage, furnace_temperature_outage, &
      & pressure_steam_outage, pressure_furnace_outage, &
      & time_full2outage, time_outage2full

       oxide_location = 'inner_metal_surface'
       ! for steiner problem oxide_location = 'outer_metal_surface'

       ! daily cycle by default
       no_pulse_per_cycle = 1

       ! default for pulse
       temp_steam_pulse = 0.0
       temp_gas_pulse = 0.0
       press_inner_pulse = 0.0
       press_outer_pulse = 0.0

       ! do not use heat flux frations
       heat_flux_fraction_idle = -1.0

       ! no sliding
       if_sliding_load = .false.

       ! event outage default
       no_event_full2low_output = -1

       ! default for simplified input
       time_boiler_last_idle = 0.0
       time_boiler_first_idle = 0.0

       ! simplified input
      time_low_load = -1.0
      time_full2low_load = -1.0
      time_full_load = -1.0
      time_low2full_load = -1.0

      steam_temperature_low_load = -1.0
      heat_flux_fraction_low_load = -1.0
      furnace_temperature_low_load = -1.0
      pressure_steam_low_load = -1.0
      pressure_furnace_low_load = -1.0

      number_load_cycles = -1
      cycle_number_output_full2low = -1

      steam_temperature_full_load = -1.0
      furnace_temperature_full_load = -1.0
      pressure_steam_full_load = -1.0
      pressure_furnace_full_load = -1.0
      
      ! simplified input for outage events
      outage_time_interval = -1.0
      time_outage_duration = -1.0
      steam_temperature_outage = -1.0
      furnace_temperature_outage = -1.0
      pressure_steam_outage = -1.0
      pressure_furnace_outage = -1.0
      time_full2outage = -1.0
      time_outage2full = -1.0
      
      end_state_boiler = 'none'

       ! Find namelist
       no_compr_namelist = .false.
       call SEARCH_NML (inp_lun, no_compr_namelist, 'boiler', 'BOILER')
       compr_namelist = .NOT. no_compr_namelist
       ! Read namelist if found
       if (compr_namelist) then
          read (inp_lun, NML= boiler, IOSTAT=ioerror)
          fatal = .not. (ioerror == 0) ! If read error, then didn't read namelist
          if (fatal)  then
            write (tty_lun, 25)
25          format (/' STOP: CANNOT Reading Boiler Namelist ...')
            stop
          else if (.not. fatal)  then
            write (tty_lun, 15) 
            write (out_lun, 15)
15          format (/' DONE Reading Boiler Namelist ...')
          endif

      tube_inner_radius = tube_outer_radius - tube_thickness
     
      ! simplified input
      if (time_low_load(1) > 0.0)  then
        time_boiler_idle = time_low_load
      endif

      if (time_full2low_load > 0.0)  then
        time_pulse_idle = time_full2low_load
      endif

      if (time_full_load(1) > 0.0)  then
        time_boiler_pulse = time_full_load
      endif

      if (time_low2full_load > 0.0)  then
        time_idle_pulse = time_low2full_load
      endif

      if (steam_temperature_low_load(1) > 0.0)  then
        temp_steam_idle = steam_temperature_low_load
      endif

      if (heat_flux_fraction_low_load(1) > 0.0)  then
        heat_flux_fraction_idle = heat_flux_fraction_low_load
      endif

      if (furnace_temperature_low_load(1) > 0.0)  then
        temp_gas_idle = furnace_temperature_low_load
      endif

      if (pressure_steam_low_load(1) > 0.0)  then
        press_inner_idle = pressure_steam_low_load
      endif

      if (pressure_furnace_low_load(1) > 0.0)  then
        press_outer_idle = pressure_furnace_low_load
      endif

      if (number_load_cycles > 0)  then
        no_total_pulse = number_load_cycles
      endif

      if (cycle_number_output_full2low(1) > 0)  then
        no_pulse_full2low_output = cycle_number_output_full2low
      endif

      if (steam_temperature_full_load(1) > 0.0)  then
        ! temp_steam_pulse(1:2) = steam_temperature_full_load(1:2)
        temp_steam_pulse = steam_temperature_full_load
      endif

      if (furnace_temperature_full_load(1) > 0.0)  then
        temp_gas_pulse = furnace_temperature_full_load
      endif

      if (pressure_steam_full_load(1) > 0.0)  then
        press_inner_pulse= pressure_steam_full_load
      endif

      if (pressure_furnace_full_load(1) > 0.0)  then
        press_outer_pulse = pressure_furnace_full_load
      endif        

      ! simplified input for outage events
      if (outage_time_interval > 0.0)  then
        event_time_interval = outage_time_interval
      endif

      if (time_outage_duration > 0.0)  then
        time_boiler_event = time_outage_duration
      endif

      if (steam_temperature_outage > 0.0)  then
        temp_steam_event = steam_temperature_outage
      endif

      if (furnace_temperature_outage > 0.0)  then
        temp_gas_event = furnace_temperature_outage
      endif

      if (pressure_steam_outage > 0.0)  then
        press_inner_event = pressure_steam_outage
      endif

      if (pressure_furnace_outage > 0.0)  then
        press_outer_event = pressure_furnace_outage
      endif

      if (time_full2outage > 0.0)  then
        time_pulse_event = time_full2outage
      endif

      if (time_outage2full > 0.0)  then
        time_event_pulse = time_outage2full
      endif

      if (time_boiler_first_idle <= 1.0e-6)  then
        ! in the old input file, the initial point was prescribed;
        ! for simplfied output, leave it out
        ! assign the start-up state to be same as to the partial load state
        time_boiler_first_idle = time_pulse_idle
        temp_steam_first_idle = temp_steam_idle(1)
        temp_gas_first_idle = temp_gas_idle(1)
        press_inner_first_idle = press_inner_idle(1)
        press_outer_first_idle = press_outer_idle(1)

      endif

      if (time_boiler_last_idle <= 1.0e-6)  then

        ! in the old input file, the last point was prescribed;
        if (TRIM(end_state_boiler) == 'partial_load')  then
          ! assign the last state to be same as to the partial load state

          time_boiler_last_idle = time_pulse_idle
          temp_steam_last_idle = temp_steam_idle(1)
          temp_gas_last_idle = temp_gas_idle(1)
          press_inner_last_idle = press_inner_idle(1)
          press_outer_last_idle = press_outer_idle(1)

        else if (TRIM(end_state_boiler) == 'outage')  then
          ! assign the last state to be same as to the outage state

          time_boiler_last_idle = 3.0
          temp_steam_last_idle = 25.0
          temp_gas_last_idle = 25.0
          press_inner_last_idle = 1.0e-1
          press_outer_last_idle = 1.0e-1

        else if (TRIM(end_state_boiler) == 'none')  then
          write(6, *) 'ERROR: end_state_boiler must be either partial_load or outage'
          STOP
        endif

      endif

      ! override the ampl data if lamp operates in pulse mode
      if (no_total_pulse == 0)  then
      ! count the number of data sets for AMPL 
        do i=100, 1, -1
          if (TEMP_BOILER_data_no == 0 .and. TEMP_STEAM(i) /= zero) &
            & TEMP_BOILER_data_no = i
        end do

      if (TEMP_BOILER_data_no == 0)  then
        write(out_lun,*) 'ERROR no table lookup for AMPL'
        stop
      else
       
        write(out_lun,*) 'allocate boiler cycle ', TEMP_BOILER_data_no
        ALLOCATE(TEMP_STEAM_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(Time_Boiler_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(TEMP_STEAM_Slope_p(1:TEMP_BOILER_data_no), STAT = Status)

        ! pressure cycle
        ALLOCATE(press_inner_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(press_outer_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(press_inner_slope_p(1:TEMP_BOILER_data_no), STAT = Status)
        ALLOCATE(press_outer_slope_p(1:TEMP_BOILER_data_no), STAT = Status)

        ! get the slope of the data
        do i= 1, TEMP_BOILER_data_no   ! no of data intervals
          TEMP_STEAM_p(i) = TEMP_STEAM(i)
          Time_Boiler_p(i) = Time_Boiler(i)
          ! 
          press_inner_p(i) = 1.0e+5
          press_outer_p(i) = 1.0e+5

        enddo

      endif

      else if (no_total_pulse > 0)  then
        
        do i = 1, no_total_pulse

          ! fill in with first pulse data
          if (temp_steam_pulse(i) <= 0.0)  then
            temp_steam_pulse(i) = temp_steam_pulse(1)
          endif

          if (temp_gas_pulse(i) <= 0.0)  then
            temp_gas_pulse(i) = temp_gas_pulse(1)
          endif
    
          if (press_inner_pulse(i) <= 0.0)  then
            press_inner_pulse(i) = press_inner_pulse(1)
          endif

          if (press_outer_pulse(i) <= 0.0)  then
            press_outer_pulse(i) = press_outer_pulse(1)
          endif

        end do

        ! set point in cycle at which printout will occur
        if (no_pulse_per_cycle == 1)  then
          high_load_id_history = 3
          low_load_id_history = 4
        else if (no_pulse_per_cycle == 1)  then
          ! they are the same
          high_load_id_history = 3
          low_load_id_history = 4
        endif

        ! moved before setting TEMP_GAS_p
        call CHECK_PARAM()
        call BOILER_SCHEDULE_TP_CYCLES(.not. if_sliding_load)

       end if

     endif

        ! get the slope of the data
        do i=1, TEMP_BOILER_data_no - 1  ! no of data intervals

          TEMP_STEAM_Slope_p(i) = (TEMP_STEAM_p(i+1) - TEMP_STEAM_p(i))/ &
            & (Time_Boiler_p(i+1) - Time_Boiler_p(i))
          TEMP_GAS_Slope_p(i) = (TEMP_GAS_p(i+1) - TEMP_GAS_p(i))/ &
            & (Time_Boiler_p(i+1) - Time_Boiler_p(i))

          press_inner_slope_p(i) = (press_inner_p(i+1) - press_inner_p(i))/ &
            & (Time_Boiler_p(i+1) - Time_Boiler_p(i))
          press_outer_slope_p(i) = (press_outer_p(i+1) - press_outer_p(i))/ &
            & (Time_Boiler_p(i+1) - Time_Boiler_p(i))

          write(out_lun, 209) i, Time_Boiler_p(i), TEMP_STEAM_p(i), &
           & TEMP_GAS_p(i), press_inner_p(i)
        enddo
       !  write(out_lun,209) i, Time_Boiler_p(i), &
       ! & TEMP_STEAM_p(i), TEMP_GAS_p(i)
 208    format('time_ampl ', 10(1pe13.6, 1x))
 209    format('time_ampl ', i4, 1x, 10(1pe13.6, 1x))

    ! get the average steam temperature
    ! get ramp-up temperature effect
    

  return

  END SUBROUTINE BOILER_INPUT

  SUBROUTINE OXIDE_INPUT(fatal)
  
  use parameter_module,     only : moxide_layer, moxide, mrate, zero, &
       & gas_constant_kjmk
  use input_utilities_module,   only: SEARCH_NML
  use oxide_data_module,       only: activ_energy_oxide, const_energy_oxide,  &
      & oxide_rate_value, oxide_rate_temp, nrate, &
      & oxide_rate_log, oxide_rate_1temp, oxide_thick_uf, &
      & oxide_thick_units, oxide_rate_units, temp_kelvin_unit, &
      & const_energy_uf, const_energy_units, &
      & flaw_defect_factor, flaw_oxide_size_factor, &
      & pilling_bedworth_ratio, pbr_bernstein_factor, &
      & if_temp_steam_or_oxide, lateral_ox_strain_length, &
      & temp_measure_ox_thickness
  use oxide_data_module,       only: log_const_energy_oxide, &
      & activ_energy_oxide_r, if_function_or_table_oxide, &
      & ox_hoop_strain_factor, ox_axial_strain_factor
  use oxide_stress_data_module,  only: if_th_oxide_growth, if_th_tube, &
      & if_exclude_th_strain_trace, if_growth_induce_stress

      ! oxide_fraction, &
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
 
    ! Argument List
    logical, intent(INOUT) :: fatal
    character(LEN = 80)    :: growth_induce_stress

    ! Local Variables
    integer :: ioerror, l, i
    logical :: no_solv_namelist, solv_namelist

    namelist /oxide_growth/  activ_energy_oxide, const_energy_oxide,  &
      & oxide_rate_value, oxide_rate_temp, &
      & oxide_thick_units, oxide_rate_units, const_energy_units, &
      & flaw_defect_factor, flaw_oxide_size_factor, &
      & pilling_bedworth_ratio, pbr_bernstein_factor, &
      & if_temp_steam_or_oxide, &
      & if_function_or_table_oxide, if_th_oxide_growth, if_th_tube, &
      & if_exclude_th_strain_trace, if_growth_induce_stress, &
      & lateral_ox_strain_length, growth_induce_stress, &
      & ox_hoop_strain_factor, ox_axial_strain_factor, &
      & temp_measure_ox_thickness

       ! default
      temp_kelvin_unit = 273.0
      flaw_defect_factor = 1.0
      ! by default grow oxide at steam temperature
      if_temp_steam_or_oxide = .true.
      if_function_or_table_oxide = .false.

      pbr_bernstein_factor = 1.0

      if_th_oxide_growth = .false.  ! no thermal expansion on h^2=Aexp(Q/RT)
      if_th_tube = .false.    ! no thermal expansion for the tube
      if_exclude_th_strain_trace = .true.
      if_growth_induce_stress = .false.
      growth_induce_stress = 'none'  ! used only in simplfied input

  ! for the oxidation/spallation, use k jules/(mol K)

       ! Find namelist
       no_solv_namelist = .false.
       call SEARCH_NML (inp_lun, no_solv_namelist, 'oxide_growth', 'OXIDE_GROWTH')
       solv_namelist = .NOT. no_solv_namelist

       ! Read namelist if found
       if (solv_namelist) then
          read (inp_lun, NML= oxide_growth, IOSTAT=ioerror)
          fatal = .not. (ioerror == 0) ! If read error, then didn't read namelist
          if (.not. fatal)  then
            write (tty_lun, 15) 
            write (out_lun, 15)
15          format (/' DONE Reading OXIDE_GROWTH Namelist ...')
          else
            write(tty_lun, *) 'STOP: OXIDE_GROWTH Input'
            stop            
          endif
       end if

       if (TRIM(growth_induce_stress) == 'YES' .or. &
         & TRIM(growth_induce_stress) == 'yes' .or. &
         & TRIM(growth_induce_stress) == 'Yes')  then
          
         if (.not. if_growth_induce_stress)  then
           if_growth_induce_stress = .true.
         endif
       endif

       if (TRIM(growth_induce_stress) == 'NO' .or. &
         & TRIM(growth_induce_stress) == 'no' .or. &
         & TRIM(growth_induce_stress) == 'No')  then
          
         if (if_growth_induce_stress)  then
           if_growth_induce_stress = .false.
         endif
       endif

       write(tty_lun, *) 'growth_induced ', TRIM(growth_induce_stress), &
          & if_growth_induce_stress

       do i = 1, 2
         ! write(aux_lun, *) 'oxlay ', oxide_fraction(i, 1:2)
       end do

       if (TRIM(const_energy_units) == 'kJ/mole')  then
         const_energy_uf = 1.0
       else
         write(out_lun, *) 'ERROR: activ energy units must be kJ/mole'
         stop
       endif

       if (TRIM(oxide_thick_units) == 'micron')  then

         if (TRIM(oxide_rate_units) == 'micron2_hour')  then
           oxide_thick_uf = 1.0
         else
           write(out_lun, *) 'ERROR oxide_rate_units should be micron2_hour'
           stop
         endif

       else if (TRIM(oxide_thick_units)== 'm')  then

         if (TRIM(oxide_rate_units) == 'm2_hour')  then
           oxide_thick_uf = 1.0e-6  ! convert to microns
         else
           write(out_lun, *) 'ERROR oxide_rate_units should be m2_hour'
           stop
         endif

       endif

       if (activ_energy_oxide /= zero)  then

         log_const_energy_oxide = log(const_energy_oxide)
         activ_energy_oxide_r = activ_energy_oxide / &
              & gas_constant_kjmk

         nrate = 35
         oxide_rate_1temp(1) = 1.0 / (temp_kelvin_unit + 20.0)
         ! oxide_rate_1temp(nrate) = 1.0 / (temp_kelvin_unit + 1220.0); used to be 720
         oxide_rate_1temp(nrate) = 1.0 / (temp_kelvin_unit + 1220.0)
         do i = 2, nrate - 1
           oxide_rate_1temp(i) = (oxide_rate_1temp(1) * (nrate - i) + &
              & oxide_rate_1temp(nrate)  * (i - 1)) / (nrate - 1)
         end do

         oxide_rate_log(1:nrate) = LOG(const_energy_oxide) - &
              & activ_energy_oxide * oxide_rate_1temp(1:nrate) / &
              & gas_constant_kjmk

         oxide_rate_value(1:nrate) = EXP(oxide_rate_log(1:nrate))

         oxide_rate_temp(1:nrate) = 1.0 / oxide_rate_1temp(1:nrate) - temp_kelvin_unit

         do i = 1, nrate
           ! write(6, *) 'exp ', oxide_rate_log(i), EXP(oxide_rate_log(i))
           write(aux_lun, 21) i, LOG(const_energy_oxide), &
              & activ_energy_oxide, oxide_rate_1temp(i), oxide_rate_temp(i), &
              & gas_constant_kjmk, &
              & LOG(const_energy_oxide) - &
              & activ_energy_oxide * oxide_rate_1temp(i) / &
              & gas_constant_kjmk, &
              & oxide_rate_value(i), &
              & 2.0*sqrt(oxide_rate_value(i))
   21  format('exp1 ', i2, 1x, 10(1pe13.6, 1x))
         end do

       else

         nrate = 0
         do i=mrate, 1, -1
           if (nrate == 0 .and. oxide_rate_value(i) /= zero) &
             & nrate = i-1
         end do

         if (nrate == 0)  then
           write(out_lun,*) 'ERROR no table lookup for OXIDATION RATE'
           stop
         else
           oxide_rate_log(1:nrate) = LOG(oxide_rate_value(1:nrate))
           oxide_rate_1temp(1:nrate) = 1.0 / &
              & (temp_kelvin_unit + oxide_rate_temp(1:nrate))
           ! use  linear interpolation between these two quantities to get 
           ! the oxidation rate at intermediate temperatures
         endif

      endif

      ! conversion to micron
      oxide_rate_value = oxide_thick_uf * oxide_rate_value
      oxide_rate_log(1:nrate) = LOG(oxide_rate_value(1:nrate))

      do i = 1, nrate
        write(out_lun, 209) i, oxide_rate_temp(i), oxide_rate_value(i), &
           & oxide_rate_1temp(i), oxide_rate_log(i)
 209    format('oxide_rate_val ', i4, 1x, 10(1pe13.6, 1x))
      end do

  return

  END SUBROUTINE OXIDE_INPUT

  SUBROUTINE DAMAGE_INPUT(fatal)
  
  use parameter_module,     only : moxide_layer, moxide, zero, &
       & gas_constant_kjmk, mdamage
  use input_utilities_module,   only: SEARCH_NML
  use oxide_data_module,       only: flaw_defect_factor, &
      & flaw_oxide_size_factor, amplitude_roughness, &
      & wavelength_roughness, delamination_size_existent, &
      & separated_area, porosity_fraction, &
      & number_cracks_volumetric, ndamage, &
      & delamination_size_factor
  use oxide_data_module,       only: flaw_size_init, flaw_name, nflaws, &
      & flaw_id_burried_crack, flaw_id_pore, flaw_id_interf_crack, &
      & flaw_id_surf_crack, flaw_id_delamin

      ! oxide_fraction, &
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
 
    ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    integer :: ioerror, l, i
    logical :: no_solv_namelist, solv_namelist

    namelist /damage/  flaw_defect_factor, flaw_oxide_size_factor, &
      & amplitude_roughness, &
      & wavelength_roughness, delamination_size_existent, &
      & separated_area, porosity_fraction, &
      & number_cracks_volumetric, delamination_size_factor, &
      & flaw_size_init, flaw_name

    ! defect size = flaw_oxide_size_factor * oxide_thickness
    ! amplitude_roughness == r (Schutze, 2005)
    ! wavelength_roughness == lambda 
    ! delamination_size_existent == R; the actual size is 2 * R 
    !   as in Schutze (2005)
    ! separated_area = A_sep
    !    A0 - entire surface area
    ! number_cracks_volumetric = N

      flaw_size_init = 0.0
      flaw_name = 'none'
      nflaws = 0

        ! Find namelist
       no_solv_namelist = .false.
       call SEARCH_NML (inp_lun, no_solv_namelist, 'damage', 'DAMAGE')
       solv_namelist = .NOT. no_solv_namelist

       ! Read namelist if found
       if (solv_namelist) then

          read (inp_lun, NML= damage, IOSTAT=ioerror)
          fatal = .not. (ioerror == 0) ! If read error, then didn't read namelist
          if (.not. fatal)  then
            write (tty_lun, 15) 
            write (out_lun, 15)
15          format (/' DONE Reading DAMAGE Namelist ...')
          else
            write(tty_lun, *) 'STOP: DAMAGE Input'
            stop            
          endif

          ndamage = mdamage

          LOOP1:  do i = 10, 1, -1

            if (nflaws == 0 .and. TRIM(flaw_name(i)) /= 'none')  then
              nflaws = i
              exit LOOP1
            endif

          end do LOOP1

          write(2, *) 'number of flaws ', nflaws
 
          do i = 1, nflaws

            if (TRIM(flaw_name(i)) == 'burried_crack')  &
               & flaw_id_burried_crack = i
            if (TRIM(flaw_name(i)) == 'surface_crack')  &
               & flaw_id_surf_crack = i
            if (TRIM(flaw_name(i)) == 'pore')  flaw_id_pore = i
            if (TRIM(flaw_name(i)) == 'interface_crack')  flaw_id_interf_crack = i

          end do

          write(2, *) 'defect types ', 'burried ', flaw_id_burried_crack
          write(2, *) 'defect types ', 'surface ', flaw_id_surf_crack
          write(2, *) 'defect types ', 'interface ', flaw_id_interf_crack
          write(2, *) 'defect types ', 'pore ', flaw_id_pore

       else

         write(2, *) 'damage namelist not found ' 
         STOP

       end if

  return

  END SUBROUTINE DAMAGE_INPUT

  SUBROUTINE OUTPUT_INPUT(fatal)
  
  use parameter_module,     only : moxide_layer
  use input_utilities_module,   only: SEARCH_NML
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun, &
     & strn_lun, enrg_lun, strs_lun, f2l_lun, gen_lun, f2l_strn_lun, &
     & enrg_dist_lun, enrg_dil_lun, strn_gen_lun, f2l_pl_lun
 
    ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    integer :: ioerror, l, i
    logical :: no_solv_namelist, solv_namelist
    character(LEN = 80)    :: strain_output_file, energy_output_file, &
      & outage_output_file, stress_output_file, general_data_file, &
      & outage_strain_output_file

    namelist /output/  strain_output_file, energy_output_file, &
      & outage_output_file, stress_output_file, general_data_file, &
      & outage_strain_output_file

       strain_output_file = 'none'
       energy_output_file = 'none'
       outage_output_file = 'none'
       stress_output_file = 'none'
       general_data_file = 'none'
       outage_strain_output_file = 'none' ! write the strain at outage

        ! Find namelist
       no_solv_namelist = .false.
       call SEARCH_NML (inp_lun, no_solv_namelist, 'output', 'OUTPUT')
       solv_namelist = .NOT. no_solv_namelist

       ! Read namelist if found
       if (solv_namelist) then

          read (inp_lun, NML= output, IOSTAT=ioerror)
          fatal = .not. (ioerror == 0) ! If read error, then didn't read namelist
          if (fatal)  then
            write (tty_lun, 25)
25          format (/' WARNING: CANNOT Reading OUTPUT Namelist ...')
            ! stop
          else if (.not. fatal)  then
            write (tty_lun, 15) 
            write (out_lun, 15)
15          format (/' DONE Reading OUTPUT Namelist ...')
          endif

         write(out_lun, *) 'read output files: ', &
           & TRIM(strain_output_file)//' ', &
           & TRIM(energy_output_file)//' ', TRIM(outage_output_file)//' ', &
           & TRIM(stress_output_file)//' ', &
           & TRIM(outage_strain_output_file)

         if (TRIM(strain_output_file) == 'none' .or. &
           & TRIM(strain_output_file) == '')  then
           write(out_lun, *) 'strain_output_file not open; use grep to get data'
           strn_lun = -1
         else

           l = LEN_TRIM(strain_output_file)
           OPEN (UNIT = strn_lun, FILE = strain_output_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'strain_output_file opened'

           l = LEN_TRIM(strain_output_file)
           OPEN (UNIT = strn_gen_lun, FILE = 'gen_'//strain_output_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'generating_strain_output_file opened'

         endif

         if (TRIM(energy_output_file) == 'none' .or. &
           & TRIM(energy_output_file) == '')  then
           write(out_lun, *) 'energy_output_file not open; use grep to get data'
           enrg_lun = -1
         else

           l = LEN_TRIM(energy_output_file)
           OPEN (UNIT = enrg_lun, FILE = energy_output_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'energy_output_file opened'

           l = LEN_TRIM(energy_output_file)
           OPEN (UNIT = enrg_dist_lun, FILE = 'dis_'//energy_output_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'energy_output_file opened'

           l = LEN_TRIM(energy_output_file)
           OPEN (UNIT = enrg_dil_lun, FILE = 'dil_'//energy_output_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'energy_output_file opened'

         endif

         if (TRIM(outage_output_file) == 'none' .or. &
           & TRIM(outage_output_file) == '')  then
           write(out_lun, *) 'outage_output_file not open; use grep to get data'
           f2l_lun = -1
           f2l_pl_lun = -1
         else

           l = LEN_TRIM(outage_output_file)
           OPEN (UNIT = f2l_lun, FILE = outage_output_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'outage_output_file opened'

           OPEN (UNIT = f2l_pl_lun, FILE = 'pl_'//outage_output_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'outage_output_file for flat_plate opened'

         endif

         if (TRIM(outage_strain_output_file) == 'none' .or. &
           & TRIM(outage_strain_output_file) == '')  then
           write(out_lun, *) 'outage_strain_output_file not open; use grep to get data'
           f2l_strn_lun = -1
         else

           l = LEN_TRIM(outage_strain_output_file)
           OPEN (UNIT = f2l_strn_lun, FILE = outage_strain_output_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'outage_strain_output_file opened'

         endif

         if (TRIM(stress_output_file) == 'none' .or. &
           & TRIM(stress_output_file) == '')  then
           write(out_lun, *) 'stress_output_file not open; use grep to get data'
           strs_lun = -1
         else

           l = LEN_TRIM(stress_output_file)
           OPEN (UNIT = strs_lun, FILE = stress_output_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'stress_output_file opened'

         endif

         if (TRIM(general_data_file) == 'none' .or. &
           & TRIM(general_data_file) == '')  then
           write(out_lun, *) 'general_data_file not open; use grep to get data'
           gen_lun = -1
         else

           l = LEN_TRIM(general_data_file)
           OPEN (UNIT = gen_lun, FILE = general_data_file(:l), &
             & STATUS = 'unknown', POSITION = 'rewind')
           write(out_lun, *) 'general_data_file opened'

         endif

       else

         write(2, *) 'output namelist not found ' 
           strn_lun = -1
           enrg_lun = -1
           f2l_lun = -1
           strs_lun = -1
           gen_lun = -1
         ! STOP

       end if

  return

  END SUBROUTINE OUTPUT_INPUT

  SUBROUTINE FILENAME_INPUT(fatal)
  
  use parameter_module,     only : moxide_layer
  use input_utilities_module,   only: SEARCH_NML
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun, &
     & adv_lun, metal_property_file, magnetite_property_file, &
     & hematite_property_file, spinel_property_file, &
     & oxide_growth_file, damage_map_file, oxide_layers_file, &
     & solver_options_file
 
    ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    integer :: ioerror, l, i
    logical :: no_solv_namelist, solv_namelist

    namelist /filenames/  metal_property_file, magnetite_property_file, &
     & hematite_property_file, spinel_property_file, &
     & oxide_growth_file, damage_map_file, oxide_layers_file, &
     & solver_options_file

     ! default options
     metal_property_file = 't22_input.txt'
     magnetite_property_file = 'fe3o4_input.txt'
     hematite_property_file = 'fe2o3_input.txt'
     spinel_property_file = 'fecr3o4_input.txt'
     oxide_growth_file = 'oxide_growth_input.txt'
     damage_map_file = 'damage_input.txt'
     oxide_layers_file = 'oxide_layer_input.txt'
     solver_options_file = 'solver_input.txt'

     ! Find namelist
     no_solv_namelist = .false.
     call SEARCH_NML (inp_lun, no_solv_namelist, 'filenames', 'FILENAMES')
          solv_namelist = .NOT. no_solv_namelist

     ! Read namelist if found
     if (solv_namelist) then

       read (inp_lun, NML= filenames, IOSTAT=ioerror)
       fatal = .not. (ioerror == 0) ! If read error, then didn't read namelist
       if (fatal)  then
         write (tty_lun, 25)
25       format (/' WARNING: CANNOT Reading FILENAMES Namelist ...')
            ! stop
       else if (.not. fatal)  then
         write (tty_lun, 15) 
         write (out_lun, 15)
15       format (/' DONE Reading FILENAMES Namelist ...')
       endif

         write(out_lun, *) 'read filenames files: ', &
           & TRIM(metal_property_file)//' ', &
           & TRIM(magnetite_property_file)//' ', &
           & TRIM(hematite_property_file)//' ', &
           & TRIM(spinel_property_file)//' ', &
           & TRIM(oxide_growth_file)//' ', &
           & TRIM(damage_map_file)//' ', &
           & TRIM(oxide_layers_file)//' ', &
           & TRIM(solver_options_file)

       else

         ! STOP

       end if

  return

  END SUBROUTINE FILENAME_INPUT

  SUBROUTINE SOLVER_INPUT(fatal)
  
  use parameter_module,     only : mstep, zero
  use input_utilities_module,   only: SEARCH_NML
  use solver_data_module,       only: max_iterat_solver, dTemp_eps, &
           & dt, time, dtime_save, time_save, ntime_save, &
           & no_thick_interval, no_temp_strain, temp_check_strain, &
           & Fe2O3_pct_check_strain, ntemp_check_strain, &
           & nFe2O3_pct_check_strain, thickness_spallation_map, &
           & no_thick_spall_map, interval_type_thick_spall_map, &
           & map_level, small_oxide_thickness, &
           & if_cylindrical, internal_htc, external_htc, &
           & no_temp_rad_points_tube, no_temp_rad_points_oxide, &
           & no_st_rad_points_tube, no_st_rad_points_oxide, &
           & Temp_reference_tube, Temp_reference_oxide, &
           & if_different_ref_temp, if_oxide_and_tube_temp, if_average_oxide
  use solver_data_module, only  : generalized_plane_strain, &
           & metal_recession_type, ratio_time_first_pulse, &
           & no_first_pulse_creep_int, start_creep_cycles, &
           & no_trans_intervals, no_intervals_f2l, &
           & no_intervals_l2f, update_sola_type
  use output_module,            only: tty_lun, out_lun, inp_lun, &
           & cea_lun, can_lun
  use creep_data_module, only: if_creep_ox_scale, if_creep_metal, &
           & stress_eq_epsilon, eps_oxide_grid, new_oxide_creep_option, &
           & no_intervals_now

    ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    integer :: ioerror, l, i
    logical :: no_solv_namelist, solv_namelist
    real    :: factor1, factor2, min1, max1
    character(LEN = 80)    :: creep_oxide_scale, creep_metal

    namelist /solver/  dTemp_eps, max_iterat_solver, dt, time, &
           & dtime_save, ntime_save, no_thick_interval, no_temp_strain, &
           & temp_check_strain, Fe2O3_pct_check_strain, &
           & thickness_spallation_map, update_sola_type, &
           & no_thick_spall_map, interval_type_thick_spall_map, &
           & map_level, small_oxide_thickness, &
           & if_cylindrical, internal_htc, external_htc, &
           & no_st_rad_points_tube, no_st_rad_points_oxide, &
           & generalized_plane_strain, Temp_reference_tube, &
           & Temp_reference_oxide, if_average_oxide, &
           & metal_recession_type, if_creep_ox_scale, if_creep_metal, &
           & stress_eq_epsilon, eps_oxide_grid, new_oxide_creep_option, &
           & ratio_time_first_pulse, creep_oxide_scale, creep_metal, &
           & no_first_pulse_creep_int, start_creep_cycles
 
       ! default
       if_average_oxide = .true.
       if_cylindrical = .true.
       internal_htc = .false.
       external_htc = .false.  ! since we do not have enough data
       ! axial strain is zero by default
       generalized_plane_strain = .false.
       Temp_reference_tube = -20.0
       Temp_reference_oxide = -20.0

       ! metal inner dimension is constant
       metal_recession_type = 'constant'

       ! no creep by default
       if_creep_ox_scale = .false.
       if_creep_metal = .false.
       creep_oxide_scale = 'none'
       creep_metal = 'none'

       stress_eq_epsilon = 0.1 ! 10% creep convergence; this is very large

       eps_oxide_grid = 0.01  ! this is 1% of grid size; could be too large
       new_oxide_creep_option = 'none'

       start_creep_cycles = 1
       ratio_time_first_pulse = 1.0
       no_first_pulse_creep_int = 1

       no_temp_strain = 0

       no_trans_intervals = 1

       ! update_sola_type = 0 by default
       update_sola_type = 0

       ! Find namelist
       no_solv_namelist = .false.
       call SEARCH_NML (inp_lun, no_solv_namelist, 'solver', 'SOLVER')
       solv_namelist = .NOT. no_solv_namelist

       ! Read namelist if found
       if (solv_namelist) then
          read (inp_lun, NML= solver, IOSTAT=ioerror)
          fatal = .not. (ioerror == 0) ! If read error, then didn't read namelist
          if (.not. fatal)  then
            write (tty_lun, 15) 
            write (out_lun, 15)
15          format (/' DONE Reading SOLVER Namelist ...')
          else
            write(tty_lun, *) 'STOP: SOLVER Input'
            stop            
          endif
       end if

       if (TRIM(creep_oxide_scale) == 'YES' .or. &
         & TRIM(creep_oxide_scale) == 'yes' .or. &
         & TRIM(creep_oxide_scale) == 'Yes')  then
          
         if (.not. if_creep_ox_scale)  then
           if_creep_ox_scale = .true.
         endif
       endif

       if (TRIM(creep_oxide_scale) == 'NO' .or. &
         & TRIM(creep_oxide_scale) == 'no' .or. &
         & TRIM(creep_oxide_scale) == 'No')  then
          
         if (if_creep_ox_scale)  then
           if_creep_ox_scale = .false.
         endif
       endif

       write(tty_lun, *) 'creep_oxide ', TRIM(creep_oxide_scale), &
         & if_creep_ox_scale

       if (TRIM(creep_metal) == 'YES' .or. &
         & TRIM(creep_metal) == 'yes' .or. &
         & TRIM(creep_metal) == 'Yes')  then
          
         if (.not. if_creep_metal)  then
           if_creep_metal = .true.
         endif
       endif

       if (TRIM(creep_metal) == 'NO' .or. &
         & TRIM(creep_metal) == 'no' .or. &
         & TRIM(creep_metal) == 'No')  then
          
         if (if_creep_metal)  then
           if_creep_metal = .false.
         endif
       endif

       write(tty_lun, *) 'creep_metal ', TRIM(creep_metal), &
         & if_creep_metal

       do i = 1, ntime_save
         time_save = dtime_save * (i - 1)
       end do

       ntemp_check_strain = 0
       do i = mstep, 1, -1
         if (ntemp_check_strain == 0 .and. temp_check_strain(i) /= zero)  &
           & ntemp_check_strain = i
       end do

       nFe2O3_pct_check_strain = 0
       do i = mstep, 1, -1
         if (nFe2O3_pct_check_strain == 0 .and. Fe2O3_pct_check_strain(i) /= zero)  &
           & nFe2O3_pct_check_strain = i
       end do

       if (no_thick_spall_map > 0)  then

         if (COUNT(thickness_spallation_map <= 1.0 .or. &
              & thickness_spallation_map >= 2000.0) == mstep)  then

     write(tty_lun, *) 'ERROR: input thickness_spallation_map values [2:999] microns'
           stop

         endif

         min1 = MINVAL(thickness_spallation_map, &
              & thickness_spallation_map > 1.0 .and. &
              & thickness_spallation_map < 2000.0)
         max1 = MAXVAL(thickness_spallation_map, &
              & thickness_spallation_map > 1.0 .and. &
              & thickness_spallation_map < 2000.0)

         if (interval_type_thick_spall_map == 'linear' .or. &
           & interval_type_thick_spall_map == 'LINEAR')  then

           do i = 1, no_thick_spall_map + 1
             thickness_spallation_map(i) = min1 + (max1 - min1) * (i - 1) / &
                & no_thick_spall_map
           end do

         else if (interval_type_thick_spall_map == 'log' .or. &
           & interval_type_thick_spall_map == 'LOG' .or. &
           & interval_type_thick_spall_map == 'logarithmic' .or. &
           & interval_type_thick_spall_map == 'LOGARITHMIC')  then

           min1 = log(min1)
           max1 = log(max1)
           do i = 1, no_thick_spall_map + 1
        thickness_spallation_map(i) = exp(min1 + (max1 - min1) * (i - 1) / &
                & no_thick_spall_map)
           end do

         else 

           write(tty_lun, *) 'ERROR: interval type must be linear or log'
           stop

         endif

       endif

      ! no_thick_interval - number of time steps per each cycle interval
      ! to get the oxide
      if (no_thick_interval <=0)  then
        write(tty_lun, *) 'ERROR: no_thick_interval must be > 0'
        stop
      endif

      ! no_temp_strain - number of temperature intervals to get the
      ! thermal strain, i.e., integrate alpha * dT

      ! no_temp_strain = number of points at which thermal expansion is checked
      if (ntemp_check_strain == 0)  then
        ntemp_check_strain = 2.0
        temp_check_strain(1) = 150.0
        temp_check_strain(2) = 550.0
      endif

      if (nFe2O3_pct_check_strain == 0) then
        Fe2O3_pct_check_strain(1) = 0.0
        Fe2O3_pct_check_strain(2) = 100.0
        nFe2O3_pct_check_strain = 2
      endif

      no_temp_rad_points_tube = 2 * no_st_rad_points_tube - 1
      no_temp_rad_points_oxide = 2 * no_st_rad_points_oxide - 1

      if (Temp_reference_tube > 0.0)  then 
        if_different_ref_temp = .true.
      else
        if_different_ref_temp = .false.
      endif

      if (Temp_reference_oxide > 0.0 .and. if_different_ref_temp)  then
        if_oxide_and_tube_temp = .true.
      else
        if_oxide_and_tube_temp = .false.
      endif

      if (if_creep_ox_scale .and. if_average_oxide)  then
        write(2, *) 'ERROR if_creep_ox_scale = if_average_oxide = .true.'
        write(2, *) 'both cannot be true'
        STOP
      endif

      if (if_creep_ox_scale .or. if_creep_metal)  then

        no_intervals_f2l = 1
        no_intervals_l2f = 0

        no_intervals_f2l = no_trans_intervals
        no_intervals_l2f = no_trans_intervals

      else

        no_intervals_f2l = 1
        no_intervals_l2f = 0

      endif

  ! define the number of total number of intervals for each segment type
  no_intervals_now = 1  ! full and partial load are done in one segment
  no_intervals_now(2) = no_intervals_l2f  ! 1-2 segment
  no_intervals_now(4) = no_intervals_f2l  ! 3-4 segment
  no_intervals_now(6) = no_intervals_l2f  ! 5-6 
  no_intervals_now(8) = no_intervals_f2l  ! 7-8
  
  if (no_intervals_f2l > 0 .and. no_intervals_l2f > 0)  then
    if (ABS(no_intervals_f2l - no_intervals_l2f)  > 0)  then
      write(tty_lun, *) 'ERROR; no_intervals_f2l and no_intervals_l2f must be the same'
      stop
    endif
  endif

  if (update_sola_type < 0 .or. update_sola_type > 11) then
    write(tty_lun, *) 'ERROR; update_sola_type must be 0, 1, 2, or 3'
    stop

  endif

  return

  END SUBROUTINE SOLVER_INPUT

  SUBROUTINE PROPERTY_INPUT (oxide_namelist)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read property namelist, put read data into place. At each call, 
    !                       data for oxide no_oxide is read.
    !   On succesful read, then no_oxide -> no_oxide + 1.
    !
    !=======================================================================
 
    use parameter_module,   only: moxide, mcp, zero, &
          & fract_tough2si, surf_energ2si
    use input_utilities_module,   only: SEARCH_NML
    use oxide_data_module,   only: material_name, cp_value_mat, cp_temp_mat, &
      & cond_value_mat, cond_temp_mat, rho_mat, rho_value_mat, rho_temp_mat, &
      & Youngs_modul_mat, Youngs_temp_mat, isubstrate, &
      & th_exp_coeff_mat, th_exp_temp_mat, no_oxide, &
      & nth_exp_mat, nyoungs_mat, ncond_mat, ncp_prop_mat, id_mat, &
      & poisson_ratio_mat, &
      & surf_fracture_energy_mat, fracture_toughness_mat
    use creep_data_module,  only: creep_const_mat, creep_exponent_mat, &
       & creep_activation_energy_mat, creep_activ_energy_stress_mat, &
       & Temp_creep_onset_mat
    use creep_data_module,  only: rupture_time_c1_mat, rupture_time_t1_mat, &
       & rupture_time_t1s1_mat, rupture_time_t1s2_mat, &
       & if_rupture_time_creep_mat, rupture_time_exp_mat
    use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

    ! Argument List
    logical, intent(INOUT)  :: oxide_namelist

    ! Local Variables
    integer :: ioerror, i
    logical :: no_oxide_namelist, fatal
    logical, save :: tube_found = .false.
    real, dimension(0:moxide, 0:moxide) :: dummy1

  ! START

  ! note; 
  ! fracture_toughness = sqrt(2.0 * surf_fracture_energy * Youngs_modul)

  namelist /property/  material_name, cp_value_mat, cp_temp_mat, &
      & cond_value_mat, cond_temp_mat, rho_mat, &
      & Youngs_modul_mat, Youngs_temp_mat, &
      & th_exp_coeff_mat, th_exp_temp_mat, &
      & poisson_ratio_mat, &
      & surf_fracture_energy_mat, fracture_toughness_mat, &
      & creep_const_mat, creep_exponent_mat, &
      & creep_activation_energy_mat, creep_activ_energy_stress_mat, &
      & Temp_creep_onset_mat, &
      & rupture_time_c1_mat, rupture_time_t1_mat, &
      & rupture_time_t1s1_mat, rupture_time_t1s2_mat, &
      & rupture_time_exp_mat

       ! initialize the variables at input
       material_name(0) = 'none'
       cp_value_mat(:, 0) = 0.0
       cp_temp_mat(:, 0) = 0.0
       cond_value_mat(:, 0) = 0.0
       cond_temp_mat(:, 0) = 0.0
       Youngs_modul_mat(:, 0) = 0.0
       Youngs_temp_mat(:, 0) = 0.0
       th_exp_coeff_mat(:, 0) = 0.0
       th_exp_temp_mat(:, 0) = 0.0
       rho_mat(0) = 0.0
       poisson_ratio_mat(0) = 0.0
       surf_fracture_energy_mat(0) = 0.0
       fracture_toughness_mat(0) = 0.0

       ! creep
       creep_const_mat(0) = 0.0
       creep_exponent_mat(0) = 0.0
       creep_activation_energy_mat(0) = 0.0
       creep_activ_energy_stress_mat(0) = 0.0
       Temp_creep_onset_mat(0) = 2000.0

       ! creep data using rupture time 
       if_rupture_time_creep_mat(0) = .false.
       rupture_time_c1_mat(0) = 0.0
       rupture_time_t1_mat(0) = 0.0
       rupture_time_t1s1_mat(0) = 0.0
       rupture_time_t1s2_mat(0) = 0.0
       rupture_time_exp_mat(0) = 0.0

       ! Find namelist
       no_oxide_namelist = .false.
       call SEARCH_NML (inp_lun, no_oxide_namelist, 'property', 'PROPERTY')
       oxide_namelist = .NOT. no_oxide_namelist

       ! Read namelist if found one
       if (oxide_namelist) then
          read (inp_lun, NML= property, IOSTAT=ioerror)
          oxide_namelist = (ioerror == 0) ! If read error, then didn't read namelist
       end if

    OXIDE_NML: if (oxide_namelist) then
       no_oxide = no_oxide + 1
       ! Read notice
       write (tty_lun, 15) no_oxide
       write (out_lun, 15) no_oxide
15     format (/' Reading Property Namelist number ',i2,' ...')
       ! Too many oxides gives a fatal error
       if (no_oxide > moxide) then
          write (tty_lun, 20) moxide
          write (out_lun, 20) moxide
20        format(/,9x,'FATAL: Exceeded maximum number: moxide = ',i5,'!',/)
          stop
       end if

       material_name(no_oxide) = material_name(0)
       rho_mat(no_oxide) = rho_mat(0)
       poisson_ratio_mat(no_oxide) = poisson_ratio_mat(0)

       Youngs_modul_mat(no_oxide, 1:moxide) = Youngs_modul_mat(0:moxide - 1, 0)
       Youngs_temp_mat(no_oxide, 1:moxide) = Youngs_temp_mat(0:moxide - 1, 0)
       nyoungs_mat(no_oxide) = 0
       do i = moxide, 1, -1
         if (nyoungs_mat(no_oxide) == 0 .and. Youngs_modul_mat(no_oxide, i) /= zero) &
             & nyoungs_mat(no_oxide) = i-1
       end do
       if (nyoungs_mat(no_oxide) == 0)  then
         write(out_lun, *) 'ERROR nyoungs = 0'
         stop
       endif

       if (surf_fracture_energy_mat(0) <= 0.0 .and. &
         & fracture_toughness_mat(0) <= 0.0)  then

         if (no_oxide == 1 .or. (TRIM(material_name(no_oxide)) == 'tube' .or. &
            & TRIM(material_name(no_oxide)) == 'TUBE' .or. &
            & TRIM(material_name(no_oxide)) == 'metal' .or. &
            & TRIM(material_name(no_oxide)) == 'METAL'))  then
            write(out_lun, *)  ' no surf_fracture_energy for metal '
         else
            write(out_lun, *) 'ERROR: surf_fracture_energy < 0 and ', &
               & ' fracture_toughness < 0 '
            STOP
         endif
       else if (surf_fracture_energy_mat(0) > 0.0 .and. &
         & fracture_toughness_mat(0) <= 0.0)  then
         fracture_toughness_mat(0) = sqrt(2.0 * &
            & surf_fracture_energy_mat(0) * Youngs_modul_mat(no_oxide, 1))
       else if (surf_fracture_energy_mat(0) <= 0.0 .and. &
         & fracture_toughness_mat(0) > 0.0)  then
         surf_fracture_energy_mat(0) = fracture_toughness_mat(0)**2 / &
         & (2.0 * Youngs_modul_mat(no_oxide, 1))
       endif

       surf_fracture_energy_mat(no_oxide) = surf_fracture_energy_mat(0)
       fracture_toughness_mat(no_oxide) = fracture_toughness_mat(0)

       ! convert to other units in order to accomodate 
       ! oxide thickness in microns
       surf_fracture_energy_mat(no_oxide) = surf_energ2si * &
           & surf_fracture_energy_mat(no_oxide) 
       fracture_toughness_mat(no_oxide) = fract_tough2si * &
           & fracture_toughness_mat(no_oxide)

       write(aux_lun, *) 'xc_mat ', TRIM(material_name(no_oxide)), &
          & poisson_ratio_mat(no_oxide), &
          & surf_fracture_energy_mat(no_oxide), &
          & fracture_toughness_mat(no_oxide)

       cp_value_mat(no_oxide, 1:moxide) = cp_value_mat(0:moxide - 1, 0)
       cp_temp_mat(no_oxide, 1:moxide) = cp_temp_mat(0:moxide - 1, 0)
       ncp_prop_mat(no_oxide) = 0
       do i = moxide, 1, -1
         if (ncp_prop_mat(no_oxide) == 0 .and. cp_value_mat(no_oxide, i) /= zero) &
             & ncp_prop_mat(no_oxide) = i-1
       end do
       if (ncp_prop_mat(no_oxide) == 0)  then
         write(out_lun, *) 'ERROR ncp_prop = 0'
         write(out_lun, *) 'cp not considered'
         ! stop
       endif
 
       cond_value_mat(no_oxide, 1:moxide) = cond_value_mat(0:moxide - 1, 0)
       cond_temp_mat(no_oxide, 1:moxide) = cond_temp_mat(0:moxide - 1, 0)
       ncond_mat(no_oxide) = 0
       do i = moxide, 1, -1
         if (ncond_mat(no_oxide) == 0 .and. cond_value_mat(no_oxide, i) /= zero) &
             & ncond_mat(no_oxide) = i-1
       end do
       if (ncond_mat(no_oxide) == 0)  then
         write(out_lun, *) 'ERROR ncond = 0'
         stop
       endif

       th_exp_coeff_mat(no_oxide, 1:moxide) = th_exp_coeff_mat(0:moxide - 1, 0)
       th_exp_temp_mat(no_oxide, 1:moxide) = th_exp_temp_mat(0:moxide - 1, 0)
       nth_exp_mat(no_oxide) = 0
       do i = moxide, 1, -1
         if (nth_exp_mat(no_oxide) == 0 .and. Th_exp_temp_mat(no_oxide, i) /= zero) &
             & nth_exp_mat(no_oxide) = i-1
       end do
       if (nth_exp_mat(no_oxide) == 0)  then
         write(out_lun, *) 'ERROR nth_exp = 0'
         stop
       endif

         write(aux_lun, *) 'Data for OXIDE ', TRIM(material_name(no_oxide)), &
       & cp_value_mat(no_oxide, 1), cp_value_mat(no_oxide, 2)
         do i = 1, no_oxide -1
           write(6, *) 'Data for OXIDE ', TRIM(material_name(no_oxide - i)), &
       & cp_value_mat(no_oxide -i, 1), cp_value_mat(no_oxide - i, 2)
         end do

       creep_const_mat(no_oxide) = creep_const_mat(0)
       creep_exponent_mat(no_oxide) = creep_exponent_mat(0)
       creep_activation_energy_mat(no_oxide) = creep_activation_energy_mat(0)
       creep_activ_energy_stress_mat(no_oxide) = &
           & creep_activ_energy_stress_mat(0)
       Temp_creep_onset_mat(no_oxide) = Temp_creep_onset_mat(0)

       if (ABS(rupture_time_exp_mat(0)) > 0.01)  then

         if_rupture_time_creep_mat(no_oxide) = .true.

         rupture_time_c1_mat(no_oxide) = rupture_time_c1_mat(0) 
         rupture_time_t1_mat(no_oxide) = rupture_time_t1_mat(0)
         rupture_time_t1s1_mat(no_oxide) = rupture_time_t1s1_mat(0)
         rupture_time_t1s2_mat(no_oxide) = rupture_time_t1s2_mat(0)
         rupture_time_exp_mat(no_oxide) = rupture_time_exp_mat(0)

       else
  
         if_rupture_time_creep_mat(no_oxide) = .false.

       endif

       write(aux_lun, *) 

    end if OXIDE_NML

    return

  END SUBROUTINE PROPERTY_INPUT 

  SUBROUTINE LAYER_INPUT (oxide_namelist)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read property namelist, put read data into place. At each call, 
    !                       data for oxide no_oxide is read.
    !   On succesful read, then no_layer -> no_layer + 1.
    !
    !=======================================================================
 
    use parameter_module,   only: moxide, mcp, zero, &
          & fract_tough2si, surf_energ2si
    use input_utilities_module,   only: SEARCH_NML
    use oxide_data_module,   only: layer_name, layer_materials_name, &
      & layer_mat_id, profile_function_ave, cp_value, cp_temp, &
      & cond_value, cond_temp, rho, rho_value, rho_temp, &
      & Youngs_modul, Youngs_temp, thickness_fr, isubstrate, &
      & th_exp_coeff, th_exp_temp, no_layer, no_layer_mat, &
      & nth_exp, nyoungs, ncond, ncp_prop, id_mat, &
      & poisson_ratio, fraction_material, thickness_fr_test, &
      & surf_fracture_energy, fracture_toughness
    use oxide_data_module,   only: material_name, cp_value_mat, cp_temp_mat, &
      & cond_value_mat, cond_temp_mat, rho_mat, rho_value_mat, rho_temp_mat, &
      & Youngs_modul_mat, Youngs_temp_mat, isubstrate, &
      & th_exp_coeff_mat, th_exp_temp_mat, no_oxide, &
      & nth_exp_mat, nyoungs_mat, ncond_mat, ncp_prop_mat, id_mat, &
      & poisson_ratio_mat, growth_ox_temperature, thickness_growth_fr, &
      & surf_fracture_energy_mat, fracture_toughness_mat, &
      & first_oxide_layer, no_metal_layers, radius_layer_inside, &
      & radius_layer_outside
    use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
    use oxide_stress_data_module,  only: oxidation_strain_mode, &
      & growth_ox_location, displacement_jump_type, &
      & displacement_jump_fr
    use creep_data_module, only: creep_const_mat, creep_exponent_mat, &
       & creep_activation_energy_mat, creep_activ_energy_stress_mat, &
       & creep_activ_energy_stress, Temp_creep_onset_mat, &
       & creep_const, creep_const_ln, creep_exponent, &
       & creep_activation_energy, Temp_creep_onset
    use creep_data_module,  only: rupture_time_c1_mat, rupture_time_t1_mat, &
       & rupture_time_t1s1_mat, rupture_time_t1s2_mat, &
       & if_rupture_time_creep_mat, rupture_time_c1, &
       & rupture_time_t1, rupture_time_exp, rupture_time_exp_mat, &
       & rupture_time_t1s1, rupture_time_t1s2, &
       & if_rupture_time_creep

    ! Argument List
    logical, intent(INOUT)  :: oxide_namelist

    ! Local Variables
    integer :: ioerror, i, j, ilay, k, k_tr
    logical :: no_oxide_namelist, fatal
    logical, save :: tube_found = .false.
    real, dimension(0:moxide, 0:moxide) :: dummy1

  ! START

  ! note; 
  ! fracture_toughness = sqrt(2.0 * surf_fracture_energy * Youngs_modul)

  namelist /oxide_layer/  layer_name, thickness_fr, &
      & profile_function_ave, layer_materials_name, fraction_material, &
      & oxidation_strain_mode, growth_ox_location, displacement_jump_type, &
      & displacement_jump_fr, growth_ox_temperature, thickness_growth_fr, &
      & radius_layer_inside

  ! it would be ideal for the growth_temperature to be based on growth_ox_location
      ! fraction_mat_dmin, fraction_mat_dmax

       ! initialize the variables at input
       layer_name(0) = 'none'
       ! profile_function_ave(0) = 'uniform'  ! constant properties
       profile_function_ave(0) = 'linear'  ! constant properties
       layer_materials_name = 'none'
       oxidation_strain_mode(0) = 'none'
       growth_ox_location(0) = 'none'
       growth_ox_temperature(0) = 'oxide_surface'

       displacement_jump_type(0) = 'none'

       thickness_fr(0) = 0.0
       thickness_growth_fr(0) = 0.0
       fraction_material(:, 0) = 0.0
       displacement_jump_fr(0) = 0.0

       radius_layer_inside(0) = -1.0

       ! Find namelist

       no_oxide_namelist = .false.
       call SEARCH_NML (inp_lun, no_oxide_namelist, &
             & 'oxide_layer', 'OXIDE_LAYER')
       oxide_namelist = .NOT. no_oxide_namelist

       ! Read namelist if found one
       if (oxide_namelist) then
          read (inp_lun, NML= oxide_layer, IOSTAT=ioerror)
          oxide_namelist = (ioerror == 0) ! If read error, then didn't read namelist
       end if

    OXIDE_NML: if (oxide_namelist) then
       no_layer = no_layer + 1
       ! Read notice
       write (tty_lun, 15) no_layer
       write (out_lun, 15) no_layer
15     format (/' Reading OXIDE_LAYER Namelist number ',i2,' ...')
       ! Too many oxides gives a fatal error
       if (no_layer > moxide) then
          write (tty_lun, 20) moxide
          write (out_lun, 20) moxide
20        format(/,9x,'FATAL: Exceeded maximum number: moxide = ',i5,'!',/)
          stop
       end if

       ! write(out_lun, *) 'start read metal_layers ', no_layer, &
       !     & no_metal_layers, first_oxide_layer, &
       !     & radius_layer_inside(0)

       ! initialize the number of metal layers
       if (no_layer == 1)  no_metal_layers = 1

       ! if there is only one metal layer assign the radius
       if (no_layer == 2 .and. radius_layer_inside(0) <= 0.0)  then
         first_oxide_layer = 2
       endif

       if (no_layer > 1 .and. radius_layer_inside(0) > 0.0) then
         ! this is an additional metal layer
         no_metal_layers = no_metal_layers + 1
         radius_layer_outside(no_metal_layers) = &
            & radius_layer_inside(no_metal_layers - 1)
         ! increment the id of the first oxide layer
         first_oxide_layer = no_metal_layers + 1
         write(out_lun, *) 'read metal_layers ', no_metal_layers, &
            & first_oxide_layer
       endif

       if (radius_layer_inside(0) > 0.0)  then
         ! no_metal_layers = no_layer at this point
         radius_layer_inside(no_layer) = radius_layer_inside(0)
       endif

       ! 2009 prescribe thickness_fr = thickness_growth_fr if missing
       if (thickness_fr(0) == 0.0)  then
         thickness_fr(0) = thickness_growth_fr(0)
         write(out_lun, *) 'thickness_fr(0) = thickness_growth_fr(0) 2009 change'
       endif

       write(aux_lun, *) 'th1 ', no_layer, thickness_fr(1:no_layer)

       layer_name(no_layer) = layer_name(0)
       layer_materials_name(1:moxide) = layer_materials_name(0:moxide-1)

       oxidation_strain_mode(no_layer) = oxidation_strain_mode(0)
       growth_ox_location(no_layer) = growth_ox_location(0)
       growth_ox_temperature(no_layer) = growth_ox_temperature(0)

       displacement_jump_type(no_layer) = displacement_jump_type(0)
       displacement_jump_fr(no_layer) = displacement_jump_fr(0)

       no_layer_mat(no_layer) = -1
       do j = moxide, 1, -1
         LOOP_MAT: do i = 1, no_oxide
           if (TRIM(layer_materials_name(j)) == TRIM(material_name(i)) .and. &
             no_layer_mat(no_layer) == -1)  then
             no_layer_mat(no_layer) = j
             exit LOOP_MAT
           endif
           ! write(out_lun, *) j, i, TRIM(layer_materials_name(j)), &
           !  & TRIM(material_name(i))
 ! 31        format('layer_mat ', 2(i2, 1x), 2('(a)'))
         end do LOOP_MAT

       end do

       if (no_layer_mat(no_layer) == -1)  then

         write(out_lun, *) 'ERROR no compounds for this layer = 0'
         stop
       endif

       ! assign new material id; each layer is made now of a new material 
       id_mat(no_layer) = no_layer

       ! identify materials that make up this layer
       do j = 1, no_layer_mat(no_layer)
         layer_mat_id(no_layer, j) = 0
         do i = 1, no_oxide
           if (TRIM(layer_materials_name(j)) == TRIM(material_name(i)))  then
             layer_mat_id(no_layer, j) = i
           endif
         end do

         if (layer_mat_id(no_layer, j) == 0)  then
           write(out_lun, *) 'ERROR no compounds name for this layer'
           stop
         endif

       end do
         
       thickness_fr(no_layer) = thickness_fr(0)
       thickness_fr_test(no_layer) = thickness_fr(no_layer)
       thickness_growth_fr(no_layer) = thickness_growth_fr(0)
       fraction_material(no_layer, 1:moxide) = fraction_material(0:moxide - 1, 0)

       ! start with uniform properties
       ! as compound fraction changes in the layer, the layer properties
       ! are location dependent; the methdology is for constant layer
       ! properties; use different layers with different void fraction
       ! thus, at this time accept only constant layer properties
       ! FE or FV needs to be used when location dependent properties are used

       ! th_exp_temp(no_oxide, 1:moxide) = th_exp_temp(0:moxide - 1, 0)

       ! temperature intervals do not change
       nth_exp(no_layer) = nth_exp_mat(layer_mat_id(no_layer, 1))
       nyoungs(no_layer) = nyoungs_mat(layer_mat_id(no_layer, 1))
       ncond(no_layer) = ncond_mat(layer_mat_id(no_layer, 1))
       ncp_prop(no_layer) = ncp_prop_mat(layer_mat_id(no_layer, 1))

       cp_temp(no_layer, 1:moxide) = cp_temp_mat( &
                 & layer_mat_id(no_layer, 1), 1:moxide)
       cond_temp(no_layer, 1:moxide) = cond_temp_mat( &
                 & layer_mat_id(no_layer, 1), 1:moxide)
       Youngs_temp(no_layer, 1:moxide) = Youngs_temp_mat( &
                 & layer_mat_id(no_layer, 1), 1:moxide)
       th_exp_temp(no_layer, 1:moxide) = th_exp_temp_mat( &
                 & layer_mat_id(no_layer, 1), 1:moxide)

       ! do only one average for the time being; the radial effect 
       ! the volumetric fraction may vary with the radius

       ilay = no_layer

       if_rupture_time_creep(no_layer) = .false.

       if (no_layer_mat(no_layer) == 1)  then

         k = layer_mat_id(no_layer, 1)  ! compound id as read in property

         ! no need for averaging when only one compound
         rho(ilay) = rho_mat(k)
         poisson_ratio(ilay) = poisson_ratio_mat(k) 
         surf_fracture_energy(ilay) = surf_fracture_energy_mat(k)
         fracture_toughness(ilay) = fracture_toughness_mat(k)

         cp_value(ilay, 1:ncp_prop(ilay)) = cp_value_mat(k, 1:ncp_prop(ilay))
         cond_value(ilay, 1:ncond(ilay)) = cond_value_mat(k, 1:ncond(ilay))
         Youngs_modul(ilay, 1:nyoungs(ilay)) = Youngs_modul_mat(k, 1:nyoungs(ilay))
         th_exp_coeff(ilay, 1:nth_exp(ilay)) = th_exp_coeff_mat(k, 1:nth_exp(ilay))

         creep_const(ilay) = creep_const_mat(k)
         creep_exponent(ilay) = creep_exponent_mat(k)
         creep_activation_energy(ilay) = creep_activation_energy_mat(k)
         creep_activ_energy_stress(ilay) = creep_activ_energy_stress_mat(k)
         Temp_creep_onset(ilay) = Temp_creep_onset_mat(k)

           ! setup creep data as a function of rupture time
           if (if_rupture_time_creep_mat(k))  then
             
             rupture_time_c1(ilay) = rupture_time_c1_mat(k)
             rupture_time_t1(ilay) = rupture_time_t1_mat(k)
             rupture_time_t1s1(ilay) = rupture_time_t1s1_mat(k)
             rupture_time_t1s2(ilay) = rupture_time_t1s2_mat(k)
             rupture_time_exp(ilay) = rupture_time_exp_mat(k)

             if_rupture_time_creep(no_layer) = .true.
             write(aux_lun, *) 'rupture_time_creep considered for layer ', &
               & no_layer

           endif

       else if (no_layer_mat(no_layer) > 1)  then

         ! initialize
         rho(ilay) = 0.0
         poisson_ratio(ilay) = 0.0
         surf_fracture_energy(ilay) = 0.0
         fracture_toughness(ilay) = 0.0

         creep_const(ilay) = 0.0
         creep_exponent(ilay) = 0.0
         creep_activation_energy(ilay) = 0.0
         creep_activ_energy_stress(ilay) = 0.0
         Temp_creep_onset(ilay) = 0.0

         cp_value(ilay, 1:ncp_prop(ilay)) = 0.0
         cond_value(ilay, 1:ncond(ilay)) = 0.0
         Youngs_modul(ilay, 1:nyoungs(ilay)) = 0.0
         th_exp_coeff(ilay, 1:nth_exp(ilay)) = 0.0

         k_tr = 0  ! counts how many material have rupture time creep functions

         do j = 1, no_layer_mat(no_layer)

           k = layer_mat_id(no_layer, j)  ! compound id as read in property
 
           rho(ilay) = rho(ilay) + rho_mat(k) * &
             & fraction_material(ilay, j)
           poisson_ratio(ilay) = poisson_ratio(ilay) + &
             & poisson_ratio_mat(k) * &
             & fraction_material(ilay, j)
           surf_fracture_energy(ilay) = surf_fracture_energy(ilay) + &
             & surf_fracture_energy_mat(k) * &
             & fraction_material(ilay, j)
           fracture_toughness(ilay) = fracture_toughness(ilay) + &
             & fracture_toughness_mat(k) * &
             & fraction_material(ilay, j)

           creep_const(ilay) = creep_const(ilay) + creep_const_mat(k) * &
             & fraction_material(ilay, j)
           creep_exponent(ilay) = creep_exponent(ilay) + &
             & creep_exponent_mat(k) * &
             & fraction_material(ilay, j)
           creep_activation_energy(ilay) = creep_activation_energy(ilay) + &
             & creep_activation_energy_mat(k) * &
             & fraction_material(ilay, j)
           creep_activ_energy_stress(ilay) = creep_activ_energy_stress(ilay) + &
             & creep_activ_energy_stress_mat(k) * &
             & fraction_material(ilay, j)
           Temp_creep_onset(ilay) = Temp_creep_onset(ilay) + &
             & Temp_creep_onset_mat(k) * &
             & fraction_material(ilay, j)

           ! store the data for the mixture, creep rate parameters are not individually averaged
           ! not needed, just use the k = layer_mat_id(no_layer, j) in DCREEP
           ! creep_exponent_mix(ilay, k) = creep_exponent_mat(k)

           ! setup creep data as a function of rupture time
           if (if_rupture_time_creep_mat(k))  then
             
             rupture_time_c1(ilay) = rupture_time_c1(ilay) + &
               & rupture_time_c1_mat(k) * &
               & fraction_material(ilay, j)
             rupture_time_t1(ilay) = rupture_time_t1(ilay) + &
               & rupture_time_t1_mat(k) * &
               & fraction_material(ilay, j)
             rupture_time_t1s1(ilay) = rupture_time_t1s1(ilay) + &
               & rupture_time_t1s1_mat(k) * &
               & fraction_material(ilay, j)
             rupture_time_t1s2(ilay) = rupture_time_t1s2(ilay) + &
               & rupture_time_t1s2_mat(k) * &
               & fraction_material(ilay, j)
             rupture_time_exp(ilay) = rupture_time_exp(ilay) + &
               & rupture_time_exp_mat(k) * &
               & fraction_material(ilay, j)

             k_tr = k_tr + 1

           endif

           cp_value(ilay, 1:ncp_prop(ilay)) = &
             & cp_value(ilay, 1:ncp_prop(ilay)) + &
             & cp_value_mat(k, 1:ncp_prop(ilay)) * &
             & fraction_material(ilay, j)
           cond_value(ilay, 1:ncond(ilay)) = &
             & cond_value(ilay, 1:ncond(ilay)) + &
             & cond_value_mat(k, 1:ncond(ilay)) * &
             & fraction_material(ilay, j)
           Youngs_modul(ilay, 1:nyoungs(ilay)) = &
             & Youngs_modul(ilay, 1:nyoungs(ilay)) + &
             & Youngs_modul_mat(k, 1:nyoungs(ilay)) * &
             & fraction_material(ilay, j)
           th_exp_coeff(ilay, 1:nth_exp(ilay)) = &
             & th_exp_coeff(ilay, 1:nth_exp(ilay)) + &
             & th_exp_coeff_mat(k, 1:nth_exp(ilay)) * &
             & fraction_material(ilay, j)

         end do

         if (k_tr == no_layer_mat(no_layer))  then
           if_rupture_time_creep(no_layer) = .true.
           write(aux_lun, *) 'rupture_time_creep considered for layer ', no_layer
         endif

       endif

       ! print out the bulk modulus
       if (ilay == 1)  then

         write(aux_lun, 16) Youngs_modul(ilay, 1:nyoungs(ilay)) / &
            & (2.0 * (1.0 + poisson_ratio(ilay)))
 16      format('bulk_modulus ', 20(1pe13.6, 1x))
       endif

       where (creep_const > 0.0)
         creep_const_ln = LOG(creep_const)
       end where

       write(aux_lun, 14) TRIM(layer_name(no_layer)), &
          & thickness_fr(no_layer), poisson_ratio(no_layer), &
          & surf_fracture_energy(no_layer), &
          & fracture_toughness(no_layer), creep_const(no_layer), &
          & creep_exponent(no_layer), creep_activation_energy(no_layer), &
          & Temp_creep_onset(no_layer)
 14    format('xc_layer ', a20, 20(1pe13.6, 1x))

       write(aux_lun, *) 'th3 ', no_layer, thickness_fr(1:no_layer)
       write(aux_lun, *) 

       if (no_layer_mat(no_layer) > 1)  then

         write(aux_lun, *) '  material_name = '//TRIM(layer_name(no_layer))
         write(aux_lun, *) '  Youngs_modul_mat = ', Youngs_modul(no_layer, 1:nyoungs(no_layer))
         write(aux_lun, *) '  Youngs_temp_mat = 10.0, 1000.0,'
         write(aux_lun, *) '  poisson_ratio_mat = ', poisson_ratio(no_layer)
         write(aux_lun, *) '  surf_fracture_energy_mat = ', surf_fracture_energy(no_layer)
         write(aux_lun, *) '  fracture_toughness_mat = ', fracture_toughness(no_layer)
         write(aux_lun, *) '  th_exp_coeff_mat = ', th_exp_coeff(no_layer, 1:nth_exp(no_layer))
         write(aux_lun, *) '  th_exp_temp_mat = 0.0,      50.0,    100.0,    150.0,    200.0,', &
            & '    250.0,    300.0,    350.0,    400.0,    450.0,    500.0,', &
            & '    550.0,    600.0,   650.0,  700.0,'
         write(aux_lun, *) '  rho_mat = ', rho
         write(aux_lun, *) '  cond_value_mat = ', cond_value
         write(aux_lun, *) '  cond_temp_mat = 20.0, 1300.0,'
         write(aux_lun, *) '  creep_const_mat = ', creep_const(no_layer)
         write(aux_lun, *) '  creep_activation_energy_mat = ', creep_activation_energy(no_layer)
         write(aux_lun, *) '  creep_exponent_mat = ', creep_exponent(no_layer)

         do i = 1, 11
           ! temp = 500.0+1.0*(i-1)
           ! stress = 100.0 ! MPa
           ! cr_out = 
           ! write(out_lun, ) 
         end do

       endif

    end if OXIDE_NML

    return

  END SUBROUTINE LAYER_INPUT 

END MODULE READ_INPUT_MODULE
