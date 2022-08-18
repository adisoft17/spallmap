MODULE OUTPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define quantities and procedures for formatted output.
  !
  !   Public Interface:
  !
  !     * call INITIALIZE_IO ()
  !
  !         Open the input file and all output files.
  !      
  !     * call WRITE_INFO_STRING (raw_string, unit1, unit2, unit3, &
  !                               unit4, unit5, formatted_output)
  !
  !         Write an informational string in the standard info format.
  !      
  !     * call WRITE_STRING (string, unit1, unit2, unit3, unit4, unit5)
  !
  !         Write a string to unit1, unit2, ..., and/or unit5.
  !
  ! Contains: INITIALIZE_IO
  !           WRITE_INFO_STRING
  !           WRITE_STRING
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !          Subroutines obtained from Telluride
  !=======================================================================
  use parameter_module, only: string_dim, string_len

  implicit none

  ! Private module
  private

  ! Public Procedures
  public :: INITIALIZE_IO, WRITE_INFO_STRING, WRITE_STRING, WRITE_ONELINE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Variables for input/output files
  character(LEN = string_len), save, public :: input_file, prefix, title

  ! Variable for specie database
  ! species_thermo_file = contains the CEA-NASA input file for Cp, h, s
  ! species_id_thermo_file = contains an index file that points to the 
  !    species record in the species_thermo_file
  ! species_transport_file = contains the CEA-NASA input file for transport
  !    properties, such as thermal conductivity, viscosity
  ! species_id_transport_file = contains an index file that points to the 
  !    species record in the species_transport_file
  character(LEN = string_len), save, public :: species_thermo_file, &
        & species_id_thermo_file, &
        & species_transport_file, species_id_transport_file

  ! input files for simplified output
  character(LEN = string_len), save, public :: metal_property_file, &
     & magnetite_property_file, &
     & hematite_property_file, spinel_property_file, &
     & oxide_growth_file, damage_map_file, oxide_layers_file, &
     & solver_options_file

  ! Logical unit numbers for all output
  integer, save, public :: inp_lun, out_lun, cea_lun, &
        & aux_lun, adv_lun, tty_lun, can_lun, &
        & spe_th_lun, spe_tr_lun, spe_id_th_lun, spe_id_tr_lun, &
        & strn_lun, enrg_lun, strs_lun, f2l_lun, tmp_lun, gen_lun, &
        & f2l_strn_lun, enrg_dist_lun, enrg_dil_lun, strn_gen_lun, &
        & f2l_pl_lun

  ! Blank line character string default
  character(LEN = string_len), save, public :: blank_line = 'Ex Libris'
  
  ! Character string to be passed into the output routines
  character(LEN = string_len), dimension(string_dim), save, public :: Output_String

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

Contains

  SUBROUTINE INITIALIZE_IO ()
    !=======================================================================
    ! Purpose(s):
    !   Open I/O files and initialize various I/O-related quantities.
    !=======================================================================

    implicit none

    ! Local Variables
    integer :: l

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Assign logical unit numbers for all I/O.
    inp_lun =  1
    out_lun =  2
    ! err_lun =  3
    aux_lun =  4
    ! stdin  = 5
    tty_lun =  6
    ! stderr = 7
    adv_lun =  8  ! advanced input file; actual input file trc
    ! int_lun =  9
    ! gra_lun = 10
    strn_lun = 7  ! strain output file
    enrg_lun = 9  ! energy ouput file
    strs_lun = 10  ! stress output file
    f2l_lun = 3   ! outage output file
    tmp_lun = 18  ! temporary file open for simplified input
    gen_lun = 19  ! general output information data file 
    f2l_strn_lun = 20   ! outage strain output file
    enrg_dist_lun = 21  ! distortion elastic energy - shear 
    enrg_dil_lun = 22   ! dilatational elastic energy - tensile/compression
    strn_gen_lun = 23
    f2l_pl_lun = 25

    ! input file for CEA - NASA combustion
    cea_lun = 11
    ! input file for cantera - combustion
    can_lun = 12  
   
    ! input file for database
    spe_th_lun = 13
    spe_id_th_lun = 14
    spe_tr_lun = 15
    spe_id_tr_lun = 16

       ! Open the input file
       l = LEN_TRIM(input_file)
       OPEN (UNIT = inp_lun, FILE = input_file(:l), STATUS = 'old', &
             POSITION = 'rewind')

       ! Open the output edit file
       input_file = Trim(prefix) // '.out'
       l = LEN_TRIM(input_file)
       OPEN (UNIT = out_lun, FILE = input_file(:l), STATUS = 'unknown', &
             POSITION = 'rewind')

       ! Open the output auxiliary file
       input_file = Trim(prefix) // '.aux'
       l = LEN_TRIM(input_file)
       OPEN (UNIT = aux_lun, FILE = input_file(:l), STATUS = 'unknown', &
             POSITION = 'rewind')

    return

  END SUBROUTINE INITIALIZE_IO

  SUBROUTINE WRITE_INFO_STRING (raw_string, unit1, unit2, unit3, unit4, &
                                unit5, formatted_output)
    !======================================================================
    ! Purpose(s):
    !   Write an informational string in the Telluride info format.
    !   Optionally, return a formatted string with the info format.
    !======================================================================
    use parameter_module, only: string_dim

    implicit none

    ! Argument List
    character (LEN=*),                         Intent(IN)  :: raw_string
    character (LEN=*), Optional, Dimension(:), Intent(OUT) :: formatted_output
    integer, Optional, Intent(IN) :: unit1
    integer, Optional, Intent(IN) :: unit2
    integer, Optional, Intent(IN) :: unit3
    integer, Optional, Intent(IN) :: unit4
    integer, Optional, Intent(IN) :: unit5

    ! Local Variables
    character (LEN=256), dimension(string_dim) :: formatted_string
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    formatted_string = blank_line
    Write(formatted_string, 10) Trim(raw_string)
10  Format(/,9x, A,/)

    if (PRESENT(formatted_output)) formatted_output = formatted_string

    if (PRESENT(unit5)) Then
       call WRITE_STRING (formatted_string, unit1, unit2, unit3, unit4, unit5)
    else if (PRESENT(unit4)) Then
       call WRITE_STRING (formatted_string, unit1, unit2, unit3, unit4)
    else if (PRESENT(unit3)) Then
       call WRITE_STRING (formatted_string, unit1, unit2, unit3)
    else if (PRESENT(unit2)) Then
       call WRITE_STRING (formatted_string, unit1, unit2)
    else if (PRESENT(unit1)) Then
       call WRITE_STRING (formatted_string, unit1)
    end if
    
    return

  END SUBROUTINE WRITE_INFO_STRING

  SUBROUTINE WRITE_STRING (String, unit1, unit2, unit3, unit4, unit5)
    !=======================================================================
    ! Purpose(s):
    !   Write a string to unit1, unit2, unit3, unit4, and unit5.
    !=======================================================================

    implicit none

    ! Argument List
    character(LEN = *),   dimension(:), intent(IN) :: String
    integer, optional, intent(IN) :: unit1
    integer, optional, intent(IN) :: unit2
    integer, optional, intent(IN) :: unit3
    integer, optional, intent(IN) :: unit4
    integer, optional, intent(IN) :: unit5
    
    ! Local Variables
    character(LEN=256)       :: Out_String
    logical :: write1
    logical :: write2
    logical :: write3
    logical :: write4
    logical :: write5
    logical :: write_tty
    integer :: line, str_end, str_size

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get the string size
    str_size = SIZE(String)

    ! Loop in reverse order over the str_size, looking for the first
    ! string that is none-null (not empty and not equal to blank_line)
    if (str_size >= 3) then
       do str_end = str_size,1,-1
          if (TRIM(String(str_end)) /= TRIM(blank_line) .and. &
              LEN_TRIM(String(str_end)) /= 0) exit
       end do
    else
       str_end = str_size
    end if

    ! Check to see which logical unit numbers are present, and
    ! which units, if any, are the tty unit (stdout)
    write1 = PRESENT(unit1)
    write2 = PRESENT(unit2)
    write3 = PRESENT(unit3)
    write4 = PRESENT(unit4)
    write5 = PRESENT(unit5)
    write_tty = .False.

    ! Set the unit1 write flag; see if it is the tty unit number
    if (write1) then
       write_tty = unit1 == tty_lun
       write1 = unit1 /= tty_lun
    end if
    ! Set the unit2 write flag; see if it is the tty unit number
    if (write2) then
       if (.not.write_tty) write_tty = unit2 == tty_lun
       write2 = unit2 /= tty_lun
    end if
    ! Set the unit3 write flag; see if it is the tty unit number
    if (write3) then
       if (.not.write_tty) write_tty = unit3 == tty_lun
       write3 = unit3 /= tty_lun
    end if
    ! Set the unit4 write flag; see if it is the tty unit number
    if (write4) then
       if (.not.write_tty) write_tty = unit4 == tty_lun
       write4 = unit4 /= tty_lun
    end if
    ! Set the unit5 write flag; see if it is the tty unit number
    if (write5) then
       if (.not.write_tty) write_tty = unit1 == tty_lun
       write5 = unit5 /= tty_lun
    end if

    ! Loop over the lines in the string and write them out
    STRING_ELEMENTS: do line = 1, str_end

       if (TRIM(String(line)) == TRIM(blank_line)) then
          Out_String = ''
       else
          Out_String = TRIM(String(line))
       end if

          if (write_tty) write (*, 1)  TRIM(Out_String)
          if (write1) write (unit1, 1) TRIM(Out_String)
          if (write2) write (unit2, 1) TRIM(Out_String)
          if (write3) write (unit3, 1) TRIM(Out_String)
          if (write4) write (unit4, 1) TRIM(Out_String)
          if (write5) write (unit5, 1) TRIM(Out_String)
1         format (a)

    end do STRING_ELEMENTS

    return

  END SUBROUTINE WRITE_STRING

  SUBROUTINE WRITE_ONELINE (string)
    !=======================================================================
    ! Purpose(s):
    !   Write one line to tty_lun and out_lun
    !=======================================================================
    implicit none

    ! Argument List
    character (LEN=*) :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Output_String    = blank_line
    Output_String(1) = string

    call WRITE_STRING (Output_String, tty_lun, out_lun)

    return

  END SUBROUTINE WRITE_ONELINE

END MODULE OUTPUT_MODULE
