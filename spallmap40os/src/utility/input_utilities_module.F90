MODULE INPUT_UTILITIES_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define the utility routines that assist in the reading of all input.
  !
  !   Public Interface:
  !
  !     * call SEARCH_NML (unit, fatal, name1, name2)
  !         Search unit for a namelist name of name1 or name2.
  !   Obtained from Telluride
  !=======================================================================

  Implicit None
  Private

  ! public procedures
  public :: GET_ARGUMENTS, SEARCH_NML

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
CONTAINS

  SUBROUTINE GET_ARGUMENTS ()
    !=======================================================================
    ! Purpose(s):
    !   Get arguments (input file names) from the command line.
    !=======================================================================
    use output_module,        only: blank_line, input_file, &
                                    Output_String, prefix, tty_lun, WRITE_STRING

    implicit none

    ! Local Variables
    logical :: file_exist = .true.
    integer :: iargc
    integer :: l
    integer :: n
    logical :: stop_flag

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize the stop flag.
    stop_flag = .false.

       prefix = ' '
1      continue

       ! The prefix must be set to "blanks" the first time this routine
       ! is called to get the prefix argument from the execute line.
       ! If the prefix is set to "stop", the call to getarg is skipped
       ! and the usr is asked to input a prefix.  Thus, programmers
       ! should set prefix = 'stop' prior to the second and subsequent
       ! calls to this subroutine.

       if (prefix == ' ') then
          n = IARGC()
          if (n > 0) call GETARG (1, prefix)

          if (prefix  ==  ' ') prefix = 'stop'
       end if

       ! If the prefix is "blank" or "stop", then
       ! ask the usr to input a prefix string.
       if (prefix == ' ' .or. prefix == 'stop') then
          Output_String = blank_line
          Output_String(1) = 'Enter a prefix name for all I/O files:'
          call WRITE_STRING (Output_String, tty_lun)
          read (*, FMT='(a)') prefix
       end if

       ! Check the prefix, abort job if the prefix is "blank"
       ! or if the usr requests it.
       if (prefix == 'stop' .or. prefix == 'STOP' .or. &
            prefix == 'quit' .or. prefix == 'QUIT' .or. &
            prefix == 'end'  .or. prefix == 'END'  .or. &
            prefix == ' ') then
          stop_flag = .true.
          go to 999
       end if

       ! If the prefix has a trailing ".inp", delete it.
       n = INDEX(prefix,'.inp')
       if (n /= 0) prefix(n:) = ' '

       input_file = TRIM(prefix) // '.inp'
       inquire (FILE=TRIM(input_file) , EXIST=file_exist)

       if (.not.file_exist) then
          Output_String = blank_line
          write (Output_String, 2) TRIM(input_file)
2         format (/' Error opening file ',a,/, &
                   ' Re-enter new prefix else stop, quit or end',/)
          call WRITE_STRING (Output_String, tty_lun)
          prefix = 'stop'
          go to 1
       end if

999 continue

    If (stop_flag)  then
      write(tty_lun, *)  'STOP - input file was not read'
      stop
    endif 

    return

  END SUBROUTINE GET_ARGUMENTS

  SUBROUTINE SEARCH_NML (unit, fatal, namelist_name1, namelist_name2)
    !=======================================================================
    ! Purpose(s):
    !   Search the input for a namelist names namelist_name1 and
    !   namelist_name2. If found, backspace one record to be correctly
    !   positioned for the next namelist read. If the namelist names 
    !   are not found (from the current position), then set fatal flag
    !   to .true.  A '$' or '&' namelist name indicator is allowed.
    !=======================================================================
    implicit none

    ! Argument List
    character(LEN = *),       intent(IN)    :: namelist_name1
    character(LEN = *),       intent(IN)    :: namelist_name2
    logical, intent(INOUT) :: fatal
    integer, intent(IN)    :: unit

    ! Local Variables
    character(LEN = 32)      :: ampersand, dollar
    character(LEN = 150)     :: line, string
    integer :: c1, c2, c3, c4, l, n, names

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Number of possible namelist name variations
    names = 2

    ! Preset fatal flag to true
    fatal = .true.

    ! Exit if any of the namelist names are blank
    if (LEN_TRIM(namelist_name1) == 0 .or. &
        LEN_TRIM(namelist_name2) == 0) go to 999

    ! Read next line of input file into a character string
1   continue
    read (unit, '(a)', END = 999) line

    NAME_LOOP: do n = 1,names

       ! Copy the Namelist name
       select case (n)
          case (1)
             string = namelist_name1
          case (2)
             string = namelist_name2
       end select

       ! Get the string length
       l = LEN_TRIM(string)

       dollar    = '$'//string(:l)
       ampersand = '&'//string(:l)

       ! Search for namelist name
       c1 = INDEX(line,dollar(:l))
       c2 = c1 + l
       c3 = INDEX(line,ampersand(:l))
       c4 = c3 + l

       if (c1 == 0 .and. c3 == 0) cycle  ! No match found

       CHECK_NAME: if (c1 > 0) then

          if (line(c1:c2) == dollar(:l+1)) then
             ! Namelist found; backspace to start of namelist record
             backspace unit
             fatal = .false.
             go to 999
          else
             ! No match found
             cycle
          end if

       else if (c3 > 0) then

          if (line(c3:c4) == ampersand(:l+1)) then
             ! Namelist found; backspace to start of namelist record
             backspace unit
             fatal = .false.
             go to 999
          else
             ! No match found
             cycle
          end if

       end if CHECK_NAME

    end do NAME_LOOP

    ! Go back to read another line
    go to 1

999 continue

    return

  END SUBROUTINE SEARCH_NML

END MODULE INPUT_UTILITIES_MODULE
