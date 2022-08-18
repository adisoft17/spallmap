MODULE UTILITIES_MODULE
  !======================================================================
  ! Purpose(s):
  !
  !   Define general utility routines which are used throughout Telluride.
  !
  !   Public Interface(s):
  !
  !     * call STRING_COMPARE (string1, string2, string_match)
  !
  !       Performs a string comparison.
  !
  !     * call TIMESTAMP (date_time)
  !
  !       Returns the date and time in string date_time.
  !
  ! Contains: STRING_COMPARE
  !           TIMESTAMP
  !
  ! Author(s): Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Anand V. Reddy, Cateripillar (reddy_anand_v@cat.com)
  ! Subroutine obtained from Telluride
  !=======================================================================
  implicit none

  ! Private Module
  private

  ! Public Subroutines
  public :: STRING_COMPARE, TIMESTAMP

  ! File Version

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE STRING_COMPARE (string1, string2, strings_match)
    !=======================================================================
    ! Purpose(s):
    !
    !   Performs a case insensitive comparison of two character strings,
    !   and returns true or false for the comparison result.
    !
    !=======================================================================

    implicit none

    ! Argument List
    character(LEN = *)       :: string1, string2
    logical :: strings_match

    ! Local Variables
    character(len = 26) :: Upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', &
                           Lower_case = 'abcdefghijklmnopqrstuvwxyz'

    integer :: i, j, letter_no1, letter_no2, string_length

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize string match flag 
    strings_match = .true.

    ! Compare character strings
    string_length = LEN_TRIM(string1)
    if (string_length /= LEN_TRIM(string2)) then

       strings_match = .false.

    else

       ! Loop over the string length
       STRING_LOOP: do i = 1,string_length
          letter_no1 = 0
          letter_no2 = 0
          ALPHABET: do j = 1,26
             if (string1(i:i) == Upper_case(j:j) .or. &
                 string1(i:i) == Lower_case(j:j)) letter_no1 = j
             if (string2(i:i) == Upper_case(j:j) .or. &
                 string2(i:i) == Lower_case(j:j)) letter_no2 = j
          end do ALPHABET
          if (letter_no1 /= letter_no2) strings_match = .false.
          if (letter_no1 == 0) then
             if (string1(i:i) /= string2(i:i)) strings_match = .false.
          end if
          if (.not. strings_match) exit STRING_LOOP
       end do STRING_LOOP

    end if

    return

  END SUBROUTINE STRING_COMPARE

  SUBROUTINE TIMESTAMP (date_time)
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutine to build a 26-character string containing current
    !   date and time using DATE_AND_TIME intrinsic.
    !
    !                                    12345678901234567890123456
    !   String returned is of the form:  Fri 20 Aug 93 09:33:35.500
    !
    !   Initial routine obtained from ftp.ora.com:/pub/nutshell/fortran90
    !   as part of the examples for the O'Reilly book "Migrating to
    !   Fortran 90" by James F. Kerrigan.
    !
    !=======================================================================

    implicit none

    ! Argument List
    character(LEN = 26), intent(OUT) :: date_time

    ! Local Variables

    character(LEN = 3), dimension(0:6) :: days
    character(LEN = 3), dimension(12)  :: months

    integer, dimension(8) :: elements
    integer               :: m, y, w

    ! Define Days and Months
    data days   / 'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat' /  
    data months / 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'  /

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! System date and time
    call DATE_AND_TIME (VALUES = elements)

    ! Format date and time.
    INVALID: if (elements(1) /= -HUGE(0)) then
       y = elements(1)
       m = elements(2)

       if (m < 3) then
          m = m + 12
          y = y -  1
       end if

       w = MOD((elements(3) + (13*m-27)/5 + y + y/4 - y/100 + y/400 ), 7)  
       CENTURY: if (elements(1) < 2000 ) then
          elements(1) = elements(1) - 1900
       else
          elements(1) = elements(1) - 2000
       end if CENTURY

       write (date_time,10) days(w), elements(3), months (elements(2)), &
                            elements(1), elements(5), elements(6), &
                            elements(7), elements(8)
10     format (a3, 1x, i2.2, 1x, a3, 1x, i2.2, 1x,   &
            i2.2, ':', i2.2, ':', i2.2, '.', i3.3 ) 
    else INVALID
       date_time = ' '
    end if INVALID

    return

  END SUBROUTINE TIMESTAMP

END MODULE UTILITIES_MODULE
