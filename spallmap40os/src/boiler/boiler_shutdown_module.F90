MODULE BOILER_SHUTDOWN_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Describe the boiler shut down, temperatures, pressures, flow rates
    !   htc 
    !  sabaua@ornl.gov read_input_module.F90
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: Temp_Shut_Down

 CONTAINS

  REAL FUNCTION Temp_Shut_Down(Function_Id, Temp0, t)

    REAL, Intent(IN)    :: t      ! time
    REAL, Intent(IN)    :: Temp0  ! Temperature at t=0
    Integer, Intent(IN) :: Function_Id

    ! Local Variables
    real        ::  m1 = 2.7880322e+2, m2=1.63327619e-2, m3=3.13957295e+1, &
          & m0 = 50.0
    real, save  ::  m1e, m3e, m0a
    logical, save  :: if_first_log = .true.

    SELECT CASE (Function_Id)

    CASE (1)

      ! superheater
! y = 50 + m1*exp(-m2*x)+ m3*exp(-10.0*m2*x)
!        Value   Error
!m1      2.7880322e+2    1.8960784e+0
!m2      1.63327619e-2   1.5279008e-4
!m3      3.13957295e+1   3.9594671e+0

      if (if_first_log)  then

        if (ABS(m0+m1+m3 - Temp0) > 0.5)  then
          m0a = Temp0 - m1 - m3
          write(2, *) 'exp_decay adjust ', m0, ' to ', m0a
        else
          m0a = m0
        endif

        if (m1 > 0) m1e = LOG(m1)
        if (m3 > 0) m3e = LOG(m3)

        if_first_log = .false.

      endif

      Temp_Shut_Down = m0a + exp(m1e - m2 * t)+ &
        & exp(m3e - 10.0 * m2 * t)

    CASE DEFAULT

       write(6, *) 'ERROR:Temp_Shut_Down not ok'
       STOP

    END SELECT

    RETURN

  END FUNCTION Temp_Shut_Down

END MODULE BOILER_SHUTDOWN_MODULE
