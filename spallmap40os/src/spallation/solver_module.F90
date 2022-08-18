MODULE SOLVER_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   solver interface  
    !
    !  Author: Adrian S. Sabau, sabaua@ornl.gov
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: NEWTON_RAPHSON, SOLVER_SIMPLE

 CONTAINS

 SUBROUTINE NEWTON_RAPHSON(iter, convergent)
    use parameter_module, only: nmax
    use solver_data_module, only       : neq_solver, fval, var, &
        & delta_var, dTemp_eps
    use output_module,        only: aux_lun, blank_line,   &
                                  & out_lun, Output_String, prefix, &
                                  & tty_lun, WRITE_STRING

    implicit none

    ! arguments
    integer,    intent(IN)    :: iter
    logical,    intent(INOUT) :: convergent

    ! local variables
    logical, dimension(nmax)                 :: Mask
    integer            :: i
    integer, dimension(1)           :: func_id_max, delta_id_max
    real               :: func_ave, func_max, delta_ave, delta_max, &
                & dvar_over_var
    
    convergent = .false.

    Mask = .false.
    do i = 1, neq_solver
      Mask(i) = .true.
    end do

    ! obtain the average and max of the absolute value for the function
    func_ave = SUM(abs(fval), Mask)
    func_max = MAXVAL(abs(fval), Mask)
    func_id_max = MAXLOC(abs(fval), Mask)

    if (iter > 1)  then
      dvar_over_var = MAXVAL(abs(delta_var / var), Mask .and. abs(var) > 1.0e-12)
    else
      dvar_over_var = 1.0
    endif

    ! solver returns delta_var; the increment
    call SOLVER()

    ! update solution
    var = var + delta_var

    ! check for convergence
    delta_ave = SUM(abs(delta_var), Mask)
    delta_max = MAXVAL(abs(delta_var), Mask)
    delta_id_max = MAXLOC(abs(fval), Mask)
    
    if (func_max < dTemp_eps .or. dvar_over_var < 1.0e-4)  convergent = .true.

      Output_String = blank_line
      write (Output_String, 61) iter, func_id_max, delta_id_max, &
       & func_max, delta_max, dvar_over_var
 61   format('it= ', i3, ' fid= ', i3, ' did= ', i3, &
      & ' fmax= ', 1pe13.6, ' dx_max=', 1pe13.6, ' max(dx/x)=', 1pe13.6)
      call WRITE_STRING (Output_String, tty_lun, out_lun)

    return

 END SUBROUTINE NEWTON_RAPHSON

 SUBROUTINE SOLVER_SIMPLE (neq_solver, matrix, solution)

    use solver_data_module, only : nrhs

    implicit none

    external                           DGETRS
    external                           DGETRF

    ! arguments
    integer, intent(IN)  :: neq_solver
    real(kind = 8), dimension(neq_solver, neq_solver), intent(INOUT) :: matrix
    real(kind = 8), dimension(neq_solver), intent(INOUT)      :: solution

    ! local variables
    integer, dimension(neq_solver)                 :: ipiv
    integer                                        :: solver_error

    integer            :: i, j, prepare_in = 1, prepare_out = 2

    do i=1, neq_solver
      ! write(2, 3) (matrix(i, j), j=1, neq_solver)
 3    format(50(1pe9.2, 1x))
    end do
!
!  Factor the matrix.
!
    call DGETRF (neq_solver, neq_solver, matrix, neq_solver, IPIV, solver_error)

  if ( solver_error /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  DGETRF returned INFO = ', solver_error
    write ( *, '(a)' ) '  The matrix is numerically singular.'
    stop
  end if
!
!  Solve the linear system; 
! 
    call DGETRS('N', neq_solver, NRHS, matrix, neq_solver, IPIV, solution, &
      & neq_solver, solver_error)

  if ( solver_error /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGETRS; Solution procedure failed!'
    write ( *, '(a,i6)' ) '  INFO = ', solver_error
    stop
  end if

    return

 END SUBROUTINE SOLVER_SIMPLE

 SUBROUTINE SOLVER ()
    use solver_data_module, only       : neq_solver, fval, jacf, var, nrhs
    implicit none

    external                           DGETRS
    external                           DGETRF

    ! arguments

    ! local variables
    real(kind = 8), dimension(neq_solver, neq_solver)   :: matrix
    real(kind = 8), dimension(neq_solver)          :: solution
    integer, dimension(neq_solver)                 :: ipiv
    integer                                        :: solver_error

    integer            :: i, j, prepare_in = 1, prepare_out = 2

    ! allocate matrix and vector for the right size
    call SET_MATRIX(prepare_in, neq_solver, matrix, solution)

    do i=1, neq_solver
      ! write(2, 3) (matrix(i, j), j=1, neq_solver)
 3    format(50(1pe9.2, 1x))
    end do
!
!  Factor the matrix.
!
    call DGETRF (neq_solver, neq_solver, matrix, neq_solver, IPIV, solver_error)

  if ( solver_error /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  DGETRF returned INFO = ', solver_error
    write ( *, '(a)' ) '  The matrix is numerically singular.'
    stop
  end if
!
!  Solve the linear system; jacf * dvar = - function
! 
    call DGETRS('N', neq_solver, NRHS, matrix, neq_solver, IPIV, solution, &
      & neq_solver, solver_error)

  if ( solver_error /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGETRS; Solution procedure failed!'
    write ( *, '(a,i6)' ) '  INFO = ', solver_error
    stop
  end if

    ! update the actual vector solution
    call SET_MATRIX(prepare_out, neq_solver, matrix, solution)

    return

 END SUBROUTINE SOLVER

 SUBROUTINE SET_MATRIX(prepare, nsize, matrix, solution)

    use solver_data_module, only       : neq_solver, fval, jacf, delta_var
    ! arguments
    integer,    intent(IN)            :: prepare
    integer,    intent(IN)            :: nsize
    real(kind = 8), dimension(nsize, nsize),  intent(INOUT)  :: matrix
    real(kind = 8), dimension(nsize),  intent(INOUT)  :: solution

    ! local variables
    integer   :: i, j
    ! allocate matrix and vector for the right size

    if (prepare == 1)  then
      do i = 1, neq_solver
        solution(i) = -fval(i)
        do j = 1, neq_solver
          matrix(i, j) = jacf(i, j)
        end do
      end do
    else if (prepare == 2)  then
      do i = 1, neq_solver
        delta_var(i) = solution(i)
      end do
    else 
      write(6, *) 'SET_MATRIX prepare has wrong option'
      stop
    endif

    return

 END SUBROUTINE SET_MATRIX

 END MODULE SOLVER_MODULE
