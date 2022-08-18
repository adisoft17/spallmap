  Program SPALLMAP
  !===============================================================================
  !                                  SPALLMAP
  !                    Spallation MAP  
  !	                    Computer Program
  !
  !        
  !  SpallMap solves for the stress-strain equations in an axisymmetric geometry, 
  !      tracking the stress/strain evolution during boiler operation including 
  !      outages at one-location along a boiler tube and compares it with scale 
  !      damage criteria represented by Armitt diagram.
  !        
  !       UT-BATTELLE, LLC AND THE GOVERNMENT MAKE NO REPRESENTATIONS AND 
  !       DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.  THERE ARE NO 
  !       EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A 
  !       PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE 
  !       ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT 
  !       THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE 
  !       OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.  THE USER ASSUMES 
  !       RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, 
  !       CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR 
  !       ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE 
  !       SOFTWARE.
  !
  !  Developed by: Adrian S. Sabau, sabaua@ornl.gov and Ian G. Wright 
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !
  !================================================================================
  use read_input_module,         only: READ_INPUT
  use output_module,             only: INITIALIZE_IO, out_lun
  use input_utilities_module,    only: GET_ARGUMENTS
  use solver_module,             only: NEWTON_RAPHSON
  use solver_data_module,        only: max_iterat_solver, no_problem, &
             & if_average_oxide 
  use solution_data_module,      only: no_oxide_layers
  use mesh_module,               only: SOLVER_STRESS_VARIABLES
  use oxide_data_module,         only: no_oxide, no_layer
  use oxide_growth_module,       only: OXIDE_GROWTH, OXIDE_GROWTH_LAYER
  use oxide_module,              only: CHECK_CYCLE_STRAIN, &
             & CYCLE_STRAIN_CYL_ALL
  use property_module,           only: PROPERTY_SCALE_AVE
  use spallation_module,         only: CHECK_STRAIN, CHECK_MAP_EPRI
  use creep_module,              only: INIT_CREEP

  implicit none
  ! Local variables
  integer   :: iter, nproblem
  logical   :: convergent

  ! read and parse the command line
  call GET_ARGUMENTS ()

  ! initalized the I/O 
  call INITIALIZE_IO ()

  ! get the credit
  call CODE_INIT ()
  
  ! read input
  call READ_INPUT ()
  write(6, *)  'DONE input'

  call PROPERTY_SCALE_AVE()

  ! store the number of layers
  if (if_average_oxide)  then
    no_oxide_layers = 2 
  else
    no_oxide_layers = no_layer  ! check this
  endif

  ! initialize creep data; 
  ! even for cases when there is not creep some data must be set
  call INIT_CREEP(no_oxide_layers)

  call SOLVER_STRESS_VARIABLES(no_oxide_layers)

  ! get oxide thickness evolution
  call OXIDE_GROWTH_LAYER ()

  ! check the strain evolution
  ! call CHECK_STRAIN()   ! plane geometry
  ! CHECK_CYCLE_STRAIN_CYL only 2-3 and 3-4 are considered; exact without creep
  ! call CHECK_CYCLE_STRAIN_CYL()  ! circular temp profile
  call CYCLE_STRAIN_CYL_ALL()  ! circular temp profile 1-2 and 4-1 included

  ! check the strain evolution
  call CHECK_CYCLE_STRAIN()

  ! 2009 moved subroutine here needed for average properties
  ! call PROPERTY_SCALE_AVE()

  ! check the original EPRI map
  call CHECK_MAP_EPRI(.true.)

  if (.not. if_average_oxide)  then
    ! check the original EPRI map for each layer
    call CHECK_MAP_EPRI(if_average_oxide)
  endif


  LOOP_SOLVER: do iter= 1, max_iterat_solver

    convergent = .true.

    ! call NEWTON_RAPHSON(iter, convergent)

    if (convergent)  then

      write(out_lun, *) 'problem= ', nproblem, 'convergence at iter = ', iter
      exit LOOP_SOLVER

    endif

  end do LOOP_SOLVER

  write(6,*) 'end of program'
  stop

  End Program SPALLMAP

  SUBROUTINE CODE_INIT ()
   !=======================================================================
  ! Purpose:
  !
  !   Determine and print the problem information
  !
  !=======================================================================
  use credit_module,          only: architecture, code_date, code_name, &
                                  code_version, host_name, libraries, &
                                  run_date, Copyright
  use output_module,        only: aux_lun, blank_line, inp_lun, &
                                  input_file,            &
                                  out_lun, Output_String, prefix, &
                                  title, adv_lun, tty_lun,        &
                                  WRITE_STRING
  use utilities_module,     only: TIMESTAMP

  implicit none
  ! Local Variables
  integer  ::  i, k

  code_name = 'spall'

  ! Get the code version.
  code_version = '3.0.0'

  ! Get the build date.
  code_date = 'Thu Mar 19 11:31:23 EDT 2009'

  ! Get the architecture.
  architecture = 'i686'

  ! Get the host name.
  host_name = 'cast'

  ! Get the libraries.
  libraries = '-'
  libraries = TRIM(ADJUSTL(libraries)) // ', ' // '-'
  ! Write out code/problem info (including Copyright notice).
  Output_String = blank_line
  write (Output_String, '(/,36x,a)') TRIM(code_name)
  ! call WRITE_STRING (Output_String, tty_lun, out_lun, err_lun, aux_lun, trc_lun)
  call WRITE_STRING (Output_String, tty_lun, out_lun, aux_lun)
  Output_String = blank_line
  do i = 1,SIZE(Copyright)
     write (Output_String(i), '(9x,a)') TRIM(Copyright(i))
  end do
  call WRITE_STRING (Output_String, tty_lun, out_lun, aux_lun)
  Output_String = blank_line
  write (Output_String, 3) TRIM(code_version), TRIM(architecture), TRIM(libraries)
3 format(/,1x,80('='),//,32x,'Problem/Code Specs',//, &
         ' Version: ',a,//,' Architecture: ',a,//,' Libraries: ',a)
  call WRITE_STRING (Output_String, tty_lun, out_lun, aux_lun)

  ! Print compilation date, run date, and problem title.

  read (inp_lun, '(a)') title

  call TIMESTAMP (run_date)
  k = LEN_TRIM(code_date) - 9
  Output_String = blank_line
  write (Output_String, 4) code_date(5:k), run_date(5:22), TRIM(title)
4 format(/,' Build Date/Time: ',a,//,' Execution Date/Time: ',a, &
       //,' Problem: ',a,//,1x,80('='))
  call WRITE_STRING (Output_String, tty_lun, out_lun, aux_lun)

  ! Print problem initialization header.
  input_file = TRIM(prefix) // '.inp'
  Output_String = blank_line
  write (Output_String, 5) TRIM(input_file)
5 format(/,35x,'Input Phase',//,' Parsing input file ',a)
  call WRITE_STRING (Output_String, tty_lun, out_lun)

  return
END SUBROUTINE CODE_INIT
