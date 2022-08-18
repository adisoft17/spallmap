MODULE CREDIT_MODULE
  !=======================================================================
  !   Define variables that describe the name, version,
  !   and credit notice for the code
  !
  ! Contains: None
  !
  ! Author: Adrian S. Sabau, sabaua@ornl.gov
  !
  !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================
  use parameter_module, only: string_len
  implicit none

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Character strings defining the code name, version, compilation
  ! date, and host name on which calculation is performed
  character(LEN = string_len), public, save :: code_name, code_version, &
                                               code_date,  host_name,   &
                                               run_date, architecture,  &
                                               libraries
  ! Copyright credit
  character(LEN = string_len), dimension(13), public, save ::             &
     Copyright =                                                          &
     (/                                                                   &
       '                                                               ', &
       'UT-BATTELLE, LLC AND THE GOVERNMENT MAKE NO REPRESENTATIONS AND', &
       'DISCLAIM ALL WARRANTIES,BOTH EXPRESSED AND IMPLIED. THERE ARE  ', &
       'NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS ', &  
       'FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE      ', &
       'WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER  ', &
       'ROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE   ', &
       'INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT     ', &
       'RESULT IN INJURY OR DAMAGE. THE USER ASSUMES RESPONSIBILITY   ', &
       'FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION,', &
       'AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING  ', &
       'OUT OF IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL       ', & 
       'OF THE SOFTWARE. Copyright 2018 UT-Battelle LLC.  All rights reserved.'  &
      /)
END MODULE CREDIT_MODULE
