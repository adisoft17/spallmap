MODULE PARAMETER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all integer and scalar parameters.
  !
  ! Author(s): Adrian S. Sabau, sabaua@ornl.gov
  !
  !=======================================================================

  implicit none
  save
  ! Public Module
  public

  ! string sizes for I/O 
  integer, parameter :: string_dim = 20, string_len = 256

  ! mf stands for molar fraction
  ! o2_air_mf - molar fraction of oxigen in air
  real, parameter ::  o2_air_mf0 = 0.21
  real, parameter ::  o2_air_mf = 0.209476

  ! n2_air_mf - molar fraction of nitrogen in air
  real, parameter ::  n2_air_mf0 = 0.79

  real, parameter ::  n2_air_mf = 0.78084

  ! Argone in air
  real, parameter :: ar_air_mf = 0.009365

  ! CO2 in air
  real, parameter :: co2_air_mf = 0.000319

  ! useful constants
  real, parameter ::    one_over_o2_mf = 1.0/o2_air_mf0
  real, parameter ::    n2_over_o2_mf = n2_air_mf/o2_air_mf0

  ! universal gas constant; jules/(mol K)
  real, parameter ::    gas_constant_jmk = 8.3144

  ! universal gas constant; k jules/(mol K)
  real, parameter ::      gas_constant_kjmk = 8.3144e-3 

  real, parameter ::   zero = 0.0

  real, parameter ::   ten_p3 = 1.0e+3

  integer, parameter :: npulse = 2000

  ! maximum number of points during the ramp down and ramp up full-low load
  integer, parameter :: mramp = 20
  ! maximum number of cycles at which results will be printed out at
  ! the ramp down full-low load 
  integer, parameter :: mramp_out = 20

  ! conversion K_Ic [MPa (micron)1/2] = fract_tough2si * K_Ic [MPa (m)1/2] 
  real, parameter ::  fract_tough2si = 1.0e+3

  ! conversion gamma [N micron/m2] = surf_energ2si * gamma [J/m2]
  real, parameter :: surf_energ2si = 1.0e+6

  ! conversion from microns to SI 
  ! ox_thickness [m] = ox_thickness [micron] * oxide_thick2si
  real, parameter :: oxide_thick2si = 1.0e-6

  ! conversion constants; from NASA to SI units
  real, parameter ::   visc_nasa2si = 1.0e-7  ! viscosity
  real, parameter ::   cond_nasa2si = 1.0e-4  ! thermal conductivity

  ! cp is given in J/(mole K)
  ! mw is given in g/mole
  ! cp/mw will be in J/(g K); it needs a conversion factor to si units
  ! real, parameter ::   cp_nasa2si = 1.0e+3 when MW was in g/mole
  real, parameter ::   cp_nasa2si = 1.0   ! MW was converted to Kg/mole

  ! conversion unit for mw to SI units
  ! do the MW conversion once and do not wory about it
  real, parameter :: mw_2si = 1.0e-3


  ! maximum number of data points

  ! maximum number of oxide layers, maximum number of oxide species; 
  ! the first species is the metal by default
  ! maximum number of times at which the most of variables are stored
  integer, parameter :: moxide_layer = 20, moxide = 50, mtime = 100

  ! mgrid = maximum number of grids within each layer
  integer, parameter :: mgrid = 100

  ! mave: maximum number of averaging quantities
  ! 1 - max
  ! 2 - min
  ! 3 - sum average
  integer, parameter  :: mave = 5

  ! maximum numbers of data for the time step
  ! sometimes, a smaller time step may be required
  integer, parameter :: mstep = 100

  ! maximum number of data sets for the oxidation rate
  integer, parameter :: mrate = 40

  integer, parameter   :: id_low = 2, id_high = 1

  ! maximum number of damage mechanisms
  integer, parameter :: mdamage = 6

  ! remove the data below this lane later on

  ! maximum number of fuel compounds, maximum number of combustion compounds
  integer, parameter :: mfuel = 10, mcomb= 20, mdiluent = 10

  ! maximum number of species read from the property database file
  integer, parameter :: mspecie = 201

  ! The following parameters are for NASA - CEA model of dealing with 
  ! chemical equilibrium reactions

  ! ntref - number of temperature points
  integer, parameter ::   mtref = 3

  ! mtrans = number of constants in NASA-CEA transport model 
  ! for thermal conductivity and viscosity
  integer, parameter ::   mtrans = 4
  ! mtint is maximum number of temperature intervals
  integer, parameter :: mtint = mtref - 1

  ! mcp - maximum number of constants in CEA-NASA model of specific heat 
  integer, parameter ::   mcp = moxide

  ! mh - maximum numbers of constans in CEA-NASA  model of enthalpy
  integer, parameter ::   mh = 8
 
  ! ms - maximum numbers of constans in CEA-NASA  model of entropy
  integer, parameter ::   ms = 8

  ! mgas - maximum numbers of gas components
  ! after combustion there will be air and gas (but in variable proportion)
  integer, parameter ::   mgas = 4

  ! mstates - maximum numbers of thermodynamic states
  ! this includes the intermediate states, isentropic or adiabatic states
  integer, parameter ::   mstates = 20

  ! maximum number of species in the oxidant
  integer, parameter ::   moxidant = 10

  ! maximum number of coolant flows from the compressor to the turbine
  integer, parameter ::   mcool = 5

  ! maximum number of mixing points
  integer, parameter ::   max_mix = 5

   ! maximum numbers of functions and unknowns for the system
   integer, parameter :: nmax = 5 * moxide_layer

   ! maximum number of compressor stages
   integer, parameter :: mcompr = 10

   ! maximum number of gas generator turbine stages
   integer, parameter :: mggt = 10

   ! maximum number of rotor expansions
   integer, parameter  :: mrotor_expansion = 10

   ! maximum number of rows for the gas generator turbine
   ! each row consists of a stator (vane) and rotor (blade)
   integer, parameter :: mrow_ggt = 2 * mggt

   ! maximum number of inlet and outlet for the gas generator turbine
   integer, parameter :: m_in_out_ggt = 4 * mggt

   ! maximum number of power generator turbine stages
   integer, parameter :: mpgt = 10

   ! maximum number of rows for the power generator turbine
   integer, parameter :: mrow_pgt = 2 * mpgt

   ! maximum number of inlet and outlet for the power generator turbine
   integer, parameter :: m_in_out_pgt = 4 * mpgt


   ! chord intervals, maximum number
   ! span intervals, maximum number; use it only at the reference, median 
   ! radius
   integer, parameter  :: mchord_cool = 10, mspan_cool = 1

   real, parameter     :: pi = 3.1415

   ! maximum number of problems to be run;
   ! when the coolant is constrained: problem 1 - small coolant, problem 2,
   ! the actual coolant
   integer, parameter  :: mcool_problem = 20

END MODULE PARAMETER_MODULE
