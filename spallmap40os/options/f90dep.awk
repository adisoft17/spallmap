################################################################################
#
#     Author: John A. Turner
#             Los Alamos National Laboratory
#             P.O. Box 1663, MS B226, Los Alamos, NM  87545
#             Group: XTM (Transport Methods)
#             Location: TA-3, Bldg. 43, Rm. D150
#             Phone: 505-665-1303
#             e-mail: turner@lanl.gov
#     $Id: f90dep.awk,v 1.1.1.1 2002/04/22 15:39:44 uxn Exp $
#
#     [g|n]awk script to generate Fortran 90 module dependencies for
#     Telluride (modified from similar script used by JTpack90).
#
#     NOTE: Either nawk or gawk is required because of the need for gsub.
#
#     USAGE:
#
#       Assuming this file is named f90dep.awk, the input file is named
#       named source.F90, the following command could be used:
#
#            [g|n]awk -f f90dep.awk  \
#                     -v outfile=source  \
#                     -v lib=libfoo.a  \
#                        source.F90
#
#       If the code in source.F90 used module types_module, contained in
#       file types_module.F90, the following would be sent to stdout:
#
#libfoo.a(source.o): libfoo.a(types_module.o)
#
#     VARIABLES:
#
#       Required:
#
#         outfile - basename of file being processed
#         lib - library in which the object files are archived
#
################################################################################

# A line to be processed:
#  1) has "use" as the first field, and
#  2) whose second field contains neither PGSLib nor JTpack.
(tolower($1) == "use" && 
 tolower($2) !~ /pgslib/ && 
 tolower($2) !~ /jtpack/ &&
 tolower($2) !~ /jt_/) {

  # The object file associated with the file we're processing, as an archive
  # member.
  outfile_obj = lib"("outfile".o)"

  # The object file associated with the dependency, as an archive member.
  module_obj = lib"("$2".o)"

  # Get rid of trailing comma (which might be there because of specifiers).
  gsub( ",", "", module_obj )

  # Write dependency to stdout.
  print outfile_obj": "module_obj
}
