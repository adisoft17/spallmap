# spallmap/GNUmakefile
#
#                         SPALLMAP Makefile
#
#             * * * * * * * * IMPORTANT * * * * * * * * *
#             *  This makefile is designed specifically *
#             *  for GNU Make Version 3.70 or higher.   *
#             *  GNU make is free software obtainable   *
#             *  from the Free Software Foundation via  *
#             *  anonymous ftp to prep.ai.mit.edu.      *
#             * * * * * * * * IMPORTANT * * * * * * * * *
#
#
# documentation is in spallmap/doc/development/build.txt
#
# ================================== targets ===================================
#
# to compile spallmap only
#    % gmake spallmap
#
# to clean up your file space
#    % gmake clean              -> removes temporary files, leaves the binary,
#                                  library, debugger and cross reference files
#    % gmake cleaner            -> clean + removes the binary, library,
#                                  debugger, cores and cross reference files
#    % gmake cleanest           -> cleaner + emacs ~ files
#
# to tar the spallmap source code for distribution
#    % gmake tar                -> tars source into ../spallmap.tar.gz
#
# To run Spallmap on one (or all) of the standard test problems:
#    % gmake run_all              -> Run all *** problems
#    % gmake syngas1     -> Run syngas1 problem
#
# To run Spallmap on one (or all) of the standard test problems
#   and diff the file.out output files with the corresponding repository versions:
#
#    add "compare=yes" to above test targets
#
# To rn the standard test suite and save the diffs
#    % gmake test               -> Run the tests, save the diffs
#
# To count the number of lines of code and comments in the entire source:
#    % gmake loc
#
# ==============================================================================

# get the spallmap version number (returned in VERSION)
include options/spallmap.tag

# build the spallmap executable 
spallmap:
	@cd src/build; $(MAKE) DEPENDS=yes spallmap

# count the lines of (spallmap) source code
loc:
	@cd src/build; $(MAKE) loc

# remove most of the leftovers
myclean:
	-@cd src/build;  $(MAKE) clean
clean:
	-@cd src/build;  $(MAKE) clean
#	-@if [ -d problems/tests ] ; then cd problems/tests; $(MAKE) clean; fi

# remove more of the leftovers
cleaner:
	-@cd src/build;  $(MAKE) cleaner
#	-@if [ -d problems/tests ] ; then cd problems/tests; $(MAKE) cleaner; fi

# remove even the emacs ~ files
cleanest distclean realclean superclean:
	@echo "scrubbing the directory structure"
	-@cd src/build;  $(MAKE) cleanest
#	-@if [ -d problems/tests ] ; then cd problems/tests; $(MAKE) cleanest; fi
	find . \( -name \*~ -o -name .\*~ -o -name core \) -exec $(RM) $(RMFLAGS) {} \;

# create a tar archive of the current tree, after cleaning
tar: cleanest
	@echo "writing tar file"
	@rm -f ../spallmap-$(VERSION) spallmap-$(VERSION).tar*
	@ln -s `pwd` ../spallmap-$(VERSION)
	@cd ..; gtar --create --dereference --exclude-from=spallmap-$(VERSION)/src/options/tar.exclude --file=- spallmap-$(VERSION) | gzip > spallmap-$(VERSION)/spallmap-$(VERSION).tar.gz
	@rm -f ../spallmap-$(VERSION)
	@echo "left distribution in `pwd`/spallmap-$(VERSION).tar.gz"

# ==============================================================================

# targets which involve running problems
# can be invoked with "compare=yes" to run diff on the output

# run_all: problem1 problem2 problem3


# ==============================================================================

# uncomment as needed
diagnostic:
	echo "version: $(VERSION)"
	cd src/build; $(MAKE) diagnostic

# ==============================================================================

