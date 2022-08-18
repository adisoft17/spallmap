#!/opt/bin/perl
#
# Author: John A. Turner
#         Los Alamos National Laboratory
#         Phone: 505-665-1303
#         e-mail: turner@lanl.gov
#
# Script to convert geom_module, mesh_module, and vertex_module use
# statements into single use mesh_module statement and eliminate
# Mesh, Cell, and Vertex from arg lists.
#
# Usage:
#
#  Assuming this file is named mesh_hack.pl, cd to a directory
#  containing files to be converted, then:
#
#  % mkdir new
#  % mesh_hack.pl
#
#  New versions will be in subdirectory new.
#
# Assumptions (some possible because of Telluride coding conventions):
#
#  o every routine has an implicit none
#  o current working directory contains all files under consideration
#  o source files have suffix .F90
#  o only one module per use statement
#  o lowercase used for module names
#  o declarations for Mesh, Cell, and Vertex are always single lines
#  o lots of other stuff I eventually lost track of
#
# Notes:
#
#  o noticed inconsistent capitalization of "use" (sometimes "use",
#    sometimes "Use")
#  o noticed inconsistent capitalization of "only" (sometimes "only",
#    sometimes "Only")
#  o order of use statements (supposedly alphabetical) is only maintained
#    if a "use" statement for mesh_module already exists in a routine
#    needing something from mesh_module.  otherwise the "use" statement for
#    mesh_module appears immediately before the "implicit none"
#  o "only" in "use" for mesh_module isn't lined up with others
#  o simply eliminate MESH_CONNECTIVITY, VERTEX_DATA, and CELL_GEOMETRY
#    from "use" statement for mesh_module, since they are seldom used, and
#    when they are the compiler will complain and it can be fixed by hand

# Set output field separator.
$, = '';

# Set output record separator.
$\ = "\n";

# Create array of source files.
@source = glob("*.F90");

# Loop through all source files.
foreach $file (@source) {

  # Open input file.
  open(INFILE, $file);

  # Initialize hashes.
  %mesh = ();    # entities needed from mesh_module

  # First read through file and create hashes containing info on entities needed
  # from mesh_module.  Keys are routine names.
  while (<INFILE>) {
    chop;

    if (/^\s*(function|Function|FUNCTION|subroutine|Subroutine|SUBROUTINE)\s+(\w+)\s*/) {
	    
      # Initialize key and arrays for this routine.
      $routine = $2;
      
    } elsif (/^\s*(use|Use)\s+vertex_module/) {

      # If there's a use statement for vertex module, it means it's providing at
      # least VERTEX_DATA, which implies that Vertex was in the arg list.  But
      # now we provide Vertex via mesh_module, so add it to the list.
      $mesh{$routine} .= "Vertex ";
      
    } elsif (/^\s*(use|Use)\s+geom_module/) {

      # If there's a use statement for geom module, it means it's providing at
      # least CELL_GEOMETRY, which implies that Cell was in the arg list.  But
      # now we provide Cell via mesh_module, so add it to the list.
      $mesh{$routine} .= "Cell ";

    } elsif (/^\s*(use|Use)\s+mesh_module.*:\s+(.*)&?$/) {

      # Use statement for mesh_module may be continued over several lines.
      @entities = split(/,/,$2);
      foreach (@entities) {

	s/MESH_CONNECTIVITY/Mesh/;
	s/&//;  # Not sure why this is needed, since parens excluded &
	$mesh{$routine} .= $_ . " ";
      }
      if (/&$/) {

	while (<INFILE>) {

	  # As long as lines end in a continuation character, continue to
	  # add the comma-delimited entities to the list.
	  if (/(.*)&$/) {

	    @entities = split(/,/,$1);
	    foreach (@entities) {

	      s/MESH_CONNECTIVITY/Mesh/;
	      $mesh{$routine} .= $_ . " ";
	    }
	  } else {

	    # Found a line that didn't end in a continuation character, so
	    # pop out.  Not that this throws away a line, but we'll live with
	    # that since it's unlikely a use statement for vertex_module will
	    # ever directly follow a use statement for mesh_module that extended
	    # over several lines.  If it does, the compiler will save us.
	    last;
	  }
	}
      }
    }
  }

  # Fix list of entities needed from mesh_module, and echo some stuff to STDOUT.
  foreach $routine (keys %mesh) {
    $mesh{$routine} =~ tr/ /,/s;
    $mesh{$routine} =~ s/,$//;
    $mesh{$routine} =~ s/,/, /g;
    print STDOUT $routine;
    print STDOUT "   " . $mesh{$routine};
  }
  
  # Close, then reopen input file.
  close(INFILE);
  open(INFILE, $file);
  
  # Open output file.
  open(OUTFILE, ">new/$file");
  select(OUTFILE);

  # Read through file again and create new version.
  while (<INFILE>) {
    chop;

    if (/^\s*(call|Call|CALL)/) {

      # Eliminate Mesh, Cell, and Vertex from arg lists.
      # NOTE: this doesn't deal with continued lines.

      s/(\(|\s|,)Mesh(,|\))\s*/\1/;
      s/(\(|\s|,)Cell(,|\))\s*/\1/;
      s/(\(|\s|,)Vertex(,|\))\s*/\1/;

      # Check for bare ( or ,
      s/\($/\(\)/;
      s/,\s*$/\)/;

      print;

    } elsif (/^\s*(subroutine|Subroutine|SUBROUTINE)\s+(\w+)\s*/) {

      # Save routine name since it is used as a key into the hash.
      $routine = $2;

      # Determine whether a use statement for mesh_module will be
      # needed for this routine.
      $mesh_needed = 0;
      if ($mesh{$routine} ne "") {
	$mesh_needed = 1;
      }

      # Eliminate Mesh, Cell, and Vertex from arg lists.
      # NOTE: this doesn't deal with continued lines.

      s/(\(|\s|,)Mesh(,|\))\s*/\1/;
      s/(\(|\s|,)Cell(,|\))\s*/\1/;
      s/(\(|\s|,)Vertex(,|\))\s*/\1/;

      # Check for bare ( or ,
      s/\($/\(\)/;
      s/,\s*$/\)/;

      print;

    } elsif (/^\s*(use|Use)\s+(geom|vertex)_module/) {

      # Nothing.

    } elsif (/^\s*(use|Use)\s+mesh_module/) {

      # Print the use statement for mesh_module.
      print "    use mesh_module, only: " . $mesh{$routine};

      # Reset flag.
      $mesh_needed = 0;

    } elsif (/^\s*implicit/) {

      if ($mesh_needed == 1) {

	print "    use mesh_module, only: " . $mesh{$routine} . "\n";
      }
      print;

    } elsif (/(type|Type)\s*\((CELL_GEOMETRY|MESH_CONNECTIVITY|VERTEX_DATA)\).*(Cell|Mesh|Vertex)$/) {

      # Nothing.

    } else {
      
      # Just echo line.

      print;
    }
  }
  
  # Close input file.
  close(INFILE);
}

# Close output file.
close(OUTFILE);

