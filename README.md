# spallmap
The SPALLation MAP (SPALLMAP) code contains modules developed for oxide growth and oxide exfoliation.  Copyright 2018 UT-Battelle LLC.  All rights reserved.

SpallMap solves for the stress-strain equations in an axisymmetric geometry, tracking the stress/strain evolution during boiler operation including outages at one-location along a boiler tube and compares it with scale damage criteria represented by Armitt diagram. 

Problem Statement: Boiler tubes in steam power plants experience exfoliation of oxide grown on the inner side of the tubes. In extreme cases, the exfoliation cause significant tube blockages that lead to forced power plant outages. It is thus desired to predict through modeling the propensity of exfoliation events in order to inform power plant operators of possible tube blockages.

Short summary of software's functionality: The SPALLMAP code contains modules developed for oxide growth, stress analysis, and classical fracture mechanics correlations by taking into account the following phenomena and features, (a) Non-uniform thermal expansion coefficient of oxides and metal substrates, (b) Plant operation schedule with periodic alternate full-load and partial-load regimes,(c) axisymmetric formulation for cylindrical tubes, (d) Multiple oxide layers, (e) oxide-growth induced stresses, and (f) damage criteria from classical fracture mechanics. The computer program is written in FORTRAN90. Its modular structure was sought for allowing the best flexibility in updating the program by implementing new constitutive equations due to availability of new material property data and/or new physical phenomena.

A preliminary software manual is provided in the spallmap_manual17.pdf file.  The manual include instructions for installation, a description of main parameters used in the input file and a description of the output, including some instructions (line commands in UNIX and UNIX scripts) on how to extrat the data from the output file.
