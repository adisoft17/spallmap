# -*-Mode: Makefile-*-
#
#   Code name and version definitions used in the INTERFACE_POLYGONS GNUmakefile
# 
# $Id: interface_polygons.tag,v 1.1.1.1 2002/04/22 15:39:44 uxn Exp $

CODE    = interface_polygons
VERSION = 0.1

# define LIB name to be the same as the CODE name
LIB = $(CODE)

# define symbols to be #ifdef'ed into the code
BUILD_DATE   = $(shell date)
ARCHITECTURE = $(shell uname -a)
HOST_NAME    = $(shell uname -n)

# auxiliary library names and versions (none)
