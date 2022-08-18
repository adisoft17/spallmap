# -*-Mode: Makefile-*-
#
#   Code name and version definitions used in the SYNGAS GNUmakefile
# 

#CODE    = SPALLMAP 
CODE    = spallmap
VERSION = 1.0.0

# define LIB name to be the same as the CODE name
LIB = $(CODE)

# define symbols to be #ifdef'ed into the code
BUILD_DATE   = $(shell date)
ARCHITECTURE = $(shell uname -a)
HOST_NAME    = $(shell uname -n)

