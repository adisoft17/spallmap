# -*-Mode: Makefile-*-
#
#    C pre-processor flags used in the TELLURIDE GNUmakefile
#    that are specific to TELLURIDE (machine-independent)
# 
#
# This listing of TELLURIDE C pre-processor flags is included in the 
# GNUmakefile.  It defines the list of C-preprocessor parameters that
# are to be enabled upon compilation.  These are NOT machine-specific
# and are ALWAYS used in the compilation.
#
# To toggle any of the options below, simply add or remove the comment
# character (#) in column 1.  Add new entries for general use as needed.

# =-=-=-=-=-=-=-=-=-=-=-=-=-= CPP Flags =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

CPPFLAGS = -P

# enable a more verbose debug mode
CPPFLAGS += -DDEBUG

# enable timing
CPPFLAGS += -DTIMERS

# define the source tags
CPPFLAGS += -DCODE_NAME=\'"$(CODE)"\'
CPPFLAGS += -DVERSION=\'"$(VERSION)"\'
CPPFLAGS += -DBUILD_DATE=\'"$(BUILD_DATE)"\'
CPPFLAGS += -DHOST_NAME=\'"$(HOST_NAME)"\'
CPPFLAGS += -DARCHITECTURE=\'"$(ARCHITECTURE)"\'

## use Augustus
#ifeq ($(use_Augustus),yes)
#  CPPFLAGS += -DUSE_AUGUSTUS
#  CPPFLAGS += -DAUGUSTUS=\'"$(AUGUSTUS)"\'
#endif

# compiler and operating system specific cpp flags, for trouble
# turn them on if needed, turn off every once in a while to see
# if the trouble has been fixed

#ifeq ($(os),dunix)
#  CPPFLAGS += -DDUNIX_WORKAROUND
#endif

#ifeq ($(os),unicos)
#  CPPFLAGS += -DUNICOS_WORKAROUND
#endif

#ifeq ($(os),irix)
#  CPPFLAGS += -DIRIX_WORKAROUND
#endif

#ifeq ($(os),sunos)
#  CPPFLAGS += -DSUNOS_WORKAROUND
#endif

#ifeq ($(os),solaris)
#  CPPFLAGS += -DSOLARIS_WORKAROUND
#endif

ifeq ($(os),aix)
  CPPFLAGS += -DAIX_WORKAROUND
endif

#ifeq ($(os),hpux)
#  CPPFLAGS += -DHPUX_WORKAROUND
#endif

#ifeq ($(compiler),dec)
#  CPPFLAGS += -DDEC_COMPILER_WORKAROUND
#endif

#ifeq ($(compiler),cray)
#  CPPFLAGS += -DCRAY_COMPILER_WORKAROUND
#endif

#ifeq ($(compiler),sgi)
#  CPPFLAGS += -DSGI_COMPILER_WORKAROUND
#endif

ifeq ($(compiler),sun)
  CPPFLAGS += -DSUN_COMPILER_WORKAROUND
endif

#ifeq ($(compiler),fujitsu)
#  CPPFLAGS += -DFUJITSU_COMPILER_WORKAROUND
#endif

#ifeq ($(compiler),ibm)
#  CPPFLAGS += -DIBM_COMPILER_WORKAROUND
#endif

#ifeq ($(compiler),hp)
#  CPPFLAGS += -DHP_COMPILER_WORKAROUND
#endif
