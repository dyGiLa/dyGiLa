# set vpath for make to searching sources and headers
vpath %.cpp glsol/src pario/src pario/src/utilities matep/src

# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative.
# If one want use cosfithe HILA instance as library, one could pass
# HILA_DIR=cosfithe option in command line
ifndef HILA_DIR
HILA_DIR:= /projappl/project_462000465/insHILA
else
  ifeq (${HILA_DIR}, cosfithe)
    HILA_DIR := /projappl/cosfithe/HILA
  endif
endif

# default ARCH if no ARCH is provided form shell
ifndef ARCH
ARCH := vanilla
endif

# Absolute UNIX path of parallel io library Ascent
ifeq ($(ARCH), lumi)
 ASCENT_DIR := /projappl/project_462000465/ascent/install/ascent-v0.9.0
else ($(ARCH), mahti)
 ASCENT_DIR := /projappl/project_2006478/a-2/install/ascent-v0.9.0
endif

APP_OPTS := -DNDIM=3
#-DEVEN_SITES_FIRST=0

# add headers searching directories
APP_OPTS += -I /projappl/project_462000465/dyGiLa/glsol/inc \
            -I /projappl/project_462000465/dyGiLa/pario/inc \
            -I /projappl/project_462000465/dyGiLa/matep/inc

# add headers path of ffw
APP_OPTS += -I /projappl/project_462000465/lib/fftw-3.3.10-fftw3f/include \
            -I /projappl/project_462000465/lib/fftw-3.3.10-fftw3/include        

# Set default goal and arch
.DEFAULT_GOAL := dyGiLa

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

USE_PARIO := ON
ifeq ($(USE_PARIO), ON)
  include $(ASCENT_DIR)/share/ascent/ascent_config.mk
  APP_OPTS += $(ASCENT_INCLUDE_FLAGS)
  LDFLAGS  += $(ASCENT_LINK_RPATH)
  LDLIBS   += $(ASCENT_MPI_LIB_FLAGS)
endif

LDFLAGS += -L/projappl/project_462000465/lib/fftw-3.3.10-fftw3f/lib \
           -L/projappl/project_462000465/lib/fftw-3.3.10-fftw3/lib

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir
dyGiLa: build/dyGiLa ; @:

# Now the linking step for each target executable
build/dyGiLa: Makefile build/glsol.o \
	      build/xdmf.o build/pstream.o build/init.o build/shutdown.o build/mesh.o build/actions.o \
	      build/matep.o \
              build/main.o \
              $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/glsol.o \
              build/xdmf.o build/pstream.o build/init.o build/shutdown.o build/mesh.o build/actions.o \
              build/matep.o \
              build/main.o \
	      $(HILA_OBJECTS) \
	      $(LDFLAGS) $(LDLIBS)
