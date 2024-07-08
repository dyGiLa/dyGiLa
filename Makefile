# set vpath for make to searching sources and headers
vpath %.cpp glsol/src pario/src matep/src

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
ASCENT_DIR := /projappl/project_462000465/spack/23.03/0.20.0/ascent-0.9.1-wsk2gp7
else ($(ARCH), mahti)
ASCENT_DIR := /projappl/project_2006478/a-2/install/ascent-v0.9.0
endif

APP_OPTS := -DNDIM=3
#-DEVEN_SITES_FIRST=0

# add headers searching directories
APP_OPTS += -I glsol/inc -I pario/inc -I matep/inc

USE_ASCENT := ON
ifeq ($(USE_ASCENT), ON)
  # See $(ASCENT_DIR)/share/ascent/ascent_config.mk for detailed linking info
  include $(ASCENT_DIR)/share/ascent/ascent_config.mk
  APP_OPTS += $(ASCENT_INCLUDE_FLAGS)
endif

# Set default goal and arch
.DEFAULT_GOAL := dygila

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir
dygila: build/dygila ; @:

# Now the linking step for each target executable
build/dygila: Makefile build/glsol.o build/pario.o build/matep.o build/main.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/glsol.o build/pario.o build/matep.o build/main.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) \
                                                                            $(ASCENT_LINK_RPATH) \
                                                                            $(ASCENT_MPI_LIB_FLAGS)


