# This is root makefile of dyGiLa

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

# Absolute UNIX path of parallel io library Ascent
ifeq ($(ARCH), lumi)
ASCENT_DIR := /projappl/project_462000465/ascent/install/ascent-v0.9.0
else
  ifeq ($(ARCH), mahti)
   ASCENT_DIR := /projappl/project_2006478/a-2/install/ascent-v0.9.0
  endif
endif

# default ARCH if no ARCH is provided form shell
ifndef ARCH
 ARCH := vanilla
endif

# absolute UNIX path of dyGiLa folder
ifeq ($(ARCH), lumi)
 DYGILA_DIR := /projappl/project_462000465/dyGiLa
endif

APP_OPTS := -DNDIM=3
#-DEVEN_SITES_FIRST=0


# Set default goal and arch
.DEFAULT_GOAL := dyGiLa

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk     \
        $(DYGILA_DIR)/glsol/glsol_conf.mk \
        $(DYGILA_DIR)/pario/pario_conf.mk \
        $(DYGILA_DIR)/matep/matep_conf.mk

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir
dyGiLa: build/dyGiLa ; @:

# Now the linking step for each target executable
build/dyGiLa: Makefile $(GLSOL_OBJECTS) $(PARIO_OBJECTS) $(MATEP_OBJECTS) \
              build/main.o \
              $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ $(GLSOL_OBJECTS) $(PARIO_OBJECTS) $(MATEP_OBJECTS) \
              build/main.o \
	      $(HILA_OBJECTS) \
	      $(LDFLAGS) $(LDLIBS)
