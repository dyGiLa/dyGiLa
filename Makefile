# This is root makefile of dyGiLa

# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative.
# If one want use cosfithe HILA instance as library, one could pass
# HILA_DIR=cosfithe option in command line

# Absolute UNIX path of library HILA
# Absolute UNIX path of parallel io library Ascent
ifeq ($(ARCH), lumi)
HILA_DIR:= /projappl/project_462000465/insHILA-II
ASCENT_DIR := /projappl/project_462000465/ascent/install/ascent-v0.9.0
#ASCENT_DIR := /projappl/project_462000465/ascent-0.9.3/scripts/build_ascent/install/ascent-develop
else
  ifeq ($(ARCH), lumi-hip-CC)
   HILA_DIR:= /projappl/project_462000809/insHILA
   ASCENT_DIR := /projappl/project_462000809/ascent/scripts/build_ascent/install/ascent-develop
  endif
  ifeq ($(ARCH), mahti)
   HILA_DIR:= /projappl/project_2006478/insHILA
   ASCENT_DIR := /projappl/project_2006478/a-2/install/ascent-v0.9.0
  endif
  ifeq ($(ARCH), mahti-cuda)
   HILA_DIR:= /projappl/project_2006478/insHILA
#   ASCENT_DIR := /projappl/project_2006478/a-2/install/ascent-v0.9.0
  endif
endif

# default ARCH if no ARCH is provided form shell
ifndef ARCH
 ARCH := vanilla
endif

# absolute UNIX path of dyGiLa folder
ifeq ($(ARCH), lumi)
 DYGILA_DIR := /projappl/project_462000465/dyGiLa-blob
else
  ifeq ($(ARCH), lumi-hip-CC)
   DYGILA_DIR := /projappl/project_462000809/dyGiLa
  endif
  ifeq ($(ARCH), mahti)
   DYGILA_DIR := /projappl/project_2006478/dyGiLa-blob
  endif
  ifeq ($(ARCH), mahti-cuda)
   DYGILA_DIR := /projappl/project_2006478/dyGiLa-GPU
  endif
endif

APP_OPTS := -DNDIM=3
#-DEVEN_SITES_FIRST=0


# Set default goal and arch
.DEFAULT_GOAL := dyGiLa

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk     \
        $(DYGILA_DIR)/glsol/glsol_conf.mk \
	$(DYGILA_DIR)/matep/matep_conf.mk \
        $(DYGILA_DIR)/pario/pario_conf.mk

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

# build/dyGiLa: Makefile $(GLSOL_OBJECTS) $(MATEP_OBJECTS) \
#               build/main.o \
#               $(HILA_OBJECTS) $(HEADERS)
# 	$(LD) -o $@ $(GLSOL_OBJECTS) $(MATEP_OBJECTS) \
#               build/main.o \
# 	      $(HILA_OBJECTS) \
# 	      $(LDFLAGS) $(LDLIBS)
