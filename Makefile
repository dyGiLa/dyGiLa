# This is root makefile of dyGiLa

# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative.

# Absolute UNIX path of library HILA
# Absolute UNIX path of parallel io library Ascent
ifeq ($(ARCH), lumi)
HILA_DIR:= /projappl/project_462000836/insHILA-v0.0.2
else
  ifeq ($(ARCH), lumi-hip-CC)
   HILA_DIR:= /projappl/project_462000960/insHILA-main
  endif
  ifeq ($(ARCH), mahti)
   HILA_DIR:= /projappl/project_2006478/insHILA
  endif
  ifeq ($(ARCH), mahti-cuda)
   # HILA_DIR:= /projappl/project_2014552/insHILA-v0.0.2
   HILA_DIR:= /projappl/project_2014552/insHILA-main
  endif
endif

# default ARCH if no ARCH is provided form shell
ifndef ARCH
 ARCH := vanilla
endif

# absolute UNIX path of dyGiLa folder
ifeq ($(ARCH), lumi)
 DYGILA_DIR := /projappl/project_462000836/dyGiLa-develop
else
  ifeq ($(ARCH), lumi-hip-CC)
   DYGILA_DIR := /projappl/project_462000960/dyGiLa-develop-lite
  endif
  ifeq ($(ARCH), mahti)
   DYGILA_DIR := /projappl/project_2006478/dyGiLa-blob
  endif
  ifeq ($(ARCH), mahti-cuda)
   DYGILA_DIR := /projappl/project_2014552/dyGiLa-develop-lite
  endif
endif

APP_OPTS := -DNDIM=3
#-DEVEN_SITES_FIRST=0

# Set default goal and arch
.DEFAULT_GOAL := dyGiLa

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk  \
        $(DYGILA_DIR)/dyGiLa/dyGiLa_conf.mk \
        $(DYGILA_DIR)/glsol/glsol_conf.mk \
	$(DYGILA_DIR)/matep/matep_conf.mk 

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir
dyGiLa: build/dyGiLa ; @:

# Now the linking step for each target executable
build/dyGiLa: Makefile $(DYGILAAPIs_OBJECTS) $(GLSOL_OBJECTS) $(MATEP_OBJECTS) \
                build/main.o \
                $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ $(DYGILAAPIs_OBJECTS) $(GLSOL_OBJECTS) $(MATEP_OBJECTS) \
                build/main.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
