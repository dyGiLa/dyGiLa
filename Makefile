# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
HILA_DIR := /projappl/cosfithe/HILA

APP_OPTS := -DNDIM=3
#-DEVEN_SITES_FIRST=0

# Set default goal and arch
.DEFAULT_GOAL := he3sim

# ARCH is important on mahti, 
# make -j8(24) ARCH=mahti
ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

he3sim: build/he3sim ; @:

# Now the linking step for each target executable
build/he3sim: Makefile build/he3sim.o build/matep.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/he3sim.o build/matep.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)



