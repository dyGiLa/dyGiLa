# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
HILA_DIR := /projappl/cosfithe/HILA

ASCENT_DIR := /projappl/project_2006478/ascent/install-debug

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

# See $(ASCENT_DIR)/share/ascent/ascent_config.mk for detailed linking info
include $(ASCENT_DIR)/share/ascent/ascent_config.mk

# make sure to enable c++11 support (conduit's interface now requires it)
CXX_FLAGS = -std=c++11                                               
INC_FLAGS = $(ASCENT_INCLUDE_FLAGS)                                  
LNK_FLAGS = $(ASCENT_LINK_RPATH) $(ASCENT_MPI_LIB_FLAGS)    

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

he3sim: build/he3sim ; @:

# Now the linking step for each target executable
build/he3sim: Makefile build/he3sim.o build/matep.o $(HILA_OBJECTS) $(HEADERS)
	$(CXX) $(CXX_FLAGS) $(INC_FLAGS) he3sim.cpp 
	$(LD) $(LNK_FLAGS) -o $@ build/he3sim.o build/matep.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)



