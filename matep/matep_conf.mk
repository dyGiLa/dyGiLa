# This is Makefile of matep
# It sets up the target matep

# *.cpp files searching path
vpath %.cpp matep/src

# add headers searching directories
APP_OPTS += -I $(DYGILA_DIR)/matep/inc

# matep object, built by HILA pattern rules
MATEP_OBJECTS = build/matep.o \
	        build/matep_utils.o \
	        build/init_global_wrapper_mp.o

.PHONY: matep
matep: $(MATEP_OBJECTS)
