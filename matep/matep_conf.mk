# This is Makefile of matep
# It sets up the target matep

# *.cpp files searching path
vpath %.cpp matep/src

# add headers searching directories
APP_OPTS += -I $(DYGILA_DIR)/matep/inc

# Object Path
MATEP_PATH = build/Targets

# matep object, built by HILA pattern rules
MATEP_OBJECTS = $(MATEP_PATH)/matep.o \
	        $(MATEP_PATH)/matep_utils.o \
	        $(MATEP_PATH)/init_global_wrapper_mp.o

.PHONY: matep
matep: $(MATEP_OBJECTS)
