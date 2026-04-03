# This is Makefile of irtb
# It sets up the target irtb

# *.cpp files searching path
vpath %.cpp irtb/src

# add headers searching directories
APP_OPTS += -I $(DYGILA_DIR)/irtb/inc

# Object Path
IRTB_PATH = build/Targets

# IRTB object, built by HILA pattern rules
IRTB_OBJECTS = $(IRTB_PATH)/init_pWaveIRTB.o

.PHONY: irtb
irtb: $(IRTB_OBJECTS)
