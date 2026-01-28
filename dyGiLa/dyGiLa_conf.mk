# Makefie of dyGiLa APIs

# *.cpp files searching path
vpath %.cpp dyGiLa/src dyGiLa/src/utils

# add headers searching directories
APP_OPTS += -I $(DYGILA_DIR)/dyGiLa/inc

# Objects path
DYGILAAPIs_PATH = build/Targets

# APIs objects, built by HILA pattern rules
DYGILAAPIs_OBJECTS = $(DYGILAAPIs_PATH)/nextBlocks.o \
                     $(DYGILAAPIs_PATH)/pstreaming.o \
                     $(DYGILAAPIs_PATH)/gammaEvolve.o \
                     $(DYGILAAPIs_PATH)/dyGiLaInit.o \
                     $(DYGILAAPIs_PATH)/heterogeneousQuench.o \
                     $(DYGILAAPIs_PATH)/homogenousQuench.o \
                     $(DYGILAAPIs_PATH)/dyGiLaPhaseMarking.o

.PHONY: dyGiLaAPIs
dyGiLaAPIs: $(DYGIALAAPIs_OBJECTS)
