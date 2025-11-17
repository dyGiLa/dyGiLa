# Makefie of dyGiLa APIs

# *.cpp files searching path
vpath %.cpp dyGiLa/src dyGiLa/src/utils

# add headers searching directories
APP_OPTS += -I $(DYGILA_DIR)/dyGiLa/inc 

# APIs objects, built by HILA pattern rules
DYGILAAPIs_OBJECTS = build/nextBlocks.o \
                     build/pstreaming.o \
                     build/gammaEvolve.o \
                     build/dyGiLaInit.o \
                     build/heterogeneousQuench.o \
                     build/homogenousQuench.o \
                     build/dyGiLaPhaseMarking.o
.PHONY: dyGiLaAPIs
dyGiLaAPIs: $(DYGIALAAPIs_OBJECTS)
