# Makefie of dyGiLa orchestrator

# *.cpp files searching path
vpath %.cpp orch/src orch/src/utils

# add headers searching directories
APP_OPTS += -I $(DYGILA_DIR)/orch/inc 

# orchestrator objects, built by HILA pattern rules
ORCH_OBJECTS = build/writeHDF5_xmls.o \
               build/nextBlocks.o \
               build/pstreaming.o \
               build/gammaEvolve.o \
               build/dyGiLaInit.o \
               build/heterogeneousQuench.o \
               build/homogenousQuench.o
.PHONY: orch
orch: $(ORCH_OBJECTS)
