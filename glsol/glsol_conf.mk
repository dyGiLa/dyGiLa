# Makefie of TDGL-Langvian eqns solver glsol

# *.cpp files searching path
vpath %.cpp glsol/src glsol/src/utilities glsol/src/next glsol/src/initialize \
$(DYGILA_DIR)/utils/GLfeContri $(DYGILA_DIR)/utils/ABOBA $(DYGILA_DIR)/utils/confCatchnRelax $(DYGILA_DIR)/utils/confInitialize

# Include path, linder flags
# of path and binary libs into
# prerequisites & recipes of dyGiLa
# building and linking rule

# add headers searching directories 
APP_OPTS += -I $(DYGILA_DIR)/glsol/inc -I $(DYGILA_DIR)/utils/AdGR

# ifeq ($(ARCH), lumi)
# # add headers path of ffw
# APP_OPTS += -I /projappl/project_462000465/lib/fftw-3.3.10-fftw3f/include \
#             -I /projappl/project_462000465/lib/fftw-3.3.10-fftw3/include

# # fft, fftwl libraries binary path
# LDFLAGS += -L/projappl/project_462000465/lib/fftw-3.3.10-fftw3f/lib \
#            -L/projappl/project_462000465/lib/fftw-3.3.10-fftw3/lib
# endif

# Objects path
GLSOL_PATH = build/Targets

# pario objects, built by HILA pattern rules
GLSOL_OBJECTS = $(GLSOL_PATH)/configure.o             \
                $(GLSOL_PATH)/glsol.o                 \
                $(GLSOL_PATH)/fstreams.o              \
                $(GLSOL_PATH)/write_energies.o        \
                $(GLSOL_PATH)/write_positions.o       \
                $(GLSOL_PATH)/write_phases.o          \
                $(GLSOL_PATH)/gaussianLP_matrix.o     \
                $(GLSOL_PATH)/phaseMarking.o          \
                $(GLSOL_PATH)/phaseCounting.o         \
                $(GLSOL_PATH)/next.o                  \
                $(GLSOL_PATH)/next_bath.o             \
                $(GLSOL_PATH)/next_AdGRz_bath.o       \
                $(GLSOL_PATH)/next_bath_UniT_quench.o \
                $(GLSOL_PATH)/next_bath_UniT_quench_Hfield.o \
                $(GLSOL_PATH)/next_bath_UniT_quench_AdGRz_Hfield.o \
                $(GLSOL_PATH)/next_bath_hotblob_quench_Hfield.o \
                $(GLSOL_PATH)/next_bath_hotblob_quench_Hfield_confCatch.o \
                $(GLSOL_PATH)/next_bath_UniT_quench_AdGRz_Hfield_confCatch.o \
                $(GLSOL_PATH)/glsol_initialize.o      \
                $(GLSOL_PATH)/glsol_initialize_T.o    \
                $(GLSOL_PATH)/glsol_initialize_H.o    \
                $(GLSOL_PATH)/case_0.o $(GLSOL_PATH)/case_1.o $(GLSOL_PATH)/case_2.o \
                $(GLSOL_PATH)/case_3.o $(GLSOL_PATH)/case_4.o $(GLSOL_PATH)/case_5.o \
                $(GLSOL_PATH)/case_6.o $(GLSOL_PATH)/case_7.o $(GLSOL_PATH)/case_8.o \
                $(GLSOL_PATH)/case_9.o $(GLSOL_PATH)/case_10.o \
                $(GLSOL_PATH)/dPiGLfe_AdGRz.o $(GLSOL_PATH)/dPiGLfe.o \
                $(GLSOL_PATH)/ABOBA.o $(GLSOL_PATH)/ABOBA_gBranch.o \
                $(GLSOL_PATH)/dampAndRelax.o $(GLSOL_PATH)/relax.o


.PHONY: glsol
glsol: $(GLSOL_OBJECTS)
