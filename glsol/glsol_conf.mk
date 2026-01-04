# Makefie of TDGL-Langvian eqns solver glsol

# *.cpp files searching path
vpath %.cpp glsol/src glsol/src/utilities glsol/src/next glsol/src/initialize \
$(DYGILA_DIR)/utils/GLfeContri $(DYGILA_DIR)/utils/ABOBA $(DYGILA_DIR)/utils/confCatchnRelax $(DYGILA_DIR)/utils/confInitialize 

# Include path, linder flags
# of path and binary libs into
# prerequisites & recipes of dyGiLa
# building and linking rule

# add headers searching directories 
APP_OPTS += -I $(DYGILA_DIR)/glsol/inc

# ifeq ($(ARCH), lumi)
# # add headers path of ffw
# APP_OPTS += -I /projappl/project_462000465/lib/fftw-3.3.10-fftw3f/include \
#             -I /projappl/project_462000465/lib/fftw-3.3.10-fftw3/include

# # fft, fftwl libraries binary path
# LDFLAGS += -L/projappl/project_462000465/lib/fftw-3.3.10-fftw3f/lib \
#            -L/projappl/project_462000465/lib/fftw-3.3.10-fftw3/lib
# endif

# glsol objects, built by HILA pattern rules
GLSOL_OBJECTS = build/configure.o             \
                build/glsol.o                 \
                build/fstreams.o              \
                build/write_energies.o        \
                build/write_positions.o       \
                build/write_phases.o          \
                build/gaussianLP_matrix.o     \
                build/phaseMarking.o          \
                build/phaseCounting.o         \
                build/next.o                  \
                build/next_bath.o             \
                build/next_AdGRz_bath.o       \
                build/next_bath_UniT_quench.o \
                build/next_bath_UniT_quench_Hfield.o \
                build/next_bath_UniT_quench_AdGRz_Hfield.o \
                build/next_bath_hotblob_quench_Hfield.o \
                build/next_bath_hotblob_quench_Hfield_confCatch.o \
                build/next_bath_UniT_quench_AdGRz_Hfield_confCatch.o \
                build/glsol_initialize.o      \
                build/glsol_initialize_T.o    \
                build/glsol_initialize_H.o    \
                build/dPiGLfe.o               \
                build/ABOBA.o                 \
                build/dampAndRelax.o          \
                build/relax.o                 \
                build/case_0.o build/case_1.o build/case_2.o \
                build/case_3.o build/case_4.o build/case_5.o \
                build/case_6.o build/case_7.o build/case_8.o \
                build/case_9.o 

.PHONY: glsol
glsol: $(GLSOL_OBJECTS)
