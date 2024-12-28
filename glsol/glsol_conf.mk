# Makefie of TDGL-Langvian eqns solver glsol

# *.cpp files searching path
vpath %.cpp glsol/src glsol/src/utilities glsol/src/next glsol/src/initialize

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

# pario objects, built by HILA pattern rules
GLSOL_OBJECTS = build/allocate.o              \
                build/write_energies.o        \
                build/write_positions.o       \
                build/write_phases.o          \
                build/write_moduli.o          \
                build/next.o                  \
                build/next_bath.o             \
                build/next_bath_UniT_quench.o \
                build/next_bath_UniT_quench_Hfield.o \
                build/next_bath_hotblob_quench_Hfield.o \
                build/next_T.o                \
                build/glsol_initialize.o      \
                build/glsol_initialize_T.o    \
                build/glsol_initialize_H.o    \
                build/point_params.o          

.PHONY: glsol
glsol: $(GLSOL_OBJECTS)
