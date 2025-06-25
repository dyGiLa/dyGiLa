# Makefie of parallel IO engine pario

# *.cpp files searching path
vpath %.cpp pario/src pario/src/utilities pario/src/xml

# Add Ascent include path, linder flags
# of path and binary libs into
# prerequisites & recipes of dyGiLa
# building and linking rule
USE_PARIO := ON
ifeq ($(USE_PARIO), ON)
  include $(ASCENT_DIR)/share/ascent/ascent_config.mk
  APP_OPTS += $(ASCENT_INCLUDE_FLAGS)
  LDFLAGS  += $(ASCENT_LINK_RPATH)
  LDLIBS   += $(ASCENT_MPI_CUDA_LIB_FLAGS)
endif

# add headers searching directories
APP_OPTS += -I $(DYGILA_DIR)/pario/inc

# pario objects, built by HILA pattern rules
PARIO_OBJECTS = build/xdmf.o     \
                build/xml_Amatrix.o \
                build/xml_pMarker.o \
                build/xml_massCurrent.o \
                build/xml_spinCurrent.o \
                build/pstream.o  \
                build/init.o     \
                build/shutdown.o \
                build/mesh.o     \
                build/mesh_gapA_FEDensity.o   \
                build/mesh_insitu_Temperature.o \
                build/mesh_phaseMarker.o \
                build/mesh_massCurrent.o      \
                build/mesh_spinCurrent.o      \
                build/mesh_AMatrix.o          \
                build/mesh_addGhost_verify.o  \
                build/actions_insitu.o  \
                build/actions_massCurrent.o \
                build/actions_spinCurrent.o \
                build/actions_AMatrix.o     \
                build/actions_phaseMarker.o \
                build/actions_printTree.o

.PHONY: pario
pario: $(PARIO_OBJECTS)
