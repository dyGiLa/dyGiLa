# Makefie of parallel IO engine pario

# *.cpp files searching path
vpath %.cpp pario/src pario/src/utilities

# Add Ascent include path, linder flags
# of path and binary libs into
# prerequisites & recipes of dyGiLa
# building and linking rule
USE_PARIO := ON
ifeq ($(USE_PARIO), ON)
  include $(ASCENT_DIR)/share/ascent/ascent_config.mk
  APP_OPTS += $(ASCENT_INCLUDE_FLAGS)
  LDFLAGS  += $(ASCENT_LINK_RPATH)
  LDLIBS   += $(ASCENT_MPI_LIB_FLAGS)
endif

# add headers searching directories
APP_OPTS += -I $(DYGILA_DIR)/pario/inc

# pario objects, built by HILA pattern rules
PARIO_OBJECTS = build/xdmf.o     \
                build/xml.o      \
                build/pstream.o  \
                build/init.o     \
                build/shutdown.o \
                build/mesh.o     \
                build/actions.o

.PHONY: pario
pario: $(PARIO_OBJECTS)
