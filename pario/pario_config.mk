# Makefie of parallel IO engine pario

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

# pario objects, built by HILA pattern rules
PARIO_OBJECTS = build/xdmf.o \
                build/pstream.o \
                build/init.o \
                build/shutdown.o \
                build/mesh.o \
                build/actions.o

.PHONY: pario
pario: $(PARIO_OBJECTS)
