bin_PROGRAMS = solver
solver_SOURCES = main.c laplace.h function.h input.c initialize.c build_linear_system_1d.c build_linear_system_2d.c solve_linear_system.c solve_linear_system_PETSC.c output.c res_norm.c error_l2.c show_linear_system.c 

AM_CFLAGS = $(GRVY_CFLAGS) $(MASA_CXXFLAGS) $(HDF5_CFLAGS)
LIBS = $(GRVY_LIBS) $(MASA_LIBS) $(HDF5_LIBS)

if PETSC_ENABLED
include $(PETSC_DIR)/lib/petsc/conf/variables
AM_CFLAGS += -DINCLUDE_PETSC $(PETSC_CC_INCLUDES)
LIBS += ${PETSC_LIB}
endif
