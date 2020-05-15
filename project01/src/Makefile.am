bin_PROGRAMS = solver
solver_SOURCES = main.c laplace.h function.h input.c initialize.c build_linear_system_1d.c build_linear_system_2d.c solve_linear_system.c output.c res_norm.c error_l2.c show_linear_system.c

AM_LDFLAGS = -lm
AM_CFLAGS = $(GRVY_CFLAGS) $(MASA_CXXFLAGS)
LIBS = $(GRVY_LIBS) $(MASA_LIBS)
