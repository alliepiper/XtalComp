set(spglib_SRCS
  arithmetic.c
  arithmetic.h
  cell.c
  cell.h
  debug.c
  debug.h
  delaunay.c
  delaunay.h
  hall_symbol.c
  hall_symbol.h
  kgrid.c
  kgrid.h
  kpoint.c
  kpoint.h
  mathfunc.c
  mathfunc.h
  niggli.c
  niggli.h
  pointgroup.c
  pointgroup.h
  primitive.c
  primitive.h
  refinement.c
  refinement.h
  sitesym_database.c
  sitesym_database.h
  site_symmetry.c
  site_symmetry.h
  spacegroup.c
  spacegroup.h
  spg_database.c
  spg_database.h
  spglib.c
  spglib_f.c
  spglib.h
  spin.c
  spin.h
  symmetry.c
  symmetry.h
  test.c
  tetrahedron_method.c
  tetrahedron_method.h
  version.h
)

add_library(spglib STATIC ${spglib_SRCS})
# Suppress spglib warnings
set_target_properties(spglib PROPERTIES
  COMPILE_FLAGS "-w"
)
