Large Test Suite: Exhaustive Lists of Normal Surfaces
-----------------------------------------------------

This directory tree contains complete lists of normal surfaces for
manifolds in the closed orientable and non-orientable censuses.

The main reason for these lists is so that, when the normal surface
enumeration algorithm is changed, the author can verify that the new code
produces the same results as the old code.

The files in this directory were created using Regina 4.5.  The
subdirectories "standard", "quad" and "an-standard" contain surfaces in
tri-quad, quad and tri-quad-oct coordinates respectively.  The data
files (stored as gzipped plain text) are named "or-<n>" and "nor-<n>",
which indicates that they run through all manifolds from the n-tetrahedron
closed orientable or non-orientable census respectively.

The regina-python scripts that produce the surfaces are stored alongside
the data files in each subdirectory.  Each script requires a census file
as an argument (either closed-or-census-large.rga or closed-nor-census.rga,
both of which are shipped with regina).

WARNING: These scripts do not delete normal surface lists, which means
they can consume a *very* large amount of memory when run over a large
census.  This is due to ownership constraints in regina-python, and will
be fixed if the need becomes pressing.

