# foreign

# Files to compile
SET ( FILES
  csvsurfacelist
  dehydration
  isosig
  orb
  pdf
  pyforeign
  recogniser
  snappea
  )

# Prepend folder name
FOREACH ( SOURCE_FILE ${FILES} )
  SET ( SOURCES ${SOURCES} foreign/${SOURCE_FILE})
ENDFOREACH(SOURCE_FILE)
SET(SOURCES ${SOURCES} PARENT_SCOPE)
