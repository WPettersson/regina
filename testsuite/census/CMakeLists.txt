# census

# Files to compile
SET ( FILES
  dim2census
  dim2edgepairing
  dim4census
  dim4facetpairing
  ncensus
  nfacepairing
  )

# Prepend folder name
FOREACH ( SOURCE_FILE ${FILES} )
  SET ( SOURCES ${SOURCES} census/${SOURCE_FILE})
ENDFOREACH(SOURCE_FILE)
SET(SOURCES ${SOURCES} PARENT_SCOPE)

