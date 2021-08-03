#
# Copyright Projet Lagadic / IRISA-INRIA Rennes, 2006
# www: http://www.irisa.fr/lagadic
#
# Author: Fabien Spindler email:Fabien.Spindler@irisa.fr
#
# SOURCE variable corresponds to the list of all the sources to build binaries.
# The generate binary comes by removing the .cpp extension to
# the source name.
#
# If you want to add/remove a source, modify here
SET (EXE_SOURCE
		DDVS.cpp
)

SET (SRC_SOURCES
		src/C_UR.h
		src/C_UR.cpp

    src/CamFlir.hpp
)

# Libraries
# ORFCV
#FIND_PACKAGE(ORFCV REQUIRED)
#IF(ORFCV_FOUND)
#  INCLUDE(${ORFCV_USE_FILE})
#ENDIF(ORFCV_FOUND)
#LINK_DIRECTORIES(${ORFCV_LINK_DIRECTORIES})

#MESSAGE("ORFCV LINK DIRECTORIES : ${ORFCV_LINK_DIRECTORIES}")

#INCLUDE_DIRECTORIES(
#	${ORFCV_INCLUDE_DIR}
#	${ORFCV_EXTERN_INCLUDE_DIR}
#  )

# Libraries
#TARGET_LINK_LIBRARIES(${binary} 
#${ORFCV_LIBRARIES}
#)

#ur_rtde
find_package(ur_rtde REQUIRED)

INCLUDE_DIRECTORIES(/usr/include/flycapture)

# rule for binary build

FOREACH(source ${EXE_SOURCE})
  # Compute the name of the binary to create
  GET_FILENAME_COMPONENT(binary ${source} NAME_WE)

  # From source compile the binary and add link rules
  ADD_EXECUTABLE(${binary} ${source} ${SRC_SOURCES})

  TARGET_LINK_LIBRARIES(${binary} ${PhOVS_INTERN_LIBS} ${PhOVS_EXTERN_LIBS} -lflycapture ur_rtde::rtde)

  # Add test
  ADD_TEST(${binary} ${binary})

  MESSAGE(${binary})

ENDFOREACH(source)

# customize clean target 
SET_DIRECTORY_PROPERTIES(PROPERTIES 
  ADDITIONAL_MAKE_CLEAN_FILES "core*;*~;gmon.out;DartTestfile.txt"
)

