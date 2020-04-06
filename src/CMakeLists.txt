#----------------------------------------------------------------------
# Include project files 
#----------------------------------------------------------------------
# include all the Example project .cpp files
INCLUDE(${PhOVS_SOURCE_DIR}/CMakeSourceFileList.cmake)
# include all the Example project .h files
INCLUDE(${PhOVS_SOURCE_DIR}/CMakeHeaderFileList.cmake)


#----------------------------------------------------------------------
# Create rule to copy all the headers from src to include/PhOVS
#----------------------------------------------------------------------
# For each header, we create a rule
SET(HEADER_IN_INCLUDE_DIR "")
SET(ADDED_SRC_DEPEND_PROPERTIES "")
FOREACH(header ${HEADER_ALL})
  GET_FILENAME_COMPONENT(headerName ${header} NAME) 
  ADD_CUSTOM_COMMAND(
    OUTPUT ${PhOVS_BINARY_DIR}/include/PhOVS/${headerName}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${header} ${PhOVS_BINARY_DIR}/include/PhOVS/${headerName}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${header}
    )
  SET(HEADER_IN_INCLUDE_DIR ${HEADER_IN_INCLUDE_DIR} ${PhOVS_BINARY_DIR}/include/PhOVS/${headerName}
    )
ENDFOREACH(header)


#----------------------------------------------------------------------
# Create a top level rule to copy all the headers from src to include/PhOVS
#----------------------------------------------------------------------
ADD_CUSTOM_TARGET(header ALL
  DEPENDS ${HEADER_IN_INCLUDE_DIR}
  )


#----------------------------------------------------------------------
# build rule for the library
#----------------------------------------------------------------------
ADD_LIBRARY (${PhOVS_INTERN_LIBS}
  ${SRC_ALL} 
  ${HEADER_ALL}
  ) 

# create the headers in include/PhOVS before compiling the lib
ADD_DEPENDENCIES(${PhOVS_INTERN_LIBS} header)

# Append the library version information to the library target
# properties. 
SET_TARGET_PROPERTIES(${PhOVS_INTERN_LIBS}
  PROPERTIES
  VERSION ${PhOVS_VERSION_MAJOR}.${PhOVS_VERSION_MINOR}
  SOVERSION ${PhOVS_VERSION_PATCH}
  )

#----------------------------------------------------------------------
# customize install target 
#----------------------------------------------------------------------
# install rule for the library
INSTALL_TARGETS(/lib ${PhOVS_INTERN_LIBS})
# install rule for all the headers
INSTALL_FILES(/include/PhOVS FILES ${HEADER_ALL})
