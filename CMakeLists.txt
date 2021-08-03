CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)
PROJECT(PhOVS)

#ADD_DEFINITIONS(-D_BMT_SAFE -DREENTRANT -D_BNOLICENSING -D_BGL_SUPPORT -D_BX11_SUPPORT -D_BCITYVIP)

#-----------------------------------------------------------------------------
# OMNIRENDERING version number.  An even minor number corresponds to releases.
SET(PhOVS_VERSION_MAJOR 1)
SET(PhOVS_VERSION_MINOR 0)
SET(PhOVS_VERSION_PATCH 0)
SET(PhOVS_VERSION "${PhOVS_VERSION_MAJOR}.${PhOVS_VERSION_MINOR}")
SET(PhOVS_VERSION_FULL "${PhOVS_VERSION}.${PhOVS_VERSION_PATCH}")      

SET(PhOVS_INCLUDE_DIR "${PhOVS_BINARY_DIR}/include/PhOVS")
SET(PhOVS_DOC_DIR "${PhOVS_BINARY_DIR}/doc")

#-----------------------------------------------------------------------------

# The location in which to install Example libraries.
SET(LIBRARY_OUTPUT_PATH ${PhOVS_BINARY_DIR}/lib)

# The location in which to install PhOVS binaries.
SET(BINARY_OUTPUT_PATH ${PhOVS_BINARY_DIR}/bin)

# the include directory we depend on
SET(PhOVS_INTERN_INCLUDE_DIR ${PhOVS_BINARY_DIR}/include)
# the third party include directories we depend on
SET(PhOVS_EXTERN_INCLUDE_DIR "")
# the internal visp libs we want to createn
SET(PhOVS_INTERN_LIBS PhOVS)
# the third party libraries we depend on
SET(PhOVS_EXTERN_LIBS "")
# the third party libraries we depend on
SET(PhOVS_DEFS "")

# Get the OS
SET(OS ${CMAKE_SYSTEM_NAME})

## where are user-specific cmake modules
## make_dist.sh relies on $ENV{CMAKE_MODULES_PATH} (TODO, FIXME)
SET(CMAKE_MODULE_PATH ${PhOVS_SOURCE_DIR}/CMakeModules)

INCLUDE(CheckLibraryExists)

# Create include directory which will contain all the headers
MAKE_DIRECTORY(${PhOVS_INCLUDE_DIR})



#---------------------------------------------------------------------
# Platform dependent 
#----------------------------------------------------------------------
IF(WIN32 AND NOT CYGWIN)
  SET(PhOVS_DEFS "${PhOVS_DEFS} -DWIN32")
  SET(PhOVS_DEFS "${PhOVS_DEFS} -DGX") # To support exceptions
  CHECK_LIBRARY_EXISTS("winmm"  getch "" HAVE_LIBWINMM) # for timeGetTime()
  IF(HAVE_LIBWINMM)
    SET(PhOVS_EXTERN_LIBS ${PhOVS_EXTERN_LIBS} "winmm")
  ENDIF(HAVE_LIBWINMM)
  IF(MSVC) #check if MSVC is the compiler and set ART_USE_MSVC accordingly
    SET(ART_USE_MSVC   TRUE)  # for header mbtConfig.h
  ENDIF(MSVC)
 IF(MSVC) #check if MSVC is the compiler and set ART_USE_MSVC accordingly
    SET(MBT_USE_MSVC   TRUE)  # for header mbtConfig.h
  ENDIF(MSVC)
ENDIF(WIN32 AND NOT CYGWIN)
IF(UNIX)
  SET(PhOVS_DEFS "${PhOVS_DEFS} -DUNIX")
  IF(APPLE)
    # We will not use rpath support.  Tell CMake to skip it.
    SET(CMAKE_SKIP_RPATH 1 CACHE INTERNAL "Whether to build with rpath." FORCE)
    SET(PhOVS_DEFS "${PhOVS_DEFS} -DAPPLE")
  ENDIF(APPLE)
ENDIF(UNIX)

#----------------------------------------------------------------------
# Third party libraries
#----------------------------------------------------------------------

# ViSP-2
FIND_PACKAGE(VISP REQUIRED visp_core visp_io visp_robot visp_visual_features visp_vs visp_vision)
IF(VISP_FOUND)
  INCLUDE(${VISP_USE_FILE})
ENDIF(VISP_FOUND)
LINK_DIRECTORIES(${VISP_LINK_DIRECTORIES})

set(CMAKE_MODULE_PATH 
# 	/usr/local/lib/OGRE/cmake
	${CMAKE_SOURCE_DIR}/CMakeModules
  )

  ADD_DEFINITIONS(  -msse -msse2 -ffast-math -funroll-loops -mfpmath=sse -fomit-frame-pointer) #-Wall


#FIND_PACKAGE(OpenGL REQUIRED)
#FIND_PACKAGE(OGRE REQUIRED)
#FIND_PACKAGE(OIS REQUIRED)
 
#LINK_DIRECTORIES(${OGRE_LINK_DIRECTORIES})
 
INCLUDE_DIRECTORIES(
# 	${OpenGL_INCLUDE_DIR}
# 	${OGRE_INCLUDE_DIRS}
# 	${OIS_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/src 
	${VISP_INCLUDE_DIR}
  )

#------------------------------------------------------------
# OpenMP
#------------------------------------------------------------

OPTION(USE_OMP "Use OpenMP" OFF)

OPTION(USE_ICC "Use Intel Compiler" OFF)

IF(USE_OMP)
    IF(UNIX)
        IF(USE_ICC)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dnoomp -openmp")
            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Dnoomp -openmp")
        ELSE(USE_ICC)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
        ENDIF(USE_ICC)
    ELSE(UNIX)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp")
    ENDIF(UNIX)
    ADD_DEFINITIONS(-D USE_OMP)
ENDIF(USE_OMP)

#----------------------------------------------------------------------
# Concatenate intern and third party include dirs
#----------------------------------------------------------------------
INCLUDE_DIRECTORIES(${PhOVS_INTERN_INCLUDE_DIR} ${PhOVS_EXTERN_INCLUDE_DIR})

#----------------------------------------------------------------------
# Add library path
#----------------------------------------------------------------------
LINK_DIRECTORIES(${PhOVS_SOURCE_DIR}/lib)


#----------------------------------------------------------------------
# add html-doc target for doxygen documentation building
#----------------------------------------------------------------------
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)  
  ## we need latex for doxygen because of the formulas
  FIND_PACKAGE(LATEX)
  IF    (NOT LATEX_COMPILER)
    MESSAGE(STATUS "latex command LATEX_COMPILER not found but usually required. You will probably get warnings and user inetraction on doxy run.")
  ENDIF (NOT LATEX_COMPILER)
  IF    (NOT MAKEINDEX_COMPILER)
    MESSAGE(STATUS "makeindex command MAKEINDEX_COMPILER not found but usually required.")
  ENDIF (NOT MAKEINDEX_COMPILER)
  IF    (NOT DVIPS_CONVERTER)
    MESSAGE(STATUS "dvips command DVIPS_CONVERTER not found but usually required.")
  ENDIF (NOT DVIPS_CONVERTER)
   
  CONFIGURE_FILE(${PhOVS_SOURCE_DIR}/doc/config-doxygen.in
    ${PhOVS_SOURCE_DIR}/doc/config-doxygen
    @ONLY )
  ADD_CUSTOM_TARGET(html-doc ${DOXYGEN_EXECUTABLE} ${PhOVS_SOURCE_DIR}/doc/config-doxygen) 
ENDIF(DOXYGEN_FOUND)

#----------------------------------------------------------------------
# Add definitions
#----------------------------------------------------------------------
ADD_DEFINITIONS(${PhOVS_DEFS})

#----------------------------------------------------------------------
# Propagation in sub dirs
#----------------------------------------------------------------------
SUBDIRS(src)
#IF (EXISTS PhotometricOmniVirtualVS)
  SUBDIRS(PhotometricOmniVirtualVS)
#ENDIF(EXISTS PhotometricOmniVirtualVS)
#IF (EXISTS PhotometricOmniTrack)
SUBDIRS(PhotometricOmniTrack)
#ENDIF(EXISTS PhotometricOmniTrack)
#IF (EXISTS PhotometricOmniRegistration)
SUBDIRS(PhotometricOmniRegistration)
#ENDIF(EXISTS PhotometricOmniRegistration)
#IF (EXISTS PhotometricOmniRegistration_scaleSpace)
SUBDIRS(PhotometricOmniRegistration_scaleSpace)
#ENDIF(EXISTS PhotometricOmniRegistration_scaleSpace)
#IF (EXISTS ShallowDepthOfFieldPhotometricVS)
SUBDIRS(ShallowDepthOfFieldPhotometricVS)
#ENDIF(EXISTS ShallowDepthOfFieldPhotometricVS)
#IF (EXISTS DirectionalBlurPhotometricVS)
SUBDIRS(DirectionalBlurPhotometricVS)
#ENDIF(EXISTS DirectionalBlurPhotometricVS)
#IF (EXISTS "./Defocus_based_DVS")
SUBDIRS("./Defocus_based_DVS")
#ENDIF(EXISTS "./Defocus_based_DVS")
SUBDIRS(PhotometricVS_preCond)

