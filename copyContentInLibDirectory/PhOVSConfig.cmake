#############################################################################
#
# PhOVSConfig.cmake, 2020/02/14, G. Caron, largely inspired from VISP
#
# $Id: VISPConfig.cmake.in,v 1.16 2008/11/07 09:56:04 fspindle Exp $
#
# Copyright (C) 1998-2006 Inria. All rights reserved.
#
# This software was developed at:
# IRISA/INRIA Rennes
# Projet Lagadic
# Campus Universitaire de Beaulieu
# 35042 Rennes Cedex
# http://www.irisa.fr/lagadic
#
# This file is part of the ViSP toolkit.
#
# This file may be distributed under the terms of the Q Public License
# as defined by Trolltech AS of Norway and appearing in the file
# LICENSE included in the packaging of this file.
#
# Licensees holding valid ViSP Professional Edition licenses may
# use this file in accordance with the ViSP Commercial License
# Agreement provided with the Software.
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# Contact visp@irisa.fr if any conditions of this licensing are
# not clear to you.
#
# Description:
# cmake PackageConfig file for VISP
#
# Typical usage in user project:
#   FIND_PACKAGE(VISP)
#   IF (VISP_FOUND)
#     INCLUDE(${VISP_USE_FILE})
#   ENDIF (VISP_FOUND)
#
# Authors:
# Fabien Spindler
#
#############################################################################

# Tells the user project ViSP library name
SET(PhOVS_LIBRARY "PhOVS")

# Tells the user project where to find PhOVS headers
SET(PhOVS_INCLUDE_DIR "${PhOVS_DIR}/../include" CACHE FILEPATH "Location of PhOVS includes")
#MESSAGE("PhOVS_INCLUDE_DIR: ${PhOVS_INCLUDE_DIR}")

SET(PhOVS_EXTERN_INCLUDE_DIR "/usr/include;/usr/local/include")

# Tells the user project where to find PhOVS library
SET(PhOVS_LINK_DIRECTORIES "${PhOVS_DIR}" CACHE FILEPATH "Location of PhOVS library")
#LIST(APPEND PhOVS_LINK_DIRECTORIES "/usr/X11/lib")

# Tells the user project ViSP library name
SET(PhOVS_LIBRARIES "PhOVS")

# Tells the user project where to find PhOVS build settings
SET(PhOVS_BUILD_SETTINGS_FILE "${PhOVS_DIR}/PhOVSBuildSettings.cmake")

# Tell the user project where to find PhOVS dependencies
INCLUDE("${PhOVS_DIR}/PhOVSLibraryDepends.cmake")

# where to find the USE file to be used by user project
SET(PhOVS_USE_FILE "${PhOVS_DIR}/PhOVSUse.cmake")

IF(BUILD_TEST_COVERAGE)
  # Add build options for test coverage. Currently coverage is only supported
  # on gcc compiler 
  # Because using -fprofile-arcs with shared lib can cause problems like:
  # hidden symbol `__bb_init_func', we add this option only for static 
  # library build
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ftest-coverage -fprofile-arcs")
ENDIF(BUILD_TEST_COVERAGE)

