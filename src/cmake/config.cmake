#
#  $Id$
#
#  Copyright (C) 2007   Troy D. Straszheim  <troy@icecube.umd.edu>
#  and the IceCube Collaboration <http://www.icecube.wisc.edu>
#
#  This file is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>
#
include(numeric_version)
colormsg("")
colormsg(_HIBLUE_ "Configuration starting")
colormsg("")

## set some cmake flags - see cmake_variables(7)
set(BUILD_SHARED_LIBS SHARED)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
## this is our flag, but in the above style
set(DOXYGEN_OUTPUT_PATH ${CMAKE_BINARY_DIR}/docs/doxygen)
set(CRA_SRC ${CMAKE_SOURCE_DIR})
set(CRA_BUILD ${CMAKE_BINARY_DIR})

#
#  Check build sanity
#
## ensure that we're not cmake'ing under an env-shell from another workspace

## default the option SYSTEM_PACKAGES to ON
option(SYSTEM_PACKAGES "Use tools provided by the operating system" ON)

#
# Find CVMFS software
# (designed for py2-v2 port-less software)
#
if($ENV{SROOT} MATCHES "^/cvmfs/icecube")
  if($ENV{SROOTBASE} MATCHES "py2-v1$")
    set(SYSTEM_PACKAGES FALSE)
  else($ENV{SROOTBASE} MATCHES "py2-v1$")
    set(SYSTEM_PACKAGES TRUE)
  endif()
  set(CMAKE_PREFIX_PATH $ENV{SROOT})
  set(USE_CVMFS TRUE CACHE BOOL "Are we using CVMFS?")
  set(CVMFS_SROOTBASE "$ENV{SROOTBASE}" CACHE STRING "CVMFS toolset path" )
else($ENV{SROOT} MATCHES "^/cvmfs/icecube")
  set(CVMFS_SROOTBASE /cvmfs/icecube.opensciencegrid.org/py3-v4.2.0 CACHE STRING "default CVMFS toolset path" )
endif()




#
# Create various info/debug files
#
execute_process(COMMAND /usr/bin/env
  OUTPUT_FILE ${NOTES_DIR}/env.txt)

execute_process(COMMAND uname -a
  OUTPUT_FILE ${NOTES_DIR}/uname.txt)

execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version
  OUTPUT_FILE ${NOTES_DIR}/compiler-version.txt)

#
#  COMPILER_ID_TAG
#
execute_process(COMMAND ${CMAKE_CXX_COMPILER} -v ERROR_VARIABLE COMPILER_ID_TAG)
set(COMPILER_ID_TAG "REGEXPS IN CMAKE SUCK\n${COMPILER_ID_TAG}")
string(REGEX REPLACE "^.*(g(cc|[+][+])|clang|Apple LLVM)[ -][Vv]ers(ion|ión|io|ão) ([^ \n]+).*"
                     "\\1-\\4" COMPILER_ID_TAG ${COMPILER_ID_TAG})

#
# Get system info
#
set(OSTYPE ${CMAKE_SYSTEM_NAME})
boost_report_value(OSTYPE)

set(OSVERSION ${CMAKE_SYSTEM_VERSION})
boost_report_value(OSVERSION)

set(ARCH ${CMAKE_SYSTEM_PROCESSOR})
boost_report_value(ARCH)

set(BUILDNAME "${OSTYPE}-${OSVERSION}/${ARCH}/${COMPILER_ID_TAG}" CACHE INTERNAL "buildname")
boost_report_value(BUILDNAME)

set(TOOLSET "${COMPILER_ID_TAG}/${ARCH}/${CMAKE_BUILD_TYPE}" CACHE INTERNAL "toolset")
boost_report_value(TOOLSET)

execute_process(COMMAND hostname
  COMMAND tr -d \\n
  OUTPUT_VARIABLE HOSTNAME)
boost_report_value(HOSTNAME)
set(SITE ${HOSTNAME})

#
# Show cmake path and version
#
boost_report_pretty("CMake path" CMAKE_COMMAND)
if(NOT CMAKE_VERSION)
  set(CMAKE_VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION})
endif(NOT CMAKE_VERSION)
math(EXPR CMAKE_VERSION_INT "${CMAKE_MAJOR_VERSION} * 10000 + ${CMAKE_MINOR_VERSION} * 100 + ${CMAKE_PATCH_VERSION}")
boost_report_pretty("CMake version" CMAKE_VERSION)
