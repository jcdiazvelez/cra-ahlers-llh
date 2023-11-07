#
#  $Id$
#
#  Copyright (C) 2007-9  Troy D. Straszheim  <troy@icecube.umd.edu>
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
colormsg("")
colormsg(_HIBLUE_ "Configuring tools...")
colormsg("")

add_custom_target(install_tool_libs)

include(tooldef)

#
#  By default, use /usr/share/fizzicks/cmake as I3_SITE_CMAKE_DIR
#
if (NOT IS_DIRECTORY $ENV{I3_SITE_CMAKE_DIR})
  set (I3_SITE_CMAKE_DIR "/usr/share/fizzicks/cmake"
    CACHE PATH "Path to site-specific cmake files")
  message(STATUS "Using default site cmake dir of ${I3_SITE_CMAKE_DIR}")
else()
  set (I3_SITE_CMAKE_DIR $ENV{I3_SITE_CMAKE_DIR}
    CACHE PATH "Path to site-specific cmake files")
  message(STATUS "Using user-configured I3_SITE_CMAKE_DIR=${I3_SITE_CMAKE_DIR}")
endif()

