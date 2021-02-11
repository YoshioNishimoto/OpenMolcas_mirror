#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2020, Ignacio Fdez. Galván                             *
#***********************************************************************

# Template for building an OpenMolcas program
# * Creates:
#   - An executable ${prog}.exe
#   - An object library ${prog} with all sources except main.*
#     (except if there are no such sources)
# * Supports:
#   - ${prog}_deplibs with program-specific libraries to link
#   - ${prog}_defs with program-specific compile definitions
#   - ${prog}_incs with program-specific include directories
# * Defines:
#   - ${prog}_src with the path for the source directory
#   - ${prog}_sources with all the source files (including main.*)

get_filename_component (prog ${CMAKE_CURRENT_SOURCE_DIR} NAME)

if (DEFINED sources)
  set_absolute_paths (sources ${CMAKE_CURRENT_SOURCE_DIR} ${sources})
else ()
  file (GLOB sources *.f *.f90 *.F *.F90 *.c)
endif ()

# This identifies any main.* file (there should be one), which is needed
# in order to be able to make a library out of the program sources. This
# is required for supermodules.
#-----------------------------------------------------------------------
set (main "")
foreach (fname ${sources})
  if (${fname} MATCHES "\/main\.(f|f90|F|F90|c)$")
    list (REMOVE_ITEM sources ${fname})
    set (main ${fname})
    break ()
  endif ()
endforeach ()

# There should be a main.* file, create an *.exe from it
#-------------------------------------------------------
add_Molcas_executable (${prog}.exe ${main})
# program-specific dependencies (supermodules)
if (DEFINED ${prog}_deplibs)
  target_link_libraries (${prog}.exe ${${prog}_deplibs})
endif ()
# program-specific compile definitions
if (DEFINED ${prog}_defs)
  target_compile_definitions (${prog}.exe PRIVATE "${${prog}_defs}")
endif ()
# program-specific include directories
list (APPEND ${prog}_incs ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (${prog}.exe PRIVATE "${${prog}_incs}")

# If there are any other files, create an object library
# for *.exe to link with
#--------------------------------------------------------
if (sources)
  add_Molcas_library (${prog} OBJECT ${sources})
  target_compile_definitions (${prog} PRIVATE $<TARGET_PROPERTY:${prog}.exe,COMPILE_DEFINITIONS>)
  target_include_directories (${prog} PRIVATE $<TARGET_PROPERTY:${prog}.exe,INCLUDE_DIRECTORIES>)
  target_link_libraries (${prog}.exe ${prog})
endif ()

# Info for the main project
list (APPEND sources ${main})
set (${prog}_src     ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
set (${prog}_sources ${sources}                  PARENT_SCOPE)
