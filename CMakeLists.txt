# Copyright 2015 Dominic Deuber
#
# This file is part of VCFvifia.
#
# VCFvifia is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# VCFvifia is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with VCFvifia. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required (VERSION 2.8.8)
project (VCFvifia)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

enable_testing()

find_package(Boost 1.36 COMPONENTS program_options iostreams system filesystem unit_test_framework REQUIRED)

set(GNUPLOTINSTALLED false)
find_package(Gnuplot)
if (GNUPLOT_FOUND)
  set(GNUPLOTINSTALLED true)
endif (GNUPLOT_FOUND)
if (NOT GNUPLOT_FOUND)
  message("Please install gnuplot if you would like to use the plotting feature of VCFvifia")
endif (NOT GNUPLOT_FOUND)


SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -fPIC -Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -g -ggdb -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -fPIC -Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -O3")

SET(VCFLIB_PATH "${CMAKE_BINARY_DIR}/3rd-party/vcflib")
SET(GNUPLOT_IOSTREAM "${CMAKE_BINARY_DIR}/3rd-party/gnuplot-iostream")
include_directories(${Boost_INCLUDE_DIR} ${VCFLIB_PATH}  ${VCFLIB_PATH}/src ${VCFLIB_PATH}/tabixpp ${VCFLIB_PATH}/tabixpp/htslib ${GNUPLOT_IOSTREAM} src)
configure_file( config.hpp.in src/config.hpp)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  message("--Found Doxygen (use \"make doc\" to create documentation)")
  configure_file(${CMAKE_BINARY_DIR}/Doxyfile.in
    ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)
if(NOT DOXYGEN_FOUND)
  message("Please install Doxygen to create documentation")
endif(NOT DOXYGEN_FOUND)
ADD_CUSTOM_TARGET(
  libvcflib.a ALL
  COMMAND make libvcflib.a
  WORKING_DIRECTORY ${VCFLIB_PATH}
  COMMENT "Original libvcf makefile target"
)

ADD_CUSTOM_TARGET(
  libhts.a ALL
  COMMAND make lib-static
  WORKING_DIRECTORY ${VCFLIB_PATH}/tabixpp/htslib
  COMMENT "Original libhts makefile target"
)

add_executable(vcfvifia
  src/main.cpp
  src/BasicStatistic.cpp
  src/Validator.cpp
  src/Plotter.cpp
  src/Filter.cpp
  src/utils.cpp
)

TARGET_LINK_LIBRARIES(vcfvifia ${VCFLIB_PATH}/libvcflib.a ${VCFLIB_PATH}/tabixpp/htslib/libhts.a ${Boost_LIBRARIES} pthread z m)

add_subdirectory(tests)
