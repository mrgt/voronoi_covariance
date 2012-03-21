# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project(Voronoi_covariance_3)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

IF(UNIX)
  IF(CMAKE_COMPILER_IS_GNUCC)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -frounding-math")
    ADD_DEFINITIONS(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
  ENDIF(CMAKE_COMPILER_IS_GNUCC)
ENDIF(UNIX)

find_package(CGAL COMPONENTS Qt4 required )
find_package(Boost COMPONENTS filesystem system program_options FATAL_ERROR)

set( QT_USE_QTXML    FALSE )
set( QT_USE_QTMAIN   FALSE )
set( QT_USE_QTSCRIPT  FALSE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)


include(${CGAL_USE_FILE})
include(${QT_USE_FILE})

include_directories (BEFORE .)

link_libraries(${QT_LIBRARIES}  ${Boost_SERIALIZATION_LIBRARY})
link_libraries(${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
link_libraries(${Boost_PROGRAM_OPTIONS_LIBRARY})
link_libraries(${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

add_executable  (vcm vcm.cpp)
add_executable  (sphere sphere.cpp)
