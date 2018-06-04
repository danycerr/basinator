# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.
# cmake -DCGAL_DIR=/u/archive/agip/cerroni/software/cgal/CGAL-4.12-beta2/ 
#       -DBOOST_ROOT=$mkBoostPrefix  
#       -DCMAKE_BUILD_TYPE=Release ../


project( generator )


cmake_minimum_required(VERSION 2.8.11)
include_directories(/u/archive/agip/cerroni/software/cgal/include)

set(CMAKE_CXX_STANDARD 14)
add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)
# CGAL and its components
find_package(CGAL)


if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package


# Creating entries for all C++ files with "main" routine
# ##########################################################

# include( CGAL_CreateSingleSourceCGALProgram )

# create_single_source_cgal_program( "main.cpp" )
# create_single_source_cgal_program( "main.cpp" CXX_FEATURES cxx_range_for cxx_auto_type )


file(GLOB SOURCES "*.cpp")
add_executable( generator  ${SOURCES} )
target_include_directories(generator PRIVATE ./)
# Link the executable to CGAL and third-party libraries
target_link_libraries(generator   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

