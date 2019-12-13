# CMake project file for FEST-3D

CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
PROJECT(FEST3D Fortran)
enable_language (C Fortran)
enable_testing ()

# Set the FOO version
SET(VERSION 1.1.1)

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" ON)
#OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

function(LinkFortranLibraries Target)
  target_link_libraries(Target ${ARGN})
  foreach(Lib IN LISTS ARGN)
    target_include_directories(Simulation PUBLIC $<TARGET_PROPERTY:${Lib},Fortran_MODULE_DIRECTORY>)
  endforeach()
endfunction()

# Define the executable name
SET(EXE FEST3D)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(MOD ${CMAKE_SOURCE_DIR}/mod)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(TST ${CMAKE_SOURCE_DIR}/tests)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})

# The source for the FEST3D binary and have it placed in bin folder
Add_subdirectory(${SRC} ${BIN})

# Add ctest function
INCLUDE(${CMAKE_MODULE_PATH}/CTestFunctions.cmake)
Add_subdirectory(${TST})

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall")
# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
