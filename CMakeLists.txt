# Top Level CMakeLists.txt

cmake_minimum_required( VERSION 3.15...3.17 )

# set( CMAKE_PROJECT_INCLUDE_BEFORE
#   "${CMAKE_CURRENT_LIST_DIR}/common-meta-info.in" )

# enable_language(Fortran)
# project(quantumScattering
#         DESCRIPTION ${project_description}
#         VERSION ${project_version}
#         HOMEPAGE_URL ${project_homepage}
#         LANGUAGES CXX)
project(quantumScattering
        LANGUAGES CXX Fortran)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

  # Optionally set things like CMAKE_CXX_STANDARD, CMAKE_POSITION_INDEPENDENT_CODE here

  # Let's ensure -std=c++xx instead of -std=g++xx
  set(CMAKE_CXX_EXTENSIONS OFF)

  # Let's nicely support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # Testing only available if this is the main app
  # Note this needs to be done in the main CMakeLists
  # since it calls enable_testing, which must be in the
  # main CMakeLists.
  # include(CTest)


  # Docs only available if this is the main app
  find_package(Doxygen)
      if(Doxygen_FOUND)
          add_subdirectory(docs)
      else()
          message(STATUS "Doxygen not found, not building docs")
  endif()
endif()

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# find_package(PETSc REQUIRED)

# find_package(PkgConfig REQUIRED)

#This is the logic to use MKL if MKL_ROOT is set in environment...
## check if MKL.
# if((MKL_ROOT) OR (DEFINED ENV{MKLROOT}))

#     find_package(MKL REQUIRED)
#     include_directories(${MKL_INCLUDE_DIRS})
#     set(LIBS ${LIBS} ${MKL_LIBRARIES})

# else()
    # look for FFTW3
    # find_package(FFTW REQUIRED)
    # include_directories(${FFTW_INCLUDES})
    # set(LIBS ${LIBS} ${FFTW_LIBRARIES})
# endif()

# find_package(LAPACK REQUIRED COMPONENTS lapacke)

# find_package(LAPACK REQUIRED)
find_package(MKL REQUIRED)
if(MKL_FOUND)
  include_directories(${MKL_INCLUDE_DIRS})
ENDIF()

find_package(MPI REQUIRED)
find_package(GSL REQUIRED)


# include_directories(${MPI_INCLUDE_PATH})

# pkg_check_modules(PETSC_PKG REQUIRED IMPORTED_TARGET PETSc)
# pkg_check_modules(SLEPC_PKG REQUIRED IMPORTED_TARGET SLEPc)

set( Boost_USE_STATIC_LIBS      FALSE )
set( Boost_USE_MULTITHREADED    TRUE )
set( Boost_USE_STATIC_RUNTIME   TRUE )
set( Boost_COMPILER             "-clang" )

find_package(Boost 1.50 REQUIRED COMPONENTS filesystem regex)

include_directories(library/include)

add_subdirectory(library)

option(PACKAGE_TESTS "Build the tests" ON)

if (PACKAGE_TESTS)
  enable_testing()
  include(GoogleTest)
  add_subdirectory(tests)
endif()
