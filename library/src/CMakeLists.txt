# Optionally glob, but only for CMake 3.12 or later:
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/library/include/*.hpp")
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/library/src/*.cpp")
file(GLOB FORTRAN_SOURCE_LIST CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/library/src/fortran/*.f90")
list(APPEND SOURCE_LIST ${FORTRAN_SOURCE_LIST})

# include(CMakePrintHelpers)
# cmake_print_variables(SOURCE_LIST)

# Make an automatic library - will be static or dynamic based on user setting
add_library(scatter_library ${SOURCE_LIST} ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(scatter_library PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Linking boost regex and filesystem
target_link_libraries(scatter_library PRIVATE Boost::boost Boost::regex Boost::filesystem)

# Linking GSL
target_link_libraries(scatter_library PRIVATE GSL::gsl)

# All users of this library will need at least C++17
target_compile_features(scatter_library PUBLIC cxx_std_17)

# IDEs should put the headers in a nice place
source_group(TREE "${CMAKE_SOURCE_DIR}/library/include" PREFIX "Header Files" FILES ${HEADER_LIST})
