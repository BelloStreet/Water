add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")

# Keeps CACHE cleaner
mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

macro(package_add_test TESTNAME)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${ARGN})
    # link the Google test infrastructure, mocking library, and a default main fuction to
    # the test executable.  Remove g_test_main if writing your own main function.
    # include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main
      ${Boost_LIBRARIES} ${MKL_LIBRARIES} scatter_library)
      # ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} scatter_library)
    # gtest_discover_tests replaces gtest_add_tests,
    # see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
    gtest_discover_tests(${TESTNAME}
        # set a working directory so your project root so that you can find test data via paths relative to the project root
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

# Quadrature Test
package_add_test(testQuadrature testQuad/testQuad.cpp)

# FEMDVR Test
package_add_test(testFEMDVR testFEMDVR/testFEMDVR.cpp)

# Molden Parser Test
package_add_test(testParseMolden testParseMolden/testParseMolden.cpp)
configure_file(testParseMolden/H2O_scf.molden H2O_scf.molden COPYONLY)
configure_file(testParseMolden/h2o.molden h2o.molden COPYONLY)

# Simple Lebedev Test
# package_add_test(testLebedev testLebedev/testLebedev.cpp)

# MOPartialWaveRepresentation Test
package_add_test(testMOPartialWaveRepresentation testMOPartialWaveRepresentation/testMOPartialWaveRepresentation.cpp)

# Laplacian + Centrifugal Test
package_add_test(testToperator
  testToperator/testToperator.cpp)

# J operator Test
package_add_test(testJoperator testJoperator/testJoperator.cpp)

# K operator Test
package_add_test(testKoperator testKoperator/testKoperator.cpp)

# add_executable(mpi main.cpp mpi.cpp common.h)
# add_executable(testParseMolden testParseMolden.cpp)
# target_link_libraries(mpi PRIVATE MPI::MPI_CXX)
# target_link_libraries(testParseMolden PRIVATE Boost::filesystem Boos::regex)
# Copy the molden file
# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
# add_test(NAME testMoldenParser COMMAND testParseMolden) # Command can be a target

