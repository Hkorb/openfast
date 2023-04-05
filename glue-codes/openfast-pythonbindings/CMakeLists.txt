cmake_minimum_required(VERSION 3.4)
project(${SKBUILD_PROJECT_NAME}  LANGUAGES CXX)

find_package(pybind11 CONFIG REQUIRED)
find_package(MPI REQUIRED)

pybind11_add_module(openfast_pythonbinding MODULE src/openfast.cpp)

include_directories(${CMAKE_SOURCE_DIR}/glue-codes/openfast-cpp/src/)

include_directories(${CMAKE_SOURCE_DIR}/modules/openfast-library/src/)
include_directories(${CMAKE_BINARY_DIR}/modules/openfoam/)
include_directories(${CMAKE_BINARY_DIR}/modules/supercontroller/)


target_include_directories(openfast_pythonbinding PRIVATE pybind11_INCLUDE_DIR ${MPI_INCLUDE_PATH})


target_link_libraries(openfast_pythonbinding PRIVATE
        openfastcpplib
        ${MPI_LIBRARIES})

set_target_properties(openfast_pythonbinding PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
set_target_properties(openfast_pythonbinding PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")

set_target_properties(  openfast_pythonbinding PROPERTIES 
                            OUTPUT_NAME "bindings")
set_property(TARGET openfast_pythonbinding PROPERTY LINKER_LANGUAGE CXX)
install(TARGETS openfast_pythonbinding EXPORT OpenFASTTargets LIBRARY DESTINATION .)
foreach(_cmake_expected_target IN ITEMS nwtclibs inflowwind_driver ifwlib ifw_c_binding unsteadyaero_driver aerodyn_driver adilib aerodynlib fvwlib uaaerolib afinfolib aeroacoustics aerodyn_inflow_c_binding aerodyn14lib servodynlib servodyn_driver elastodynlib beamdynlib subdynlib subdyn_driver hydrodynlib hydrodyn_driver hydrodyn_c_binding orcaflexlib orca_driver extptfm_mckflib openfoamtypeslib foamfastlib sctypeslib scfastlib scdataextypeslib scdataexlib sclib openfastlib openfast_prelib openfast_postlib versioninfolib feamlib feam_driver moordynlib icedynlib icefloelib maplib mapcpplib wdlib awaelib openfastcpplib)
        install(TARGETS ${_cmake_expected_target} LIBRARY DESTINATION ${Python_ROOT_DIR}/lib)
endforeach()

