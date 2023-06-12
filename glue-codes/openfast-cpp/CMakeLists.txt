#
# Copyright 2016 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  cmake_policy(SET CMP0074 NEW)
endif()

enable_language(CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# find_package(MPI REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(yaml-cpp REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/modules/openfast-library/src/)
include_directories(${CMAKE_BINARY_DIR}/modules/openfoam/)
include_directories(${CMAKE_BINARY_DIR}/modules/supercontroller/)
# include_directories(${MPI_INCLUDE_PATH})

add_library(openfastcpplib
  src/OpenFAST.cpp 
  src/SC.cpp
  $<TARGET_OBJECTS:openfastlib_common_obj>
  $<TARGET_OBJECTS:openfast_prelib_obj>
  $<TARGET_OBJECTS:maplib_obj>
  $<TARGET_OBJECTS:mapcpplib_obj>
  $<TARGET_OBJECTS:nwtcbaselib_obj> 
  $<TARGET_OBJECTS:nwtcsyslib_obj>
  $<TARGET_OBJECTS:nwtclibs_obj>
  $<TARGET_OBJECTS:ifwlib_obj>
  $<TARGET_OBJECTS:fvwlib_obj>
  $<TARGET_OBJECTS:uaaerolib_obj>
  $<TARGET_OBJECTS:afinfolib_obj>
  $<TARGET_OBJECTS:aeroacoustics_obj>
  $<TARGET_OBJECTS:aerodynlib_obj>
  $<TARGET_OBJECTS:aerodyn14lib_obj>
  $<TARGET_OBJECTS:servodynlib_obj>
  $<TARGET_OBJECTS:elastodynlib_obj>
  $<TARGET_OBJECTS:beamdynlib_obj>
  $<TARGET_OBJECTS:subdynlib_obj>
  $<TARGET_OBJECTS:hydrodynlib_obj>
  $<TARGET_OBJECTS:orcaflexlib_obj>
  $<TARGET_OBJECTS:extptfm_mckflib_obj>
  $<TARGET_OBJECTS:openfoamtypeslib_obj>
  $<TARGET_OBJECTS:foamfastlib_obj>
  $<TARGET_OBJECTS:scdataextypeslib_obj>
  $<TARGET_OBJECTS:scdataexlib_obj>
  $<TARGET_OBJECTS:feamlib_obj>
  $<TARGET_OBJECTS:moordynlib_obj>
  $<TARGET_OBJECTS:icedynlib_obj>
  $<TARGET_OBJECTS:icefloelib_obj>
  $<TARGET_OBJECTS:openfast_postlib_obj>
  $<TARGET_OBJECTS:versioninfolib_obj>
)

set_property(TARGET openfastcpplib PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(openfastcpplib
  openfastlib
  ${HDF5_C_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${LIBXML2_LIBRARIES}
  # ${MPI_LIBRARIES}
)
target_include_directories(openfastcpplib PUBLIC 
  ${HDF5_INCLUDES}
  ${HDF5_INCLUDE_DIR}
  ${ZLIB_INCLUDES}
  ${LIBXML2_INCLUDE_DIR}
  ${MPI_INCLUDE_PATH}
)
set_target_properties(openfastcpplib PROPERTIES PUBLIC_HEADER "src/OpenFAST.H;src/SC.h")

add_executable(openfastcpp src/FAST_Prog.cpp)
target_link_libraries(openfastcpp 
  openfastcpplib
  # ${MPI_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  openfastcpplib 
)

if (OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(openfastcpplib OpenMP::OpenMP_CXX)
    target_link_libraries(openfastcpp OpenMP::OpenMP_CXX)
  endif()
endif(OPENMP)

# if(MPI_COMPILE_FLAGS)
#   set_target_properties(openfastcpp PROPERTIES
#     COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
# endif(MPI_COMPILE_FLAGS)


set_property(TARGET openfastcpp PROPERTY LINKER_LANGUAGE CXX)

install(TARGETS openfastcpplib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)
