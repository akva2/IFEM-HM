project(HeatMassTransfer)

cmake_minimum_required(VERSION 3.0)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

# Required packages
IF (NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
ENDIF(NOT IFEM_CONFIGURED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR})

set(HM_SOURCES HeatTransfer.C
               HeatTransfer.h
               HMArgs.C
               HMArgs.h
               HMProperties.C
               HMProperties.h
               MassTransfer.C
               MassTransfer.h
               SIMHeatTransfer.h
               SIMHM.h
               SIMMassTransfer.h)

add_library(CommonHM ${HM_SOURCES})

add_executable(HeatMassTransfer main_HM.C)
target_link_libraries(HeatMassTransfer CommonHM IFEMAppCommon ${IFEM_LIBRARIES})
add_executable(HeatTransfer main_HeatTransfer.C)
target_link_libraries(HeatTransfer CommonHM IFEMAppCommon ${IFEM_LIBRARIES})

list(APPEND CHECK_SOURCES ${HM_SOURCES} main_HM.C main_HeatTransfer.C)

# Installation
install(TARGETS HeatMassTransfer DESTINATION bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}")
add_doc_target(HeatMassTransfer HeatMassTransfer)

# Testing
enable_testing()
include(IFEMTesting)

if(MPI_FOUND)
endif()
if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  ifem_add_test(Slab.reg HeatMassTransfer)
  ifem_add_vtf_test(Slab.vreg HeatMassTransfer)
  ifem_add_hdf5_test(Slab.hreg HeatMassTransfer)
endif()
list(APPEND TEST_APPS HeatMassTransfer HeatTransfer)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  IFEM-THM CommonHM
                  IFEMAppCommon ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
