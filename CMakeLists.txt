cmake_minimum_required(VERSION 3.20)

project(HeatMassTransfer LANGUAGES C CXX)

if(IFEM_AS_SUBMODULE)
  if(NOT TARGET IFEM)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  endif()
else()
  get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
  set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
  find_package(IFEM REQUIRED)
endif()

ifem_add_library(
  NAME
    CommonHM
  SOURCES
    HeatTransfer.C
    HMArgs.C
    HMProperties.C
    MassTransfer.C
  HEADERS
    HeatTransfer.h
    HMArgs.h
    HMProperties.h
    MassTransfer.h
    SIMHeatTransfer.h
    SIMHM.h
    SIMMassTransfer.h
  LIBRARIES
    IFEMAppCommon IFEM
)

ifem_add_application(
  NAME
    HeatMassTransfer
  SOURCES
    main_HM.C
  LIBRARIES
    CommonHM
)

ifem_add_application(
  NAME
    HeatTransfer
  SOURCES
    main_HeatTransfer.C
  LIBRARIES
    CommonHM
)

# Testing
enable_testing()

if(MPI_FOUND)
endif()
if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  ifem_add_test(Slab.reg HeatMassTransfer)
  ifem_add_vtf_test(Slab.vreg HeatMassTransfer)
  ifem_add_hdf5_test(Slab.hreg HeatMassTransfer)
endif()
list(APPEND TEST_APPS HeatMassTransfer HeatTransfer)

# Unit tests
ifem_add_test_app(SOURCES ${PROJECT_SOURCE_DIR}/Test/TestHMProperties.C
                  WORKDIR ${PROJECT_SOURCE_DIR}/Test
                  NAME IFEM-HM
                  LIBRARIES CommonHM)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(APP HeatMassTransfer)

# Installation
install(TARGETS HeatMassTransfer DESTINATION bin)
