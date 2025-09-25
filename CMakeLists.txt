cmake_minimum_required(VERSION 3.22)

project(HeatMassTransfer LANGUAGES C CXX)

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
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

ifem_add_hdf5_test(
  TARGET
    HeatMassTransfer
  TEST_FILES
    Slab.hreg
)

ifem_add_vtf_test(
  TARGET
    HeatMassTransfer
  TEST_FILES
    Slab.vreg
)

ifem_add_regression_test(
  TARGET
    HeatMassTransfer
  TEST_FILES
    Slab.reg
)

# Unit tests
ifem_add_test_app(
  NAME
    HM
  SOURCES
    ${PROJECT_SOURCE_DIR}/Test/TestHMProperties.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/Test
  LIBRARIES
    CommonHM
)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(TARGET HeatMassTransfer)

# Installation
install(TARGETS HeatMassTransfer DESTINATION bin)
