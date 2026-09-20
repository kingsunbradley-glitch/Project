cmake_minimum_required(VERSION 3.20)
project(SuperconductingSolenoidGeant4 VERSION 1.0 LANGUAGES CXX)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(SOLENOID_BUILD_TESTS "Build the unit tests" ON)

find_package(Geant4 11.4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

set(SOLENOID_HEADERS
    include/SimulationConfig.hh
    include/FiniteSolenoidField.hh
    include/ReactionKinematics.hh
    include/EventInformation.hh
    include/DetectorConstruction.hh
    include/PrimaryGeneratorAction.hh
    include/ActionInitialization.hh
    include/RunAction.hh
    include/EventAction.hh
    include/SteppingAction.hh)

set(SOLENOID_SOURCES
    src/SimulationConfig.cc
    src/FiniteSolenoidField.cc
    src/ReactionKinematics.cc
    src/DetectorConstruction.cc
    src/PrimaryGeneratorAction.cc
    src/ActionInitialization.cc
    src/RunAction.cc
    src/EventAction.cc
    src/SteppingAction.cc)

add_library(solenoid_core ${SOLENOID_HEADERS} ${SOLENOID_SOURCES})
target_include_directories(solenoid_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(solenoid_core PUBLIC ${Geant4_LIBRARIES})
target_compile_options(solenoid_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(solenoid_transport src/main.cc)
target_link_libraries(solenoid_transport PRIVATE solenoid_core)

file(COPY macros DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(SOLENOID_BUILD_TESTS)
  enable_testing()
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_executable(solenoid_unit_tests tests/unit_tests.cc)
  target_link_libraries(solenoid_unit_tests PRIVATE solenoid_core)
  add_test(NAME solenoid_unit_tests COMMAND solenoid_unit_tests)
  add_test(NAME solenoid_integration_checks
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration_checks.py
                   --exe $<TARGET_FILE:solenoid_transport>)
  set_tests_properties(solenoid_integration_checks PROPERTIES TIMEOUT 300)
endif()
