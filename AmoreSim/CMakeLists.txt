#----------------------------------------------------------------------------
# Project Name: AmoreSim
# Description: Simulation library & executable for AmoreSim package
# Author: BaseHardware(basehw); 2018/05/31
#----------------------------------------------------------------------------
# Setup the project and require CMake version above 3.5
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(AmoreSim)
if(EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    message("Build information file for conan framework has been detected. Switching to conan mode.")
    include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    conan_basic_setup()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

#----------------------------------------------------------------------------
# Add variables for the git revision
#----------------------------------------------------------------------------
# Get the current working branch
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE AmoreSim_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

# Get the latest abbreviated commit hash of the working branch
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE AmoreSim_GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

#----------------------------------------------------------------------------
# Find ROOT and setup it
#----------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Find Geant4 and setup it
#----------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
    add_definitions(-DG4VIS_USE)
else()
    find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})
#For backward compatibility
set(CMAKE_CXX_STANDARD 17)

#----------------------------------------------------------------------------
# Check environment variables for run scripts and macros
#----------------------------------------------------------------------------
if(NOT DEFINED ENV{MUON_PATH})
    message("Environment variable MUON_PATH is not defined. Output path will be set to the build directory.")
endif()

if(NOT DEFINED ENV{NEUTRON_PATH})
    message("Environment variable NEUTRON_PATH is not defined. Output path will be set to the build directory.")
endif()

if(NOT DEFINED ENV{SIMOUT_PATH})
    message("Environment variable SIMOUT_PATH is not defined. Output path will be set to the build directory.")
endif()


#----------------------------------------------------------------------------
# Add variables for run scripts and macros
#----------------------------------------------------------------------------
set(AMORESIM_YEMI_MUON_FLUX_NAME "atHandeok_BIGNUM" CACHE STRING "Muon flux set name to use for muon simulation at YEMI (AMoRE-II/200)")
set(AMORESIM_Y2L_MUON_FLUX_NAME "atYangYangNew_BIGNUM" CACHE STRING "Muon flux set name to use for muon simulation at Y2L (AMoRE-Pilot/I)")
set(AMORESIM_NEUTRON_FLUX_NAME "th_neut_bin1_uniform_1000" CACHE STRING "Neutron flux set name to use for this simulation")

set(AMORESIM_VETO_CONF "WithBeam" CACHE STRING "Veto geometry configuration for this simulation")
set(AMORESIM_JOB_NAME "scratch" CACHE STRING "Representative name for this compilation")
if("${AMORESIM_JOB_NAME}" STREQUAL "")
    set(AMORESIM_JOB_NAME "scratch" FORCE)
endif()

set(AMORESIM_NEUTRON_MODE ON CACHE BOOL "Enable neutron mode for AMoRE-II, Pilot simulation")

if(${AMORESIM_NEUTRON_MODE})
    set(AMORESIM_NMODE "true")
else()
    set(AMORESIM_NMODE "false")
endif()

set(AMORESIM_DECAY_SOURCE_NAME "U238 Th232 K40 U235 Pb210 Na22 I125 I126 Te121 Te121m Te123m Te125m Te127m H3" CACHE STRING "Target isotopes list for decay simulation (Seperated by space)")
set(AMORESIM_DECAY_RATE "4.9E-18 1.6E-18 1.8E-17 3.1E-17 1E-9 8.5E-9 1.4E-7 6.2E-7 4.2E-7 4.9E-8 6.7E-8 1.4E-7 7.6E-8 1.8E-9" CACHE STRING "Half life list for target isotopes (Seperated by space, Unit in Hz)")
set(AMORESIM_NCORE "100" CACHE STRING "Total number of jobs for simulation of each isotope in the list")
set(AMORESIM_SOURCE_CATEGORY "internal" CACHE STRING "Simulation category for decay simulation (internal or external-pmt)")
set(AMORESIM_DECAY_EVTNUM "1000" CACHE STRING "The number of events each job will simulate")

#----------------------------------------------------------------------------
# Set advanced properties to variables
#----------------------------------------------------------------------------
mark_as_advanced(AMORESIM_Y2L_MUON_FLUX_NAME)
mark_as_advanced(AMORESIM_YEMI_MUON_FLUX_NAME)
mark_as_advanced(AMORESIM_NEUTRON_FLUX_NAME)
mark_as_advanced(AMORESIM_VETO_CONF)
mark_as_advanced(AMORESIM_NEUTRON_MODE)
mark_as_advanced(AMORESIM_DECAY_SOURCE_NAME)
mark_as_advanced(AMORESIM_DECAY_RATE)
mark_as_advanced(AMORESIM_NCORE)
mark_as_advanced(AMORESIM_SOURCE_CATEGORY)
mark_as_advanced(AMORESIM_DECAY_EVTNUM)

#----------------------------------------------------------------------------
# Set variables for run scripts and macros
#----------------------------------------------------------------------------
set(MUON_Y2L_PATH $ENV{MUON_PATH}/${AMORESIM_Y2L_MUON_FLUX_NAME})
set(MUON_YEMI_PATH $ENV{MUON_PATH}/${AMORESIM_YEMI_MUON_FLUX_NAME})
set(NEUTRON_PATH $ENV{NEUTRON_PATH})

if(NOT DEFINED ENV{SIMOUT_PATH})
    set(SIMOUT_PATH ${PROJECT_BINARY_DIR}/output)
else ()
    set(SIMOUT_PATH $ENV{SIMOUT_PATH})
endif()

set(AMORESIM_WORK_DIR ${CMAKE_BINARY_DIR})
set(RUN_SCRIPTS
		olaf_jobsubmit.sh
    run_II_muon.sh
    run_II_neut.sh
    run_II_decay.sh
    run_I_muon.sh
    run_I_neut.sh
    run_Pilot_muon.sh
    run_Pilot_neut.sh
    run_Pilot_decay.sh
    run_II_muon_GV.sh
    run_II_neut_GV.sh
    run_I_muon_GV.sh
    run_I_neut_GV.sh
    run_Pilot_muon_GV.sh
    run_Pilot_neut_GV.sh
    run_Pilot_decay_GV.sh
    rerun.sh
    )
set(SIM_MACROS
    II_muonbckg.mac
    II_neutbckg.mac
    II_dc_internal.mac
    I_muonbckg.mac
    I_neutbckg.mac
    Pilot_muonbckg.mac
    Pilot_neutbckg.mac
    Pilot_dc_external-pmt.mac
    Pilot_dc_internal.mac
    geom_validation.mac
    )
set(SESSION_MACROS
    init_vis.mac
    vis.mac
    gui.mac
    )
set(MCObjs_ENVIRONMENT_SCRIPTS environment.sh)

#----------------------------------------------------------------------------
# Configure the scripts and macros
#----------------------------------------------------------------------------
foreach(_SCRIPT ${RUN_SCRIPTS})
    configure_file (
        ${PROJECT_SOURCE_DIR}/run/${_SCRIPT}.in
        ${PROJECT_BINARY_DIR}/run/${_SCRIPT}
        @ONLY
        )
endforeach()
foreach(_SCRIPT ${MCObjs_ENVIRONMENT_SCRIPTS})
    configure_file (
        ${PROJECT_SOURCE_DIR}/${_SCRIPT}.in
        ${PROJECT_BINARY_DIR}/${_SCRIPT}
        @ONLY
        )
endforeach()

foreach(_MACRO ${SIM_MACROS})
    configure_file (
        ${PROJECT_SOURCE_DIR}/mac/${_MACRO}.in
        ${PROJECT_BINARY_DIR}/mac/${_MACRO}
        @ONLY
        )
endforeach()
foreach(_MACRO ${SESSION_MACROS})
    configure_file (
        ${PROJECT_SOURCE_DIR}/${_MACRO}
        ${PROJECT_BINARY_DIR}/${_MACRO}
        @ONLY
        )
endforeach()

#----------------------------------------------------------------------------
# Copy the data files to build directory
#----------------------------------------------------------------------------
file(COPY ${PROJECT_SOURCE_DIR}/mac ${PROJECT_SOURCE_DIR}/run ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR}
    PATTERN *.in EXCLUDE)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output/${AMORESIM_JOB_NAME})
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output/${AMORESIM_JOB_NAME}/log)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output/${AMORESIM_JOB_NAME}/mac)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output/${AMORESIM_JOB_NAME}/root)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
file(GLOB AmoreSim_EXEC_SOURCE ${PROJECT_SOURCE_DIR}/test/amoresim.cc)
file(GLOB AmoreSim_GitRev_SOURCE ${PROJECT_SOURCE_DIR}/src/AmoreSimGitRevision.cc)
file(GLOB AmoreSim_LIB_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB AmoreSim_LIB_HEADERS ${PROJECT_SOURCE_DIR}/AmoreSim/*.hh)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Add definitions for git commit hash
#----------------------------------------------------------------------------
set_source_files_properties(${AmoreSim_EXEC_SOURCE} PROPERTIES
    COMPILE_DEFINITIONS "AmoreSim_GIT_BRANCH=${AmoreSim_GIT_BRANCH};AmoreSim_GIT_COMMIT_HASH=${AmoreSim_GIT_COMMIT_HASH}" )
set_source_files_properties(${AmoreSim_GitRev_SOURCE} PROPERTIES
    COMPILE_DEFINITIONS "AmoreSim_GIT_BRANCH=${AmoreSim_GIT_BRANCH};AmoreSim_GIT_COMMIT_HASH=${AmoreSim_GIT_COMMIT_HASH}" )

#----------------------------------------------------------------------------
# Add libraries and executables, and link it to the ROOT and Geant4 framework library
#----------------------------------------------------------------------------
add_library(AmoreSimL SHARED ${AmoreSim_LIB_SOURCES})
target_link_libraries(AmoreSimL CupSimL ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
add_executable(amoresim ${AmoreSim_EXEC_SOURCE})
target_link_libraries(amoresim AmoreSimL CupSimL ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Check dependencies for this project and set include directories and libraries
#----------------------------------------------------------------------------
if(DEFINED CupSim)
    include_directories(${CupSim_INCLUDE_DIR})
    target_include_directories(amoresim PUBLIC ${CupSim_INCLUDE_DIR})
    target_include_directories(AmoreSimL PUBLIC ${CupSim_INCLUDE_DIR})
    target_link_libraries(amoresim CupSimL)
elseif(DEFINED CONAN_INCLUDE_DIRS_CUPSIM)
    include_directories(${CONAN_INCLUDE_DIRS_CUPSIM})
    target_include_directories(AmoreSimL PUBLIC ${CONAN_INCLUDE_DIRS_CUPSIM})
    target_include_directories(amoresim PUBLIC ${CONAN_INCLUDE_DIRS_CUPSIM})
    target_link_libraries(amoresim ${CONAN_LIBS_CUPSIM})
else()
    message(SEND_ERROR "This project cannot be built without CupSim. Please check about this.")
endif()

if(DEFINED MCObjs)
    include_directories(${MCObjs_INCLUDE_DIR})
    target_include_directories(amoresim PUBLIC ${MCObjs_INCLUDE_DIR})
    target_include_directories(AmoreSimL PUBLIC ${MCObjs_INCLUDE_DIR})
    target_link_libraries(amoresim MCObjs)
elseif(DEFINED CONAN_INCLUDE_DIRS_MCOBJS)
    include_directories(${CONAN_INCLUDE_DIRS_MCOBJS})
    target_include_directories(AmoreSimL PUBLIC ${CONAN_INCLUDE_DIRS_MCOBJS})
    target_include_directories(amoresim PUBLIC ${CONAN_INCLUDE_DIRS_MCOBJS})
    target_link_libraries(amoresim ${CONAN_LIBS_MCOBJS})
else()
    message(SEND_ERROR "This project cannot be built without MCObjs. Please check about this.")
endif()

#----------------------------------------------------------------------------
# Setup global and local include directories for this project
#----------------------------------------------------------------------------
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${Geant4_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIRS})
target_include_directories(amoresim PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(amoresim PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(amoresim PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories(AmoreSimL PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(AmoreSimL PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(AmoreSimL PUBLIC ${ROOT_INCLUDE_DIRS})
if(CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(AmoreSimL PUBLIC -fdiagnostics-color=always)
    target_compile_options(amoresim PUBLIC -fdiagnostics-color=always)
endif()

#----------------------------------------------------------------------------
# Expose this public includes and library to other subprojects through cache
# variable.
#----------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIR ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
set(${PROJECT_NAME}_LIBRARY ${PROJECT_BINARY_DIR}/lib${PROJECT_NAME}.so
    CACHE INTERNAL "${PROJECT_NAME}: Path to library" FORCE)
set(${PROJECT_NAME} "Exist!"
    CACHE INTERNAL "${PROJECT_NAME}: Variable for checking dependencies" FORCE)
