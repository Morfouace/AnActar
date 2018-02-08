cmake_minimum_required (VERSION 2.8)
include(CheckCXXCompilerFlag)
project(actar_reco)
set(CMAKE_BUILD_TYPE Release)

# Setting the policy to match Cmake version
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})

# include the actar standard CMake preamble
include("ressources/CMake/ACTAR_CMake_Preamble.cmake")

#find_package(PCL 1.2 REQUIRED)
#find_package(Boost)

#if(PCL_FOUND)
#  MESSAGE("-- PCL Support Found")
#  MESSAGE(${PCL_DEFINITIONS})
#  MESSAGE(${PCL_INCLUDE_DIRS})
#endif()

#if(VTK_FOUND)
#  MESSAGE("-- VTK Support Found")
#endif()

set(CMAKE_BINARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set(CMAKE_INCLUDE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/include )

#set(MFM_DIR /home/morfouace/Physics/Actar/Analysis/actar_analysis/libMFM)
#set(CCFG_DIR /home/morfouace/Physics/Actar/Analysis/CompoundConfig/src/CCfg)

#set(INCLUDE_DIRECTORIES
#${ROOT_INCLUDE_DIR}
#${Boost_INCLUDE_DIRS}
#${BASE_INCLUDE_DIRECTORIES}
#${PCL_INCLUDE_DIRS}
#${MFM_DIR}
#${CCFG_DIR}
#)
#include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
#${ROOT_LIBRARY_DIR}
#${Boost_LIBRARY_DIRS}
#${PCL_LIBRARY_DIRS}
#${MFM_DIR}
#${CCFG_DIR}
${CMAKE_SOURCE_DIR}/core
${CMAKE_SOURCE_DIR}/reco
${CMAKE_SOURCE_DIR}/utils
)
link_directories(${LINK_DIRECTORIES})

#Set(DEPENDENCIES
#    Base
#    ${Boost_LOG_LIBRARY}
#    ${Boost_LOG_SETUP_LIBRARY}
#    ${Boost_THREAD_LIBRARY}
#    ${Boost_DATE_TIME_LIBRARY}
#    ${Boost_FILESYSTEM_LIBRARY}
#    ${Boost_SYSTEM_LIBRARY}
#    ${PCL_LIBRARIES}
#)
    

add_subdirectory (core)
add_subdirectory (reco)
add_subdirectory (utils)


export(PACKAGE ACTARLib)
file(RELATIVE_PATH REL_INCLUDE_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" "${CMAKE_INCLUDE_OUTPUT_DIRECTORY}")
set(CONF_INCLUDE_DIRS "./" "./")
configure_file(ressources/CMake/ACTARLibConfig.cmake.in "${PROJECT_BINARY_DIR}/ACTARLibConfig.cmake" @ONLY)









