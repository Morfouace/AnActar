cmake_minimum_required (VERSION 2.8) 
project (ACTARAnalysis)
set(CMAKE_BUILD_TYPE Release)  
# Setting the policy to match Cmake version
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})

#Finding ACTAR
find_package(ACTARLib)
#include(${ACTARLib_USE_FILE})

#Finding NPTool
find_package(NPLib)
include(${NPLib_USE_FILE})

#include the ACTAR standard CMake preamble 
include("${ACTARLIB}/ressources/CMake/ACTAR_CMake_Preamble.cmake")

#find_package(PCL 1.2 REQUIRED)

#set(MFM_DIR /home/morfouace/Physics/Actar/Analysis/actar_analysis/libMFM)
#set(CCFG_DIR /home/morfouace/Physics/Actar/Analysis/CompoundConfig/src/CCfg)

set(INCLUDE_DIRECTORIES
#${ROOT_INCLUDE_DIR}
#${MFM_DIR}
#${CCFG_DIR}
#${Boost_INCLUDE_DIRS}
#${BASE_INCLUDE_DIRECTORIES}
#${PCL_INCLUDE_DIRS}
)


include_directories(${INCLUDE_DIRECTORIES})

# allow link to not care for undefined reference
if(${CMAKE_CXX_COMPILER_ID} MATCHES ".*Clang.*")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments -fcolor-diagnostics -undefined dynamic_lookup")
endif()

add_executable(actaranalysis Analysis.cxx)
target_link_libraries(actaranalysis ${ROOT_LIBRARIES} -L${ACTARLIB}/lib -lACTARreco -lACTARcore -L${NPLIB}/lib -lNPCore -lNPPhysics)
#target_link_libraries(actaranalysis ${ROOT_LIBRARIES} ${PCL_LIBRARIES} -L${ACTARLIB}/lib)
