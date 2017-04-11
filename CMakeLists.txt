cmake_minimum_required(VERSION 2.6)
project(eema)
enable_language(CXX)

# Eigen support
option (USE_EIGEN  "Build library with Eigen support" ON)
set (EIGEN_INCLUDE_DIR "INCLUDE directory for Eigen library" third-party-libs/eigen3)
include_directories(${EIGEN_INCLUDE_DIR})

#===============================================================================
# Used to link multiple static libraries
#===============================================================================
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")


#Bring the headers, such as *.h into the project, should be renamed to include
include_directories(include)

# All the source code we care about
add_subdirectory(src)

file(COPY templates DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY third-party-libs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY include DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#===============================================================================
# Examples and tests
#===============================================================================
option(EEMA_EXAMPLES "Activate examples" OFF)
if(EEMA_EXAMPLES)
  add_subdirectory(examples)
endif()

#===============================================================================
# Documentation
#===============================================================================
option(EEMA_DOCUMENTATION "Build source documentation using Doxygen." OFF)
if(EEMA_DOCUMENTATION)
  add_subdirectory(doc)
endif()
