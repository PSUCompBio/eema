cmake_minimum_required(VERSION 2.6)
project(eema)
enable_language(CXX)

# Eigen support
option (USE_EIGEN  "Build library with Eigen support" ON)
set (EIGEN_INCLUDE_DIR "INCLUDE directory for Eigen library" third-party-libs/eigen3)
include_directories(${EIGEN_INCLUDE_DIR})

#Bring the headers, such as *.h into the project, should be renamed to include
include_directories(include)

# All the source code we care about
#file(GLOB SRC_CORE "src/*/*.cpp")
add_subdirectory(src) 

# this creates a static library (i.e. libeem.a)
#add_library(eema STATIC ${SRC_CORE})
#install(TARGETS eema DESTINATION ./lib)

#add_executable(eem  ${SOURCES})
#===============================================================================
# Examples and tests
#===============================================================================
option(EEMA_EXAMPLES "Activate examples" OFF)

if(EEMA_EXAMPLES)
  add_subdirectory(examples)
endif()


