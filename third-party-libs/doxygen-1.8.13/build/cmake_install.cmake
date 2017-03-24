# Install script for directory: /Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/libmd5/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/qtools/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/vhdlparser/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/src/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/examples/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/doc/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/addon/doxmlparser/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/addon/doxyapp/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/addon/doxysearch/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/addon/doxywizard/cmake_install.cmake")
  include("/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/testing/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/vsg111/Dropbox/Work/Papers/Paper_EEM_Computational/EEM_Dynamic/third-party-libs/doxygen-1.8.13/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
