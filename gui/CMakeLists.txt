# configure Qt

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Core Widgets REQUIRED)
#find_package(Qt5Widgets REQUIRED)
#find_package(Qt5Test REQUIRED)

# Testing configuration
enable_testing()
set(TEST_LINK_LIBRARIES Qt5::Test)



configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defines.h.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/defines.h)

file(GLOB_RECURSE UI_FILES *.ui)
file(GLOB_RECURSE CODE_FILES *.cpp)


qt5_wrap_ui(UI_HEADERS ${UI_FILES})
qt5_add_resources(RESOURCE_FILES ${CMAKE_SOURCE_DIR}/gui/resources/resources.qrc)

# Windows application icon
if (WIN32)
  set(WINDOWS_RES_FILE ${CMAKE_CURRENT_BINARY_DIR}/resources.obj)
  if (MSVC)
    add_custom_command(OUTPUT ${WINDOWS_RES_FILE}
      COMMAND rc.exe /fo ${WINDOWS_RES_FILE} resources.rc
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/win
    )
  else()
    add_custom_command(OUTPUT ${WINDOWS_RES_FILE}
      COMMAND windres.exe resources.rc ${WINDOWS_RES_FILE}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/win
    )
  endif()
endif()

add_executable(eema_gui WIN32
  ${UI_HEADERS}
  ${CODE_FILES}
  ${RESOURCE_FILES}
  ${WINDOWS_RES_FILE}
)

target_link_libraries(eema_gui
  Qt5::Widgets
)
