SET(example1_src main.cpp)

ADD_EXECUTABLE(eema_example1 ${example1_src})

TARGET_LINK_LIBRARIES(eema_example1 eema)

file(COPY input_load.inp DESTINATION .)
file(COPY input_displacement.inp DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)
