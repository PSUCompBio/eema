SET(Input_Output_src 
    Input/fe_text2matrix.cpp
    Input/fe_mainRead.cpp
    Input/Mesh.cpp
    Input/Materials.cpp
    Output/fe_write.cpp
    Output/fe_vtu.cpp
    Output/fe_vtk.cpp
    Output/fe_vector2text.cpp
    Output/fe_simulationTime.cpp
    Output/fe_matrix2text.cpp
    Output/fe_display.cpp
     )
ADD_LIBRARY(Input-Output STATIC ${Input_Output_src})
target_link_libraries(Input-Output Materials)
