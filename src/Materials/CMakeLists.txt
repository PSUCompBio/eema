SET(Materials_src 
    fe_get_mats.cpp
    StressUpdate/fe_stressUpdate_1d.cpp
    StressUpdate/fe_stressUpdate.cpp
    Saint-Venant/fe_saintvenant_elastic.cpp
    Ogden-Hyperelastic/fe_ogden_hyperelastic.cpp
    Mooney-Rivlin-Hyperelastic/fe_mooneyrivlin_hyperelastic.cpp
    Elastic/fe_simple_elastic.cpp
     )
ADD_LIBRARY(Materials STATIC ${Materials_src})
target_link_libraries(Materials Math)
