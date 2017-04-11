SET(FEM_src 
    RHSNodalForce/fe_update.cpp
    RHSNodalForce/fe_getforce.cpp
    Explicit/fe_mainEXPLICIT.cpp
    Electromechanics/fe_electromechanics.cpp
    BoundaryConditions/fe_apply_bc.cpp
    Assembly-Global-Local/fe_scatter.cpp
    Assembly-Global-Local/fe_gather.cpp
    Assembly-Global-Local/fe_find_index.cpp
    Assembly-Global-Local/fe_assemble_mass.cpp
     )
ADD_LIBRARY(FEM OBJECT ${FEM_src})
