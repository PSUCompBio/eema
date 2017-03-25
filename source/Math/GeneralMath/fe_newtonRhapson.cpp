#include "functions.h"

using namespace Eigen;

VectorXd fe_newtonRhapson(VectorXd nat_coord, double xcoord[], double ycoord[], double zcoord[], int nnel){

    VectorXd coord = VectorXd::Zero(nat_coord.size()*nnel);
    double eps = 1e6;

    for (int i=0;i<nnel;i++){
        coord((i*3)) = xcoord[i];
        coord((i*3)+1) = ycoord[i];
        coord((i*3)+2) = zcoord[i];
    }

    double edof = coord.size();
    VectorXd iso_coord = VectorXd::Zero(3); /** Vector showing the coordinates in iso-parametric system */

    if(nnel==8){
        VectorXd shapes = fe_shapes_8(iso_coord(0),iso_coord(1),iso_coord(2));
    }

    MatrixXd shape_mat = fe_shapeMatix(edof,nnel,shapes);
    VectorXd fx = VectorXd::Zero(3);


	dndr = fe_dndr_8(iso_coord(0),iso_coord(1),iso_coord(2));
	dnds = fe_dnds_8(iso_coord(0),iso_coord(1),iso_coord(2));
	dndt = fe_dndt_8(iso_coord(0),iso_coord(1),iso_coord(2));
	jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
    MatrixXd transJacobian = jacobian.transpose();
    MatrixXd it_Jacobian = transJacobian.inverse();
    
    
    while(eps >= eps_nr){

        iso_coord_new = iso_coord - (it_jacobian*(fx))
        iso_coord = iso_coord_new;
        
        shapes = fe_shapes_8(iso_coord(0),iso_coord(1),iso_coord(2));
        shape_mat = fe_shapeMatix(edof,nnel,shapes);
        fx = (nat_coord - (shape_mat*coord)); 

        eps = fx.norm();

        dndr = fe_dndr_8(iso_coord(0),iso_coord(1),iso_coord(2));
	    dnds = fe_dnds_8(iso_coord(0),iso_coord(1),iso_coord(2));
	    dndt = fe_dndt_8(iso_coord(0),iso_coord(1),iso_coord(2));
        jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
        transJacobian = jacobian.transpose();
        it_Jacobian = transJacobian.inverse(); 

        

    }



    return iso_coord;
}