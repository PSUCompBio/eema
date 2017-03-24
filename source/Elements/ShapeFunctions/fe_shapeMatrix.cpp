

#include"functions.h"

using namespace Eigen;

MatrixXd fe_shapeMatrix(int edof, int nnel, VectorXd shapes){
	
	MatrixXd shape_mat(6,edof);
	shape_mat = MatrixXd::Zero(6,edof);
	int i1,i2,i3;

	for(int i=0;i<nnel;i++){
		
		i1 = 0 + i*3;
		i2 = i1 + 1;
		i3 = i2 + 1;
		
		// Shape Matrix for Constructing Mass Matrix
		shape_mat(0,i1) = shapes(i);
		shape_mat(1,i2) = shapes(i);
		shape_mat(2,i3) = shapes(i);

	}
	//std::cout<<"Exiting Shape Function Matrix Generator\n";
	return shape_mat;
}
