#include"functions.h"

using namespace Eigen;

VectorXd fe_apply_bc_displacement(VectorXd U,double time){

	if(bcdof.size()!=0){
		for(int i=0;i<bcdof.size();i++){
			int c = bcdof(i);
			U(c) = bcval(i);
		}
	}

	if(dcdof.size()!=0){
		for(int i=0;i<dcdof.size();i++){
			int c = dcdof(i);
			U(c) = fe_function(input_disp_amp,dc_curve,time);
		}
	}

	return U;
}

VectorXd fe_apply_bc_load(VectorXd fe,double time){

	if(fcdof.size()!=0){
		for(int i=0;i<fcdof.size();i++){
			int c = fcdof(i);
			fe(c) = fe_function(input_load_amp,fc_curve,time);
		}
	}

	return fe;
}

MatrixXd fe_apply_bc_stiffness(MatrixXd kk,VectorXi bcdof,VectorXd bcval){

	int n = bcdof.size();
	int sdof = kk.rows();

	for(int i=0;i<n;i++){
		int c = bcdof(i);
		for(int j=0;j<sdof;j++){
			kk(c,j) = 0;
		}
		kk(c,c) = 1;
	}
	return kk;
}
