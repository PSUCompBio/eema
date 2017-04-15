#include "functions.h"

using namespace Eigen;

VectorXd fe_getForce_3d(VectorXd U, VectorXd V, VectorXd fe, int time_step_counter){

	VectorXd f_int_e;


	MatrixXd nodes_host = mesh[0].getNewNodes();
	MatrixXi elements_host = mesh[0].getNewElements();
	MatrixXd nodes_embed = mesh[1].getNewNodes();
	MatrixXi elements_embed = mesh[1].getNewElements();

	return f_int_e;
}