/*
 * Mesh.h
 *
 *  Created on: Dec 15, 2016
 *      Author: vsg111
 */

#ifndef HEADERS_MESH_H_
#define HEADERS_MESH_H_

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace Eigen;

class Mesh{
	MatrixXd nodes;
	MatrixXd elements;
	MatrixXd nodes_new;
	MatrixXd elements_new;
public:
	void readMesh(MatrixXd n, MatrixXd e);
	MatrixXd getNodes();
	MatrixXd getElements();
	MatrixXd getNewNodes();
	MatrixXd getNewElements();
	void preprocessMesh();
};

#endif /* HEADERS_MESH_H_ */
