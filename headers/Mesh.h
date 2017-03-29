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
	MatrixXi elements;
	MatrixXd nodes_new;
	MatrixXi elements_new;
public:
	void readMesh(MatrixXd n, MatrixXi e);
	MatrixXd getNodes();
	MatrixXi getElements();
	MatrixXd getNewNodes();
	MatrixXi getNewElements();
	void preprocessMesh();
};

#endif /* HEADERS_MESH_H_ */
