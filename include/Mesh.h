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
	/** This contains the nodal data read from the mesh */
	MatrixXd nodes;
	/** This contains the element data read from the mesh */
	MatrixXi elements;
	/** This contains the nodal data after preprocessing or any other operation */
	MatrixXd nodes_new;
	/** This contains the element data after preprocessing or any other operation */
	MatrixXi elements_new;
public:
	void readMesh(MatrixXd n, MatrixXi e);
	MatrixXd getNodes();
	MatrixXi getElements();
	MatrixXd getNewNodes();
	MatrixXi getNewElements();
	VectorXd getMaxCharLength(std::string choice);
	VectorXd getMinCharLength(std::string choice);
	void preprocessMesh();
	void replaceNodes(MatrixXd A, std::string B);
	void replaceElements(MatrixXi A, std::string B);
	void checkMesh();
};

#endif /* HEADERS_MESH_H_ */
