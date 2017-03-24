#include "functions.h"
#include "Mesh.h"

using namespace Eigen;

void Mesh::readMesh(MatrixXd n, MatrixXd e){
	nodes = n;
	elements = e;
}

MatrixXd Mesh::getNodes(void){
	return nodes;
}

MatrixXd Mesh::getElements(void){
	return elements;
}

MatrixXd Mesh::getNewNodes(void){
	return nodes_new;
}

MatrixXd Mesh::getNewElements(void){
	return elements_new;
}

void Mesh::preprocessMesh(void){

	nodes_new = nodes;
	elements_new = elements;
	/** Nodes Preprocessing - Putting the numbering in order */
	for(int i=0;i<nodes.rows();i++){
		nodes_new(i,0) = i;
	}
	/** Elements Preprocessing - Correcting the element definitions */
	for(int i=0;i<elements.rows();i++){
		elements_new(i,0) = i;
		elements_new(i,1) = elements(i,1);
		for(int j=2;j<elements.cols();j++){
			elements_new(i,j) = fe_find(nodes.col(0),elements(i,j)); /** find_in_matrix gives the row number of node in nodes matrix */
		}
	}
}
