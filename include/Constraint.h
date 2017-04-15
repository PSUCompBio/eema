#ifndef HEADERS_CONSTRAINT_H_
#define HEADERS_CONSTRAINT_H_

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace Eigen;

class Constraint{

private:
	std::string constraint_name;
	int constraint_id;
	int master_name;
	int slave_name;
public:
	void readConstraints();
	void printInfo();
	std::string getName();
	std::string getMaster();
	std::string getSlave();
}

#endif /* HEADERS_CONSTRAINT_H_ */