#include "functions.h"
#include "Constraint.h"

using namespace Eigen;

void Constraint::readConstraints(std::string name, int id, std::string master, std::string slave){
	constraint_name = name;
	constraint_id = id;
	master_name = master;
	slave_name = slave;
}

void Constraint::printInfo(){
	std::cout << "Constraint Type: " << constraint_name << "\n";
	std::cout << "Constraint ID: " << constraint_id << "\n";
	std::cout << "According to this constraint definition" << "\n" ;
	std::cout << "Master Mesh - " << master_name << "\n";
	std::cout << "Slave Mesh - " << slave_name << "\n";
}

std::string Constraint::getName(){
	return constraint_name;
}

std::string Constraint::getMaster(){
	return master_name;
}

std::string Constraint::getSlave(){
	return slave_name;
}