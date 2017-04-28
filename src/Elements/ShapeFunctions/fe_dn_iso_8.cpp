
#include"functions.h"

using namespace Eigen;

VectorXd fe_dndr_8(double rvalue, double svalue, double tvalue) {

	VectorXd dndr(8);
	int i = 0;
	double tmp = 0.125 * (-1) * (1 - svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1) * (1 - svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1) * (1 + svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (-1) * (1 + svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (-1) * (1 - svalue) * (1 + tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (+1) * (1 - svalue) * (1 + tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (+1) * (1 + svalue) * (1 + tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (-1) * (1 + svalue) * (1 + tvalue);
	dndr(i) = tmp;

	return dndr;
}

VectorXd fe_dnds_8(double rvalue, double svalue, double tvalue) {

	VectorXd dnds(8);
	int i = 0;
	double tmp = 0.125 * (1 - rvalue) * (-1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (-1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (+1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (+1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (-1) * (1 + tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (-1) * (1 + tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (+1) * (1 + tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (+1) * (1 + tvalue);
	dnds(i) = tmp;

	return dnds;
}

VectorXd fe_dndt_8(double rvalue, double svalue, double tvalue) {

	VectorXd dndt(8);
	int i = 0;
	double tmp = 0.125 * (1 - rvalue) * (1 - svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 - svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 + svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (1 + svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (1 - svalue) * (+1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 - svalue) * (+1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 + svalue) * (+1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (1 + svalue) * (+1);
	dndt(i) = tmp;

	return dndt;
}

/****************************************/

void fe_dndr_8_pbr(VectorXd& dndr, double rvalue, double svalue, double tvalue) {


	int i = 0;
	double tmp = 0.125 * (-1) * (1 - svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1) * (1 - svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1) * (1 + svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (-1) * (1 + svalue) * (1 - tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (-1) * (1 - svalue) * (1 + tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (+1) * (1 - svalue) * (1 + tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (+1) * (1 + svalue) * (1 + tvalue);
	dndr(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (-1) * (1 + svalue) * (1 + tvalue);
	dndr(i) = tmp;


}

void fe_dnds_8_pbr(VectorXd& dnds, double rvalue, double svalue, double tvalue) {


	int i = 0;
	double tmp = 0.125 * (1 - rvalue) * (-1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (-1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (+1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (+1) * (1 - tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (-1) * (1 + tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (-1) * (1 + tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (+1) * (1 + tvalue);
	dnds(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (+1) * (1 + tvalue);
	dnds(i) = tmp;

}

void fe_dndt_8_pbr(VectorXd& dndt, double rvalue, double svalue, double tvalue) {


	int i = 0;
	double tmp = 0.125 * (1 - rvalue) * (1 - svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 - svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 + svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (1 + svalue) * (-1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (1 - svalue) * (+1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 - svalue) * (+1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 + rvalue) * (1 + svalue) * (+1);
	dndt(i) = tmp;

	i = i + 1;
	tmp = 0.125 * (1 - rvalue) * (1 + svalue) * (+1);
	dndt(i) = tmp;


}


