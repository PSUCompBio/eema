#include "functions.h"

using namespace Eigen;

/** This function calculates the wave speed for an element based on its material properties */
double fe_calWaveSpeed(double E, double nu, double rho){

	double c_wave;

	if(nu==0.5){
		c_wave = sqrt(E/rho);
	} else{
		c_wave = sqrt((E*(1-nu))/(rho*(1+nu)*(1-(2*nu))));
	}

	return c_wave;

}
