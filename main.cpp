/*! \mainpage
 *
 * This is the property of The Penn State Computational Biomechanics Group. \n
 * This code is developed by Harsha Teja Garimella under the supervision of Dr. Reuben H Kraft.
 *
 * ## Motivation:
 * Computational Brain Biomechanics
 *
 * ## Acknowledgements:
 * Funding from CFDRC and Funding from ARL
 *
 * ## Contact:
 * Harsha T Garimella,
 * \n Ph.D. Candidate, Mechanical Engineering,
 * \n The Pennsylvania State University,
 * \n University Park, Pennsylvania, USA.
 * \n Email: harshatejagarimella@gmail.com
 *
 * Reuben H. Kraft, Ph.D.
 * \n Shuman Asst. Professor,
 * \n Department of Mechanical Engineering,
 * \n Department of Biomedical Engineering,
 * \n The Pennsylvania State University,
 * \n University Park, Pennsylvania, USA.
 * \n Email: reuben.kraft@psu.edu
 */

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <fstream>
#include "time.h"

#include <Eigen/Dense>

#include "functions.h"
#include "Mesh.h"
#include "Materials.h"

using namespace Eigen;

/*
 * Global Variables
 */

std::string home_path; /** Job folder name */
std::string job_file; /** Job input file name */

/** Dimension of the Simulation (3 for 3D, 2 for 2D and 1 for 1D) */
int ndof;

/** Input meshes - number of meshes and type of meshes */
int num_meshes;
Mesh *mesh;

/* Time Related Initialized Variables Stored Here */
double t_start; /** Start Time */
double t_end; /** End Time */
int output_frequency;
double dt_initial; /** Initial time step */
double dt_min; /** Minimum time step */
double reduction; /** Reduces the timestep value by this amount */

/** Material Properties are stored in here */
int material_types;
Materials *mat;

int num_load_constraints;
double input_load_amp; /** Loading Amplitude */
VectorXi fcdof; /** Degrees of freedom -- corresponding to the applied load */
std::string fc_curve; /** Type of loading curve */

int num_disp_constraints;
double input_disp_amp;
VectorXi dcdof;
std::string dc_curve; /** Type of displacement curve */

/* Boundary Conditions */
/** Number of boundary conditions = nodes x dof */
int num_constraints;
/** DOF numbers for boundary conditions */
VectorXi bcdof;
/** Displacement of the constrained DOF's */
VectorXd bcval;

/* Results Stored in these variables */
MatrixXd element_stress_host; /** Stores the stress tensor of an element */
MatrixXd element_strain_host; /** Stores the strain tensor of an element */
MatrixXd element_stress_truss; /** Stores the stress of a truss element */
MatrixXd element_strain_truss; /** Stores the strain of a truss element */

/* Internal nodal force at previous and current timesteps.
   Used in the calculation of internal energy of the system. */
VectorXd fi_prev; /** Internal nodal force vector at previous time step */
VectorXd fi_curr; /** Internal nodal force vector at current time step */
VectorXd W_int; /** Internal energy at different time steps */

VectorXd U_host; /** Variable storing the nodal displacements of the host mesh */
VectorXd V_host; /** Variable storing the nodal velocities of the host mesh */
VectorXd A_host; /** Variable storing the nodal accelerations of the host mesh */

VectorXd U_truss; /** Variable storing the nodal displacements of the embedded mesh */
VectorXd V_truss; /** Variable storing the nodal velocities of the embedded mesh */
VectorXd A_truss; /** Variable storing the nodal accelerations of the embedded mesh */

/**
 * This is the main file. If you want to submit a new job -- this is where you do it
 */

int main(){

	clock_t t;
	t = clock();

	/** Enter the path address for your job folder */
	home_path = "./../Jobs/job-3/";
	job_file = "input.inp";

	/** What is the Problem Dimension - Degrees of Freedom per node? */
	ndof = 3; // DOFs per node

	fe_mainRead(home_path+job_file);

	// Initialize Time Variables
	t_start = 0; /** Simulation start time */
	t_end = 0.001; /** Simulation end time */
	output_frequency = 1; /** Output frequency -- result ouputs for every 100 timesteps */
	dt_initial = 1*pow(10,-6); /** Maximum time step */
	dt_min = dt_initial;
	reduction = 0.5; /** Reduction factor */

	/** Type of Loading
	   1 - Force based
	   2 - Displacement based
	 */
	int type_of_loading = 1;

	if(type_of_loading == 1) {
		// Input the loading conditions
		input_load_amp  = 0.08;/** Loading Amplitude */
		num_load_constraints = 4; /** number of nodes * number of directions */
		fcdof = VectorXi::Zero(num_load_constraints);
		fcdof << 14, 17, 20, 23; /** Degrees of freedom index data */
		fc_curve = "RAMP";
	}

	if(type_of_loading == 2) {
		// Nodes upon which displacement is forced
		input_disp_amp = 0.001;
		num_disp_constraints = 4;
		dcdof = VectorXi::Zero(num_disp_constraints);
		dcdof << 14, 17, 20, 23;
		dc_curve = "RAMP";
	}

	// Boundary Conditions
	num_constraints = 12; /** Number of boundary condition constraints */
	bcdof = VectorXi::Zero(num_constraints); /** Which nodes are constrained ? */
	bcdof << 0, 1, 2, 4, 5, 8, 9, 11, 12, 13, 16, 21; /** index corresponding to the constrained degrees of freedom */
	bcval = VectorXd::Zero(num_constraints);
	bcval << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; /** Displacements of the constrained degrees of freedom */

	/** Enter your choice for type of simulations below, based on the following options:
	 * 1 - Explicit Dynamic
	 * 2 - Implicit Dynamic
	 * 3 - Static Analysis
	 */

	int choice_simulation = 1;
	switch(choice_simulation) {
	case 1:
		fe_mainEXPLICIT();
		break;

	default: // Currently default is set to explicit dynamic simulations
		fe_mainEXPLICIT();
		break;
	}

	t = clock()-t;

	std::cout<<"--------------------------------------"<<"\n";
	std::cout<<"Total Simulation CPU Time: "<<(((float)t)/CLOCKS_PER_SEC)<<"s \n";
	std::cout<<"Simulation Completed."<<"\n";
	std::cout<<"--------------------------------------"<<"\n";
	return 0;
}
