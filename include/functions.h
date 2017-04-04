#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>
#include<iomanip>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "unsupported/Eigen/MatrixFunctions"

#include "Mesh.h"
#include "Materials.h"

using namespace Eigen;

/*
 * Global Variables
 */

/* Job Folder */
extern std::string home_path;
extern std::string job_file;

/* Dimension of the Simulation (3 for 3D, 2 for 2D and 1 for 1D) */
extern int ndof;

/* Input Host and Truss Meshes */
// extern Mesh host, truss;
extern int num_meshes;
extern Mesh *mesh;

/* Time Related Initialized Variables Stored Here */
extern double t_start;
extern double t_end;
extern int output_frequency;
extern double dt_initial;
extern double dt_min;
extern double reduction;

/* Material Properties are stored in here */
extern int material_types;
extern Materials *mat;

/* Loading Conditions */
extern int num_load_constraints;
extern double input_load_amp;
extern VectorXi fcdof;
extern std::string fc_curve;

/** Displacement Loading Conditions */
extern int num_disp_constraints;
extern double input_disp_amp;
extern VectorXi dcdof;
extern std::string dc_curve;

/* Boundary Conditions */
extern int num_constraints;
extern VectorXi bcdof;
extern VectorXd bcval;

/* Results Stored in these variables */
extern MatrixXd element_stress_host;
extern MatrixXd element_strain_host;
extern MatrixXd element_stress_truss;
extern MatrixXd element_strain_truss;

/* Internal nodal force at previous and current timesteps.
Used in the calculation of internal energy of the system. */
extern VectorXd fi_prev;
extern VectorXd fi_curr;
extern VectorXd W_int;

/* System External Energy */
extern VectorXd fe_prev; 
extern VectorXd fe_curr;
extern VectorXd W_ext;

/* System Kinetic Energy */
extern VectorXd W_kin;

/* System Total Energy */
extern VectorXd W_tot;

/* System Max Energy */
extern VectorXd W_max;

extern VectorXd U_host;
extern VectorXd V_host;
extern VectorXd A_host;

extern VectorXd U_truss;
extern VectorXd V_truss;
extern VectorXd A_truss;

extern double eps_nr; /** Convergence Criteria for Newton Rhapson Method */
extern double iterations_nr;
extern double eps_energy;

/*
 * Functions
 */

/** Reads a text file into matrix */
MatrixXd text2matrix(std::string name, int cols);

/** Creates the shape functions for an 8 noded element */
VectorXd fe_shapes_8(double rvalue, double svalue, double tvalue);

/** dndr of isoparametric element calculated for particular r, s, and t */
VectorXd fe_dndr_8(double rvalue, double svalue, double tvalue);

/** dnds of isoparametric element calculated for particular r, s, and t */
VectorXd fe_dnds_8(double rvalue, double svalue, double tvalue);

/** dndt of isoparametric element calculated for particular r, s, and t */
VectorXd fe_dndt_8(double rvalue, double svalue, double tvalue);

/** Create a guass_point vector of n values */
VectorXd guass_points(int n);

/** Creates a guass_weight vector of n values */
VectorXd guass_weights(int n);

/** Creates a guass point matrix in 3D */
MatrixXd guass_points_3d(int nx, int ny, int nz);

/** Creates a guass weight matrix in 3D */
MatrixXd guass_weights_3d(int ndof, int nx, int ny, int nz);

/** Create material matrix for isotropic elastic case */
MatrixXd fe_calculate_matlmat(int n, double E, double nu);

/** Calculates the jacobian -- using the derivates of shape functions */
MatrixXd fe_calJacobian(int dim, int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);

/** dndx of actual element calculates using jacobian and shape function derivates calculated in the isoparametric element */
VectorXd fe_dndx_8(int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, MatrixXd invJacobian);

/** dndy of actual element calculates using jacobian and shape function derivates calculated in the isoparametric element */
VectorXd fe_dndy_8(int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, MatrixXd invJacobian);

/** dndz of actual element calculates using jacobian and shape function derivates calculated in the isoparametric element */
VectorXd fe_dndz_8(int nnel, VectorXd dndr, VectorXd dnds, VectorXd dndt, MatrixXd invJacobian);

/** Strain displacement matrix B */
MatrixXd fe_strDispMatrix(int edof, int nnel, VectorXd dndx, VectorXd dndy, VectorXd dndz);

MatrixXd fe_apply_bc_stiffness(MatrixXd kk, VectorXi bcdof, VectorXd bcval);

/** Writes a matrix into a text file */
void matrix2text(std::string name,MatrixXd new_slave_master,int width);

/** Writes a vector into a text file */
void vector2text(std::string name, VectorXd vector, int width);

/** Calculates the transformation matrix - transformation from local (truss) coordinate system to global (3d hex) coordinate system */
MatrixXd fe_calTransformation(MatrixXd truss_nodes, int choice);

MatrixXd fe_calSimpTransformation(MatrixXd truss_nodes);

/** Internal nodal force vector for a hexahedral element */
MatrixXd fe_stiffness_hex(double E, double nu, int ndof, int nnel, int edof, double xcoord[], double ycoord[], double zcoord[]);

/** Internal nodal force vector for a truss (1D) element */
MatrixXd fe_stiffness_embed_truss(MatrixXd nodes_truss, MatrixXd elements_truss, double E_truss, double A_truss, int ndof, int nnel, int edof, double xcoord[], double ycoord[], double zcoord[]);

/** Outputs the shape function matrix for an element */
MatrixXd fe_shapeMatrix(int edof, int nnel, VectorXd shapes);

/** Calculates the mass matrix for a hex element */
MatrixXd fe_mass_hex(double rho, int ndof, int nnel, int edof, VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);

/** Updates the stress at each time step based on the material model */
VectorXd fe_stressUpdate(VectorXd dndx, VectorXd dndy, VectorXd dndz, MatrixXd disp_mat, VectorXd u, int opt, int return_opt);

/** Calculates the resultant force vector - Box 6.1 of Belytschko */
VectorXd fe_getforce(MatrixXd nodes, MatrixXi elements, int ndof, VectorXd u, VectorXd v, VectorXd fext, int size_counter, MatrixXd nodes_truss, MatrixXi elements_truss);

/** Outputs the critical time step based on all the elements in a FE analysis */
double fe_getTimeStep(MatrixXd nodes, MatrixXi elements, int ndof, VectorXd u, VectorXd v, VectorXd fext);

/** Calculates the time step for a single element based on its dimensions and material model */
double fe_calTimeStep(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord, double E, double nu, double rho);

/** Calculates the area of a face with 4 vertices */
double fe_calArea_4(double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double c1, double c2, double c3, double c4);

/** Calculates the volume of a 3d element */
double fe_calVolume(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);

/** Prints out a vector */
void fe_display_vector(VectorXd A);

/** Prints out a matrix */
void fe_display_matrix(MatrixXd A);

/** Converts a normal mass matrix into a lumped mass matrix */
MatrixXd fe_transformMass(MatrixXd m, int opt);

/** Extracts the material parameter values based on the material id */
double fe_get_mats(int matl_code, int obj_interest);

std::string fe_get_model(int matl_code);

/** Calculates the strain displacement matrix in total lagrangian system */
MatrixXd fe_strDispMatrix_totalLagrangian(int edof, int nnel, VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u);

/** Calculates the deformation gradient */
MatrixXd fe_calDefGrad(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u);

/** Calculates the wavespeed for a particular material model */
double fe_calWaveSpeed(double E, double nu, double rho);

void fe_writeElementStress(MatrixXd sigma_all, double time);

/** Calculates the mass of a truss element */
MatrixXd fe_mass_truss(double rho,double A_truss,int edof,MatrixXd nodes,MatrixXd elements);

/** Updates the stress at each time step based on the material model for a 1d element */
VectorXd fe_stressUpdate_1d(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u_e, int opt, MatrixXd nodes);

/** Writes VTK files for host mesh */
void fe_vtkWrite_host(std::string output,int format_choice,int mesh_choice,int time_step,MatrixXd nodes, MatrixXi elements);

/** Writes VTK files for truss mesh */
void fe_vtkWrite_truss(std::string output,int format_choice,int mesh_choice,int time_step,MatrixXd nodes, MatrixXi elements);

/** Updates the nodal coordinates based on the displacements */
MatrixXd fe_updateNodes(MatrixXd nodes,VectorXd displacements);

/** Run the finite element analysis using an explicit dynamic method */
void fe_mainEXPLICIT();

/* Date: March 15, 2017 */
/** Find the index based on the DOF of a particular node */
VectorXi fe_find_index(VectorXi node_list);

/** Assembles the global mass matrix */
MatrixXd fe_assemble_mass(MatrixXd mm, MatrixXd m, VectorXi node_list, int sdof);

/** Creates element level vector (displacement, velocity, acceleration etc.) from a system level vector */
VectorXd fe_gather(VectorXd global_vec, VectorXd local_vec, VectorXi node_list, int sdof);

/** Updates a system level vector based on the element level vector */
VectorXd fe_scatter(VectorXd global_vec, VectorXd local_vec, VectorXi node_list, int sdof);

/* Date: March 20, 2017 */

/** Writes a symmetric tensor into a Voigt vector form */
VectorXd fe_tensor2voigt(MatrixXd A);

/** Writes a Voigt vector into a symmetric tensor form */
MatrixXd fe_voigt2tensor(VectorXd B);

/** Concatenate a vector to a matrix -- rowwise or coloumn wise */
MatrixXd fe_concatenate_vector2matrix(MatrixXd A, VectorXd B,int opt);

/** find the poistion index of a value in a vector -- analogous to 'find' function in MATLAB */
int fe_find(VectorXd A, double a);

int fe_find(VectorXd A, int a);

/** Read the input text file -- for a particular job */
void fe_mainRead(std::string file);

/* Date: March 21, 2017 */
/** Function outputs the standard curve values */
double fe_function(double a, std::string b, double time);

VectorXd fe_apply_bc_displacement(VectorXd U,double time);

VectorXd fe_apply_bc_load(VectorXd fe,double time);

/* Date: March 22, 2017 */
VectorXd fe_mooneyrivlin_hyperelastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u, int opt, int return_opt);

VectorXd fe_ogden_hyperelastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u, int opt, int return_opt);

VectorXd fe_simple_elastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, MatrixXd disp_mat, VectorXd u, int opt, int return_opt);

VectorXd fe_saintvenant_elastic(VectorXd dndx, VectorXd dndy, VectorXd dndz, VectorXd u, int opt, int return_opt);

void fe_vtuWrite(std::string output, int time_step, Mesh mesh1);

/* Date: March 26, 2017 */

VectorXd fe_newtonRhapson(VectorXd nat_coord, VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);

double fe_minElementLength(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);

double fe_maxElementLength(VectorXd xcoord, VectorXd ycoord, VectorXd zcoord);

/* Embedded Element Constraint */

MatrixXd fe_insert_vector2matrix(MatrixXd A, VectorXd B, int num, int opt);

VectorXd fe_embed_preprocessing_mapping(Mesh host, Mesh embed);

VectorXd fe_embed_preprocessing(Mesh host,Mesh embed);

void fe_embed_preprocessing_length(Mesh host, Mesh embed);

int fe_compute_host(VectorXd A, MatrixXd nodes_host, MatrixXd elements_host_tmp);

MatrixXd fe_create_bbox(VectorXd A, MatrixXd nodes_host, MatrixXd elements_host, double length);

#endif
