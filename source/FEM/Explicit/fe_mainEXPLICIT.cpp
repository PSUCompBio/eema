#include "functions.h"
using namespace Eigen;

/*! \brief
 * This function carries out the explicit dynamic analysis of the FEM problem.
 */

void fe_mainEXPLICIT(){

	for(int i=0;i<num_meshes;i++){
		mesh[i].preprocessMesh();
	}

	MatrixXd nodes = mesh[0].getNewNodes();
	MatrixXi elements = mesh[0].getNewElements();
	MatrixXd nodes_truss = mesh[1].getNewNodes();
	MatrixXi elements_truss = mesh[1].getNewElements();

	// Following variables - Only for Hex Element
	int nel = elements.rows(); /*! number of elements */
	int nnel = (elements.cols()-2); // number of nodes per element
	int nnode = nodes.rows(); // number of nodes
	int sdof = nnode * ndof; // system degrees of freedom
	int edof = nnel * ndof; // element degrees of freedom

	// Element Stress and Strain
	element_stress_host = MatrixXd::Zero(nel,(3*ndof));
	element_strain_host = MatrixXd::Zero(nel,(3*ndof));
	element_stress_truss = MatrixXd::Zero(elements_truss.rows(),(3*ndof));
	element_strain_truss = MatrixXd::Zero(elements_truss.rows(),(3*ndof));
	U_host = VectorXd::Zero(sdof);

	// Updated Nodes and Elements
	MatrixXd updated_nodes = nodes;

	// Writing the Read Mesh
	//fe_vtkWrite_host("eem_matrix",1,5,0,updated_nodes,elements);
	//return;

	// Initialization
	double dT = dt_initial;
	double t_half = 0; // t(n + 1/2)
	double t=t_start;
	int size_counter = 0; // time step count
	int time_temp_1 = 1;
	double output_temp_1 = ((double)(t_end/output_frequency)); 

	VectorXd A = VectorXd::Zero(sdof); // Acceleration Vector
	MatrixXd AA = MatrixXd::Zero(A.size(),1); // Matrix storing nodal accln at different time steps
	VectorXd V = VectorXd::Zero(sdof); // Velocity Vector
	MatrixXd VV = MatrixXd::Zero(V.size(),1); // Matrix storing nodal velocities at different time steps
	VectorXd V_half = VectorXd::Zero(sdof); // Velocity Vector at n+1/2
	VectorXd U = VectorXd::Zero(sdof); // Displacement Vector
	MatrixXd UU = MatrixXd::Zero(U.size(),1); // Matrix storing nodal displacements at different time steps

	VectorXd fe = VectorXd::Zero(sdof); // External Nodal force vector

	fi_prev = VectorXd::Zero(sdof); // Internal nodal force vector at previous timestep
	fi_curr = VectorXd::Zero(sdof); // Internal Nodal force vector at current timestep
	W_int = VectorXd::Zero(1);
	VectorXd U_prev = VectorXd::Zero(sdof); // Nodal displacements at previous time
	VectorXd U_curr = VectorXd::Zero(sdof); // Nodal displacements at current time

	fe_prev = VectorXd::Zero(sdof);
	fe_curr = VectorXd::Zero(sdof);
	W_ext = VectorXd::Zero(1);

	W_kin = VectorXd::Zero(1);
	W_tot = VectorXd::Zero(1);


	VectorXd F_net = VectorXd::Zero(sdof); // Total Nodal force vector


	MatrixXd mm = MatrixXd::Zero(sdof,sdof); // System Mass Matrix

	// Loading Conditions
	fe = fe_apply_bc_load(fe,0);

/* The algorithm was developed based on the Belytschko Textbook */
// ----------------------------------------------------------------------------
	//Step-1: Initialize and Calculate the Mass Matrix of the system
	VectorXi nodes_local = VectorXi::Zero(nnel);
	VectorXd xcoord = VectorXd::Zero(nnel);
	VectorXd ycoord = VectorXd::Zero(nnel);
	VectorXd zcoord = VectorXd::Zero(nnel);

	for(int i=0;i<nel;i++){

		double rho = fe_get_mats(elements(i,1),0); // material density

		for(int j=0;j<nnel;j++){
			int g = -1;
			for(int f=0;f<nnode;f++){
				if(elements(i,j+2)==nodes(f,0)){
					g = f;
					break;
				}
			}
			nodes_local(j) = elements(i,j+2);
			xcoord(j) = nodes(g,1);
			ycoord(j) = nodes(g,2);
			zcoord(j) = nodes(g,3);
		}

		MatrixXd m_hex = MatrixXd::Zero(edof,edof); // mass of hex elements
		m_hex = fe_mass_hex(rho, ndof, nnel, edof, xcoord, ycoord, zcoord);

		/* TRUSS ANALYSIS - WORKS FOR ONLY SINGLE ELEMENT
		// Truss elements
		// by default, mass of the truss is added to the hex element in a lumped form
		double rho_truss = fe_get_mats(elements_truss(i,1),0);
		double A_truss = 0.5;
		MatrixXd m_truss = MatrixXd::Zero(edof,edof);
		m_truss = fe_mass_truss(rho_truss,A_truss,edof,nodes_truss,elements_truss);

		MatrixXd m_correction = MatrixXd::Zero(edof,edof);
		m_correction = fe_mass_truss(rho,A_truss,edof,nodes_truss,elements_truss); // This is the redundant volume

		// mm = fe_assemble_mass(m_hex+m_truss); // This is where the mass magic happens (m_hex + m_truss - m_correction)*/

		mm = fe_assemble_mass(mm,m_hex,nodes_local,sdof); // With Volume Redundancy
	}

	mm = fe_transformMass(mm,2); // Transforming consistent mass to lumped mass

	std::string mass = home_path+"results/system_mass.txt";
	matrix2text(mass.c_str(),mm,mm.cols());

    //std::string mass_inverse = home_path+"results/system_mass_inverse.txt";
    //matrix2text(mass_inverse.c_str(),mm.inverse(),24);

// ----------------------------------------------------------------------------
	//Step-2: getforce step from Belytschko
	F_net = fe_getforce(nodes,elements,ndof,U,V,fe,size_counter,nodes_truss,elements_truss);
	fe_vtuWrite("eem_matrix",size_counter,mesh[0]);
	fe_vtuWrite("eem_truss",size_counter,mesh[1]);
	

	dT = reduction * fe_getTimeStep(nodes,elements,ndof,U,V,fe);
	if(dT>dt_min){
		dT = dt_min;
	}

// ----------------------------------------------------------------------------

	//Step-3: Calculate accelerations
	A  = mm.inverse()*(F_net);

	UU.col(0) = U;
	VV.col(0) = V;
	AA.col(0) = A;
	U_host = U;
	V_host = V;
	A_host = A;
	U_prev = U;

// ----------------------------------------------------------------------------

	//Step-4: Time loop starts....
    size_counter = size_counter+1;
	clock_t s,s_prev,ds;
	s = clock();

    while(t<=t_end){

		fe = fe_apply_bc_load(fe,t);

		t_half = 0.5*(t+t+dT); // Here t = t+dT so 2t-dT = 2t + 2dT - dT
		t = t+dT;

		V_half = V + (t_half - (t-dT))*A; // Updating node velocities
		// velocity boundary conditions should be enforced here.

		// Displacement Calculations
		U = U + dT*(V_half); // Calculating the new displacements from nodal velocites.
		U = fe_apply_bc_displacement(U,t); // Enforcing displacement BCs.
		UU = fe_concatenate_vector2matrix(UU,U,2);

		F_net = fe_getforce(nodes,elements,ndof,U,V,fe,size_counter,nodes_truss,elements_truss); // Calculating the force term.
		dT = reduction * fe_getTimeStep(nodes, elements, ndof, U, V, fe);
		if(dT>dt_min){
			dT = dt_min;
		}


		// Acceleration calculations
		A = mm.inverse()*(F_net); // Calculating the new accelerations from total nodal forces.
		AA = fe_concatenate_vector2matrix(AA,A,2);

		// Velocity calculations
		V = V_half + (t - t_half)*A; // Calculating the new velocities.
		VV = fe_concatenate_vector2matrix(VV,V,2);

		/* Calculating the internal energy terms */
		double old_int_energy = W_int(size_counter-1);
		U_curr = U;
		VectorXd del_U = U_curr - U_prev;
		double new_int_energy = old_int_energy + 0.5*(del_U.dot(fi_prev + fi_curr));
		fi_prev = fi_curr;
		U_prev = U_curr;
		W_int.conservativeResize(W_int.size()+1);
		W_int(size_counter) = new_int_energy;

		/* Calculating the external energy terms */
		fe_curr = fe ;
		double old_ext_energy = W_ext(size_counter-1);
		double new_ext_energy = old_ext_energy + 0.5*(del_U.dot(fe_prev + fe_curr));
		fe_prev = fe_curr;
		W_ext.conservativeResize(W_ext.size()+1);
		W_ext(size_counter) = new_ext_energy;

		/* Calculating the kinetic energy */
		double kin_energy = 0.5*(V.transpose())*(mm)*(V);
		W_kin.conservativeResize(W_kin.size()+1);
		W_kin(size_counter) = kin_energy;

		/* Calculating the total energy of the system */
		W_tot.conservativeResize(W_tot.size()+1);
		W_tot(size_counter) = std::abs((W_kin(size_counter) + W_int(size_counter) - W_ext(size_counter)));

		if(W_tot(size_counter)>eps_energy){
			std::cout << "**********************************************" << std::endl;
			std::cout << "ALERT: INSTABILITIES IN THE SYSTEM DETECTED \n BASED ON THE ENERGY BALANCE CHECK \n";
			std::cout << "**********************************************" << std::endl;
			// std::exit(-1);
		}


		/** Writing the output to VTK files */
		// updated_nodes = fe_updateNodes(nodes,U);

		// std::cout<<"Z Strain: "<<(U(26)/2)<<"\n";
		// std::cout << "Z displacement: "<<(U(14))<<"\n";

    	U_host = U;
		V_host = V;
		A_host = A;

		if(t >= (time_temp_1 * (output_temp_1))){
			fe_vtuWrite("eem_matrix",size_counter,mesh[0]);
			fe_vtuWrite("eem_truss",size_counter,mesh[1]);
			time_temp_1 = time_temp_1 + 1;
			//fe_vtkWrite_host("eem_matrix",1,5,size_counter,nodes,elements);
			//fe_vtkWrite_truss("eem_truss",1,5,size_counter,nodes_truss,elements_truss);
		}

		s_prev = s;
		s = clock();
		ds = s - s_prev;

		std::cout << "Timestep Value = "<<dT<<" ; Current Time = "<< t << " ; Timestep Number = " << (size_counter) << "; CPU Time (Step) = " << ((float)ds/CLOCKS_PER_SEC) << "s \n";

		size_counter = size_counter+1;

	}


	std::string displacements = home_path+"results/nodal_displacements.txt";
	matrix2text(displacements.c_str(),UU,UU.cols());

	/*std::string velocities = home_path+"results/nodal_velocities.txt";
	matrix2text(velocities.c_str(),VV,VV.cols());

	std::string accelerations = home_path+"results/nodal_accelerations.txt";
	matrix2text(accelerations.c_str(),AA,AA.cols());
	*/

	std::string internal_energy = home_path + "results/internal_energy_system.txt";
	vector2text(internal_energy,W_int,10);

	std::string external_energy = home_path + "results/external_energy_system.txt";
	vector2text(external_energy,W_ext,10);

	std::string kinetic_energy = home_path + "results/kinetic_energy_system.txt";
	vector2text(kinetic_energy,W_kin,10);

	std::string total_energy = home_path + "results/total_energy_system.txt";
	vector2text(total_energy,W_tot,10);

}
