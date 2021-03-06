#include "functions.h"
using namespace Eigen;

double eps_energy = 0.01;
double area_truss = 5e-6;
double failure_time_step = 1e-8;

/*! \brief
 * This function carries out the explicit dynamic analysis of the FEM problem.
 */

void
fe_mainEXPLICIT()
{

    if (num_meshes == 0) {
        std::cout << "No meshes included - Simulation is not possible !! " << "\n";
        std::exit(-1);
    }

    for (int i = 0; i < num_meshes; i++) {
        mesh[i].preprocessMesh();
    }

    // Following variables - Only for Hex Element
    int nnode = mesh[0].getNumNodes();          // number of nodes
    int sdof  = nnode * ndof;          // system degrees of freedom

    // Initialization
    double dT             = dt_initial;
    double t              = t_start;
    int time_step_counter = 0; // time step count
    int plot_state_counter = 1;
    double output_temp_1  = ((double) (t_end / output_frequency));
    VectorXd A            = VectorXd::Zero(sdof); // Acceleration Vector
    VectorXd V            = VectorXd::Zero(sdof); // Velocity Vector
    VectorXd V_half       = VectorXd::Zero(sdof); // Velocity Vector at n+1/2
    VectorXd U            = VectorXd::Zero(sdof); // Displacement Vector
    VectorXd F_net        = VectorXd::Zero(sdof); // Total Nodal force vector
    VectorXd fe           = VectorXd::Zero(sdof); // External Nodal force vector
    VectorXd fe_prev      = VectorXd::Zero(sdof);

    VectorXd fr_prev      = VectorXd::Zero(sdof);
    VectorXd fr_curr      = VectorXd::Zero(sdof);
    VectorXd fi_prev      = VectorXd::Zero(sdof); // Internal nodal force vector at previous timestep
    VectorXd fi_curr      = VectorXd::Zero(sdof); // Internal Nodal force vector at current timestep
    VectorXd U_prev       = VectorXd::Zero(sdof);

    double energy_int_old = 0;
    double energy_int_new = 0;
    double energy_ext_old = 0;
    double energy_ext_new = 0;
    double energy_kin     = 0;
    double energy_total   = 0;
    double energy_max     = 0;

    std::string internal_energy = home_path + "/" + "results/internal_energy_system.txt";
    std::string external_energy = home_path + "/" + "results/external_energy_system.txt";
    std::string kinetic_energy = home_path + "/" + "results/kinetic_energy_system.txt";
    std::string total_energy = home_path + "/" + "results/total_energy_system.txt";
    fe_energyWrite_new(internal_energy, external_energy, kinetic_energy, total_energy, plot_state_counter, t, energy_int_new, energy_ext_new, energy_kin, energy_total);

    // Loading Conditions
    fe_apply_bc_load(fe, t_start);

    // ----------------------------------------------------------------------------
    // Step-1: Calculate the mass matrix similar to that of belytschko.
    VectorXd m_system = VectorXd::Zero(sdof);
    fe_calculateMass(m_system, "direct_lumped");

    std::string mass = home_path + "/" + "results/system_mass.txt";
    new_vector2text(mass.c_str(), m_system, m_system.cols());

    // ----------------------------------------------------------------------------
    // Step-2: getforce step from Belytschko
    fe_getforce(F_net, ndof, U, fe, time_step_counter);

    mesh[0].readNodalKinematics(U, V, A);

    for (int i = 0; i < num_meshes; i++) {
        fe_vtuWrite(plot_state_counter - 1, t, mesh[i]);
    }

    dT = fe_getTimeStep();

    // ----------------------------------------------------------------------------
    // Step-3: Calculate accelerations
    fe_calculateAccln(A, m_system, F_net);
    U_prev = U;

    // ----------------------------------------------------------------------------
    // Step-4: Time loop starts....
    time_step_counter = time_step_counter + 1;
    clock_t s, s_prev, ds;
    s = clock();

    while (t < t_end) {


        if (((t + dT) >= t_end) && (t != t_end)) {
            dT = t_end - t;
            if (dT <= 0) {
                break;
            }
        }


        /** Update Loading Conditions - time dependent loading conditions */
        fe_apply_bc_load(fe, t);

        /** Steps - 4,5,6 and 7 from Belytschko Box 6.1 - Update time, velocity and displacements */
        fe_timeUpdate(U, V, V_half, A, t, dT, "newmark-beta-central-difference");

        /** Step - 8 from Belytschko Box 6.1 - Calculate net nodal force*/
        fe_getforce(F_net, ndof, U, fe, time_step_counter); // Calculating the force term.



        /** Step - 9 from Belytschko Box 6.1 - Calculate Accelerations */
        fe_calculateAccln(A, m_system, F_net); // Calculating the new accelerations from total nodal forces.
        fe_apply_bc_acceleration(A, t);

        /** Step- 10 from Belytschko Box 6.1 - Second Partial Update of Nodal Velocities */
        fe_timeUpdate_velocity(V, V_half, A, t, dT, "newmark-beta-central-difference");


        fi_curr = fe - F_net;
        fe_calculateFR(fr_curr, sdof, fi_curr, m_system, A);

        std::cout << "****************************************************\n";
        std::cout << "current time = " << t << std::endl;
        std::cout << "RF3 node #4 = " << fr_curr[14] << std::endl;
        std::cout << "RF3 node #5 = " << fr_curr[17] << std::endl;
        std::cout << "RF3 node #6 = " << fr_curr[20] << std::endl;
        std::cout << "RF3 node #7 = " << fr_curr[23] << std::endl;
        std::cout << "****************************************************\n";



        /** Step - 11 from Belytschko Box 6.1 - Calculating energies and Checking Energy Balance */
        fe_checkEnergies(U_prev, U, fi_prev, fi_curr, fe_prev, fe, fr_prev, fr_curr, m_system, V, energy_int_old, energy_int_new, energy_ext_old, energy_ext_new, energy_kin, energy_total, energy_max);

        mesh[0].readNodalKinematics(U, V, A);

        if (t >= (plot_state_counter * (output_temp_1))) {

            for (int i = 0; i < num_meshes; i++) {
                fe_vtuWrite(plot_state_counter, t, mesh[i]);
            }

            plot_state_counter = plot_state_counter + 1;

            std::cout << "-----------------------------------" << "\n";
            std::cout << " Timestep Value = " << std::setw(5) << std::scientific << std::setprecision(1) << dT
                      << "\n Current Time = " << std::setw(5) << std::setprecision(1) << t
                      << "\n Timestep Number = " << (time_step_counter)
                      << "\n Plot State = " << (plot_state_counter - 1)
                      << "\n CPU Time = " << std::setw(5) << std::setprecision(1)
                      << ((float) ds / CLOCKS_PER_SEC) << "s \n";
            std::cout << std::setw(5) << std::scientific << std::setprecision(5) << "Z Displacement: " << U(20) << "\n";
            // std::cout << std::setw(5)<<std::scientific<<std::setprecision(5) <<"Current Precise Time: " << t << "\n";
            // std::cout << std::setw(5)<<std::scientific<<std::setprecision(5) <<"Internal Energy: " << energy_int_new << "\n";
            // std::cout << std::setw(5)<<std::scientific<<std::setprecision(5) <<"External Work: " << energy_ext_new << "\n";
            // std::cout << std::setw(5)<<std::scientific<<std::setprecision(5) <<"Kinetic Energy: " << energy_kin << "\n";

            /*print current frame, current time, and energy to individual .txt files*/
            fe_energyWrite_append(internal_energy, external_energy, kinetic_energy, total_energy, plot_state_counter, t, energy_int_new, energy_ext_new, energy_kin, energy_total);
        }

        s_prev = s;
        s      = clock();
        ds     = s - s_prev;
        time_step_counter = time_step_counter + 1;

        dT = fe_getTimeStep();


    }
} // fe_mainEXPLICIT
