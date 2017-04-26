#include "functions.h"

using namespace Eigen;

VectorXd
fe_getForce_3d_normal(VectorXd u, VectorXd fext, int time_step_counter, int host_id)
{

    MatrixXd nodes_host    = mesh[host_id].getNewNodes();
    MatrixXi elements_host = mesh[host_id].getNewElements();

    // Variables - Mesh Details
    int nel   = elements_host.rows();
    int nnel  = (elements_host.cols() - 2);
    int nnode = nodes_host.rows();
    int sdof  = nnode * ndof;
    int edof  = nnel * ndof;

    VectorXd element_stress_host_local = VectorXd::Zero(nel * 9);
    VectorXd element_strain_host_local = VectorXd::Zero(nel * 9);

    // Element Data
    VectorXi nodes_local = VectorXi::Zero(nnel);
    VectorXd xcoord      = VectorXd::Zero(nnel);
    VectorXd ycoord      = VectorXd::Zero(nnel);
    VectorXd zcoord      = VectorXd::Zero(nnel);

    // Total nodal force vector for the system.
    VectorXd f_tot(sdof);

    f_tot = VectorXd::Zero(sdof);

    for (int i = 0; i < nel; i++) {
        for (int j = 0; j < nnel; j++) {
            //  int g = -1;
            //  for(int f=0;f<nnode;f++){
            //      if(elements(i,j+2)==nodes(f,0)){
            //          g = f;
            //          break;
            //      }
            //  }

            int g = elements_host(i, j + 2);
            nodes_local(j) = g;
            xcoord(j)      = nodes_host(g, 1);
            ycoord(j)      = nodes_host(g, 2);
            zcoord(j)      = nodes_host(g, 3);
        }

        VectorXd u_e = VectorXd::Zero(edof); // element displacements
        u_e = fe_gather(u, u_e, nodes_local, sdof);

        VectorXd f_ext_e = VectorXd::Zero(edof);
        f_ext_e = fe_gather(fext, f_ext_e, nodes_local, sdof); // element external nodal forces

        VectorXd f_int_e = VectorXd::Zero(edof);
        f_int_e = VectorXd::Zero(edof); // element internal nodal force

        VectorXd f_tot_e = VectorXd::Zero(edof);
        f_tot_e = VectorXd::Zero(edof); // element total nodal force

        int nglx = 2;
        int ngly = 2;
        int nglz = 2;

        MatrixXd disp_mat(6, edof);

        VectorXd dndr(nnel);
        VectorXd dnds(nnel);
        VectorXd dndt(nnel);
        VectorXd dndx(nnel);
        VectorXd dndy(nnel);
        VectorXd dndz(nnel);
        MatrixXd jacobian(ndof, ndof);
        MatrixXd invJacobian(ndof, ndof);

        // MatrixXd points_3d = guass_points_3d(nglx,ngly,nglz);
        // MatrixXd weights_3d = guass_weights_3d(ndof,nglx,ngly,nglz);

        VectorXd points  = guass_points(nglx);
        VectorXd weights = guass_weights(nglx);

        int node_counter = 0;

        if (time_step_counter != 0) { // if this is not the first time step the go into the loop
            for (int intx = 0; intx < nglx; intx++) {
                double x   = points(intx);
                double wtx = weights(intx);
                for (int inty = 0; inty < ngly; inty++) {
                    double y   = points(inty);
                    double wty = weights(inty);
                    for (int intz = 0; intz < nglz; intz++) {
                        double z   = points(intz);
                        double wtz = weights(intz);

                        // VectorXd dndr(edof);
                        dndr = fe_dndr_8(x, y, z);
                        // VectorXd dnds(edof);
                        dnds = fe_dnds_8(x, y, z);
                        // VectorXd dndt(edof);
                        dndt = fe_dndt_8(x, y, z);

                        jacobian = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
                        double detJacobian = jacobian.determinant();
                        invJacobian = jacobian.inverse();

                        dndx = fe_dndx_8(nnel, dndr, dnds, dndt, invJacobian);
                        dndy = fe_dndy_8(nnel, dndr, dnds, dndt, invJacobian);
                        dndz = fe_dndz_8(nnel, dndr, dnds, dndt, invJacobian);

                        disp_mat = fe_strDispMatrix_totalLagrangian(edof, nnel, dndx, dndy, dndz, u_e);
                        // disp_mat = fe_strDispMatrix(edof, nnel, dndx, dndy, dndz);

                        VectorXd sigma_e = VectorXd::Zero(6);
                        sigma_e = fe_stressUpdate(dndx, dndy, dndz, disp_mat, u_e, elements_host(i, 1), 0);

                        // f_int_e = f_int_e + ((disp_mat.transpose())*sigma_e*wtx*wty*wtz*detJacobian); (previous correct)
                        // std::cout<<k<<std::endl;

                        f_int_e = f_int_e + ((disp_mat.transpose()) * sigma_e * wtx * wty * wtz * detJacobian);
                    }
                }
            }

            element_stress_host_local.segment<9>(i) = fe_calCentroidStress_3d(nnel, xcoord, ycoord, zcoord, u_e, elements_host(
                        i,
                        1));
            element_strain_host_local.segment<9>(i) = fe_calCentroidStrain_3d(nnel, xcoord, ycoord, zcoord, u_e);

            // Host element internal nodal forces

            // TRUSS ANALYSIS - WORKS FOR ONLY SINGLE ELEMENT PROBLEMS
            // truss element nodal force

            /*int truss_intg = 2;
             * VectorXd truss_intg_points = guass_points(2);
             * VectorXd truss_intg_weights = guass_weights(2);
             *
             * for(int q=0;q<truss_intg;q++){
             * double t = truss_intg_points(q);
             * double wtt = truss_intg_weights(q);
             * dndr = fe_dndr_8(0,0,t);
             * dnds = fe_dnds_8(0,0,t);
             * dndt = fe_dndt_8(0,0,t);
             * jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
             * invJacobian = jacobian.inverse();
             * dndx = fe_dndx_8(nnel, dndr, dnds, dndt, invJacobian);
             * dndy = fe_dndy_8(nnel, dndr, dnds, dndt, invJacobian);
             * dndz = fe_dndz_8(nnel, dndr, dnds, dndt, invJacobian);
             * disp_mat = fe_strDispMatrix_totalLagrangian(edof,nnel,dndx,dndy,dndz,u_e);
             * VectorXd sigma_truss = VectorXd::Zero(6);
             * sigma_truss = fe_stressUpdate_1d(dndx,dndy,dndz,u_e,elements_truss(0,1),nodes_truss);
             * VectorXd f_int_truss = (disp_mat.transpose()*sigma_truss*wtt*0.5*0.5);
             *
             * VectorXd sigma_correction = fe_stressUpdate_1d(dndx,dndy,dndz,u_e,elements(i,1),nodes_truss);
             * VectorXd f_int_correction = (disp_mat.transpose()*sigma_correction*wtt*0.5*0.5);
             *
             * f_int_e = f_int_e + f_int_truss - f_int_correction; // no Volume Redundancy
             * // f_int_e = (0.5*f_int_e) + f_int_truss; // no Volume Redundancy using volume fractions
             *
             * // f_int_e = f_int_e + f_int_truss ; // With Volume Redundancy
             * }
             *
             * // Calculating Truss Element Stress and Strains
             *
             * dndr = fe_dndr_8(0,0,0);
             * dnds = fe_dnds_8(0,0,0);
             * dndt = fe_dndt_8(0,0,0);
             * jacobian = fe_calJacobian(ndof,nnel,dndr,dnds,dndt,xcoord,ycoord,zcoord);
             * invJacobian = jacobian.inverse();
             * dndx = fe_dndx_8(nnel, dndr, dnds, dndt, invJacobian);
             * dndy = fe_dndy_8(nnel, dndr, dnds, dndt, invJacobian);
             * dndz = fe_dndz_8(nnel, dndr, dnds, dndt, invJacobian);
             * disp_mat = fe_strDispMatrix_totalLagrangian(edof,nnel,dndx,dndy,dndz,u_e);
             *
             * // Strain - Calculation
             * F = fe_calDefGrad(dndx,dndy,dndz,u_e);
             * MatrixXd C = F.transpose()*F;
             * VectorXd dir_truss = VectorXd::Zero(3);
             * dir_truss << (nodes_truss(1,1)-nodes_truss(0,1)) , (nodes_truss(1,2)-nodes_truss(0,2)), (nodes_truss(1,3)-nodes_truss(0,3));
             * VectorXd deformed_truss = C*dir_truss;
             * double lambda_tmp = dir_truss.dot(deformed_truss);
             * double lambda = sqrt(lambda_tmp);
             * element_strain_truss(i,0) = lambda-1;
             * element_strain_truss(i,1) = 0;
             * element_strain_truss(i,2) = 0;
             * element_strain_truss(i,3) = 0;
             * element_strain_truss(i,4) = 0;
             * element_strain_truss(i,5) = 0;
             * element_strain_truss(i,6) = 0;
             * element_strain_truss(i,7) = 0;
             * element_strain_truss(i,8) = 0;
             *
             * // Stress- Calculation
             * std::string model;
             * model = fe_get_model(elements_truss(0,1));
             * double sigma = 0;
             * if(model=="simple_elastic"){
             *  double E = fe_get_mats(elements_truss(0,1),1);
             *  sigma = E*log(lambda);
             * }
             * element_stress_truss(i,0) = sigma;
             * element_stress_truss(i,1) = 0;
             * element_stress_truss(i,2) = 0;
             * element_stress_truss(i,3) = 0;
             * element_stress_truss(i,4) = 0;
             * element_stress_truss(i,5) = 0;
             * element_stress_truss(i,6) = 0;
             * element_stress_truss(i,7) = 0;
             * element_stress_truss(i,8) = 0;*/
        }
        f_tot_e = f_ext_e - f_int_e;
        f_tot   = fe_scatter(f_tot, f_tot_e, nodes_local, sdof);
    }

    mesh[host_id].readElementStressStrain(element_stress_host_local, element_strain_host_local);
    return f_tot;
} // fe_getForce_3d_normal
