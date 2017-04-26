#include "functions.h"

using namespace Eigen;

VectorXd
fe_getForce_3d_embed(VectorXd u, VectorXd fext, int time_step_counter, int host_id, int embed_id, bool address_vr)
{

    MatrixXd nodes_host     = mesh[host_id].getNewNodes();
    MatrixXi elements_host  = mesh[host_id].getNewElements();
    MatrixXd nodes_embed    = mesh[embed_id].getNewNodes();
    MatrixXi elements_embed = mesh[embed_id].getNewElements();

    // Variables - Mesh Details
    int nel   = elements_host.rows();
    int nnel  = (elements_host.cols() - 2);
    int nnode = nodes_host.rows();
    int sdof  = nnode * ndof;
    int edof  = nnel * ndof;

    VectorXd element_stress_host_local  = VectorXd::Zero(nel * 9);
    VectorXd element_strain_host_local  = VectorXd::Zero(nel * 9);
    VectorXd element_stress_embed_local = VectorXd::Zero(elements_embed.rows() * 9);
    VectorXd element_strain_embed_local = VectorXd::Zero(elements_embed.rows() * 9);

    // Element Data
    VectorXi nodes_local = VectorXi::Zero(nnel);
    VectorXd xcoord      = VectorXd::Zero(nnel);
    VectorXd ycoord      = VectorXd::Zero(nnel);
    VectorXd zcoord      = VectorXd::Zero(nnel);

    VectorXd nodes_local_embed = VectorXd::Zero(elements_embed.cols() - 2);
    VectorXd xcoord_embed      = VectorXd::Zero(elements_embed.cols() - 2);
    VectorXd ycoord_embed      = VectorXd::Zero(elements_embed.cols() - 2);
    VectorXd zcoord_embed      = VectorXd::Zero(elements_embed.cols() - 2);
    VectorXd u_embed       = VectorXd::Zero((nodes_embed.rows()) * ndof);
    VectorXd u_embed_local = VectorXd::Zero((elements_embed.cols() - 2) * ndof);
    VectorXd v_embed       = VectorXd::Zero((nodes_embed.rows()) * ndof);
    VectorXd a_embed       = VectorXd::Zero((nodes_embed.rows()) * ndof);

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

        /*int counter = 0;
         * for (int var1 = 0; var1 < xcoord.size(); var1++) {
         *  xcoord(var1) = xcoord(var1) + u_e(counter);
         *  ycoord(var1) = ycoord(var1) + u_e(counter + 1);
         *  zcoord(var1) = zcoord(var1) + u_e(counter + 2);
         *  counter      = counter + 3;
         * }*/

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

            element_stress_host_local.segment<9>(i * 9) = fe_calCentroidStress_3d(nnel, xcoord, ycoord, zcoord, u_e, elements_host(
                        i,
                        1));
            element_strain_host_local.segment<9>(i * 9) = fe_calCentroidStrain_3d(nnel, xcoord, ycoord, zcoord, u_e);

            /*
             *  EMBEDDED FIBER ANALYSIS - STARTS FROM HERE
             *
             *  Step-1: Find the number of truss elements present this single host element using the mapping vector.
             *  Step-2: For that number of truss elements and using the mapping vector, one should perform the analysis and update the local element stress/strain information.
             */

            int ngl_embed          = 2;
            VectorXd points_embed  = guass_points(ngl_embed);
            VectorXd weights_embed = guass_weights(ngl_embed);

            for (int fib = 0; fib < elements_embed.rows(); fib++) {
                for (int j = 0; j < (elements_embed.cols() - 2); j++) {
                    int g = elements_embed(fib, j + 2);
                    nodes_local_embed(j) = g;
                    xcoord_embed(j)      = nodes_embed(g, 1);
                    ycoord_embed(j)      = nodes_embed(g, 2);
                    zcoord_embed(j)      = nodes_embed(g, 3);

                    VectorXd nat_cooord_nodes = VectorXd::Zero(3);
                    nat_cooord_nodes(0) = xcoord_embed(j);
                    nat_cooord_nodes(1) = ycoord_embed(j);
                    nat_cooord_nodes(2) = zcoord_embed(j);
                    VectorXd iso_coord_nodes = fe_newtonRhapson(nat_cooord_nodes, xcoord, ycoord, zcoord);
                    VectorXd shapes = fe_shapes_8(iso_coord_nodes(0), iso_coord_nodes(1), iso_coord_nodes(2));
                    MatrixXd shape_mat_embed = fe_shapeMatrix(edof, nnel, shapes);
                    u_embed.segment<3>(g * ndof)   = (shape_mat_embed * u_e);
                    u_embed_local.segment<3>(j * ndof) = u_embed.segment<3>(g * ndof);
                }

                //std::cout << "Embedded Displacements: \n" << u_embed << "\n";
                //std::cout << "Local Embed Displacement: \n" << u_embed_local << "\n";


                double length_embed = fe_calVolume(xcoord_embed, ycoord_embed, zcoord_embed);

                for (int embed_intg = 0; embed_intg < ngl_embed; embed_intg++) {

                    VectorXd local_intg_points = fe_findIntgPoints_1d(xcoord_embed, ycoord_embed, zcoord_embed, points_embed(embed_intg), length_embed);
                    VectorXd global_intg_poins = fe_newtonRhapson(local_intg_points, xcoord, ycoord, zcoord);

                    VectorXd global_intg_points(3);
                    global_intg_points(0) = 0;
                    global_intg_points(1) = 0;
                    global_intg_points(2) = points_embed(embed_intg);

                    double wtt = weights_embed(embed_intg);

                    dndr        = fe_dndr_8(global_intg_points(0), global_intg_points(1), global_intg_points(2));
                    dnds        = fe_dnds_8(global_intg_points(0), global_intg_points(1), global_intg_points(2));
                    dndt        = fe_dndt_8(global_intg_points(0), global_intg_points(1), global_intg_points(2));
                    jacobian    = fe_calJacobian(ndof, nnel, dndr, dnds, dndt, xcoord, ycoord, zcoord);
                    invJacobian = jacobian.inverse();
                    dndx        = fe_dndx_8(nnel, dndr, dnds, dndt, invJacobian);
                    dndy        = fe_dndy_8(nnel, dndr, dnds, dndt, invJacobian);
                    dndz        = fe_dndz_8(nnel, dndr, dnds, dndt, invJacobian);
                    disp_mat    = fe_strDispMatrix_totalLagrangian(edof, nnel, dndx, dndy, dndz, u_e);

                    VectorXd sigma_truss = VectorXd::Zero(6);

                    // Procedure - 1: (Same Deformation Gradient - Because No Slip)
                    sigma_truss = fe_stressUpdate(dndx, dndy, dndz, disp_mat, u_e, elements_embed(fib, 1), 0);
                    VectorXd f_int_truss = (disp_mat.transpose() * sigma_truss * wtt * (length_embed / 2) * area_truss);

                    // Procedure - 2: (Same Displacements - Same Deformation Gradient - Transformation Matrix Inside)
                    /* sigma_truss = fe_stressUpdate_1d(elements_embed(fib, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, dndx, dndy, dndz, u_e, 0);
                    VectorXd f_int_truss = (disp_mat.transpose() * sigma_truss * wtt * (length_embed / 2) * area_truss); */

                    // Procedure - 3: (Same Displacements - Same Deformation Gradient - Transformation Matrix Outside)
                    /* MatrixXd stress_transformation_mat = fe_calTransformation(xcoord_embed, ycoord_embed, zcoord_embed, 1);
                    sigma_truss = fe_stressUpdate_1d(elements_embed(fib, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, 0);
                    VectorXd f_int_truss = (disp_mat.transpose() * stress_transformation_mat.transpose() * sigma_truss * wtt * (length_embed / 2) * area_truss);*/

                    f_int_e = f_int_e + f_int_truss;

                    if (address_vr == true) {
                        VectorXd sigma_correction = VectorXd::Zero(6);

                        sigma_correction = fe_stressUpdate(dndx, dndy, dndz, disp_mat, u_e, elements_embed(fib, 1), 0);
                        VectorXd f_int_correction = (disp_mat.transpose() * sigma_correction * wtt * (length_embed / 2) * area_truss);

                        /*sigma_correction = fe_stressUpdate_1d(elements_host(i, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, 0);
                        VectorXd f_int_correction = (disp_mat.transpose() * stress_transformation_mat.transpose() * sigma_correction * wtt * (length_embed / 2) * area_truss );*/

                        f_int_e = f_int_e - f_int_correction;
                    }

                    // f_int_e = f_int_e + f_int_truss - f_int_correction; // no Volume Redundancy
                    // f_int_e = (0.5*f_int_e) + f_int_truss; // no Volume Redundancy using volume fractions
                    // With Volume Redundancy
                }

                element_strain_embed_local.segment<9>(fib * 9) = fe_calCentroidStrain_embed_3d(u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed);
                element_stress_embed_local.segment<9>(fib * 9) = fe_calCentroidStress_embed_3d(elements_embed(fib, 1), u_embed_local, xcoord_embed, ycoord_embed, zcoord_embed, length_embed, xcoord, ycoord, zcoord);
            }
        }

        //std::cout << "Host Final - nodal force: " << f_int_e.maxCoeff() << "\n";

        f_tot_e = f_ext_e - f_int_e;
        f_tot   = fe_scatter(f_tot, f_tot_e, nodes_local, sdof);
    }

    mesh[host_id].readElementStressStrain(element_stress_host_local, element_strain_host_local);
    mesh[embed_id].readElementStressStrain(element_stress_embed_local, element_strain_embed_local);
    mesh[embed_id].readNodalKinematics(u_embed, v_embed, a_embed);

    return f_tot;
} // fe_getForce_3d_embed
