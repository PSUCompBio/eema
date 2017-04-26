#include "functions.h"

using namespace Eigen;

void fe_calculateMass(std::string type) {

  if (type == "direct_lumped" && embedded_constraint == true) { // Embedded Mass
    for (int i = 0; i < num_constraints; ++i) {
      if (cons[i].getName() == "embedded") {
        std::string host = cons[i].get_EmbedMaster();
        std::string slave = cons[i].get_EmbedSlave();
        bool correct_vr = cons[i].get_EmbedAddressVR();
        int host_id, embed_id;
        for (int i = 0; i < num_meshes; ++i) {
          std::string name = mesh[i].getName();
          if (name == host) {
            host_id = i;
          }
          if (name == slave) {
            embed_id = i;
          }
        }
        fe_calculateMassDirectLumped_embed(host_id, embed_id, correct_vr);
      }
    }
  }
  else if (type == "direct_lumped" && embedded_constraint != true) {
    fe_calculateMassDirectLumped(0);
  }
  else {
    fe_calculateMassDirectLumped(0);
  }

}

void fe_calculateMassDirectLumped(int mesh_id) {

  MatrixXd nodes = mesh[mesh_id].getNewNodes();
  MatrixXi elements = mesh[mesh_id].getNewElements();

  int nel = elements.rows(); /*! number of elements */
  int nnel = (elements.cols() - 2); // number of nodes per element
  int nnode = nodes.rows(); // number of nodes
  int sdof = nnode * ndof; // system degrees of freedom
  int edof = nnel * ndof; // element degrees of freedom

  m_system = VectorXd::Zero(sdof);

  for (int i = 0; i < nel; i++) {

    VectorXd m_element = VectorXd::Zero(edof); // mass of hex elements
    m_element = fe_massLumped(nodes, elements.row(i));

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

    fe_scatterMass(m_element, elements.block<1, 8>(i, 2), sdof);
  }

}

void fe_calculateMassDirectLumped_embed(int host_id, int embed_id, bool address_vr) {

  std::cout << "Entered Mass - Embedded" << "\n";

  MatrixXd nodes_host = mesh[host_id].getNewNodes();
  MatrixXi elements_host = mesh[host_id].getNewElements();
  MatrixXd nodes_embed = mesh[embed_id].getNewNodes();
  MatrixXi elements_embed = mesh[embed_id].getNewElements();

  int nel = elements_host.rows(); /*! number of elements */
  int nnel = (elements_host.cols() - 2); // number of nodes per element
  int nnode = nodes_host.rows(); // number of nodes
  int sdof = nnode * ndof; // system degrees of freedom
  int edof = nnel * ndof; // element degrees of freedom

  m_system = VectorXd::Zero(sdof);

  VectorXd xcoord_embed(elements_embed.cols() - 2);
  VectorXd ycoord_embed(elements_embed.cols() - 2);
  VectorXd zcoord_embed(elements_embed.cols() - 2);

  for (int i = 0; i < nel; i++) {

    VectorXd m_element = VectorXd::Zero(edof); // mass of hex elements
    m_element = fe_massLumped(nodes_host, elements_host.row(i));

    for (int fib = 0; fib < elements_embed.rows(); fib++) {
      for (int j = 0; j < elements_embed.cols() - 2; j++) {
        int g = elements_embed(fib, j + 2);
        xcoord_embed(j) = nodes_embed(g, 1);
        ycoord_embed(j) = nodes_embed(g, 2);
        zcoord_embed(j) = nodes_embed(g, 3);
      }

      double rho_truss = fe_get_mats(elements_embed(fib, 1), 0);
      double rho_host = fe_get_mats(elements_host(i, 1), 0);

      double m_truss = fe_calVolume(xcoord_embed, ycoord_embed, zcoord_embed) * area_truss * rho_truss;
      double m_correction = fe_calVolume(xcoord_embed, ycoord_embed, zcoord_embed) * area_truss * rho_host;

      for (int k = 0; k < edof; k++) {
        m_element(k) = m_element(k) + (m_truss / nnel);
        if (address_vr == true) {
          m_element(k) = m_element(k) - (m_correction / nnel);
        }
      }
    }

    fe_scatterMass(m_element, elements_host.block<1, 8>(i, 2), sdof);
  }

}
