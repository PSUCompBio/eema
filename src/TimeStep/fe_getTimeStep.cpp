#include "functions.h"

using namespace Eigen;

/** For all elements -- this function calculates the minimum critical timestep */
double
fe_getTimeStep(MatrixXd nodes, MatrixXi elements, int ndof, VectorXd u, VectorXd v, VectorXd fext)
{
    double deltaT_crit = 0;

    int nel         = elements.rows();
    int nnel        = (elements.cols() - 2);
    int nnode       = nodes.rows();
    int sdof        = nnode * ndof;
    int edof        = nnel * ndof;
    VectorXd xcoord = VectorXd::Zero(nnel);
    VectorXd ycoord = VectorXd::Zero(nnel);
    VectorXd zcoord = VectorXd::Zero(nnel);

    VectorXd deltaT_element = VectorXd::Zero(nel);

    for (int i = 0; i < nel; i++) {
        VectorXi nodes_local = VectorXi::Zero(nnel);

        for (int j = 0; j < nnel; j++) {
            int g = elements(i, j + 2);
            xcoord(j) = nodes(g, 1);
            ycoord(j) = nodes(g, 2);
            zcoord(j) = nodes(g, 3);
        }

        VectorXd u_e(edof);
        u_e = fe_gather(u, u_e, elements.block<1, 8>(i, 2), sdof);

        deltaT_element(i) = fe_calTimeStep(xcoord, ycoord, zcoord, elements(i, 1), u_e); // reduction factor for time step added.
    }

    deltaT_crit = deltaT_element.minCoeff();

    return deltaT_crit;
} // fe_getTimeStep
