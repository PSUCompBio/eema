#include "functions.h"

using namespace Eigen;

void fe_vtuWrite(std::string output, int time_step, Mesh mesh1){

    MatrixXd nodes = mesh1.getNewNodes();
    MatrixXi elements = mesh1.getNewElements();

    std::string name;

    std::ostringstream ss;
    ss << time_step;

    name = home_path + "results/vtu/" + output + "_" + ss.str() + ".vtu";

    std::ofstream myfile(name.c_str());

    myfile << "<?xml version=\"1.0\"?>\n";
    myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    myfile << "\t<UnstructuredGrid>\n";
    myfile << "\t\t<Piece NumberOfPoints=\"" << nodes.rows() << "\" NumberOfCells=\"" << elements.rows() << "\">\n";

    /** Points Data */
    myfile << "\t\t\t<Points>\n";
    myfile << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i=0;i<nodes.rows();i++){
        	for(int j=1;j<nodes.cols();j++){
				myfile<<"\t\t\t\t\t"<<std::scientific<<std::setprecision(10)<<nodes(i,j)<<" ";
        	}
        	myfile<<"\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t</Points>\n";

    /** Cell Data */
    myfile << "\t\t\t<Cells>\n";
    myfile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(int i=0;i<elements.rows();i++){
                myfile << "\t\t\t\t\t" ;
        	    for(int j=2;j<elements.cols();j++){
				    myfile<<elements(i,j)<<"\t";
                }
    }
    myfile<<"\n";
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for(int i=0;i<elements.rows();i++){
            myfile << "\t\t\t\t\t" << (elements.row(i).cols()-2) << "\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";

    myfile << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for(int i=0;i<elements.rows();i++){
        if((elements.row(i).cols()-2)==8) myfile << "\t\t\t\t\t" << "12" << "\n";
        if((elements.row(i).cols()-2)==2) myfile << "\t\t\t\t\t" << "3" << "\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";

    myfile << "\t\t\t</Cells>\n";

    /** Point Vector Data - Displacements */
    myfile << "\t\t\t<PointData>\n";
    myfile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Displacement\" NumberOfComponents=\"3\" ComponentName0=\"X\" ComponentName1=\"Y\" ComponentName2=\"Z\" format=\"ascii\">\n";
    int num = 0;
    for(int i=0;i<nodes.rows();i++){
        	myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << U_host(num) << " " << U_host(num+1) << " " << U_host(num+2) << " \n";
        	num = num+3;
    }
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t</PointData>\n";

    /** Cell Data - Stresses and Strains */
    myfile << "\t\t\t<CellData>\n";
    myfile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Stress (Cauchy)\" NumberOfComponents=\"9\" ComponentName0=\"XX\" ComponentName1=\"XY\" ComponentName2=\"XZ\" ComponentName3=\"YX\" ComponentName4=\"YY\" ComponentName5=\"YZ\" ComponentName6=\"ZX\" ComponentName7=\"ZY\" ComponentName8=\"ZZ\" format=\"ascii\">\n";
    for(int i=0;i<element_stress_host.rows();i++){
       	    myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_stress_host(i,0) << " " << element_stress_host(i,1) << " " << element_stress_host(i,2) << "\n";
       	    myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_stress_host(i,3) << " " << element_stress_host(i,4) << " " << element_stress_host(i,5) << "\n";
        	myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_stress_host(i,6) << " " << element_stress_host(i,7) << " " << element_stress_host(i,8) << "\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Strain (Logarithmic)\" NumberOfComponents=\"9\" ComponentName0=\"XX\" ComponentName1=\"XY\" ComponentName2=\"XZ\" ComponentName3=\"YX\" ComponentName4=\"YY\" ComponentName5=\"YZ\" ComponentName6=\"ZX\" ComponentName7=\"ZY\" ComponentName8=\"ZZ\" format=\"ascii\">\n";
    for(int i=0;i<element_stress_host.rows();i++){
       	    myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_strain_host(i,0) << " " << element_strain_host(i,1) << " " << element_strain_host(i,2) << "\n";
       	    myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_strain_host(i,3) << " " << element_strain_host(i,4) << " " << element_strain_host(i,5) << "\n";
        	myfile << "\t\t\t\t\t" << std::scientific << std::setprecision(10) << element_strain_host(i,6) << " " << element_strain_host(i,7) << " " << element_strain_host(i,8) << "\n";
    }
    myfile << "\t\t\t\t</DataArray>\n";
    myfile << "\t\t\t</CellData>\n"; 

    myfile << "\t\t</Piece>\n";
    myfile << "\t</UnstructuredGrid>\n";
    myfile << "</VTKFile>\n";
}

