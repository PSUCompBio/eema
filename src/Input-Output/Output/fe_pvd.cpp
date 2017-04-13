#include "functions.h"

using namespace Eigen;

void fe_pvdNew(std::string output, int time_step, double time){

    /** Output File Name */
    std::string name;
    name = home_path + "/" + "results/vtu/" + output + ".pvd";
    std::ofstream myfile(name.c_str());

    /** Write Header */
    myfile << "<?xml version=\"1.0\"?>\n";
    myfile << "<VTKFile type=\"Collection\" version=\"0.1\"\n";
    myfile << "\t\t\t\tbyte_order=\"LittleEndian\"\n";
    myfile << "\t\t\t\tcompressor=\"vtkZLibDataCompressor\">\n";
    myfile << "\t\t<collection>\n";

    /** Write First Line of Data */
    myfile << "\t\t\t<DataSet timestep= " << time << "\" group=\"\" part=\"0\"\n";
    std::ostringstream ss;
    ss << time_step;
    std::string vtuFileName;
    vtuFileName = home_path + "/" + "results/vtu/" + output + "_" + ss.str() + ".vtu";
    myfile << "\t\t\t\t\t\t<file= \"" << name << "\"/>\n";
}
