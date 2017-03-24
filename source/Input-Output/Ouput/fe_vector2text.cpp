#include<iomanip>
#include"functions.h"
using namespace Eigen;

void vector2text(std::string name, VectorXd vector, int width){
	if(width==0){
			width=3;
		}
	std::ofstream myfile(name.c_str());

	for(int i=0;i<vector.size();i++){
		myfile << std::setw(width) << vector(i) << "\n";
	}
	myfile.close();
}
