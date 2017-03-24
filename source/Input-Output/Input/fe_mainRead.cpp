#include<iomanip>
#include "functions.h"
using namespace Eigen;

void fe_mainRead(std::string file){

  num_meshes = 0;
  int num_meshes_counter = 0;

  material_types = 0;
  int material_types_counter = 0;

  std::ifstream myfile(file.c_str());
  std::string line;
  std::getline(myfile,line);

  while(line!="*END"){
    if(line[0]=='*'){
        if(line=="*MESH"){
          num_meshes = num_meshes+1;
        }
        if(line=="*MATERIAL"){
          material_types = material_types+1;
        }
    }
    std::getline(myfile,line);
  }

  mesh = new Mesh[num_meshes];
  mat = new Materials[material_types];

  std::ifstream myfile1(file.c_str());
  std::getline(myfile1,line);



  while(line!="*END"){

    if(line[0]=='*'){

        if(line=="*MESH"){
          std::getline(myfile1,line);

          while(line!="*END_MESH"){
            MatrixXd nodes;
            MatrixXd elements;

            if(line=="*NODES"){
              int rows = 0;
              int cols = 0;
              myfile1 >> rows;
              myfile1 >> cols;
              nodes = MatrixXd::Zero(rows,cols+1);
              for(int i=0;i<rows;i++){
                for(int j=0;j<(cols+1);j++){
                  myfile1 >> nodes(i,j);
                }
              }
              myfile1 >> line;
            }

            if(line=="*ELEMENTS"){
              int rows = 0;
              int cols = 0;
              myfile1 >> rows;
              myfile1 >> cols;
              elements = MatrixXd::Zero(rows,cols+2);
              for(int i=0;i<rows;i++){
                for(int j=0;j<(cols+2);j++){
                  myfile1 >> elements(i,j);
                }
              }
              myfile1 >> line;
            }

            if(nodes.rows()!=0 && elements.rows()!=0){
              mesh[num_meshes_counter].readMesh(nodes,elements);
              num_meshes_counter = num_meshes_counter + 1;
            }
        }
      }

        if(line=="*MATERIAL"){
            int mat_id;
            std::string mat_model;
            VectorXd mat_properties = VectorXd::Zero(100);
            myfile1 >> mat_id;
            myfile1 >> mat_model;
            int i = 0;
            std::getline(myfile1,line);
          while(line!="*END_MATERIAL"){
            myfile1.seekg(-(line.length()),std::ios::cur);
            myfile1 >> mat_properties(i);
            myfile1 >> line;
            i = i+1;
          }
          mat[material_types_counter].readMats(mat_id,mat_model,mat_properties);
          material_types_counter = material_types_counter + 1;
        }
    }
    std::getline(myfile1,line);
  }

}
