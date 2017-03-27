#include "functions.h"

using namespace Eigen;

/** This function carries out all the required pre-processing steps for embedded element constraint */

VectorXd fe_embed_preprocessing(Mesh host,Mesh embed){

    fe_embed_preprocessing_length(host,embed);

}

VectorXd fe_embed_preprocessing_length(Mesh host, Mesh embed){

    MatrixXd nodes_host = host.getNewNodes();
    MatrixXd elements_host = host.getNewElements();
    MatrixXd nodes_embed = embed.getNewNodes();
    MatrixXd elements_embed = embed.getNewElements();
    double xcoord[(elements_host.cols()-2)];
    double ycoord[(elements_host.cols()-2)];
    double zcoord[(elements_host.cols()-2)];

    double min_length = 10000000;
    double lc_host;
    double lc_embed;

    /* Calculate the min length of the host system */
    for(int i=0;i<elements_host.rows();i++){
        for(int j=0;j<(elements_host.cols()-2);j++){
			int g = -1;
			for(int f=0;f<nodes_host.rows();f++){
				if(elements_host(i,j+2)==nodes_host(f,0)){
					g = f;
					break;
				}
			}
			xcoord[j] = nodes_host(g,1);
			ycoord[j] = nodes_host(g,2);
			zcoord[j] = nodes_host(g,3);
		}

        lc_host = fe_minElementLength(xcoord,ycoord,zcoord);

        if(lc_host<min_length){
            min_length = lc_host;
        }

        if(lc_host == 0){
        std::cout << "**********************************************" << std::endl;
        std::cout << "ZERO LENGTH HOST ELEMENT IN THE SYSTEM \n FIBER LENGTH PREPROCESSING NOT POSSIBLE" << etd::endl;
        std::cout << "**********************************************" << std::endl;
        std::exit(-1);
        }

    }

    /* Modify the fiber length so that every fiber's length is less than that of the system min length */
    for(int i=0;i<elements_embed.rows();i++){
        for(int j=0;j<(elements_embed.cols()-2);j++){
			int g = -1;
			for(int f=0;f<nodes_embed.rows();f++){
				if(elements_embed(i,j+2)==nodes_embed(f,0)){
					g = f;
					break;
				}
			}
			xcoord[j] = nodes_embed(g,1);
			ycoord[j] = nodes_embed(g,2);
			zcoord[j] = nodes_embed(g,3);
		}

        lc_embed = fe_maxElementLength(xcoord,ycoord,zcoord);

        if(lc_embed<min_length){
            min_length = lc_embed;
        }

        if(lc_embed == 0){
        std::cout << "**********************************************" << std::endl;
        std::cout << "ZERO LENGTH HOST ELEMENT IN THE SYSTEM \n FIBER LENGTH PREPROCESSING NOT POSSIBLE" << etd::endl;
        std::cout << "**********************************************" << std::endl;
        std::exit(-1);
        }
        
    }

    while (lc_embed > lc_host){

    }
    




    



}

VectorXd fe_embed_preprocessing_mapping(){
    
   // VectorXd embed_vector;
//
//    MatrixXd nodes1 = host.getNewNodes();
//    MatrixXd elements1 = host.getNewElements();
//    MatrixXd nodes2 = embed.getNewNodes();
//    MatrixXd elements2 = embed.getNewElements();
//
//    MatrixXd nodes_embed = MatrixXd::Zero(nodes2.rows(),(nodes2.cols()+1));
//
//    return embed_vector;
}