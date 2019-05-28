#include <iostream>
#include "Grid.hpp"

int test_Grid_1d(){
    Grid_1d g1(-1,1,10);
    for(int i=0;i<g1.size();i++){
	std::cout<<g1[i][0]<<"\t";
    }
    std::cout<<std::endl;
    std::vector<double> vd;
    for(int i=0;i<10;i++) vd.push_back(i);
    Grid_1d g2(vd);
    for(int i=0;i<g2.size();i++){
	std::cout<<g2[i][0]<<"\t";
    }
    std::cout<<std::endl;
    MeshGrid g3(g1,g2);
    for(int i=0;i<g3.lengths()[0];i++){
	for(int j=0;j<g3.lengths()[1];j++)
	    std::cout<<g3[i*g3.lengths()[1]+j][0]<<","<<g3[i*g3.lengths()[1]+j][1]<<"\t";
	std::cout<<std::endl;
    }
    std::cout<<std::endl;
    
}

int main(){
    test_Grid_1d();
    return 0;
}
