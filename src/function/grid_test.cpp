# include <iostream>

# include "grid.hpp"

int test_grid(){
    std::vector<double> w(10,0);
    std::vector<double> q(20,0);
    for(int i=0;i<w.size();i++) w[i]=i*i;
    for(int i=0;i<q.size();i++) q[i]=i;
    std::vector<std::vector<double>> g;
    g.push_back(w);    g.push_back(q);
    Grid g1(g);
    std::cout<<g1.size()<<"\t"<<g1.dimension()<<std::endl;
    for(int i=0;i<g1.size();i++)
	std::cout<<i<<"\t"<<g1.coordinate(i,0)<<"\t"<<g1.coordinate(i,1)
		 <<"\t"<<g1.point(i)[0]<<"\t"<<g1.point(i)[1]<<std::endl;
    return 1;
}

int main(){
    test_grid();
    return 0;
}
