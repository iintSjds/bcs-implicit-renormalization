# include <iostream>

# include "function.hpp"


int test_function(){
    std::vector<double> w(5,0);
    std::vector<double> q(5,0);
    for(int i=0;i<w.size();i++) w[i]=i*0.5;
    for(int i=0;i<q.size();i++) q[i]=i*0.8;
    std::vector<std::vector<double>> g;
    g.push_back(w);    g.push_back(q);
    Grid g1(g);
    Function f1(g1);
    std::cout<<f1.size()<<"\t"<<f1.dimension()<<std::endl;
    for(int i=0;i<f1.size();i++)
	std::cout<<i<<"\t"<<f1[i]
		 <<"\t"<<f1.point(i)[0]<<"\t"<<f1.point(i)[1]<<std::endl;
    return 1;
}

int main(){
    test_function();
    return 0;
}
