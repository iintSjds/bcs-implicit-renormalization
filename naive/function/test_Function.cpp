#include <iostream>
#include "Function.hpp"


double f(Point p){
    double result=0;
    for(int i=0;i<p.dimension();i++)
	result+=p[i];
    return result;
}

int test_Function(){
    Grid_1d g(0,10,10);
    TabFunc tf(g,f);
    std::cout<<"f="<<tf[1]<<std::endl;
    tf.at(1).print();
    BareFunc bf(g,f);
    std::cout<<"f="<<tf[1]<<std::endl;
    bf.at(1).print();
    return 1;
}

int main(){
    test_Function();
    return 0;
}
