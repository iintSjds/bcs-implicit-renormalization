#include "GG.hpp"
#include <iostream>

int test(){
    GG GG_test(1,1,3.14,0.0005,0.0314,100,2000);
    for(int i=0;i<GG_test.N0;i++) std::cout<<GG_test.freq0[i]<<"\t"<<GG_test.result0(i)<<"\t"<<GG_test.freq0[i]*GG_test.result0(i)<<std::endl;
    //for(int i=0;i<GG_test.N1;i++) std::cout<<GG_test.freq1[i]<<"\t"<<GG_test.result1(i)<<"\t"<<GG_test.freq1[i]*GG_test.result1(i)<<std::endl;    
    return 1;
}

int main(){
    test();
    return 0;
}
