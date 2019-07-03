#include <iostream>
#include "Point.hpp"

int test_Point(){
    std::vector<double> co(2,1);
    Point p1(co,2);
    Point p2(2,3);
    std::cout<<p1[0]<<"\t"<<p1[1]<<"\t"<<p1.area()<<std::endl;
    std::cout<<p2[0]<<"\t"<<p2.area()<<std::endl;
    p2.set_a(10);
    p2.set_co(0,5);
    std::cout<<p2[0]<<"\t"<<p2.area()<<std::endl;
    (p1*p2).print();
}

int main(){
    test_Point();
    return 1;
}
