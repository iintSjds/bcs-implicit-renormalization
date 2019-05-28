#include "Point.hpp"
#include <iostream>

Point::Point(std::vector<double> co,double ar)
    :coor(co),a(ar)
{
}

Point::Point(std::vector<double> co)
    :coor(co),a(1)
{
}

Point::Point(int d)
{
    std::vector<double> co(d,0);
    coor=co;
    a=1;
}

Point::Point(){
    std::vector<double> co(1,0);
    coor=co;
    a=1;    
}

Point::Point(double co0,double area){
    std::vector<double> co(1,co0);
    coor=co;
    a=area;
}

Point operator*(const Point& p1,const Point& p2){
    std::vector<double> co;
    for(int i=0;i<p1.dimension();i++) co.push_back(p1[i]);
    for(int i=0;i<p2.dimension();i++) co.push_back(p2[i]);    
    return Point(co,p1.area()*p2.area());
}

void Point::print(){
    std::cout<<dimension()<<"d"<<std::endl;
    for(int i=0;i<coor.size();i++) std::cout<<coor[i]<<"\t";
    std::cout<<std::endl<<"area:"<<area()<<std::endl;
}
