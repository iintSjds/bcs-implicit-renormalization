#ifndef _GENERAL_MULTIDIMENSION_POINT_
#define _GENERAL_MULTIDIMENSION_POINT_


#include <vector>

class Point{
    // points of a grid
    // contains the location(a set of coordinates)
    // and the area around the location, useful for integration
public:
    Point();
    Point(int);
    Point(std::vector<double> co);
    Point(std::vector<double> co, double ar);
    Point(double co0,double area=0);
    
    //return coordinates of the point
    std::vector<double>& coordinates() {return coor;};
    //return dimension of the point
    int dimension() const{return coor.size();};
    //return i-th coordinate of the point
    double operator[](int i) const{return coor[i];};
    //return area of the region around the point
    double area() const{return a;};
    double set_a(double area) {return a=area;};
    double set_co(int i,double co) {return coor[i]=co;};
    void print();// for debug
private:
    std::vector<double> coor;
    double a;
};

Point operator*(const Point& p1,const Point& p2);


#endif
