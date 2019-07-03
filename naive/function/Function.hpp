#ifndef _GENERAL_FUNCTION_ON_GRID_
#define _GENERAL_FUNCTION_ON_GRID_

#include <vector>
#include "Grid.hpp"
#include "Point.hpp"

class Function{
    // abstract class of function
public:
    virtual double operator[](int i)const=0;  // return func value of site i
    virtual Point& at(int i)=0;  // return Point of site i

    virtual int size()const=0;
    virtual int dimension()const=0;
    virtual std::vector<int> lengths()const=0;
    
};

class TabFunc:public Function{
public:
    TabFunc(Grid& g, double(& f)(Point p));

    double operator[](int i)const{return funcVal[i];};
    Point& at(int i){return grid[i];};
    
    int size()const{return funcVal.size();};
    int dimension()const{return grid.dimension();};
    std::vector<int> lengths()const {return grid.lengths();};
private:
    std::vector<double> funcVal;
    Grid& grid;
};

class BareFunc:public Function{
public:
    BareFunc(Grid& g, double(&f)(Point p)):grid(g),func(f){};

    double operator[](int i)const{return func(grid[i]);}
    Point& at(int i){return grid[i];}

    int size()const{return grid.size();}
    int dimension()const{return grid.dimension();}
    std::vector<int> lengths()const{return grid.lengths();}
private:
    Grid& grid;
    double(&func)(Point p);
};

#endif
