#ifndef _GENERAL_MULTIDIMENSION_GRID_
#define _GENERAL_MULTIDIMENSION_GRID_

// a general purpose multidimensional grid
// provide conversion between index and coordinates
// use std vector to store information
// only involves grids for continuous variables, lattice not included
// only an interface defined here, implimentations should be derived

#include <vector>
#include "Point.hpp"

class Grid{
    // this should be an abstract class grid
    // for both 1d and multidimensional grid
    // which should be derived from this class
    // WHEN USED, ONE GRID SHOULD ONLY BE GENERATED ONCE
public:
    //grid information
    virtual int dimension()const=0;
    //basics of a container
    virtual int size()const=0;
    virtual Point& operator[](int i)=0;
    virtual std::vector<int> lengths()const=0;
};

class Grid_1d:public Grid{
public:
    //constructor for all purposes
    Grid_1d(double init, double fin, int size);  //evenly from init to fin, size points
    Grid_1d(std::vector<double> p);
    inline Point& operator[](int i) {return points[i];};
    inline int dimension() const{return dim;};
    inline int size() const{return points.size();};
    inline std::vector<int> lengths() const {return lens;};
private:
    int dim=1;
    std::vector<Point> points;
    std::vector<int> lens;
};

class MeshGrid:public Grid{
public:
    MeshGrid(Grid& g1,Grid& g2); // made from 2 grids

    inline Point& operator[](int i) {return points[i];};
    inline int dimension() const{return dim;};
    inline int size() const{return points.size();};
    inline std::vector<int> lengths()const{return lens;};
    
private:
    int dim;
    std::vector<int> lens;
    std::vector<Point> points;
};

#endif
