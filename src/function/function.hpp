# ifndef _FUNCTION_ON_GRID_
# define _FUNCTION_ON_GRID_

# include "grid.hpp"

class Function{
public:
    //constructors
    Function(Grid g_in);
    Function(std::vector<double> g1,std::vector<double> g2);
    
    //basic information
    int size() const {return g.size();};
    int dimension() const {return g.dimension();};

    //return value
    double& operator[](int n) {return val[n];};
    std::vector<double> point(int n) const{return g[n];};
    double coordinate(int n,int d) const{return g.coordinate(n,d);};
    Grid grid() {return g;}
    std::vector<double> value() const{return val;}
private:
    std::vector<double> val;
    Grid g;
};

# endif
