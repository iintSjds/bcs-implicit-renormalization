# ifndef _MULTIDIMENSION_GRID_
# define _MULTIDIMENSION_GRID_

# include <vector>

class Grid{
public:
    // constructors
    Grid(std::vector<std::vector<double> > g_in);
    Grid(std::vector<double> g1,std::vector<double> g2);
    
    // basic information
    int dimension() const {return g.size();};
    int size() const;
    int lengths(int i) const {return g[i].size();};
    // point coordinates
    std::vector<double> point(int n) const;
    std::vector<double> operator[](int n) const;
    double coordinate(int n,int d) const;
    // su information
    std::vector<double> gg(int i) const{return g[i];};
    

private:
    std::vector<std::vector<double> > g;
};


# endif
