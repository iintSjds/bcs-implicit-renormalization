# ifndef _MULTIDIMENSION_GRID_
# define _MULTIDIMENSION_GRID_

# include <vector>

class Grid{
public:
    // constructors
    Grid(std::vector<std::vector<double>> g_in);

    // basic information
    int dimension() const {return g.size();};
    int size() const;

    // point coordinates
    std::vector<double> point(int n) const;
    std::vector<double> operator[](int n) const;
    double coordinate(int n,int d) const;
    

private:
    std::vector<std::vector<double>> g;
};


# endif
