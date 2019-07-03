# include "grid.hpp"

Grid::Grid(std::vector<std::vector<double>> g_in)
    :g(g_in)
{
}

int Grid::size() const{
    int N=1;
    for(int i=0;i<g.size();i++){
	N*=g[i].size();
    }
    return N;
}

std::vector<double> Grid::point(int n) const{
    // n=nd+Ld*(n(d-1)+L(d-1)*... ))
    std::vector<double> p(g.size(),0);
    for(int i=g.size()-1;i>-1;i--){
	p[i]=g[i][n%g[i].size()];
	n=n/g[i].size();
    }
    return p;
}


std::vector<double> Grid::operator[](int n) const{
    // n=nd+Ld*(n(d-1)+L(d-1)*... ))
    std::vector<double> p(g.size(),0);
    for(int i=g.size()-1;i>-1;i--){
	p[i]=g[i][n%g[i].size()];
	n=n/g[i].size();
    }
    return p;
}

double Grid::coordinate(int n,int d) const{
    for(int i=g.size()-1;i>d;i--){
	n=n/g[i].size();
    }
    return g[d][n%g[d].size()];
}

