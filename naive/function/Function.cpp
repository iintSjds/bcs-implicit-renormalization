#include "Function.hpp"

TabFunc::TabFunc(Grid& g,double(&f)(Point p))
    :grid(g)
{
    for(int i=0;i<g.size();i++)
	funcVal.push_back(f(g[i]));
}
