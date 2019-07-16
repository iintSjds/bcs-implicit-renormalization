# include "function.hpp"

Function::Function(Grid g_in)
    :g(g_in),val(g_in.size(),0)
{
}

Function::Function(std::vector<double> g1,std::vector<double> g2)
    :g(g1,g2),val(g1.size()*g2.size(),0)
{
}
