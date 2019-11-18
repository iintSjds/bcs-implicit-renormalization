# include "abstract_function.hpp"
# include <iostream>

double f(std::vector<double> x){return x[0]*x[0];}

int main(){
  std::vector<double> domainX;
  domainX.push_back(0);domainX.push_back(1);
  std::vector<std::vector<double> > domains;
  domains.push_back(domainX);

  AnalyticalFunction Func(f,domains);

  std::cout<<Func.domain(0)[0]<<"<"<<Func(std::vector<double>(1,0.5))<<std::endl;

  std::vector<std::vector<double> > grids;
  for(int i=0;i<3;i++){
    std::vector<double> g(5,0);
    for(int j=0;j<g.size();j++) g[j]=j*(1+0.1*i);
    grids.push_back(g);
  }
  std::vector<double> arg(3,0.8);
  TabFunction TFunc(grids);
  std::cout<<TFunc[3]<<TFunc(arg)<<std::endl;
  for(int i=0;i<TFunc.size();i++) TFunc[i]=i;
  std::cout<<TFunc[3]<<TFunc(arg)<<std::endl;
  return 0;
}
