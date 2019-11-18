# include "abstract_function.hpp"
# include <stdexcept>
# include <algorithm>

// AnalyticalFunction 

int AnalyticalFunction::dimension() const {return domains.size();}
std::vector<double> AnalyticalFunction::domain(int d)const {return domains[d];}

double AnalyticalFunction::operator()(std::vector<double> arguments){
  // chech validity of arguments
  if(arguments.size()!=dimension())
    throw std::invalid_argument{"AnalyticalFunction::operator()"};
  for(int i=0;i<dimension();i++){
    if(domain(i)[0]>arguments[i] || domain(i)[1]<=arguments[i])
      throw std::domain_error{"AnalyticalFunction::operator()"};
  }
  return func( arguments);
}


AnalyticalFunction::AnalyticalFunction(
                                       double (&f)(std::vector<double> arguments),
                                       std::vector<std::vector<double> > d)
  :
  func(f),
  domains(d)
{
}

// TabFunction

TabFunction::TabFunction(std::vector<std::vector<double> > grids_)
  :grids(grids_)
{
  // check validity
  int N=1;
  for(auto i : grids) {
    if(!std::is_sorted(i.begin(),i.end(),std::less_equal<double>()))
       throw std::invalid_argument{"TabFunction:grids should be sorted"};
    N*=i.size();
  }
  if(N==0)
    throw std::invalid_argument{"TabFunction:All grids should have non-zero components"};
  // generate val with zero init value.
  val.assign(N,0);
}

TabFunction::TabFunction(std::vector<std::vector<double> > grids_, std::vector<double> val_)
  :grids(grids_),val(val_)
{
  // check validity
  int N=1;
  for(auto i : grids) {
    if(!std::is_sorted(i.begin(),i.end(),std::less_equal<double>()))
      throw std::invalid_argument{"TabFunction:grids should be sorted"};
    N*=i.size();
  }
  if(N==0)
    throw std::invalid_argument{"TabFunction:All grids should have non-zero components"};
  if(N!=val.size())
    throw std::logic_error{"TabFunction:Total Num. of points should match size of val."};
}

double TabFunction::operator() (std::vector<double> arguments){
  // first check validity
  if(arguments.size()!=dimension())
    throw std::invalid_argument{"TabFunction::operator()"};

  // find nearest point, check range at the same time
  std::vector<int> indices;
  int index=0;
  for(int i=0;i<dimension();i++){
    if(arguments[i]<domain(i)[0] || arguments[i]>=domain(i)[1])
      throw std::domain_error{"TabFunction::operator(): arguments out of boundary"};
    indices.push_back(
                      std::lower_bound(grids[i].begin(),grids[i].end(),arguments[i],std::less_equal<double>())
                      - grids[i].begin()
                      );
    indices[i]=(arguments[i]-grids[i][indices[i]-1]<grids[i][indices[i]]-arguments[i])?indices[i]-1:indices[i];
    index=index*grids[i].size()+indices[i];
  }
  return val[index];
}

double TabFunction::operator() (std::vector<double> arguments, bool isExtrapolate){
  if(!isExtrapolate) return operator()(arguments);
  else{
    //TBA
    return operator()(arguments);
  }
}

int TabFunction::dimension() const{
  return grids.size();
}

std::vector<double> TabFunction::domain(int d) const {
  std::vector<double> result;
  result.push_back(grids[d].front());
  result.push_back(grids[d].back());
  return result;
}



