# ifndef _GENERAL_ABSTRACT_FUNCTION_
# define _GENERAL_ABSTRACT_FUNCTION_

# include <vector>

class Function{
  // an abstract class for functions
  // realization could be tabulated function, analytical function, etc.
  // require all functions to provide such interface:
public:
  // receive arguments from a vector of double, return result in double
  virtual double operator() (std::vector<double> arguments) =0;
  // return dimension of the arguments
  virtual int dimension() const =0;
  // return range of specific dimension
  virtual std::vector<double> domain(int) const =0;

  virtual ~Function(){}; // destructor
};

class AnalyticalFunction :public Function
{
private:
  std::vector<std::vector<double> > domains;
  double (&func)(std::vector<double> arguments);
public:
  double operator() (std::vector<double> arguments) override;
  int dimension() const override;
  std::vector<double> domain(int d) const override;

  AnalyticalFunction(double (&f)(std::vector<double> arguments), std::vector<std::vector<double> > domains);
  ~AnalyticalFunction(){};
};

class TabFunction :public Function
{
private:
  std::vector<double> val;
  std::vector<std::vector<double> > grids;
public:
  double & operator[](int i) {return val[i];}
  unsigned size() {return val.size();}

  double operator() (std::vector<double> arguments) override;
  double operator() (std::vector<double> arguments, bool isExtrapolate);

  int dimension() const override;
  std::vector<double> domain(int d) const override;

  TabFunction(std::vector<std::vector<double> > grids_);
  TabFunction(std::vector<std::vector<double> > grids_, std::vector<double> val_);
  ~TabFunction(){};
};

#endif
