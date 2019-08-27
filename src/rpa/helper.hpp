# ifndef _HELPER_FUNCTION_
# define _HELPER_FUNCTION_

# include "../function/function.hpp"
# include <iostream>
# include <ostream>
# include <fstream>
# include <iomanip>
# include <cstdlib>
# include <cmath>
# include <omp.h>
# include <string>
# include "H5Cpp.h"

class Helper_Func{
public:
    Helper_Func(Function pf,double e2_,int n);
    Helper_Func(std::string filename,double e2_,int n);
    int order() const{return H.size();}

    Function& operator[](int i){return H[i];}
    void save(std::string filename);
    Grid grid(){return H[0].grid();};
private:
    double e2;
    std::vector<Function> H;
};

# endif
