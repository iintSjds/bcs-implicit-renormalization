#ifndef _DOUBLE_GREENS_FUNCTION_
#define _DOUBLE_GREENS_FUNCTION_

#include <vector>

//  generate GG function

class GG{
public:
    double omega_c;// cut off freq
    double Ec;// upper bound of freq
    double T;// temperature
    double mu;// chemical potential
    double t;// kinetic 
    
    int dim=2;//dimension of momentum space
    int L;//No. momentum points

    std::vector<double> freq0,freq1;
    std::vector<double> gg0,gg1;
    double gg(double freq,std::vector<double> k);
    int N0,N1;//No. of frequency point in/out of cut off
    double interval1;
    GG(double mu_p,double t_p,double Ec_p,double T_p,double omega_c_p,int N1_p,int L_p);
    double result0(int i);
    double result1(int i);
};

#endif
