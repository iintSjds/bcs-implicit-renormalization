#include "GG.hpp"
#include <cmath>
#define PI 3.1415926535384626

double GG::gg(double freq,std::vector<double> k){
    double e=4*t;
    for(int i=0;i<dim;i++){
	e-=2*t*std::cos(k[i]);
    }
    return 1/(freq*freq+(e-mu)*(e-mu));
}

GG::GG(double mu_p,double t_p,double Ec_p,double T_p,double omega_c_p,int N1_p,int L_p)
    :mu(mu_p),t(t_p),Ec(Ec_p),T(T_p),omega_c(omega_c_p),N1(N1_p),L(L_p)
{
    //if(DEBUG) std::cout<<"start init iterator"<<std::endl;
    N0=std::floor(omega_c/PI/T/2)*2; // even number, accurate below omega_c
    //N1=std::floor(Ec/PI/T/2)*2-N0; // all the rest
    //initialize freqs
    for(int i=0;i<N0;i++){
	freq0.push_back(-(N0-1)*PI*T+i*2*PI*T);
    }
    
    for(int i=0;i<N1/2;i++){
	//freq1.push_back(i*2*PI*T-(N0-1)*PI*T-N1*PI*T);
	freq1.push_back(-(Ec-(Ec-omega_c)*2.0/N1*(i+0.5)));	
    }
    for(int i=0;i<N1/2;i++){
	freq1.push_back(-freq1[N1/2-1-i]);
    }
    N1=freq1.size();
    interval1=(Ec-omega_c)/N1/PI/T;   
    //if(DEBUG) std::cout<<"end init iterator"<<std::endl;
    for(int i=0;i<N0;i++){
	double temp=0;
	for(int px=0;px<L;px++){
	    for(int py=0;py<L;py++){
		std::vector<double> k;
		k.push_back(-PI+px*2*PI/(L));
		k.push_back(-PI+py*2*PI/(L));
		temp+=gg(freq0[i],k)/(L)/(L);
	    }
	}
	gg0.push_back(temp);
    }
    for(int i=0;i<N1;i++){
	double temp=0;
	for(int px=0;px<L;px++){
	    for(int py=0;py<L;py++){
		std::vector<double> k;
		k.push_back(-PI+px*2*PI/(L));
		k.push_back(-PI+py*2*PI/(L));
		temp+=gg(freq1[i],k)/(L)/(L);
	    }
	}
	gg1.push_back(temp);
    }

}

double GG::result0(int i){
    return gg0[i];
}
double GG::result1(int i){
    return gg1[i];
}
