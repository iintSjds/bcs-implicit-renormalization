#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

//constants

#define PI 3.1415926535897932
#define ALPHA 0.8819384944323748

double compress_func(double frac){
    //return 16*std::ceil(frac*20)-15;
    return 1;
}



class Iterator{
private:
    double g;
    double omega;
    double Ec;
    double T;
    int L;
    std::vector<long double> delta;
    std::vector<long double> temp;
    std::vector<long double> freq;
    std::vector<std::vector<long double>> gamma;

    double Gamma(int i, int j);
    double GG(int i);
public:
    Iterator(double,double,double,double,int);
    double update(double shift);
    double measure();
    int print();
};

Iterator::Iterator(double G, double O, double E, double t, int l)
    :g(G),omega(O),Ec(E),T(t),L(l)
{
    int i=0;
    double f=-E;
    while(f<0){
	freq.push_back(f);
	f+=(E-PI*T)/(L-1)*compress_func(std::fabs(f)/omega);
	//if(i<L) freq.push_back(-E+(E-PI*T)/(L-1)*i);
	//else freq.push_back(E-(E-PI*T)/(L-1)*(2*L-i-1));
    }
    for(int i=freq.size()-1;i>-1;i--){
	freq.push_back(-freq[i]);
    }
    L=freq.size()/2;
    for(int i=0;i<2*L;i++){	
	//delta.push_back(1-1.5*std::abs(i-L)*1.0/L);
	delta.push_back((std::fabs(freq[i])<omega)?1:-1);
	temp.push_back(0);
    }
    for(int i=0;i<2*L;i++){
	std::vector<long double> g;
	for(int j=0;j<2*L;j++) g.push_back(Gamma(i,j));

	gamma.push_back(g);
    }
}

double Iterator::Gamma(int i,int j){
    return (freq[i]-freq[j])*(freq[i]-freq[j])/(omega*omega+(freq[i]-freq[j])*(freq[i]-freq[j]));
}

double Iterator::GG(int i){
    return 1/std::fabs(freq[i])*compress_func(std::fabs(freq[i])/omega);
}

double Iterator::update(double shift=0){
    for(int i=0;i<2*L;i++){
	temp[i]=shift*delta[i];
	for(int j=0;j<2*L;j++){
	    temp[i] += -PI*T*g*delta[j]*GG(j)*gamma[i][j];
	}
    }
    //double max=*std::max_element(temp.begin(),temp.end());
    //if(max<0) max=1;
    double max=0;
    //for (int i=0;i<2*L;i++) max+=temp[i]*temp[i];
    //max=std::sqrt(max/2/L);
    max=temp[L];
    //std::cout<<max<<std::endl;
    for(int i=0;i<2*L;i++){
	delta[i]=temp[i]/max;//renormalize
    }
    return max;
}

double Iterator::measure(){
    for(int i=0;i<2*L;i++){
	temp[i]=0;
	for(int j=0;j<2*L;j++){
	    temp[i] += -PI*T*g*delta[j]*GG(j)*Gamma(i,j);
	}
    }
    return std::accumulate(temp.begin(),temp.end(),0.0)/std::accumulate(delta.begin(),delta.end(),0.0);
}

int Iterator::print(){
   for(int i=0;i<2*L;i++){
       std::cout<<freq[i]<<"\t"<<delta[i]<<std::endl;
   }
}


double test(double l){
    int N=31;
    double g=2.0/PI,o=0.5,e=15,t=e/(l-1)/PI;
    Iterator it(g,o,e,t,l/2);
    for(int i=0;i<N;i++) std::cout<<it.update(2)<<std::endl;
    it.print();
    return it.measure();
}

//  _ __ ___   __ _(_)_ __  
// | '_ ` _ \ / _` | | '_ \ 
// | | | | | | (_| | | | | |
// |_| |_| |_|\__,_|_|_| |_|


int main(){
    // for(int l=16;l<8000;l*=2){
    // 	std::cout<<1.0/(l-1)<<"\t"<<test(l)<<std::endl;
    // 	}
    std::cout<<test(4774)<<std::endl;
    return 0;

}

//              _
//  ___ _ _  __| |
// / -_) ' \/ _` |
// \___|_||_\__,_|
