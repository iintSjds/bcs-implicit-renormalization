#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<numeric>

#define PI 3.1415926535897932

#include <random>

#define DEBUG 0

//random generator setup
std::mt19937_64 generator;

std::uniform_real_distribution<double> dist(0.0,1.0);

double rndm(){
    return dist(generator);
}

class Iterator{
private:
    double g;
    double omega;
    double Ec;
    double T;
    double omega_c;
    double interval1;
    int N0,N1;
    std::vector<double> delta0,delta1;
    std::vector<double> freq0,freq1;

    double Gamma(int i, int j ,int k);//all factors combined into gamma function
    // k=0,1,2,3; stand for 00,01,10,11
public:
    Iterator(double,double,double,double,double, int);
    int update0(double shift);
    int update1(int M);
    double measure();
    int print();
};

Iterator::Iterator(double G, double O, double E, double t,double O_c, int n1)
    :g(G),omega(O),Ec(E),T(t),omega_c(O_c),N1(n1)
{
    if(DEBUG) std::cout<<"start init iterator"<<std::endl;
    N0=std::floor(omega_c/PI/T/2)*2; // even number, accurate below omega_c
    //N1=std::floor(Ec/PI/T/2)*2-N0; // all the rest
    //initialize freqs and deltas
    for(int i=0;i<N0;i++){
	freq0.push_back(-(N0-1)*PI*T+i*2*PI*T);
	delta0.push_back(1);
    }
    
    for(int i=0;i<N1/2;i++){
	//freq1.push_back(i*2*PI*T-(N0-1)*PI*T-N1*PI*T);
	freq1.push_back(-(Ec-(Ec-omega_c)*2.0/N1*(i+0.5)));	
	delta1.push_back(-1);
    }
    for(int i=0;i<N1/2;i++){
	freq1.push_back(-freq1[N1/2-1-i]);
	delta1.push_back(-1);
    }
    interval1=(Ec-omega_c)/N1/PI/T;   
    if(DEBUG) std::cout<<"end init iterator"<<std::endl;
}

double Iterator::Gamma(int i,int j, int k){
    //if(DEBUG) std::cout<<"start gamma"<<std::endl;
    std::vector<double> freqi,freqj;
    if(k%2==0) {freqj=freq0;}
    else {freqj=freq1;}
    if(k/2==0) freqi=freq0;
    else freqi=freq1;
    //if(DEBUG) std::cout<<"end gamma"<<std::endl;    
    return (freqi[i]-freqj[j])*(freqi[i]-freqj[j])
	/(omega*omega+(freqi[i]-freqj[j])*(freqi[i]-freqj[j]))
	*PI*T*g/std::fabs(freqj[j]);
}

int Iterator::update1(int M){
    // update delta1 with given delta0 by iterate methods
    if(DEBUG) std::cout<<"start update1"<<std::endl;    
    std::vector<double> sum_delta;
    //define Gamma_10 * delta0 and init
    std::vector<double> gamma10delta0;
    for(int i=0;i<N1;i++){
	double d1=0;
	for(int j=0;j<N0;j++){
	    d1+=Gamma(i,j,2)*delta0[j];
	}
	gamma10delta0.push_back(d1);
	sum_delta.push_back(delta1[i]);
    }

    for(int i=1;i<M+1;i++){
	for(int j=0;j<N1;j++){
	    for(int k=0;k<N1;k++)
		sum_delta[j]+= -interval1*Gamma(j,k,3)*sum_delta[k]/i;
	    sum_delta[j]+= -gamma10delta0[j];
	}
    }
    //update result to delta1

    for(int i=0;i<N1;i++){
	double d1=0;
	for(int j=0;j<N1;j++){
	    d1-=interval1*Gamma(i,j,3)*sum_delta[j]/M;
	}
	d1-=gamma10delta0[i];
	delta1[i]=d1;
    }
    if(DEBUG) std::cout<<"end update1"<<std::endl;    
    return 1;
}

int Iterator::update0(double shift=1){
    std::vector<double> temp;
    if(DEBUG) std::cout<<"start update0"<<std::endl;    
    for(int i=0;i<N0;i++){
	temp.push_back(shift*delta0[i]);
	for(int j=0;j<N0;j++){
	    temp[i] += -Gamma(i,j,0)*delta0[j];
	}
	for(int j=0;j<N1;j++) temp[i] += -interval1*Gamma(i,j,1)*delta1[j];
    }
    //double max=*std::max_element(temp.begin(),temp.end());
    //if(max<0) max=1;
    double max=temp[N0/2];
    //for (int i=0;i<2*L;i++) max+=temp[i]*temp[i];
    //max=std::sqrt(max/2/L);
    //std::cout<<max<<std::endl;
    for(int i=0;i<N0;i++){
	delta0[i]=temp[i]/max;//renormalize
    }
    if(DEBUG) std::cout<<"end update0"<<std::endl;    
    return 1;    
}

double Iterator::measure(){
    std::vector<double> temp;
    for(int i=0;i<N0;i++){
	temp.push_back(0);
	for(int j=0;j<N0;j++){
	    temp[i] += -Gamma(i,j,0)*delta0[j];
	}
	for(int j=0;j<N1;j++) temp[i] += -interval1*Gamma(i,j,1)*delta1[j];
    }

    return std::accumulate(temp.begin(),temp.end(),0.0)/std::accumulate(delta0.begin(),delta0.end(),0.0);
}

int Iterator::print(){
    for(int i=0;i<N1/2;i++)
	std::cout<<freq1[i]<<"\t"<<delta1[i]<<std::endl;
    for(int i=0;i<N0;i++)
	std::cout<<freq0[i]<<"\t"<<delta0[i]<<std::endl;
    for(int i=N1/2;i<N1;i++)
	std::cout<<freq1[i]<<"\t"<<delta1[i]<<std::endl;

    return 1;
}

double test(double t,double o_c){
    int N1=400;
    int N=20,M=50;
    double g=0.7/PI,o=0.1,e=PI;
    Iterator it(g,o,e,t,o_c,N1);
    for(int i=0;i<N;i++){
	it.update1(M/2*((i*2)/N+1));
	for(int j=0;j<2;j++)
	    it.update0(1);
	//std::cout<<"\t"<<it.measure()<<std::endl;
    }
    //it.print();
    return it.measure();
}

int main(){
    double t=1.0/1023.0;
    for(int l=2000;l<40000/8;l*=2){
	t=1.0/(l);
	std::cout<<t<<"\t"<<test(t,0.1)<<std::endl;
    }
    return 0;
}
