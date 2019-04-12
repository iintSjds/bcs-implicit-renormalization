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
    int N0,N1;
    std::vector<double> delta0,delta1;
    std::vector<double> freq0,freq1;

    double Gamma(int i, int j ,int k);//all factors combined into gamma function
    // k=0,1,2,3; stand for 00,01,10,11
public:
    Iterator(double,double,double,double,double,int);
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
    //initialize freqs and deltas
    for(int i=0;i<N0;i++){
	freq0.push_back(-(N0-1)*PI*T+i*2*PI*T);
	delta0.push_back(1);
    }
    
    for(int i=0;i<N1/2;i++){
	freq1.push_back(-(Ec-(Ec-omega_c)*2.0/N1*i));
	delta1.push_back(-1);
    }
    for(int i=0;i<N1/2;i++){
	freq1.push_back(-freq1[N1/2-1-i]);
	delta1.push_back(-1);
    }
   
    if(DEBUG) std::cout<<"end init iterator"<<std::endl;
}

double Iterator::Gamma(int i,int j, int k){
    //if(DEBUG) std::cout<<"start gamma"<<std::endl;
    std::vector<double> freqi,freqj;
    double interval;
    if(k%2==0) {freqj=freq0; interval=1;}
    else {freqj=freq1; interval=(Ec-omega_c)/N1/PI/T;}
    if(k/2==0) freqi=freq0;
    else freqi=freq1;
    //if(DEBUG) std::cout<<"end gamma"<<std::endl;    
    return interval*(freqi[i]-freqj[j])*(freqi[i]-freqj[j])
	/(omega*omega+(freqi[i]-freqj[j])*(freqi[i]-freqj[j]))
	*PI*T*g/std::fabs(freqj[j]);
}

int Iterator::update1(int M){
    // update delta1 with given delta0 by mc methods
    if(DEBUG) std::cout<<"start update1"<<std::endl;    
    std::vector<double> delta_count;
    long int count=1;
    int n=0;
    int index1=0;
    double f=1;
    //define Gamma_10 * delta0 and init
    std::vector<double> gamma10delta0;
    for(int i=0;i<N1;i++){
	double d1=0;
	for(int j=0;j<N0;j++){
	    d1+=Gamma(i,j,2)*delta0[j];
	}
	gamma10delta0.push_back(d1);
	delta_count.push_back(0);
    }

    for(int i=0;i<M;i++){
	if(DEBUG) std::cout<<"\tstart mc loop"<<std::endl;	
	//update the config
	int new_n=3*rndm(),new_index1=0;
	if (new_n!=0)
	    new_index1=N1*rndm();

	//calculate acceptance ratio
	if(DEBUG) std::cout<<"\tcalculate R"<<std::endl;
	double R;
	double new_f;
	if(new_n==0) new_f=1;
	else if (new_n==2) new_f=gamma10delta0[new_index1];
	else{
	    new_f=0;
	    for(int j=0;j<N1;j++)
		new_f+=Gamma(new_index1,j,3)*delta_count[j]/(i+1)/count;
	}
	if(std::fabs(new_f)>std::fabs(f)) R=1;
	else R=std::fabs(new_f)/std::fabs(f);

	//test, R===1
	R=1;
	//update if accepted
	if(rndm()<R){
	    n=new_n;
	    index1=new_index1;
	    f=new_f;
	}
	if(DEBUG) std::cout<<"\tcount"<<std::endl;
	//count
	if(n==0) count++;
	else{
	    //delta_count[index1]-=((f>0)?1:-1);
	    delta_count[index1]-=f;	    
	}
    }
    //update result to delta1

    for(int i=0;i<N1;i++){
	double d1=0;
	for(int j=0;j<N1;j++){
	    d1-=Gamma(i,j,3)*delta_count[j]/N1/count;
	}
	d1-=gamma10delta0[i];
	delta1[i]=d1;
    }
    if(DEBUG) std::cout<<"end update1"<<std::endl;    
    return 1;
}

int Iterator::update0(double shift=3){
    std::vector<double> temp;
    if(DEBUG) std::cout<<"start update0"<<std::endl;    
    for(int i=0;i<N0;i++){
	temp.push_back(shift*delta0[i]);
	for(int j=0;j<N0;j++){
	    temp[i] += -Gamma(i,j,0)*delta0[j];
	}
	for(int j=0;j<N1;j++) temp[i] += -Gamma(i,j,1)*delta1[j];
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
	for(int j=0;j<N1;j++) temp[i] += -Gamma(i,j,1)*delta1[j];
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
    int N1=248;
    int N=30,M=N1*100;
    double g=2/PI,o=0.1,e=PI;
    Iterator it(g,o,e,t,o_c,N1);
    for(int i=0;i<N;i++){
	it.update1(M);
	for(int j=0;j<30;j++)it.update0(1);
    }
    it.print();
    return it.measure();
}

int main(){
    double t=0.00392157;
    std::cout<<test(t,0.1)<<std::endl;
    return 0;
}
