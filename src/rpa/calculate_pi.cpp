# include "../function/function.hpp"
# include <iostream>
# include <cmath>

# include <trng/yarn2.hpp>
# include <trng/uniform_dist.hpp>
# include <trng/exponential_dist.hpp>

trng::yarn2 gen;


// calculate bold interaction w
// produce a bare w result file and a helper result file for further use

const int M=1000000;
const double pi=3.141592653;

class Pi_Function{
public:
    Pi_Function(Grid& g_in,double T_in,double mu_in,double m_in,int N);

    int size() const {return func.size();}
    int dimension() const {return func.dimension();}

    double operator[](int n) {return func[n];}
    std::vector<double> point(int n) const{return func.point(n);}
    double coordinate(int n,int d) const{return func.coordinate(n,d);}
private:
    Function func;
    double T;
    double mu;
    double m;

    double f(double p,double q,double w);
    double calculate(int n,int N); // calculate func value of point n with mc
};

Pi_Function::Pi_Function(Grid& g_in,double T_in,double mu_in,double m_in,int N)
    :func(g_in),T(T_in),mu(mu_in),m(m_in)
{
    for(int i=0;i<func.size();i++){
	func[i]=calculate(i,N);
    }
}

double Pi_Function::f(double p,double q,double w){
    return p*m/2/pi/pi/q/(std::exp((p*p/2/m-mu)/T)+1)*std::log((w*w+(q*q+2*p*q)*(q*q+2*p*q)/4/m/m)/(w*w+(q*q-2*p*q)*(q*q-2*p*q)/4/m/m));
}

double Pi_Function::calculate(int n,int M){
    double w=coordinate(n,0),q=coordinate(n,1);
    double kf=std::sqrt(2*mu*m);
    trng::uniform_dist<> uniform(1.0,0);
    trng::exponential_dist<> exponential(kf);
    double I0=0,I1=0,D0=0.1;
    double P0=0.5;
    double p=uniform(gen),newp=0.1;
    int index=1,newindex=0;
    double R=1.0;

    for(int i=0;i<M/10;i++){
	newindex=(uniform(gen)<P0)?0:1;
	if (newindex==1){
	    newp=exponential(gen);
	    R=(index==1)?
		(f(newp,q,w)/f(p,q,w)/std::exp((p-newp)/kf))
		:(f(newp,q,w)*std::exp(newp/kf)/D0);
	}
	else{
	    R=(index==1)?
		(D0/f(p,q,w)/exp(p/kf))
		:1;
	}
	if(uniform(gen)<R){
	    index=newindex;
	    p=newp;
	}
	if(index==0) I0++;
	else I1++;
    }

    D0=I1/I0*D0+0.001;

    for(int i=0;i<M;i++){
	newindex=(uniform(gen)<P0)?0:1;
	if (newindex==1){
	    newp=exponential(gen);
	    R=(index==1)?
		(f(newp,q,w)/f(p,q,w)/std::exp((p-newp)/kf))
		:(f(newp,q,w)*std::exp(newp/kf)/D0);
	}
	else{
	    R=(index==1)?
		(D0/f(p,q,w)/exp(p/kf))
		:1;
	}
	if(uniform(gen)<R){
	    index=newindex;
	    p=newp;
	}
	if(index==0) I0++;
	else I1++;
    }

    return I1/I0*D0;
    
}

void test(){
    double T=0.01,mu=1.0,m=0.5;
    std::vector<double> freq(1,0);
    for(int i=0;i<freq.size();i++) freq[i]=2*pi*T*(i+1);
    std::vector<double> mmt(50,0);
    for(int i=0;i<mmt.size();i++) mmt[i]=(i<15)?(0.00001*std::pow(2,i)):(0.1*i-1.3);
    std::vector<std::vector<double>> g_in;
    g_in.push_back(freq);g_in.push_back(mmt);
    Grid g(g_in);
    Pi_Function pf(g,T,mu,m,M);

    for(int i=0;i<pf.size();i++){
	std::cout<<std::setw(10)<<pf.coordinate(i,0)<<"\t"
		 <<std::setw(10)<<pf.coordinate(i,1)<<"\t"
		 <<pf[i]<<"\t"
		 <<1/(pf.coordinate(i,1)*pf.coordinate(i,1)+pf[i])<<std::endl;//"\t"
	    //<<1/(1+pf[i]/pf.coordinate(i,1)/pf.coordinate(i,1))<<std::endl;
    }
}

int main(){
    test();
    return 0;
}
