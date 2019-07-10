# include "../function/function.hpp"
# include <iostream>
# include <ostream>
# include <fstream>
# include <cmath>

# include <string>
# include "H5Cpp.h"

# include <trng/yarn2.hpp>
# include <trng/uniform_dist.hpp>
# include <trng/exponential_dist.hpp>

trng::yarn2 gen;


// calculate bold interaction w
// produce a bare w result file and a helper result file for further use

const int M=100000;
const double pi=3.141592653;

class Pi_Function{
public:
    Pi_Function(Grid& g_in,double T_in,double mu_in,double m_in,int N);

    int size() const {return func.size();}
    int dimension() const {return func.dimension();}

    double operator[](int n) {return func[n];}
    std::vector<double> point(int n) const{return func.point(n);}
    double coordinate(int n,int d) const{return func.coordinate(n,d);}
    Grid grid(){return func.grid();}
    std::vector<double> f(){return func.value();}
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
    double I0=0,I1=0,D0=0.001;
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

    D0=I1/I0*D0+0.0000001;

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

Pi_Function test_pi(){
    double T=0.5,mu=1.0,m=0.5;
    std::vector<double> freq(10,0);
    for(int i=0;i<freq.size();i++) freq[i]=2*pi*T*(i);
    std::vector<double> mmt(50,0);
    for(int i=0;i<mmt.size();i++) mmt[i]=(i<13)?(0.0001*std::pow(2,i)):(0.2*i-2.0);
    std::vector<std::vector<double>> g_in;
    g_in.push_back(freq);g_in.push_back(mmt);
    Grid g(g_in);
    Pi_Function pf(g,T,mu,m,M);

    // std::ofstream out;
    // out.open("pi.txt");
    // for(int i=0;i<pf.size();i++){
    // 	out<<std::setw(10)<<pf.coordinate(i,0)<<"\t"
    // 	   <<std::setw(10)<<pf.coordinate(i,1)<<"\t"
    // 	   <<pf[i]<<std::endl;//<<"\t"
    // 	    //<<1/(pf.coordinate(i,1)*pf.coordinate(i,1)+pf[i])<<std::endl;//"\t"
    // 	    //<<1/(1+pf[i]/pf.coordinate(i,1)/pf.coordinate(i,1))<<std::endl;
    // }
    return pf;
}

int main(){
    double e2=1.0;
    Pi_Function pf=test_pi();

    H5::H5File file("pi.h5",H5F_ACC_TRUNC);
    H5::Group g1(file.createGroup("/pi"));
    hsize_t dims[1];
    //store freq grid
    dims[0]=pf.grid().lengths(0);
    H5::DataSpace dataspace=H5::DataSpace(1,dims);
    H5::DataSet dataset(file.createDataSet("/pi/w",
					   H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(pf.grid().gg(0)[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    //store momentum grid
    dims[0]=pf.grid().lengths(1);
    dataspace=H5::DataSpace(1,dims);
    dataset=H5::DataSet(file.createDataSet("/pi/q",
					   H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(pf.grid().gg(1)[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    //store val grid
    dims[0]=pf.f().size();
    dataspace=H5::DataSpace(1,dims);
    dataset=H5::DataSet(file.createDataSet("/pi/pi",
					   H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(pf.f()[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    // Function H1(pf.grid());
    // for(int i=0;i<H1.size();i++){
    // 	double val=0;
    // 	for(int j=(i/H1.grid().lengths(1))*H1.grid().lengths(1);j<i;j++){
    // 	    val+=(1/(pf.coordinate(j,1)/4/pi/e2+pf[j]/pf.coordinate(j,1))
    // 		+1/(pf.coordinate(j+1,1)/4/pi/e2+pf[j+1]/pf.coordinate(j+1,1)))
    // 		//e2*(pf.coordinate(j,1)+pf.coordinate(j+1,1))
    // 		*(pf.coordinate(j+1,1)-pf.coordinate(j,1))/2;
    // 	}
    // 	H1[i]=val;
    // }
    // std::ofstream hout;
    // hout.open("h1.txt");
    // for(int i=0;i<pf.size();i++){
    // 	hout<<std::setw(10)<<H1.coordinate(i,0)<<"\t"
    // 	   <<std::setw(10)<<H1.coordinate(i,1)<<"\t"
    // 	   <<H1[i]<<std::endl;//<<"\t"
    // 	    //<<1/(pf.coordinate(i,1)*pf.coordinate(i,1)+pf[i])<<std::endl;//"\t"
    // 	    //<<1/(1+pf[i]/pf.coordinate(i,1)/pf.coordinate(i,1))<<std::endl;
    // }

    // std::vector<std::vector<double>> gwin;
    // gwin.push_back(pf.grid().gg(0));
    // std::vector<double> mmt(20,0);
    // for(int i=0;i<mmt.size();i++) mmt[i]=0.1+0.2*i;
    // gwin.push_back(mmt);gwin.push_back(mmt);

    // Function w0(gwin);
    // for(int i=0;i<w0.size();i++){
    // 	double h1=0,h2=0;
    // 	double w=w0.coordinate(i,0);
    // 	double p1=w0.coordinate(i,1)+w0.coordinate(i,2);
    // 	double p2=std::abs(w0.coordinate(i,1)-w0.coordinate(i,2));
    // 	int m=0,n1=0,n2=0;
    // 	for(int j=0;j<H1.size();j+=H1.grid().gg(1).size())
    // 	    if(std::abs(H1.coordinate(j,0)-w)<1e-9)
    // 		m=j;
    // 	for(int j=0;j<H1.grid().gg(1).size();j++){
    // 	    if(H1.coordinate(m+j,1)<p1) n1++;
    // 	    if(H1.coordinate(m+j,1)<p2) n2++;	    
    // 	}
    // 	if(n1>0)
    // 	    h1=( (H1[m+n1]-H1[m+n1-1]) *p1
    // 		 +H1[m+n1-1]*H1.coordinate(m+n1,1)-H1[m+n1]*H1.coordinate(m+n1-1,1) )
    // 		/(H1.coordinate(m+n1,1)-H1.coordinate(m+n1-1,1));
    // 	if(n2>0)
    // 	    h2=( (H1[m+n2]-H1[m+n2-1]) *p2
    // 		 +H1[m+n2-1]*H1.coordinate(m+n2,1)-H1[m+n2]*H1.coordinate(m+n2-1,1) )
    // 		/(H1.coordinate(m+n2,1)-H1.coordinate(m+n2-1,1));
    // 	w0[i]=(h1-h2)/w0.coordinate(i,1)/w0.coordinate(i,2);
    // 	// std::cout<<std::setw(8)<<w<<"\t"
    // 	// 	 <<std::setw(5)<<n1<<"\t"
    // 	// 	 <<std::setw(8)<<h1<<"\t"
    // 	// 	 <<std::setw(5)<<n2<<"\t"
    // 	// 	 <<std::setw(8)<<h2<<"\t"
    // 	// 	 <<w0[i]<<std::endl;//<<"\t"
       
    // }
    // std::ofstream wout;
    // wout.open("w0.txt");
    // for(int i=0;i<w0.size();i++){
    // 	wout<<std::setw(10)<<w0.coordinate(i,0)<<"\t"
    // 	   <<std::setw(10)<<w0.coordinate(i,1)<<"\t"
    // 	   <<std::setw(10)<<w0.coordinate(i,2)<<"\t"	    
    // 	   <<w0[i]<<std::endl;//<<"\t"
    // 	    //<<1/(pf.coordinate(i,1)*pf.coordinate(i,1)+pf[i])<<std::endl;//"\t"
    // 	    //<<1/(1+pf[i]/pf.coordinate(i,1)/pf.coordinate(i,1))<<std::endl;
    // }
    
    
    return 0;
}
