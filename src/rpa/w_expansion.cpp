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
# include "./helper.hpp"

const double pi=3.141592653;

class W_Func{
public:
    W_Func(Grid gwin,Helper_Func H);
    int order() const{return W.size();}
    Function& operator[](int i){return W[i];}
private:
    std::vector<Function> W;
};

W_Func::W_Func(Grid gwin,Helper_Func H){
    for(int m=0;m<H.order();m++){
	Function w0(gwin);
#pragma omp parallel num_threads(omp_get_max_threads()-2)
	{
#pragma omp for
	    for(int i=0;i<w0.size();i++){
		double h1=0,h2=0;
		double w=w0.coordinate(i,0);
		double p1=w0.coordinate(i,1)+w0.coordinate(i,2);
		double p2=std::abs(w0.coordinate(i,1)-w0.coordinate(i,2));
		double k=w0.coordinate(i,1),p=w0.coordinate(i,2);
		int m=0,n1=0,n2=0;
		for(int j=0;j<H[0].size();j+=H[0].grid().gg(1).size())
		    if(std::abs(H[0].coordinate(j,0)-w)<1e-9)
			m=j;
		for(int j=0;j<H[0].grid().gg(1).size();j++){
		    if(H[0].coordinate(m+j,1)<p1) n1++;
		    if(H[0].coordinate(m+j,1)<p2) n2++;	    
		}
		double hp=0,hm=0;
		if(n1>0){
		    if(m==0){
			hp=H[0][m+n1];
			hm=H[0][m+n1-1];
		    }
		    if(m==1){
			hp=((k*k+p*p)*H[0][m+n1]-H[1][m+n1])/2/k/p;
			hm=((k*k+p*p)*H[0][m+n1-1]-H[1][m+n1-1])/2/k/p;
		    }
		    if(m==2){
			hp=(((k*k+p*p)*(k*k+p*p)-4*k*k*p*p)*H[0][m+n1]
			    -6*(k*k+p*p)*H[1][m+n1]
			    +3*H[2][m+n1]
			    )
			    /8/k/k/p/p;
			hm=(
			    ((k*k+p*p)*(k*k+p*p)-4*k*k*p*p)*H[0][m+n1-1]
			    -6*(k*k+p*p)*H[1][m+n1-1]
			    +3*H[2][m+n1-1]
			    )
			    /8/k/k/p/p;
		    }
		    h1=( (hp-hm) *p1
			 +hm*H[0].coordinate(m+n1,1)-hp*H[0].coordinate(m+n1-1,1) )
			/(H[0].coordinate(m+n1,1)-H[0].coordinate(m+n1-1,1));
		}
		if(n2>0){
		    if(m==0){
			hp=H[0][m+n2];
			hm=H[0][m+n2-1];
		    }
		    if(m==1){
			hp=((k*k+p*p)*H[0][m+n2]-H[1][m+n2])/2/k/p;
			hm=((k*k+p*p)*H[0][m+n2-1]-H[1][m+n2-1])/2/k/p;
		    }
		    if(m==2){
			hp=(((k*k+p*p)*(k*k+p*p)-4*k*k*p*p)*H[0][m+n2]
			    -6*(k*k+p*p)*H[1][m+n2]
			    +3*H[2][m+n2]
			    )
			    /8/k/k/p/p;
			hm=(
			    ((k*k+p*p)*(k*k+p*p)-4*k*k*p*p)*H[0][m+n2-1]
			    -6*(k*k+p*p)*H[1][m+n2-1]
			    +3*H[2][m+n2-1]
			    )
			    /8/k/k/p/p;
		    }
		    h2=( (hp-hm) *p2
			 +hm*H[0].coordinate(m+n2,1)-hp*H[0].coordinate(m+n2-1,1) )
			/(H[0].coordinate(m+n2,1)-H[0].coordinate(m+n2-1,1));
		}
		w0[i]=(h1-h2)/k/p;       
	    }
	}
	W.push_back(w0);
    }
}

int main(){
    std::string line;
    std::ifstream infile("w.in");
    double T=0.5,mu=1.0,m=0.5,rs=1.0;
    std::vector<double> mmt;
    if (infile.is_open()){
	while(std::getline(infile,line)){
	    std::size_t eq=line.find("=");
	    std::string cmd,val;
	    if(eq!=std::string::npos){
		cmd=line.substr(0,eq);
		val=line.substr(eq+1);
	    }
	    if(cmd.compare("set T")==0) T=std::strtod(val.c_str(),NULL);
	    if(cmd.compare("set mu")==0) mu=std::strtod(val.c_str(),NULL);
	    if(cmd.compare("set m")==0) m=std::strtod(val.c_str(),NULL);
	    if(cmd.compare("set rs")==0) rs=std::strtod(val.c_str(),NULL);
	    if(cmd.compare("set mmt")==0){
		while(std::getline(infile,line)&&line.compare("end grid")!=0){
		    mmt.push_back(std::strtod(line.c_str(),NULL));
		    eq=line.find("=");
		    if(eq!=std::string::npos){
			cmd=line.substr(0,eq);
			val=line.substr(eq+1);
		    }
		}
	    }
	}
	infile.close();
    }
    double e2=rs*1.0421235224;
    
    H5::H5File file;
    H5::DataSet dataset;
    file.openFile("pi.h5",H5F_ACC_RDWR);
    // read freq grid
    dataset=file.openDataSet("/pi/w");
    hsize_t dims_out[1];
    dataset.getSpace().getSimpleExtentDims(dims_out,NULL);

    std::vector<double> freq(dims_out[0],0);
    dataset.read(&(freq[0]),H5::PredType::IEEE_F64LE);

    //read momentum grid
    dataset=file.openDataSet("/pi/q");
    dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
    std::vector<double> q(dims_out[0],0);
    dataset.read(&(q[0]),H5::PredType::IEEE_F64LE);

    //read value of pi
    dataset=file.openDataSet("/pi/pi");
    dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
    std::vector<double> pval(dims_out[0],0);
    dataset.read(&(pval[0]),H5::PredType::IEEE_F64LE);

    //restore pf
    std::vector<std::vector<double> > g_in;
    g_in.push_back(freq);g_in.push_back(q);
    Function pf(g_in);
    for(int i=0;i<pf.size();i++) pf[i]=pval[i];

    Function H1(pf.grid());
#pragma omp parallel num_threads(omp_get_max_threads()-2)
    {
#pragma omp for
    for(int i=0;i<H1.size();i++){
    	double val=0;
    	for(int j=(i/H1.grid().lengths(1))*H1.grid().lengths(1);j<i;j++){
    	    val+=(1/(pf.coordinate(j,1)/4/pi/e2+pf[j]/pf.coordinate(j,1))
    		+1/(pf.coordinate(j+1,1)/4/pi/e2+pf[j+1]/pf.coordinate(j+1,1)))
    		//e2*(pf.coordinate(j,1)+pf.coordinate(j+1,1))
    		*(pf.coordinate(j+1,1)-pf.coordinate(j,1))/2;
    	}
    	H1[i]=val;
    }
    }
    std::ofstream hout;
    hout.open("h1.txt");
    for(int i=0;i<pf.size();i++){
    	hout<<std::setw(10)<<H1.coordinate(i,0)<<"\t"
    	   <<std::setw(10)<<H1.coordinate(i,1)<<"\t"
    	   <<H1[i]<<std::endl;//<<"\t"
    	    //<<1/(pf.coordinate(i,1)*pf.coordinate(i,1)+pf[i])<<std::endl;//"\t"
    	    //<<1/(1+pf[i]/pf.coordinate(i,1)/pf.coordinate(i,1))<<std::endl;
    }

    Helper_Func Hs(pf,e2,3);
    Hs.save("h.h5");
    
    std::vector<std::vector<double> > gwin;
    gwin.push_back(pf.grid().gg(0));
    //    std::vector<double> mmt(20,0);
    //    for(int i=0;i<mmt.size();i++) mmt[i]=0.1+0.2*i;
    gwin.push_back(mmt);gwin.push_back(mmt);

    Function w0(gwin);
#pragma omp parallel num_threads(omp_get_max_threads()-2)
    {
#pragma omp for
    for(int i=0;i<w0.size();i++){
    	double h1=0,h2=0;
    	double w=w0.coordinate(i,0);
    	double p1=w0.coordinate(i,1)+w0.coordinate(i,2);
    	double p2=std::abs(w0.coordinate(i,1)-w0.coordinate(i,2));
    	int m=0,n1=0,n2=0;
    	for(int j=0;j<H1.size();j+=H1.grid().gg(1).size())
    	    if(std::abs(H1.coordinate(j,0)-w)<1e-14)
    		m=j;
    	for(int j=0;j<H1.grid().gg(1).size();j++){
    	    if(H1.coordinate(m+j,1)<p1) n1++;
    	    if(H1.coordinate(m+j,1)<p2) n2++;	    
    	}
    	if(n1>0)
    	    h1=( (H1[m+n1]-H1[m+n1-1]) *p1
    		 +H1[m+n1-1]*H1.coordinate(m+n1,1)-H1[m+n1]*H1.coordinate(m+n1-1,1) )
    		/(H1.coordinate(m+n1,1)-H1.coordinate(m+n1-1,1));
    	if(n2>0)
    	    h2=( (H1[m+n2]-H1[m+n2-1]) *p2
    		 +H1[m+n2-1]*H1.coordinate(m+n2,1)-H1[m+n2]*H1.coordinate(m+n2-1,1) )
    		/(H1.coordinate(m+n2,1)-H1.coordinate(m+n2-1,1));
	if(p2<1e-14){
	    // correction for log integral
	    double pp=0,pm=0;
	    if(i%w0.grid().gg(2).size()!=0) pm=w0.coordinate(i, 2)-w0.coordinate(i-1, 2);
	    if((i+1)%w0.grid().gg(2).size()!=0) pp=w0.coordinate(i+1, 2)-w0.coordinate(i, 2);
	    for(int j=1;j<H1.grid().gg(1).size()&&H1.coordinate(m+j, 1)<pp;j++){
		h2+=H1[m+j]*(H1.coordinate(m+j, 1)-H1.coordinate(m+j-1, 1));
	    }
	    for(int j=1;j<H1.grid().gg(1).size()&&H1.coordinate(m+j, 1)<pm;j++){
		h2+=H1[m+j]*(H1.coordinate(m+j, 1)-H1.coordinate(m+j-1, 1));
	    }
	    h2/=(pp+pm);
	}
	w0[i]=(h1-h2)/w0.coordinate(i,1)/w0.coordinate(i,2);
    	// std::cout<<std::setw(8)<<w<<"\t"
    	// 	 <<std::setw(5)<<n1<<"\t"
    	// 	 <<std::setw(8)<<h1<<"\t"
    	// 	 <<std::setw(5)<<n2<<"\t"
    	// 	 <<std::setw(8)<<h2<<"\t"
    	// 	 <<w0[i]<<std::endl;//<<"\t"
       
    }
    }
    std::ofstream wout;
    wout.open("w0.txt");
    for(int i=0;i<w0.size();i++){
    	wout<<std::setw(10)<<w0.coordinate(i,0)<<"\t"
    	   <<std::setw(10)<<w0.coordinate(i,1)<<"\t"
    	   <<std::setw(10)<<w0.coordinate(i,2)<<"\t"	    
    	   <<w0[i]<<std::endl;//<<"\t"
    	    //<<1/(pf.coordinate(i,1)*pf.coordinate(i,1)+pf[i])<<std::endl;//"\t"
    	    //<<1/(1+pf[i]/pf.coordinate(i,1)/pf.coordinate(i,1))<<std::endl;
    }

    file=H5::H5File("w.h5",H5F_ACC_TRUNC);
    H5::Group g2(file.createGroup("/w0"));
    hsize_t dims[1];
    //store freq grid
    dims[0]=w0.grid().lengths(0);
    H5::DataSpace dataspace=H5::DataSpace(1,dims);
    dataset=H5::DataSet(file.createDataSet("/w0/w",
					   H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(w0.grid().gg(0)[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    //store momentum grid
    dims[0]=w0.grid().lengths(1);
    dataspace=H5::DataSpace(1,dims);
    dataset=H5::DataSet(file.createDataSet("/w0/q0",
					   H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(w0.grid().gg(1)[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    dims[0]=w0.grid().lengths(2);
    dataspace=H5::DataSpace(1,dims);
    dataset=H5::DataSet(file.createDataSet("/w0/q1",
					   H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(w0.grid().gg(2)[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();

    //store val grid
    dims[0]=w0.value().size();
    dataspace=H5::DataSpace(1,dims);
    dataset=H5::DataSet(file.createDataSet("/w0/w0",
					   H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(w0.value()[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();

    
    return 0;

}
