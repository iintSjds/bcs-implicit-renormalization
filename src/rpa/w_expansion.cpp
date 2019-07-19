# include "../function/function.hpp"
# include <iostream>
# include <ostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <omp.h>
# include <string>
# include "H5Cpp.h"

const double pi=3.141592653;

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
    	    if(std::abs(H1.coordinate(j,0)-w)<1e-9)
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
