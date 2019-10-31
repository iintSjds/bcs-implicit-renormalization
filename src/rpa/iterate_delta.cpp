# include "../function/function.hpp"
# include <iostream>
# include <ostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <omp.h>

# include <algorithm>

# include <string>
# include "H5Cpp.h"

#include <csignal>
#include <unistd.h>
#include <exception>

#pragma omp declare reduction(vecdplus: std::vector<double>:            \
                              std::transform(omp_out.begin(),omp_out.end(),omp_in.begin(),omp_out.begin(), \
                                             std::plus<double>())) initializer(omp_priv=omp_orig)

# include <trng/yarn2.hpp>
# include <trng/uniform01_dist.hpp>
trng::yarn2 gen;
const double pi=3.141592653;

class sigint:public std::exception{
public:
  const char* what() const throw(){
    return "keyboard interrupt";
  }
};

void signalHandler( int signum ) {
  std::cout << "Interrupt signal (" << signum << ") received.\n";
  // cleanup and close up stuff here
  // terminate program
  throw sigint();
}

class Extrapolation{
  // provide extrapolate for a point,
  // f=p1*f(n1)+p2*f(n2)
public:
  Extrapolation(int i,int j,double a,double b){n1=i;n2=j;p1=a;p2=b;};
  int operator[](int i)const{return (i==0)?n1:n2;}
  double operator()(int i)const{return (i==0)?p1:p2;}
private:
  int n1,n2;
  double p1,p2;
};

inline bool exist_file(const std::string& name){
  std::ifstream f(name.c_str());
  return f.good();
}

class Freq_Extrapolator{
  //provide a table of extrapolation for all frequency points
  //from a pair of freq of function delta (w1,w2)
  //to a point of freq of function W0, (w1+w2) or (w1-w2)
public:
  Freq_Extrapolator(std::vector<double> w,std::vector<double> v,bool is_plus);
  Extrapolation operator()(int i,int j)const{return Exp[v_size*i+j];}
private:
  int w_size,v_size;
  std::vector<Extrapolation> Exp;
};

Freq_Extrapolator::Freq_Extrapolator(std::vector<double> w,std::vector<double> v,bool is_plus){
  for(int i=0;i<v.size();i++){
    // for all points of w
    for(int j=0;j<v.size();j++){
	    int m=0,n=w.size()-1;
	    double a,b;
	    // find v1+v2 or |v1-v2|
	    double freq=(is_plus)?(v[i]+v[j]):(std::abs(v[i]-v[j]));
	    for(int k=0;k<w.size();k++){
        // find two nearest points of freq
        if(w[k]<=freq) m=k; 
        if(w[w.size()-k-1]>freq) n=w.size()-k-1;
	    }
	    if(n==0) n++;
	    if(m==w.size()-1) m--;

	    // calculate extrapolation coefficients and store
	    a=(freq-w[n])/(w[m]-w[n]);
	    b=(freq-w[m])/(w[n]-w[m]);
	    Exp.push_back(Extrapolation(m,n,a,b));
    }
  }
  w_size=w.size();
  v_size=v.size();
}

class Iterator{
public:
  Iterator(double T_,double mu_,double m_, Function w0_ ,
           std::vector<double> v,double wc,double kc);
  double func(int m,int n,int k,int p);
  double update0(double shift,int N);
  void update1(int M);
  void update1_mc(int M,int N);

  Function get_delta(){return delta;};
  Function get_area(){return area;}

  void save_delta(std::string filename);
  bool load_delta(std::string filename);
private:
  double T,mu,mass;
  double wc,kc,kf;//running cut off for freq and mmt,0~wc and kf-kc~kf+kc
  Function delta;
  Function w0;
  Freq_Extrapolator exp_plus,exp_minus;
  Function area;//store integrate area for delta
  std::vector<double> sum1;
  double count;
};

void Iterator::save_delta(std::string filename){
  H5::H5File file;
  //if(!exist_file(filename)){
  file=H5::H5File(filename,H5F_ACC_TRUNC);
  H5::Group g1(file.createGroup("/delta"));
	//}
	//else{
	//file.openFile(filename,H5F_ACC_RDWR);
	//}
  hsize_t dims[1];
  dims[0]=delta.grid().lengths(0);
  H5::DataSpace dataspace=H5::DataSpace(1,dims);
  H5::DataSet dataset(file.createDataSet("/delta/w",
                                         H5::PredType::IEEE_F64LE,dataspace));
  dataset.write(&(delta.grid().gg(0)[0]),H5::PredType::IEEE_F64LE);
  dataspace.close();
  dataset.close();
  //store momentum grid
  dims[0]=delta.grid().lengths(1);
  dataspace=H5::DataSpace(1,dims);
  dataset=H5::DataSet(file.createDataSet("/delta/q",
                                         H5::PredType::IEEE_F64LE,dataspace));
  dataset.write(&(delta.grid().gg(1)[0]),H5::PredType::IEEE_F64LE);
  dataspace.close();
  dataset.close();
  //store val grid
  dims[0]=delta.size();
  dataspace=H5::DataSpace(1,dims);
  dataset=H5::DataSet(file.createDataSet("/delta/delta",
                                         H5::PredType::IEEE_F64LE,dataspace));
  dataset.write(&(delta[0]),H5::PredType::IEEE_F64LE);
  dataspace.close();
  dataset.close();
  file.close();
}

bool Iterator::load_delta(std::string filename){
  H5::H5File file;
  if(!exist_file(filename)) return false;
  file.openFile(filename,H5F_ACC_RDWR);
  H5::DataSet dataset;

  //read value of delta
  hsize_t dims_out[1];
  dataset=file.openDataSet("/delta/delta");
  dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
  std::vector<double> dval(dims_out[0],0);
  dataset.read(&(dval[0]),H5::PredType::IEEE_F64LE);

  if(delta.size()!=dval.size()) return false;
  for(int i=0;i<delta.size();i++){
    delta[i]=dval[i];
    sum1[i]=dval[i];
  }
  count=0;
  return true;
}

Iterator::Iterator(double T_,double mu_,double m_,
                   Function w0_ ,std::vector<double> v,
                   double wc_,double kc_)
  :w0(w0_),T(T_),mu(mu_),mass(m_),wc(wc_),kc(kc_),kf(std::sqrt(mu_/2/m_)),
   exp_plus(w0_.grid().gg(0),v,true),exp_minus(w0_.grid().gg(0),v,false)
  ,delta(v,w0_.grid().gg(1)),area(v,w0_.grid().gg(1)),
   count(0),sum1(v.size()*w0_.grid().gg(1).size(),0)

{
  //make a table for extrapolating freq
  exp_plus=Freq_Extrapolator(w0.grid().gg(0),v,true);
  exp_minus=Freq_Extrapolator(w0.grid().gg(0),v,false);    
  // delta init value, could be random
  // random init
  trng::uniform01_dist<> uniform;
  for(int i=0;i<delta.size();i++){
    delta[i]=1;//2*uniform(gen)-1;
    sum1[i]=delta[i];
  }
    
  for(int i=0;i<area.size();i++){
    int m=i/area.grid().gg(1).size(),k=i%area.grid().gg(1).size();
    double freq_area=0,mmt_area=0;
    if(m!=area.grid().gg(0).size()-1)
	    freq_area=area.grid().gg(0)[m+1]-area.grid().gg(0)[m];
    else
	    freq_area=0;//(area.grid().gg(0)[m]-area.grid().gg(0)[m-1])/2;
    if(k==0)
	    mmt_area=(area.grid().gg(1)[k]+area.grid().gg(1)[k+1])/2;
    else if(k==area.grid().gg(1).size()-1)
	    mmt_area=(area.grid().gg(1)[k]-area.grid().gg(1)[k-1])/2;
    else
	    mmt_area=(area.grid().gg(1)[k+1]-area.grid().gg(1)[k-1])/2;
    //std::cout<<freq_area<<"\t"<<mmt_area<<std::endl;
    area[i]=freq_area/2/pi/T*1;
  }
  std::cout<<"cut off:"<<std::endl
           <<"freq:"<<wc<<std::endl
           <<"mmt:"<<kf+kc<<std::endl;
}


double Iterator::func(int m,int n,int k,int p){
  // return T*p^2/(2*pi)^2*GG(wn,p)*W0(wm,wn,k,p)
  int psize=w0.grid().gg(1).size();
  int psize_sq=psize*psize;
  double mmt=delta.grid().gg(1)[p];
  double mmt_sq=mmt*mmt;
  double freq=delta.grid().gg(0)[n];
  double freq0=delta.grid().gg(0)[m];
  return -(exp_plus(m,n)(0)*w0[exp_plus(m,n)[0]*psize_sq+k*psize+p]
    	     +exp_plus(m,n)(1)*w0[exp_plus(m,n)[1]*psize_sq+k*psize+p]
    	     +exp_minus(m,n)(0)*w0[exp_minus(m,n)[0]*psize_sq+k*psize+p]
    	     +exp_minus(m,n)(1)*w0[exp_minus(m,n)[1]*psize_sq+k*psize+p]
    	     )
    /(freq*freq+(mmt_sq/2/mass-mu)*(mmt_sq/2/mass-mu))
    *T*mmt_sq/4/pi/pi;
  //return -0.7*((freq0+freq)*(freq0+freq)/((freq0+freq)*(freq0+freq)+0.25)
  // 	     +(freq0-freq)*(freq0-freq)/((freq0-freq)*(freq0-freq)+0.25))
  // 	/(freq*freq+(mmt*mmt/2/mass-mu)*(mmt*mmt/2/mass-mu))
  // 	*T*(mmt/2/pi)*(mmt/2/pi);
  // return -2.0*((freq0+freq)*(freq0+freq)/((freq0+freq)*(freq0+freq)+0.01)
  // 	     +(freq0-freq)*(freq0-freq)/((freq0-freq)*(freq0-freq)+0.01))
  // 	/freq/psize
  // 	*T;
}

double Iterator::update0(double shift,int N){
  Function newdelta0(delta.grid());
  std::vector<double> WGGdelta1(delta.size(),0);
  int psize=delta.grid().gg(1).size();
  int fsize=delta.grid().gg(0).size();
  double lambda=0;
  double norm=0,kn=0;
#pragma omp for
  for(int m=0;m<fsize;m++){
    for(int k=0;k<psize;k++){
	    if((delta.grid().gg(1)[k]<=kf+kc)&&(delta.grid().gg(0)[m]<=wc)){
        for(int n=0;n<fsize;n++){
          for(int p=0;p<psize-1;p++){
            if((delta.grid().gg(1)[p]>kf+kc)
               ||(delta.grid().gg(0)[n]>wc)){
              WGGdelta1[m*psize+k]
                +=(func(m,n,k,p)
                   *delta[n*psize+p]*area[n*psize+p]
                   +func(m,n,k,p+1)
                   *delta[n*psize+p+1]*area[n*psize+p+1])
                *(delta.coordinate(n*psize+p+1,1)
                  -delta.coordinate(n*psize+p,1))/2.0;
            }
          }
        }
	    }
    }
  }
  for(int i=0;i<N;i++){
    norm=0;
    kn=0;
#pragma omp parallel num_threads(omp_get_max_threads()-2)
    {
#pragma omp for
	    for(int m=0;m<fsize;m++){
        for(int k=0;k<psize;k++){
          if((delta.grid().gg(1)[k]<=kf+kc)
             &&(delta.grid().gg(0)[m]<=wc)){
            newdelta0[m*psize+k]=0;
            for(int n=0;n<fsize;n++){
              for(int p=0;p<psize-1;p++){
                if((delta.grid().gg(1)[p]<=kf+kc)
                   &&(delta.grid().gg(0)[n]<=wc)){

                  newdelta0[m*psize+k]
                    +=(func(m,n,k,p)
                       *delta[n*psize+p]*area[n*psize+p]
                       +func(m,n,k,p+1)
                       *delta[n*psize+p+1]*area[n*psize+p+1])
                    *(delta.coordinate(n*psize+p+1,1)
                      -delta.coordinate(n*psize+p,1))/2.0;
                }
              }
            }
            newdelta0[m*psize+k]+=
              shift*delta[m*psize+k]+WGGdelta1[m*psize+k];
          }
        }
	    }
    }
	
    // for(int m=0;m<fsize;m++){
    //     for(int k=0;k<psize;k++){
    // 	norm+=newdelta0[m*psize+k]
    // 	    *area[m*psize+k]
    // 	    *(delta.coordinate(m*psize+((k==psize-1)?(k-1):k)+1,1)
    // 	      -delta.coordinate(m*psize+k,1));
    //     }
    // }
    norm=newdelta0[0];
    for(int m=0;m<fsize;m++){
	    for(int k=0;k<psize;k++){
        if((delta.grid().gg(1)[k]<=kf+kc)&&(delta.grid().gg(0)[m]<=wc)){
          delta[m*psize+k]=newdelta0[m*psize+k]
            /norm;
        }
	    }
    }
    //std::cout<<"sample delta:"<<delta[fsize*psize-100]<<std::endl;
    lambda=norm-shift;
  }

  return lambda;
    
}

void Iterator::update1(int M){
  std::vector<double> sum1(delta.size(),0);
  std::vector<double> newsum1(delta.size(),0);
  std::vector<double> WGGdelta0(delta.size(),0);
  int psize=delta.grid().gg(1).size();
  int fsize=delta.grid().gg(0).size();
#pragma omp parallel num_threads(omp_get_max_threads()-2)
  // init
#pragma omp for
  for(int m=0;m<fsize;m++){
    for(int k=0;k<psize;k++){
	    if((delta.grid().gg(1)[k]>kf+kc)||(delta.grid().gg(0)[m]>wc)){
        for(int n=0;n<fsize;n++){
          for(int p=0;p<psize-1;p++){
            if((delta.grid().gg(1)[p]<=kf+kc)
               &&(delta.grid().gg(0)[n]<=wc)){
              WGGdelta0[m*psize+k]
                +=(func(m,n,k,p)
                   *delta[n*psize+p]*area[n*psize+p]
                   +func(m,n,k,p+1)
                   *delta[n*psize+p+1]*area[n*psize+p+1])
                *(delta.coordinate(n*psize+p+1,1)
                  -delta.coordinate(n*psize+p,1))/2.0;
            }
          }
        }
        newsum1[m*psize+k]=sum1[m*psize+k]=delta[m*psize+k];
	    }
    }
  }
  // calculate with accumulation
#pragma omp for
  for(int i=1;i<M+1;i++){
    for(int m=0;m<fsize;m++){
	    for(int k=0;k<psize;k++){
        if((delta.grid().gg(1)[k]>kf+kc)||(delta.grid().gg(0)[m]>wc)){
          //newdelta1[m*psize+k]=0;
          for(int n=0;n<fsize;n++){
            for(int p=0;p<psize-1;p++){
              if((delta.grid().gg(1)[p]>kf+kc)
                 ||(delta.grid().gg(0)[n]>wc)){
                newsum1[m*psize+k]
                  +=(func(m,n,k,p)
                     *sum1[n*psize+p]*area[n*psize+p]
                     +func(m,n,k,p+1)
                     *sum1[n*psize+p+1]*area[n*psize+p+1])
                  *(delta.coordinate(n*psize+p+1,1)
                    -delta.coordinate(n*psize+p,1))/2.0/i;
              }
            }
          }
          sum1[m*psize+k]=newsum1[m*psize+k]+WGGdelta0[m*psize+k];
          newsum1[m*psize+k]=sum1[m*psize+k];
        }
	    }
    }
    std::cout<<"sample sum:"<<sum1[fsize*psize-100]<<std::endl;
  }
  // update delta1
#pragma omp for
  for(int m=0;m<fsize;m++){
    for(int k=0;k<psize;k++){
	    if((delta.grid().gg(1)[k]>kf+kc)||(delta.grid().gg(0)[m]>wc)){
        delta[m*psize+k]=0;
        for(int n=0;n<fsize;n++){
          for(int p=0;p<psize-1;p++){
            if((delta.grid().gg(1)[p]>kf+kc)
               ||(delta.grid().gg(0)[n]>wc)){
              delta[m*psize+k]
                +=(func(m,n,k,p)
                   *sum1[n*psize+p]*area[n*psize+p]
                   +func(m,n,k,p+1)
                   *sum1[n*psize+p+1]*area[n*psize+p+1])
                *(delta.coordinate(n*psize+p+1,1)
                  -delta.coordinate(n*psize+p,1))/2.0/(M+1);
            }
          }
        }
        delta[m*psize+k]+=WGGdelta0[m*psize+k];
	    }
    }
  }
  //std::cout<<"sample delta:"<<delta[fsize*psize-100]<<std::endl;
  // for(int m=0;m<fsize;m++){
  // 	for(int k=0;k<psize;k++){
  // 	    if((delta.grid().gg(1)[k]>kf+kc)||(delta.grid().gg(0)[m]>wc)){
  // 		delta[m*psize+k]=newdelta1[m*psize+k];
  // 	    }
  // 	}
  // }
}

void Iterator::update1_mc(int M, int N){
  //    std::vector<double> sum1(delta.value());
  //    double count=0;
  trng::yarn2 mc;
  trng::uniform01_dist<> uniform;
  int psize=delta.grid().gg(1).size();
  int fsize=delta.grid().gg(0).size();
  double all=delta.size()*fsize*(psize-1);
  for(int i=0;i<M;i++){
    double multi=1;
    double subcount=0;
    std::vector<double> subsum1(delta.size(),0);
#pragma omp parallel num_threads(omp_get_max_threads()-2)	\
  reduction(+:count) reduction(vecdplus:subsum1)
    for(int j=0;j<N;j++){
	    mc.split(omp_get_num_threads(),omp_get_thread_num());
	    // propose a new point
	    int m=uniform(mc)*fsize;
	    int n=uniform(mc)*fsize;
	    int k=uniform(mc)*psize;
	    int p=uniform(mc)*(psize-1);
	    subcount+=multi;
	    if((delta.grid().gg(1)[k]>kf+kc)||(delta.grid().gg(0)[m]>wc)){
        if((delta.grid().gg(1)[p]>kf+kc)||(delta.grid().gg(0)[n]>wc)){
          subsum1[m*psize+k]+=multi*(func(m,n,k,p)
                                     *(sum1[n*psize+p]+subsum1[n*psize+p])
                                     *area[n*psize+p]
                                     +func(m,n,k,p+1)
                                     *(sum1[n*psize+p+1]+subsum1[n*psize+p])
                                     *area[n*psize+p+1])
            *(delta.coordinate(n*psize+p+1,1)
              -delta.coordinate(n*psize+p,1))/2.0/
            (count+subcount+all)*all;
        }
        else{
          subsum1[m*psize+k]+=multi*(func(m,n,k,p)
                                     *delta[n*psize+p]*area[n*psize+p]
                                     +func(m,n,k,p+1)
                                     *delta[n*psize+p+1]*area[n*psize+p+1])
            *(delta.coordinate(n*psize+p+1,1)
              -delta.coordinate(n*psize+p,1))/2.0;
        }
	    }
    }
#pragma omp critical
    count+=subcount;
    std::transform(sum1.begin(),sum1.end(),
                   subsum1.begin(),sum1.begin(),std::plus<double>());
  }
  // update delta1
#pragma omp for
  for(int m=0;m<fsize;m++){
    for(int k=0;k<psize;k++){
	    if((delta.grid().gg(1)[k]>kf+kc)||(delta.grid().gg(0)[m]>wc)){
        delta[m*psize+k]=sum1[m*psize+k]/(count+all)*all;
	    }
    }
  }

}

int main(){
  std::string line;
  std::ifstream infile("delta.in");
  double T=0.5,mu=1.0,m=0.5,rs=1.0;
  std::vector<double> v;
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
	    if(cmd.compare("set v")==0){
        while(std::getline(infile,line)&&line.compare("end grid")!=0){
          v.push_back(std::strtod(line.c_str(),NULL));
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
  file.openFile("w.h5",H5F_ACC_RDWR);
  // read freq grid
  dataset=file.openDataSet("/w0/w");
  hsize_t dims_out[1];
  dataset.getSpace().getSimpleExtentDims(dims_out,NULL);

  std::vector<double> freq(dims_out[0],0);
  dataset.read(&(freq[0]),H5::PredType::IEEE_F64LE);

  //read momentum grid
  dataset=file.openDataSet("/w0/q0");
  dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
  std::vector<double> q(dims_out[0],0);
  dataset.read(&(q[0]),H5::PredType::IEEE_F64LE);

  //read value of w0
  dataset=file.openDataSet("/w0/w0");
  dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
  std::vector<double> wval(dims_out[0],0);
  dataset.read(&(wval[0]),H5::PredType::IEEE_F64LE);

  //restore w0
  std::vector<std::vector<double> > g_in;
  g_in.push_back(freq);g_in.push_back(q);g_in.push_back(q);
  Function w0(g_in);
  for(int i=0;i<w0.size();i++) w0[i]=wval[i];

  // Freq_Extrapolator exp_plus(freq,v,true);
  // std::cout<<v[v.size()-1]<<"\t"<<v[v.size()-1]<<"\t"
  // 	     <<freq[exp_plus(v.size()-1,v.size()-1)[0]]<<"\t"
  // 	     <<exp_plus(v.size()-1,v.size()-1)(0)<<std::endl;

  signal(SIGINT, signalHandler);

    
  Iterator it(T,mu,m,w0,v,100,100);
  try{
    if(exist_file("delta.h5")){
	    bool suc=it.load_delta("delta.h5"); 
	    std::cout<<std::string("load ")+std::string(suc?"success":"fail")
               <<std::endl;
    }
    for(int i=0;i<50;i++) {
	    //it.update1_mc(100,100000);
	    std::cout<<it.update0(6.0,3)<<std::endl;
	    it.save_delta("delta.h5");
    }
    std::cout<<it.update0(6.0,1)<<std::endl;
    it.save_delta("delta.h5");
    Function delta(it.get_delta()),area(it.get_area());
    std::ofstream dout;
    dout.open("delta.txt");
    for(int i=0;i<delta.size();i++){
	    dout<<std::setw(10)<<delta.coordinate(i,0)<<"\t"
          <<std::setw(10)<<delta.coordinate(i,1)<<"\t"
          <<delta[i]<<std::endl;//<<"\t"
      //<<1/(pf.coordinate(i,1)*pf.coordinate(i,1)+pf[i])<<std::endl;//"\t"
      //<<1/(1+pf[i]/pf.coordinate(i,1)/pf.coordinate(i,1))<<std::endl;
    }
    std::cout<<"done"<<std::endl;
  }
  catch(sigint& e){
    std::cout<<"saving files"<<std::endl;
    it.save_delta("delta.h5");
    std::cout<<e.what()<<std::endl;
    exit(SIGINT);
  }
    
  return 0;
}
