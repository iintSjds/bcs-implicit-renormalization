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

#include <boost/math/special_functions/lambert_w.hpp>

#include "./helper.hpp"

#pragma omp declare reduction(vecdplus: std::vector<double>:            \
                              std::transform(omp_out.begin(),omp_out.end(),omp_in.begin(),omp_out.begin(), \
                                             std::plus<double>())) initializer(omp_priv=omp_orig)

# include <trng/lcg64_shift.hpp>
# include <trng/uniform01_dist.hpp>
# include <trng/uniform_int_dist.hpp>
trng::lcg64_shift mc;

const double pi=3.141592653;
const double k_lower_cut=1e-30;


class sigint:public std::exception{
public:
  const char* what() const throw(){
    return "keyboard interrupt";
  }
};

class gotnan:public std::exception{
public:
  const char* what() const throw(){
    return "nan appear in calculation";
  }
};

void signalHandler( int signum ) {
  std::cout << "Interrupt signal (" << signum << ") received.\n";
  // cleanup and close up stuff here
  // terminate program
  throw sigint();
}

inline bool exist_file(const std::string& name){
  std::ifstream f(name.c_str());
  return f.good();
}

class K_Distributor{
public:
  K_Distributor(double kf_,double kc_):kf(kf_),kc(kc_),kmin(1e-9) {};
  K_Distributor(double kf_,double kc_,double kmin_):kf(kf_),kc(kc_),kmin(kmin_){};
  double dist(double p);// receive a number from 0 to 1, propose a k with dist
  double prob(double k);// return probability of proposing k
private:
  double kf;
  double kc;
  double kmin;
};

double K_Distributor::dist(double p){
  int region=p*3;
  double pp=p*3-region;
  if(region==0) return kmin*std::pow((0.5*kf+kmin)/kmin,pp)-kmin;
  else if(region==1) return kf+kmin-kmin*std::pow((0.5*kf+kmin)/kmin,pp);
  else return kmin*std::pow((kmin+kc-kf)/kmin,pp)+kf-kmin;
  // double k=-pp*(1+std::log(2))/2.0
  // 	/boost::math::lambert_wm1(-pp*(1+std::log(2))/2/2.718281828);
  // if(region==0) return k*kf;
  // else if(region==1) return kf-k*kf;
  // else return kf+k*2*(kc-kf);
}

double K_Distributor::prob(double k){
  if(k<kf){
    double kk=kf/2.0-std::abs(k-kf/2.0);
    //return (std::log(kf/kk))*2.0/kf/(1+std::log(2));
    return 1.0/(kk+kmin)/std::log(0.5*kf/kmin+1);
  }
  else{
    //    	double kk=k-kf;
    //    	double kff=2*(kc-kf);
    //    	return (std::log(kff/kk))*2.0/kff/(1+std::log(2));
    return 1.0/(k-kf+kmin)/std::log(1+kc/kmin+kf/kmin);
  }
  //    return 1.0/kc;
}

class W_Distributor{
public:
  W_Distributor(double mu_,double wc_,double T_):mu(mu_),wc(wc_),T(T_),wm(pi*T){};
  W_Distributor(double mu_,double wc_,double T_,double wm_):mu(mu_),wc(wc_),T(T_),wm(wm_){};
  long long int dist(double p);// receive a number from 0 to 1, propose a w with dist
  double prob(double w);// return probability of proposing w
  double prob(double w1,double w2);// return p of proposing between w1<w2
private:
  double mu;
  double wc;
  double T;
  double wm;
};

long long int W_Distributor::dist(double p){
  return (std::pow(wc/wm,p)-1.0)/2.0;
}

double W_Distributor::prob(double w){
  return std::log((w+2*pi*T)/(w))/std::log(wc/wm);
}

double W_Distributor::prob(double w1,double w2){
  return std::abs(std::log((w2)/(w1))/std::log(wc/wm));
}

class Iterator{
public:
  Iterator(double T_,double mu_,double m_,
           Helper_Func H_,std::vector<double> v,std::vector<double> k,
           double wc,double kc);
  double func(double w1,double w2,double k1,double k2);
  double delta(double w,double k);

  double update0(double shift,int N,unsigned long seed=0);
  void update1(int N);

  void save_delta(std::string filename);
  bool load_delta(std::string filename);
  double get_count0() const{return count0;}
  void print(std::string filename);
private:
  int order;
  double T,mu,mass;
  double wc,kc;
  Function delta_count;
  Function area;
  Helper_Func H;
  double count0;
  double count1;
};

void Iterator::print(std::string filename){
  std::ofstream dout;
  dout.open(filename);
  for(int i=0;i<delta_count.size();i++){
    double w=delta_count.coordinate(i,0),k=delta_count.coordinate(i,1);
    dout<<std::setw(10)<<w<<"\t"
        <<std::setw(10)<<k<<"\t"
        <<std::setw(10)<<delta(w,k)<<"\t"
        <<std::setw(10)<<area[i]<<"\t"
        <<std::endl;//<<"\t"
    //<<1/(pf.coordinate(i,1)*pf.coordinate(i,1)+pf[i])<<std::endl;//"\t"
    //<<1/(1+pf[i]/pf.coordinate(i,1)/pf.coordinate(i,1))<<std::endl;
  }
}

void Iterator::save_delta(std::string filename){
  H5::H5File file;
  //if(!exist_file(filename)){
  file=H5::H5File(filename,H5F_ACC_TRUNC);
  try{
    H5::Group g1(file.createGroup("/delta"));
    //}
    //else{
    //file.openFile(filename,H5F_ACC_RDWR);
    //}
    hsize_t dims[1];
    dims[0]=delta_count.grid().lengths(0);
    H5::DataSpace dataspace=H5::DataSpace(1,dims);
    H5::DataSet dataset(file.createDataSet("/delta/w",
                                           H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(delta_count.grid().gg(0)[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    //store momentum grid
    dims[0]=delta_count.grid().lengths(1);
    dataspace=H5::DataSpace(1,dims);
    dataset=H5::DataSet(file.createDataSet("/delta/q",
                                           H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(delta_count.grid().gg(1)[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    //store val grid
    dims[0]=delta_count.size();
    dataspace=H5::DataSpace(1,dims);
    dataset=H5::DataSet(file.createDataSet("/delta/delta",
                                           H5::PredType::IEEE_F64LE,dataspace));
    std::vector<double> delta_val(delta_count.size(),0);
    for(int i=0;i<delta_count.size();i++){
	    delta_val[i]=delta(delta_count.coordinate(i,0),delta_count.coordinate(i,1));
    }
    dataset.write(&(delta_val[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    file.close();
  }
  catch(sigint){
    //TBA
    file.close();
    throw;
  }
}

bool Iterator::load_delta(std::string filename){
  H5::H5File file;
  double N=1e9;
  if(!exist_file(filename)) return false;
  file.openFile(filename,H5F_ACC_RDWR);
  H5::DataSet dataset;

  //read value of delta
  hsize_t dims_out[1];
  dataset=file.openDataSet("/delta/delta");
  dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
  std::vector<double> dval(dims_out[0],0);
  dataset.read(&(dval[0]),H5::PredType::IEEE_F64LE);

  if(delta_count.size()!=dval.size()) {
    file.close();
    dataset.close();
    return false;
  }

  for(int i=0;i<delta_count.size();i++){
    delta_count[i]=dval[i]*N;
  }
  count0=count1=dval[0]*N;
  dataset.close();
  file.close();
  return true;
}

double Iterator::update0(double shift,int N,unsigned long seed){
  //mc.seed(int(count0));
  if(seed!=0) mc.seed((seed));
  trng::uniform01_dist<> urn;
  // trng::uniform_int_dist uin1(0,int((wc-T*pi)/2/pi/T));
  // trng::uniform_int_dist uin2(0,int((delta_count.grid().gg(0).back())/2/pi/T)+1);
  trng::uniform_int_dist uinp(0,2);
  int psize=delta_count.grid().gg(1).size();
  double pmax=delta_count.grid().gg(1).back()
    *((delta_count.grid().gg(0).back()-pi*T)/2/pi/T);
  double c0=0,c1=0;

  double temp_count0=0;
  std::vector<double> temp_delta_count(delta_count.size(),0);

  int part=1;
  double kmin=delta_count.grid().gg(1).front();
  K_Distributor kd1(std::sqrt(2*mu*mass),kc,kmin);
  K_Distributor kd2(std::sqrt(2*mu*mass),delta_count.grid().gg(1).back(),kmin);
  W_Distributor wd1(mu,wc,T);
  W_Distributor wd2(mu,delta_count.grid().gg(0).back(),T);
  double w1=pi*T*(2*wd1.dist(urn(mc))+1),w2=pi*T*(2*wd2.dist(urn(mc))+1),
    k1=kd1.dist(urn(mc)),k2=kd2.dist(urn(mc));
  double f0=0,f1=0;
  // std::ofstream fout;
  // fout.open("f.txt",std::ofstream::out | std::ofstream::app);
  int hitcount=0;
  for(int i=0;i<N;i++){

    //	part=uinp(mc);

    w1=pi*T*(2*wd1.dist(urn(mc))+1);
    w2=pi*T*(2*wd2.dist(urn(mc))+1);
    k1=kd1.dist(urn(mc))+k_lower_cut;
    k2=kd2.dist(urn(mc))+k_lower_cut;
    double prob=kd1.prob(k1)*kd2.prob(k2)*wd1.prob(w1)*wd2.prob(w2);

    f1=(func(w1,w2,k1,k2)*delta(w2,k2))/prob;
    f0=((shift*0.01)*delta(w1,k1))/pmax/prob;

    //	    std::cout<<f<<std::endl;
    // update new config to delta_count
    int nw=std::lower_bound(delta_count.grid().gg(0).begin(),
                            delta_count.grid().gg(0).end()-1,
                            w1+0.1*T) - delta_count.grid().gg(0).begin();
    int nk=std::lower_bound(delta_count.grid().gg(1).begin(),
                            delta_count.grid().gg(1).end()-1,
                            k1) - delta_count.grid().gg(1).begin();

    if(nw!=0) nw--;
    // fout<<std::setw(15)<<delta_count.grid().gg(0)[nw]<<"\t"
    //     <<std::setw(15)<<w1<<"\t"
    //     <<std::setw(15)<<delta_count.grid().gg(1)[nk]<<"\t"
    //     <<std::setw(15)<<k1<<"\t"
    //     <<std::setw(15)<<area[psize*nw+nk]//func(w1,w2,k1,k2)
    //     <<std::endl;
    c0+=f0;c1+=f1;

    temp_delta_count[psize*nw+nk]+=(f0+f1)
	    /area[psize*nw+nk];///area[psize*nw+nk];
    if(psize*nw+nk==0){
	    temp_count0+=(f0)/area[0];
	    hitcount++;
    }
 }
  for(int j=0;j<delta_count.size();j++){
    if(delta_count.coordinate(j,0)<wc && delta_count.coordinate(j,1)<kc ){
      delta_count[j] += 99*delta_count[j]/count0*temp_count0
        +temp_delta_count[j];//*temp_sign;
      temp_delta_count[j]=0;
    }
  }
  count0=delta_count[0];
  temp_count0=0;
  return 0.01*shift/c0*c1;
}

void Iterator::update1(int N){
  trng::uniform01_dist<> urn;
  trng::uniform_int_dist uinp(0,2);
  int psize=delta_count.grid().gg(1).size();
  double pmax=delta_count.grid().gg(1).back()
    *((delta_count.grid().gg(0).back()-pi*T)/2/pi/T);
  double c0=0,c1=0;

  count1=1;
  std::vector<double> temp_delta_count(delta_count.size(),0);

  int part=1;
  double kmin=delta_count.grid().gg(1).front();
  K_Distributor kd1(std::sqrt(2*mu*mass),kc,kmin);
  K_Distributor kd2(std::sqrt(2*mu*mass),delta_count.grid().gg(1).back(),kmin);
  W_Distributor wd1(mu,delta_count.grid().gg(0).back(),T,wc);
  W_Distributor wd2(mu,delta_count.grid().gg(0).back(),T);
  double w1=pi*T*(2*wd1.dist(urn(mc))+1),w2=pi*T*(2*wd2.dist(urn(mc))+1),
    k1=kd1.dist(urn(mc)),k2=kd2.dist(urn(mc));
  double f0=0,f1=0;
  int hitcount=0;
  for(int i=0;i<N;i++){
    w1=pi*T*(2*wd1.dist(urn(mc))+1);
    w2=pi*T*(2*wd2.dist(urn(mc))+1);
    k1=kd1.dist(urn(mc))+k_lower_cut;
    k2=kd2.dist(urn(mc))+k_lower_cut;
    double prob=kd1.prob(k1)*kd2.prob(k2)*wd1.prob(w1)*wd2.prob(w2);

    f1=(func(w1,w2,k1,k2)*delta(w2,k2))/prob;
    int nw=std::lower_bound(delta_count.grid().gg(0).begin(),
                            delta_count.grid().gg(0).end()-1,
                            w1+0.1*T) - delta_count.grid().gg(0).begin();
    int nk=std::lower_bound(delta_count.grid().gg(1).begin(),
                            delta_count.grid().gg(1).end()-1,
                            k1) - delta_count.grid().gg(1).begin();

    if(nw!=0) nw--;
    c1+=f1;

    delta_count[psize*nw+nk]+=(f1)
      /area[psize*nw+nk];///area[psize*nw+nk];
    count1++;
  }
  //update recalculated delta1
  for(int j=0;j<delta_count.size();j++){
    if(delta_count.coordinate(j,0)>=wc || delta_count.coordinate(j,1)>=kc){
      delta_count[j]=delta_count[j]/count1;
    }
  }
  count1=1;
  //int temp_sign=temp_delta_count[0]/std::abs(temp_delta_count[0]);
}

Iterator::Iterator(double T_,double mu_,double m_,
                   Helper_Func H_,std::vector<double> v,std::vector<double> k,
                   double wc_,double kc_)
  :T(T_),mu(mu_),mass(m_),
   H(H_),delta_count(v,k),area(v,k),
   wc((v.back()>wc_)?wc_:v.back()),kc((k.back()>(kc_))?kc_:k.back()),
   order(0),
   count0(1),count1(1)
{
  // wc should be set to minimal value exist in grid larger than wc_
  int nwc=std::lower_bound(v.begin(),
                           v.end()-1,
                           wc) - v.begin();
  wc=v[nwc];

  trng::uniform01_dist<> urn;
  W_Distributor wd1(mu,wc,T);
  for(int i=0;i<delta_count.size();i++){
    if((delta_count.coordinate(i,1)<=kc)&&(delta_count.coordinate(i,0)<=0.5)){
	    delta_count[i]=1;//urn(mc)*2-1;
    }
    else{
	    delta_count[i]=1;//urn(mc)*2-1;
    }
    delta_count[0]=1;
    int m=i/area.grid().gg(1).size(),k=i%area.grid().gg(1).size();
    double freq_area=0,mmt_area=0;
    if(m!=area.grid().gg(0).size()-1)
	    freq_area=area.grid().gg(0)[m+1]-area.grid().gg(0)[m];
    // freq_area=wd1.prob(area.grid().gg(0)[m+1],area.grid().gg(0)[m])
    // 	*(area.grid().gg(0).back()-pi*T)/2/pi/T;
    else
	    freq_area=area.grid().gg(0)[m]-area.grid().gg(0)[m-1];
    // freq_area=wd1.prob(area.grid().gg(0)[m],area.grid().gg(0)[m-1])
    // 	*(area.grid().gg(0).back()-pi*T)/2/pi/T;
    //(area.grid().gg(0)[m]-area.grid().gg(0)[m-1])/2;
    if(k==0)//area.grid().gg(1).size()-1)
	    mmt_area=area.grid().gg(1)[k];//-area.grid().gg(1)[k-1];
    else
	    mmt_area=(area.grid().gg(1)[k]-area.grid().gg(1)[k-1]);
    //std::cout<<freq_area<<"\t"<<mmt_area<<std::endl;
    area[i]=freq_area/2/pi/T*mmt_area;
  }
}

double Iterator::delta(double w,double k){
  //    std::cout<<"in delta"<<std::endl;
  int nw=std::upper_bound(delta_count.grid().gg(0).begin(),
                          delta_count.grid().gg(0).end()-1,
                          w+0.1*T) - delta_count.grid().gg(0).begin();
  int nk=std::upper_bound(delta_count.grid().gg(1).begin(),
                          delta_count.grid().gg(1).end()-1,
                          k) - delta_count.grid().gg(1).begin();
  int psize=delta_count.grid().gg(1).size();
  if(nw!=0) nw--;
  // std::cout<<w<<"\t"<<k<<"\t"
  // 	     <<delta_count.grid().gg(0)[nw]<<"\t"
  // 	     <<delta_count.grid().gg(1)[nk]<<"\t"
  // 	     <<delta_count[psize*nw+nk]<<"\t"
  // 	     <<area[psize*nw+nk]<<"\t"
  // 	     <<std::endl;
  if((k<=kc)&&(w<=wc)){
    return delta_count[psize*nw+nk]/count0;
  }
  else{
    return delta_count[psize*nw+nk]/count1;
  }
}

double Iterator::func(double w1,double w2,double k1,double k2){
  //    std::cout<<"in func"<<std::endl;
  double pp=k1+k2,pm=std::abs(k1-k2)+k_lower_cut;
  double vm=std::abs(w1-w2),vp=w1+w2;
  double h=0;
  // int kp=std::lower_bound(H[0].grid().gg(1).begin(),
  // 			    H[0].grid().gg(1).end()-1,
  // 			    pp) - H[0].grid().gg(1).begin();
  // int km=std::lower_bound(H[0].grid().gg(1).begin(),
  // 			    H[0].grid().gg(1).end()-1,
  // 			    pm) - H[0].grid().gg(1).begin();
  // int fp=std::lower_bound(H[0].grid().gg(0).begin(),
  // 			    H[0].grid().gg(0).end()-1,
  // 			    vp+0.01*T) - H[0].grid().gg(0).begin();
  // int fm=std::lower_bound(H[0].grid().gg(0).begin(),
  // 			    H[0].grid().gg(0).end()-1,
  // 			    vm+0.01*T) - H[0].grid().gg(0).begin();
  // if(fp!=0) fp--;
  // if(fm!=0) fm--;
  // int psize=H[0].grid().gg(1).size();
  // if(order==0){
  // 	h=(H[0][psize*fp+kp]+H[0][psize*fm+kp]
  // 	    -H[0][psize*fp+km]-H[0][psize*fm+km])/k1/k2;
  // }

  // double g=1;
  // h+=4*pi*g/k1/k2*(std::log(pp/pm)-
  // 		     0.5*g/(vp*vp+g)*(2*std::log(pp/pm)-
  // 				      std::log((pp*pp+vp*vp+g)/(pm*pm+vp*vp+g))));
  // h+=4*pi*g/k1/k2*(std::log(pp/pm)-
  // 		     0.5*g/(vm*vm+g)*(2*std::log(pp/pm)-
  // 				      std::log((pp*pp+vm*vm+g)/(pm*pm+vm*vm+g))));

  // double result=-h/(w2*w2+(k2*k2/2.0/mass-mu)*(k2*k2/2.0/mass-mu))
  // 	*T*k2*k2/4.0/pi/pi;

  // if(std::isnan(result)){
  // 	std::cout<<"result==nan!"<<std::endl;
  // 	std::cout<<"current state:"<<std::endl;
  // 	std::cout<<"w1="<<w1<<std::endl;
  // 	std::cout<<"w2="<<w2<<std::endl;
  // 	std::cout<<"k1="<<k1<<std::endl;
  // 	std::cout<<"k2="<<k2<<std::endl;
  // 	throw gotnan();
  // }
  // else return result;
  return -8.0*pi*((vp*vp)/(vp*vp+0.25)+(vm*vm)/(vm*vm+0.25))
    /(w2*w2+(k2*k2/2.0/mass-mu)*(k2*k2/2.0/mass-mu))
    *T*k2*k2/4.0/pi/pi;
}


int main(){
  std::string line;
  std::ifstream infile("delta.in");
  double T=0.5,mu=1.0,m=0.5,rs=1.0;
  std::vector<double> v,mmt;
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

  Helper_Func H("./h.h5",e2,3);
  Iterator it(T,mu,m,H,v,mmt,50,100);

  K_Distributor kd(1,10);
  //    std::cout<<it.func(0.1,0.1,0.9,1.1)<<"\t"<<std::endl;
  signal(SIGINT, signalHandler);
  if(it.load_delta("delta.h5")) std::cout<<"load success"<<std::endl;

  double sum_eigen=0;
  for(unsigned long long int i=0;i<200;i++)
    {
	    double eigen=it.update0(3, 10000000, 23*i+19);
      it.update1(10000);
	    sum_eigen+=eigen;
	    std::cout<<eigen<<"\t"
               <<sum_eigen/(i+1)<<"\t"
               <<it.get_count0()<<std::endl;
	    it.save_delta("delta.h5");
    }
  it.print("delta_mc.txt");
  return 0;
}
