#include<iostream>
#include<vector>
#include<cmath>

#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

trng::yarn2 generator;
trng::uniform01_dist<> dist;

double urn(){
    return dist(generator);
}

double const pi=3.14159265;

double const T=0.01;
double const mu=1;

double f(double p, double q, double w, double cost){
    return p*p/(w*w+4*p*p*q*q*cost*cost)/(std::exp((p*p+q*q/4+p*q*cost-mu)/T)+1)/2/pi/pi
}

class Iterator{
public:
    Iterator(int N,double qinit);

    void update();
    void print();

    
private:
    std::vector<double> ws;
    std::vector<double> count;
    double q;
    double total;
    double p;
    double cost;
};

Iterator::Iterator(int N, double qinit)
    :q(qinit),
     total(0),
     count(N,0),
     ws(),
     p(qinit),
     cost(1)
{
    for(int i=0;i<N;i++){
	ws.push_back(2*pi*T*(0.5+i));
    }
}

void Iterator::update(){
    
}


