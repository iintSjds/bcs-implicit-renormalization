#include "Grid.hpp"
#include "cmath"
#include "Point.hpp"

Grid_1d::Grid_1d(double init, double fin, int size){
    double step=(fin-init)/size;
    double area=std::fabs(step);
    double x=init+0.5*step;
    dim=1;
    for(int i=0;i<size;i++){
	points.push_back(Point(x,area));
	x+=step;
    }
    std::vector<int> len(1,size);
    lens=len;
}

Grid_1d::Grid_1d(std::vector<double> p){
    dim=1;
    std::vector<int> len(1,p.size());
    lens=len;
    if(p.size()>0){
	if(p.size()>1){
	    if(p.size()>2){
		points.push_back(Point(p[0],std::fabs(p[1]-p[0])));
		for(int i=1;i<p.size()-1;i++)
		    points.push_back(Point(p[i],std::fabs(p[i+1]-p[i-1])));
		points.push_back(Point(p[p.size()-1],std::fabs(p[p.size()-1]-p[p.size()-2])));		
	    }
	    else{
		//==2
		double area=std::fabs(p[0]-p[1]);
		points.push_back(Point(p[0],area));
		points.push_back(Point(p[1],area));		
	    }
	}
	else{
	    //==1
	    points.push_back(Point(p[0]));
	}
    }
}

MeshGrid::MeshGrid(Grid& g1, Grid& g2){
    dim = g1.dimension()+g2.dimension();
    for(int i=0;i<g1.dimension();i++) lens.push_back(g1.lengths()[i]);
    for(int i=0;i<g2.dimension();i++) lens.push_back(g2.lengths()[i]);
    int size=g1.size()*g2.size();
    for(int i=0;i<size;i++){
	points.push_back(g1[i/g2.size()]*g2[i%g2.size()]);//point multiply
    }
}

