# include "./helper.hpp"
const double pi=3.141592653;


Helper_Func::Helper_Func(std::string filename,double e2_,int n)
    :e2(e2_)
{
    H5::H5File file;
    file.openFile(filename,H5F_ACC_RDWR);
    H5::DataSet dataset;

    std::vector<std::vector<double> > hs;
    hsize_t dims_out[1];
    
    dataset=file.openDataSet("/h/w");
    dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
    std::vector<double> w(dims_out[0],0);
    dataset.read(&(w[0]),H5::PredType::IEEE_F64LE);
    
    dataset=file.openDataSet("/h/q");
    dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
    std::vector<double> q(dims_out[0],0);
    dataset.read(&(q[0]),H5::PredType::IEEE_F64LE);

    for(int i=0;i<n;i++){
	dataset=file.openDataSet(std::string("/h/h")+std::to_string(2*i+1));
	dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
	std::vector<double> h(dims_out[0],0);
	dataset.read(&(h[0]),H5::PredType::IEEE_F64LE);
	hs.push_back(h);
    }
    Grid g(w,q);
    for(int m=0;m<n;m++){
	Function h(g);
	for(int i=0;i<h.size();i++) h[i]=hs[m][i];
	H.push_back(h);
    }
}

Helper_Func::Helper_Func(Function pf,double e2_,int n)
    :e2(e2_)
{
    for(int m=0;m<n;m++){
	Function h(pf.grid());
#pragma omp parallel num_threads(omp_get_max_threads()-2)
	{
#pragma omp for
	    for(int i=0;i<h.size();i++){
		double val=0;
		for(int j=(i/h.grid().lengths(1))*h.grid().lengths(1);j<i;j++){
		    val+=(std::pow(pf.coordinate(j,1),2*m)
			  /(pf.coordinate(j,1)/4/pi/e2
			    +pf[j]/pf.coordinate(j,1))
			  +std::pow(pf.coordinate(j+1,1),2*m)
			  /(pf.coordinate(j+1,1)/4/pi/e2
			   +pf[j+1]/pf.coordinate(j+1,1)))
			//e2*(pf.coordinate(j,1)+pf.coordinate(j+1,1))
			*(pf.coordinate(j+1,1)-pf.coordinate(j,1))/2;
		}
		h[i]=val;
	    }
	}
	H.push_back(h);
    }
}

void Helper_Func::save(std::string filename){
        H5::H5File file;
    //if(!exist_file(filename)){
    file=H5::H5File(filename,H5F_ACC_TRUNC);
    H5::Group g1(file.createGroup("/h"));
	//}
	//else{
	//file.openFile(filename,H5F_ACC_RDWR);
	//}
    hsize_t dims[1];
    dims[0]=grid().lengths(0);
    H5::DataSpace dataspace=H5::DataSpace(1,dims);
    H5::DataSet dataset(file.createDataSet("/h/w",
					   H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(grid().gg(0)[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    //store momentum grid
    dims[0]=grid().lengths(1);
    dataspace=H5::DataSpace(1,dims);
    dataset=H5::DataSet(file.createDataSet("/h/q",
					   H5::PredType::IEEE_F64LE,dataspace));
    dataset.write(&(grid().gg(1)[0]),H5::PredType::IEEE_F64LE);
    dataspace.close();
    dataset.close();
    //store val grid
    for(int i=0;i<H.size();i++){
	dims[0]=H[i].size();
	dataspace=H5::DataSpace(1,dims);
	dataset=H5::DataSet(file.createDataSet(std::string("/h/h")
					       +std::to_string(2*i+1),
					       H5::PredType::IEEE_F64LE,dataspace));
	dataset.write(&(H[i][0]),H5::PredType::IEEE_F64LE);
	dataspace.close();
	dataset.close();
    }
    file.close();
}
