#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

#include <stdio.h>
class horizon{
	public:
		std::vector<double> x,y,z;
		double x_min_=1.e+20,y_min_=1.e+20,x_max_=1.e-20,y_max_=1.e-20;
		int id;
		horizon(int cid, std::vector<double> &points,int nc);
		void dump_off_file();
		void dump_obj_file();
	private:
};
horizon::horizon(int cid, std::vector<double> &points,int nc):id(cid){
	std::cout<< "building the horizon "<< id<<std::endl;
	int nr=points.size()/nc;
	double scale_x=2.e+1*1.e-0;
	double scale_y=2.e+1*1.e-0;
	double scale_z=1.e-0*1.e-0;

	std::cout<<"we have "<< nr<<" points"<<std::endl;
	// x.resize(nr);y.resize(nr);z.resize(nr);
	for(int i=0;i<nr;i++){
		if(((points[i*nc+0]<250 +31) && (points[i*nc+0]>250 -31)&& 
					(points[i*nc+1]<100 +31) && (points[i*nc+1]>100 -31))
		  )
		{ 
			x.push_back(points[i*nc+0]*scale_x);
			y.push_back(points[i*nc+1]*scale_y);
			z.push_back(points[i*nc+2+id]*scale_z);
			if(points[i*nc+0]<x_min_)x_min_=points[i*nc+0];
			if(points[i*nc+0]>x_max_)x_max_=points[i*nc+0];
			if(points[i*nc+1]<y_min_)y_min_=points[i*nc+1];
			if(points[i*nc+1]>y_max_)y_max_=points[i*nc+1];
		}
	}
	std::cout<< "end building the horizon "<< id<<std::endl;    
}
///////////////////////////////////////////
void horizon::dump_off_file()
{
	std::stringstream filename;
	filename<<"horizon_"<<id<<".off";
	std::ofstream ofs (filename.str().c_str(), std::ofstream::out);
	ofs << "OFF"<<std::endl;
	ofs << x.size() <<" 0 0"<<std::endl<<std::endl;
	for(int i=0;i<x.size();i++){
		ofs << x[i] << " "
			<< y[i] << " "
			<< z[i] << std::endl;
	}
	// 
	ofs.close();
}

///////////////////////////////////////////
void horizon::dump_obj_file()
{
	std::stringstream filename;
	filename<<"horizon_"<<id<<".obj";
	FILE * pFile;
	pFile = fopen (filename.str().c_str(),"w");
	fprintf (pFile, "o Plane\n");
	for(int i=0;i<x.size();i++){

		fprintf (pFile, "v %f %f %f\n",x[i],y[i],z[i]);
	}
	fclose (pFile);
}
int main(int argc, char* argv[])
{
	std::cout<<"Coverter"<<std::endl;
	std::ifstream myfile (argv[1]);
	std::string line;
	std::vector<double> points;
	int nc,nr;
	nr=0;
	int buf;
	myfile>>buf;
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			std::vector<double> intline;
			std::istringstream iss(line);
			double a;nc=0;
			while(iss>>a){
				points.push_back(a);
				nc++;
			}
			nr++;
		}
		myfile.close();
	}
	std::cout<<"there are "<< nc-2 << "horizons"
		<<" and "<< nr << " points"<<std::endl;
	std::vector<horizon> hrzs;
	for(int id=0;id<nc-2; id++){
		horizon hbuf(id,points,nc);
		hrzs.push_back(hbuf);
		hrzs[id].dump_off_file();
		hrzs[id].dump_obj_file();
	}

}
