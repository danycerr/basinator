#include <patcher.hpp>
patcher::patcher(horizon   &h1,horizon &h2, int label):name(label){
    std::cout <<"Building class patcher with "<<h1.hname
    << " and "<<h2.hname <<std::endl;
    {
    double xmin=1.e+9;double xmax=1.e-9;
    double ymin=1.e+9;double ymax=1.e-9;
    std::cout<<"Searching min max"<<std::endl;
    this->find_min_max(h1,xmin,xmax,ymin,ymax);

    std::cout<<" xmin "<<xmin
             <<" xmax "<<xmax
             <<" ymin "<<ymin
             <<" ymax "<<ymax<<std::endl;
    std::cout<<"Adding points to the class"<<std::endl;
    add_points(xmin,0,0,h1,points_xmin_);
    add_points(ymin,1,0,h1,points_ymin_);
    add_points(xmax,0,1,h1,points_xmax_);
    add_points(ymax,1,1,h1,points_ymax_);  
    }
    {
    double xmin=1.e+9;double xmax=1.e-9;
    double ymin=1.e+9;double ymax=1.e-9;
    std::cout<<"Searching min max"<<std::endl;
    this->find_min_max(h2,xmin,xmax,ymin,ymax);

    std::cout<<" xmin "<<xmin
             <<" xmax "<<xmax
             <<" ymin "<<ymin
             <<" ymax "<<ymax<<std::endl;
    std::cout<<"Adding points to the class"<<std::endl;
    add_points(xmin,0,0,h2,points_xmin_);
    add_points(ymin,1,0,h2,points_ymin_);
    add_points(xmax,0,1,h2,points_xmax_);
    add_points(ymax,1,1,h2,points_ymax_);
    }
//      for(int i =0; i<points_ymax_.size();i++)
//          std::cout<< points_ymax_[i]<<std::endl;
   build_lateral_surfaces("xmin",points_xmin_);
   build_lateral_surfaces("ymin",points_ymin_);
   build_lateral_surfaces("xmax",points_xmax_);
   build_lateral_surfaces("ymax",points_ymax_);
}
///////////////////////////////////////////////
void patcher::find_min_max(horizon   &h1,
                      double& xmin,
                      double& xmax,
                      double& ymin,
                      double& ymax                      
){
for(int ip=0;ip<h1.get_points().size();ip++){ 
        // std::cout<<h1.get_points()[ip] <<std::endl;
        double x=(h1.get_points()[ip])[0];
        if(x<xmin) xmin=x;
        else if(x>xmax) xmax=x;
        double y=(h1.get_points()[ip])[1];
        if(y<ymin) ymin=y;
        else if(y>ymax) ymax=y;
    }
}
////////////////////////////////////////////////////7
void patcher::add_points(double trashold,int dir,int method,horizon &h,std::vector<Point>& container){
    double eps=20.e+0;
    if(method==0) for(int ip=0;ip<h.get_points().size();ip++){
             if ((h.get_points()[ip])[dir] < trashold+eps )
            container.push_back(h.get_points()[ip]);
        }
    if(method==1) for(int ip=0;ip<h.get_points().size();ip++){
        if ((h.get_points()[ip])[dir] > trashold-eps )
            container.push_back(h.get_points()[ip]);
        }
}
//////////////////////////////////7
void patcher::build_lateral_surfaces(std::string hname,std::vector<Point> &points_){
             reconstruction_.reset(
             new Reconstruction(points_.begin(), points_.end()));
  Timer t;
  t.start();
  Smoother smoother( 10, 200 );
  reconstruction_->increase_scale (4,smoother);
  Mesher mesher( smoother.squared_radius(),
                 false, // Do not separate shells
                 true // Force manifold output
               );                 
  reconstruction_->reconstruct_surface(mesher);
  std::cerr << "Reconstruction done in " << t.time() << " sec." << std::endl;
  std::cout << "generating "<<reconstruction_->number_of_facets()<<std::endl;
  const char* interface_off_name_ = (hname+"."+ std::to_string(name) +".off").c_str();
  std::ofstream interface_out (interface_off_name_);
  interface_out << *reconstruction_;
  interface_out.close();    
  std::cout<< "end printing off"<<std::endl;
  std::cout<< "Start exporter with file "<< (hname+"."+ std::to_string(name) +".off").c_str()<<std::endl;
  exporter stl_exp((hname+"."+ std::to_string(name) +".off").c_str());
//   sharpner sharp(hname+"."+ std::to_string(name));
  std::cout<< "end exporting"<<std::endl;
  
}
