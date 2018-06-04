#include <reconstructurer.hpp>
reconstructurer::reconstructurer(){
    std::cout <<"Building class reconstructurer with "<<std::endl;
   
}
////////////////////////////////////////////////////7
void reconstructurer::add_points(std::vector<Point>& container){
    for(int i=0; i<container.size(); i++) points_.push_back(container[i]);

}
//////////////////////////////////7
void reconstructurer::build_lateral_surfaces(std::string hname){
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
  const char* interface_off_name_ = (hname +".off").c_str();
  std::ofstream interface_out (interface_off_name_);
  interface_out << *reconstruction_;
  interface_out.close();    
  std::cout<< "end printing off "<<interface_off_name_<<std::endl;
  std::cout<< "Start exporter with file "<< (hname+".off").c_str()<<std::endl;
  exporter stl_exp((hname + ".off").c_str());
//   sharpner sharp(hname+"."+ std::to_string(name));
  std::cout<< "end exporting"<<std::endl;
  
}
