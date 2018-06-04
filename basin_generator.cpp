#include<basin_generator.hpp>

///////////////////////////////7
horizon::horizon(char * file):hname(file){
    std::cout<<"Building basin class"<<std::endl;
    std::cout<<"horizon "<<file<<std::endl;
    hname.erase(0,hname.find_last_of("/")+1);
    hname.erase(hname.find_last_of("."));
    std::cout<<"name of horizon "<<hname<<std::endl;
    
}
///////////////////////////////////7
void horizon::read_horizons_points(char * file){
    std::cout<<"reading points clouds"<<std::endl;
    std::ifstream in(file);
    std::cerr << "Reading " << std::flush;
    if( !in || !CGAL::read_off_points( in, std::back_inserter( points_ ) ) ) {
       std::cerr << "Error: cannot read file" << std::endl;
       }
  std::cerr << "done: " << points_.size() << " points." << std::endl;
}
////////////////////////////////////////////7
void horizon::reconstruction_surface() {
         reconstruction_.reset(
             new Reconstruction(points_.begin(), points_.end()));
  Timer t;
  t.start();
  Smoother smoother( 10, 100 );
  reconstruction_->increase_scale (4,smoother);
  Mesher mesher( smoother.squared_radius(),
                 false, // Do not separate shells
                 true // Force manifold output
               );                 
  reconstruction_->reconstruct_surface(mesher);
  
//   // Jet smoother
//     reconstruction_->increase_scale<JetSmoother> (4);
//     reconstruction_->reconstruct_surface ();
  std::cerr << "Reconstruction done in " << t.time() << " sec." << std::endl;
  std::cout << "generating "<<reconstruction_->number_of_facets()<<std::endl;
  
 
std::cout<<"Looping on reconstruction points"<<std::endl;
  for (Reconstruction::Point_iterator p=reconstruction_->points_begin(); p!=reconstruction_->points_end(); p++){
rpoints_.push_back(*p);
}
  check_normals();
  print_stl();
  print_interface();
}
//////////////////////////////////////////////////7
void horizon::print_interface(){
    std::cout<< "start printing off of "<<hname <<std::endl;
    Timer t;
    t.start();
    interface_off_name_ = (hname+".off").c_str();
    interface_oriented_off_name_ = ("../interfaces_oriented.off");
    std::ofstream interface_out (interface_off_name_);
    interface_out << *reconstruction_;
    interface_out.close();    
    std::cout<< "end printing off"<<std::endl;
    {
//     sharpner sharp(hname);
    }
    std::ifstream input(interface_oriented_off_name_);
    input >> poly_;
    // Create a vector with only one element: the pointer to the polyhedron.
    std::vector<Polyhedron*> poly_ptrs_vector(1, &poly_);
    int PRINT_MEDIT=0;
    if(PRINT_MEDIT){
    // Create a polyhedral domain, with only one polyhedron,
    // and no "bounding polyhedron", so the volumetric part of the domain will be
    // empty.
    Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
    // Get sharp features
    domain.detect_features(); //includes detection of borders
    // Mesh criteria    
    Mesh_criteria criteria(
                            edge_size = 0.25,
                            facet_angle = 25,
                            facet_size = 0.1,
                            facet_distance = 0.1
                          );
    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
    // Output
    std::ofstream medit_file("../interfaces.mesh");
    c3t3.output_to_medit(medit_file);
    medit_file.close();
    }
    std::cerr << "Interface printing  done in " << t.time() << " sec." << std::endl;
}
/////////////////////////////////////////////////
// function for writing the reconstruction output in the off format
void horizon::print_interface_oriented()
{
  std::ofstream output("../interfaces_oriented.off");
  output << "OFF \n" << reconstruction_->number_of_points() << " "
         << reconstruction_->number_of_facets() << " 0\n\n";
  std::copy(reconstruction_->points_begin(),
            reconstruction_->points_end(),
            std::ostream_iterator<Point>(output,"\n"));
  for( int iface=0; iface<faces_.size();iface++)
      output << "3 " << faces_[iface] << std::endl;
  output << "\n"; 
  output.close();
}
//////////////////////////////////////////////////////////
// function for writing the reconstruction output in the geo format
void horizon::print_stl()
{
//binary file
//     std::string filename="prova";
std::cout<<"stl file printing"<<(hname + ".stl").c_str()<<std::endl;
std::string header_info = "solid " + hname + "-output";
char head[80];
std::strncpy(head,header_info.c_str(),sizeof(head)-1);
char attribute[2] = "0";
uint32_t nTriLong =(uint32_t)  faces_.size();

std::ofstream myfile;

myfile.open((hname + ".stl").c_str(),  std::ios::out | std::ios::binary);
myfile.write(head,sizeof(head));
myfile.write((char*)&nTriLong,sizeof nTriLong);
std::cout<<"size of ntri "<< sizeof nTriLong<<std::endl;
std::cout<<"size of points "<< sizeof (fl_points_[0])[0]<<std::endl;
std::cout<<"size of normals "<< sizeof(normals_[0])[0]<<std::endl;
    //write down every triangle
  for( int iface=0; iface<nTriLong;iface++)
{
    //normal vector coordinates
    // std::cout<<"writing element "<<iface<<std::endl;
    for(int idim=0;idim<3;idim++){
        float buf=(normals_[iface])[idim];
    myfile.write((char*)&buf, sizeof(float));
    // std::cout<<(normals_[iface])[idim]<<std::endl;
    }
    //p1 coordinates
    for(int ipoint =0; ipoint<3;ipoint++)
    for(int idim=0;idim<3;idim++){
        float buf=(float)(fl_points_[((faces_[iface])[ipoint])])[idim];
    myfile.write((char*)&buf, sizeof(float));
    //     std::cout<<(fl_points_[((faces_[iface])[ipoint])])[idim]<<std::endl;
    }


    myfile.write(attribute,2);
}
myfile.close();
}
/////////////////////////////////////////////////////////
void horizon::create_poligonal_domain(char * file){
    std::cout<<"Geneartion of poligonal costrained mesh"<<std::endl;
    // creation of the polygonal container
    domain_off_name_ = (file);
    std::ifstream input(domain_off_name_);
    input >> poly_bounding_;  // Create domain
    domain_costrained_.reset(new Mesh_domain (poly_, poly_bounding_));
    // Get sharp features
    domain_costrained_->detect_features();
    std::cout<<"End of geneartion of poligonal costrained mesh"<<std::endl;
    // Mesh criteria
    Mesh_criteria criteria(edge_size = 0.025,
                         facet_angle = 25, facet_size = 0.05, facet_distance = 0.005,
                         cell_radius_edge_ratio = 3, cell_size = 0.05);
  
    // Mesh generation
    std::cout<<"Geneartion of poligonal costrained mesh"<<std::endl;
    Timer t;
    t.start();
    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_costrained_, criteria,
                                        no_perturb(), no_exude());
    std::cerr << "Generation of 3d mesh costrained in "<< t.time() << " sec." << std::endl;
    // Output
    dump_c3t3(c3t3, "../3dbasin");
}
///////////////////////////////////////
void horizon::check_normals(){
    std::cout<<"Checking normals"<<std::endl;
    std::cout << "number of points" << reconstruction_->number_of_points() <<std::endl;
    int ipoint=0; 
    fl_points_.resize(reconstruction_->number_of_points());
    faces_.resize(reconstruction_->number_of_facets()); // sizing for faces
    normals_.resize(reconstruction_->number_of_facets()); // sizing for normals
    for( Point_iterator it = reconstruction_->points_begin(); 
	it != reconstruction_->points_end(); ++it )
  	{
        double a=it->x();
        (fl_points_[ipoint])[0]=(float) it->x();
        (fl_points_[ipoint])[1]=(float) it->y();
        (fl_points_[ipoint])[2]=(float) it->z();
        // std::cout<<points[ipoint]<<std::endl;
        ipoint++;
    } 
    int iface=0;
  for( Facet_iterator it = reconstruction_->facets_begin(); 
        it != reconstruction_->facets_end(); ++it ){
        // Point3(1,1,1);
        // std::cout<<"\t point  1x "<< (points[(*it)[0]])[0] << std::endl;
        // std::cout<<"\t point  1y "<< (points[(*it)[0]])[1] << std::endl;
        // std::cout<<"\t point  1z "<< (points[(*it)[0]])[2] << std::endl;
        Point3  p1=Point3((fl_points_[(*it)[0]])[0],(fl_points_[(*it)[0]])[1],(fl_points_[(*it)[0]])[2]);
        Point3  p2=Point3((fl_points_[(*it)[1]])[0],(fl_points_[(*it)[1]])[1],(fl_points_[(*it)[1]])[2]);
        Point3  p3=Point3((fl_points_[(*it)[2]])[0],(fl_points_[(*it)[2]])[1],(fl_points_[(*it)[2]])[2]);
        Vector3             n  = CGAL::cross_product(p2-p1, p3-p1);
        // std::cout << "element " << *it << std::endl;
        // std::cout << "normal " << (n/ std::sqrt(n*n))[0] << std::endl;
        (faces_[iface])[0]=(*it)[0];
        if( (n/ std::sqrt(n*n))[2]>0){
            (faces_[iface])[1]=(*it)[1];
            (faces_[iface])[2]=(*it)[2];
        }
        else{
            (faces_[iface])[1]=(*it)[2];
            (faces_[iface])[2]=(*it)[1];}
        // std::cout<<"\t point  1 "<< points[(*it)[0]] << std::endl;        
        // std::cout<<"\t point  2 "<< points[(*it)[1]] << std::endl;        
        // std::cout<<"\t point  3 "<< points[(*it)[2]] << std::endl;        
        iface++;
        }
        // second run for security ad for setting the normals
        for(int iiface=0;iiface<faces_.size();iiface++ ){
        Point3  p1=Point3((fl_points_[((faces_[iiface])[0])])[0],(fl_points_[((faces_[iiface])[0])])[1],(fl_points_[((faces_[iiface])[0])])[2]);
        Point3  p2=Point3((fl_points_[((faces_[iiface])[1])])[0],(fl_points_[((faces_[iiface])[1])])[1],(fl_points_[((faces_[iiface])[1])])[2]);
        Point3  p3=Point3((fl_points_[((faces_[iiface])[2])])[0],(fl_points_[((faces_[iiface])[2])])[1],(fl_points_[((faces_[iiface])[2])])[2]);
        Vector3             n  = CGAL::cross_product(p2-p1, p3-p1);
        for(int idim =0;idim<3;idim++)(normals_[iiface])[idim]=(float)(n/ std::sqrt(n*n))[idim];
        }
        
        
    print_interface_oriented();
}
