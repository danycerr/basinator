#include<domain.hpp>
// ///////////////////////////////
domain::domain(){
//     create_domain();
    create_multiple_domain();
}
///////////////////////////////
void domain::create_domain(){
    std::cout<<"domain::Creating domain"<<std::endl;
    
    Polyhedron p_domain;   
    std::ifstream input_xmin("./xmin.off");input_xmin>>p_domain;
    std::ifstream input_xmax("./xmax.off");input_xmax>>p_domain;
    std::ifstream input_ymin("./ymin.off");input_ymin>>p_domain;
    std::ifstream input_ymax("./ymax.off");input_ymax>>p_domain;
    
    std::ifstream input_top("./horizon_1.off");
    input_top>>p_domain;
    std::ifstream input_bot("./horizon_2.off");
    input_bot>>p_domain;
    
    
    // Create domain
    Mesh_domain domain(p_domain);
    
    // Get sharp features
    domain.detect_features();

                        

    // Mesh generation
    std::cout<<"Geneartion of poligonal costrained mesh"<<std::endl;
    Timer t;
    t.start();
    // defining meshing criteria
//     Surface_mesh s_mesh;

    std::ofstream domain_out ("domain.off");
    domain_out << p_domain;
    domain_out.close(); 
    std::ifstream domain_pt("domain.off");
    
//     domain_pt>>s_mesh ;
//     CGAL::Polygon_mesh_processing::triangulate_faces(s_mesh);
//     std::ofstream cube_off("../s_3dbasin.off");
//     cube_off << s_mesh;
//     std::cout<<"crete domain:: exporting basin surface"<<std::endl;
    // exporter exp("../s_3dbasin.off");
    
    
    
//       // surface mesh generation
//       Trs trs;            // 3D-Delaunay triangulation
//       C2t3 c2t3 (trs);   // 2D-complex in 3D-Delaunay triangulation
//     
//     
//       // defining meshing criteria
//       CGAL::Surface_mesh_default_criteria_3<Trs> criteria_surface(30.,  // angular bound
//                                                      0.1,  // radius bound
//                                                      0.1); // distance bound
//     
//       // meshing surface
//       CGAL::make_surface_mesh(c2t3, domain, criteria_surface, CGAL::Non_manifold_tag());
//       std::cout << "Final number of points: " << trs.number_of_vertices() << "\n";
    
    
    
    std::cerr << "Generation of 3d mesh costrained in "<< t.time() << " sec." << std::endl;
    
// 	Mesh_criteria criteria(
// 			edge_size = 10,
// 			facet_angle = 25
// 			, facet_size = 40
// 			, facet_distance = 0.05
// 			,cell_radius_edge_ratio = 3
//                         , cell_size = 2
// 			);
// // //   Mesh criteria
//   Mesh_criteria criteria (
//       edge_size=1 //0.1
//       ,facet_angle=30
//       ,facet_size=20
//       // ,facet_topology=0
//       //,facet_distance=4
//       //,cell_radius_edge_ratio=3
//      // cell_size=1
//                          );

    // Criteria
Edge_criteria edge_criteria(0);
Facet_criteria facet_criteria(30, 100, 10); // angle, size, approximation
Cell_criteria cell_criteria(2, 100); // radius-edge ratio, size
Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);



    C3t3 c3t3; // tridimensional mesh with costrained surfaces
    c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
            no_perturb(), no_exude());
    // Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
    std::cerr << "End of Generation of 3d mesh "<< t.time() << " sec." << std::endl;
    // CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);
    std::cerr << "End of Perturbation "<< t.time() << " sec." << std::endl;
    // Exudation
    // CGAL::exude_mesh_3(c3t3,12);
    std::cerr << "End of exudation "<< t.time() << " sec." << std::endl;
    std::cerr << "Generation of 3d mesh costrained in "<< t.time() << " sec." << std::endl;
    // Output
    dump_c3t3(c3t3, "../3dbasin_dom");
    // Polyhedron N_top(top);
    //   //Polyhedron_ei bot;
    //   //std::ifstream input_bot("./horizon_3.off");
    //   //input_bot>>bot;
}

void domain::create_multiple_domain(){
    std::cout<<"domain::Creating domain"<<std::endl;
    bool lateral_polylines=true;
    Polyhedron p_domain;   
    Polyhedron pint_domain;   
    if(1){
    int label=10;    
    std::ifstream input_xmin("./xmin." + std::to_string(label) + ".off");input_xmin>>p_domain;
    std::ifstream input_xmax("./xmax." + std::to_string(label) + ".off");input_xmax>>p_domain;
    std::ifstream input_ymin("./ymin." + std::to_string(label) + ".off");input_ymin>>p_domain;
    std::ifstream input_ymax("./ymax." + std::to_string(label) + ".off");input_ymax>>p_domain;
    }   
    if(1){
    int label=20;    
    std::ifstream input_xmin("./xmin." + std::to_string(label) + ".off");input_xmin>>p_domain;
    std::ifstream input_xmax("./xmax." + std::to_string(label) + ".off");input_xmax>>p_domain;
    std::ifstream input_ymin("./ymin." + std::to_string(label) + ".off");input_ymin>>p_domain;
    std::ifstream input_ymax("./ymax." + std::to_string(label) + ".off");input_ymax>>p_domain;
    }
    std::ifstream input_top("./horizon_1.off");
    input_top>>p_domain;
	std::ifstream input_bot("./horizon_3.off");
    input_bot>>p_domain;
     std::ifstream input_mid("./horizon_2.off");
     input_mid>>pint_domain;
    Polylines polylines;
    if(lateral_polylines){
    //Point3* p_min= new Point3(1.e+9,1.e+9,1.e+9);
    //b_box={x_min,x_max,x_min,x_max,z_min,z_max}
    std::vector<double> b_box; b_box.resize(6);
    bouding_box(0,0,p_domain,b_box);
    // Generation of all lateral polylines
    //Polylines xmin/max,ymin/max,zmin--xmin/max,ymin/max,zmax
    { 
        for(int iy=0;iy<2;iy++)
            for(int ix=0;ix<2;ix++){
                std::vector<Point3> b_points;
                //point xmin/max,ymin/max,zmin &zmax
                for(int iz=0;iz<2;iz++) 
                    b_points.push_back(Point3(b_box[ix+0*2],
                                            b_box[iy+1*2],
                                            b_box[iz+2*2]));
                //point xmin,ymin,zmax
                std::vector<Point3> rb_points;rb_points.resize(b_points.size());
                std::vector<Point3*> sb_points_buf;sb_points_buf.resize(b_points.size());
                std::vector<Point3> sb_points;
                //Search for the real points near the bounding points
                find_nearest_point(p_domain,b_points,rb_points);
//                 std::cout<<"Finded    "<<b_points[0]<<" <-> "<<rb_points[0]<<std::endl;
//                 std::cout<<"Finded    "<<b_points[1]<<" <-> "<<rb_points[1]<<std::endl;
//                 std::cout<<"z rbpoint    1 "<<(rb_points[0])[2]<<std::endl;
//                 std::cout<<"z rbpoint    2 "<<(rb_points[1])[2]<<std::endl;
//                 sb_points[0]=((rb_points[0])[2],(rb_points[0])[2],(rb_points[0])[2]);
                // scale a bit the TWO point in vertical direction
                double delta=0.e-0+0;
                double zc = delta*fabs((rb_points[1])[2]- (rb_points[0])[2]);
                double dir=-1.;
                for(int iz=0;iz<2;iz++){
                    dir*=-1;
                    sb_points_buf[iz]=new Point3 ((rb_points[iz])[0],
                                            (rb_points[iz])[1],
                                            ((rb_points[iz])[2]+dir*zc));
                
                    sb_points.push_back(*sb_points_buf[iz]);
//                     sb_points.push_back(rb_points[iz]);
                }
                //Push the polyline in the structure polylines
                create_polyline(sb_points,polylines);
            }
        
    }
    }
    std::cout<<"End of generation of lateral polylines"<<std::endl;
    insert_features_of_horizon("./horizon_2.off",polylines);
    insert_features_of_horizon("./horizon_1.off",polylines);
    insert_features_of_horizon("./horizon_3.off",polylines);
    
    
    
//     // Create domain
// 	Mesh_domain domain(p_domain);
    // Create a polyhedral domain, with only one polyhedron,
    // and no "bounding polyhedron", so the volumetric part of the domain will be
    // empty.
    std::vector<Polyhedron*> poly_ptrs_vector(1, &p_domain);
//     Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
    
	Mesh_domain domain(pint_domain,p_domain);
    // Get sharp features
    domain.add_features(polylines.begin(),polylines.end());    
    domain.detect_features();                 

    // Mesh generation
    std::cout<<"Geneartion of poligonal costrained mesh"<<std::endl;
    Timer t;
    t.start();

    
    std::cerr << "Generation of 3d mesh costrained in "<< t.time() << " sec." << std::endl;

    // Criteria
Edge_criteria edge_criteria(100);
Facet_criteria facet_criteria(30, 100, 1.); // angle, size, approximation
Cell_criteria cell_criteria(3, 100); // radius-edge ratio, size
Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

C3t3 c3t3; // tridimensional mesh with costrained surfaces
c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
            no_perturb(), no_exude());
    // Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
    std::cerr << "End of Generation of 3d mesh "<< t.time() << " sec." << std::endl;
    CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);
    std::cerr << "End of Perturbation "<< t.time() << " sec." << std::endl;
    // Exudation
    CGAL::exude_mesh_3(c3t3,12);
    std::cerr << "End of exudation "<< t.time() << " sec." << std::endl;
    std::cerr << "Generation of 3d mesh costrained in "<< t.time() << " sec." << std::endl;
// 	dump_c3t3(c3t3, "../3dbasin_dom");
    std::string name="./patch_6";
    dump_c3t3(c3t3, name.c_str());
//     std::ofstream domain_out ((name + ".off").c_str());
//     domain_out << c3t3;
//     domain_out.close();
    dump_off(c3t3,name.c_str(),false); //  planar not planar
    exporter exp((name + ".off").c_str());
    
    // Get sharp features
    // omain.add_features(polylines.begin(),polylines.end());
//     domain_2.detect_features();                     

    // Mesh generation o
//     std::cout<<"Geneartion of poligonal costrained mesh"<<std::endl; 
//     C3t3 c3t3_v2; // tridimensional mesh with costrained surfaces
//     c3t3_v2 = CGAL::make_mesh_3<C3t3>(domain_2, criteria,
//             no_perturb(), no_exude());
//     // Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
//     std::cerr << "End of Generation of 3d mesh "<< t.time() << " sec." << std::endl;
//     std::cerr << "Generation of 3d mesh costrained in "<< t.time() << " sec." << std::endl;
//     dump_c3t3(c3t3, "../3dbasin_dom");
}   

/////////////////////////////////////////////////7
void domain::bouding_box( 
                        int dir,
                        int minmax,
                        Polyhedron &poly,
                        std::vector<double>& b_box
                        //Point3* point
                        )
{   
    Point3* point=new Point3(0,0,0);
    std::vector<double> pt_min={1.e+10, 1.e+10, 1.e+10};   
    std::vector<double> pt_max={-1.e+10, -1.e+10, -1.e+10};
    
    for ( Vertex_iterator v = poly.vertices_begin(); 
                v != poly.vertices_end(); ++v){
        for (int idim=0;idim<3;idim++){
            pt_min[idim]=((v->point()[idim]<pt_min[idim])? v->point()[idim]:pt_min[idim]);
            pt_max[idim]=((v->point()[idim]>pt_max[idim])? v->point()[idim]:pt_max[idim]); 
        }
    //*point=v->point();
    }
    
    std::cout<<"Bounding points are "<<std::endl;
    std::cout<<"Min ";
    for(int i=0; i<pt_min.size(); i++) 
    {
        std::cout<<pt_min[i]<<" ";
        b_box[2*i]=pt_min[i];
        b_box[2*i+1]=pt_max[i];
    }
    std::cout<<std::endl;
    std::cout<<"Max ";
    for(int i=0; i<pt_max.size();i++) std::cout<<pt_max[i]<<" ";
    std::cout<<std::endl;
}
/////////////////////////////////////////////////////////
//Method for finding the nerast real point
////////////////////////////////////////////////////////////
    void domain::find_nearest_point(
        Polyhedron &poly,             ///Polyedron to search
        std::vector<Point3>& points,  /// searcing points
        std::vector<Point3>& d_points /// Real finded points
    )
    {
        for (int i_point=0;i_point<points.size();i_point++){
            double dist=1.e+8;
            for ( Vertex_iterator v = poly.vertices_begin(); 
                v != poly.vertices_end(); ++v){
                if(dist>CGAL::sqrt( CGAL::to_double(
                CGAL::squared_distance(points[i_point],v->point())) ))
                    {
                        d_points[i_point]=v->point();
                        dist = CGAL::sqrt( CGAL::to_double(
                            CGAL::squared_distance(
                                points[i_point],v->point())) );
                        
                    }
                    
                }     
            dist = CGAL::sqrt( CGAL::to_double(
            CGAL::squared_distance(points[i_point],d_points[i_point])) );
            std::cout<<i_point<<": Distance  "<<dist <<
            " of point "<< d_points[i_point] <<std::endl;
            
        }
    }
//////////////////////////////////////////////
//Method for polyline generation
//////////////////////////////////////////////
void domain::create_polyline(       
        std::vector<Point3>& points,
        Polylines& polylines
                            )
{
    Polyline_3 polyline;
    for(int i_point=0;i_point<points.size(); i_point++)
        polyline.push_back(points[i_point]);
    polylines.push_back(polyline);
    std::cout<<"End of Polylines creation"<<std::endl;
    
}
//////////////////////////////////////////////
//Method for inserting features of horizons
//////////////////////////////////////////////
void domain::insert_features_of_horizon(
        const char* off_filename,
        Polylines& polylines,
        int method
    )
{
    std::cout<<"Inserting features of "<<off_filename<<std::endl;
    std::ifstream input_file(off_filename);
    Polyhedron horizon; 
    input_file>>horizon;
    std::vector<double> b_box; b_box.resize(6);
    bouding_box(0,0,horizon,b_box);
    if (method==0)
    {

    for (int i_pos=0;i_pos<4;i_pos++)//xmin max ymin ymax
    {
        std::vector<Point3> rb_points;
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); ++v)
            if(fabs(CGAL::to_double(v->point()[(int)(i_pos/2)])-b_box[i_pos])<50)                
                        rb_points.push_back(v->point());
    //     for ( Vertex_iterator v = horizon.vertices_begin();
    //          v != horizon.vertices_end(); ++v)
    //          if(fabs(CGAL::to_double(v->point()[1])-b_box[3])<50)    //ymax            
    //                     rb_points.push_back(v->point());
        //for (int i=0; i<rb_points.size();i++) std::cout<<rb_points[i]<<std::endl;
            //erase first and last component already used for vertical
        rb_points.erase(rb_points.begin()+rb_points.size()-1);
        rb_points.erase(rb_points.begin());
        create_polyline(rb_points,polylines);
    }
    }
    if (method==1) //complete horizon
    {
    	// Declaring the type of Predicate that accepts 2 pairs and return a bool
	typedef std::function<bool(std::pair<int, double>, std::pair<int, double>)> Comparator;
    	// Defining a lambda function to compare two pairs. 
//     It will compare two pairs using second field
	Comparator compFunctor =
			[](std::pair<int, double> elem1 ,std::pair<int, double> elem2)
			{
				return elem1.second < elem2.second;
			};
        int norm_dir; double delta=1.e+9;
    for (int idir=0; idir<3;idir++)
        if(fabs(b_box[2*idir] - b_box[2*idir+1])<delta){
            norm_dir=idir;
            delta=fabs(b_box[2*idir] - b_box[2*idir+1]);
        }
//         Searching center of the bounding box
        double pt_c[3]={0,0,0};
        for (int idir=0; idir<3;idir++)
        pt_c[idir]=(b_box[2*idir+1] + b_box[2*idir])/2;
        Point3 box_center(pt_c[0],pt_c[1],pt_c[2]);
    std::cout<<"Normal of "<<off_filename<<" is  "<<norm_dir <<std::endl;
    std::cout<<"And center of the bounding box is "
            <<pt_c[0]<<" "<<pt_c[1]<<" "<<pt_c[2]<<std::endl;
//     std::cout<<"and in direction "<<norm_dir<< " we use "<< minmax<<std::endl;
    int dir1=(norm_dir+1)%3;
    int dir2=(norm_dir+2)%3;
    std::cout<<" dir 1 "<< dir1<<" dir 2 "<< dir2<<std::endl;
//     for ( Vertex_iterator v = horizon.vertices_begin();
//             v != horizon.vertices_end(); ++v)
//         std::cout<<v->point()<<std::endl;
    double eps=150;
    std::vector<Point3> rrb_points;
    std::vector<Point3> srb_points;
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); ++v)
            if(fabs(CGAL::to_double(v->point()[dir1])-b_box[2*dir1])<eps)                
            { rrb_points.push_back(v->point());
//                 std::cout<<"adding 1"<<std::endl;
            }
            
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); ++v)
            if(fabs(CGAL::to_double(v->point()[dir2])-b_box[2*dir2+1])<eps)    //ymax            
            { rrb_points.push_back(v->point());
//                 std::cout<<"adding 2"<<std::endl;
            }
            
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); v++)
            if(fabs(CGAL::to_double(v->point()[dir1])-b_box[2*dir1+1])<eps)                
            { rrb_points.push_back(v->point());
//                 std::cout<<"adding 3"<<std::endl;
            }
            
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); v++)
            if(fabs(CGAL::to_double(v->point()[dir2])-b_box[2*dir2])<eps)    //ymax            
            { rrb_points.push_back(v->point());
//                 std::cout<<"adding 4 "<< v->point() <<std::endl;
            }
    
        double cum_dist=0;//cumulate distance
        double radius=0;
        std::map <int, double> ponint_angle;
        for(int i=0; i<rrb_points.size();i++) {
            double angle=  asin( fabs( (rrb_points[i])[1] - box_center[1] ) / 
                        CGAL::sqrt( CGAL::to_double(
                CGAL::squared_distance(rrb_points[i],box_center)) ) );
            if((rrb_points[i])[1] - box_center[1]>0 && 
                (rrb_points[i])[0] - box_center[0]>0 ){
                    // first quadrant nothing to do
                }
                else if((rrb_points[i])[1] - box_center[1]>0 && 
                (rrb_points[i])[0] - box_center[0]<0 ){
                    // second quadrant
                    angle=3.14-angle;
                }
                else if((rrb_points[i])[1] - box_center[1]<0 && 
                (rrb_points[i])[0] - box_center[0]<0 ){
                    // second quadrant
                    angle=3.14+angle;
                }
                else if((rrb_points[i])[1] - box_center[1]<0 && 
                (rrb_points[i])[0] - box_center[0]>0 ){
                    // second quadrant
                    angle=2*3.14-angle;
                }
                ponint_angle.insert(std::pair <int, double> (i, angle));
//             std::cout<<rrb_points[i]<< " angle " << angle << " and "<<
//                                                     atan2(
//                                                         (rrb_points[i])[0] - box_center[0],
//                                                         (rrb_points[i])[1] - box_center[1])
//                                                      <<std::endl;
        }
            // printing map gquiz1
    // Declaring a set that will store the pairs using above comparision logic
	std::set<std::pair<int, double>, Comparator> pnt_ang_sorted(
			ponint_angle.begin(), ponint_angle.end(), compFunctor);
//     std::map <int, double> :: iterator itr;
//     for (itr = ponint_angle.begin(); itr != ponint_angle.end(); ++itr)
//     {
//         std::cout  <<  '\t' << itr->first 
//               <<  '\t' << itr->second << '\n';
//     }
    	for (std::pair<int , double> element : pnt_ang_sorted){
		std::cout << element.first << " :: " << element.second <<" Point "<<rrb_points[element.first] << std::endl;
        srb_points.push_back(rrb_points[element.first]);
        }
        srb_points.push_back(srb_points[0]);
        create_polyline(srb_points,polylines);
    }
    std::cout<<"End of Inserting features of "<<off_filename<<std::endl;
}

////////////////////////////////////////////////////////
// function for writing the reconstruction output in the off format
void domain::dump_off(const C3t3& c3t3, std::string name,bool planar_face)
{
  typedef typename C3t3::Triangulation Triangulation;
  typedef typename Triangulation::Vertex_handle Vertex_handle;

  std::map<Vertex_handle, std::size_t> V;
  
  std::size_t inum = 0; 
  std::size_t nfacets = 0;
  std::array<std::size_t,3> indices={{0,0,0}};
  std::stringstream facet_buffer,vertex_buffer;
  std::vector<Point3> points_buffer;
  std::vector<int> faces_buffer;
  int node_per_facet=3;
  for(typename C3t3::Facets_in_complex_iterator 
        fit = c3t3.facets_in_complex_begin(),
        end = c3t3.facets_in_complex_end();
      fit != end; ++fit) 
  {
    typename C3t3::Subdomain_index cell_sd=c3t3.subdomain_index(fit->first);
    typename C3t3::Subdomain_index opp_sd=c3t3.subdomain_index(fit->first->neighbor(fit->second));
//     
//     if (cell_sd!=sd_index && opp_sd!=sd_index) continue;
// 
    ++nfacets;
    int j=-1;
//     
//     
    for (int i = 0; i < 4; ++i)
      if (i != fit->second)
          indices[++j]=mesh_3_export::get_vertex_index((*fit).first->vertex(i), V, inum,vertex_buffer,points_buffer);
//     if ( ( (cell_sd==sd_index) == (fit->second%2 == 1) ) == normals_point_outside_of_the_subdomain )
//       std::swap(indices[0],indices[1]);
//     std::cout << "3" << " " << indices[0] <<" " << indices[1] <<" " << indices[2] << "\n";
      for (int i_idx = 0; i_idx < node_per_facet; ++i_idx) faces_buffer.push_back(indices[i_idx]);
  }
  
//   std::cout<<"Pointssss "<<std::endl;
//   for (int i=0; i<points_buffer.size(); i++)std::cout<<points_buffer[i]<<std::endl;
  
  std::cout<<"Num of facet "<< nfacets<<std::endl;
 
  // Reorient faces
    std::vector<double> b_box; b_box.resize(6);
    bouding_box_of_point_cloud(points_buffer,b_box);
    if(planar_face){
      int min=-1;double delta=1.e+20;
      for(int idir=0;idir<3;idir++)
	if(fabs(b_box[idir*2+1]-b_box[idir*2]) < delta){
	  min=idir;
	  delta=fabs(b_box[idir*2+1]-b_box[idir*2]);
	}
      b_box[2*min]=0;b_box[2*min+1]=0;
      std::cout<<"Found planar face in "<<min<<std::endl;
    }
  //         Searching center of the bounding box
    double pt_c[3]={0,0,0};
    for (int idir=0; idir<3;idir++)
        pt_c[idir]=(b_box[2*idir+1] + b_box[2*idir])/2;
    Point3 box_center(pt_c[0],pt_c[1],pt_c[2]);
    for (int i_face=0; i_face<faces_buffer.size()/3;i_face++){
      Point3  p1=Point3( points_buffer[faces_buffer[3*i_face+0]][0],
			 points_buffer[faces_buffer[3*i_face+0]][1],
			 points_buffer[faces_buffer[3*i_face+0]][2]
 		      );
      Point3  p2=Point3( points_buffer[faces_buffer[3*i_face+1]][0],
			 points_buffer[faces_buffer[3*i_face+1]][1],
			 points_buffer[faces_buffer[3*i_face+1]][2]
 		      );
      Point3  p3=Point3( points_buffer[faces_buffer[3*i_face+2]][0],
			 points_buffer[faces_buffer[3*i_face+2]][1],
			 points_buffer[faces_buffer[3*i_face+2]][2]
 		      );
      Vector3             n  = CGAL::cross_product(p2-p1, p3-p1);
      //Center of element
      double pt_ce[3]={0,0,0};
      for (int idir=0; idir<3;idir++)
	for (int idx=0; idx<3;idx++){
	  double check=points_buffer[faces_buffer[3*i_face+idx]][idir];
	  pt_ce[idir]+=points_buffer[faces_buffer[3*i_face+idx]][idir]/3.;
	}
      Point3 cent_el(pt_ce[0],pt_ce[1],pt_ce[2]);
      Point3 cent_el_translated(pt_ce[0],pt_ce[1],pt_ce[2]);
      Transformation scale(CGAL::SCALING,1.e+0/ std::sqrt(n*n));
      n=scale.transform(n);
      Transformation translate(CGAL::TRANSLATION, n);
      cent_el_translated = translate.transform(cent_el);
      double proj=CGAL::scalar_product(n, cent_el - box_center);
      double dist_tra = CGAL::sqrt(
	                CGAL::to_double(
		        CGAL::squared_distance(cent_el_translated,box_center)
		    	   ));
      double dist_nat = CGAL::sqrt(
	                 CGAL::to_double(
		         CGAL::squared_distance(cent_el,box_center)
			   ));
      double dist_tx=(cent_el_translated[0]-box_center[0])*(cent_el_translated[0]-box_center[0]);
      double dist_ty=(cent_el_translated[1]-box_center[1])*(cent_el_translated[1]-box_center[1]);
      double dist_tz=(cent_el_translated[2]-box_center[2])*(cent_el_translated[2]-box_center[2]);
      double dist_nx=(cent_el[0]-box_center[0])*(cent_el[0]-box_center[0]);
      double dist_ny=(cent_el[1]-box_center[1])*(cent_el[1]-box_center[1]);
      double dist_nz=(cent_el[2]-box_center[2])*(cent_el[2]-box_center[2]);
//       if(dist_nat>dist_tra){
//       if(  (dist_tx-dist_nx)*fabs(n[0])
//           +(dist_tx-dist_nx)*fabs(n[1])
//           +(dist_tx-dist_nx)*fabs(n[2])<0 )
      if(proj<0){
// 	std::cout<<scale.transform(n) <<std::endl;
	std::swap(faces_buffer[3*i_face+0],
	          faces_buffer[3*i_face+1]);
// 	std::cout<<"swapped "<< 3*i_face<<std::endl;
      }
//       else std::cout<<"Not Printing off!!! Swapping"<<std::endl;
	
    }// end element face
    std::cout<<"Dump off :: End face reorientation"<<std::endl;
std::ofstream output((name+".off").c_str());
//  c3t3.output_boundary_to_off(output);
// output.close();
  output << "OFF \n" << points_buffer.size() << " "
         << nfacets << " 0\n";

  output << vertex_buffer.str();
  for( int iface=0; iface<faces_buffer.size()/3;iface++){
      output << "3 ";
      for( int ipt=0; ipt<3;ipt++)
          output << faces_buffer[3*iface+ipt]<<" ";
      output << "\n";
  }
  output << "\n"; 
  output.close();
}
//////////////////////////////////////////////
//////////////////////////////////////////7
    void domain::bouding_box_of_point_cloud(
        std::vector<Point3>& pt_cloud,
        std::vector<double>& b_box
    )
    {
    std::vector<double> pt_min={1.e+10, 1.e+10, 1.e+10};   
    std::vector<double> pt_max={-1.e+10, -1.e+10, -1.e+10};
    
    for ( int v=0; v< pt_cloud.size(); v++){
        for (int idim=0;idim<3;idim++){
            pt_min[idim]=(((pt_cloud[v])[idim]<pt_min[idim])? (pt_cloud[v])[idim]:pt_min[idim]);
            pt_max[idim]=(((pt_cloud[v])[idim]>pt_max[idim])? (pt_cloud[v])[idim]:pt_max[idim]); 
        }
    //*point=v->point();
    }
    std::cout<<"Bounding points of points clouds are "<<std::endl;
    std::cout<<"Min ";
    for(int i=0; i<pt_min.size(); i++) 
    {
        std::cout<<pt_min[i]<<" ";
        b_box[2*i]=pt_min[i];
        b_box[2*i+1]=pt_max[i];
    }
    std::cout<<std::endl;
    std::cout<<"Max ";
    for(int i=0; i<pt_max.size();i++) std::cout<<pt_max[i]<<" ";
    std::cout<<std::endl;
    }
    
    
////////////////////////////////////////////////////////
void domain::create_single_domain(){
      std::cout<<"Second round "<<std::endl;
    Polyhedron internal_surf, p2_domain;
    std::cout<<"domain::Creating domain from shell"<<std::endl;
//     std::ifstream input_top_p("./patch_1.off");
//     input_top_p>>p2_domain;
//     std::ifstream input_bot_p("./patch_3.off");
//     input_bot_p>>p2_domain;    
//     std::ifstream input_mid_p("./patch_2.off");
//     input_mid_p>>internal_surf;
        std::ifstream input_mid_p("./patch_4.off");
    input_mid_p>>internal_surf;
	// Create domain
//     Mesh_domain domain(internal_surf,p2_domain);
    Mesh_domain domain(p2_domain);
    std::cout<<"End second reading"<<std::endl;
//     Polylines polylines;
//     insert_features_of_horizon("./patch_2.off",polylines);
//     domain.add_features(polylines.begin(),polylines.end());    
        domain.detect_features();
        // Criteria
    Edge_criteria edge_criteria(100);
    Facet_criteria facet_criteria(30, 100, 10); // angle, size, approximation
    Cell_criteria cell_criteria(2, 100); // radius-edge ratio, size
    Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

    C3t3 c3t3; // tridimensional mesh with costrained surfaces
    c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                no_perturb(), no_exude());
    std::cerr << "Generation of 3d mesh costrained "<< std::endl;
	dump_c3t3(c3t3, "../3dbasin_single_dom");
}

void domain::create_horizon(const char * file){     
    Polyhedron p_domain;   
    std::ifstream input(file);
    input>>p_domain;
    std::vector<Polyhedron*> poly_ptrs_vector(1, &p_domain);
    Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
// 	Mesh_domain domain(p_domain);
    Polylines polylines;
    insert_features_of_horizon(file,polylines, 1); // 1 for complete horizon
    // Get sharp features
   domain.add_features(polylines.begin(),polylines.end());    
//     domain.detect_features();           
    // Mesh generation
    std::cout<<"Geneartion of poligonal costrained mesh"<<std::endl;
    // Criteria
    Edge_criteria edge_criteria(100);
    Facet_criteria facet_criteria(30, 100, 1.); // angle, size, approximation
    Cell_criteria cell_criteria(3, 100); // radius-edge ratio, size
    Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

    C3t3 c3t3; // tridimensional mesh with costrained surfaces
    c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
            no_perturb(), no_exude());
    // Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
    std::cerr << "End of Generation of 3d mesh " << std::endl;
//     std::string name=;
    std::string hname=file;
    hname.erase(0,hname.find_last_of("/")+1);
    hname.erase(hname.find_last_of("."));
    dump_c3t3(c3t3, ("mesh_" + hname).c_str());
    dump_off(c3t3,("mesh_" + hname).c_str(),false); //  planar not planar
    exporter exp(("mesh_" + hname + ".off").c_str());
    store_boundary(("mesh_" + hname + ".off").c_str());
}
////////////////////////////////
///Method for store bounding of horizon
////////////////7
void domain::store_boundary(
        const char* off_filename){
    std::cout<<"Storing features of "<<off_filename<<std::endl;
    std::ifstream input_file(off_filename);
    Polyhedron horizon; 
    input_file>>horizon;
    std::vector<double> b_box; b_box.resize(6);
    bouding_box(0,0,horizon,b_box);
   	// Declaring the type of Predicate that accepts 2 pairs and return a bool
	typedef std::function<bool(std::pair<int, double>, std::pair<int, double>)> Comparator;
    // Defining a lambda function to compare two pairs. 
    // It will compare two pairs using second field
	Comparator compFunctor =
			[](std::pair<int, double> elem1 ,std::pair<int, double> elem2)
			{
				return elem1.second < elem2.second;
			};
    int norm_dir; double delta=1.e+9;
    for (int idir=0; idir<3;idir++)
        if(fabs(b_box[2*idir] - b_box[2*idir+1])<delta){
            norm_dir=idir;
            delta=fabs(b_box[2*idir] - b_box[2*idir+1]);
        }
//         Searching center of the bounding box
        double pt_c[3]={0,0,0};
        for (int idir=0; idir<3;idir++)
        pt_c[idir]=(b_box[2*idir+1] + b_box[2*idir])/2;
        Point3 box_center(pt_c[0],pt_c[1],pt_c[2]);
    std::cout<<"Normal of "<<off_filename<<" is  "<<norm_dir <<std::endl;
    std::cout<<"And center of the bounding box is "
            <<pt_c[0]<<" "<<pt_c[1]<<" "<<pt_c[2]<<std::endl;
//     std::cout<<"and in direction "<<norm_dir<< " we use "<< minmax<<std::endl;
    int dir1=(norm_dir+1)%3;
    int dir2=(norm_dir+2)%3;
    std::cout<<" dir 1 "<< dir1<<" dir 2 "<< dir2<<std::endl;
//     for ( Vertex_iterator v = horizon.vertices_begin();
//             v != horizon.vertices_end(); ++v)
//         std::cout<<v->point()<<std::endl;
    double eps=5;
    std::vector<Point3> rrb_points;
    std::vector<Point3> srb_points;
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); ++v)
            if(fabs(CGAL::to_double(v->point()[dir1])-b_box[2*dir1])<eps)                
            { rrb_points.push_back(v->point());
//                 std::cout<<"adding 1"<<std::endl;
            }
            
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); ++v)
            if(fabs(CGAL::to_double(v->point()[dir2])-b_box[2*dir2+1])<eps)    //ymax            
            { rrb_points.push_back(v->point());
//                 std::cout<<"adding 2"<<std::endl;
            }
            
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); v++)
            if(fabs(CGAL::to_double(v->point()[dir1])-b_box[2*dir1+1])<eps)                
            { rrb_points.push_back(v->point());
//                 std::cout<<"adding 3"<<std::endl;
            }
            
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); v++)
            if(fabs(CGAL::to_double(v->point()[dir2])-b_box[2*dir2])<eps)    //ymax            
            { rrb_points.push_back(v->point());
//                 std::cout<<"adding 4 "<< v->point() <<std::endl;
            }
    
        double cum_dist=0;//cumulate distance
        double radius=0;
        std::map <int, double> ponint_angle;
        for(int i=0; i<rrb_points.size();i++) {
            double angle=  asin( fabs( (rrb_points[i])[1] - box_center[1] ) / 
                        CGAL::sqrt( CGAL::to_double(
                CGAL::squared_distance(rrb_points[i],box_center)) ) );
            if((rrb_points[i])[1] - box_center[1]>0 && 
                (rrb_points[i])[0] - box_center[0]>0 ){
                    // first quadrant nothing to do
                }
                else if((rrb_points[i])[1] - box_center[1]>0 && 
                (rrb_points[i])[0] - box_center[0]<0 ){
                    // second quadrant
                    angle=3.14-angle;
                }
                else if((rrb_points[i])[1] - box_center[1]<0 && 
                (rrb_points[i])[0] - box_center[0]<0 ){
                    // second quadrant
                    angle=3.14+angle;
                }
                else if((rrb_points[i])[1] - box_center[1]<0 && 
                (rrb_points[i])[0] - box_center[0]>0 ){
                    // second quadrant
                    angle=2*3.14-angle;
                }
                ponint_angle.insert(std::pair <int, double> (i, angle));
        }
            // printing map gquiz1
    // Declaring a set that will store the pairs using above comparision logic
	std::set<std::pair<int, double>, Comparator> pnt_ang_sorted(
			ponint_angle.begin(), ponint_angle.end(), compFunctor);
   	for (std::pair<int , double> element : pnt_ang_sorted){
// 	std::cout << element.first << " :: " << element.second <<" Point "<<rrb_points[element.first] << std::endl;
        srb_points.push_back(rrb_points[element.first]);
        }
        srb_points.push_back(srb_points[0]);
        bounding_point_.push_back(srb_points);
//         create_polyline(srb_points,polylines);
    std::cout<<"End of Storing features of "<<off_filename<<std::endl;
    
}
//////////////////////////////////////////////////
void domain::lateral_face_building(std::vector<int>& index)
{
    reconstructurer r;
    for(int i=0;i<index.size();i++) r.add_points(bounding_point_[index[i]]);
    r.build_lateral_surfaces("lateral");
    int n_ref=18; double w=1./(n_ref+1.);
    std::cout<<"Number points in horizons "<<(bounding_point_[index[0]]).size()
        <<" " << (bounding_point_[index[1]]).size()<<std::endl ;
        // pairing points
    //
        std::vector<double> pared_z;
    for (int i_pt=0;i_pt<(bounding_point_[index[0]]).size(); i_pt++)
    {
        double xp=((bounding_point_[index[0]])[i_pt])[0];
        double yp=((bounding_point_[index[0]])[i_pt])[1];
        double dist=1.e+9;
        int min_dist=-1;
        for (int j_pt=0;j_pt<(bounding_point_[index[1]]).size(); j_pt++){
            double xp_j=((bounding_point_[index[1]])[j_pt])[0];
            double yp_j=((bounding_point_[index[1]])[j_pt])[1];
            double dist_buf=(xp-xp_j)*(xp-xp_j) + (yp-yp_j)*(yp-yp_j);
            if (dist_buf<dist){
                min_dist=j_pt;
                dist=dist_buf;
            }
        }
        std::cout<<"Nearest node of "<<i_pt<<" is "<<min_dist<<std::endl;
        pared_z.push_back(((bounding_point_[index[1]])[min_dist])[2]);
    }
    for (int iref=0;iref<n_ref; iref++){
        std::vector<Point3> pt_buf;
        std::vector<double> pt_buf_d={0,0,0};
//         std::cout<<"Refinig adding points at "<<zm << " w "<< w<<std::endl;
        for (int i_pt=0;i_pt<(bounding_point_[index[0]]).size(); i_pt++)
        {
            double zm=  (1-w)*((bounding_point_[index[0]])[i_pt])[2] +w*(pared_z[i_pt]);
            for(int idir=0; idir<2;idir++){
                pt_buf_d[idir]=((bounding_point_[index[0]])[i_pt])[idir];
//                 pt_buf_d[idir]+=(1-w)*((bounding_point_[index[1]])[i_pt])[idir];
        
            }
            pt_buf.push_back(Point3 (pt_buf_d[0],pt_buf_d[1],zm));
//             std::cout<<pt_buf_d[0]<<" "<< pt_buf_d[1]<<" "<<pt_buf_d[2]<<std::endl;
        }
// //             pt_buf.push_back(bounding_point_[index[0]])
         r.add_points(pt_buf);
         w+= 1./(n_ref+1.);
    }
    r.build_lateral_surfaces("lateral");
    std::cout<<"reading lateral surface "<<std::endl; 
    std::ifstream input_file("lateral.off");
    Polyhedron lateral; 
    input_file>>lateral;input_file.close();
    
//     std::cout<<"numb facet"<<lateral.size_of_facets() <<std::endl;
//     typedef typename Polyhedron::Facet_iterator   	       poly_facet_iterator;
//     for(Polyhedron::Vertex_iterator vit; vit!=lateral.vertices_end();
//         ++vit)
//         {
//            std::cout<<"looping verticies"<<std::endl;
//         }
    
    //     for(poly_facet_iterator fit; fit!=lateral.facets_end() ; ++fit) 
//     {
// //         std::cout<<"looping faces"<<std::endl;
//     }
//     std::cout<<lateral<<std::endl;
    
    
    
    
    std::vector<Polyhedron*> poly_ptrs_vector(1, &lateral);
    Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
// 	Mesh_domain domain(p_domain);
    Polylines polylines;
    create_polyline(bounding_point_[index[0]],polylines);
    create_polyline(bounding_point_[index[1]],polylines);
//     insert_features_of_horizon(file,polylines);
    // Get sharp features
   domain.add_features(polylines.begin(),polylines.end());    
    domain.detect_features();           
    // Mesh generation
    std::cout<<"Geneartion of poligonal costrained mesh"<<std::endl;
    // Criteria
    Edge_criteria edge_criteria(100);
    Facet_criteria facet_criteria(30, 100, 1.); // angle, size, approximation
    Cell_criteria cell_criteria(3, 100); // radius-edge ratio, size
    Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

    C3t3 c3t3; // tridimensional mesh with costrained surfaces
    c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
            no_perturb(), no_exude());
    // Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
    std::cerr << "End of Generation of 3d mesh " << std::endl;
//     std::string name=;
    std::string hname="lateral";
    dump_c3t3(c3t3, ("mesh_" + hname).c_str());
    dump_off(c3t3,("mesh_" + hname).c_str(),false); //  planar not planar
    exporter exp((("mesh_" + hname) + ".off").c_str());
}
