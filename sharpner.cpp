#include<sharpner.hpp>

sharpner::sharpner(std::string hname){
    std::cout<<"Welcome to sharpner "<<hname<<std::endl;
    face_off_name_<< path_ << hname << ".off";
    face_sharp_off_name_<<path_ << hname << ".sharp.off";
    face_sharp_mesh_name_<<path_ << hname << ".sharp";
//     std::cout<<"Sharpenig quadrilateral face "
//     << (path_+ hname + ".off").c_str() <<std::endl;
    create_domain();
    exporter exp(face_sharp_off_name_.str().c_str());
}
/////////////////////////////////////
void sharpner::create_domain(){
	std::cout<<"domain::Creating domain form "<< face_off_name_.str().c_str()<<std::endl;
	Polyhedron p_domain;   
    std::ifstream off_input(face_off_name_.str().c_str());
    off_input>>p_domain;
    std::vector<Polyhedron*> poly_ptrs_vector(1, &p_domain);

//     //Point3* p_min= new Point3(1.e+9,1.e+9,1.e+9);
//     //b_box={x_min,x_max,x_min,x_max,z_min,z_max}
//     std::vector<double> b_box; b_box.resize(6);
// 	bouding_box(p_domain,b_box);
//     // Generation of all lateral polylines
    Polylines polylines;
    insert_features_of_quadrilateral_face(
        face_off_name_.str().c_str(), ///off file
        polylines                     /// polyline to protect
    );
//     insert_features_of_horizon("./horizon_2.off",polylines);
//     //insert_features_of_horizon("./horizon_3.off",polylines);
    
    
//     // Create domain
     Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
    
    // Get sharp features
     domain.add_features(polylines.begin(),polylines.end());
     domain.detect_features();                     

	// Mesh generation
	Timer t;
	t.start();
    std::cerr << "Generation of 3d mesh costrained in "<< t.time() << " sec." << std::endl;
  
    // Criteria
    Edge_criteria edge_criteria(10);
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
	dump_c3t3(c3t3, face_sharp_mesh_name_.str().c_str());
    std::ofstream domain_out (face_sharp_off_name_.str().c_str());
    domain_out << c3t3;
    domain_out.close(); 
}   
/////////////////////////////////////////////////7
// Create the bounding box of the domain
///////////////////////////////////////////////////
void sharpner::bouding_box(
                           Polyhedron &poly, /// <- polyedron to build the bounding box
                           std::vector<double>& b_box /// -> vector bounding box points
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
//////////////////////////////////////////////
//Method for inserting features of horizons
//////////////////////////////////////////////
void sharpner::insert_features_of_quadrilateral_face(
        const char* off_filename, /// off file name
        Polylines& polylines      /// polyline to protect
    )
{
    std::cout<<"sharpner:: Inserting features of "<<off_filename<<std::endl;
	std::ifstream input_file(off_filename);
    Polyhedron horizon; 
    input_file>>horizon;
    std::vector<double> b_box; b_box.resize(6);
	bouding_box(horizon,b_box);
    int norm_dir; double delta=1.e+9;
    for (int idir=0; idir<3;idir++)
        if(fabs(b_box[2*idir] - b_box[2*idir+1])<delta){
            norm_dir=idir;
            delta=fabs(b_box[2*idir] - b_box[2*idir+1]);
        }
//     int minmax=(
//         (
//         (horizon.vertices_begin()->point()[norm_dir]-b_box[2*norm_dir])
//         <
//          (horizon.vertices_begin()->point()[norm_dir]-b_box[2*norm_dir+1])
//         )? 0:1);
    std::cout<<"Normal of "<<off_filename<<" is  "<<norm_dir<<std::endl;
//     std::cout<<"and in direction "<<norm_dir<< " we use "<< minmax<<std::endl;
    int dir1=(norm_dir+1)%3;
    int dir2=(norm_dir+2)%3;
    std::cout<<" dir 1 "<< dir1<<" dir 2 "<< dir2<<std::endl;
    double eps=50;
//     for (int i_pos=0;i_pos<4;i_pos++)//xmin max ymin ymax
//     {
        std::vector<Point3> rb_points;
        for ( Vertex_iterator v = horizon.vertices_begin();
            v != horizon.vertices_end(); ++v)
            if(fabs(CGAL::to_double(v->point()[dir1])-b_box[2*dir1])<eps)                
             { rb_points.push_back(v->point());
                 std::cout<<"adding 1"<<std::endl;
            }
            
        for ( Vertex_iterator v = horizon.vertices_begin();
             v != horizon.vertices_end(); ++v)
             if(fabs(CGAL::to_double(v->point()[dir2])-b_box[2*dir2+1])<eps)    //ymax            
             { rb_points.push_back(v->point());
                 std::cout<<"adding 2"<<std::endl;
            }
             
        for ( Vertex_iterator v = horizon.vertices_end();
            v != horizon.vertices_begin(); --v)
            if(fabs(CGAL::to_double(v->point()[dir1])-b_box[2*dir1+1])<eps)                
            { rb_points.push_back(v->point());
                 std::cout<<"adding 3"<<std::endl;
            }
            
        for ( Vertex_iterator v = horizon.vertices_end();
            v != horizon.vertices_begin(); --v)
             if(fabs(CGAL::to_double(v->point()[dir2])-b_box[2*dir2])<eps)    //ymax            
             { rb_points.push_back(v->point());
                 std::cout<<"adding 4"<<std::endl;
            }
            
        for (int i=0; i<rb_points.size();i++) std::cout<<rb_points[i]<<std::endl;
//             //erase first and last component already used for vertical
//         rb_points.erase(rb_points.begin()+rb_points.size()-1);
//         rb_points.erase(rb_points.begin());
        create_polyline(rb_points,polylines);
//     }
    std::cout<<"End of Inserting features of "<<off_filename<<std::endl;
}
//////////////////////////////////////////////
//Method for polyline generation//////////////
//////////////////////////////////////////////
void sharpner::create_polyline(       
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
