#ifndef DOMAINHS
#define DOMAINHS
#include <iostream>
// cgal
// polyedron must be on top for parameters configuration
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
// Triangulation
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Timer.h>

#include <CGAL/aff_transformation_tags.h>
#include <CGAL/IO/read_off_points.h>
// Surface mesh generation
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
// 
#include <exporter.hpp>
#include <reconstructurer.hpp>

#include <vector>
#include <math.h>      
using namespace CGAL::parameters;
namespace mesh_3_export{
template <class Vertex_handle,class Point3>
std::size_t get_vertex_index(Vertex_handle v,std::map<Vertex_handle, std::size_t>& V,std::size_t& inum,
			     std::stringstream& vertex_buffer,std::vector<Point3> &points_buffer){
  std::pair<typename std::map<Vertex_handle, std::size_t>::iterator,bool> res=
    V.insert(std::make_pair(v,inum));
  if (res.second){
    ++inum;
    vertex_buffer <<   res.first->first->point().point() <<"\n"; //point is weighted!
    points_buffer.push_back(res.first->first->point().point());
  }
  return res.first->second;
}
} // end of namespace mesh_3_export

class domain{
typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
//Polyedron
typedef CGAL::Mesh_polyhedron_3<Kernel>::type                               Polyhedron;
typedef Polyhedron::Vertex_iterator                                         Vertex_iterator;
typedef Kernel::Point_3                                                     Point3;
typedef Kernel::Vector_3                                                    Vector3;
// Triangulation
typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel>                Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type                       Tr;
typedef typename Tr::Finite_vertices_iterator                               Finite_vertices_iterator;
typedef typename Tr::Vertex_handle                                          Vertex_handle;
typedef typename Tr::Weighted_point                                         Weighted_point;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index>            C3t3;
  
typedef typename C3t3::Facets_in_complex_iterator                          Facet_iterator;
//Features
typedef std::vector<Point3>                                                Polyline_3;
typedef std::list<Polyline_3>                                              Polylines;
typedef CGAL::Aff_transformation_3<Kernel>  Transformation;

  
  
  
// Criteria
typedef CGAL::Mesh_criteria_3<Tr>     Mesh_criteria;
typedef Mesh_criteria::Edge_criteria  Edge_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria  Cell_criteria;
// default triangulation for Surface_mesher

// typedef CGAL::Surface_mesh<Point3>          Surface_mesh;

// default triangulation for Surface_mesher
// typedef CGAL::Surface_mesh_default_triangulation_3 Trs;
// c2t3
// typedef CGAL::Complex_2_in_triangulation_3<Trs> C2t3;

typedef CGAL::Timer Timer;
public:
    domain();
    //////////////////////////////////////////////
    //Method for creating domain from shell
    //////////////////////////////////////////////
    void create_single_domain();
    
    //////////////////////////////////////////////
    //Method from a reconstruction to a mesh
    //////////////////////////////////////////////
    void create_horizon(const char * file);
    
    //////////////////////////////////////////////////////
    //Method for bulding lateral face between two horizons
    ////////////////////////////////////////////////////
    void lateral_face_building(std::vector<int>& index);
private:
    std::vector<std::vector<Point3>> bounding_point_;
    void create_domain();
    void create_multiple_domain();
    ///////////////////////////////////////
    void bouding_box(
        int dir,
        int minmax,
        Polyhedron &poly,
        std::vector<double>& b_box
    );
    //////////////////////////////////////////////
    ///////////////////////////////////////
    void bouding_box_of_point_cloud(
        std::vector<Point3>& point_cloud,
        std::vector<double>& b_box
    );
    //////////////////////////////////////////////
    //Method for finding nearest point of domain
    //////////////////////////////////////////////
    void find_nearest_point(
        Polyhedron &poly,             ///Polyedron to search
        std::vector<Point3>& points, 
        std::vector<Point3>& d_points);
    //////////////////////////////////////////////
    //Method for polyline generation
    //////////////////////////////////////////////
    void create_polyline(       
        std::vector<Point3>& points, //point that will form polyline
        Polylines& polylines
    );   
    //////////////////////////////////////////////
    //Method for inserting features of horizons
    //////////////////////////////////////////////
    void insert_features_of_horizon(
        const char* off_filename,
        Polylines& polylines,
        int method=0
    );
    //////////////////////////////////////////////
    //Method for storing a boundary
    //////////////////////////////////////////////
    void store_boundary(
        const char* off_filename
    );
    //////////////////////////////////////////////
    //Method for dumping oriented off file
    //////////////////////////////////////////////
    void dump_off(const C3t3& reconstruct, std::string name,bool planar=false);
};
#endif
