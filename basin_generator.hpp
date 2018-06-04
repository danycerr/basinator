#ifndef BASIN_GENERATORHS
#define BASIN_GENERATORHS

#include <iostream>
#include <memory>

// polyedron must be on top for parameters configuration
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
// Triangulation
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
// face reconstruction
#include <fstream>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Timer.h>

#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Scale_space_reconstruction_3/Jet_smoother.h>
#include <CGAL/Scale_space_reconstruction_3/Advancing_front_mesher.h>

#include<sharpner.hpp>

using namespace CGAL::parameters;
class horizon{
typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >                Reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother< Kernel > Smoother;
typedef CGAL::Scale_space_reconstruction_3::Jet_smoother<Kernel>            JetSmoother;
typedef CGAL::Scale_space_reconstruction_3::Advancing_front_mesher<Kernel>  Adv_Mesher;

typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher< Kernel >    Mesher;
typedef Reconstruction::Point                                               Point;
typedef Reconstruction::Facet_const_iterator                                Facet_iterator;
typedef Reconstruction::Point_const_iterator                                Point_iterator;
typedef Mesher::Facet_const_iterator                                        Mesher_iterator;
//Polyedron
typedef CGAL::Mesh_polyhedron_3<Kernel>::type                               Polyhedron;
typedef Kernel::Point_3                                                     Point3;
typedef Kernel::Vector_3                                                    Vector3;
// Triangulation
typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel>                Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type                       Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index>            C3t3;
  
// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                                           Mesh_criteria;
typedef CGAL::Timer Timer;
public:
    std::string hname;
    horizon(char* file);
    void read_horizons_points(char* file);
    void reconstruction_surface();
    Polyhedron poly_,poly_bounding_; //poliedron for the global mesh
    const char* interface_off_name_; // name for the off file of interface
    const char* interface_oriented_off_name_; // name for the off file of interface oriented
    const char* domain_off_name_; // name for the off file of interface
    void create_poligonal_domain(char * file);
    std::vector<Point> get_points(){return rpoints_;}
    void create_domain();
    
private:
    std::unique_ptr<Mesh_domain> domain_costrained_;
    C3t3 c3t3; // tridimensional mesh with costrained surfaces
    std::vector<Point> points_; // original points
    std::vector<Point> rpoints_; // reconstruction points
    
    typedef std::array<int, 3>    mface;
    typedef std::array<float, 3> mvector;
    typedef std::array<float, 3>  mpoint;
    std::vector<mpoint> fl_points_;
    std::vector<mface> faces_;
    std::vector<mvector> normals_;
    std::unique_ptr<Reconstruction> reconstruction_;
    void print_interface();
    void print_stl();
    void print_interface_oriented(); // print off and medit oriented surface
    void check_normals();
};
#endif
