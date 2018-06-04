#ifndef HPATCHER_HS
#define HPATCHER_HS
#include <basin_generator.hpp>
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

//exporting 
#include <exporter.hpp>
#include<sharpner.hpp>
class horizon;
class patcher{
typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >                Reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother< Kernel > Smoother;
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
    patcher(horizon &h1, horizon &h2, int label);
private:
    void find_min_max(horizon   &h1,
                      double& xmin,
                      double& xmax,
                      double& ymin,
                      double& ymax);
    //////////////////////////////////////////////
    void add_points   (double trashold, // trashold value
                       int dir,         //direction
                       int method,      // 0 min 1 max
                       horizon &h,      //cloud points
                       std::vector<Point>& container);
    /////////////////////////////////////////////////
    void build_lateral_surfaces(std::string hname,
                                 std::vector<Point> &points_);
    std::unique_ptr<Reconstruction> reconstruction_;
    std::vector<Point> points_xmin_;
    std::vector<Point> points_ymin_;
    std::vector<Point> points_xmax_;
    std::vector<Point> points_ymax_;
    int name;
    // horizon top,bot;
};
#endif
