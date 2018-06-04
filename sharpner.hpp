#ifndef SHARPNERHS
#define SHARPNERHS

#include <iostream>
#include <fstream>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::ostringstream

//Cgal includes
// polyedron must be on top for parameters configuration
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Polyhedron_3.h>
// Triangulation
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/read_off_points.h>
#include <CGAL/Timer.h>

#include <exporter.hpp>
using namespace CGAL::parameters;
class sharpner{
    // cgal typedef
    typedef CGAL::Exact_predicates_inexact_constructions_kernel         Kernel;
    //Polyedron
    typedef CGAL::Mesh_polyhedron_3<Kernel>::type                       Polyhedron;
    typedef Kernel::Point_3                                             Point3;
    typedef Kernel::Vector_3                                            Vector3;
    // Triangulation
    typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel>        Mesh_domain;
    typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type               Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<
      Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index>    C3t3;
    // Criteria
    typedef CGAL::Mesh_criteria_3<Tr>                                   Mesh_criteria;
    typedef Mesh_criteria::Edge_criteria                                Edge_criteria;
    typedef Mesh_criteria::Facet_criteria                               Facet_criteria;
    typedef Mesh_criteria::Cell_criteria                                Cell_criteria;
    typedef CGAL::Timer                                                 Timer;
    typedef std::vector<Point3>                                         Polyline_3;
    typedef std::list<Polyline_3>                                       Polylines;
    typedef Polyhedron::Vertex_iterator                                 Vertex_iterator;
public:
    sharpner(std::string hname);
private:
    std::ostringstream face_off_name_; // name for the off file of interface
    std::ostringstream face_sharp_off_name_; // name for the off file of interface oriented
    std::ostringstream face_sharp_mesh_name_; // name for the off file of interface oriented
    std::string path_="./";
    void create_domain();
    void bouding_box(
        Polyhedron &poly,          /// <- polyedron to build the bounding box
        std::vector<double>& b_box /// -> vector bounding box points
                        );
    ////////////////////////////////////////////
    //creation of general quadrilateral features
    ////////////////////////////////////////////
    void insert_features_of_quadrilateral_face(
        const char* off_filename, /// off file name
        Polylines& polylines      /// polyline to protect
    );
    ////////////////////////////////////////////
    //creation of general quadrilateral features
    ////////////////////////////////////////////
    void create_polyline(       
        std::vector<Point3>& points,
        Polylines& polylines
                            );
};

#endif
