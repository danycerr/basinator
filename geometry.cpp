// #include<basin_generator.hpp>
// 
// void basin::create_domain(){
//   std::cout<<"Reading sufraces"<<std::endl;
//   Polyhedron_ei top;
//   std::ifstream input_top("./horizon_1.off");
//   input_top>>top;
//   Nef_polyhedron N_top(top);
//   //Polyhedron_ei bot;
//   //std::ifstream input_bot("./horizon_3.off");
//   //input_bot>>bot;
//   //Nef_polyhedron N_bot(top);
//   std::cout<<"---geometry---"<<std::endl;
//   Nef_polyhedron N1(Plane_3( 1, 0, 0,-280000));
//   Nef_polyhedron N2(Plane_3(-1, 0, 0, 220000));
//   Nef_polyhedron N3(Plane_3( 0, 1, 0,-130000));
//   Nef_polyhedron N4(Plane_3( 0,-1, 0, -80000));
//   Nef_polyhedron N5(Plane_3( 0, 0, 1,-5000));
//   Nef_polyhedron N6(Plane_3( 0, 0,-1,-1000)); 
//   // Nef_polyhedron I1(!N1 + !N2);  // open slice in yz-plane
//   // Nef_polyhedron I2(N3 - !N4);   // closed slice in xz-plane
//   // Nef_polyhedron I3(N5 ^ N6);    // open slice in yz-plane
//   // Nef_polyhedron Cube1(I2 * !I1);
//   // Cube1 *= !I3;
//   // Nef_polyhedron Cube2 = N1 * N2 * N3 * N4 * N5 * N6;
//    Nef_polyhedron Cube2 = N1 * N2 * N3 * N4 * N_top * N6;
//   //  Nef_polyhedron Cube3 = N_top ;
//   std::cout<<"conversion"<<std::endl;
//   Polyhedron_ei P;
//   Cube2.convert_to_polyhedron(P);
//   std::cout<<P<<std::endl;
//   std::ofstream out_poly;
//   out_poly.open("nef_poly.off");
//    out_poly << P;
//   out_poly.close();
//   // Mesh_domain_ek domain(P);
//   // output the result into a Surface_mesh
// //   Surface_mesh output;
// //   CGAL::convert_nef_polyhedron_to_polygon_mesh(Cube2, output);
// //   std::ofstream out;
// //   out.open("nef.off");
// //   out << output;
// //   out.close();
//   
//   std::ifstream input("nef_poly.off");
//   Polyhedron poly_bounding_buf;
//   input>>poly_bounding_buf;
//   
//   
//   Mesh_domain domain_costrained_buf (poly_bounding_buf);
//   domain_costrained_buf.detect_features();
//     std::cout<<"End of geneartion of poligonal costrained mesh"<<std::endl;
//     // Mesh criteria
//     Mesh_criteria criteria(//edge_size = 0.25,
//                          facet_angle = 25
//                          // , facet_size = 0.05, facet_distance = 0.005,
//                          // cell_radius_edge_ratio = 3, cell_size = 0.05
//                           );
//   
//     // Mesh generation
//     std::cout<<"Geneartion of poligonal costrained mesh"<<std::endl;
//     Timer t;
//     t.start();
//     c3t3 = CGAL::make_mesh_3<C3t3>(domain_costrained_buf, criteria,
//                                         no_perturb(), no_exude());
//     std::cerr << "Generation of 3d mesh costrained in "<< t.time() << " sec." << std::endl;
//     // Output
//     dump_c3t3(c3t3, "../3dbasin_nef");
//   //CGAL_assertion(Cube1 == Cube2);  // both are closed cube
//   //CGAL_assertion(Cube1 == Cube1.closure());
//   //CGAL_assertion(Cube1 == Cube1.regularization());
//   //CGAL_assertion((N1 - N1.boundary()) == N1.interior());
//   //CGAL_assertion(I1.closure() == I1.complement().interior().complement());
//   //CGAL_assertion(I1.regularization() == I1.interior().closure());
// }
