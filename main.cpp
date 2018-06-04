#include <basin_generator.hpp>
#include <patcher.hpp>
#include <domain.hpp>

int main(int argc, char* argv[]){
    std::cout<<"Welcome to the basin generator"<<std::endl;
    std::cout<<"we have "<<argc<<" horizons"<<std::endl;
    std::vector<horizon> h;
    for(int ih=0;ih<argc-1;ih++){
        h.push_back(horizon (argv[ih+1]));
    }
    for(int ih=0; ih<h.size();ih++){
    std::cout<<"*** Building horizon "<<ih<<" ***"<<std::endl;
    h[ih].read_horizons_points(argv[ih+1]);
    h[ih].reconstruction_surface();
    }
    
    patcher p1(h[0],h[1], 10);
    patcher p2(h[1],h[2], 20);
    domain d;
//     d.create_horizon((h[0].hname + ".off").c_str());
    d.create_horizon((h[1].hname + ".off").c_str());
//     std::vector<int> hs={0,1};
//     d.lateral_face_building(hs);
//     d.create_single_domain(); // create single domain shell patsces
//     b.create_poligonal_domain(argv[2]);
 return 1;   
}
