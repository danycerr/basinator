#ifndef EXPORTERHS
#define EXPORTERHS
#include <iostream>
#include <vector>
#include <array>

#include <math.h>       /* sqrt */
#include<fstream>
// #include <string.h>
#include <cstring>
class exporter{
    typedef std::array<int, 3>    mface; ///type for int
    typedef std::array<float, 3> mvector;
    typedef std::array<float, 3>  mpoint;
    std::vector<mpoint>  fl_points_;
    std::vector<mface>   faces_;
    std::vector<mvector> normals_;
public:
    exporter(const char* off_filename);
private:
    std::string outname_;
    void print_stl();
};
#endif
