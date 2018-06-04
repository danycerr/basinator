#include "exporter.hpp"

exporter::exporter(const char* off_filename):outname_(off_filename){
    std::cout<<"Exporter::building and exporter  "<<std::endl;
    std::cout<<off_filename <<std::endl;
    // Check if the last three characters match the ext.
   const std::string ext(".off");
   if ( outname_ != ext &&
     outname_.size() > ext.size() &&
     outname_.substr(outname_.size() - ext.size()) == ".off" )
   {
   // if so then strip them off
   outname_ = outname_.substr(0, outname_.size() - ext.size());
   }
   std::cout<<"Exporter::outname is "<< outname_<<std::endl;
    std::cout<<"Exporter::reading off file"<<std::endl;
    std::ifstream off_file;
    off_file.open (off_filename);
    char buff[100];
    off_file>>buff; // getting rid of OFF initial line
    int nb_points=0;int nb_facet=0;
    off_file >> nb_points;
    off_file >> nb_facet;
     off_file>>buff; // getting rid of a dummy number
    std::cout<<"Exporter::Number of points of "<< off_filename <<" is "<<nb_points<<std::endl;
    std::cout<<"Exporter::Number of facets of "<< off_filename <<" is "<<nb_facet<<std::endl;
    fl_points_.resize(nb_points);
    faces_.resize(nb_facet); // sizing for faces
    normals_.resize(nb_facet); // sizing for faces
    for( int it=0;  it<nb_points; it++ )
  	{
        off_file >> (fl_points_[it])[0];
        off_file >> (fl_points_[it])[1];
        off_file >> (fl_points_[it])[2];
    }
    
    for( int it=0;  it<nb_facet; it++ )
    { 
        off_file>>buff; // getting rid of a dummy number
        off_file >> (faces_[it])[0]; 
        off_file >> (faces_[it])[1]; 
        off_file >> (faces_[it])[2]; 
    }
    for( int it=0;  it<nb_facet; it++ ) 
    {
    double a[3],b[3],n[3];
       for( int idim=0;  idim<3; idim++ )
       { 
          a[idim]=(fl_points_[(faces_[it])[2]])[idim]-(fl_points_[(faces_[it])[0]])[idim];
          b[idim]=(fl_points_[(faces_[it])[1]])[idim]-(fl_points_[(faces_[it])[0]])[idim];
       }
       double mod=sqrt((a[1]*b[2]-a[2]*b[1])*(a[1]*b[2]-a[2]*b[1])+
                  (a[0]*b[2]-a[2]*b[0])*(a[0]*b[2]-a[2]*b[0])+
                  (a[0]*b[1]-a[1]*b[0])*(a[0]*b[1]-a[1]*b[0]));
       (normals_[it])[0]=(float)(a[1]*b[2]-a[2]*b[1])/mod;
       (normals_[it])[0]=(float)(a[0]*b[2]-a[2]*b[0])/mod;
       (normals_[it])[0]=(float)(a[0]*b[2]-a[2]*b[0])/mod;
    }
    print_stl();
    std::cout<<"Exporter end printing stl"<<std::endl;    
}



void exporter::print_stl(){
//binary file
//     std::string filename="prova";
std::cout<<"stl file printing"<<(outname_ + ".stl").c_str()<<std::endl;
std::string header_info = "solid " + outname_ + "-output";
char head[80];
std::strncpy(head,header_info.c_str(),sizeof(head)-1);
char attribute[2] = "0";
uint32_t nTriLong =(uint32_t)  faces_.size();

std::ofstream myfile;

myfile.open((outname_ + ".stl").c_str(),  std::ios::out | std::ios::binary);
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
