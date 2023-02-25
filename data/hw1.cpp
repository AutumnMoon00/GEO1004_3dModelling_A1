#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
struct FaceInfo {
  bool interior;
  FaceInfo() {
    interior = false;
  }
};
typedef CGAL::Triangulation_vertex_base_2<Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> Triangulation;

const std::string input_file = "D:/Geomatics/Q3/GEO1004 3D Modelling/assignments/1/data/station.hw1";
const std::string output_file = "D:/Geomatics/Q3/GEO1004 3D Modelling/assignments/1/output/station.obj";

struct Vertex {
  int id;
  double x, y, z;
};

struct Face {
  std::list<int> outer_ring;
  std::list<std::list<int>> inner_rings;
  Kernel::Plane_3 best_plane;
  Triangulation triangulation;
};

int main(int argc, const char * argv[]) {
  
  std::map<int, Vertex> vertices;
  std::map<int, Face> faces;
  
  // Read file
  std::ifstream input_stream;
  input_stream.open(input_file);
  if (input_stream.is_open()) {
    std::string line;
    
    // Read vertex header
    getline(input_stream, line);
    std::cout << "Vertex header: " << line << std::endl;
    std::istringstream vertex_header_stream(line);
    int number_of_vertices;
    vertex_header_stream >> number_of_vertices;
    std::cout << "Parsing " << number_of_vertices << " vertices..." << std::endl;
    
    // Read vertices
    for (int i = 0; i < number_of_vertices; ++i) {
      getline(input_stream, line);
      std::istringstream line_stream(line);
      int id;
      double x, y, z;
      line_stream >> id >> x >> y >> z;
      std::cout << "Vertex " << id << ": (" << x << ", " << y << ", " << z << ")" << std::endl;
      vertices[id].x = x;
      vertices[id].y = y;
      vertices[id].z = z;
    }
    
    // TO DO
  }
  
  return 0;
}
