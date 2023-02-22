#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
//#include <string>

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

const std::string input_file = "station.hw1";
//const std::string input_file = "D:\\Geomatics\\Q3\\GEO1004 3D Modelling\\assignments\\1\\station.hw1";
const std::string output_file = "D:/Geomatics/Q3/GEO1004 3D Modelling/assignments/1/output/station.obj";

struct Vertex {
    int id;
    double x, y, z;
};

struct Face {
    int fid;  // ADDED BY ME
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
        std::getline(input_stream, line);
        std::cout << "Vertex header: " << line << std::endl;

        std::istringstream vertex_header_stream(line);
        int number_of_vertices;
        vertex_header_stream >> number_of_vertices;

        std::cout << "Parsing " << number_of_vertices << " vertices..." << std::endl;
        // Read Vertices
        for (int i = 0; i < number_of_vertices; ++i) {
            std::getline(input_stream, line);
            std::istringstream line_stream(line);
            int id;
            double x, y, z;
            line_stream >> id >> x >> y >> z;
//            std::cout << "Vertex " << id << ": (" << x << ", " << y << ", " << z << ")" << std::endl;
            vertices[id].id = id;  // ADDED BY ME
            vertices[id].x = x;
            vertices[id].y = y;
            vertices[id].z = z;
        }

        // Read face header
        std::getline(input_stream, line);
        std::cout << "Face header: " << line << std::endl;

        std::istringstream face_header_stream(line);
        int number_of_faces;
        face_header_stream >> number_of_faces;
        std::cout << "Parsing " << number_of_faces << " faces..." << std::endl;

        // Read Faces
        for (int j = 0; j < number_of_faces; j++) {
            faces[j].fid = j;  // ADDED BY ME
            std::getline(input_stream, line);
            std::istringstream line_stream(line);

            int outer_ring {1}, inner_rings {}, total_rings {};
            line_stream >> outer_ring >> inner_rings;
            total_rings = outer_ring + inner_rings;

            for (int k = 0; k < total_rings; k++) {
                std::getline(input_stream, line);
                std::istringstream line_stream(line);
                int num_vertices {};
                line_stream >> num_vertices;
                if (k == 0) {
                    // outer ring
                    // adding vertices to outer ring list
                    for (int l = 0; l < num_vertices; l++) {
                        int vertex_id;
                        line_stream >> vertex_id;
                        faces[j].outer_ring.emplace_back(vertex_id);
                    }
                }
                else {
                    // inner ring
                    std::list<int> inner_ring_vertices = {};
                    // adding vertices to inner ring list
                    for (int l = 0; l < num_vertices; l++) {
                        int vertex_id;
                        line_stream >> vertex_id;
                        inner_ring_vertices.emplace_back(vertex_id);
                    }
                    faces[j].inner_rings.emplace_back(inner_ring_vertices);
                }
            }
        }
    }

    for (const auto& [key, face]: faces) {
        std::cout << "\nFace id: " << face.fid << std::endl << "\touter ring vertices: ";
        for (const auto& vertex_outer: face.outer_ring) {
            std::cout << vertex_outer << " ";
        }
        for (const auto& inner_ring: face.inner_rings) {
            std::cout << "\n\tinner ring vertices: ";
            for (const auto& vertex_inner: inner_ring) {
                std::cout << vertex_inner << " ";
            }
        }
    }
    return 0;
}
