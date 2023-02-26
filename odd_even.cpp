#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>


#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Exact_predicates_tag Tag;
struct FaceInfo {
    bool processed;
    bool interior;
    FaceInfo() {
        processed = false;
        interior = false;
//        interior = true;
    }
};
typedef CGAL::Triangulation_vertex_base_2<Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> Triangulation;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Point_2 Point_2;

typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
typedef Triangulation::Finite_faces_iterator FaceIterator;


struct Vertex {
    int id;
    double x, y, z;
};

struct Face {
    int fid;  // ADDED BY ME
    std::list<int> outer_ring;
    std::list<std::list<int>> inner_rings;
    Plane_3 best_plane;
    Triangulation triangulation;
};

typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;

const std::string input_file = "station.hw1";
const std::string output_file = "D:/Geomatics/Q3/GEO1004 3D Modelling/assignments/1/output/station.obj";

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

    // make the best fitting plane using all the points in both the inner and outer ring vertices
    int total_triangles = 0;
    for (auto& [key, face]: faces) {
        std::vector<Point_3> face_vertices{};

        std::cout << "\nFace id: " << face.fid << std::endl << "\touter ring vertices: ";
        for (const auto &vertex_outer: face.outer_ring) {
            std::cout << vertex_outer << " ";
            face_vertices.emplace_back(
                    Point_3(vertices[vertex_outer].x, vertices[vertex_outer].y, vertices[vertex_outer].z));
        }
        for (const auto &inner_ring: face.inner_rings) {
            std::cout << "\n\tinner ring vertices: ";
            for (const auto &vertex_inner: inner_ring) {
                std::cout << vertex_inner << " ";
                face_vertices.emplace_back(
                        Point_3(vertices[vertex_inner].x, vertices[vertex_inner].y, vertices[vertex_inner].z));
            }
        }

        // best fitting plane - assigned to the face
        FT rms = CGAL::linear_least_squares_fitting_3(face_vertices.begin(), face_vertices.end(), faces[key].best_plane,
                                                      CGAL::Dimension_tag<0>());
        // copying a plane to a new plane to not fuck up by mistake the original plane
        Plane_3 plane = faces[key].best_plane;

        // Output the plane equation and RMS error
        std::cout << "\n\tnumber of vertices: " << face_vertices.size();
        std::cout << "\n\tPlane equation: " << faces[key].best_plane << std::endl;
        std::cout << "\tRMS error: " << rms << std::endl;

        // projecting points onto best_fitting plane
        std::vector<Point_2> proj_points {};
        for (const auto& pt_3: face_vertices) {
            proj_points.emplace_back(plane.to_2d(pt_3));
//            std::cout << plane.to_2d(pt_3) << std::endl;
        }

        std::cout << "\t==================" << std::endl;
        std::vector<Point_2> outer_points_2 {};  // all the vertices on the outer ring
        for (const auto& vertex_outer: face.outer_ring) {
            outer_points_2.emplace_back(plane.to_2d(Point_3(vertices[vertex_outer].x, vertices[vertex_outer].y, vertices[vertex_outer].z)));
//            std::cout << "\n\t outer vertex: ";
//            std::cout << plane.to_2d(Point_3(vertices[vertex_outer].x, vertices[vertex_outer].y, vertices[vertex_outer].z)) << std::endl;
        }
        face.triangulation.insert_constraint(outer_points_2.begin(), outer_points_2.end(), true);
        face.triangulation.insert(outer_points_2.begin(), outer_points_2.end());

        for (const auto& inner_ring: face.inner_rings) {
            std::vector<Point_2> inner_points_2 {};
            for (const auto& vertex_inner: inner_ring) {
                inner_points_2.emplace_back(plane.to_2d(Point_3(vertices[vertex_inner].x, vertices[vertex_inner].y, vertices[vertex_inner].z)));
                // printing the vertices
//                std::cout << "\n\t inner vertex: ";
//                std::cout << plane.to_2d(Point_3(vertices[vertex_inner].x, vertices[vertex_inner].y, vertices[vertex_inner].z)) << std::endl;
            }
            face.triangulation.insert_constraint(inner_points_2.begin(), inner_points_2.end(), true);
            face.triangulation.insert(inner_points_2.begin(), inner_points_2.end());
        }


        // looping through the vertices to see whether all the vertices are loaded
//        for (Finite_vertices_iterator vit = face.triangulation.finite_vertices_begin(); vit != face.triangulation.finite_vertices_end(); ++vit)
//            std::cout << "\tVertex " << ": " << vit->point() << std::endl;

        int num_triangles = 0;
        for (auto f = face.triangulation.finite_faces_begin(); f != face.triangulation.finite_faces_end(); ++f) {
            num_triangles++;
            total_triangles++;
        }
        std::cout << "\tNumber of triangles: " << num_triangles << std::endl;
        std::cout << "\ttotal number of triangles: " << total_triangles;


        // depth first search to lable the faces as interior and exterior
        // code from Ken
        std::list<Triangulation::Face_handle> to_check {};
        face.triangulation.infinite_face()->info().processed = true;

        CGAL_assertion(face.triangulation.infinite_face()->info().processed == true);
        CGAL_assertion(face.triangulation.infinite_face()->info().interior == false);
        to_check.push_back(face.triangulation.infinite_face());
        while (!to_check.empty()) {
            CGAL_assertion(to_check.front()->info().processed == true);
            for (int neighbour = 0; neighbour < 3; ++neighbour) {
                if (to_check.front()->neighbor(neighbour)->info().processed == true) {
                    // Note: validation code.
                    if (face.triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour)))
                        CGAL_assertion(to_check.front()->neighbor(neighbour)->info().interior != to_check.front()->info().interior);
                    else
                        CGAL_assertion(to_check.front()->neighbor(neighbour)->info().interior == to_check.front()->info().interior);
                }
                else {
                    to_check.front()->neighbor(neighbour)->info().processed = true;
                    CGAL_assertion(to_check.front()->neighbor(neighbour)->info().processed == true);
                    // to see triangle 1
//                    std::cout << "\n\t\ttriangle 1 spotted" << std::endl;
                    if (face.triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) {
//                        std::cout << "\n\t\ttriangle 1 spotted" << std::endl;
                        to_check.front()->neighbor(neighbour)->info().interior = !to_check.front()->info().interior;
                        to_check.push_back(to_check.front()->neighbor(neighbour));
                    }
                    else {
                        to_check.front()->neighbor(neighbour)->info().interior = to_check.front()->info().interior;
                        to_check.push_back(to_check.front()->neighbor(neighbour));
                    }
                }
            } to_check.pop_front();
        }


        int tri_count {0};
        for (FaceIterator fit = face.triangulation.finite_faces_begin(); fit != face.triangulation.finite_faces_end(); ++fit)
        {
            std::cout << "\n\ttringle: " << tri_count << "\tinterior: " << fit->info().interior;
            tri_count++;
        }


    }
    std::cout << "\ntotal number of triangles: " << total_triangles;
    std::cout << "\n=====================================" << std::endl;

    return 0;
}
