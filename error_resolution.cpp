#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <vector>
#include <exception>  // for std::exception
#include <CGAL/exceptions.h>
#include <typeinfo>
#include <cmath>
#include <algorithm>
#include <iomanip>


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
typedef Kernel::Segment_2 Segment_2;

typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
typedef Triangulation::Finite_faces_iterator FaceIterator;
typedef Triangulation::Finite_edges_iterator EdgeIterator;
typedef Triangulation::Face_handle Face_handle;


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
    std::vector<std::vector<int>> interior_triangles;
};

typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;

const std::string input_file = "stilted_house.hw1";
const std::string output_file = "stilted_house_cpp.obj";

void read_vertices_faces(std::map<int, Vertex>& vertices, std::map<int, Face>& faces) {
//    std::map<int, Vertex> vertices;
//    std::map<int, Face> faces;

    struct retVals{
        std::map<int, Vertex> ret_v;
        std::map<int, Face> ret_f;
    };

    // Read file
    std::ifstream input_stream;
    input_stream.open(input_file);

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
        std::cout << "Vertex " << id << ": (" << x << ", " << y << ", " << z << ")" << std::endl;
//            std::cout << "v " << " " << x << " " << y << " " << z << std::endl;
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

        std::cout << "\nface: " << j << " Outer rings: " << outer_ring << " - Inner rings: " << inner_rings << std::endl;


        for (int k = 0; k < total_rings; k++) {
            std::getline(input_stream, line);
            std::istringstream line_stream(line);
            int num_vertices {};
            line_stream >> num_vertices;
            if (k == 0) {
                // outer ring
                // adding vertices to outer ring list
                std::cout << "    outer vertices: ";
                for (int l = 0; l < num_vertices; l++) {
                    int vertex_id;
                    line_stream >> vertex_id;
                    faces[j].outer_ring.emplace_back(vertex_id);
                    std::cout << vertex_id << " ";

                }
            }
            else {
                // inner ring
                std::list<int> inner_ring_vertices = {};
                // adding vertices to inner ring list
                std::cout << "\n\tInner ring: " << k << ": ";
                for (int l = 0; l < num_vertices; l++) {
                    int vertex_id;
                    line_stream >> vertex_id;
                    inner_ring_vertices.emplace_back(vertex_id);
                    std::cout << vertex_id << " ";
                }
                faces[j].inner_rings.emplace_back(inner_ring_vertices);
            }
        }
    }
//    return retVals{vertices, faces};
}

void make_best_fitting_plane(std::map<int, Vertex>& vertices, std::map<int, Face>& faces, std::vector<int>& invalid_face_indices) {
    // make the best fitting plane using all the points in both the inner and outer ring vertices
    int total_triangles {}, interior_triangles {}, exterior_triangles {};
    for (auto &[key, face]: faces) {
        std::vector<Point_3> face_vertices{};
        for (const auto &vertex_outer: face.outer_ring) {
            face_vertices.emplace_back(
                    Point_3(vertices[vertex_outer].x, vertices[vertex_outer].y, vertices[vertex_outer].z));
        }
        for (const auto &inner_ring: face.inner_rings) {
            for (const auto &vertex_inner: inner_ring) {
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
        std::cout << "\nface id: " << key ;
        std::cout << "\n\tnumber of vertices: " << face_vertices.size();
        if (std::isnan(rms)) {
            std::cout << "\n\tFACE IS INVALID";
            invalid_face_indices.emplace_back(key);
            continue;
        }
        std::cout << "\n\tPlane equation: " << faces[key].best_plane;
        std::cout << "\n\tRMS error: " << rms;
        // /////////////////////////////////////////////////////////

        // projecting points onto best_fitting plane
        std::vector<Point_2> proj_points {};
        for (const auto& pt_3: face_vertices)
            proj_points.emplace_back(plane.to_2d(pt_3));

        Point_2 pt_o, pt_i;
        std::map<Point_2, int> proj_2d_to_vertex;

        // inserting outer ring vertices and constraints into triangulation
        std::vector<Point_2> outer_points_2 {};  // all the vertices on the outer ring
        for (const auto& vertex_outer: face.outer_ring) {
            pt_o = plane.to_2d(Point_3(vertices[vertex_outer].x, vertices[vertex_outer].y, vertices[vertex_outer].z));
            outer_points_2.emplace_back(pt_o);
            proj_2d_to_vertex[pt_o] = vertex_outer;
        }
//        face.triangulation.set_max_number_of_steiner_points(0);
        face.triangulation.insert(outer_points_2.begin(), outer_points_2.end());
        face.triangulation.insert_constraint(outer_points_2.begin(), outer_points_2.end(), true);

        // inserting inner ring vertices and constraints into triangulation
        for (const auto& inner_ring: face.inner_rings) {
            std::vector<Point_2> inner_points_2 {};
            for (const auto& vertex_inner: inner_ring) {
                pt_i = plane.to_2d(Point_3(vertices[vertex_inner].x, vertices[vertex_inner].y, vertices[vertex_inner].z));
                inner_points_2.emplace_back(pt_i);
                proj_2d_to_vertex[pt_i] = vertex_inner;
            }
            face.triangulation.insert(inner_points_2.begin(), inner_points_2.end());
            face.triangulation.insert_constraint(inner_points_2.begin(), inner_points_2.end(), true);
        }

        // looping through the vertices to see whether all the vertices are loaded
        int i = 1;
        Point_2 point_error_catch {};
        Point_3 error_pt3 {};
        std::cout << "\n\tNumber of vertices CDT: " << face.triangulation.number_of_vertices() << std::endl;

        int additional_vertex = vertices.size()-1;
        for (Finite_vertices_iterator vit = face.triangulation.finite_vertices_begin(); vit != face.triangulation.finite_vertices_end(); ++vit){
            std::cout << "\n\tVertex " << i << ": " << vit->point();
            point_error_catch = vit->point();
            try{
                std::cout << " - " << proj_2d_to_vertex.at(vit->point());
//                std::cout << "\n\t heehehhehhehehe";
                std::cout << "\n\tprojected to 3d: " << plane.to_3d(point_error_catch);
//                std::cout << "\n\tblahprojected to 3d: " << plane.to_3d(vit->point());
                std::cout << "\n\tI'm here - fucked up";
            }
            catch (const std::out_of_range& e) {
                additional_vertex++;
                proj_2d_to_vertex.emplace(point_error_catch, additional_vertex);
                Vertex new_vertex ;
                new_vertex.id = additional_vertex;
                new_vertex.x = plane.to_3d(point_error_catch).x();
                new_vertex.y = plane.to_3d(point_error_catch).y();
                new_vertex.z = plane.to_3d(point_error_catch).z();
                vertices.emplace(additional_vertex, new_vertex);
                std::cout << " - " << proj_2d_to_vertex.at(vit->point());
                std::cout << "\n\tprojected to 3d: " << point_error_catch;
            }

            i++;
        }


        // odd-even labelling
        std::list<Triangulation::Face_handle> to_check{};
        face.triangulation.infinite_face()->info().processed = true;
        to_check.push_back(face.triangulation.infinite_face());
        while (!to_check.empty()) {
            try {
                CGAL_assertion(to_check.front()->info().processed == true);
                for (int neighbour = 0; neighbour < 3; ++neighbour) {
                    if (to_check.front()->neighbor(neighbour)->info().processed) {
                        if (face.triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) {
                            CGAL_assertion(to_check.front()->neighbor(neighbour)->info().interior !=
                                           to_check.front()->info().interior);
                        } else {
                            CGAL_assertion(to_check.front()->neighbor(neighbour)->info().interior ==
                                           to_check.front()->info().interior);
                        }
                    } else {
                        to_check.front()->neighbor(neighbour)->info().processed = true;
                        CGAL_assertion(to_check.front()->neighbor(neighbour)->info().processed == true);
                        if (face.triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) {
                            to_check.front()->neighbor(neighbour)->info().interior = !to_check.front()->info().interior;
                            to_check.push_back(to_check.front()->neighbor(neighbour));
                        } else {
                            to_check.front()->neighbor(neighbour)->info().interior = to_check.front()->info().interior;
                            to_check.push_back(to_check.front()->neighbor(neighbour));
                        }
                    }
                }
            }
            catch (...) {
                continue;
            }
            to_check.pop_front();
        }




//        Point_2 point_outside(1000, 1000);
//        for (FaceIterator fit = face.triangulation.finite_faces_begin(); fit != face.triangulation.finite_faces_end(); ++fit) {
//            int num_of_intersections {0};
//            Point_2 p1 = fit->vertex(0)->point();
//            Point_2 p2 = fit->vertex(1)->point();
//            Point_2 p3 = fit->vertex(2)->point();
//            Point_2 centroid((p1.x() + p2.x() + p3.x()) / 3.0, (p1.y() + p2.y() + p3.y()) / 3.0);
//            std::cout << "Centroid of face: " << centroid << std::endl;
//            Segment_2 line(centroid, point_outside);
//
//            for (EdgeIterator eit = face.triangulation.finite_edges_begin(); eit != face.triangulation.finite_edges_end(); ++eit) {
//                Kernel::Segment_2 edge = face.triangulation.segment(eit);
//                if (face.triangulation.is_constrained(*eit)) {
//                    if (CGAL::do_intersect(line, edge)) {
////                        std::cout << "Line intersects edge!" << std::endl;
//                        num_of_intersections++;
//                    }
//                }
//            }
//
//            if (num_of_intersections % 2 == 0)
//                fit->info().interior = false;
//            else
//                fit->info().interior = true;
//
//        }


        std::cout << "\n\t=======================";
        std::cout << "\n\tnumber of triangles: " << face.triangulation.number_of_faces();
        int tri_count {0};
        Point_2 tri_vertex_2d;
        for (FaceIterator fit = face.triangulation.finite_faces_begin(); fit != face.triangulation.finite_faces_end(); ++fit)
        {
            std::cout << "\n\ttringle: " << tri_count << "\tinterior: " << fit->info().interior;
            std::vector<int> triangle_vertices_numbers {};
            if (fit->info().interior) {
                interior_triangles++;
                Face_handle face_handle = fit;
                for (int k = 0; k < 3; ++k) {
                    tri_vertex_2d = face_handle->vertex(k)->point();
                    triangle_vertices_numbers.emplace_back(proj_2d_to_vertex.at(tri_vertex_2d));
                }
                face.interior_triangles.emplace_back(triangle_vertices_numbers);

            }

            else
                exterior_triangles++;
            tri_count++;
            total_triangles++;
        }


    }

    std::cout << "\n============================";
    std::cout << "\ninvalid faces: ";
    for (const auto& i: invalid_face_indices) {
        std::cout << i << " ";
    }
    std::cout << "\n============================";
    std::cout << "\ntotal triangles: " << total_triangles;
    std::cout << "\ninterior triangles: " << interior_triangles;
    std::cout << "\nexterior triangles: " << exterior_triangles;
}


int main (int argc, const char * argv[]) {

    std::map<int, Vertex> vertices;
    std::map<int, Face> faces;

    // Read file
    std::ifstream input_stream;
    input_stream.open(input_file);

    if (input_stream.is_open()) {
        read_vertices_faces(vertices, faces);

        std::cout << "\n===========================";
        std::cout << "\n\nlength of vertices: " << vertices.size();
        std::cout << "\nlength of faces: " << faces.size();
        std::cout << "\n===========================";

        std::vector<int> invalid_face_indices{};
        make_best_fitting_plane(vertices, faces, invalid_face_indices);
        std::cout << "\n===========================" << std::endl;
        for (const auto &i: invalid_face_indices)
            std::cout << i << " ";

        std::cout << "\n";
        std::ofstream fw(output_file, std::ofstream::out);
        if (fw.is_open()) {
            std::map<Point_3, int> xyz_to_index;
            for (const auto &[i, vertex]: vertices) {
                std::cout << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
                fw << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
                xyz_to_index[Point_3(vertex.x, vertex.y, vertex.z)] = i;
            }
            for (const auto [k, face]: faces) {
                for (const auto vertex_vector: face.interior_triangles) {
                    std::cout << "\nf ";
                    fw << "\nf ";
                    for (auto vertex: vertex_vector) {
                        std::cout << vertex + 1 << " ";
                        fw << vertex + 1 << " ";
                    }
                }
            }

            fw.close();
        }
    }
}