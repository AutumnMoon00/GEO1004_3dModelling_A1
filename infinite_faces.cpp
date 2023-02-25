//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Triangulation_2.h>
//
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Triangulation_2<K> Delaunay;
//
//int main()
//{
//    // Define a set of 2D points
//    std::vector<K::Point_2> points = {
//            K::Point_2(0, 0),
//            K::Point_2(1, 0),
//            K::Point_2(0, 1),
//            K::Point_2(1, 1),
//            K::Point_2(0.5, 0.5)
//    };
//
//    // Create a Delaunay triangulation from the input points
//    Delaunay dt;
//    dt.insert(points.begin(), points.end());
//
//    // Iterate over the faces and count the number of infinite faces
//    int num_infinite_faces = 0;
//    int count {0};
//    for (auto it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it) {
//        count++;
//        if (dt.is_infinite(it)) {
//            ++num_infinite_faces;
//        }
//    }
//
//    std::cout << "Number of infinite faces: " << num_infinite_faces << std::endl;
//    std::cout << "count: " << count << std::endl;
//    return 0;
//}



// //////////////////////
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

int main()
{
    // Define a set of 2D points
    std::vector<K::Point_2> points = {
            K::Point_2(0, 0),
            K::Point_2(1, 0),
            K::Point_2(0, 1),
            K::Point_2(1, 1),
            K::Point_2(0.5, 0.5)
    };

    // Create a Delaunay triangulation from the input points
    Delaunay dt;
    dt.insert(points.begin(), points.end());

    // Loop through the faces of the triangulation
    for (auto it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it) {
        // Get the vertices of the current face
        auto v0 = it->vertex(0)->point();
        auto v1 = it->vertex(1)->point();
        auto v2 = it->vertex(2)->point();

        // Do something with the vertices, e.g. print them to the console
        std::cout << "Face: " << v0 << ", " << v1 << ", " << v2 << std::endl;
    }

    std::cout << "==============================" << std::endl;
    for (auto it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it) {
        if (dt.is_infinite(it)) {
            // Get the vertices of the current face
            auto v0 = it->vertex(0)->point();
            auto v1 = it->vertex(1)->point();
            auto v2 = it->vertex(2)->point();

            // Do something with the vertices, e.g. print them to the console
            std::cout << "Infinite face: " << v0 << ", " << v1 << ", " << v2 << std::endl;
        }
    }

    auto infinite_face = dt.infinite_face();
    for (int i = 0; i < 3; ++i) {
        auto vertex = infinite_face->vertex(i);
        std::cout << "Vertex " << i << ": " << vertex->point() << std::endl;
    }

     return 0;
}


void label_triangles() {
    std::list<Triangulation::Face_handle> to_check;
    triangulation.infinite_face()->info().processed = true;
    CGAL_assertion(triangulation.infinite_face()->info().processed == true);
    CGAL_assertion(triangulation.infinite_face()->info().interior == false);
    to_check.push_back(triangulation.infinite_face());
    while (!to_check.empty()) {
        CGAL_assertion(to_check.front()->info().processed == true);  // useless
        for (int neighbour = 0; neighbour < 3; ++neighbour) {
            if (to_check.front()->neighbor(neighbour)->info().processed == true)  {
                // Note: validation code.
//          if (triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour)))
//              CGAL_assertion(to_check.front()->neighbor(neighbour)->info().interior != to_check.front()->info().interior);
//          else
//              CGAL_assertion(to_check.front()->neighbor(neighbour)->info().interior == to_check.front()->info().interior);
            } else {
                to_check.front()->neighbor(neighbour)->info().processed = true;
                CGAL_assertion(to_check.front()->neighbor(neighbour)->info().processed == true);
                if (triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) {

                    to_check.front()->neighbor(neighbour)->info().interior = !to_check.front()->info().interior;
                    to_check.push_back(to_check.front()->neighbor(neighbour));
                } else {
                    to_check.front()->neighbor(neighbour)->info().interior = to_check.front()->info().interior;
                    to_check.push_back(to_check.front()->neighbor(neighbour));
                }
            }
        } to_check.pop_front();
    }
}