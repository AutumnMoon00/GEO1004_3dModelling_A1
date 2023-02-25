#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;

int main() {
    // Define the input points and constraints
    std::vector<Point> points = {{0, 0}, {1, 0}, {0, 1}, {1, 1}, {0.5, 0.5}};
    std::vector<std::pair<Point, Point>> segments = {{{0, 0}, {1, 1}}, {{0, 1}, {1, 0}}};

    // Construct a constrained Delaunay triangulation
    CDT cdt;
    cdt.insert_constraint(segments.begin(), segments.end());
    for (auto p : points) {
        cdt.insert(p);
    }

    // Output the resulting triangles
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        std::cout << f->vertex(0)->point() << " "
                  << f->vertex(1)->point() << " "
                  << f->vertex(2)->point() << std::endl;
    }

    return 0;
}
