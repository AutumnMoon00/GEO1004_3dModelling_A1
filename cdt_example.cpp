#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_3<K> PT;
typedef CGAL::Plane_3<K> Plane;
typedef CGAL::Point_3<K> Point_3;
typedef CGAL::Point_2<K> Point_2;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<PT, CGAL::Default, Itag> CDT;

int main()
{
    // Define the plane using a, b, c, and d values
    K::FT a = 1, b = 2, c = -1, d = 3;
    Plane plane(a, b, c, d);

    // Define a normal vector to the plane
    K::Vector_3 normal(a, b, c);

    // Create a point on the plane
    K::Point_3 p(0, 0, 0);
    K::Point_3 q = p + normal;

    // Create the projection traits object
    PT pt(q);


    // Define the four points
    std::vector<Point_3> ps(4);
    ps[0] = Point_3(0, 0, 0);
    ps[1] = Point_3(3, 1, -1);
    ps[2] = Point_3(-1, 3, -3);
    ps[3] = Point_3(1, 0.5, -0.5);

    // Project the points onto the plane
    std::vector<Point_2> projected_ps(4);
    for(int i = 0; i < 4; ++i)
        projected_ps[i] = plane.to_2d(ps[i]);

//    PT pt {plane};
    CDT cdt(pt);

    // Insert the points and constrained edges into the CDT
    for(int i = 0; i < 4; ++i)
        cdt.insert(projected_ps[i]);
    for(int i = 1; i < 3; ++i)
        cdt.insert_constraint(projected_ps[i], projected_ps[i+1]);

    // Print the triangles in the CDT
    for(CDT::Face_handle f : cdt.all_face_handles())
    {
        for(int i = 0; i < 3; ++i)
            std::cout << plane.to_3d(f->vertex(i)->point()) << " ";
        std::cout << std::endl;
    }

    return 0;
}
