#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Projection_traits_xy_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef CGAL::Projection_traits_xy_3<K> Projection_traits;

int main()
{
    // Define the coefficients of the plane.
    double a = 1.0, b = 2.0, c = 3.0, d = 4.0;

    // Define the 3D point to be projected.
    Point_3 point(1.0, 2.0, 3.0);

    // Define the plane.
    Plane_3 plane(a, b, c, d);

    // Project the point onto the plane.
    Point_3 projected_point = plane.projection(point);

    // Print out the original and projected coordinates.
    std::cout << "Original point: (" << point.x() << ", " << point.y() << ", " << point.z() << ")" << std::endl;
    std::cout << "Projected point: (" << projected_point.x() << ", " << projected_point.y() << ", " << projected_point.z() << ")" << std::endl;

//     Convert the projected point to a 2D point in the plane's coordinate system.
    Projection_traits traits;
//    Projection_traits::Point_2 projected_point_2 = traits(projected_point);
    Projection_traits::Point_2 projected_point_2 = traits(projected_point);
//
//    // Print out the projected point in the plane's coordinate system.
//    std::cout << "Projected point in the plane's coordinate system: (" << projected_point_2.x() << ", " << projected_point_2.y() << ")" << std::endl;

    return 0;
}
