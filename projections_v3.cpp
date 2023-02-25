#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Plane_3.h>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Point_2 Point_2;

int main()
{
    // Define the plane using a,b,c,d values.
    double a = 0.0, b = 0.0, c = 1.0, d = 0.0;
    Plane_3 plane(a, b, c, d);

    // Define the 2D point.
    Point_2 projection(1, 2);

    // Convert the 2D point to a 3D point on the plane.
    Point_3 point = plane.to_3d(projection);

    // Print the 3D point in CGAL representation.
    std::cout << "3D point: " << point << std::endl;

    return 0;
}
