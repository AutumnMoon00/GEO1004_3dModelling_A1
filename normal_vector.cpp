#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Vector_3 Vector_3;

int main()
{
    // Define three points in 3D space
    Point_3 p1(0, 0, 0);
    Point_3 p2(1, 0, 0);
    Point_3 p3(0, 1, 0);

    // Create a plane from the three points
    Plane_3 plane(p1, p2, p3);

    // Calculate the normal vector of the plane using the normal() function
//    Vector_3 normal_vector = CGAL::normal(plane);
    Vector_3 normal_vector = CGAL::normal(p1, p2, p3);
    std::cout << "The normal vector of the plane is: " << normal_vector << std::endl;
    return 0;
}
