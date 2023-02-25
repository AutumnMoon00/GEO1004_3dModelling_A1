#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

struct My_face_base : public CGAL::Triangulation_face_base_2<CGAL::Exact_predicates_inexact_constructions_kernel> {
    int my_index;

    My_face_base() : my_index(-1) {}
    My_face_base(const My_face_base& other) : my_index(other.my_index) {}
    My_face_base& operator=(const My_face_base& other) {
        my_index = other.my_index;
        return *this;
    }

    int index() const { return my_index; }
};

using My_triangulation = CGAL::Triangulation_2<CGAL::Exact_predicates_inexact_constructions_kernel, My_face_base>;

int main() {
    My_triangulation triangulation;
    // add vertices and faces
    for (auto face_it = triangulation.finite_faces_begin(); face_it != triangulation.finite_faces_end(); ++face_it) {
        int face_index = face_it->index();
        // do something with the face index
    }
    return 0;
}
