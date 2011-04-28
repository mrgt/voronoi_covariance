#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <fstream>
#include <vector>
#include <CGAL/voronoi_covariance_sphere_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

int main()
{
  typedef K::Point_3 Point;
  std::vector <Point> p;
  CGAL::internal::make_icosahedron<Point>(1.0, std::back_inserter(p));
  //make_dodecahedron<Point>(1.0, std::back_inserter(p));

  std::ofstream op("/tmp/sphere.cloud");
  for (size_t i = 0; i < p.size(); ++i)
    op << p[i] << "\n";

  std::ofstream ow("/tmp/sphere.w");
  for (size_t i = 0; i < p.size(); ++i)
    ow << "1\n";
}
