#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <fstream>
#include "sphere.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

int main()
{
  std::vector <K::Point_3> p0;
  make_unit_sphere<K>(p0, 20);

  std::ofstream op("/tmp/sphere.cloud");
  for (size_t i = 0; i < p0.size(); ++i)
    op << p0[i] << "\n";

  std::ofstream ow("/tmp/sphere.w");
  for (size_t i = 0; i < p0.size(); ++i)
    ow << "1\n";
}
