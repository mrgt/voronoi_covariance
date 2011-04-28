#ifndef CGAL_VORONOI_COVARIANCE_SPHERE_3_HPP
#define CGAL_VORONOI_COVARIANCE_SPHERE_3_HPP

#include <CGAL/point_generators_3.h>

CGAL_BEGIN_NAMESPACE

namespace internal
{
  template <class Point, class OutputIterator>
  void
  make_icosahedron (double R,
		    OutputIterator out)
  {
    const double phi = (1.0 + ::sqrt(5))/2.0;
    const double s = R / ::sqrt(phi + 2);
    
    *out ++ = Point(0, +s, +s*phi);
    *out ++ = Point(0, -s, +s*phi);
    *out ++ = Point(0, +s, -s*phi);
    *out ++ = Point(0, -s, -s*phi);
    
    *out ++ = Point(+s, +s*phi, 0);
    *out ++ = Point(+s, -s*phi, 0);
    *out ++ = Point(-s, +s*phi, 0);
    *out ++ = Point(-s, -s*phi, 0);
    
    *out ++ = Point(+s*phi, 0, +s);
    *out ++ = Point(-s*phi, 0, +s);
    *out ++ = Point(+s*phi, 0, -s);
    *out ++ = Point(-s*phi, 0, -s);
  }
  
  template <class Point, class OutputIterator>
  void
  make_dodecahedron (double R,
		     OutputIterator out)
  {
    const double phi = (1.0 + ::sqrt(5))/2.0;
    const double one_phi = 1.0/phi;
    const double s = R / ::sqrt(3.0);

    *out ++ = Point(+s, +s, +s);
    *out ++ = Point(-s, +s, +s);
    *out ++ = Point(+s, -s, +s);
    *out ++ = Point(-s, -s, +s);
    *out ++ = Point(+s, +s, -s);
    *out ++ = Point(-s, +s, -s);
    *out ++ = Point(+s, -s, -s);
    *out ++ = Point(-s, -s, -s);
    
    *out ++ = Point(0, +s*one_phi, +s*phi);
    *out ++ = Point(0, -s*one_phi, +s*phi);
    *out ++ = Point(0, +s*one_phi, -s*phi);
    *out ++ = Point(0, -s*one_phi, -s*phi);
    
    *out ++ = Point(+s*one_phi, +s*phi, 0);
    *out ++ = Point(-s*one_phi, +s*phi, 0);
    *out ++ = Point(+s*one_phi, -s*phi, 0);
    *out ++ = Point(-s*one_phi, -s*phi, 0);
    
    *out ++ = Point(+s*phi, 0, +s*one_phi);
    *out ++ = Point(-s*phi, 0, +s*one_phi);
    *out ++ = Point(+s*phi, 0, -s*one_phi);
    *out ++ = Point(-s*phi, 0, -s*one_phi);
  }

  template <class Vector>
  void
  perturb_points(Vector &p, double eps)
  {
    typedef typename Vector::value_type Point;
    typename CGAL::Random_points_in_sphere_3<Point> r;

    for (size_t i = 0; i < p.size(); ++i)
      p[i] = p[i] + eps * (*r++ - CGAL::ORIGIN);
  }
}

template <class DT>
void
make_icosahedron_dt(DT &dt, double R,
		    double perturbation = 1e-3)
{
  typedef typename DT::Point Point;
  typename std::vector<typename DT::Point> p;

  internal::make_icosahedron<Point> (2*R, std::back_inserter(p));
  internal::perturb_points(p, perturbation);
  dt.insert(p.begin(), p.end());
}

template <class DT>
void
make_dodecahedron_dt(DT &dt, double R,
		     double perturbation = 1e-3)
{
  typedef typename DT::Point Point;
  typename std::vector<typename DT::Point> p;

  internal::make_dodecahedron<Point> (2*R, std::back_inserter(p));
  internal::perturb_points(p, perturbation);
  dt.insert(p.begin(), p.end());
}


// void
// make_truncated_icosahedron_dt(DT &dt, double R,
// 			      double perturbation = 1e-3)
// {
//   typedef typename DT::Point Point;
  
//   typename std::vector<Point> p;
//   internal::make_icosahedron<Point>  (2*R, p);
//   internal::make_dodecahedron<Point> (2*R, p);

//   typename CGAL::Random_points_in_sphere_3<Point> r;
//   for (size_t i = 0; i < p.size(); ++i)
//     p[i] = p[i] + perturbation * (*r++ - CGAL::ORIGIN);

//   dt.insert(p.begin(), p.end());
// }


CGAL_END_NAMESPACE

#endif

