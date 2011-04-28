#ifndef CGAL_VORONOI_COVARIANCE_SPHERE_3_HPP
#define CGAL_VORONOI_COVARIANCE_SPHERE_3_HPP

#include <CGAL/point_generators_3.h>

CGAL_BEGIN_NAMESPACE

namespace internal {
  template <class Point>
  double get_coulomb_energy (const std::vector<Point> &p)
  {
    double e = 0;

    size_t N = p.size();
    for(size_t i = 0; i != N; ++i)
      for(size_t j = i+1; j != N; ++j )
	{
	  e += 1 / CGAL::squared_distance(p[i], p[j]);
	}

    return e;
  }

  template <class Point, class Vector>
  void get_forces (const std::vector<Point> &p,
		   std::vector<Vector> &f)
  {
    size_t N = p.size();

    f.resize(N);
    for (size_t i = 0; i < N; i++)
      f[i] = Vector(0,0,0);
  
    double rr, ff, l;
    for (size_t i = 0; i < N; i++)
      {
	for(size_t j = i+1; j<N; j++ )
	  {
	    Vector r = p[i] - p[j];
	    double ll = ::sqrt(r.squared_length());
	    double l = 1.0/(ll*ll*ll);
	    Vector ff = l * r;
	    f[i] = f[i] +  ff;
	    f[j] = f[i] - ff;
	  }
      }

    for (size_t i = 0; i < N; i++)
      {
	const Vector pp = (p[i] - CGAL::ORIGIN);
	f[i] = f[i] - (f[i] * pp) * pp;
      }

    return;
  }

  template <class K>
  void
  make_unit_sphere (std::vector<typename K::Point_3> &p0,
		    size_t N = 30,
		    size_t Nsteps = 1000,
		    double minimal_step = 1e-10,
		    double perturbation = 1e-3)
  {
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typename CGAL::Random_points_in_sphere_3<Point> r;

    p0.resize (N);
    for (size_t i = 0; i < N; ++i)
      p0[i] = *r++;

    double e0 = get_coulomb_energy(p0);
    double step = 0.01;

    for(size_t k = 0; k < Nsteps; ++k)
      {
	double t = exp(- 20.0 * double(k+3)/double(Nsteps));

	std::vector<Vector> f;
	get_forces(p0, f);

	std::vector<Point> p1(N);
	for (size_t i = 0; i < N; ++i)
	  {
	    Vector r = (p0[i] - CGAL::ORIGIN) + step * f[i];
	    typename K::FT rl = 1.0 / ::sqrt(r.squared_length());
	    p1[i] = CGAL::ORIGIN + (rl * r);
	  }
	double e = get_coulomb_energy(p1);

	if(e >= (1+t) * e0)
	  {
	    step /= 2;
	    if (step < minimal_step)
	      break;
	    continue;
	  }
	else
	  {
	    p0 = p1;
	    e0 = e;
	    step *= 2;
	  }
      }
    
    for (size_t i = 0; i < N; ++i)
      p0[i] = p0[i] + perturbation * (*r++ - CGAL::ORIGIN);
  }

}

template <class DT>
void
make_sphere_dt (DT &sphere, 
		double R,
		size_t N = 30,
		size_t Nsteps = 1000,
		double minimal_step = 1e-10)
{
  typedef typename DT::Geom_traits::Kernel K;
  std::vector <typename DT::Point> p0;
  
  internal::make_unit_sphere<K>(p0, N, Nsteps, minimal_step);
  
  for (size_t i = 0; i < N; ++i)
    {
      p0[i] = typename DT::Point(R * p0[i].x(),
				 R * p0[i].y(),
				 R * p0[i].z());
    }
  
  sphere.insert(p0.begin(), p0.end());
}

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

