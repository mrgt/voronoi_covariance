#ifndef VCM_SPHERE_HPP
#define VCM_SPHERE_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>

namespace VCM
{
  namespace details
  {
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
	      double ll = sqrt(r.squared_length());
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

    double e0 = details::get_coulomb_energy(p0);
    double step = 0.01;

    for(size_t k = 0; k < Nsteps; ++k)
      {
	double t = exp(- 20.0 * double(k+3)/double(Nsteps));

	std::vector<Vector> f;
	details::get_forces(p0, f);

	std::vector<Point> p1(N);
	for (size_t i = 0; i < N; ++i)
	  {
	    Vector r = (p0[i] - CGAL::ORIGIN) + step * f[i];
	    p1[i] = CGAL::ORIGIN + (r / sqrt(r.squared_length()));
	  }
	double e = details::get_coulomb_energy(p1);

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

  template <class DT>
  void
  make_sphere (DT &sphere, 
	       double R,
	       size_t N = 30,
	       size_t Nsteps = 1000,
	       double minimal_step = 1e-10)
  {
    typedef typename DT::Geom_traits::Kernel K;
    std::vector <typename DT::Point> p0;

    make_unit_sphere<K>(p0, N, Nsteps, minimal_step);

    for (size_t i = 0; i < N; ++i)
      {
	p0[i] = typename DT::Point(R * p0[i].x(),
				   R * p0[i].y(),
				   R * p0[i].z());
	std::cerr << "{ "
		  <<  p0[i].x() << ", "
		  <<  p0[i].y() << ", "
		  <<  p0[i].z()
		  << " }";
	if (i < (N-1))
	  std::cerr << ",";
	std::cerr << "\n";
      }

    //sphere.clear();
    sphere.insert(p0.begin(), p0.end());
  }
}

#endif

