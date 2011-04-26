#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>

#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <iostream>
#include <fstream>
#include <cassert>
#include <iterator>
#include <list>
#include <vector>
#include "sphere.hpp"

namespace ublas = boost::numeric::ublas;
typedef double scalar;
typedef ublas::range range;
typedef ublas::vector<scalar> uvector;
typedef ublas::matrix<scalar, ublas::row_major> umatrix;

template <class Matrix>
class Covariance_accumulator
{
public:
  typedef Matrix result_type;
  result_type _result;
  size_t _ntri;
  
public:
  Covariance_accumulator() : 
    _result(3, 3, 0.0),
    _ntri(0)
  {}
    
  template <class Vector>
  inline void operator () (const Vector &a,
			   const Vector &b,
			   const Vector &c)
  {
    // assume center is the origin
    const double m11 = a[0], m12 = b[0], m13 = c[0];
    const double m21 = a[1], m22 = b[1], m23 = c[1];
    const double m31 = a[2], m32 = b[2], m33 = c[2];

    // std::cerr << "a = " << a << "\n";
    // std::cerr << "b = " << b << "\n";
    // std::cerr << "c = " << c << "\n";
    
    const double det = (m11*m33*m22 - m11*m32*m23 - m21*m12*m33 +
			m21*m13*m32 + m31*m12*m23 - m31*m13*m22);
    const double det60 = fabs(det)/60.0;
    
    umatrix R(3,3);
    R(0,0) = (m11*m11 + m11*m12 + m11*m13 +
	      m12*m12 + m12*m13 + m13*m13) * det60;
    R(0,1) = (m11*m21 + m11*m22/2.0 + m11*m23/2.0 +
	      m12*m21/2.0 + m12*m22 + m12*m23/2.0 +
	      m13*m21/2.0 + m13*m22/2.0 + m13*m23) * det60;
    R(0,2) = (m11*m31 + m11*m32/2.0 + m11*m33/2.0 + 
	      m12*m31/2.0 + m12*m32 + m12*m33/2.0 + 
	      m13*m31/2.0 + m13*m32/2.0 + m13*m33) * det60;
    R(1,0) = R(0,1);
    R(1,1) = (m21*m21 + m21*m22 + m21*m23 + 
	      m22*m22 + m22*m23 + m23*m23) * det60;
    R(1,2) = (m31*m21 + m31*m22/2.0 + m31*m23/2.0 +
	      m32*m21/2.0 + m32*m22 + m32*m23/2.0 +
	      m33*m21/2.0 + m33*m22/2.0 + m33*m23) * det60;
    R(2,0) = R(0,2);
    R(2,1) = R(1,2);
    R(2,2) = (m31*m31 + m31*m32 + m31*m33 +
	      m32*m32 + m32*m33 + m33*m33) * det60;
    
    _result += R;
    _ntri++;
  }
  
  const result_type &result() const
  {
    //std::cerr << "ntri = " << _ntri << "\n";
    return _result;
  }
};


template <class DT, class F>
F &tessellate (const DT &dt,
	       typename DT::Vertex_handle v,
	       F &f)
{
  typedef typename DT::Point Point;
  typedef typename DT::Geom_traits::Kernel K;
  typedef typename K::Vector_3 Vector;
  typedef typename DT::Cell_handle Cell_handle;
  typedef typename DT::Vertex_handle Vertex_handle;
  
  const Point A (v->point());

  // get all vertices incident to v
  typename std::list<Vertex_handle> vertices;

  dt.incident_vertices(v,std::back_inserter(vertices));
  for(typename std::list<Vertex_handle>::iterator it = vertices.begin();
      it != vertices.end(); it++)
    {
      // build edge from two vertices
      typename DT::Cell_handle cell;
      int i1,i2;
      
      if(!dt.is_edge(v, *it, cell, i1, i2))
	continue;
      
      // tessellate the polygon around its center
      typename DT::Cell_circulator c = dt.incident_cells(cell, i1, i2);
      typename DT::Cell_circulator begin = c; 
      typename DT::Cell_circulator done = c; 
      
      Vector C(dt.dual(c) - CGAL::ORIGIN);  c++;
      size_t count = 1;
      while (c != done)
	{	   
	  C = C + (dt.dual(c) - CGAL::ORIGIN);
	  count++;
	  c++;
	}
      Point center = CGAL::ORIGIN + (1.0/count) * C;
      
      c = begin;	
      const Point u = dt.dual(c); c++;
      const Point v = dt.dual(c);
      f(center-A, u-A, v-A);
      
      while (c != done)
	{
	  const Point u = dt.dual(c); c++;
	  const Point v = dt.dual(c);
	  f(center-A, u-A, v-A);
	}
    }
  return f;
}


template <class DT, class F>
F& tessellate_and_intersect(const DT &model,
			    typename DT::Vertex_handle v,
			    const DT &sphere, 
			    F &f)
{
  typedef typename DT::Vertex_handle Vertex_handle;
  typedef typename DT::Point Point;

  DT local = sphere;

  std::list<Vertex_handle> vertices;
  model.incident_vertices(v,std::back_inserter(vertices));
      
  std::vector<Point> points;
  for(typename std::list<Vertex_handle>::iterator
	it = vertices.begin();
        it != vertices.end(); ++it)
    {
      points.push_back( CGAL::ORIGIN + ((*it)->point() - v->point()) );
    }
  points.push_back(CGAL::ORIGIN);

  //std::cerr << points.size() << " Delaunay neighbors\n";

  local.insert(points.begin(), points.end());
  Vertex_handle nv = local.nearest_vertex(CGAL::ORIGIN);

  return tessellate (local, nv, f);
}

template <class DT>
umatrix
voronoi_covariance_3 (const DT &dt,
		      typename DT::Vertex_handle v,
		      const DT &sphere)
{
  Covariance_accumulator<umatrix> accum;
  return tessellate_and_intersect(dt, v, sphere, accum).result();
}

template <class DT>
void load_xyz(DT &dt, const std::string &s,
	      std::vector<typename DT::Vertex_handle> &vh)
{
  typedef typename DT::Point Point;

  std::ifstream is (s.c_str());
  double x, y, z;

  typename std::vector<Point> points;
  while (is >> x >> y >> z)
    points.push_back(Point(x,y,z));

  dt.insert(points.begin(), points.end());

  vh.clear();
  for (size_t i = 0; i < points.size(); ++i)
    vh.push_back(dt.nearest_vertex(points[i]));
}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K> DT;

std::ostream &
operator << (std::ostream &os, const umatrix &m)
{
  os << m(0,0) << " " << m(1,0) << " " << m(2, 0) << " "
     << m(1,1) << " " << m(2,1) << " "
     << m(2,2);
  return os;
}

template <class DT>
void
make_sphere (DT &sphere, 
             double R,
             size_t N = 30)
{
  typedef typename DT::Point Point;
  std::vector<Point> points;

#if 0
  CGAL::Random_points_in_sphere_3<Point> r (R);

  for (size_t i = 0; i < N; ++i)
    points.push_back(*r++);
  sphere.insert(points.begin(), points.end());
#elif 0
  static const double sampling[][3] = 
    {
      { 0.0863709168,    0.1854629722,    0.9788480733},
      { 0.1752854391,   -0.8379426457,   -0.5168434360},
      {-0.5623793283,    0.2791448137,   -0.7783364722},
      { 0.7585690449,    0.5989658402,   -0.2565403018},
      { 0.1607583477,    0.5927997765,   -0.7891420523},
      {-0.6979759978,    0.7044759450,   -0.1286201751},
      {-0.5526722392,   -0.5508598110,   -0.6253853729},
      { 0.7883873022,   -0.1027831797,   -0.6065320104},
      { 0.0716498794,   -0.9639822164,    0.2561339126},
      {-0.7885431056,    0.1022925191,    0.6064124102},
      { 0.9451052504,    0.0737868052,    0.3183262053},
      {-0.3178613514,   -0.4869599454,    0.8135319126},
      {-0.9898952969,   -0.0422659821,   -0.1353546748},
      { 0.1046132024,   -0.1846039662,   -0.9772294785},
      { 0.4959735376,    0.6891424646,    0.5282924506},
      {-0.6724384833,   -0.7276983541,    0.1352094363},
      { 0.5135779664,   -0.4328307126,    0.7408746498},
      { 0.7766471502,   -0.6292979298,   -0.0283428929},
      {-0.3330976774,    0.7409342586,    0.5831486618},
      { 0.0379306227,    0.9922356610,   -0.1184468693},
    };
  static size_t Nsampling = sizeof(sampling)/sizeof(sampling[0]);
  for (size_t i = 0; i < Nsampling; ++i)
    {
      points.push_back(Point (R*sampling[i][0],
			      R*sampling[i][1],
			      R*sampling[i][2]));
      std::cerr << sampling[i][0] << "\t"
		<< sampling[i][1] << "\t" 
		<< sampling[i][2] << "\n";
      double d = (pow(sampling[i][0], 2.0) + 
		  pow(sampling[i][1], 2.0) +
		  pow(sampling[i][2], 2.0));
      std::cerr << d << "\n";
    }
#else 
  static const double sampling[][3] = 
    {
{ -0.122809, -0.157329, -0.012865 },
{ -0.161661, -0.0139973, 0.116918 },
{ -0.145037, 0.00473613, -0.137629 },
{ 0.0295226, 0.156518, -0.120957 },
{ 0.145654, 0.00783818, 0.136834 },
{ 0.110176, -0.108739, -0.126637 },
{ -0.13169, 0.111213, 0.101437 },
{ 0.02939, -0.120378, 0.156989 },
{ 0.12702, -0.125138, 0.09059 },
{ 0.141686, 0.0796305, -0.116551 },
{ -0.0950974, -0.129045, 0.119599 },
{ -0.106145, 0.130281, -0.108444 },
{ -0.103153, -0.113623, -0.128255 },
{ 0.103153, 0.113623, 0.128255 },
{ 0.0126554, 0.0102644, -0.199335 },
{ -0.0126554, -0.0102644, 0.199335 },
{ -0.00811274, 0.199597, 0.00976275 },
{ 0.00811274, -0.199597, -0.00976278 },
{ -0.199434, -0.00746801, -0.0130463 },
{ 0.199434, 0.00746801, 0.0130463 }
    };
  static size_t Nsampling = sizeof(sampling)/sizeof(sampling[0]);
  for (size_t i = 0; i < Nsampling; ++i)
    {
      points.push_back(Point (sampling[i][0],
			      sampling[i][1],
			      sampling[i][2]));
    }
  std::cerr << "npoints " << points.size() << "\n";
#endif
  sphere.insert(points.begin(), points.end());
}


int main()
{
  boost::timer t;
  double R = 0.1;
  size_t N = 20;

  DT dt, sphere;
  std::vector<DT::Vertex_handle> vh;

  load_xyz(dt, "fandisk.cloud", vh);
  VCM::make_sphere (sphere, 2*R, N);
  //make_sphere (sphere, 2*R, N);
  std::cerr <<  sphere.number_of_vertices() << "\n";

  std::cerr << "built Delaunay triangulation in " << t.elapsed() << "s\n";
  std::cerr << "number of vertices: " << dt.number_of_vertices() << "\n";  

  boost::progress_display p(vh.size(),
			    std::cerr);
  t.restart();
  for (size_t i = 0; i < vh.size(); ++i)
    {
      umatrix m = voronoi_covariance_3(dt, vh[i], sphere);
      std::cout << m << "\n";
      ++p;
    }  
  std::cerr << "time: " <<  t.elapsed() << "s\n";    

  return 0;
}

