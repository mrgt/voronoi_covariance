#ifndef CGAL_VORONOI_COVARIANCE_3_HPP

#include <list>
#include <CGAL/array.h>
#include <CGAL/voronoi_covariance_sphere_3.h>

CGAL_BEGIN_NAMESPACE

namespace internal {

   template <class FT, class Array>
   inline void
   covariance_matrix_tetrahedron (FT ax, FT ay, FT az,
				  FT bx, FT by, FT bz,
				  FT cx, FT cy, FT cz,
				  Array &R)
   { 
      const FT det = (ax*cz*by - ax*bz*cy - ay*bx*cz +
		      ay*cx*bz + az*bx*cy - az*cx*by);
      const double det60 = fabs(det)/60.0;
	
      R[0] += (ax*ax + ax*bx + ax*cx +
	       bx*bx + bx*cx + cx*cx) * det60;
      R[1] += (ax*ay + ax*by/2.0 + ax*cy/2.0 +
	       bx*ay/2.0 + bx*by + bx*cy/2.0 +
	       cx*ay/2.0 + cx*by/2.0 + cx*cy) * det60;
      R[2] += (ax*az + ax*bz/2.0 + ax*cz/2.0 + 
	       bx*az/2.0 + bx*bz + bx*cz/2.0 + 
	       cx*az/2.0 + cx*bz/2.0 + cx*cz) * det60;
      R[3] += (ay*ay + ay*by + ay*cy + 
	       by*by + by*cy + cy*cy) * det60;
      R[4] += (az*ay + az*by/2.0 + az*cy/2.0 +
	       bz*ay/2.0 + bz*by + bz*cy/2.0 +
	       cz*ay/2.0 + cz*by/2.0 + cz*cy) * det60;
      R[5] += (az*az + az*bz + az*cz +
	       bz*bz + bz*cz + cz*cz) * det60;
   }
   

   template <class FT>
   class Covariance_accumulator_3
   {
   public:
      typedef array<FT, 6> Result_type;
      
   private:
      Result_type _result;
      size_t _ntri;
      
   public:
      Covariance_accumulator_3() : _ntri(0)
      {
	 std::fill (_result.begin(), _result.end(), FT(0));
      }
      
      template <class Vector>
      inline void operator () (const Vector &a,
			       const Vector &b,
			       const Vector &c)
      {
	internal::covariance_matrix_tetrahedron (a[0], a[1], a[2],
						 b[0], b[1], b[2],
						 c[0], c[1], c[2],
						 _result);
      }
       
      const Result_type &result() const
      {
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
   F& tessellate_and_intersect(const DT &dt,
			       typename DT::Vertex_handle v,
			       const DT &sphere, 
			       F &f)
   {
      typedef typename DT::Vertex_handle Vertex_handle;
      typedef typename DT::Point Point;
	
      DT local = sphere;
	
      std::list<Vertex_handle> vertices;
      dt.incident_vertices(v,std::back_inserter(vertices));
	
      std::vector<Point> points;
      for(typename std::list<Vertex_handle>::iterator
	     it = vertices.begin();
	  it != vertices.end(); ++it)
      {
	 points.push_back( CGAL::ORIGIN + ((*it)->point() - v->point()) );
      }
      points.push_back(CGAL::ORIGIN);
	
      local.insert(points.begin(), points.end());
      Vertex_handle nv = local.nearest_vertex(CGAL::ORIGIN);
	
      return tessellate (local, nv, f);
   }
}
   
template <class DT, class FT>
void
voronoi_covariance_3 (const DT &dt,
		      typename DT::Vertex_handle v,
		      const DT &sphere,
		      FT covariance[6])
{
  typename internal::Covariance_accumulator_3<FT> ca;
  internal::tessellate_and_intersect(dt, v, sphere, ca);  
  std::copy (ca.result().begin(), ca.result().end(), covariance);
}

template<class FT>
class Voronoi_covariance_3 : public CGAL::array<FT,6>
{
  typedef typename CGAL::array<FT,6> Parent;
  typedef typename Parent::iterator iterator;
  typedef typename Parent::const_iterator const_iterator;

public:
  Voronoi_covariance_3 (const Parent &p) : Parent(p)
  {}

  Voronoi_covariance_3 (FT m[6])
  {
    std::copy (m, m + 6, Parent::begin());
  }
  
  Voronoi_covariance_3 ()
  {
    std::fill(Parent::begin(), Parent::end(), FT(0));
  }

  Voronoi_covariance_3<FT> &
  operator += (const Voronoi_covariance_3<FT>& c)
  {
    for (size_t i = 0; i < 6; ++i)
      Parent::operator[] (i) += c[i];
  }
};

template <class DT>
array<typename DT::Geom_traits::FT, 6>
voronoi_covariance_3 (const DT &dt,
		      typename DT::Vertex_handle v,
		      const DT &sphere)
{
  typedef typename DT::Geom_traits::FT FT;
  typename internal::Covariance_accumulator_3<FT> ca;

  return internal::tessellate_and_intersect(dt, v, sphere, ca).result();
}

template <class FT>
std::ostream &
operator << (std::ostream &os, const CGAL::Voronoi_covariance_3<FT> &cov)
{
  return os << cov[0] << " " << cov[1] << " " << cov[2] << " "
	    << cov[3] << " " << cov[4] << " " 
	    << cov[5] << "\n";
}

template <class FT>
std::istream &
operator >> (std::istream &is, CGAL::Voronoi_covariance_3<FT> &cov)
{
  return  is >> cov[0] >> cov[1] >> cov[2] 
	     >> cov[3] >> cov[4]
	     >> cov[5];
}


/*
template <class DT>
CGAL::array<typename DT::FT,6>
voronoi_covariance_3 (const DT &dt,
		      typename DT::Vertex_handle v,
		      const DT &sphere)
{
  typename internal::Covariance_accumulator_3<typename DT::FT> ca;
  return internal::tessellate_and_intersect(dt, v, sphere, ca).result();
}
*/

CGAL_END_NAMESPACE

#endif
