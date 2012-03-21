#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/halfspaces_intersection.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iterator>
#include <fstream>
#include <stdlib.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Convex_hull_traits_3<K> Traits;
typedef Traits::Polyhedron_3 Polyhedron;
typedef Polyhedron::Facet_iterator                   Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Plane_3 Plane;


double mu = 0.4;

Point
project_back(const Point & p)
{
  double sk = p.x() * p.x() + p.y() * p.y() + (p.z() * p.z() / (mu*mu));
  return Point(p.x()/sqrt(sk), p.y()/sqrt(sk), p.z()/sqrt(sk));
}

Plane
tangent_plane (Point &p)
{
  Vector v (p.x(), p.y(), p.z()/(mu*mu));
  v = v/sqrt(v.squared_length());
  return Plane(v.x(), v.y(), v.z(), - (p - CGAL::ORIGIN) * v);
}

 Point
 rescale (Point p)
 {
   return Point(p.x(), p.y(), mu * p.z());
 }


 double lloyd_step (std::vector<Point> &points)
 {
   std::vector<Plane> planes;
   for (size_t i = 0; i < points.size(); ++i)
     planes.push_back(tangent_plane(points[i]));

   Polyhedron P;

   CGAL::internal::halfspaces_intersection(planes.begin(),
					   planes.end(), P, K());

   size_t N = points.size();
   points.clear();
   double totarea = 0.0;
   double vararea = 0.0;

   for (Polyhedron::Facet_iterator it = P.facets_begin();
	it != P.facets_end(); ++it)
     {
       Polyhedron::Halfedge_around_facet_circulator
	 h0 = it->facet_begin(), hf = h0--, hs = hf;
       hs ++;

       double farea = 0.0;
       Vector fcentroid (CGAL::NULL_VECTOR);
       while(1)
	 {
	   const Point &a  = h0->vertex()->point(),
	     &b = hf->vertex()->point(),
	     &c = hs->vertex()->point();
	   Vector v = 
	     CGAL::cross_product(b -a, c - a);
	   double tarea = .5 * sqrt(v.squared_length());
	   farea += tarea;
	   fcentroid = fcentroid +  (tarea/3.0) * ((a - CGAL::ORIGIN) + 
						   (b - CGAL::ORIGIN) + 
						   (c - CGAL::ORIGIN));
	   if (hs == h0)
	     break;
	   ++hs; ++hf;
	 }
       totarea += farea;
       vararea += pow(4.0 * M_PI / N - farea, 2.0);
       fcentroid = fcentroid/farea;
       points.push_back(project_back(CGAL::ORIGIN + fcentroid));
     }
   std::cerr << "area = " << totarea << ", "
	     << "vararea = " << vararea << "\n";
   return vararea/totarea;
 }

 int main(int argc, const char **argv )
 {
   std::vector<Point> points;
   CGAL::Random_points_on_sphere_3<Point> g;

   size_t N = 0;
   if (argc > 1)
     N = atof(argv[1]);
   N = std::max(size_t(100), N);

   for (size_t i = 0; i < N; ++i)
     points.push_back(rescale(*g++));

   for (size_t n = 0; n < 100; ++n)
     {
       std::cerr << "step " << n << ":\n\t";
       lloyd_step(points);
     }

   Polyhedron P;
   CGAL::convex_hull_3(points.begin(), points.end(), P);

   CGAL::set_ascii_mode( std::cout);
   std::cout << "OFF" << std::endl << P.size_of_vertices() << ' '
	     << P.size_of_facets() << " 0" << std::endl;
   std::copy( P.points_begin(), P.points_end(),
	      std::ostream_iterator<Point>( std::cout, "\n"));
   for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
     Halfedge_facet_circulator j = i->facet_begin();
     // Facets in polyhedral surfaces are at least triangles.
     CGAL_assertion( CGAL::circulator_size(j) >= 3);
     std::cout << CGAL::circulator_size(j) << ' ';
     do {
       std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());
     } while ( ++j != i->facet_begin());
     std::cout << std::endl;
   }

   std::ofstream os ("test.cloud");
   std::copy(points.begin(), points.end(),
	     std::ostream_iterator<Point>(os, "\n"));
}
