voronoi_covariance
==================

The purpose of this software is to compute the Voronoi covariance measure of a point set, a notion introduced
for the purpose of geometric inference from point cloud data in the following article 

/Voronoi-based curvature and feature estimation from point clouds,/
Quentin Mérigot, Maks Ovsjanikov, and Leonidas Guibas,
IEEE Transactions on Visualization and Computer Graphics


Build
=====

You need a decently recent version of CGAL.

cd .. 
mkdir voronoi_covariance-build
cd voronoi_covariance-build
ccmake ../voronoi_covariance
make

Usage
=====

The executable vcm takes as input a 3D point cloud. This is a text file with one line per point, with the three
coordinates of the point. The VCM is output as a text file containing 6 floats per line, corresponding to the entries of a 3x3 symmetric matrix. This symmetric is the convolved (r,R)-VCM of the corresponding point in the cloud.

Example:
./vcm --input=fandisk.cloud --R=0.1 --r=0.01 # R is the offset radius, r is the convolution radius

