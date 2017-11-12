#ifndef CENTROID_VORONOI_TESSELLATION_H_INCLUDED
#define CENTROID_VORONOI_TESSELLATION_H_INCLUDED

#include "geometry.h"
#include "Voronoi_Tessellation.h"

//	Model Class
class Model
{
private:
	VoronoiSubdivision *vs;
	vector<Point> v;
public:
	double tol;
	int iter, maxiter;
	Model();
	~Model();
	void input_generator(string filename);
	void output(string filename);
	void Lloyd();
};

#endif // CENTROID_VORONOI_TESSELLATION_H_INCLUDED
