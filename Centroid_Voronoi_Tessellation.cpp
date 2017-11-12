#include "geometry.h"
#include "Voronoi_Tessellation.h"
#include "Centroid_Voronoi_Tessellation.h"
using namespace std;

Model::Model()
{
	vs = NULL;
	v.clear();
	iter = 0;
	//	default setting
	tol = 1e-6;
	maxiter = 500;
}

Model::~Model()
{
	delete vs;
}

void Model::input_generator(string filename)
{
	ifstream file;
	file.open(filename);
	int n;
	file >> n;
	for (int i = 0; i < n; i++)
	{
		double x, y;
		file >> x >> y;
		v.push_back(Point(x, y));
	}
	file.close();
}

void Model::output(string filename)
{
	ofstream file;
	file.open(filename);
	vs->print(file);
	file.close();
}

void Model::Lloyd()
{
	//	v has been initialized by input_generator
	iter = 0;
	while(true)
	{
		iter++;
		vs = new VoronoiSubdivision;
		vs->get_sites(v);
		vs->do_voronoi_tessellation();
		vector<Point> v_new;
		v_new.clear();
		vs->region_centroids(v_new);
		//	calculate change of generators(L2 norm)
		double change = 0, vabs = 0;
		int n = v.size();
		for (int i = 0; i < n; i++)
		{
			change += (v[i] - v_new[i]).dot(v[i] - v_new[i]);
			vabs += v[i].dot(v[i]);
		}
		change = sqrt(change);
		vabs = sqrt(vabs);
		//	if the change or relative change is less than
		//	the tolerance or iter>maxiter then exit
		if (change < tol || change / vabs < tol || iter > maxiter)
		{
			return;
		}
		else
		{
			v = v_new;
			delete vs;
		}
	}
}
