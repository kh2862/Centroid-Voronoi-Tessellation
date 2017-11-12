#ifndef VORONOI_TESSELLATION_H_INCLUDED
#define VORONOI_TESSELLATION_H_INCLUDED

#include "geometry.h"
#include <vector>
#include <queue>
#include <deque>
#include <string>
#include <fstream>
#include <cmath>
#include <ext/rope>

using namespace std;

class Event;
class Arc;
struct LessThanEvent
{
    bool operator()(Event *lhs, Event *rhs) const;
};

//	Data Structure of the Voronoi Subdivision

class HalfEdge;

class Site
{
public:
	Point pos;
	HalfEdge *edge;
	Site(Point p);
	Polygon make_polygon();
};

class Vertex
{
public:
	Point pos;
	HalfEdge *edge;
	Vertex(Point p);
};

class HalfEdge
{
public:
	Vertex *orig, *dest;
	Site *site;
	HalfEdge *prev, *next, *twin;
	HalfEdge();
};

class VoronoiSubdivision
{
private:
	vector<Site*> sites;
	vector<Vertex*> vertices;
	vector<HalfEdge*> edges;
	__gnu_cxx::rope<Arc*> beachline;
	priority_queue<Event*, vector<Event*>, LessThanEvent> Q;
public:
	VoronoiSubdivision();
	~VoronoiSubdivision();
	void get_sites(vector<Point> &v);
	void do_voronoi_tessellation();
	void region_centroids(vector<Point> &v);
	void print(std::ofstream &file);
	friend class SiteEvent;
	friend class CircleEvent;
};


//	Data Structure of event and arc

class Event
{
protected:
	double height;
public:
	double get_height();
	virtual string event_type() = 0;
    virtual void action(VoronoiSubdivision *vs) = 0;
    virtual ~Event(){}
};

class SiteEvent : public Event
{
private:
	Site *site;
public:
	SiteEvent(Site *_site);
	virtual string event_type();
	virtual void action(VoronoiSubdivision *vs);
	virtual ~SiteEvent(){}
	void insert_to_arc(int ind, VoronoiSubdivision *vs);
};

class CircleEvent : public Event
{
public:
    Arc *arc;
	Site *sitel, *sitem, *siter;
	bool false_alarm;

    Point center;
	CircleEvent(Arc *_arc, Site *_sitel, Site *_sitem, Site *_siter);
	virtual string event_type();
    virtual void action(VoronoiSubdivision *vs);
    virtual ~CircleEvent(){}
    void set_false_alarm();
};

class Arc
{
public:
	Site *site;
	CircleEvent* circ_event;
	HalfEdge *left_edge, *right_edge;
	Arc *leftarc, *rightarc;
	Arc(Site *_site);
	Arc();
};


typedef priority_queue<Event*, vector<Event*>, LessThanEvent> EventQueue;
double calc_right(__gnu_cxx::rope<Arc*> &beachline, int ind, double h);
int bin_search(__gnu_cxx::rope<Arc*> &beachline, double x, double h);
void insert_circle_event(Site *sitel, Site *sitem, Site *siter,
                          Arc *arc, double h, EventQueue &Q);

#endif // VORONOI_TESSELLATION_H_INCLUDED
