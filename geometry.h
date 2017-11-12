#ifndef GEOMETRY_H_INCLUDED
#define GEOMETRY_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

//	Bounding box [-L_box,L_box]*[-L_box,L_box]
const double L_box = 0.5;
//	precision for geometric computing
const double precision = 1e-8;

class Point
{
public:
	double x, y;
	Point();
	Point(double x_, double y_);
	void print(std::ofstream &file);
	Point operator+(const Point &p_);
	Point operator-(const Point &p_);
	//	cross and dot products
	double cross(const Point &p_);
	double dot(const Point &p_);
	double norm();
};

class Polygon
{
private:
	std::vector<Point> vertices;
public:
	Polygon();
	void add_point(const Point &p);
	void print(std::ofstream &file);
	Point centroid();
};

class Circle
{
public:
    Point center;
    double r;
    //  Given center and radius
    Circle(Point _center, double _r);
    //  Given a,b,c not on a line,
    //  return circle going through these points
    Circle(Point a, Point b, Point c);
    //  Bottom point of the circle
    //  only return the y-coordinate
    double bottom();
};

//  compute the parabola intersection
//  parabola 1: d(r,point p)=d(r,line y=h)
//  parabola 2: d(r,point q)=d(r,line y=h)
//  only returns x-coordinate of the intersection point
//  assume -----------
double parabola_intersect(Point p, Point q, double h);

//  cross square--------
void cross_square(Point st, Point d, int &side, Point &p);

#endif // GEOMETRY_H_INCLUDED
