#include "geometry.h"
using namespace  std;

Point::Point()
{
	x = 0;
	y = 0;
}

Point::Point(double x_, double y_)
{
	x = x_;
	y = y_;
}

void Point::print(ofstream &file)
{
	file << x << " " << y << endl;
}

Point Point::operator+(const Point &p_)
{
	return Point(x + p_.x, y + p_.y);
}

Point Point::operator-(const Point &p_)
{
	return Point(x - p_.x, y - p_.y);
}

double Point::cross(const Point &p_)
{
	return x * p_.y - y * p_.x;
}

double Point::dot(const Point &p_)
{
	return x * p_.x + y * p_.y;
}

double Point::norm()
{
	return sqrt(x * x + y * y);
}

Polygon::Polygon()
{
	vertices.clear();
}

void Polygon::add_point(const Point &p)
{
	//Points should be added anti-clockwisely
	vertices.push_back(p);
}

void Polygon::print(ofstream &file)
{
	for (auto iter = vertices.begin(); iter != vertices.end(); iter++)
	{
		(*iter).print(file);
	}
}

Point Polygon::centroid()
{
	double cx = 0, cy = 0;
	double S = 0;
	Point p_last = vertices.back();
	for (auto iter = vertices.begin(); iter != vertices.end(); iter++)
	{
		Point p = *iter;
		double S_tri = p_last.cross(p) / 2;
		double cx_tri = (p_last + p).x / 3;
		double cy_tri = (p_last + p).y / 3;
		S += S_tri;
		cx += cx_tri * S_tri;
		cy += cy_tri * S_tri;
		p_last = p;
	}
	return Point(cx / S, cy / S);
}

Circle::Circle(Point _center, double _r)
{
    center = _center;
    r = _r;
}

Circle::Circle(Point a, Point b, Point c)
{
    double a11 = (a-b).x, a12 = (a-b).y;
	double a21 = (a-c).x, a22 = (a-c).y;
	double b1 = a.dot(a) - b.dot(b);
	double b2 = a.dot(a) - c.dot(c);
	// Solve equations Ax=b
	// assume solvable, no coline case----
	double d = a11 * a22 - a12 * a21;
	double sx = (a22*b1 - a12*b2) / d / 2;
	double sy = (a11*b2 - a21*b1) / d / 2;
	//	center of circle
	center = Point(sx, sy);
	//	radius of circle
	r = (center-a).norm();
}

double Circle::bottom()
{
	return center.y - r;
}

double parabola_intersect(Point p, Point q, double h)
{
    double A = p.y - q.y;
    double B = 2 * (p.x*(q.y-h) - q.x*(p.y-h));
    double C = (p.y-h)*q.x*q.x - (q.y-h)*p.x*p.x
             - (p.y-q.y)*(p.y-h)*(q.y-h);
    //  Solve Ax^2+Bx+C=0
    if(abs(A)<precision)
    {
        if(abs(B)<precision){cout<<p.x<<" "<<p.y<<" "<<q.x<<" "<<q.y<<" "<<h<<endl;cout<<"Exception!B=0"<<endl;}//debug
        return -C/B;
    }
    else
    {
        double D = B*B - 4*A*C;
        //----------
        if(abs(D)<precision){cout<<"Exception!D=0"<<endl;}//debug
        //----------
        double x1 = (-B+sqrt(D)) / (2*A);
        double x2 = (-B-sqrt(D)) / (2*A);
        //  -----------------
        if(A>0)return min(x1,x2);
        else return max(x1,x2);
    }
}

void cross_square(Point st, Point d, int &side, Point &p)
{
    //  side 0
    if(d.y>precision)
    {
        double u = d.x/d.y*(L_box-st.y)+st.x;
        if(u>=-L_box-precision&&u<=L_box+precision)
        {
            if(abs(u+L_box)<=precision)side = 1;
            else side = 0;
            p = Point(u,L_box);
            return;
        }
    }
    //  side 2
    if(d.y<-precision)
    {
        double u = d.x/d.y*(-L_box-st.y)+st.x;
        if(u>=-L_box-precision&&u<=L_box+precision)
        {
            if(abs(u-L_box)<=precision)side = 3;
            else side = 2;
            p = Point(u, -L_box);
            return;
        }
    }
    //  side 1
    if(d.x<-precision)
    {
        double u = d.y/d.x*(-L_box-st.x)+st.y;
        if(u>=-L_box-precision&&u<=L_box+precision)
        {
            if(abs(u+L_box)<=precision)side = 2;
            else side = 1;
            p = Point(-L_box, u);
            return;
        }
    }
    //  side 3
    if(d.x>precision)
    {
        double u = d.y/d.x*(L_box-st.x)+st.y;
        if(u>=-L_box-precision&&u<=L_box+precision)
        {
            if(abs(u-L_box)<=precision)side = 0;
            else side = 3;
            p = Point(L_box, u);
            return;
        }
    }
}
