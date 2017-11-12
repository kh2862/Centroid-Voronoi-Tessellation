#include "geometry.h"
#include "Voronoi_Tessellation.h"
using namespace std;
using namespace __gnu_cxx;

Site::Site(Point p)
{
	pos = p;
	edge = NULL;
}

Polygon Site::make_polygon()
{
	deque<Point> pts;
	pts.clear();
	bool attainbound = true;
	HalfEdge *ed = NULL;
	for (ed = edge; ed->next != NULL; ed = ed->next)
	{
		pts.push_back(ed->dest->pos);
		if(ed == edge)
        {
            attainbound = false;
            break;
        }
	}
	if(attainbound)
    {
        Point pl = pos;
        Point pr = ed->twin->site->pos;
        int sidebg;
        Point bg;
        cross_square(Point((pl+pr).x/2, (pl+pr).y/2), Point(-(pr-pl).y, (pr-pl).x), sidebg, bg);
        for(ed = edge; ed->prev != NULL; ed = ed->prev)
        {
            pts.push_front(ed->orig->pos);
        }
        pl = pos;
        pr = ed->twin->site->pos;
        int sidefn;
        Point fn;
        cross_square(Point((pl+pr).x/2, (pl+pr).y/2), Point((pr-pl).y, -(pr-pl).x), sidefn, fn);
        if(abs(bg.x-fn.x)<precision&&abs(bg.y-fn.y)<precision)
        {
            pts.push_back(bg);
        }
        else
        {
            //  wrap along square
            pts.push_back(bg);
            pts.push_front(fn);
            while(sidebg != sidefn)
            {
                switch(sidebg)
                {
                    case 0: bg = Point(-L_box,  L_box); break;
                    case 1: bg = Point(-L_box, -L_box); break;
                    case 2: bg = Point( L_box, -L_box); break;
                    case 3: bg = Point( L_box,  L_box); break;
                }
                pts.push_back(bg);
                sidebg = (sidebg + 1) % 4;
            }
        }
    }
    //  make polygon, ignore coincide points
    Polygon poly;
    int n = pts.size();
    while(abs(pts[n-1].x-pts[0].x)<precision&&abs(pts[n-1].y-pts[0].y)<precision)
    {
        n--;
    }
    poly.add_point(pts[0]);
    for(int i = 1 ; i < n; i++)
    {
        if(abs(pts[i].x-pts[i-1].x)>precision ||
           abs(pts[i].y-pts[i-1].y)>precision)
        {
            poly.add_point(pts[i]);
        }
    }
	return poly;
}


Vertex::Vertex(Point p)
{
	pos = p;
	edge = NULL;
}


HalfEdge::HalfEdge()
{
	orig = NULL;
	dest = NULL;
	site = NULL;
	prev = NULL;
	next = NULL;
	twin = NULL;
}


VoronoiSubdivision::VoronoiSubdivision()
{
	sites.clear();
	vertices.clear();
	edges.clear();
	beachline.clear();
	while(!Q.empty())Q.pop();
}

VoronoiSubdivision::~VoronoiSubdivision()
{
	for (auto iter = sites.begin(); iter != sites.end(); iter++)
	{
		Site *site = *iter;
		delete site;
	}
	for (auto iter = vertices.begin(); iter != vertices.end(); iter++)
	{
		Vertex *vertex = *iter;
		delete vertex;
	}
	for (auto iter = edges.begin(); iter != edges.end(); iter++)
	{
		HalfEdge *edge = *iter;
		delete edge;
	}
	for(auto iter = beachline.begin(); iter != beachline.end(); iter++)
    {
        Arc *arc = *iter;
        delete arc;
    }
}

void VoronoiSubdivision::get_sites(vector<Point> &v)
{
    for(auto iter = v.begin(); iter != v.end(); iter++)
    {
        Site *site = new Site(*iter);
        sites.push_back(site);
    }
}

//	Fortune's algorithm
void VoronoiSubdivision::do_voronoi_tessellation()
{
    //  initialize
    beachline.clear();
    while(!Q.empty())Q.pop();
    for(auto iter = sites.begin(); iter != sites.end(); iter++)
    {
        Event *e = new SiteEvent(*iter);
        Q.push(e);
    }
    //  handle events
    while(!Q.empty())
    {
        Event *cur = Q.top();
        Q.pop();
        cur->action(this);
        //debug
//        cout << cur->event_type() << " "<<cur->get_height() << endl;
//        for(int i =0;i<beachline.size();i++)
//        {
//            Arc *arc = beachline[i];
//            cout<<arc->site->pos.x<<" "<<arc->site->pos.y<<endl;
//            if(arc->circ_event!=NULL)cout<<arc->circ_event->sitem->pos.x<<" "<<arc->circ_event->sitem->pos.y<<endl;
//        }
        //debug
        delete cur;
    }
}

void VoronoiSubdivision::region_centroids(vector<Point> &v)
{	//	assume v is cleared
    //debug
//    ofstream f;f.open("tmp.txt");
//    f<<"L"<<endl;
//    for(auto iter = vertices.begin(); iter!=vertices.end();iter++)
//    {
//        f <<(*iter)->pos.x<<(*iter)->pos.y<<endl;
//    }
	for (auto iter = sites.begin(); iter != sites.end(); iter++)
	{
		Site *site = *iter;
		Polygon poly = site->make_polygon();
		//poly.print(f);
		Point p = poly.centroid();
		//f<<"J "<<p.x<<" "<<p.y<<endl;
		v.push_back(p);
	}
	//f.close();
}

void VoronoiSubdivision::print(ofstream &file)
{
    cout << sites.size() << endl;
	for (auto iter = sites.begin(); iter != sites.end(); iter++)
	{
		Site *site = *iter;
		Polygon poly = site->make_polygon();
		poly.print(file);
	}
}


double Event::get_height()
{
	return height;
}

bool LessThanEvent::operator()(Event *lhs, Event *rhs) const
{
    return lhs->get_height() < rhs->get_height();
}

SiteEvent::SiteEvent(Site *_site)
{
	site = _site;
	height = site->pos.y;
}

string SiteEvent::event_type()
{
	return "Site Event";
}

double calc_right(rope<Arc*> &beachline, int ind, double h)
{
    int n = beachline.size();
    //debug
    //cout<<ind<<" "<<h<<" "<<n<<endl;
    //debug
    if(ind==n-1)
    {
        return(1.0);
    }
    else
    {
        Point pl = beachline[ind]->site->pos;
        Point pr = beachline[ind+1]->site->pos;
        return parabola_intersect(pl, pr, h);
    }
}

int bin_search(rope<Arc*> &beachline, double x, double h)
{
    int left = 0, right = beachline.size()-1;
    while(left<=right)
    {
        int mid = (left+right) / 2;
        double xn = calc_right(beachline, mid, h);
        if(xn>=x-precision)right=mid-1;
        else left=mid+1;
    }
    return left;
}

void SiteEvent::action(VoronoiSubdivision *vs)
{
    //  highest sites should be specially treated-----
    if(vs->beachline.empty())
    {
        Arc *arc = new Arc(this->site);
        vs->beachline.push_back(arc);
        return;
    }

    int ind = bin_search(vs->beachline, site->pos.x, height);
    insert_to_arc(ind, vs);
    double x = calc_right(vs->beachline, ind+2, height);
    if(abs(x-site->pos.x)<precision)
    {
        //  ind+1 position should be treated----
        insert_to_arc(ind+3, vs);
    }
}

void SiteEvent::insert_to_arc(int ind, VoronoiSubdivision *vs)
{
    Arc *arc = vs->beachline[ind];
    if(arc->circ_event!=NULL)
    {
        arc->circ_event->set_false_alarm();
    }
    HalfEdge *edge1 = new HalfEdge();
    HalfEdge *edge2 = new HalfEdge();
    edge1->site = this->site;
    edge2->site = arc->site;
    edge1->twin = edge2;
    edge2->twin = edge1;
    this->site->edge = edge1;
    arc->site->edge = edge2;
    //-----------
    Arc *arcn2 = new Arc(*arc);
    Arc *arcn = new Arc(this->site);
    arcn->left_edge = edge1;
    arcn->right_edge = edge1;
    arc->left_edge = edge2;
    arcn2->right_edge = edge2;
    arcn->leftarc = arcn2;
    arcn->rightarc = arc;
    arc->leftarc = arcn;
    arcn2->rightarc = arcn;
    if(arcn2->leftarc!=NULL)arcn2->leftarc->rightarc = arcn2;
    vs->beachline.insert(ind, arcn);
    vs->beachline.insert(ind, arcn2);
    // debug
    //cout<<"sddsdsdds "<<arc.site->pos.x<<" "<<arc.site->pos.y<<endl;
    //cout<<"sddsdsdds "<<this->site->pos.x<<" "<<this->site->pos.y<<endl;
    // debug
    //-------
    if(ind>0)
    {
        Arc *arcl = vs->beachline[ind-1];
        insert_circle_event(arcl->site, arc->site, this->site,
                            arcn2, this->height, vs->Q);
    }
    int n = vs->beachline.size();
    if(ind<n-3)
    {
        Arc *arcr = vs->beachline[ind+3];
        insert_circle_event(this->site, arc->site, arcr->site,
                            arc, this->height, vs->Q);
    }
}

void insert_circle_event(Site *sitel, Site *sitem, Site *siter,
                        Arc *arc, double h, EventQueue &Q)
{
    Point pl = sitel->pos;
    Point pm = sitem->pos;
    Point pr = siter->pos;
    //  colinear case
    double v = (pl-pm).cross(pr-pm);
    if(abs(v)<precision)return;
    //  ------
    Circle circ(pl, pm, pr);
    double x1, y1, x2, y2;
    if(abs(pl.y-h)>precision)
    {
        x1 = parabola_intersect(pl, pm, h);
        y1 = ((x1-pl.x)*(x1-pl.x)-h*h+pl.y*pl.y)/2/(pl.y-h);
    }
    else
    {

    }
    if(abs(pr.y-h)>precision)
    {
        x2 = parabola_intersect(pm, pr, h);
        y2 = ((x2-pr.x)*(x2-pr.x)-h*h+pr.y*pr.y)/2/(pr.y-h);
    }
    else
    {

    }
    if(circ.bottom()>h+precision||circ.center.y<-L_box-precision
        ||circ.center.y>y1-precision||circ.center.y>y2-precision)
    {
        return;
    }
    else
    {
        //debug
        //cout<<"LOL "<<pl.x<<" "<<pl.y<<" "<<pm.x<<" "<<pm.y<<" "<<pr.x<<" "<<pr.y<<" "<<circ.bottom()<<endl;
        Event *e = new CircleEvent(arc, sitel, sitem, siter);
        Q.push(e);
    }
}


CircleEvent::CircleEvent(Arc *_arc, Site *_sitel, Site *_sitem, Site *_siter)
{
    arc = _arc;
	sitel = _sitel;
	sitem = _sitem;
	siter = _siter;
	Point pl = sitel->pos;
	Point pm = sitem->pos;
	Point pr = siter->pos;
	Circle circ(pl, pm, pr);
	center = circ.center;
	height = circ.bottom();
	false_alarm = false;
}

string CircleEvent::event_type()
{
	return "Circle Event";
}

void CircleEvent::action(VoronoiSubdivision *vs)
{
    if(false_alarm)return;
    CircleEvent *cel = arc->leftarc->circ_event;
    if(cel!=NULL)
    {
        cel->set_false_alarm();
    }
    CircleEvent *cer = arc->rightarc->circ_event;
    if(cer!=NULL)
    {
        cer->set_false_alarm();
    }
    //--------
    Vertex *v = new Vertex(center);
    v->edge = arc->right_edge;
    HalfEdge *edge1 = new HalfEdge();
    HalfEdge *edge2 = new HalfEdge();
    edge1->twin = edge2;
    edge2->twin = edge1;
    edge1->orig = v;
    edge2->dest = v;
    edge2->next = arc->left_edge->twin;
    edge1->prev = arc->right_edge->twin;
    edge1->site = siter;
    edge2->site = sitel;
    arc->left_edge->dest = v;
    arc->left_edge->next = arc->right_edge;
    arc->left_edge->twin->orig = v;
    arc->left_edge->twin->prev = edge2;
    arc->right_edge->orig = v;
    arc->right_edge->prev = arc->left_edge;
    arc->right_edge->twin->dest = v;
    arc->right_edge->twin->next = edge1;
    //---------
    Arc *arcl = arc->leftarc;
    Arc *arcr = arc->rightarc;
    arcl->right_edge = edge2;
    arcr->left_edge = edge1;
    arcl->rightarc = arcr;
    arcr->leftarc = arcl;
    //--------
    if(arcl->leftarc!=NULL)
    {
        Site *sitell = arcl->leftarc->site;
        insert_circle_event(sitell, sitel, siter, arcl,
                            this->height, vs->Q);
    }
    if(arcr->rightarc!=NULL)
    {
        Site *siterr = arcr->rightarc->site;
        insert_circle_event(sitel, siter, siterr, arcr,
                            this->height, vs->Q);
    }
    //--------------
    int n = vs->beachline.size();
    for(int i = 0; i < n; i++)
    {
        if(vs->beachline[i]==arc)
        {
            vs->beachline.erase(i,1);
            break;
        }
    }
}

void CircleEvent::set_false_alarm()
{
    false_alarm = true;
}

Arc::Arc(Site *_site)
{
	site = _site;
	circ_event = NULL;
	left_edge = NULL;
	right_edge = NULL;
	leftarc = NULL;
	rightarc = NULL;
}

Arc::Arc()
{
    site=NULL;
    circ_event=NULL;
    left_edge=NULL;
    right_edge=NULL;
    leftarc = NULL;
	rightarc = NULL;
}
