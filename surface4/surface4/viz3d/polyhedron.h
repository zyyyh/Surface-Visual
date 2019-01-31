#ifndef POLY
#define POLY

#include <sys/time.h>
#include <set>
#include <map>
#include <iomanip>
#include <fstream>

#include "object.h"
#include "pv.h"
#include "rbtree.h"
#include "octree.h"

using namespace std;
using namespace acp;

double getTime ();

typedef vector<PV3> PV3s;

Parameter cross (const PV3 &a, const PV3 &b, const int coord);

typedef unsigned int ID;

typedef vector<ID> IVector;

typedef vector<IVector> IVectors;

typedef set<ID> IDSet;

IDSet intersection (const IDSet &a, const IDSet &b);

class PlaneData {
 public:
  PlaneData () {}
  PlaneData (const PV3 &n, const Parameter &k) : n(n), k(k) {}
  int size () const { return 4; }
  Parameter & operator[](int i) { return i < 3 ? n[i] : k; }
  PV3 n;
  Parameter k;
};

class Plane : public Object<PlaneData> {
  friend class Point;
  friend class PPoint;
  friend class Face;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  static ID planeid;
  ID id;
 public:
  PV3 getN () { return get().n; }
  Parameter getK () { return get().k; }
  ID getid () const { return id; }
};

typedef vector<Plane *> Planes;

Primitive2(SamePlane, Plane *, p, Plane *, q);

Primitive1(ProjectionCoordinate, Plane *, p);

class Point : public Object<PV3> {
  friend class Plane;
  friend class TrianglePlane;
  friend class EPPoint;
  friend class EdgeProjectionPoint;
  friend class MidPoint;
  friend class CentroidPoint;
  friend class Vertex;
  friend class Edge;
  friend class Face;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  IDSet ps;
 public:
  PV3 getP () { return get(); }
  IDSet getps () const { return ps; }
  void getBBox (double *bbox);
  bool identical (Point *p);
  int order (Point *p);
  bool identicalI (Point *p);
  bool onLine (Point *a, Point *b);
  int side (Plane *a);
};

typedef vector<PTR<Point> > Points;

Primitive3(TripleProduct, Point *, a, Point *, b, Point *, c);

Primitive3(OnLine, Point *, a, Point *, t, Point *, h);

Primitive3(Order, Point *, a, Point *, t, Point *, h);

bool onEdge (Point *a, Point *t, Point *h, bool strict);

Primitive2(Side, Plane *, p, Point *, a);

Primitive3(PointOrder, Point *, a, Point *, b, Point *, r);

Primitive2(PointOrderR, Point *, a, Point *, b);

class LeftTurn : public Primitive {
  Point *a, *b, *c;
  int pc;
  int sign () {
    return cross(c->getP() - b->getP(), a->getP() - b->getP(), pc).sign();
  }
 public:
  LeftTurn (Point *a, Point *b, Point *c, int pc) 
    : a(a), b(b), c(c), pc(pc) {}
};

Primitive4(Orientation, Point *, a, Point *, b, Point *, c, Point *, d);

Primitive4(CloserPair, Point *, a, Point *, b, Point *, c, Point *, d);

Primitive2(PlaneRayAlignment, Plane *, p, Point *, u);

void copyBBox (const double *bbf, double *bbt);

void mergeBBox (const double *bbf, double *bbt);

bool bboxOverlap (const double *a, const double *b, double s = 0.0);

bool bboxOverlap (Point *a, const double *bbox);

class TrianglePlane : public Plane {
  friend class EPPoint;
  friend class RayPlanePoint;
 protected:
  PTR<Point> a, b, c;
  PlaneData calculate () {
    PV3 ap = a->getP(), bp = b->getP(), cp = c->getP(),
      n = (cp - bp).cross(ap - bp);
    Parameter k = - n.dot(ap);
    return PlaneData(n, k);
  }  
 public:
  TrianglePlane () : a(0), b(0), c(0) {}
  TrianglePlane (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {
    IDSet ps = intersection(a->ps, intersection(b->ps, c->ps));
    if (ps.empty()) {
      id = planeid++;
      a->ps.insert(id); b->ps.insert(id); c->ps.insert(id);
    }
    else
      id = *ps.begin();
  }
  
  Point * getA () { return a; }
  Point * getB () { return b; }
  Point * getC () { return c; }
};

class InputPoint : public Point {
 public:
  InputPoint (double x, double y, double z, bool perturb = true) {
    set(perturb ? PV3::input(x, y, z) : PV3::constant(x, y, z));
  }
};

extern InputPoint Rdir;

class NegPoint : public Point {
  friend class Point;
  PTR<Point> a;
  PV3 calculate () { return - a->getP(); }
 public:
  NegPoint (Point *a) : a(a) {}
};

class SumPoint : public Point {
  friend class Point;
  PTR<Point> a, b;
  PV3 calculate () { return a->getP() + b->getP(); }
 public:
  SumPoint (Point *a, Point *b) : a(a), b(b) {}
};

class DiffPoint : public Point {
  friend class Point;
  PTR<Point> a, b;
  PV3 calculate () { return a->getP() - b->getP(); }
 public:
  DiffPoint (Point *a, Point *b) : a(a), b(b) {}
};

class MidPoint : public Point {
  friend class Point;
  PTR<Point> a, b;
  PV3 calculate () { return (a->getP() + b->getP())/2; }
 public:
  MidPoint (Point *a, Point *b) : a(a), b(b) {
    ps = intersection(a->ps, b->ps);
  }
};

class CentroidPoint : public Point {
  friend class Point;
  PTR<Point> a, b, c;
  PV3 calculate () { return (a->getP() + b->getP() + c->getP())/3; }
 public:
  CentroidPoint (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {
    ps = intersection(a->ps, b->ps);
    ps = intersection(ps, c->ps);
  }
};

class ScalePoint : public Point {
  PTR<Point> p;
  double unit;
  PV3 calculate () { return unit * p->getP(); }
 public:
  ScalePoint(PTR<Point> p, double unit) : p(p), unit(unit) {}
};

class EPPoint : public Point {
  PTR<Point> t, h, a, b, c;
  PV3 calculate () {
    PV3 ap = a->getP(), n = (b->getP() - ap).cross(c->getP() - ap),
      tp = t->getP(), u = h->getP() - tp;
    Parameter k = n.dot(ap - tp)/n.dot(u);
    return tp + k*u;
  }
 public:
  EPPoint (Point *t, Point *h, TrianglePlane *q)
    : t(t), h(h), a(q->a), b(q->b), c(q->c) {
    IDSet psth = intersection(t->ps, h->ps);
    for (IDSet::iterator i = psth.begin(); i != psth.end(); ++i)
      ps.insert(*i);
    ps.insert(q->id);
  }
};

class RayPlanePoint : public Point {
  PTR<Point> t, r;
  TrianglePlane p;
  PV3 calculate () {
    PV3 a = t->getP(), u = r->getP(), n = p.getN();
    Parameter k = -(n.dot(a) + p.getK())/n.dot(u);
    return a + k*u;
  }
 public:
  RayPlanePoint (Point *t, Point *r, TrianglePlane *q)
    : t(t), r(r), p(q->a, q->b, q->c) {}
};

class RayZPlanePoint : public Point {
  PTR<Point> t, r;
  double z;
  PV3 calculate () {
    PV3 a = t->getP(), u = r->getP();
    Parameter k = (z - a.z)/u.z;
    return a + k*u;
  }
 public:
  RayZPlanePoint (Point *t, Point *r, double z)
    : t(t), r(r), z(z) {}
};

class Edge;

class HEdge;

class Face;

class Vertex {
  friend class Edge;
  friend class HEdge;
  friend class Face;
  friend class Shell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  PTR<Point> p;
  vector<Edge *> edges;
  double rint[2], bbox[6];
  RBTree<Vertex*>::Node *node;
 public:
  Vertex (Point *p, bool perturbed);
  Point * getP () { return p; }
  int EdgesN () const { return edges.size(); }
  Edge * getEdge (int i) const { return edges[i]; }
  double * getBBox () { return bbox; }
  void outgoingHEdges (vector<HEdge *> &ed) const;
  vector<Face *> incidentFaces () const;
  HEdge * connected (Vertex *a) const;
  int order (Vertex *v) const;
};

typedef vector<Vertex *> Vertices;

typedef vector<const Vertices *> VVertices;

typedef set<Vertex *> VertexSet;

class Triangle {
 public:
  Vertex *a, *b, *c;
  Triangle (Vertex *a, Vertex *b, Vertex *c) : a(a), b(b), c(c) {}
};

typedef vector<Triangle> Triangles;

class Vertex;

class HEdge {
  friend class Edge;
  friend class Face;
  friend class HFace;
  friend class Shell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  Edge *e;
  bool forward, flag, iflag;
  Face *f;
  HEdge *next;
 public:
  HEdge () : e(0), forward(false), flag(false), iflag(false), next(0), f(0) {}
  HEdge (Edge *e, bool forward) 
    : e(e), forward(forward), flag(false), iflag(true), next(0), f(0) {}
  Vertex * tail () const;
  Vertex * head () const;
  Edge * getE () const { return e; }
  bool getForward () const { return forward; }
  bool getFlag () const { return flag; }
  void setFlag (bool f) { flag = f; }
  Face * getF () const { return f; }
  HEdge * getNext () const { return next; }
  void setNext (HEdge *h) { next = h; }
  PV3 getU ();
  PV3 getN ();
  HEdge * cw () const;
  HEdge * ccw () const;
  void loop (Vertices &ve);
  HEdge * findLoop (Vertices &ve);
  void loop (vector<HEdge *> &ed);
};

typedef vector<HEdge *> HEdges;

typedef set<HEdge *> HEdgeSet;

class Edge {
  friend class EEPoint;
  friend class Vertex;
  friend class HEdge;
  friend class Face;
  friend class Shell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  Vertex *t, *h;
  double bbox[6];
  HEdge hei[2];		 
  vector<HEdge *> hedges;
  Vertices vertices;
  vector <Edge *> dedges;
 public:
  Edge (Vertex *t, Vertex *h);
  ~Edge ();
  void setBBox ();
  HEdge * addHEdge (bool forward);
  void removeHEdge (HEdge *e);
  Vertex * getT () const { return t; }
  Vertex * getH () const { return h; }
  PV3 getU () { return h->p->getP() - t->p->getP(); }
  int HEdgesN () const { return hedges.size(); }
  HEdge * getHEdge (int i) const { return i < hedges.size() ? hedges[i] : 0; }
  double * getBBox () { return bbox; }
  IDSet ps () const { return intersection(t->p->ps, h->p->ps); }
  VertexSet vertexSet () const { 
    return VertexSet(vertices.begin(), vertices.end()); 
  }
  void edgeVertices (Vertices &ve) const;
  Face * otherFace (Face *f) const;
  void sortVertices ();
  void sortHEdges ();
};

typedef vector<Edge *> Edges;

typedef set<Edge *> EdgeSet;

Primitive3(EdgeVertexOrderP, Edge *, e, Vertex *, v, Vertex *, w);

Primitive3(EdgeOrderP, Edge *, e, HEdge *, f, HEdge *, g);

class EdgeOrder {
 public:
  EdgeOrder (Edge *e) : e(e) {}
  bool operator() (HEdge *f, HEdge *g) const {
    return f != g && EdgeOrderP(e, f, g) == 1;
  }

  Edge *e;
};

class EdgeVertexOrder {
 public:
  EdgeVertexOrder (Edge *e) : e(e) {}
  bool operator() (Vertex *v, Vertex *w) const {
    return v != w && EdgeVertexOrderP(e, v, w) == 1;
  }

  Edge *e;
};

class EEPoint : public Point {
  PTR<Point> et, eh, ft, fh;
  int pc;
  PV3 calculate () {
    PV3 a = et->getP(), u = eh->getP() - a, b = ft->getP(), v = fh->getP() - b;
    Parameter k = cross(b - a, v, pc)/cross(u, v, pc);
    return a + k*u;
  }
 public:
  EEPoint (Edge *e, Edge *f, int pc) : et(e->t->getP()), eh(e->h->getP()), 
    ft(f->t->getP()), fh(f->h->getP()), pc(pc) {
    IDSet pe = e->ps(), pf = f->ps(), pef;
    set_union(pe.begin(), pe.end(), pf.begin(), pf.end(),
	      inserter(pef, pef.begin()));
    for (IDSet::iterator i = pef.begin(); i != pef.end(); ++i)
      ps.insert(*i);
  }

  EEPoint (Point *et, Point *eh, Point *ft, Point *fh) : et(et), eh(eh), ft(ft), fh(fh) {
    TrianglePlane p(et, eh, ft);
    pc = ProjectionCoordinate(&p);
  }
};

class Shell;

class Polyhedron;

typedef pair<Point *, Vertex*> PVPair;

typedef map<Point *, Vertex *> PVMap;

class HFace {
  friend class Face;
  friend class Shell;
  friend class Cell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  Face *f;
  Shell *s;
 public:
  HFace () : f(0), s(0) {}
  Face * getF () const { return f; }
  Shell * getS () const { return s; }
  bool pos () const;
  HFace * twin () const;
  PV3 getN () const;
  void neighbors (vector<HFace *> &hfaces) const;
  HFace * neighbor (HEdge *e) const;
};

typedef vector<HFace *> HFaces;

class HFaceNormal : public Point {
  HFace *f;
  PV3 calculate () { return f->getN(); }
 public:
  HFaceNormal (HFace *f) : f(f) {}
};

void triangulate (const VVertices &reg, int coord, Triangles &tr);

class Face {
  friend class HEdge;
  friend class HFace;
  friend class Shell;
  friend class Cell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  TrianglePlane p;
  int pc;
  HEdges boundary;
  HFace hfaces[2];
  double bbox[6];
  Edges edges;
 public:
  Face (const Vertices &b, Vertex *u, Vertex *v, Vertex *w, int pc, bool flag);
  void update ();
  void addBoundary (HEdge *e);
  TrianglePlane * getP () { return &p; }
  HEdge * getBoundary (int i) const { return boundary[i]; }
  const HEdges &getBoundary () { return boundary; } 
  double * getBBox () { return bbox; }
  HFace * getHFace (int i) { return hfaces + i; }
  const Edges & getEdges () const { return edges; }
  int getPC ();
  bool boundaryVertex (Point *a) const;
  bool boundaryVertex (Vertex *v) const;
  void boundaryVertices (Vertices &ve) const;
  bool boundaryEdge (Edge *e) const;
  void boundaryHEdges (HEdges &ed) const;
  void sharedBoundaryVertices (Face *f, Vertices &vfg) const;
  void sharedBoundaryVertices1 (const Face *f, Vertices &vfg) const;
  void containedBoundaryVertices (Face *f, Vertices &vfg);
  void containedBoundaryVertices1 (const Face *f, Vertices &vfg);
  bool sharedBoundaryEdge (Face *f) const;
  void edgeVertices (VertexSet &vs) const;
  bool coplanar (Face *f);
  PTR<Point> rayIntersection (Point *a, Point *r);
  bool contains (Point *a, bool strict);
  bool triangle () const;
  bool containsConvex (Point *a, bool strict);
  bool boundaryContains (Point *a, int i);
  bool intersects (Edge *e, bool strict);
  bool intersectsEE (Edge *e, Edge *f);
  PTR<Point> centroid () const;
  void triangulate (Triangles &tr);
};

typedef vector<Face *> Faces;

typedef set<Face *> FaceSet;

Face * newFace (const Vertices &ve, int pc, bool flag);

void extremalR (const Vertices &ve, Vertex *&u, Vertex *&v, Vertex *&w);

class RayEdgeIntersection : public Primitive {
  Point *a, *r, *t, *h;
  int c;
  int sign ();
 public:
  RayEdgeIntersection (Point *a, Point *r, Point *t, Point *h, int c)
    : a(a), r(r), t(t), h(h), c(c) {}
};

void deleteRegion (const VVertices &reg);

typedef Octree<Face *> FOctree;

class Cell;

class Shell {
  friend class Cell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  HFaces hfaces;
  double bbox[6];
  FOctree *octreef;
  Cell *c;
  void setBBox ();
  void setOctree ();
 public:
  Shell () : octreef(0) {};
  ~Shell ();
  Shell (const HFaces &hf);
  HFaces & getHFaces () { return hfaces; }
  double * getBBox () { return bbox; }
  Cell * getC() const { return c; }
  bool outer () const;
  Vertex * vmax (Point *r) const;
  Vertex * vmax (Face *f, Point *r) const;
  bool contains (Shell *s) const;
  int contains (Point *a) const;
  void rayBBox (Point *a, Point *r, double *rb) const;
  int euler () const;
};

typedef vector<Shell *> Shells;

void deleteShells (const Shells &sh);

Primitive3(SlopeOrder, Edge *, e, Edge *, f, Point *, x);

Primitive2(Convex, HEdge *, e, HFace *, f);

class Cell {
  friend class Polyhedron;
  friend class Convolution;
 protected:
  Shell *outer;
  Shells inner;
  int wn;
 public:
  Cell (Shell *outer) : outer(outer), wn(0) { if (outer) outer->c = this; }
  ~Cell () { delete outer; deleteShells(inner); }
  int nShells () const { return inner.size() + (outer ? 1 : 0); }
  
  Shell * getShell (int i) const {
    if (outer)
      return i == 0 ? outer : inner[i-1];
    return inner[i];
  }
    
  double * getBBox () { return outer ? outer->getBBox() : 0; }
  void addInner (Shell *s) { inner.push_back(s); s->c = this; }
  bool contains (Point *p) const;
  PTR<Point> interiorPoint () const;
  int getWN () const { return wn; }
};

typedef vector<Cell *> Cells;

typedef set<Cell *> CellSet;

class Attribute {
 public:
  Attribute (ID id, string name, string value) 
    : id(id), name(name), value(value) {}

  ID id;
  string name, value;
};

typedef vector<Attribute> Attributes;

typedef pair<Edge *, Face *> EFPair;

typedef pair<EFPair, Vertex *> EFVPair;

typedef map<EFPair, Vertex *> EFVMap;

typedef pair<Face *, Face *> FFPair;

FFPair ffpair (Face *f, Face *g);

typedef map<Face *, Face *> FFMap;

typedef pair<Vertex *, Vertex *> VVPair;

typedef set<VVPair> VVPairSet;

class VertexHHEdgeOrder : public Primitive {
  Vertex *v;
  HEdge *e, *f;
  int c;
  int sign ();
 public:
  VertexHHEdgeOrder (Vertex *v, HEdge *e, HEdge *f, int c) 
    : v(v), e(e), f(f), c(c) {}
};

class HHEdge {
 public:
  HEdge *e;
  bool f;
  HHEdge (HEdge *e, bool f) : e(e), f(f) {}
  Vertex * tail () const { return f ? e->tail() : e->head(); }
  Vertex * head () const { return f ? e->head() : e->tail(); }
};

typedef vector<HHEdge> HHEdges;

class HHEdgeOrder {
 public:
  HHEdgeOrder (int c) : c(c) {}
  bool operator() (const HHEdge &e, const HHEdge &f) const {
    Vertex *u = e.tail(), *v = f.tail();
    if (u != v)
      return u < v;
    if (e.head() == f.head())
      return f.f < e.f;
    return VertexHHEdgeOrder(u, e.e, f.e, c) == 1;
  }

  int c;
};

class HEdgeOrder {
 public:
  bool operator() (HEdge *e, HEdge *f) const {
    if (e->getE() == f->getE())
      return false;
    Point *et = e->tail()->getP(), *eh = e->head()->getP(),
      *ft = f->tail()->getP(), *fh = f->head()->getP(),
      *pe = PointOrderR(et, eh) == 1 ? et : eh,
      *pf = PointOrderR(ft, fh) == 1 ? ft : fh;
    return pe != pf && PointOrderR(pe, pf) == 1;
  }
};

Primitive4(FFOrderP, Face *, f, Face *, g, Vertex *, v, Vertex *, w);

class FFOrder {
  Face *f, *g;
 public:
  FFOrder (Face *f, Face *g) : f(f), g(g) {}
  bool operator() (Vertex *v, Vertex *w) const {
    return v != w && FFOrderP(f, g, v, w) == 1;
  }
};

class FFF {
 public:
  FFF (Face *f1, Face *f2, Face *f3) {
    Face *fa[] = {f1, f2, f3};
    sort(fa, fa + 3);
    a = fa[0]; b = fa[1]; c = fa[2];
  }

  bool operator< (const FFF &x) const {
    return a < x.a || (a == x.a && b < x.b) || (a == x.a && b == x.b && c < x.c);
  }

  Face *a, *b, *c;
};

typedef pair<FFF, Vertex *> FFFVPair;

typedef map<FFF, Vertex *> FFFVMap;

typedef pair<Edge *, Edge *> EEPair;

typedef pair<EEPair, Vertex *> EEVPair;

typedef map<EEPair, Vertex *> EEVMap;

typedef pair<Point *, Point *> PPPair;

typedef pair<PPPair, EdgeSet> PPEPair;

typedef map<PPPair, EdgeSet> PPEMap;

enum SetOp { Union, Intersection, Complement };

class FaceDesc {
  VertexSet vs;
 public:
  FaceDesc () {}
  FaceDesc (Vertex **ve, int n) { vs.insert(ve, ve + n); }
  FaceDesc (const Vertices &ve) { vs.insert(ve.begin(), ve.end()); }

  bool operator< (const FaceDesc &f) const {
    if (vs.size() < f.vs.size())
      return true;
    if (f.vs.size() < vs.size())
      return false;
    for (VertexSet::const_iterator i = vs.begin(), j = f.vs.begin(); 
	 i != vs.end(); ++i, ++j)
      if (*i < *j)
	return true;
      else if (*j < *i)
	return false;
    return false;
  }
};

typedef set<FaceDesc> FaceDescSet;

class Polyhedron {
 public:
  bool perturbed, oneWayInt;
  RBTree<Vertex *> vtree;
  Vertices vertices;
  Edges edges;
  Faces faces;
  Cells cells;
  Attributes attributes;
  double bbox[6];
  EEVMap eevmap;
  EdgeSet iedges;
  PPEMap ppemap;

  Polyhedron (bool perturbed = true) : perturbed(perturbed), oneWayInt(false) {} 
  ~Polyhedron ();
  Vertex * getVertex (Point *p);
  Vertex * getVertex (double x, double y, double z) {
    Point *p = new InputPoint(x, y, z, perturbed);
    return getVertex(p);
  }
  Vertex * getVertex (Vertex *v, PVMap &pvmap);
  HEdge * addHEdge (Vertex *a, Vertex *b);
  Edge * getEdge (Vertex *a, Vertex *b);
  HEdge * getHEdge (Vertex *a, Vertex *b);
  HEdge * findHEdge (Vertex *a, Vertex *b);
  void addVertex (Edge *e, Vertex *v);
  Face * addFace (const VVertices &reg, int pc = 0, bool flag = false);
  HEdge * addLoop (const Vertices &ve);
  HEdge * addLoop (Vertex **v, int n);
  Face * addFace (const Vertices &ve, int pc = 0, bool flag = false);
  Face * addFace (HEdge *e, int pc = 0, bool flag = false);
  Face * addTriangle (Vertex *a, Vertex *b, Vertex *c, int pc = 0, bool flag = false);
  Face * addRectangle (Vertex *a, Vertex *b, Vertex *c, Vertex *d);
  void formCells ();
  void formCellsAux (const Shells &shells);
  void formShells (Shells &shells);
  Shell * formShell (HFace *f) const;
  Cell * enclosingCell (Shell *s, Octree<Cell *> *octreec) const;
  Polyhedron * triangulate () const;
  Face * getTriangle (Vertex *ta, Vertex *tb, Vertex *tc, PVMap &pvmap, int pc = 0);
  Face * getTriangle (Face *f, PVMap &pvmap);
  Polyhedron * subdivide ();
  void intersectFF ();
  void intersectFF (Face *f, Face *g, EFVMap &efvmap);
  void intersectFFP (Face *f, Face *g, EFVMap &efvmap);
  void intersectFEP (Face *f, Edge *e, EFVMap &efvmap);
  Vertex * intersectEE (Edge *e, Edge *f, int pc);
  bool intersectEV (Edge *e, Vertex *v);
  void intersectEEL (Edge *e, Edge *f);
  bool intersectPE (Face *f, Face *g, EFVMap &efvmap, Edges &iedges, Faces &ifaces);
  void intersectFV (Face *f, Vertex *v);
  void intersectEF (Edge *e, Face *f, EFVMap &efvmap, Vertices &vfg);
  void formFF (Face *f, Face *g, Vertices &vfg);
  Polyhedron * subdivideAux ();
  void intersectFFF (Face *f, FFFVMap &fffvmap);
  void subfaces (Face *f, Polyhedron *a, PVMap &pvmap);
  void subedges (Face *f, Polyhedron *a, PVMap &pvmap, HEdges &he);
  void subedges (HEdge *e, VVPairSet &vvps);
  void subedge (Edge *e, Vertex *t, Vertex *h, VVPairSet &vvps);
  void setNext (const HEdges &he, int c) const;
  Face * enclosingFace (const Faces &fa, Point *p) const;
  void removeLoops (const HEdges &ed);
  Polyhedron * copy () const;
  Polyhedron * scale (double unit) const;
  Polyhedron * negative () const;
  Polyhedron * translate (Point *t) const;
  Polyhedron * negativeTranslate (Point *t) const;
  bool intersects (Polyhedron *a, bool strict) const;
  bool contains (Point *p) const;
  int containingCell (Point *p) const;
  bool intersectsEdges (const Polyhedron *a, bool strict) const;
  Polyhedron * boolean (Polyhedron *a, SetOp op);
  Polyhedron * cellPolyhedron (int i) const;
  void addHFaces (const HFaces &hf, PVMap &pvmap);
  void replaceVertex (Face *f, Vertex *v, Vertex *w);
  void removeLoop (HEdge *e);
  void removeHEdge (HEdge *e);
  void moveVertex (Vertex *v, Point *p);
  void removeNullFaces ();
  void updateCells ();
  Cell * updateCell (Cell *c) const;
  Shell * updateShell (Shell *s) const;
  FOctree * faceOctree (double s = 0.0) const;
  Octree<Cell *> * cellOctree () const;
  Polyhedron * round (double d);
  Polyhedron * selfUnion ();
  void computeWindingNumbers ();
  void updateWN (Cell *c, HFaces &st) const;
  Polyhedron * encaseNonManifold (double d);
  bool manifold (Vertex *v) const;
  void addTetrahedron (Point *o, double r);
  void describe () const;
};

typedef vector<Polyhedron *> Polyhedrons;

Plane * facePlane (const VVertices &reg, Polyhedron *a, bool check);

bool outerLoop (const Vertices &ve, int pc);

bool inSet (bool ina, bool inb, SetOp op);

Polyhedron * overlay (Polyhedron **poly, int n);

Polyhedron * multiUnion (Polyhedron **poly, int n);

Polyhedron * box (double *b, bool perturb = true);

Polyhedron * sphere (double ox, double oy, double oz, double r, double err);

Polyhedron * sphere (double err);

Polyhedron * octohedron ();

double sphereError (Polyhedron *a);

Polyhedron * sphereRefine (Polyhedron *a);

Polyhedron * lbox (double x1, double x2, double y1, double y2, double z1, double z2);

Polyhedron * room (double x1, double x2, double x3, double y1, double y2,
		   double y3, double y4, double y5, double z1, double z2,
		   bool perturb = true);

#endif
