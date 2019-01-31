/*
  ACP (Adaptive Controlled Precision) Library
  for robust computational geometry

  Copyright (c) 2013-07-15
  Victor Joseph Milenkovic
  University of Miami
  vjm@miami.edu
  Elisha Sacks
  Purdue University
  eps@cs.purdue.edu

  This file contains code described in

Robust Complete Path Planning in the Plane
Victor Milenkovic, Elisha Sacks, and Steven Trac
Proceedings of the Workshop on the Algorithmic Foundations of Robotics (WAFR)
pages 37-52, 2012

   This code is under development.
   It is free for use for academic and research purposes.
   Permission from the authors is required for commercial use.
*/

#ifndef POLY_H
#define POLY_H

#include "object.h"
#include "ppoly.h"

namespace acp {

typedef Object<Parameter> Knot;
typedef InputObject<Parameter> InputKnot;
//class Knot : public Object<Parameter> {
//};

//class InputKnot :public Knot {
// public:
//  InputKnot (double x) { set(Parameter(x)); }
//};
   
class SubKnot :public Knot {
  PTR<Knot> l, r;
  unsigned int i, n;
  
  Parameter calculate () { return Parameter((l->get()*(n - i) + r->get()*i)/n); }
 public:
  SubKnot (PTR<Knot> l, PTR<Knot> r, unsigned int i, unsigned int n)
    : l(l), r(r), i(i), n(n) {}
};

class DInt {
 public:
  PTR<Knot> l, u;
  unsigned int n;
   
  DInt (PTR<Knot> l, PTR<Knot> u, unsigned int n) : l(l), u(u), n(n) {}
};

class Root;

typedef vector<ObjPTR<Parameter>> Roots;

class PolySolver {
  Roots linearRoots (PTR<Knot> l, PTR<Knot> u);
  Roots quadraticRoots (PTR<Knot> l, PTR<Knot> u);
  Roots descartesRoots (PTR<Knot> l, PTR<Knot> u);
  bool descartes1 (vector<DInt> &st, const DInt &i, int v, bool lflag);
  bool descartes2 (vector<DInt> &st, const DInt &i, int v);
  int descartes2int (const DInt &i, int v);
  ObjPTR<PPoly<Parameter>> poly;
  PPoly<Parameter> get () { return poly->get(); }
  PPoly<Parameter> getApprox (double e) { return poly->getApprox(e); }
  Primitive1(Degree, Object<PPoly<Parameter>>*, poly);
 public:
  PolySolver (Object<PPoly<Parameter>> *poly) : poly(poly) {}
  int deg () { return Degree(poly); }
  Roots getRoots ();
  Roots getRoots (PTR<Knot> l, PTR<Knot> u);
  Roots getRoots (double l, double u) {
    // return getRoots(new InputKnot(l), new InputKnot(u));
    return getRoots(new InputKnot(Parameter::constant(l)),
                    new InputKnot(Parameter::constant(u)));
  }
};

/*
class InputPoly : public Poly {
 public:
  InputPoly (int d, double *ai) {
    PPoly<Parameter> p;
    for (int i = 0; i <= d; ++i)
      p.add(Parameter(ai[i]), i);
    set(p);
  }
};
*/
class Root : public Object<Parameter> {
  friend class Poly;
 protected:
  ObjPTR<PPoly<Parameter>> p;
 public:
  Root (Object<PPoly<Parameter>> *p) : p(p) {}
  PTR<Object<PPoly<Parameter>>> getPoly () { return p; }
};

class LinearRoot : public Root {
  Parameter calculate () {
    PPoly<Parameter> q = p->get();
    return - q[0]/q[1];
  }
 public:
  LinearRoot (Object<PPoly<Parameter>> *p) : Root(p) {}
};

class QuadraticRoot : public Root {
  bool flag;
  
  Parameter calculate () {
    PPoly<Parameter> q = p->get();
    const Parameter &a = q[2], &b = q[1], &c = q[0];
    Parameter d = (b*b - 4.0*a*c).sqrt();
    if (b.sign() == -1)
      return flag ? 2.0*c/(d-b) : 0.5*(d-b)/a;
    return flag ? -0.5*(d+b)/a : -2.0*c/(d+b);
  }
  
 public:
  QuadraticRoot (Object<PPoly<Parameter>> *p, bool flag) : Root(p), flag(flag) {}
};

class PolyRoot : public Root {
  PTR<Knot> l, u;
  Parameter calculate () {
    PPoly<Parameter> q = p->get();
    int slb = q.value(l->get()).sign(), sub = q.value(u->get()).sign();
    if (uninitialized())
      return q.shrinkBracket(Parameter(l->get().ubP(), u->get().lbP()), sub);
    return q.shrinkBracket(getCurrentValue(), sub);
  }
 public:
  PolyRoot (Object<PPoly<Parameter>> *p, PTR<Knot> l, PTR<Knot> u) : Root(p), l(l), u(u) {}
};   

class CauchyBound : public Knot {
  ObjPTR<PPoly<Parameter>> p;
  
  Parameter calculate () {
    PPoly<Parameter> q = p->get();
    q.checkDegree();
    Parameter b(1.0), qd = q.lc();
    for (int i = 0; i < q.deg(); ++i) {
      Parameter bi = (q.a[i]/qd).abs();
      if (b < bi)
	b = bi;
    }
    return 1.0 + b;
  }
  
 public:
  CauchyBound (Object<PPoly<Parameter>> *p) : p(p) {}
};

Primitive1(QuadraticRoots, Object<PPoly<Parameter>> *, p);

Primitive2(Order, Knot *, a, Knot *, b);

Primitive3(Descartes, Object<PPoly<Parameter>> *, p, Knot *, l, Knot *, u);

bool zerop (const Parameter &p);

class Root2;

typedef ObjPTR<Parameter> RootPTR;
/*
class RootPTR : public PTR<Root> {
public:
  RootPTR () {}
  RootPTR (Root *r) : PTR<Root>(r) {}
  operator ObjPTR<Parameter> () const;
};
*/
class Root2PTR : public PTR<Root2> {
public:
  Root2PTR () {}
  Root2PTR (Root2 *r) : PTR<Root2>(r) {}
  operator ObjPTR<PV2> () const;
};

class Poly2 : public RefCnt {
public:
  Poly2 () {}
  virtual PPoly2 getPoly() = 0;
};

class Root2 : public Object<PV2> {
  friend class Poly2;
  PTR<Poly2> f, g;
public:
  Root2(Poly2 *f, Poly2 *g, PV2 &r) : f(f), g(g) {
    set(r);
  }
  PTR<Poly2> getF () { return f; }
  PTR<Poly2> getG () { return g; }
  PV2 calculate ();
  static PV2 polish(PV2 p, PPoly2 f, PPoly2 g);
private:

  //order of roots in v: bot, right, top, left
  typedef struct {
    PV2 I;
    vector< vector<RootPTR> > v;
  } RootBoundary;

  static bool polishing;
  static std::vector<int> intersections(RootBoundary &b, PPoly2 f, PPoly2 g);
  static int parity(std::vector<int> alt);
  static void newton(RootBoundary &I, PPoly2 f, PPoly2 g);
  static void setRoots(RootBoundary &rb, PPoly2 f, PPoly2 g);
  static RootBoundary splitVert(RootBoundary &rb, Parameter c, PPoly2 f, PPoly2 g);
  static RootBoundary splitHoriz(RootBoundary &rb, Parameter c, PPoly2 f, PPoly2 g);
  static bool subdivide(RootBoundary &I, PPoly2 f, PPoly2 g);
  static bool solve(PV2, PV2, PV2, PV2 &);
};

inline Root2PTR::operator ObjPTR<PV2> () const {
  return ObjPTR<PV2>(operator Root2 *());
}

}
#endif
