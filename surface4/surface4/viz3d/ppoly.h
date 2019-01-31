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

#ifndef PPOLY
#define PPOLY

#include "pv.h"
#include <vector>
#include <unordered_map>
using namespace std;
using namespace acp;

namespace acp {

template<class T>
class PPoly {
 public:
  PPoly () {}
  PPoly (double c) { add(T(c), 0); }
  int size () const { return a.size(); }
  T &operator[] (int i) { return a[i]; }
  T lc () const { return *a.rbegin(); }
  int deg () const { return a.size() - 1; }
  bool zero () const { return a.empty(); }

  void add (const T &c, int e) {
    int d = deg();
    if (e <= d)
      a[e] = a[e] + c;
    else {
      while (d + 1 < e) {
	a.push_back(T(0.0));
	++d;
      }
      a.push_back(c);
    }
  }

  PPoly operator+ (const PPoly &p) const {
    int d = deg(), e = p.deg();
    PPoly q;
    if (d <= e) {
      for (int i = 0; i <= d; ++i)
	q.add(a[i] + p.a[i], i);
      for (int i = d + 1; i <= e; ++i)
	q.add(p.a[i], i);
    }
    else {
      for (int i = 0; i <= e; ++i)
	q.add(a[i] + p.a[i], i);
      for (int i = e + 1; i <= d; ++i)
	q.add(a[i], i);
    }
    return q;
  }

  PPoly operator- (const PPoly &p) const { return (*this) + (- p); }

  PPoly operator- () const {
    PPoly q;
    for (int i = 0; i <= deg(); ++i)
      q.add(- a[i], i);
    return q;
  }

  PPoly operator* (const PPoly &p) const {
    int d = deg(), e = p.deg();
    PPoly q;
    for (int i = 0; i <= d; ++i)
      for (int j = 0; j <= e; ++j)
	q.add(a[i]*p.a[j], i + j);
    return q;
  }
  
  PPoly rem (const PPoly &p, PPoly &q) const {
    T lp = p.lc();
    int dp = p.deg();
    PPoly r(*this);
    while (dp <= r.deg()) {
      T k = r.lc()/lp;
      int dk = r.deg() - dp;
      q.add(k, dk);
      r.a.pop_back();
      for (int i = 0; i < dp; ++i)
	r.add(- k*p.a[i], i + dk);
      r.checkDegree();
    }
    return r;
  }

  void checkDegree () {
    while (!zero() && zerop(lc()))
      a.pop_back();
  }

  PPoly gcd (const PPoly &p, PPoly &s, PPoly &t) const
  {
    s = PPoly(1);
    t = PPoly(0);
    PPoly r(*this), nr(p), ns, nt(1);
    while (!nr.zero()) {
      PPoly q, nnr = r.rem(nr, q), nns = s - q*ns, nnt = t - q*nt;
      r = nr;
      nr = nnr;
      s = ns;
      ns = nns;
      t = nt;
      nt = nnt;
    }
    return r;
  }

  T value (const T &x) const
  {
    int d = deg();
    T y = a[d];
    for (int i = d - 1; i >= 0; --i)
      y = x*y + a[i];
    return y;
  }

  T der (const T &x) const
  {
    int d = deg();
    T y = d*a[d];
    for (int i = d - 1; i > 0; --i)
      y = x*y + i*a[i];
    return y;
  }

  PPoly der () const
  {
    PPoly g;
    for (int i = 1; i <= deg(); ++i)
      g.a.push_back(i*a[i]);
    return g;
  }
  
  PPoly neg () const
  {
    PPoly p;
    for (unsigned int i = 0u; i <= deg(); ++i)
      p.a.push_back(i%2 == 0 ? a[i] : - a[i]);
    return p;
  }

  PPoly dual () const
  {
    PPoly p;
    int d = deg();
    for (unsigned int i = 0u; i <= d; ++i)
      p.a.push_back(a[d-i]);
    return p;
  }

  T shrinkBracket (T x, int sub) const {
    bool bflag = false;
    while (true) {
      bool flag = true;
      T xm = x.midP(), fm = value(xm), df = der(x);
      if (df.sign(false) == 0)
	flag = false;
      else {
	T nx = (xm - fm/df).intersect(x);
	if (nx.subset(x))
	  x = nx;
	else
	  flag = false;
      }
      if (!flag) {
	int sm = fm.sign(false);
	if (sm == sub) {
	  T nx = x.interval(xm);
	  if (!nx.subset(x))
	    break;
	  x = nx;
	}
	else if (sm == - sub) {
	  T nx = xm.interval(x);
	  if (!nx.subset(x))
	    break;
	  x = nx;
	}
	else if (!bflag) {
	  bflag = true;
	  x = shrinkBracketB(x, sub);
	}
	else
	  break;
      }
    }
    return x;
  }

  T shrinkBracketB (const T &x, int sub) const {
    T xlb = x.lbP(), xm = x.midP(), xub = x.ubP();
    while ((xlb - xm).sign(false) < 0) {
      T nx = xlb.interval(xm).midP();
      if ((nx - xlb).sign(false) == 0 || (nx - xm).sign(false) == 0)
	break;
      int sm = value(nx).sign(false);
      if (sm == 0)
	xm = nx;
      else if (sm == sub)
	xub = nx;
      else
	xlb = nx;
    }
    xm = x.midP();
    while ((xm - xub).sign(false) < 0) {
      T nx = xm.interval(xub).midP();
      if ((nx - xm).sign(false) == 0 || (nx - xub).sign(false) == 0)
	break;
      int sm = value(nx).sign(false);
      if (sm == 0)
	xm = nx;
      else if (sm == -sub)
	xlb = nx;
      else
	xub = nx;
    }
    T nx = xlb.interval(xub);
    return nx.subset(x) ? nx : x;
  }

  PPoly moebius (const T &l, const T &u) const
  {
    return shift(l).dual().shift(1.0/(u - l));
  }

  PPoly shift (const T &s) const
  {
    PPoly p = *this;
    if (s.sign() == 0)
      return p;
    int d = deg();
    for (int i = 0; i < d; ++i)
      for (int j = d - 1; j >= i; --j)
	p.a[j] = p.a[j] + s * p.a[j+1];
    return p;
  }

  PPoly shiftX (const T &s) const
  {
    PPoly p = *this;
    if (s.sign() == 0)
      return p;
    p.scale(s);
    p.shift1();
    p.scale(1.0/s);
    return p;
  }

  void scale (const T &s)
  {
    T ss = s;
    for (int i = 1; i <= deg(); ++i) {
      a[i] = ss*a[i];
      ss = ss*s;
    }
  }

  void shift1 ()
  {
    int d = deg();
    for (int i = 0; i < d; ++i)
      for (int j = d - 1; j >= i; --j)
	a[j] = a[j] + a[j+1];
  }

  void print () const {
    int d = deg();
    cerr << "( ";
    for (int i = 0; i <= d; ++i)
      cerr << a[i].mid() << " ";
    cerr << ")" << endl;
  }

  vector<T> a;
};

template<class T>
bool zerop (const PPoly<T> &p) {
  return p.zero() || (p.deg() == 0 && zerop(p.lc()));
 }

typedef std::pair<Parameter, Parameter> Int;
typedef std::vector<Int> Ints;
class Mobius;

class PPoly1 : public PPoly<Parameter> {
public:

  static bool print_all;
  int d;

  PPoly1 () : d(-1) {}
  PPoly1 (int d, const Parameter *a = 0) : d(d) {
    if (a == 0)
      add(Parameter::constant(0), d);
    else
      for (int i = 0; i <= d; i++)
        add(a[i], i);
  }
  PPoly1 (const Parameter &c, int d) : d(d) { add(c, d); }
  PPoly1 (const PPoly<Parameter> &ppoly) : PPoly<Parameter>(ppoly), d(ppoly.deg()) {}

  bool increased () const { return a[0].increased(); }
  bool decreased () const { return a[0].decreased(); }
  Parameter value (const Parameter &x) const;
  Parameter der (const Parameter &x) const;
  PPoly1 der () const;

  static PPoly1 XminusR (const Parameter &r) {
    Parameter a[2] = { -r, Parameter::constant(1) };
    return PPoly1(1, a);
  }

  std::vector<Parameter> roots (const Parameter &lb, const Parameter &ub) const
    { return rootsR(lb, ub); }
  
  Parameter polish (const Parameter &x) const;
  Parameter shrinkBracket (Parameter x, int sub) const;
  Parameter shrinkBracketB (const Parameter &x, int sub) const;
  std::vector<Parameter> rootsR (const Parameter &lb, const Parameter &ub) const;
  std::vector<Parameter> roots1 (const Parameter &lb, const Parameter &ub) const;
  std::vector<Parameter> roots2 (const Parameter &lb, const Parameter &ub) const;

  PPoly1 operator* (const PPoly1 &g) const;
  PPoly1 operator* (const Parameter &p) const;
  PPoly1 operator* (const double p) const;
  PPoly1 plusCtimes (double c, const PPoly1 &b) const;
  PPoly1 plusCtimes (const Parameter &c, const PPoly1 &b) const;
  PPoly1 operator+ (const PPoly1 &b) const { return plusCtimes(1.0, b); }
  PPoly1 operator- (const PPoly1 &b) const { return plusCtimes(-1.0, b); }

  PPoly1 karatsuba (const PPoly1 &g) const;

  PPoly1 normalized() const;

  PPoly1 quotient (const PPoly1 &b, PPoly1 &remainder);
  PPoly1 gcd (const PPoly1 &b);
  int degree () { return deg(); }

  void print();

  int size () { return a.size(); }
  Parameter &operator[] (int i) { return a[i]; }

  PPoly1 monic () const;
  PPoly1 normal () const;
  PPoly1 shift (double s) const;
  PPoly1 shift (const Parameter &s) const;
  PPoly1 neg () const;
  PPoly1 dual () const;
  PPoly1 removeZeroRoots () const;
  unsigned int descartes () const;
  Parameter ubk () const;

  std::vector<Parameter> rootsD (const Parameter &lb, const Parameter &ub) const;
  Ints rootInts (const Parameter &lb, const Parameter &ub) const;
  void vas (PPoly1 p, Mobius m, Ints &ints) const;
  Int vasInt (const PPoly1 &p, const Mobius &m) const;
};

class PPoly2 {
 public:
  PPoly2 () : nt(0), a(0), m(0), degx(-2), degy(-2) {}
  PPoly2 (int nt) : nt(nt), a(new Parameter [nt]), m(new int [2*nt]), degx(-2), degy(-2) {}
  PPoly2 (int nt, Parameter *a, int *m);
  PPoly2 (const PPoly2 &p);
  void setDegree ();
  ~PPoly2 () {delete [] a; delete [] m;}
  PPoly2 & operator= (const PPoly2 &p);
  int degree () const { return d; }
  bool increased () const { return a[0].increased(); }
  bool decreased () const { return a[0].decreased(); }
  Parameter value (Parameter *x) const;
  Parameter value (const std::initializer_list<Parameter>& x) const;

  static PPoly2 one() {
    int nt = 1;
    Parameter a[1];
    a[0] = Parameter::constant(1);
    int m[] = {0, 0};
    return PPoly2(nt, a, m);
  }

  static PPoly2 zero() {
    int nt = 1;
    Parameter a[1];
    a[0] = Parameter::constant(0);
    int m[] = {0, 0};
    return PPoly2(nt, a, m);
  }

  PPoly1 subX(Parameter x) const;
  PPoly1 subY(Parameter y) const;

  PPoly2 derX () const;
  PPoly2 derY () const;

  int degX() const;
  int degY() const;

  void print_input() const;
  void print() const;

  PV2 polish (const PV2 &x) const;

  PPoly2 operator* (const PPoly2 &g) const;
  PPoly2 operator* (const Parameter &p) const;
  PPoly2 operator* (const double p) const;
  PPoly2 plusCtimes (double c, const PPoly2 &b) const;
  PPoly2 operator+ (const PPoly2 &b) const { return plusCtimes(1.0, b); }
  PPoly2 operator- (const PPoly2 &b) const { return plusCtimes(-1.0, b); }

  int size () { return nt; }
  Parameter &operator[] (int i) { return a[i]; }

  int nt, d;
  Parameter *a;
  int *m;
  int degx, degy;

};

class PPoly3 {
public:
  PPoly3 () {}

  void add (int i, int j, int k, const Parameter& c) { 
    auto it = ia.find(Index3(i, j, k));
    if (it == ia.end()) {
      ia[Index3(i, j, k)] = a.size();
      a.push_back(c);
    }
    else
      a[it->second] += c;
  }

  Parameter value
    (const Parameter& x, const Parameter& y, const Parameter& z) const;

  PV3 gradient
    (const Parameter& x, const Parameter& y, const Parameter& z) const;

  // xyz = 0, 1, 2:  substituting for yz, xz, xy.
  PPoly<Parameter> substitute2
    (int xyz, const Parameter& x, const Parameter& y, const Parameter& z) const;

  class Index3 {
  public:
    int i[3];
    Index3 (int i0, int i1, int i2) { i[0] = i0; i[1] = i1; i[2] = i2; }
  };  

  class Hasher {
  public:
    size_t operator() (const Index3 &a) const {
      return a.i[0] * 779230947 + a.i[1] * 247091631 + a.i[2] * 1194289623;
    }
  };

  class Equals {
  public:
    bool operator() (const Index3 &a, const Index3 &b) const {
      return a.i[0] == b.i[0] && a.i[1] == b.i[1] && a.i[2] == b.i[2];
    }
  };
  
  unordered_map<Index3, int, Hasher, Equals> ia;
  vector<Parameter> a;
  int size () { return a.size(); }
  Parameter &operator[] (int i) { return a[i]; }
};

}

#endif
