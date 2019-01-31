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

#ifndef ACP_H
#define ACP_H

//#define NO_MODE

#include <gmp.h>
#include <mpfr.h>
#include <assert.h>
#include <iostream>
#include <fenv.h>
#include <exception>
#include <math.h>
#include <algorithm>
#include <float.h>

//#define ROOT_SEPARATION
#ifdef ROOT_SEPARATION
extern int nGCD106;
extern int nGCD212;
extern int nRootCmps;
extern bool inGCD;
#endif

//#define MOD_PRIME26
#ifdef MOD_PRIME26
extern int prime26;
#endif
extern bool inGetApprox;
namespace acp {

double randomNumber (double rmin, double rmax);

inline double nextD (double x) {
  static double p = 1.0/(1 << 26)/(1 << 26);
  double x2 = x + p*fabs(x);
  if (x2 == x)
    return nextafter(x, 1.0);
  double x3 = 0.5*(x + x2);
  return x3 > x ? x3 : x2;
}
 
inline double prevD (double x) {
  static double p = 1.0/(1 << 26)/(1 << 26);
  double x2 = x - p*fabs(x);
  if (x2 == x)
    return nextafter(x, -1.0);
 double x3 = 0.5*(x + x2);
 return x3 < x ? x3 : x2;
}

class SignException : public std::exception {
 public:
  virtual const char* what() const throw() {
    return "Not enough precision";
  }
};

class PrecisionException : public std::exception {
 public:
  virtual const char* what() const throw() {
    return "Maximum precision exceeded";
  }
};

class IdentityException : public std::exception {
 public:
  virtual const char* what() const throw() {
    return "Identically zero primitive";
  }
};

extern SignException signException;
extern PrecisionException precisionException;
extern IdentityException identityException;

enum RoundMode { RoundUp=1, RoundDown=-1, RoundNearest=0 };

class MValue {
 public:
  MValue (unsigned int ip) : p(ip) { mpfr_init2(m, p); }
  MValue (double x, unsigned int ip) : p(ip) { 
    mpfr_init2(m, p);
    mpfr_set_d(m, x, GMP_RNDN);
  }
  MValue (const MValue &v, unsigned int ip) : p(ip) { 
    mpfr_init2(m, ip); 
    mpfr_set(m, v.m, GMP_RNDN); 
  }
  MValue (const MValue &v) : p(v.p) { 
    mpfr_init2(m, v.p); 
    mpfr_set(m, v.m, GMP_RNDN); 
  }
  ~MValue () { mpfr_clear(m); }
  void operator= (const MValue &v) { mpfr_set(m, v.m, GMP_RNDN); }
  double value () const { return mpfr_get_d(m, GMP_RNDN); }
  MValue plus (const MValue &b, mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_add(res.m, m, b.m, round);
    return res;
  }
  MValue plus (double b, mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_add_d(res.m, m, b, round);
    return res;
  }
  MValue minus () const {
    MValue res(p);
    mpfr_neg(res.m, m, GMP_RNDN);
    return res;
  }
  MValue minus (const MValue &b, mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_sub(res.m, m, b.m, round);
    return res;
  }
  MValue times (const MValue &b, mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_mul(res.m, m, b.m, round);
    return res;
  }
  MValue times (double b, mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_mul_d(res.m, m, b, round);
    return res;
  }
  MValue divide (const MValue &b, mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_div(res.m, m, b.m, round);
    return res;
  }
  int sign () const { return mpfr_sgn(m); }
  bool operator< (const MValue &b)  const { return mpfr_less_p(m, b.m); }
  MValue sqrt (mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_sqrt(res.m, m, round);
    return res;
  }
  MValue root (unsigned long int k, mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_root(res.m, m, k, round);
    return res;
  }
  MValue cos (mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_cos(res.m, m, round);
    return res;
  }
  MValue acos (mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_acos(res.m, m, round);
    return res;
  }

  mpfr_t m;
  unsigned int p;
};

typedef unsigned int uint;

typedef unsigned long ulong;

unsigned int inverse (unsigned int a, unsigned int n);

#define ulong(a) ((unsigned long)(a))
 
class ModP {
  static const int eShift;
  static const int eMax;
  static const int eMin;

  unsigned int *pow2;

public:
  const unsigned int p;

  ModP (unsigned int p) : p(p) {
    pow2 = new unsigned int[eMax - eMin + 1] - eMin;
    pow2[0] = 1;
    unsigned long x = 1;
    for (int e = 1; e <= eMax; e++) {
      x = (2 * x) % p;
      pow2[e] = x;
    }
    x = 1;
    unsigned int i2 = (p + 1) / 2;
    for (int e = -1; e >= eMin; e--) {
      x = (i2 * x) % p;
      pow2[e] = x;
    }
  }

  ~ModP () { delete[] (pow2 + eMin); }

  unsigned int mod (double x) {
    if (x == (long) x)
      return x >= 0 ? ((long) x) % p : p + ((long) x) % p;
    int e;
    double m = frexp(x, &e) * (1l << eShift);
    assert(m == (long) m);
    e -= eShift;
    long mp = ((long) m) % p;
    if (mp < 0)
      mp += p;
    assert(eMin <= e && e <= eMax);
    return ((unsigned long) mp * pow2[e]) % p;
  }
};

class ModInt {
 public:
  static unsigned int p1, p2;
  static ModP modP1, modP2;

  ModInt () {}
  
  ModInt (double r, unsigned int ip) : p(ip) {
    if (p == p1)
      a = modP1.mod(r);
    else if (p == p2)
      a = modP2.mod(r);
    else
      assert(0);
    /*
    long l = r;
    if (l == r) {
      l = l%p;
      a = l < 0 ? l + p : l;
    }
    else {
      union {
	double r;
	unsigned long l;
      } u;
      u.r = r;
      a = u.l%p;
    }
    */
  }

  ModInt operator- () const {
    long l = (-long(a))%p;
    return ModInt(l < 0 ? l + p : l, p);
  }
  
  ModInt operator+ (const ModInt &x) const {

    return ModInt((long(a) + long(x.a))%p, p);
  }

  ModInt operator+ (double x) const {
    return *this + ModInt(x, p);
  }

  ModInt operator- (const ModInt &x) const {
    long l = (long(a) - long(x.a))%p;
    return ModInt(l < 0 ? l + p : l, p);
  }

  ModInt operator- (double x) const {
    return *this - ModInt(x, p);
  }

  ModInt operator* (const ModInt &x) const {
    return ModInt((ulong(a)*ulong(x.a))%p, p);
  }

  ModInt operator* (double x) const {
    return *this*ModInt(x, p);
  }
  
  ModInt operator/ (const ModInt &x) const {
    return ModInt((ulong(a)*ulong(inverse(x.a, p)))%p, p);
  }

  unsigned int a, p;
};

class Parameter;

class MInt {
protected:
  unsigned int refCnt;

  MInt () : refCnt(1), m1(ModInt(random(), ModInt::p1)), m2(ModInt(random(), ModInt::p2)) {}
  MInt (const ModInt &m1, const ModInt &m2) : refCnt(1), m1(m1), m2(m2) {}
  virtual ~MInt () {}

public:
  const ModInt m1, m2;

  void incRef () { refCnt++; }
  void decRef () { if (--refCnt == 0) delete this; }

  static MInt* make (int precision, const Parameter &p);
  static MInt* make (int precision, const MInt *m);

  virtual unsigned int precision () const = 0;
  virtual double intervalWidth () const = 0;
  virtual double lb () const = 0;
  virtual double ub () const = 0;
  virtual MInt* lbP () const = 0;
  virtual MInt* ubP () const = 0;
  virtual MInt* plus (const MInt* b) const = 0;
  virtual MInt* plus (double b) const = 0;
  virtual MInt* minus (const MInt* b) const = 0;
  virtual MInt* minus () const = 0;
  virtual MInt* times (const MInt* b) const = 0;
  virtual MInt* times (double b) const = 0;
  virtual MInt* divide (const MInt* b) const = 0;
  virtual int sign () const = 0;
  virtual MInt* mid () const = 0;
  virtual bool subset (const MInt* b) const = 0;
  virtual MInt* interval (const MInt* b) const = 0;
  virtual MInt* innerInterval (const MInt* b) const = 0;
  virtual bool intersects (const MInt* b) const = 0;
  virtual MInt* intersect (const MInt* b) const = 0;
  virtual MInt* sqrt () const = 0;
  virtual MInt* root (unsigned long int k) const = 0;
  virtual MInt* cos () const = 0;
  virtual MInt* acos () const = 0;
};

class AnglePoly;
class PFun;
class PPoly1;
class PPoly2;

// The number class:  an interval of double with multiple precision backup.
// Treat it like a double.
class Parameter {
  friend class MValue;
  template<class P> friend class Object;
  friend class Primitive;
  friend class Poly;
  friend class AnglePoly;
  friend class PPoly1;
  friend class PPoly2;
  friend class Encasement;
  friend class PInt;

  // static double force_rounding (double x) { volatile double x_ = x; return x_; }
  static double no_optimize (volatile double x) { return x; }

 public:

  static unsigned long int identityCount;

  int size () { return 1; }
  Parameter &operator[] (int i) { return *this; }

  // Perturbation magnitude.
  static double delta; // default is 2^{-26}

  // Maximum allowed precision.  Default is 848 bits.
  static unsigned int maxPrecision;

  // Handle sign exception.  Default is true.  If false, programmer
  // must catch signException.
  static bool handleSignException;

  // throw precision exception. Default is true.  If false, return 0
  static bool usePrecisionException;

  // sign returns 0 for identity if true (the default).  If false, sign
  // throws an IdentityException
  static bool handleIdentity;

  // Call Parameter::enable() before using ACP
  static void enable () { 
#ifndef NO_MODE
    if (!enabled) { enabled = true; fesetround(FE_UPWARD); }
#endif
  }

  // Call Parameter::disable() before non-ACP calculations.
  static void disable () { 
#ifndef NO_MODE
    if (enabled) { enabled = false; fesetround(FE_TONEAREST); }
#endif
  }

  static bool isEnabled () { return enabled; }

  Parameter () : l(0.0) { u.r = -1.0; }

  bool uninitialized () const { return l == 0.0 && u.r == -1.0; }
  
  static Parameter input (double x) { 
    double y = x + delta*(1.0 + fabs(x))*randomNumber(-1.0, 1.0);
    Parameter p(y, y);
    if (highPrecision > 53u)
      p.increasePrecision();
    return p;
  }

  static Parameter constant (double x) { 
    Parameter p(x, x); 
    if (highPrecision > 53u)
      p.increasePrecision();
    return p;
  }

  static Parameter interval (double x, double y) { 
    Parameter p(x, y);
    if (highPrecision > 53u)
      p.increasePrecision();
    return p;
  }

  // x is not perturbed
  Parameter (double x) {
    l = u.r = x;
    if (highPrecision > 53u)
      increasePrecision();
  }

  Parameter (const Parameter &p) : l(p.l) {
    if (l != sentinel)
      u.r = p.u.r;
    else {
      u.m = p.u.m;
      u.m->incRef();
    }
  }
  
  ~Parameter () { 
    if (l == sentinel)
      u.m->decRef();
  }
  
  void intersectIntervals (const Parameter &p);

  Parameter & operator= (const Parameter &p) {
    if (p.l == sentinel) {
      p.u.m->incRef();
      if (l != sentinel && l <= u.r && p.u.m->precision() == 106u)
	intersectIntervals(p);
    }
    if (l == sentinel)
      u.m->decRef();
    l = p.l;
    u = p.u;
    return *this;
  }
  
  // Determine the sign of a parameter: -1 or 1.
  // If fail==true, will return 0 on sign failure.
  // Otherwise throws SignException.
  int sign (bool fail = true) const {
    if (l != sentinel) {
      if (l > 0.0)
	return 1;
      if (u.r < 0.0)
	return -1;
      if (handleIdentity && l == 0 && u.r == 0)
	return 0;
      if (fail) {
	highPrecision = 106u;
	throw signException;
      }
      return 0;
    }
    int s = u.m->sign();
    if (s)
      return s;
    if (u.m->m1.a == 0u && u.m->m2.a == 0u) {
      ++identityCount;
      if (handleIdentity)
	return 0;
      throw identityException;
    }
    if (fail) {
      if (highPrecision <= maxPrecision) {
        highPrecision *= 2;
        throw signException;
      }
      else if (usePrecisionException)
        throw precisionException;
      else
        return 0;
    }
    return 0;
  }
  
  bool operator== (int zero) {
    assert(zero == 0);
    return sign() == 0;
  }

  // Approximate value of parameter for display purposes.
  double mid () const { return 0.5*(lb() + ub()); }

  // double precision lower bound of interval
  double lb () const { 
    return l != sentinel ? l : u.m->lb(); 
  }

  // double precision upper bound of interval
  double ub () const { 
    return l != sentinel ? u.r : u.m->ub(); 
  }

  // width of interval
  double intervalWidth () const {
    return l != sentinel ? u.r - l : u.m->intervalWidth();
  }
  
  // Parameter has unary - and binary +,-,*,/, <, >
  // If one argument of binary operation is a double instead of a
  // Parameter, it is treated as a constant, not perturbed.
  // Hence 2 * p does not perturb the 2.
  Parameter operator+ (const Parameter &b) const {
    if (l != sentinel && b.l != sentinel)
#ifndef NO_MODE
      return Parameter(- no_optimize(no_optimize(- l) - b.l), u.r + b.u.r);
#else
      return Parameter(prevD(l + b.l), nextD(u.r + b.u.r));
#endif    
    assert(l == sentinel && b.l == sentinel);
    return Parameter(u.m->plus(b.u.m));
  }

  Parameter & operator+= (const Parameter &p) {
    return *this = *this + p;
  }

  Parameter operator+ (double b) const { return *this + Parameter(b); }
  
  Parameter operator- (const Parameter &b) const {
    if (l != sentinel && b.l != sentinel)
#ifndef NO_MODE
      return Parameter(- no_optimize(b.u.r - l), u.r - b.l);
#else
      return Parameter(prevD(l - b.u.r), nextD(u.r - b.l));
#endif    
    assert(l == sentinel && b.l == sentinel);
    return Parameter(u.m->minus(b.u.m));
  }
  
  Parameter operator- (double b) const { return *this - Parameter(b); }

  Parameter operator- () const {
    if (l != sentinel)
      return Parameter(- u.r, - l);
    return Parameter(u.m->minus());
  }
  
  Parameter operator* (const Parameter &b) const {
#ifndef NO_MODE
    if (l != sentinel && b.l != sentinel) {
      Parameter s = u.r < 0.0 ? - *this : *this, t = u.r < 0.0 ? - b : b;
      if (s.l > 0.0) {
	volatile double k = t.l > 0.0 ? - s.l : - s.u.r;
	return Parameter(- no_optimize(k*t.l), t.u.r > 0.0 ? s.u.r*t.u.r : s.l*t.u.r);
      }
      if (t.l > 0.0) {
	volatile double k = - s.l;
	return Parameter(- no_optimize(k*t.u.r), s.u.r*t.u.r);
      }
      if (t.u.r < 0.0) {
	volatile double k = - s.u.r;
	return Parameter(- no_optimize(k*t.l), s.l*t.l);
      }
      volatile double k1 = - s.l, k2 = - s.u.r;
      double cl1 = no_optimize(k1*t.u.r), cl2 = no_optimize(k2*t.l), 
	cu1 = s.l*t.l, cu2 = s.u.r*t.u.r,
	cl = cl1 < cl2 ? - cl2 : - cl1, 
        cu = cu1 < cu2 ? cu2 : cu1;
      return Parameter(cl, cu);
    }
#else
    if (l != sentinel && b.l != sentinel) {
      Parameter s = u.r < 0.0 ? - *this : *this, t = u.r < 0.0 ? - b : b;
      if (s.l >= 0.0)
	if (t.l >= 0.0)
	  return Parameter(prevD(s.l*t.l), nextD(s.u.r*t.u.r));
	else if (t.u.r <= 0.0)
	  return Parameter(prevD(s.u.r*t.l), nextD(s.l*t.u.r));
	else
	  return Parameter(prevD(s.u.r*t.l), nextD(s.u.r*t.u.r));
      if (t.l >= 0.0)
	return Parameter(prevD(s.l*t.u.r), nextD(s.u.r*t.u.r));
      if (t.u.r <= 0.0)
	return Parameter(prevD(s.u.r*t.l), nextD(s.l*t.l));
      double k1 = s.l*t.u.r, k2 = s.u.r*t.l, nl = k1 < k2 ? k1 : k2,
	k3 = s.l*t.l, k4 = s.u.r*t.u.r, nu = k3 < k4 ? k4 : k3;
      return Parameter(prevD(nl), nextD(nu));
    }
#endif
    assert(l == sentinel && b.l == sentinel);
    return Parameter(u.m->times(b.u.m));
  }
  
  Parameter operator* (double b) const { return *this * Parameter(b); }

  Parameter operator/ (const Parameter &b) const {
    int bs = b.sign();
    assert(bs != 0);
#ifndef NO_MODE
    if (l != sentinel && b.l != sentinel) {
      if (bs == 1) {
	if (l >= 0.0)
	  return Parameter(- no_optimize(no_optimize(- l)/b.u.r), u.r/b.l);
	if (u.r <= 0.0)
	  return Parameter(- no_optimize(no_optimize(- l)/b.l), u.r/b.u.r);
	return Parameter(- no_optimize(no_optimize(- l)/b.l), u.r/b.l);
      }
      if (l >= 0.0)
	return Parameter(- no_optimize(no_optimize(- u.r)/b.u.r), l/b.l);
      if (u.r <= 0.0)
	return Parameter(- no_optimize(no_optimize(- u.r)/b.l), l/b.u.r);
      return Parameter(- no_optimize(no_optimize(- u.r)/b.u.r), l/b.u.r);
    }
#else
    if (l != sentinel && b.l != sentinel) {
      if (bs == 1)
	if (l >= 0.0)
	  return Parameter(prevD(l/b.u.r), nextD(u.r/b.l));
	else if (u.r <= 0.0)
	  return Parameter(prevD(l/b.l), nextD(u.r/b.u.r));
	else
	  return Parameter(prevD(l/b.l), nextD(u.r/b.l));
      if (l >= 0.0)
	return Parameter(prevD(u.r/b.u.r), nextD(l/b.l));
      if (u.r <= 0.0)
	return Parameter(prevD(u.r/b.l), nextD(l/b.u.r));
      return Parameter(prevD(u.r/b.u.r), nextD(l/b.u.r));
    }
#endif
    assert(l == sentinel && b.l == sentinel);
    return Parameter(u.m->divide(b.u.m));
  }

  Parameter operator/ (double b) const { return *this / Parameter(b);}
    
  bool operator< (const Parameter &b) const { return (b - *this).sign() == 1; }
  bool operator< (double b) const { return (*this - b).sign() == -1; }
  bool operator> (const Parameter &b) const { return (b - *this).sign() == -1; }
  bool operator> (double b) const { return (*this - b).sign() == 1; }

  static Parameter max (const Parameter &a, const Parameter &b) { return a < b ? b : a; }

  Parameter abs () const { 
    return sign() == 1 ? *this : -*this; 
  }

  Parameter sqrt () const;

  Parameter root (unsigned long int n) const;

  Parameter pow (long int n) const {
    if (n == 0)
      return Parameter(1.0);
    if (n < 0)
      return Parameter(1.0)/Parameter::pow(*this, -n);
    return Parameter::pow(*this, n);
  }

  Parameter cos () const {
    if (highPrecision == 106u) {
      highPrecision = 212u;
      throw signException;
    }
    if (l != sentinel) {
      bool save_enabled = isEnabled();
      if (save_enabled)
        disable();
      double c1 = ::cos(l), c2 = ::cos(u.r), cl, cu;
      if (c1 < c2) {
	cl = nextafter(c1, -10.0);
	cu = nextafter(c2, 10.0);
      }
      else {
	cl = nextafter(c2, -10.0);
	cu = nextafter(c1, 10.0);
      }
      if (save_enabled)
        enable();
      return Parameter(cl, cu);
    }
    return Parameter(u.m->cos());
  }

  Parameter acos () const {
    if (highPrecision == 106u) {
      highPrecision = 212u;
      throw signException;
    }
    if (l != sentinel) {
      bool save_enabled = isEnabled();
      if (save_enabled)
        disable();
      double a1 = ::acos(l < -1.0 ? -1.0 : l),
	a2 = ::acos(u.r > 1.0 ? 1.0 : u.r), al, au;
      if (a1 < a2) {
      	al = nextafter(a1, -10.0);
	au = nextafter(a2, 10.0);
      }
      else {
	al = nextafter(a2, -10.0);
	au = nextafter(a1, 10.0);
      }
      if (save_enabled)
        enable();
      return Parameter(al, au);
    }
    return Parameter(u.m->acos());
  }

  void print ();

  void print (FILE *out);

 private:
  static Parameter sqrt (double);
  static Parameter root (double a, unsigned int n);
  static Parameter pow (const Parameter &x, unsigned long int n);

  static const double sentinel;

 public:
  static unsigned int highPrecision;
 private:

  static bool enabled;

  Parameter (double il, double iu) : l(il) {
    u.r = iu;
    // increasePrecision();
#ifdef MOD_PRIME26
    if (l == (long) l)
      l = (double) (((long) l) % prime26);
    if (u.r == (long) u.r)
      u.r = (double) (((long) u.r) % prime26);
#endif
  }

  Parameter (MInt *m) : l(sentinel) { u.m = m; }

 public:
  Parameter lbP () const {
    if (l != sentinel)
      return Parameter(l, l);
    return Parameter(u.m->lbP());
  }

  Parameter ubP () const {
    if (l != sentinel)
      return Parameter(u.r, u.r);
    return Parameter(u.m->ubP());
  }

  Parameter (const Parameter &il, const Parameter &iu) {
    if (il.l != sentinel) {
      l = il.l;
      u.r = iu.u.r;
    }
    else {
      l = sentinel;
      u.m = il.u.m->interval(iu.u.m);
    }
  }

  unsigned int precision () const {
    return l != sentinel ? 53u : u.m->precision();
  }

  bool increased () const { return precision() == highPrecision; }
  bool decreased () const { return precision() == 53u; }

  void increasePrecision () {
    assert(!uninitialized());
    if (increased())
      return;
    if (l != sentinel) {
      u.m = MInt::make(highPrecision, *this);
      l = sentinel;
    }
    else {
      MInt* m = MInt::make(highPrecision, u.m);
      u.m->decRef();
      u.m = m;
    }
    if (highPrecision >= 212u)
      disable();
  }

  void decreasePrecision () {
    Parameter::highPrecision = 53u;
    if (l == sentinel) {
      enable();
      double dl = lb(), du = ub();
      u.m->decRef();
      l = dl;
      u.r = du;
    }
  }

  Parameter midP () const {
    if (l != sentinel) {
      double m = 0.5*(l + u.r);
      return Parameter(m, m);
    }
    return Parameter(u.m->mid());
  }

  bool subset (const Parameter &b) const {
    if (l != sentinel && b.l != sentinel)
      return !(l < b.l) && !(b.u.r < u.r) && (b.l < l || u.r < b.u.r);
    assert (l == sentinel && b.l == sentinel);
    return u.m->subset(b.u.m);
  }

  Parameter interval (const Parameter &b) const {
    if (l != sentinel && b.l != sentinel)
      return Parameter(l, b.u.r);
    assert (l == sentinel && b.l == sentinel);
    return Parameter(u.m->interval(b.u.m));
  }

  Parameter innerInterval (const Parameter &b) const {
    if (l != sentinel && b.l != sentinel)
      return Parameter(u.r, b.l);
    assert (l == sentinel && b.l == sentinel);
    return Parameter(u.m->innerInterval(b.u.m));
  }

  bool intersects(const Parameter &b) const {
    if(l != sentinel && b.l != sentinel) {
      return !(u.r < b.l || b.u.r < l);
    } 
    assert(l == sentinel && b.l == sentinel);
    return u.m->intersects(b.u.m);
  }

  Parameter intersect (const Parameter &b) const {
    if (l != sentinel && b.l != sentinel) {
      assert(!(u.r < b.l || b.u.r < l));
      double il = l < b.l ? b.l : l;
      double iu = u.r < b.u.r ? u.r : b.u.r;
      return Parameter(il, iu);
    }
    assert (l == sentinel && b.l == sentinel);
    return Parameter(u.m->intersect(b.u.m));
  }

  double l;
  union {
    double r;
    MInt *m;
  } u;
};

inline Parameter operator+ (double a, const Parameter &b)
{
  return b + a;
}

inline Parameter operator- (double a, const Parameter &b)
{
  return (- b) + a;
}

inline Parameter operator* (double a, const Parameter &b)
{
  return b*a;
}

inline Parameter operator/ (double a, const Parameter &b)
{
  return Parameter(a)/b;
}

inline bool operator< (double a, const Parameter &b)
{
  return (b - a).sign() == 1;
}

inline bool operator> (double a, const Parameter &b)
{
  return (b - a).sign() == -1;
}

class PInt : public MInt {
  Parameter p;
  friend class Parameter;
 public:
  PInt (const Parameter &p) { this->p = p;
    assert(p.l != Parameter::sentinel);
  }

  PInt (const Parameter &p, const ModInt &m1, const ModInt &m2) 
    : MInt(m1, m2) { this->p = p;
    assert(p.l != Parameter::sentinel);
  }
  
  unsigned int precision () const { return 106u; }

  double intervalWidth () const { return p.intervalWidth(); }
  
  double lb () const { return p.lb(); }
  
  double ub () const { return p.ub(); }
  
  MInt* lbP () const { return new PInt(p.lbP()); }

  MInt* ubP () const { return new PInt(p.ubP()); }

  MInt* plus (const MInt* that) const {
    const PInt &b = *dynamic_cast<const PInt*>(that);
    return new PInt(p + b.p, m1 + b.m1, m2 + b.m2);
  }
  
  MInt* plus (double b) const { return new PInt(p + b, m1 + b, m2 + b); }
  
  MInt* minus (const MInt* that) const {
    const PInt &b = *dynamic_cast<const PInt*>(that);
    return new PInt(p - b.p, m1 - b.m1, m2 - b.m2);
  }
  
  MInt* minus () const { return new PInt(- p, - m1, - m2); }
  
  MInt* times (const MInt* that) const {
    const PInt &b = *dynamic_cast<const PInt*>(that);
    return new PInt(p * b.p, m1*b.m1, m2*b.m2);
  }
  
  MInt* times (double b) const { return new PInt(p * b, m1*b, m2*b); }

  MInt* divide (const MInt* that) const {
    const PInt &b = *dynamic_cast<const PInt*>(that);
    return new PInt(p / b.p, m1/b.m1, m2/b.m2); 
  }

  int sign () const {
    if (p.l > 0)
      return 1;
    if (p.u.r < 0)
      return -1;
    return 0;
  }

  MInt* mid () const { return new PInt(p.midP()); }

  bool subset (const MInt* that) const {
    const PInt &b = *dynamic_cast<const PInt*>(that);
    return p.subset(b.p);
  }

  MInt* interval (const MInt* that) const {
    const PInt &b = *dynamic_cast<const PInt*>(that);
    return new PInt(p.interval(b.p));
  }

  MInt* innerInterval (const MInt* that) const {
    const PInt &b = *dynamic_cast<const PInt*>(that);
    return new PInt(p.innerInterval(b.p));
  }

  bool intersects (const MInt* that) const {
    const PInt &b = *dynamic_cast<const PInt*>(that);
    return p.intersects(b.p);
  }

  MInt* intersect (const MInt* that) const {
    const PInt &b = *dynamic_cast<const PInt*>(that);
    return new PInt(p.intersect(b.p));
  }

  MInt* sqrt () const { assert(0); return new PInt(p.sqrt()); }

  MInt* root (unsigned long int k) const { assert(0); return new PInt(p.root(k)); }

  MInt* cos () const { assert(0); }

  MInt* acos () const { assert(0); }
};

inline void Parameter::intersectIntervals (const Parameter &p) {
  Parameter &pp = dynamic_cast<PInt*>(p.u.m)->p;
  assert(!(u.r < pp.l || pp.u.r < l));
  pp.l = l > pp.l ? l : pp.l;
  pp.u.r = u.r < pp.u.r ? u.r : pp.u.r;
}

class EInt : public MInt {
 public:
  EInt (const MValue &l, const MValue &u) : lm(l), um(u) {}
  EInt (const MValue &l, const MValue &u, const ModInt &m1, const ModInt &m2) 
    : MInt(m1, m2), lm(l), um(u) {}
  ~EInt () {}

  unsigned int precision () const { return lm.p; }

  double intervalWidth () const { return um.minus(lm, GMP_RNDN).value(); }

  double lb () const { return mpfr_get_d(lm.m, GMP_RNDD); }

  double ub () const { return mpfr_get_d(um.m, GMP_RNDU); }

  MInt* lbP () const { return new EInt(lm, lm); }

  MInt* ubP () const { return new EInt(um, um); }

  MInt* plus (const MInt* that) const {
    const EInt &b = *dynamic_cast<const EInt*>(that);
    return new EInt(lm.plus(b.lm, GMP_RNDD), um.plus(b.um, GMP_RNDU),
		    m1 + b.m1, m2 + b.m2);
  }

  MInt* plus (double b) const {
    return new EInt(lm.plus(b, GMP_RNDD), um.plus(b, GMP_RNDU),
		    m1 + b, m2 + b);
  }

  MInt* minus (const MInt* that) const {
    const EInt &b = *dynamic_cast<const EInt*>(that);
    return new EInt(lm.minus(b.um, GMP_RNDD), um.minus(b.lm, GMP_RNDU),
		    m1 - b.m1, m2 - b.m2);
  }

  MInt* minus () const { 
    return new EInt(um.minus(), lm.minus(), - m1, - m2);
  }

  MInt* times (const MInt* that) const;

  MInt* times (double b) const {
    return b > 0.0 ? new EInt(lm.times(b, GMP_RNDD), um.times(b, GMP_RNDU),
			      m1*b, m2*b)
      : new EInt(um.times(b, GMP_RNDD), lm.times(b, GMP_RNDU),
		 m1*b, m2*b);
  }

  MInt* divide (const MInt* b) const;

  int sign () const {
    if (lm.sign() == 1) return 1;
    if (um.sign() == -1) return -1;
    return 0;
  }

  MInt* mid () const {
    MValue m = lm.plus(um, GMP_RNDN).times(0.5, GMP_RNDN);
    return new EInt(m, m);
  }

  bool subset (const MInt* that) const {
    const EInt &b = *dynamic_cast<const EInt*>(that);
    return !(lm < b.lm) && !(b.um < um) && (b.lm < lm || um < b.um);
  }

  MInt* interval (const MInt* that) const {
    const EInt &b = *dynamic_cast<const EInt*>(that);
    return new EInt(lm, b.um);
  }

  MInt* innerInterval (const MInt* that) const {
    const EInt &b = *dynamic_cast<const EInt*>(that);
    return new EInt(um, b.lm);
  }

  bool intersects (const MInt* that) const {
    const EInt &b = *dynamic_cast<const EInt*>(that);
    return !(um < b.lm || b.um < lm);
  }

  MInt* intersect (const MInt* that) const {
    const EInt &b = *dynamic_cast<const EInt*>(that);
    assert(!(um < b.lm || b.um < lm));
    MValue l = lm < b.lm ? b.lm : lm;
    MValue u = um < b.um ? um : b.um;
    return new EInt(l, u);
  }

  MInt* sqrt () const {
    return new EInt(lm.sqrt(GMP_RNDD), um.sqrt(GMP_RNDU));
  }

  MInt* root (unsigned long int k) const {
    return new EInt(lm.root(k, GMP_RNDD), um.root(k, GMP_RNDU));
  }

  MInt* cos () const {
    MValue v[] = {lm.cos(GMP_RNDD), lm.cos(GMP_RNDU), um.cos(GMP_RNDD),
		  um.cos(GMP_RNDU)}, vmin = v[0], vmax = v[0];
    for (int i = 1; i < 4; ++i)
      if (v[i] < vmin)
	vmin = v[i];
      else if (vmax < v[i])
	vmax = v[i];
    return new EInt(vmin, vmax);
  }

  MInt* acos () const {
    return new EInt(um.acos(GMP_RNDD), lm.acos(GMP_RNDU));
  }

  MValue lm, um;
};

inline void Parameter::print () {
  assert(intervalWidth() == 0);
  if (l == sentinel)
    mpfr_printf("%.0Rf", dynamic_cast<EInt*>(u.m)->lm.m);
}

inline void Parameter::print (FILE *out) {
  assert(intervalWidth() == 0);
  if (l == sentinel)
    mpfr_fprintf(out, "%.0Rf", dynamic_cast<EInt*>(u.m)->lm.m);
}

inline MInt* MInt::make (int precision, const Parameter &p) {
  double l = p.lb(), u = p.ub();
  ModInt m1 = l == u ? ModInt(l, ModInt::p1) : ModInt(random(), ModInt::p1);
  ModInt m2 = l == u ? ModInt(l, ModInt::p2) : ModInt(random(), ModInt::p2);
  if (precision == 106u)
    return new PInt(p, m1, m2);
  else
    return new EInt(MValue(l, precision), MValue(u, precision), m1, m2);
}

inline MInt* MInt::make (int precision, const MInt *m) {
  if (m->precision() == 106u)
    return new EInt(MValue(m->lb(), precision), MValue(m->ub(), precision), 
                    m->m1, m->m2);
  else {
    const EInt* e = dynamic_cast<const EInt*>(m);
    return new EInt(MValue(e->lm, precision), MValue(e->um, precision),
                    m->m1, m->m2);
  }
}

}
#endif
