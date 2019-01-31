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

#include "acp.h"
using namespace acp;

bool inGetApprox;

namespace acp {

double randomNumber (double rmin, double rmax)
{
  return rmin + (rmax - rmin)*random()/double(RAND_MAX);
}

unsigned int inverse (unsigned int a, unsigned int n)
{
  long t = 0, nt = 1, r = n, nr = a;
  while (nr != 0) {
    long q = r/nr, pt = t, pr = r;
    t = nt;
    nt = pt - q*nt;
    r = nr;
    nr = pr - q*nr;
  }
  return t < 0 ? t + n : t;
}

const int ModP::eShift = 53;
const int ModP::eMax = 1024 - eShift;
const int ModP::eMin = -1073 - eShift;

unsigned int ModInt::p1 = 4294967291u;
unsigned int ModInt::p2 = 4294967279u;
ModP ModInt::modP1(p1);
ModP ModInt::modP2(p2);

MInt* EInt::times (const MInt* that) const
{
  const EInt &b = *dynamic_cast<const EInt*>(that);
  bool pflag = um.sign() == -1;
  MValue sl = pflag ? um.minus() : lm, su = pflag ? lm.minus() : um,
    tl = pflag ? b.um.minus() : b.lm, tu = pflag ? b.lm.minus() : b.um;
  if (sl.sign() == 1) {
    MValue &l1 = tl.sign() == 1 ? sl : su, &u1 = tu.sign() == 1 ? su : sl;
    return new EInt(l1.times(tl, GMP_RNDD), u1.times(tu, GMP_RNDU),
		    m1*b.m1, m2*b.m2);
  }
  if (tl.sign() == 1)
    return new EInt(sl.times(tu, GMP_RNDD), su.times(tu, GMP_RNDU),
		    m1*b.m1, m2*b.m2);
  if (tu.sign() == -1)
    return new EInt(su.times(tl, GMP_RNDD), sl.times(tl, GMP_RNDU),
		    m1*b.m1, m2*b.m2);
  MValue cl1 = sl.times(tu, GMP_RNDD), cl2 = su.times(tl, GMP_RNDD),
    cu1 = sl.times(tl, GMP_RNDU), cu2 = su.times(tu, GMP_RNDU);
  return new EInt(cl1 < cl2 ? cl1 : cl2, cu1 < cu2 ? cu2 : cu1,
		  m1*b.m1, m2*b.m2);
}

MInt* EInt::divide (const MInt* that) const
{
  const EInt &b = *dynamic_cast<const EInt*>(that);
  int as = sign(), bs = b.sign();
  if (bs == 1)
    switch (as) {
    case 1:
      return new EInt(lm.divide(b.um, GMP_RNDD), um.divide(b.lm, GMP_RNDU),
		      m1/b.m1, m2/b.m2);
    case 0:
      return new EInt(lm.divide(b.lm, GMP_RNDD), um.divide(b.lm, GMP_RNDU),
		      m1/b.m1, m2/b.m2);
    case -1:
      return new EInt(lm.divide(b.lm, GMP_RNDD), um.divide(b.um, GMP_RNDU),
		      m1/b.m1, m2/b.m2);
    }
  switch (as) {
  case 1:
    return new EInt(um.divide(b.um, GMP_RNDD), lm.divide(b.lm, GMP_RNDU),
		    m1/b.m1, m2/b.m2);
  case 0:
    return new EInt(um.divide(b.um, GMP_RNDD), lm.divide(b.um, GMP_RNDU),
		    m1/b.m1, m2/b.m2);
  case -1:
    return new EInt(um.divide(b.lm, GMP_RNDD), lm.divide(b.um, GMP_RNDU),
		    m1/b.m1, m2/b.m2);
  }
  return 0;
}

unsigned long int Parameter::identityCount = 0;

double Parameter::delta = ::pow(2.0, -27);
const double Parameter::sentinel = 1e300;
unsigned int Parameter::highPrecision = 53u;
unsigned int Parameter::maxPrecision = 848u;
bool Parameter::handleSignException = true;
bool Parameter::usePrecisionException = true;
bool Parameter::handleIdentity = true;
bool Parameter::enabled = false;
SignException signException;
PrecisionException precisionException;
IdentityException identityException;

Parameter Parameter::sqrt () const {
  if (highPrecision == 106u) {
    highPrecision = 212u;
    throw signException;
  }
  int s = sign();
  assert(s > -1);
  if (l != sentinel)
    return s == 0 ? Parameter(0, 0) : 
      Parameter(Parameter::sqrt(l).l, Parameter::sqrt(u.r).u.r);
  return Parameter(u.m->sqrt());
}

Parameter Parameter::sqrt (double x) {
  assert(highPrecision <= 106u);
  double s = ::sqrt(x);
  Parameter xx(x, x);
  Parameter lr = xx / s;
  if (s < lr.l)
    return Parameter(s, lr.u.r);
  if (s > lr.u.r)
    return Parameter(lr.l, s);
  return lr;
}

Parameter Parameter::pow (const Parameter &x, unsigned long int n) {
  if (n == 1)
    return x;
  Parameter y = pow(x, n / 2);
  if (n %2 == 0)
    return y * y;
  else
    return x * y * y;
}

Parameter Parameter::root (unsigned long int n) const {
  if (highPrecision == 106u)
    throw signException;
  assert(sign() > 0);
  if (l != sentinel)
    return Parameter(Parameter::root(l, n).l, Parameter::root(u.r, n).u.r);
  return Parameter(u.m->root(n));
}

Parameter Parameter::root (double a, unsigned int n) {
  double x = a < 1 ? 1 : a;
  double x_old;

  do {
    x_old = x;
    Parameter xx(x);
    // x = o - f(o) / f'(o)
    // x = o - (o^n - a) / (n o^(n-1))
    // x = o - (o - a / o^(n-1)) / n
    x = (xx - (xx - a / xx.pow(n-1)) / n).u.r;
  } while (x < x_old);

  Parameter xx(x);
  Parameter l = a / xx.pow(n-1);
  return Parameter(l.l, x);
}
}

