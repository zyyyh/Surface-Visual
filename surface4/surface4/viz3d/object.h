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

#ifndef OBJECT_H
#define OBJECT_H

#include "acp.h"
#include <vector>

namespace acp {

class RefCnt {
  template<class T> friend class PTR;
  int refCnt;
  void incRef () { refCnt++; }
  void decRef () { if (--refCnt == 0) delete this; }
public:
  RefCnt () : refCnt(0) {}
  virtual ~RefCnt () { assert(refCnt == 0); }
};

template<class T>
class PTR {
  T *t;
public:
  PTR () : t(0) {}
  PTR (T *t) : t(t) { incRef(); }
  PTR (const PTR &p) : t(p.t) { incRef(); }
  const PTR &operator= (const PTR &p) { 
    p.incRef(); decRef(); t = p.t; return *this; 
  }
  ~PTR () { decRef(); }
  void incRef () const { if (t != 0) t->incRef(); }
  void decRef () const { if (t != 0) t->decRef(); }
  operator T* () const { return t; }
  T *operator-> () const { return t; }
};

class BaseObject {
  friend class Primitive;
  friend class Parameter;
  template<class P> friend class Object;
  template<class P> friend class InputObject;
  friend class Poly;
  friend class AnglePoly;

  static std::vector<BaseObject*> precisionIncreased;

  virtual void decreasePrecision () = 0;

public:
  static void decreaseAll () {
    for (int i = 0; i < precisionIncreased.size(); i++)
      precisionIncreased[i]->decreasePrecision();
    precisionIncreased.clear();
  }

  virtual ~BaseObject () {
    for (int i = precisionIncreased.size(); --i >= 0;)
      if (precisionIncreased[i] == this) {
        //precisionIncreased[i]->decreasePrecision();
        precisionIncreased.erase(precisionIncreased.begin()+i);
      }
  }
};

template<class P> class InputObject;

template<class P>
class Object : BaseObject, public RefCnt {
  friend class InputObject<P>;
  friend class Poly;
  friend class AnglePoly;

  P p;
  
  int size () { return p.size(); }
  Parameter &get (int i) { return p[i]; }

  // The precision of an empty object is the prevailing precision;
  // otherwise, it is the precision of the first (0th) coordinate.
  unsigned int precision () { 
    return size() == 0 ? Parameter::highPrecision : get(0).precision(); 
  }

  // "increased" means this object is at the prevailing precision.
  bool increased () { return precision() == Parameter::highPrecision; }

  void increasePrecision () {
    // std::cout << "increasing precision of " << this << std::endl;
    for (int i = 0; i < size(); i++)
      get(i).increasePrecision();
  }

  void decreasePrecision () {
    // std::cout << "decreasing precision of " << this << std::endl;
    for (int i = 0; i < size(); i++)
      get(i).decreasePrecision();
  }

  virtual P calculate () {
    if (uninitialized())
      std::cerr << "An extension of Object must call set(p) in the constructor "
                << "or override calculate()" << std::endl;
    assert(!uninitialized());
    if (size() > 0 && get(0).lb() < get(0).ub())
      std::cerr << "Calculate is not overridden but coordinates are non-trivial intervals.\n";
    assert(!(size() > 0 && get(0).lb() < get(0).ub()));
    return getCurrentValue();
  }

public:
  // An empty object is assumed to be uninitialized, although the user
  // might really want it to be an empty object.
  bool uninitialized () { return size() == 0 || get(0).uninitialized(); }

  const P &getCurrentValue () {
    assert(!uninitialized());
    int precis = precision();
    if (precis < Parameter::highPrecision) {
      increasePrecision();
      if (precis == 53u)
        precisionIncreased.push_back(this);
    }
    return p;
  }

  const P &get () {
    assert(!Parameter::handleSignException);
    int precis = uninitialized() ? 0 : precision();
    if (precis < Parameter::highPrecision) {
      p = calculate();

      if (Parameter::highPrecision > 53u && precis <= 53u)
        precisionIncreased.push_back(this);
    }
    else
      assert(precis == Parameter::highPrecision);
    return p;
  }

  void copyMod (P q) {
    assert(p.size() == q.size());
    for (int i = 0; i < p.size(); i++)
      p[i].copyMod(q[i]);
  }

  const P &getApprox (double accuracy=1e-17) {
    int precis = uninitialized() ? 0 : precision();
    if (precis < Parameter::highPrecision || 
        !checkAccuracy(accuracy, true)) {
      if (Parameter::handleSignException)
        safe_setp(accuracy);
      else {
	p = calculate();
	if (Parameter::highPrecision > 53u && precis <= 53u)
	  precisionIncreased.push_back(this);
        checkAccuracy(accuracy);
      }

      if (Parameter::highPrecision > 53u && precis <= 53u)
        precisionIncreased.push_back(this);
    }
    else
      assert(precis == Parameter::highPrecision);
    return p;
  }

  void set (const P &p) {
    assert(uninitialized());
    this->p = p;
    assert(!uninitialized());
    if (precision() > 53u && size() > 0)
      precisionIncreased.push_back(this);
  }

private:
  void safe_setp (double accuracy) {
    assert(Parameter::handleSignException == true);
    Parameter::handleSignException = false;
    int precis = uninitialized() ? 0 : precision();
    bool failed = false;
    while (true)
      try {
	if (uninitialized() || precision() < Parameter::highPrecision)
	  p = calculate();
        checkAccuracy(accuracy);
        if (failed) {
          decreasePrecision();  // ???
          decreaseAll();
          Parameter::highPrecision = 53u;
        }
        Parameter::handleSignException = true;
        return;
      } catch (SignException se) {
        failed = true;
      }
  }

  bool checkAccuracy (double accuracy, bool nothrow=false) {
    if (accuracy == 1.0) return true;
    for (int i = 0; i < p.size(); i++) {
      Parameter &x = p[i];
      if (nothrow && p[i].lb() < p[i].ub() && p[i].sign(false) == 0)
	return false;
      int s = p[i].sign();
      if (s == 0) {
        x = Parameter::constant(0);
        continue;
      }

      double l = x.lb();
      double u = x.ub();

      if ((l > 0 && u - l > accuracy * l) ||
          (u < 0 && u - l > accuracy * -u)) {
        double m = (l + u) / 2;
	double lm = (l + m) / 2;
	double mu = (m + u) / 2;
	if ((l < lm) + (lm < m) + (m < mu) + (mu < u) <= 3)
        // if (l == m || m == u)
          continue;
        if (!nothrow)
          (x - m).sign();
        return false;
      }
    }
    return true;
  }
};

template<class P>
class InputObject : public Object<P> {
public:
  InputObject (const P &q) { 
    P &p = *(P*)&q;
    for (int i = 0; i < p.size(); i++)
      if (p[i].lb() < p[i].ub()) {
        std::cerr << "Input object has non-trivial interval." << std::endl;
        assert(0);
      }
    Object<P>::set(p); 
  }

  virtual P calculate () {
    return Object<P>::getCurrentValue();
  }
};

template<class P>
class ObjPTR : public PTR< Object<P> > {
public:
  ObjPTR () {}
  ObjPTR (Object<P> *o) : PTR< Object<P> >(o) {}
  ObjPTR (const ObjPTR &p) : PTR< Object<P> >(p) {}
  const ObjPTR &operator= (const ObjPTR &p) {
    PTR< Object<P> >::operator=(p);
    return p;
  }
};

class Primitive {
  // Calculate the sign from the objects.
  virtual int sign () = 0;
  
public:
  operator int () {
    if (!Parameter::handleSignException)
      return sign();
    Parameter::handleSignException = false;
    bool failed = false;
    while (true)
      try {
        int s = sign();
        if (failed) {
          BaseObject::decreaseAll();
          Parameter::highPrecision = 53u;
        }
        Parameter::handleSignException = true;
        return s;
      } catch (SignException se) {
        failed = true;
      }
  }

  int unsafe () { return sign(); }
};

}  

// Macro to create one-argument primitive.
#define Primitive1(P, t1, v1)                                    \
  class P : public acp::Primitive {                              \
    t1 v1;                                                       \
  int sign ();                                                   \
  public:                                                        \
  P (t1 v1) : v1(v1) {}                                          \
  };

// Macro to create two-argument primitive.
#define Primitive2(P, t1, v1, t2, v2)                            \
  class P : public acp::Primitive {                              \
    t1 v1; t2 v2;                                                \
  int sign ();                                                   \
  public:                                                        \
  P (t1 v1, t2 v2) : v1(v1), v2(v2) {}                           \
  };

// Macro to create three-argument primitive.
#define Primitive3(P, t1, v1, t2, v2, t3, v3)                       \
  class P : public acp::Primitive {                                 \
    t1 v1; t2 v2; t3 v3;                                            \
  int sign ();                                                      \
  public:                                                           \
  P (t1 v1, t2 v2, t3 v3) : v1(v1), v2(v2), v3(v3) {}               \
  };

// Macro to create four-argument primitive.
#define Primitive4(P, t1, v1, t2, v2, t3, v3, t4, v4)                 \
  class P : public acp::Primitive {                                   \
    t1 v1; t2 v2; t3 v3; t4 v4;                                       \
    int sign ();                                                      \
  public:                                                             \
  P (t1 v1, t2 v2, t3 v3, t4 v4) : v1(v1), v2(v2), v3(v3), v4(v4) {}  \
  };

// Macro to create five-argument primitive.
#define Primitive5(P, t1, v1, t2, v2, t3, v3, t4, v4, t5, v5)           \
  class P : public acp::Primitive {                                     \
    t1 v1; t2 v2; t3 v3; t4 v4; t5 v5;                                  \
    int sign ();                                                        \
  public:                                                               \
  P (t1 v1, t2 v2, t3 v3, t4 v4, t5 v5)                                 \
    : v1(v1), v2(v2), v3(v3), v4(v4), v5(v5) {}                         \
  };

#endif
