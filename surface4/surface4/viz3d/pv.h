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

#ifndef PV_H
#define PV_H

#include "acp.h"

namespace acp {

class PV2 {
public:
  int size () { return 2; }
  Parameter &operator[] (int i) { return i == 0 ? x : y; }

  PV2 () {}
  static PV2 input (double x, double y) {
    return PV2(Parameter::input(x), Parameter::input(y));
  }
  static PV2 constant (double x, double y) {
    return PV2(Parameter::constant(x), Parameter::constant(y));
  }
  PV2 (const Parameter &ix, const Parameter &iy) : x(ix), y(iy) {}
  PV2 (const PV2 &p) : x(p.x), y(p.y) {}
  PV2 & operator= (const PV2 &p) { x = p.x; y = p.y; return *this; }
  bool uninitialized () { return x.uninitialized() && y.uninitialized(); }
  Parameter getX () const { return x; }
  Parameter getY () const { return y; }
  Parameter dot (const PV2 &b) const { return x*b.x + y*b.y; }
  Parameter cross (const PV2 &b) const { return x*b.y - y*b.x; }
  PV2 operator+ (const PV2 &b) const { return PV2(x + b.x, y + b.y); }
  PV2 operator- (const PV2 &b) const { return PV2(x - b.x, y - b.y); }
  PV2 operator- () const { return PV2(- x, - y); }    
  PV2 operator* (double k) const { return PV2(x*k, y*k); }
  PV2 operator* (const Parameter &k) const { return PV2(x*k, y*k); }
  PV2 operator/ (double k) const { return PV2(x/k, y/k); }
  PV2 operator/ (const Parameter &k) const { return PV2(x/k, y/k); }
  PV2 unit () const { return *this / (*this).dot(*this).sqrt(); }

  Parameter length() const { return (*this).dot(*this).sqrt(); }

  // 2 times the signed area of the triangle this,b,c.
  Parameter area (const PV2 &b, const PV2 &c) const {
    return (*this - c).cross(b - c);
  }

  PV2 mid () const { return constant(x.mid(), y.mid()); }

  Parameter x, y;
};

inline PV2 operator* (double k, PV2 v) { return v * k; }
inline PV2 operator* (const Parameter &k, PV2 v) { return v * k; }

class PV3 {
 public:
  int size () { return 3; }

  PV3 () {}
  static PV3 input (double x, double y, double z) {
    return PV3(Parameter::input(x),
               Parameter::input(y),
               Parameter::input(z));
  }
  static PV3 constant (double x, double y, double z) {
    return PV3(Parameter::constant(x),
               Parameter::constant(y),
               Parameter::constant(z));
  }
  PV3 (const Parameter &ix, const Parameter &iy, const Parameter &iz)
    : x(ix), y(iy), z(iz) {}
  PV3 (const PV3 &p) : x(p.x), y(p.y), z(p.z) {}
  PV3 & operator= (const PV3 &p) 
    { x = p.x; y = p.y; z = p.z; return *this; }

  bool uninitialized () { return x.uninitialized() && y.uninitialized() && z.uninitialized(); }

  Parameter getX () const { return x; }
  Parameter getY () const { return y; }
  Parameter getZ () const { return z; }

  Parameter &operator[] (int i) { return i == 0 ? x : i == 1 ? y : z; }

  Parameter dot (const PV3 &b) const { 
    return x*b.x + y*b.y + z*b.z; 
  }

  PV3 operator+ (const PV3 &b) const {
    return PV3(x + b.x, y + b.y, z + b.z);
  }

  PV3 operator- () const {
    return PV3(- x, - y, - z);
  }

  PV3 operator- (const PV3 &b) const {
    return PV3(x - b.x, y - b.y, z - b.z);
  }

  PV3 operator* (double k) const {
    return PV3(x*k, y*k, z*k);
  }

  PV3 operator* (const Parameter &k) const {
    return PV3(x*k, y*k, z*k);
  }

  PV3 operator/ (double k) const {
    return PV3(x/k, y/k, z/k);
  }

  PV3 operator/ (const Parameter &k) const {
    return PV3(x/k, y/k, z/k);
  }

  PV3 cross (const PV3 &b) const {
    return PV3(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);
  }

  Parameter tripleProduct (const PV3 &b, const PV3 &c) const {
    return dot(b.cross(c));
  }

  PV3 unit () const { return *this / (*this).dot(*this).sqrt(); }

  Parameter length() const { return (*this).dot(*this).sqrt(); }

  PV3 mid () const { return constant(x.mid(), y.mid(), z.mid()); }

  Parameter x, y, z;
};

inline PV3 operator* (double k, PV3 v) { return v * k; }
inline PV3 operator* (const Parameter &k, PV3 v) { return v * k; }

}

#endif
