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

#include "poly.h"
using namespace std;
using namespace acp;

namespace acp {

bool Root2::polishing = false;

int PolySolver::Degree::sign () {
  PPoly<Parameter> p = poly->get();
  int d = p.deg();
  while (d >= 0 && zerop(p.a[d]))
    d--;
  return d;
}
/*
int Poly::deg ()
{
  PPoly<Parameter> p = getApprox(1.0);
  while (!p.zero() && p.lc().sign(false) == 0)
    p.a.pop_back();
  return p.deg();
}
*/

Roots PolySolver::getRoots ()
{
  double b = CauchyBound(poly).getApprox().ub();
  return getRoots(new InputKnot(-b), new InputKnot(b));
}

Roots PolySolver::getRoots (PTR<Knot> l, PTR<Knot> u)
{
  int d = deg();
  if (d < 1)
    return Roots();
  if (d == 1)
    return linearRoots(l, u);
  if (d == 2)
    return quadraticRoots(l, u);
  return descartesRoots(l, u);
}

Roots PolySolver::linearRoots (PTR<Knot> l, PTR<Knot> u)
{
  Roots res;
  RootPTR r = new LinearRoot(poly);
  if (Order(l, r) == 1 && Order(r, u) == 1)  
    res.push_back(r);
  return res;
}

Roots PolySolver::quadraticRoots (PTR<Knot> l, PTR<Knot> u)
{
  Roots res;
  if (QuadraticRoots(poly) == 1) {
    RootPTR r1 = new QuadraticRoot(poly, false), r2 = new QuadraticRoot(poly, true);
    if (Order(l, r1) == 1 && Order(r1, u) == 1)
      res.push_back(r1);
    if (Order(l, r2) == 1 && Order(r2, u) == 1)
      res.push_back(r2);
    if (res.size() == 2 && Order(res[0], res[1]) == -1)
      reverse(res.begin(), res.end());
  }
  return res;
}

Roots PolySolver::descartesRoots (PTR<Knot> l, PTR<Knot> u)
{
  Roots res;
  vector<DInt> st;
  st.push_back(DInt(l, u, 4));
  int k = 0;
  while (!st.empty()) {
    ++k;
    DInt i = *st.rbegin();
    st.pop_back();
    int v = Descartes(poly, i.l, i.u);
    if (v == 1)
      res.push_back(new PolyRoot(poly, i.l, i.u));
    else if (v > 1 // && !descartes1(st, i, v, true) && !descartes1(st, i, v, false)
	     && !descartes2(st, i, v)) {
      PTR<Knot> m = new SubKnot(i.l, i.u, 1, 2);
      unsigned int n = max(4u, (unsigned int) sqrt(i.n));
      st.push_back(DInt(i.l, m, n));
      st.push_back(DInt(m, i.u, n));
    }
  }
  //cerr << "descartes iterations = " << k << endl;
  reverse(res.begin(), res.end());
  return res;
}

bool PolySolver::descartes1 (vector<DInt> &st, const DInt &i, int v, bool lflag)
{
  PTR<Knot> m = new SubKnot(i.l, i.u, lflag ? 1 : i.n - 1, i.n),
    nl = lflag ? i.l : m, nu = lflag ? m : i.u;
  if (Descartes(poly, nl, nu) != v)
    return false;
  unsigned long int ln = i.n, nn = min(ln*ln, 4294967295ul);
  st.push_back(DInt(nl, nu, nn));
  return true;
}

bool PolySolver::descartes2 (vector<DInt> &st, const DInt &i, int v)
{
  int k = descartes2int(i, v);
  if (k == 0)
    return false;
  PTR<Knot> nl = new SubKnot(i.l, i.u, k - 2, 4*i.n),
    nu = new SubKnot(i.l, i.u, k + 2, 4*i.n);
  if (Descartes(poly, nl, nu) != v)
    return false;
  unsigned long int ln = i.n, nn = min(ln*ln, 4294967295ul);
  st.push_back(DInt(nl, nu, nn));
  return true;
}

int PolySolver::descartes2int (const DInt &i, int v)
{
  PPoly<Parameter> p = getApprox(1.0);
  Parameter il = i.l->getApprox(1.0), iu = i.u->getApprox(1.0),
    x = 0.5*(il + iu), dp = p.der(x);
  if (dp.sign(false) == 0)
    return 0;
  Parameter nx = x - v*p.value(x)/dp, w = (nx - il)/(iu - il);
  int k = min(max(4*i.n*w.mid(), 2.0), 4*i.n - 2.0);
  return k;
}

int QuadraticRoots::sign ()
{
  PPoly<Parameter> q = p->get();
  Parameter d = q[1]*q[1] - 4.0*q[0]*q[2];
  return d.sign();
}

int Order::sign ()
{
  return (b->get() - a->get()).sign();
}

int Descartes::sign ()
{
  PPoly<Parameter> q = p->get().moebius(l->get(), u->get());
  int v = 0u, d = q.deg();
  int s = q.a[d].sign();
  for (int i = d - 1u; i >= 0; --i)
    if (q.a[i].sign() == - s) {
      ++v;
      s = - s;
    }
  return v;
}

bool zerop (const Parameter &p)
{
  return p.sign() == 0;
}


PV2 Root2::calculate () {

  return polish(getCurrentValue(), f->getPoly(), g->getPoly());

}

PV2 Root2::polish(PV2 p, PPoly2 ff, PPoly2 gg) {

  vector<PV2> b;
  b.push_back(p);
  vector< ObjPTR<PPoly2> > func;
  func.push_back(new InputObject<PPoly2>(ff));
  func.push_back(new InputObject<PPoly2>(gg));

  //debugDraw(b, func);

  polishing = true;

  //printf("------------------------------------------------------------------------------------------\n\n");

  RootBoundary rb = {.I = p, .v = vector< vector<RootPTR> >()};

  //run Newton's until it fails then try to subdivide
  do {
    newton(rb, ff, gg);
  } while(subdivide(rb, ff, gg));

  polishing = false;

  //printf("I [%.16lf, %.16lf] [%.16lf, %.16lf]\n", rb.I.x.lb(), rb.I.x.ub(), rb.I.y.lb(), rb.I.y.ub());

  //printf("------------------------------------------------------------------------------------------\n\n");

  return rb.I;
}

//sign() == -1 if a < b
Primitive2(LessThan, ObjPTR<Parameter> , a, ObjPTR<Parameter> , b);
int LessThan::sign() {
  return (a->get() - b->get()).sign();
}

vector<int> Root2::intersections(RootBoundary &rb, PPoly2 f, PPoly2 g) {
 
  PV2 b = rb.I;
  vector< vector<RootPTR> > fv;
  vector< vector<RootPTR> > gv;

  for(int i = 0; i < 4; i++) fv.push_back(rb.v[i]);
  for(int i = 4; i < 8; i++) gv.push_back(rb.v[i]);

  vector<int> alt;

  for(int k = 0; k < 4; k++) {

    vector<RootPTR> frr = fv[k];
    vector<RootPTR> grr = gv[k];

    for(int i = 0, j = 0; i < frr.size() || j < grr.size(); ) {

      RootPTR fr = (i < frr.size() ? frr[i] : nullptr);
      RootPTR gr = (j < grr.size() ? grr[j] : nullptr);

      //k < 2 is bottom and right which are ascending roots, k >= 2 are top and left which are descending
      if(fr != nullptr && (gr == nullptr || (k < 2 ? LessThan(fr, gr) == -1 : LessThan(fr, gr) == 1))) {
        alt.push_back(0);
        i++;
      } else {
        alt.push_back(1);
        j++;
      }

    }

  }

  return alt;

}

int Root2::parity(vector<int> alt) {
 
  //find the first f intersection
  int k = 0;
  int j = 0;
  while(j < alt.size()) {
    if(alt[j] == 0) break;
    j++;
  }
  
  //if f doesn't intersect this cell, keep going
  if(j == alt.size()) {
    return 0;
  }

  int count = 0;
  k = (j+1)%alt.size();

  vector<int> par;

  //compute x1, x2, x3, x4, ...
  //count the number of g intersections inbetween pairs of f intersections
  while(k != j) {
    if(alt[k] == 1) {
      count++;
    } else {
      par.push_back(count);
      count = 0;
    }

    k = (k+1)%alt.size();
  }

  //don't forget the last pair
  par.push_back(count);

  //parity of intersections == parity of sum of odd indexed xj's
  int parity = 0;
  for(j = 1; j < par.size(); j += 2) {
    parity += par[j];
  }

  return parity % 2;

}

class SubXPoly : public Object<PPoly<Parameter>> {
  ObjPTR<PPoly2> f;
  ObjPTR<Parameter> x;
  PPoly<Parameter> calculate() { return f->get().subX(x->get()); }
public:
  SubXPoly (ObjPTR<PPoly2> f, ObjPTR<Parameter> x) : f(f), x(x) {}
};

class SubYPoly : public Object<PPoly<Parameter>> {
  ObjPTR<PPoly2> f;
  ObjPTR<Parameter> y;
  PPoly<Parameter> calculate() { return f->get().subY(y->get()); }
public:
  SubYPoly (ObjPTR<PPoly2> f, ObjPTR<Parameter> y) : f(f), y(y) {}
};

void Root2::setRoots(RootBoundary &rb, PPoly2 f, PPoly2 g) {

  PV2 b = rb.I;

  rb.v = vector< vector<RootPTR> >();

  ObjPTR<Parameter> xl = new InputObject<Parameter>(b.x.lbP());
  ObjPTR<Parameter> xu = new InputObject<Parameter>(b.x.ubP());
  ObjPTR<Parameter> yl = new InputObject<Parameter>(b.y.lbP());
  ObjPTR<Parameter> yu = new InputObject<Parameter>(b.y.ubP());
  
  ObjPTR<PPoly2> ff = new InputObject<PPoly2>(f);

  ObjPTR<PPoly<Parameter>> bot_f   = new SubYPoly(ff, yl);
  ObjPTR<PPoly<Parameter>> right_f = new SubXPoly(ff, xu);
  ObjPTR<PPoly<Parameter>> top_f   = new SubYPoly(ff, yu);
  ObjPTR<PPoly<Parameter>> left_f  = new SubXPoly(ff, xl);

  rb.v.push_back(PolySolver(bot_f).getRoots(xl, xu));
  rb.v.push_back(PolySolver(right_f).getRoots(yl, yu));
  rb.v.push_back(PolySolver(top_f).getRoots(xl, xu));
  rb.v.push_back(PolySolver(left_f).getRoots(yl, yu));

  ObjPTR<PPoly2> gg = new InputObject<PPoly2>(g);

  ObjPTR<PPoly<Parameter>> bot_g   = new SubYPoly(gg, yl);
  ObjPTR<PPoly<Parameter>> right_g = new SubXPoly(gg, xu);
  ObjPTR<PPoly<Parameter>> top_g   = new SubYPoly(gg, yu);
  ObjPTR<PPoly<Parameter>> left_g  = new SubXPoly(gg, xl);

  rb.v.push_back(PolySolver(bot_g).getRoots(xl, xu));
  rb.v.push_back(PolySolver(right_g).getRoots(yl, yu));
  rb.v.push_back(PolySolver(top_g).getRoots(xl, xu));
  rb.v.push_back(PolySolver(left_g).getRoots(yl, yu));


  //top and left need to be reversed for parity check
  reverse(rb.v[2].begin(), rb.v[2].end());
  reverse(rb.v[3].begin(), rb.v[3].end());
  reverse(rb.v[6].begin(), rb.v[6].end());
  reverse(rb.v[7].begin(), rb.v[7].end());

}

Root2::RootBoundary Root2::splitHoriz(RootBoundary &rb, Parameter c, PPoly2 f, PPoly2 g) {

  ObjPTR<Parameter> s = new InputObject<Parameter>(c);
  
  vector< vector<RootPTR> > bot(8);
  vector< vector<RootPTR> > top(8);

  //bot side for bot f and g roots are same
  bot[0] = rb.v[0];
  bot[4] = rb.v[4];

  //top side for top f and g roots are the same 
  top[2] = rb.v[2];
  top[6] = rb.v[6];

  
  //Separate roots around the split point, give left to left and right to right
  //for f [0, 4) and g [4, 8)
  for(int i = 0; i < 8; i += 4) {

    int li = 0;
    while(li < rb.v[i+1].size() && LessThan(rb.v[i+1][li], s) == -1) { li++; }
    //[0, li) are < s [li, size) are >= s


    int ui = 0;
    while(ui < rb.v[i+3].size() && LessThan(rb.v[i+3][ui], s) == 1) { ui++; }
    //[0, ui) are > s [ui, size) are <= s  

    bot[i+1] = vector<RootPTR>(rb.v[i+1].begin(), rb.v[i+1].begin()+li);
    top[i+1] = vector<RootPTR>(rb.v[i+1].begin()+li, rb.v[i+1].end());

    top[i+3] = vector<RootPTR>(rb.v[i+3].begin(), rb.v[i+3].begin()+ui);
    bot[i+3] = vector<RootPTR>(rb.v[i+3].begin()+ui, rb.v[i+3].end());

  }

  ObjPTR<PPoly2> ff = new InputObject<PPoly2>(f);
  ObjPTR<PPoly2> gg = new InputObject<PPoly2>(g);

  ObjPTR<Parameter> xl = new InputObject<Parameter>(rb.I.x.lbP());
  ObjPTR<Parameter> xu = new InputObject<Parameter>(rb.I.x.ubP());
 
  ObjPTR<PPoly<Parameter>> mid_line_f = new SubYPoly(ff, s);
  ObjPTR<PPoly<Parameter>> mid_line_g = new SubYPoly(gg, s);

  vector<RootPTR> mid_roots_f_asc = PolySolver(mid_line_f).getRoots(xl, xu);
  vector<RootPTR> mid_roots_g_asc = PolySolver(mid_line_g).getRoots(xl, xu);

  vector<RootPTR> mid_roots_f_desc(mid_roots_f_asc.begin(), mid_roots_f_asc.end());
  vector<RootPTR> mid_roots_g_desc(mid_roots_g_asc.begin(), mid_roots_g_asc.end());

  reverse(mid_roots_f_desc.begin(), mid_roots_f_desc.end());
  reverse(mid_roots_g_desc.begin(), mid_roots_g_desc.end());

  top[0] = mid_roots_f_asc;
  top[4] = mid_roots_g_asc;

  bot[2] = mid_roots_f_desc;
  bot[6] = mid_roots_g_desc;

  PV2 I = rb.I;

  Parameter yl = I.y.lbP();
  Parameter yu = I.y.ubP();

  rb = {.I = PV2(I.x, Parameter(c, yu)), .v = top};

  return {.I = PV2(I.x, Parameter(yl, c)), .v = bot};

}

Root2::RootBoundary Root2::splitVert(RootBoundary &rb, Parameter c, PPoly2 f, PPoly2 g) {

  ObjPTR<Parameter> s = new InputObject<Parameter>(c);
  
  vector< vector<RootPTR> > left(8);
  vector< vector<RootPTR> > right(8);

  //left side for left f and g roots are same
  left[3] = rb.v[3];
  left[7] = rb.v[7];

  //right side for right f and g roots are the same 
  right[1] = rb.v[1];
  right[5] = rb.v[5];

  
  //Separate roots around the split point, give left to left and right to right
  //for f [0, 4) and g [4, 8)
  for(int i = 0; i < 8; i += 4) {
    int li = 0;
    while(li < rb.v[i+0].size() && LessThan(rb.v[i+0][li], s) == -1) { li++; }
    //[0, li) are < s [li, size) are >= s


    int ui = 0;
    while(ui < rb.v[i+2].size() && LessThan(rb.v[i+2][ui], s) == 1) { ui++; }
    //[0, ui) are > s [ui, size) are <= s  

    left[i+0] = vector<RootPTR>(rb.v[i+0].begin(), rb.v[i+0].begin()+li);
    right[i+0] = vector<RootPTR>(rb.v[i+0].begin()+li, rb.v[i+0].end());

    right[i+2] = vector<RootPTR>(rb.v[i+2].begin(), rb.v[i+2].begin()+ui);
    left[i+2] = vector<RootPTR>(rb.v[i+2].begin()+ui, rb.v[i+2].end());
  }

  ObjPTR<PPoly2> ff = new InputObject<PPoly2>(f);
  ObjPTR<PPoly2> gg = new InputObject<PPoly2>(g);

  ObjPTR<Parameter> yl = new InputObject<Parameter>(rb.I.y.lbP());
  ObjPTR<Parameter> yu = new InputObject<Parameter>(rb.I.y.ubP());
 
  ObjPTR<PPoly<Parameter>> mid_line_f = new SubXPoly(ff, s);
  ObjPTR<PPoly<Parameter>> mid_line_g = new SubXPoly(gg, s);

  vector<RootPTR> mid_roots_f_asc = PolySolver(mid_line_f).getRoots(yl, yu);
  vector<RootPTR> mid_roots_g_asc = PolySolver(mid_line_g).getRoots(yl, yu);

  vector<RootPTR> mid_roots_f_desc(mid_roots_f_asc.begin(), mid_roots_f_asc.end());
  vector<RootPTR> mid_roots_g_desc(mid_roots_g_asc.begin(), mid_roots_g_asc.end());

  reverse(mid_roots_f_desc.begin(), mid_roots_f_desc.end());
  reverse(mid_roots_g_desc.begin(), mid_roots_g_desc.end());

  left[1] = mid_roots_f_asc;
  left[5] = mid_roots_g_asc;

  right[3] = mid_roots_f_desc;
  right[7] = mid_roots_g_desc;

  PV2 I = rb.I;

  Parameter xl = I.x.lbP();
  Parameter xu = I.x.ubP();

  rb = {.I = PV2(Parameter(c, xu), I.y), .v = right};

  return {.I = PV2(Parameter(xl, c), I.y), .v = left};;;;

}

//subdivide the cell once along major axis
bool Root2::subdivide(RootBoundary &rb, PPoly2 f, PPoly2 g) {

  PV2 I = rb.I;

  //printf("SUBDIVIDE\n\n");

  bool hse = Parameter::handleSignException;
  Parameter::handleSignException = false;
  unsigned int hp = Parameter::highPrecision;

  int even_count = 0;
  int odd_count  = 0;
  int odd_index = -1;

  vector<RootBoundary> b;

  try {

    if(rb.v.size() == 0) {
      setRoots(rb, f, g);
    }

    if(I.x.intervalWidth() > I.y.intervalWidth()) {

      Parameter xl = I.x.lbP();
      Parameter xu = I.x.ubP();
      Parameter xl75 = 0.75*xl + 0.25*xu;
      Parameter xl25 = 0.25*xl + 0.75*xu;
   
      b.push_back(splitVert(rb, xl75, f, g));
      b.push_back(splitVert(rb, xl25, f, g));
      b.push_back(rb);

    } else {
      Parameter yl = I.y.lbP();
      Parameter yu = I.y.ubP();
      Parameter yl75 = 0.75*yl + 0.25*yu;
      Parameter yl25 = 0.25*yl + 0.75*yu;

      b.push_back(splitHoriz(rb, yl75, f, g));
      b.push_back(splitHoriz(rb, yl25, f, g));
      b.push_back(rb);

    }

    for(int it = 0; it < b.size(); it++) {

      //vector<PV2> list;
      //list.push_back(bb);

      //vector< ObjPTR<PPoly2> > fs;
      //fs.push_back(new InputObject<PPoly2>(f));
      //fs.push_back(new InputObject<PPoly2>(g));

      //debugDraw(list, fs);

      if(parity(intersections(b[it], f, g))) {
        odd_count++;
        odd_index = it;
      } else {
        even_count++;
      }

    }

  } catch(SignException e) {
    Parameter::highPrecision = hp;
    Parameter::handleSignException = hse;
    return false;
  }
  
  Parameter::handleSignException = hse;

  assert(odd_count >= 1);
  //if(odd_count < 1)
  //  return false;

  if(odd_count >= 1) {
    //I = b[odd_index];
    rb = b[odd_index];
    return true;
  } else {
    return false;
  }

}

void Root2::newton(RootBoundary &rb, PPoly2 f, PPoly2 g) {

  PV2 I = rb.I;

  PPoly2 fx = f.derX();
  PPoly2 fy = f.derY();
  PPoly2 gx = g.derX();
  PPoly2 gy = g.derY();

  int count = 0;

  do {

    PV2 y(I.x.midP(), I.y.midP());
    
    PV2 foy(f.value(&y.x), g.value(&y.x));

    PV2 fr(fx.value(&I.x), fy.value(&I.x));
    PV2 gr(gx.value(&I.x), gy.value(&I.x));

    PV2 x;


    //let this be an exception in the end code
    if(!solve(fr, gr, foy, x)) {
      //printf("NEWTON determinant 0, need to subdivide\n\n");
      break;
    }

    PV2 newI = y - x;

    //printf("Newton Result [%20.16g, %20.16g] [%20.16g, %20.16g]\n\n", newI.x.lb(), newI.x.ub(), newI.y.lb(), newI.y.ub());


    //this should just be an assertion, not exception in the end code
    assert(newI.x.intersects(I.x) && newI.y.intersects(I.y));

    if(!newI.x.intersects(I.x) || !newI.y.intersects(I.y)) {
      //printf("NEWTON result disjoint\n\n");
      break;
    }

    vector<PV2> b;
    b.push_back(I);
    b.push_back(newI);
    vector< ObjPTR<PPoly2> > fs;
    fs.push_back(new InputObject<PPoly2>(f));
    fs.push_back(new InputObject<PPoly2>(g));
    
    //debugDraw(b, fs);
 
    PV2 inter(newI.x.intersect(I.x), newI.y.intersect(I.y));

    //this is a termination condition
    //if the new interval, (result of intersecting the result and old interval, is not a proper subset of the old interval
    //if(!inter.x.subset(I.x) && !inter.y.subset(I.y))

    //this is less efficient
    //if(inter.x.intervalWidth() >= I.x.intervalWidth() && inter.y.intervalWidth() >= I.y.intervalWidth()) {
    if(!inter.x.subset(I.x) && !inter.y.subset(I.y)) {
      //printf("NEWTON no progress\n\n");
      break;
    }
 /*
    int pI = parity(intersections(I, f, g));
    int pIn = parity(intersections(inter, f, g));

    if(pI == 0 || pIn == 0) {
      printf("here\n");
    }
*/
    I = inter;
  
    rb = {.I = I, .v = vector< vector<RootPTR> >()};

    //printf("      I       [%20.16g, %20.16g] [%20.16g, %20.16g]\n\n", I.x.lb(), I.x.ub(), I.y.lb(), I.y.ub());
    
    //printf("NEWTON STEP\n\n");

  } while(true);
  //} while(++count < 100);

  //return I;

}

bool Root2::solve(PV2 fr, PV2 gr, PV2 b, PV2 &sol) {

  //pivot on the smallest element!
  //compare two things: sing(false) for comparison

  //if det contains zero, need to subdivide
  Parameter d = fr.x*gr.y - fr.y*gr.x;
  if(!d.sign(false)) return false;

  if(!fr.x.sign(false) || !gr.x.sign(false)) return false;
 
  //pivot first row 
  if( (fr.x.abs() - gr.x.abs()).sign(false) >= 0) {  

    Parameter c = gr.x / fr.x;

    gr.x = gr.x - c*fr.x;
    gr.y = gr.y - c*fr.y;
    b.y  = b.y - c*b.x;

    if(!gr.y.sign(false)) return false;
    Parameter y = b.y / gr.y;

    if(!fr.x.sign(false)) return false;
    Parameter x = (b.x - fr.y*y) / fr.x;

    sol = PV2(x, y); 

  } 
  //pivot on the second row
  else {

    Parameter c = fr.x / gr.x;

    fr.x = fr.x - c*gr.x;
    fr.y = fr.y - c*gr.y;
    b.x = b.x - c*b.y;

    if(!fr.y.sign(false)) return false;
    Parameter y = b.x / fr.y;

    if(!gr.x.sign(false)) return false;
    Parameter x = (b.y - gr.y*y) / gr.x;

    sol = PV2(x, y);

  }

  return true;
}


}

// debug

void pp (Knot *p)
{
  cerr << p->getApprox().mid() << endl;
}

void pp (Object<PPoly<Parameter>> *p)
{
  p->getApprox(1.0).print();
}

