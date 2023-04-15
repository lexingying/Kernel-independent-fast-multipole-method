/* Kernel Independent Fast Multipole Method
   Copyright (C) 2004 Lexing Ying, New York University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */
#ifndef _MATMGNT3D_HPP_
#define _MATMGNT3D_HPP_

#include "common/nummat.hpp"
#include "common/numtns.hpp"
#include "common/offtns.hpp"
#include "common/vec3t.hpp"
#include "kernel3d.hpp"
#include "comobject.hpp"

#include "rfftw.h"

using std::map;
using std::pair;

//--------------------------------------
//unique identifier: equation
class MatMgnt3d
{
public:
  enum {	 UE=0,	 UC=1,	 DE=2,	 DC=3,  };
protected:
  //PARAMS(REQ) -- has to be set by parent class
  Kernel3d _knl; //the elq used by matmagnt (provide from fmm)
  int _np;
  bool _hom;
  vector<int> _degvec;
  //COMPONENTS
  map<int, DblNumMat> _uc2ue;
  map<int, OffTns<DblNumMat> > _ue2uc;
  map<int, DblNumMat> _dc2de;
  map<int, OffTns<DblNumMat> > _de2dc;
  map<int, OffTns<DblNumMat> > _ue2dc;
  DblNumMat _splpos[4]; //sample  positions
  DblNumMat _regpos; //regular positions
  rfftwnd_plan _forplan;
  rfftwnd_plan _invplan;
  
public:
  MatMgnt3d();
  ~MatMgnt3d();
  //MEMBER ACCESS
  Kernel3d& knl() { return _knl; }
  int& np() { return _np; }
  double alt(); //TODO: decide it based on np
  //...
  int sdof() { return _knl.sdof(); }  //int tdof() { return eq().tdof(qt()); }
  int tdof() { return _knl.tdof(); }
  int dim() { return 3; }
  //SETUP AND USE
  int setup();
  int report();
  int plndatasize(int tp); //the size of pln data
  int effdatasize(int tp); //the size of eff data
  
  int UC2UE_dgemv(int level,             const DblNumVec&, DblNumVec&);
  int UE2UC_dgemv(int level, Index3 ii, const DblNumVec&, DblNumVec&);
  int DC2DE_dgemv(int level,             const DblNumVec&, DblNumVec&);
  int DE2DC_dgemv(int level, Index3 ii, const DblNumVec&, DblNumVec&);
  int UE2DC_dgemv(int level, Index3 ii, const DblNumVec& effden, DblNumVec& effval);
  
  int plnden2effden(int level, const DblNumVec&, DblNumVec&); //pln->reg->eff
  int splden2regden(const DblNumVec&, DblNumVec&);
  int effval2plnval(int level, const DblNumVec&, DblNumVec&); //eff->reg->pln
  int regval2splval(const DblNumVec&, DblNumVec&);
  
  const DblNumMat& splpos(int tp) { return _splpos[tp]; }
  const DblNumMat& regpos()       { return _regpos; }
  int locpos(int, Point3, double, DblNumMat&);
  int splposcal(int n, double R, DblNumMat& ret); //calculate various grids
  int regposcal(int n, double R, DblNumMat& ret); //calculate the fft regular position grid
  int cptwvv(int, double, fftw_complex*, int, fftw_complex*, int, fftw_complex*, int);
  
public:
  static double _wsbuf[];
  static vector<MatMgnt3d> _mmvec;
public:
  static MatMgnt3d* getmmptr(Kernel3d, int);  //static void clearmmptrs();
};




/*	 int plnnum(ItlGrdType tp) { return _plnpos[tp].n(); }
	 int effnum(ItlGrdType tp) { return (2*_np+2)*(2*_np)*(2*_np); }
	 int regnum(ItlGrdType tp) { return _regpos[tp].n(); }  */
//int SE2TC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval);
//int SE2UC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, Point3 trgctr, double trgrad, const DblNumVec& srcden, DblNumVec& trgval);
//int SE2DC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, Point3 trgctr, double trgrad, const DblNumVec& srcden, DblNumVec& trgval);
//int DE2TC_dgemv(Point3 srcctr, double srcrad, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval);
//int UE2TC_dgemv(Point3 srcctr, double srcrad, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval);




#endif


