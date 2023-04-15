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
#ifndef _FMM3D_HPP_
#define _FMM3D_HPP_

#include "knlmat3d.hpp"
#include "let3d.hpp"
#include "matmgnt3d.hpp"


//-------------------------------------------
class FMM3d: public KnlMat3d
{
public:
  typedef pair<int,int> intpair;
  enum {	 UE=0,	 UC=1,	 DE=2,	 DC=3,  };
  //------------------------------------
  class Node {
  protected:
	 int _Vinnum, _Vincnt;
	 DblNumVec _effval;
	 int _Votnum, _Votcnt; //for user
	 DblNumVec _effden;
  public:
	 Node() : _Vinnum(0), _Vincnt(0), _Votnum(0), _Votcnt(0) {;}
	 int& Vinnum() { return _Vinnum; }
	 int& Vincnt() { return _Vincnt; }
	 DblNumVec& effval() { return _effval; }
	 int& Votnum() { return _Votnum; }
	 int& Votcnt() { return _Votcnt; }
	 DblNumVec& effden() { return _effden; }
  };
  
protected:
  //PARAMS (REQ)
  Point3 _ctr;
  int    _rootlvl; //the level of the root box, radius of the box is 2^(-_rootlvl)
  //PARAMS (OPT)
  int _np;
  //COMPONENTS, local member and data  vector<int> _stphds, _evlhds;
  Let3d* _let;
  MatMgnt3d* _matmgnt;
  
  vector<Node> _ndevec;
  
  DblNumMat _sextpos, _sextnor;
  DblNumVec _sextden, _supeden, _supcval;
  DblNumMat _textpos;
  DblNumVec _textval, _tdneden, _tdncval;
  
  //IMPORTANT LEXING
  Kernel3d _knl_mm; //elq used in matmgnt
  int _mul_mm; //mul used in matmgnt  //FUNCTIONS: eq, lt, qt, sdof, tdof
public:
  FMM3d(const string& p);
  ~FMM3d();
  //MEMBER ACCESS
  Point3& ctr() { return _ctr; }
  int& rootlvl() { return _rootlvl; }
  int& np() { return _np; }  //int& usechk() { return _usechk; }  //int& numchk() { return _numchk; }
  
  //SETUP and USE
  int setup(map<string,string>& opts);
  int eval(const DblNumVec& srcden, DblNumVec& trgval);
  //   others
  int check(const DblNumVec& srcden, DblNumVec& trgval, int numchk);  //int report();
  //OTHER ACCESS
  Let3d* let() { return _let; }
  MatMgnt3d* matmgnt() { return _matmgnt; }
  vector<Node>& ndevec() { return _ndevec; }
  Node& node(int gni) { return _ndevec[gni]; }

protected:
  int datasize(int tp) { return _matmgnt->plndatasize(tp); }
  //multiplication
  int SE2TC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval);
  int SE2UC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, Point3 trgctr, double trgrad, const DblNumVec& srcden, DblNumVec& trgval);
  int SE2DC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, Point3 trgctr, double trgrad, const DblNumVec& srcden, DblNumVec& trgval);
  int DE2TC_dgemv(Point3 srcctr, double srcrad, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval);
  int UE2TC_dgemv(Point3 srcctr, double srcrad, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval);
  //contributor data
  DblNumMat sextpos(int gni);
  DblNumMat sextnor(int gni);
  DblNumVec sextden(int gni);
  DblNumVec supeden(int gni);
  DblNumVec supcval(int gni);
  //evaluator data
  DblNumMat textpos(int gni);
  DblNumVec textval(int gni);
  DblNumVec tdneden(int gni);
  DblNumVec tdncval(int gni);
  
  //auxilary functions
  int srcdata();
  int trgdata();
  //int setupGData();
  //int setupCData();
  //int setupUData();
  //int setupEData();
};


//IMPORTANT LEXING
//int reshape_UP(); //forward, elq form -> elq_mm form  //int reshape_DP(); //backward, elq_mm form -> elq form
//Equation feq() { return _feq; }  int flt() { return LYR_S; }  int      fqt() { return QNT_U; }  int      fdof() { return feq().sdof(flt()); }



#endif
