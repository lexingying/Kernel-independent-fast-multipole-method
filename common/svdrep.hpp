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
/*! \file */
#ifndef _SVDREP_HPP_
#define _SVDREP_HPP_

#include "nummat.hpp"

using std::endl;

class SVDRep
{
public:
  SVDRep()   {}
  SVDRep(const SVDRep& C): _matU(C._matU), _matS(C._matS), _matVT(C._matVT)  {}
  ~SVDRep()  {}

  SVDRep& operator=(const SVDRep& c)  { _matU = c._matU; _matS = c._matS; _matVT = c._matVT; return *this; }
  //access
  DblNumMat& U() { return _matU; }
  DblNumVec& S() { return _matS; }
  DblNumMat& VT(){ return _matVT; }
  //ops
  int construct(double epsilon, const DblNumMat& M);
  int dgemv(double alpha, const DblNumVec& X, double beta, DblNumVec& Y, double tol=0.0); // y <- a Mx + b y
  int dgemv(double alpha, double* X, double beta, double* Y, double tol=0.0);
  
  int m() const { return _matU.m(); }
  int k() const { return _matS.m(); } //length
  int n() const { return _matVT.n(); }
  
protected:
  DblNumMat _matU;
  DblNumVec _matS;
  DblNumMat _matVT;
  static int _wssize;
  static double _wsbuf[];
  //static int _wisize;
  //static int _wibuf[];
};	

inline ostream& operator<<( ostream& os, SVDRep& svdrep)
{
  os<<svdrep.U().m()<<" "<<svdrep.S().m()<<" "<<svdrep.VT().n()<<endl;
  os<<svdrep.U()<<svdrep.S()<<svdrep.VT()<<endl;
  return os;
}

//int matvec(double a, const DblNumVec& X, double b, DblNumVec& Y); // y <- a Mx + b y
//int matvec(double a, double* X, double b, double* Y);

#endif
