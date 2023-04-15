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
#ifndef _KERNEL3D_HPP_
#define _KERNEL3D_HPP_

#include "common/nummat.hpp"

using std::vector;

//eqt: 1 2 3 4 5 6
//lyr: s d r p
//qnt: u p ...

enum {
  //laplace kernels
  KNL_LAP_S_U = 111,
  KNL_LAP_D_U = 121,
  KNL_LAP_I   = 191, //identity tensor
  //stokes kernels
  KNL_STK_F_U = 301, //kernel used by FMM3d algorithm for stokes equation
  KNL_STK_S_U = 311,
  KNL_STK_S_P = 312,
  KNL_STK_D_U = 321,
  KNL_STK_D_P = 322,
  KNL_STK_R_U = 331,
  KNL_STK_R_P = 332,
  KNL_STK_I   = 391, //identity tensor
  KNL_STK_E   = 392, //levi-civita tensor
  //navier kernels  //KNL_NAV_F_U = 501, //used for fmm
  KNL_NAV_S_U = 511, //single displacement
  KNL_NAV_D_U = 521, //double displacement
  KNL_NAV_R_U = 531,
  KNL_NAV_I   = 591, //identity tensor
  KNL_NAV_E   = 592, //levi-civita tensor
  //other kernels
  //error
  KNL_ERR = -1
};

//static class kernel
class Kernel3d
{
protected:
  int _kt;
  vector<double> _coefs;
  static double _mindif; //minimal difference
public:
  Kernel3d(): _kt(KNL_ERR) {;}
  Kernel3d(int kt, const vector<double>& coefs): _kt(kt), _coefs(coefs) {;}
  Kernel3d(const Kernel3d& c): _kt(c._kt), _coefs(c._coefs) {;}
  Kernel3d& operator=(const Kernel3d& c) {
	 _kt = c._kt; _coefs = c._coefs; return *this;
  }
  int& kt() { return _kt; }
  const int& kt() const { return _kt; }
  vector<double>& coefs() { return _coefs; }
  const vector<double>& coefs() const { return _coefs; }
  int dim() { return 3; }
  int sdof() const;
  int tdof() const;
  bool hom() const; //homogeneous or not
  void homsdeg(vector<int>&) const; //homogeneous degree, vector size == sdof
  //each kt handles coef in its own way, can be empty
  int kernel(const DblNumMat& srcpos, const DblNumMat& srcnor, const DblNumMat& trgpos, DblNumMat& inter);
};

inline bool operator==(const Kernel3d& a, const Kernel3d& b) {
  return (a.kt()==b.kt() && a.coefs()==b.coefs());
}

#endif
