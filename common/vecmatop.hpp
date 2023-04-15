/*! \file */
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

#ifndef _VECMATOP_HPP_
#define _VECMATOP_HPP_

#include "nummat.hpp"

//--------------------------------------------------
//x = a x
int dscal(double alpha, DblNumVec& X);
int dscal(int n, double alpha, double* X);
//y = a x + y
int daxpy(double a, const DblNumVec& X, DblNumVec& Y);
int daxpy(double a, const DblNumMat& X, DblNumMat& Y);
int daxpy(int n, double a, double* X, double* Y);
// c = alpha*a*b + beta*c
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C);
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C);
// a = alpha* x*y' + a
int dger(double alpha, const DblNumVec& X, const DblNumVec& Y, DblNumMat& A);
int dger(int m, int n, double alpha, double* X, double* Y, double* A);
// y <= alpha A x + beta y
int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y);
int dgemv(int m, int n, double alpha, double* A, double* X, double beta, double* Y);
// R <= tran(M)
int tran(const DblNumMat& M, DblNumMat& R);
// R <= pinv(M, epsilon)
int pinv(const DblNumMat& M, double epsilon, DblNumMat& R);
// R <= inv(M);
int inv(const DblNumMat& M, DblNumMat& R);
//--------------------------------------------------
// interpolation, etc.
//evaluation flags
enum {  EVFLAG_VL = 1,  EVFLAG_FD = 2,  EVFLAG_SD = 4 };
//domain flag
enum {  DMFLAG_PERIOD = 0,  DMFLAG_CLOSED = 1 };

// cubic spline interpolation
//int spev1d( int evflag, int dmflag, const DblNumMat& M, int n,   int i,   double u,   DblNumMat& res);
//int spev2d( int evflag, int dmflag, const DblNumMat& M, int* mn, int* ij, double* uv, DblNumMat& res);
int spev1d( int evflag, int dmflag, int dof, double* M, int n,   double e,   int i,   double u,   double* res);
int spev2d( int evflag, int dmflag, int dof, double* M, int* mn, double* ef, int* ij, double* uv, double* res);
int spcoef( int evflag, double u, double* us);
// cubic lagrangian interpolation
//int lagev1d(int evflag, int pdflag, const DblNumMat& M, int n,   int i,   double u,   DblNumMat& res);
//int lagev2d(int evflag, int pdflag, const DblNumMat& M, int* mn, int* ij, double* uv, DblNumMat& res);
int lagev1d(int evflag, int pdflag, int dof, double* M, int n,   double e,   int i,   double u,   double* res);
int lagev2d(int evflag, int pdflag, int dof, double* M, int* mn, double* ef, int* ij, double* uv, double* res);
int lagev3d(int evflag, int pdflag, int dof, double* M, int* mno, double* efg, int* ijk, double* uvw, double* res);
int lagcoef(int evflag, double u, double* us);
// fft-based periodic refinement
//int fftrf1d(const DblNumMat& M, int  m,  int  ref, DblNumMat& R);
//int fftrf2d(const DblNumMat& M, int* mn, int* ref, DblNumMat& R);
int fftrf1d(int dof, double* M, int  m,  int  ref, double* res);
int fftrf2d(int dof, double* M, int* mn, int* ref, double* res);


#endif
