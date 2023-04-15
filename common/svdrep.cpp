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
#include "blas.h"
#include "lapack.h"
#include "svdrep.hpp"

using std::min;
using std::max;
using std::abs;

int    SVDRep::_wssize = 4194304;
double SVDRep::_wsbuf[4194304];

/* ********************************************************************** */
int SVDRep::construct(double epsilon, const DblNumMat& K)
{
  int m = K.m();
  int n = K.n();
  int k = min(m, n);
  
  DblNumMat tU(m, k);
  DblNumVec tS(k);
  DblNumMat tVT(k, n);
  int wssize = _wssize;

  //SVD
  int INFO;
  char JOBU  = 'S';
  char JOBVT = 'S';
  iA( wssize >= max(3*min(m,n)+max(m,n), 5*min(m,n)));
  double* wsbuf = _wsbuf;
  DGESVD(&JOBU, &JOBVT, &m, &n, K.data(), &m, tS.data(), tU.data(), &m, tVT.data(), &k, wsbuf, &wssize, &INFO);
  iA(INFO==0);
  
  //cutoff
  double cutoff = epsilon*tS(0);
  int cnt=0;
  while(cnt< k)
    if(abs(tS(cnt)) >= cutoff)
      cnt++;
    else
      break;
  
  _matU.resize(m, cnt);
  _matS.resize(cnt);	
  _matVT.resize(cnt,n);
  
  for(int i=0; i<m; i++)
    for(int j=0; j<cnt; j++)
      _matU(i,j) = tU(i,j);
  for(int i=0; i<cnt; i++)
    _matS(i) = tS(i);
  for(int i=0; i<cnt; i++)
    for(int j=0; j<n; j++)
      _matVT(i,j) = tVT(i,j);
  
  return 0;
}

/* ********************************************************************** */
int SVDRep::dgemv(double alpha, const DblNumVec& X, double beta, DblNumVec& Y, double tol)
{
  iA(Y.m() == _matU.m());
  iA(X.m() == _matVT.n());
  iC( dgemv(alpha, X.data(), beta, Y.data(), tol) );
  
  return 0;
}

/* ********************************************************************** */
int SVDRep::dgemv(double alpha, double* X, double beta, double* Y, double tol)
{
  int K = 1; //prevent matrix of zero size
  while(K<_matS.m() && _matS(K)>=tol) 
	 K++;
  //buf = VT(1:K,:) * X
  double* buf = _wsbuf;
  {
	 char TRANS = 'N';
	 int M = _matVT.m();
	 int N = _matVT.n();
	 double ALPHA = 1.0;
	 double BETA = 0.0;
	 int INC = 1;
	 DGEMV(&TRANS, &K, &N, &ALPHA, _matVT.data(), &M, X, &INC, &BETA, buf, &INC); //first K rows
  }
  // buf = S(1:K) .* buf;
  for(int i=0; i<K; i++)
	 buf[i] = buf[i] * _matS(i);
  // y = U(:,1:K) * buf
  {
	 char TRANS = 'N';
	 int M = _matU.m(); //int N = _matU.n();
	 double ALPHA = alpha;
	 double BETA = beta;
	 int INC = 1;
	 DGEMV(&TRANS, &M, &K, &ALPHA, _matU.data(), &M, buf, &INC, &BETA, Y, &INC);	
  }
  
  return 0;
}
