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
#include "blas.h"
#include "lapack.h"
#include "svdrep.hpp"
#include "numvec.hpp"
#include "nummat.hpp"
#include "numtns.hpp"

#include "vecmatop.hpp"

// ---------------------------------------------------------------------- 
int dscal(double alpha, DblNumVec& X)
{
  dscal( X.m(), alpha, X.data() );
  return 0;
}
// ---------------------------------------------------------------------- 
int dscal(int n, double alpha, double* X)
{
  int incx = 1;
  DSCAL(&n, &alpha, X, &incx);
  return 0;
}
// ---------------------------------------------------------------------- 
int daxpy(double a, const DblNumVec& X, DblNumVec& Y)
{
  assert( X.m() == Y.m() );
  daxpy(X.m(), a, X.data(), Y.data());
  return 0;
}
// ---------------------------------------------------------------------- 
int daxpy(double a, const DblNumMat& X, DblNumMat& Y)
{
  assert( X.m() == Y.m() );  assert( X.n() == Y.n() );
  iC( daxpy(X.m()*X.n(), a, X.data(), Y.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int daxpy(int n, double a, double* X, double* Y)
{
  int incx = 1;  int incy = 1;
  DAXPY(&n, &a, X, &incx, Y, &incy);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  iC( dgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C)
{
  char transa = 'N';
  char transb = 'N';
  assert(m!=0 && n!=0 && k!=0);
  DGEMM(&transa, &transb, &m, &n, &k,
		  &alpha, A, &m, B, &k, &beta, C, &m);
  return 0;
}
// ---------------------------------------------------------------------- 
int dger(double alpha, const DblNumVec& X, const DblNumVec& Y, DblNumMat& A)
{
  assert(X.m() == A.m());
  assert(Y.m() == A.n());
  iC( dger(A.m(), A.n(), alpha, X.data(), Y.data(), A.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dger(int m, int n, double alpha, double* X, double* Y, double* A)
{
  assert(m!=0 && n!=0);
  int incx = 1;  int incy = 1;
  DGER(&m, &n, &alpha, X, &incx, Y, &incy, A, &m);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  iC( dgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dgemv(int m, int n, double alpha, double* A, double* X, double beta, double* Y)
{
  char trans = 'N';
  assert(m!=0 && n!=0);
  int incx = 1;
  int incy = 1;
  DGEMV(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
  return 0;
}
// ---------------------------------------------------------------------- 
int tran(const DblNumMat& M, DblNumMat& R)
{
  assert(R.m()==M.n() && R.n()==M.m());  //R.resize(M.n(), M.m());
  for(int i=0; i<M.m(); i++)
	 for(int j=0; j<M.n(); j++)
		R(j,i) = M(i,j);
  return 0;
}
// ----------------------------------------------------------------------
int pinv(const DblNumMat& M, double epsilon, DblNumMat& R)
{
  assert(M.m() == R.n());  assert(M.n() == R.m());
  SVDRep svd;
  iC( svd.construct(epsilon, M) );
  //invert Svd
  double cutoff = epsilon * svd.S()(0);
  for(int i=0; i<svd.S().m(); i++) {
    if( svd.S()(i) >= cutoff) {
      svd.S()(i) = 1.0/(svd.S()(i));
	 } else {
		assert(0);      //svd.S()(i) = 0.0;
	 }
  }
  DblNumMat UT(svd.U().n(),  svd.U().m());
  DblNumMat V( svd.VT().n(), svd.VT().m());
  iC( tran(svd.U(), UT) );
  iC( tran(svd.VT(), V) );
  for(int i=0; i<V.m(); i++)
    for(int j=0; j<V.n(); j++) {
      V(i,j) = V(i,j) * svd.S()(j);
	 }
  
  char transa = 'N';
  char transb = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int m = V.m();
  int n = UT.n();
  int k = V.n();
  DGEMM(&transa, &transb, &m, &n, &k, &alpha,
		  V.data(), &m, UT.data(), &k, 
		  &beta, R.data(), &m);  
  
  return 0;
}

// ---------------------------------------------------------------------- 
int inv(const DblNumMat& M, DblNumMat& R) //Gaussian Elimination
{
  //OR pinv(M, 0.0, R);
  assert(M.m()==M.n() && R.m()==R.n() && M.m()==R.m());
  memcpy(R.data(), M.data(), M.m()*M.n()*sizeof(double));
  int info;
  int m = M.m();
  int* ipiv = new int[m];
  DGETRF(&m, &m, R.data(), &m, ipiv, &info); assert(info==0);
  int lwork = m;
  double* work = new double[lwork];
  DGETRI(&m, R.data(), &m, ipiv, work, &lwork, &info);  assert(info==0);
  delete [] ipiv;
  delete [] work;
  return 0;
}
// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------- 
int spev1d(int evflag, int dmflag, int dof, double* data, int m, double e, int i, double u, double* res)
{
  int is[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++) is[k]=(i+k-1);
  }
  DblNumMat M(dof,m,false,data); //assert(M.n()==n && M.m()==res.m()); //double dof = M.m();  //int cnt = 0;
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 double scl = 1.0;
	 double us[4]; iC( spcoef(EVFLAG_VL, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++) 
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 double scl = double(m) / e;
	 double us[4]; iC( spcoef(EVFLAG_FD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 double scl = double(m*m)/(e*e);
	 double us[4]; iC( spcoef(EVFLAG_SD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int spev2d(int evflag, int dmflag, int dof, double* data, int* mn, double* ef, int* ij, double* uv, double* res)
{
  int m = mn[0];  int n = mn[1];
  double e = ef[0];  double f = ef[1];
  int i = ij[0];  int j = ij[1];
  double u = uv[0];  double v = uv[1];
  
  int is[4]; int js[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
	 for(int k=0; k<4; k++) js[k]=(j+k-1 + n) % n;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++)	is[k]=(i+k-1);
	 assert(j>=1 && j<=n-3);
	 for(int k=0; k<4; k++) js[k]=(j+k-1);
  }
  DblNumMat M(dof,m*n,false,data);
  double scl;
  double us[4], vs[4];
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 scl = 1.0;
	 iC( spcoef(EVFLAG_VL, u, us) );
	 iC( spcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 scl = double(m)/e;
	 iC( spcoef(EVFLAG_FD, u, us) );
	 iC( spcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b];
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n)/f;
	 iC( spcoef(EVFLAG_VL, u, us) );
	 iC( spcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 scl = double(m*m)/(e*e);
	 iC( spcoef(EVFLAG_SD, u, us) );
	 iC( spcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(m*n)/(e*f);
	 iC( spcoef(EVFLAG_FD, u, us) );
	 iC( spcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n*n)/(f*f);
	 iC( spcoef(EVFLAG_VL, u, us) );
	 iC( spcoef(EVFLAG_SD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int spcoef(int evflag, double u, double* us)
{
  double u1 = u;
  double u2 = u*u;
  double u3 = u*u*u;
  if(       evflag==EVFLAG_VL) {
	 us[0] = (  1 - 3*u1 + 3*u2 -   u3)/6.0;
	 us[1] = (  4        - 6*u2 + 3*u3)/6.0;
	 us[2] = (  1 + 3*u1 + 3*u2 - 3*u3)/6.0;
	 us[3] = (                  +   u3)/6.0;
  } else if(evflag==EVFLAG_FD) {
	 us[0] = (- 3 + 6*u1 - 3*u2)/6.0;
	 us[1] = (    -12*u1 + 9*u2)/6.0;
	 us[2] = (  3 + 6*u1 - 9*u2)/6.0;
	 us[3] = (             3*u2)/6.0;
  } else if(evflag==EVFLAG_SD) {
	 us[0] = (  6 - 6*u1 ) / 6.0;
	 us[1] = (-12 +18*u1 ) / 6.0;
	 us[2] = (  6 -18*u1 ) / 6.0;
	 us[3] = (      6*u1 ) / 6.0;	 //assert(0); //TODO;
  }  
  return 0;
}
// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------- 
int lagev1d(int evflag, int dmflag, int dof, double* data, int m, double e, int i, double u, double* res)
{
  int is[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++) is[k]=(i+k-1);
  }
  DblNumMat M(dof,m,false,data); //assert(M.n()==n && M.m()==res.m()); //double dof = M.m();  //int cnt = 0;
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 double scl = 1.0;
	 double us[4]; iC( lagcoef(EVFLAG_VL, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++) 
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof; //cnt ++;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 double scl = double(m) / e;
	 double us[4]; iC( lagcoef(EVFLAG_FD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof; //cnt ++;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 double scl = double(m*m)/(e*e);
	 double us[4]; iC( lagcoef(EVFLAG_SD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof; //cnt ++;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int lagev2d(int evflag, int dmflag, int dof, double* data, int* mn, double* ef, int* ij, double* uv, double* res)
{
  int m = mn[0];  int n = mn[1];
  double e = ef[0];  double f = ef[1];
  int i = ij[0];  int j = ij[1];
  double u = uv[0];  double v = uv[1];
  
  int is[4]; int js[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
	 for(int k=0; k<4; k++) js[k]=(j+k-1 + n) % n;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++)	is[k]=(i+k-1);
	 assert(j>=1 && j<=n-3);
	 for(int k=0; k<4; k++) js[k]=(j+k-1);
  }
  DblNumMat M(dof,m*n,false,data);
  double scl; 
  double us[4], vs[4];
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 scl = 1.0;
	 iC( lagcoef(EVFLAG_VL, u, us) );
	 iC( lagcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 scl = double(m)/e;
	 iC( lagcoef(EVFLAG_FD, u, us) );
	 iC( lagcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b];
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n)/f;
	 iC( lagcoef(EVFLAG_VL, u, us) );
	 iC( lagcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 scl = double(m*m)/(e*e);
	 iC( lagcoef(EVFLAG_SD, u, us) );
	 iC( lagcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(m*n)/(e*f);
	 iC( lagcoef(EVFLAG_FD, u, us) );
	 iC( lagcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n*n)/(f*f);
	 iC( lagcoef(EVFLAG_VL, u, us) );
	 iC( lagcoef(EVFLAG_SD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int lagev3d(int evflag, int dmflag, int dof, double* data, int* mno, double* efg, int* ijk, double* uvw, double* res)
{
  assert(evflag==EVFLAG_VL);
  
  int m = mno[0];  int n = mno[1];  int o = mno[2];
  double e = efg[0];  double f = efg[1];  double g = efg[2];
  int i = ijk[0];  int j = ijk[1];  int k = ijk[2];
  double u = uvw[0];  double v = uvw[1];  double w = uvw[2];
  
  int is[4]; int js[4];  int ks[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int h=0; h<4; h++) is[h]=(i+h-1 + m) % m;
	 for(int h=0; h<4; h++) js[h]=(j+h-1 + n) % n;
	 for(int h=0; h<4; h++) ks[h]=(k+h-1 + o) % o;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int h=0; h<4; h++)	is[h]=(i+h-1);
	 assert(j>=1 && j<=n-3);
	 for(int h=0; h<4; h++) js[h]=(j+h-1);
	 assert(k>=1 && k<=o-3);
	 for(int h=0; h<4; h++) ks[h]=(k+h-1);
  }
  DblNumMat M(dof,m*n*o,false,data);
  double scl; 
  double us[4], vs[4], ws[4];
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 scl = 1.0;
	 iC( lagcoef(EVFLAG_VL, u, us) );
	 iC( lagcoef(EVFLAG_VL, v, vs) );
	 iC( lagcoef(EVFLAG_VL, w, ws) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++)
		  for(int c=0; c<4; c++) {
			 double coef = us[a]*vs[b]*ws[c]; 
			 for(int d=0; d<dof; d++)
				res[d] += coef * M(d, is[a]+js[b]*m+ks[c]*m*n);
		  }
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int lagcoef(int evflag, double u, double* us)
{
  //double u1 = u;
  double u2 = u*u;
  double u3 = u*u*u;
  if(       evflag==EVFLAG_VL) {
	 us[0] = (u3-3*u2+2*u)/(-6.0);
	 us[1] = (u3-2*u2-u+2)/(2.0);
	 us[2] = (u3-u2-2*u)  /(-2.0);
	 us[3] = (u3-u)       /(6.0);
  } else if(evflag==EVFLAG_FD) {
	 us[0] = (3*u2-6*u+2)/(-6.0);
	 us[1] = (3*u2-4*u-1)/(2.0);
	 us[2] = (3*u2-2*u-2)/(-2.0);
	 us[3] = (3*u2-1)    /(6.0);
  } else if(evflag==EVFLAG_SD) {
	 assert(0); //TODO
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int fftrf1d(int, double*, int, int, double*)//int fftrf1d(int dof, double* data, int n, int ref, double* res)
{
  //TODO
  return 0;
}
// ---------------------------------------------------------------------- 
int fftrf2d(int dof, double* data, int* mm, int* rr, double* res)
{
  //DblNumMat M(dof,m   * n,  false,data);
  //DblNumMat R(dof,m*r * n*s,false,data);  //assert(M.m()==R.m() && M.n()==m*n && R.n()==m*r*n*s);  //int dof = M.m();
  //int m = mn[0];  int n = mn[1];  //int r = rs[0];  int s = rs[1];  //int m1 = m;  int m2 = n;  //int n1 = m*r;  int n2 = n*s;
  int m1 = mm[0];        int m2 = mm[1]; //old matrix size
  int n1 = mm[0]*rr[0];  int n2 = mm[1]*rr[1]; //new matrix size
  int mq1 = m1/2+1;  int mq2 = m2/2+1;
  int nq1 = n1/2+1;  //int nq2 = n2/2+1;
  //allocate space
  rfftwnd_plan regfftplan = rfftw2d_create_plan(m1,m2,FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
  rfftwnd_plan rfdfftplan = rfftw2d_create_plan(n1,n2,FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
  double* _sorg = data;//M.data();
  double* _sref = res; //R.data();
  double* _forg = new double[dof * 2*mq1 * m2];  assert(_forg!=NULL);
  double* _fref = new double[dof * 2*nq1 * n2];  assert(_fref!=NULL);
  memset( _forg, 0, dof*2*mq1*m2*sizeof(double) );
  memset( _fref, 0, dof*2*nq1*n2*sizeof(double) );
  //scale  DblNumVec sorgvec(dof*m1*m2, false, _sorg);  iC( dscal(1.0/double(m1*m2), sorgvec) );
  //do fft
  rfftwnd_real_to_complex(regfftplan, dof, _sorg, dof, 1, (fftw_complex*)_forg, dof, 1);
  //rearrange
  NumMat<fftw_complex> forgmat( dof*mq1, m2, false, (fftw_complex*)_forg );
  NumMat<fftw_complex> frefmat( dof*nq1, n2, false, (fftw_complex*)_fref );
  //    j direction
  for(int i=0; i<dof*mq1; i++) {
	 for(int j=0; j<mq2; j++) {
		frefmat(i,j) = forgmat(i,j);
	 }
	 for(int j=0; j<m2-mq2; j++) {
		frefmat(i,j+mq2+n2-m2) = forgmat(i,j+mq2);
	 }
	 if(m2 % 2 == 0) {
		frefmat(i,mq2-1).re /= 2.0;
		frefmat(i,mq2-1).im /= 2.0;
		frefmat(i,mq2+n2-m2-1).re = frefmat(i,mq2-1).re;
		frefmat(i,mq2+n2-m2-1).im = frefmat(i,mq2-1).im;
	 }
  }
  //    i direction
  if(m1 % 2 == 0) {
	 for(int i=(mq1-1)*dof; i<mq1*dof; i++) {
		for(int j=0; j<n2; j++) {
		  frefmat(i,j).re /= 2.0;
		  frefmat(i,j).im /= 2.0;
		}
	 }
  }
  //do ifft
  rfftwnd_complex_to_real(rfdfftplan, dof, (fftw_complex*)_fref, dof, 1, _sref, dof, 1);
  //scale
  DblNumVec srefvec(dof*n1*n2, false, _sref);  iC( dscal(1.0/double(m1*m2), srefvec) );
  //deallocate space
  delete [] _forg;
  delete [] _fref;
  rfftwnd_destroy_plan(regfftplan);
  rfftwnd_destroy_plan(rfdfftplan);
  
  return 0;
}


