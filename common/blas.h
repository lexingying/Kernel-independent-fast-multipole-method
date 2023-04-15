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

#ifndef _BLAS_H_
#define _BLAS_H_

#include "commoninc.hpp"

#define DGEMM dgemm_
#define DGEMV dgemv_
#define DAXPY daxpy_
#define DGER  dger_
#define DSCAL dscal_

extern "C"
{
  void DAXPY(int* N, double* ALPHA, double* X, int* INCX, double* Y, int* INCY); 
  void DGEMM(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA, double* A,
				 int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);
  void DGEMV(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX,
				 double* BETA, double* Y, int* INCY);
  void DGER (int* M, int * N, double* ALPHA, double* X, int* INCX, double* Y, int* INCY,
				 double* A, int* LDA);
  void DSCAL(int* N, double* ALPHA, double* X, int* INCX);
}

#endif

