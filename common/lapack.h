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

#ifndef _LAPACK_H_
#define _LAPACK_H_

#include "commoninc.hpp"

//blas and lapack workspace for the code

#define DGESVD dgesvd_
#define DGESDD dgesdd_
#define DGETRF dgetrf_
#define DGETRI dgetri_

//EXTERN_C_BEGIN
extern "C"
{
  extern void DGESVD(char *JOBU, char *JOBVT, int *M, int *N, double *A, int *LDA, 
							double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *INFO);
  extern void DGESDD(char *jobz, int* m, int* n, double* a, int* lda,
							double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* iwork, int* info);
  extern void DGETRF(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
  extern void DGETRI(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);
}
//EXTERN_C_END

#endif
