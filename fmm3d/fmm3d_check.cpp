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
#include "fmm3d.hpp"
#include "common/vecmatop.hpp"

using std::cerr;
using std::endl;

/* ********************************************************************** */
int FMM3d::check(const DblNumVec& srcden, DblNumVec& trgval, int numchk)
{
  iA(_trgpos->n()>0); //have things to check
  int dim = this->dim();
  int sdof = this->sdof();
  int tdof = this->tdof();
  int snum = _srcpos->n();
  int tnum = _trgpos->n();   //int lclnum = this->lclnum();
  
  //1. select point to check
  vector<int> chkvec(numchk);
  for(int k=0; k<numchk; k++) {
	 chkvec[k] = int( floor(drand48()*tnum) );	 iA(chkvec[k]>=0 && chkvec[k]<tnum);
  }
  DblNumMat chkpos(dim, numchk);
  DblNumVec chkval(numchk*tdof);
  DblNumVec chkdns(numchk*tdof);
  for(int k=0; k<numchk; k++) {
	 for(int i=0; i<dim; i++)		chkpos(i, k) = (*_trgpos)(i, chkvec[k]);
	 for(int i=0; i<tdof;i++)		chkval(k*tdof+i) = trgval(chkvec[k]*tdof+i);
  }
  
  //2. compute and gather
  DblNumMat inter(tdof, snum*sdof);
  for(int i=0; i<numchk; i++) {
	 DblNumMat onechkpos(dim, 1, false, chkpos.clmdata(i));
	 DblNumVec onechkdns(tdof, false, chkdns.data()+i*tdof);
	 iC( _knl.kernel(*_srcpos, *_srcnor, onechkpos, inter) );
	 iC( dgemv(1.0, inter, srcden, 0.0, onechkdns) );
  }
  
  //3. distribute to individual
  for(int k=0; k<numchk; k++)
	 for(int i=0; i<tdof; i++)
		chkdns(k*tdof+i) -= chkval(k*tdof+i);
  
  double vn = 0;  double en = 0;
  for(int k=0; k<numchk; k++)
	 for(int i=0; i<tdof; i++) {
		vn += chkval(k*tdof+i) * chkval(k*tdof+i);
		en += chkdns(k*tdof+i) * chkdns(k*tdof+i);
	 }
  vn = sqrt(vn);
  en = sqrt(en);
  cerr<<"relative error: "<<en/vn<<endl;
  cerr<<"error norm    : "<<en<<endl;
  cerr<<"solution norm : "<<vn<<endl;
  
  return 0;
}



