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
#include "dense3d.hpp"
#include "common/vecmatop.hpp"

Dense3d::Dense3d(const string& p): KnlMat3d(p)
{
}

Dense3d::~Dense3d()
{
}

// ----------------------------------------------------------------------
int Dense3d::setup(map<string,string>&)
{
  iA(_srcpos!=NULL && _srcnor!=NULL && _trgpos!=NULL);
  iA((*_srcpos).m()==dim() && (*_trgpos).m()==dim());  //nothing to do
  return 0;
}

// ---------------------------------------------------------------------- 
int Dense3d::eval(const DblNumVec& srcden, DblNumVec& trgval) 
{
  //-----------------------------------
  iA(srcden.m()==sdof()*(*_srcpos).n());  iA(trgval.m()==tdof()*(*_trgpos).n());
  
  int dim  = this->dim();
  int sdof = this->sdof();
  int tdof = this->tdof();
  int srcnum = (*_srcpos).n();
  int trgnum = (*_trgpos).n();
  
  DblNumMat inter(tdof, srcnum*sdof);
  for(int i=0; i<trgnum; i++) {
	 DblNumMat oneposmat(dim, 1, false, (*_trgpos).clmdata(i));
	 DblNumVec onevalvec(tdof, false, trgval.data()+tdof*i);
	 iC( _knl.kernel((*_srcpos), (*_srcnor), oneposmat, inter) );
	 iC( dgemv(1.0, inter, srcden, 0.0, onevalvec) );
  }
  
  return 0;
}
