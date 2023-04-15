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
#include "common/vecmatop.hpp"
#include "fmm3d.hpp"

using std::istringstream;

FMM3d::FMM3d(const string& p):
  KnlMat3d(p), _ctr(0,0,0), _rootlvl(0),
  _np(6), _let(NULL), _matmgnt(NULL)
{
}

FMM3d::~FMM3d()
{
  if(_let!=NULL)	 delete _let;
}

// ---------------------------------------------------------------------- 
int FMM3d::SE2TC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval)
{
  int TMAX = 1024;
  if(trgpos.n()<=TMAX) {
	 int M = trgpos.n() * _knl.tdof();
	 int N = srcpos.n() * _knl.sdof();
	 DblNumMat tmp(M,N);
	 iC( _knl.kernel(srcpos, srcnor, trgpos, tmp) );
	 iC( dgemv(1.0, tmp, srcden, 1.0, trgval) );
  } else {
	 int RUNS = (trgpos.n()-1) / TMAX + 1;
	 for(int r=0; r<RUNS; r++) {
		int stt = r*TMAX;
		int end = min((r+1)*TMAX, trgpos.n());
		int num = end-stt;
		int M = num * _knl.tdof();
		int N = srcpos.n() * _knl.sdof();
		DblNumMat tps(dim(), num, false, trgpos.data() + stt*dim() );
		DblNumVec tvl(num*_knl.tdof(), false, trgval.data() + stt*_knl.tdof());
		DblNumMat tmp(M,N);
		iC( _knl.kernel(srcpos, srcnor, tps, tmp) );
		iC( dgemv(1.0, tmp, srcden, 1.0, tvl) );
	 }
  }
  return (0);
}
// ---------------------------------------------------------------------- 
int FMM3d::SE2UC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, Point3 trgctr, double trgrad, const DblNumVec& srcden, DblNumVec& trgval)
{
  DblNumMat trgpos; iC( _matmgnt->locpos(UC, trgctr, trgrad, trgpos) );
  int M = trgpos.n() * _knl.tdof();
  int N = srcpos.n() * _knl.sdof();
  DblNumMat tmp(M,N);
  iC( _knl.kernel(srcpos, srcnor, trgpos, tmp) );
  iC( dgemv(1.0, tmp, srcden, 1.0, trgval) );
  return (0);
}
// ---------------------------------------------------------------------- 
int FMM3d::SE2DC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, Point3 trgctr, double trgrad, const DblNumVec& srcden, DblNumVec& trgval)
{
  DblNumMat trgpos; iC( _matmgnt->locpos(DC, trgctr, trgrad, trgpos) );
  int M = trgpos.n() * _knl.tdof();
  int N = srcpos.n() * _knl.sdof();
  DblNumMat tmp(M,N);
  iC( _knl.kernel(srcpos, srcnor, trgpos, tmp) );
  iC( dgemv(1.0, tmp, srcden, 1.0, trgval) );
  return 0;
}
// ---------------------------------------------------------------------- 
int FMM3d::DE2TC_dgemv(Point3 srcctr, double srcrad, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval)
{
  int TMAX = 1024;
  if(trgpos.n()<=TMAX) {
	 DblNumMat srcpos; iC( _matmgnt->locpos(DE, srcctr, srcrad, srcpos) );
	 int M = trgpos.n() * _knl_mm.tdof();
	 int N = srcpos.n() * _knl_mm.sdof();
	 DblNumMat tmp(M,N);
	 iC( _knl_mm.kernel(srcpos, srcpos, trgpos, tmp) );
	 iC( dgemv(1.0, tmp, srcden, 1.0, trgval) );
  } else {
	 DblNumMat srcpos; iC( _matmgnt->locpos(DE, srcctr, srcrad, srcpos) );
	 int RUNS = (trgpos.n()-1) / TMAX + 1;
	 for(int r=0; r<RUNS; r++) {
		int stt = r*TMAX;
		int end = min((r+1)*TMAX, trgpos.n());
		int num = end-stt;
		int M = num * _knl_mm.tdof();
		int N = srcpos.n() * _knl_mm.sdof();
		DblNumMat tps(dim(), num, false, trgpos.data() + stt*dim());
		DblNumVec tvl(num*_knl_mm.tdof(), false, trgval.data() + stt*_knl_mm.tdof());
		DblNumMat tmp(M, N);
		iC( _knl_mm.kernel(srcpos, srcpos, tps, tmp) );
		iC( dgemv(1.0, tmp, srcden, 1.0, tvl) );
	 }
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int FMM3d::UE2TC_dgemv(Point3 srcctr, double srcrad, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval)
{
  int TMAX = 1024;
  if(trgpos.n()<=TMAX) {
	 DblNumMat srcpos; iC( _matmgnt->locpos(UE, srcctr, srcrad, srcpos) );
	 int M = trgpos.n() * _knl_mm.tdof();
	 int N = srcpos.n() * _knl_mm.sdof();
	 DblNumMat tmp(M,N);
	 iC( _knl_mm.kernel(srcpos, srcpos, trgpos, tmp) );
	 iC( dgemv(1.0, tmp, srcden, 1.0, trgval) );
  } else {
	 DblNumMat srcpos; iC( _matmgnt->locpos(UE, srcctr, srcrad, srcpos) );
	 int RUNS = (trgpos.n()-1) / TMAX + 1;
	 for(int r=0; r<RUNS; r++) {
		int stt = r*TMAX;
		int end = min((r+1)*TMAX, trgpos.n());
		int num = end-stt;
		int M = num * _knl_mm.tdof();
		int N = srcpos.n() * _knl_mm.sdof();
		DblNumMat tps(dim(), num, false, trgpos.data() + stt*dim());
		DblNumVec tvl(num*_knl_mm.tdof(), false, trgval.data() + stt*_knl_mm.tdof());
		DblNumMat tmp(M,N);
		iC( _knl_mm.kernel(srcpos, srcpos, tps, tmp) );
		iC( dgemv(1.0, tmp, srcden, 1.0, tvl) );
	 }
  }
  return 0;
}

// ---------------------------------------------------------------------- 
DblNumMat FMM3d::sextpos(int gni) 
{
  Let3d::Node& node=_let->node(gni);
  int beg = node.sextbeg();
  int num = node.sextnum();
  return DblNumMat(dim(), num, false, _sextpos.data()+beg*dim());
}
DblNumMat FMM3d::sextnor(int gni)
{
  Let3d::Node& node=_let->node(gni);
  int beg = node.sextbeg();
  int num = node.sextnum();
  return DblNumMat(dim(), num, false, _sextnor.data()+beg*dim());
}
DblNumVec FMM3d::sextden(int gni)
{
  Let3d::Node& node=_let->node(gni);
  int beg = node.sextbeg();
  int num = node.sextnum();
  return DblNumVec(sdof()*num, false, _sextden.data()+beg*sdof());
}
DblNumVec FMM3d::supeden(int gni)
{
  Let3d::Node& node=_let->node(gni);
  int idx = node.sndeidx();
  return DblNumVec(datasize(UE), false, _supeden.data()+idx*datasize(UE));
}
DblNumVec FMM3d::supcval(int gni)
{
  Let3d::Node& node=_let->node(gni);
  int idx = node.sndeidx();
  return DblNumVec(datasize(UC), false, _supcval.data()+idx*datasize(UC));
}
// ---------------------------------------------------------------------- 
DblNumMat FMM3d::textpos(int gni)
{
  Let3d::Node& node=_let->node(gni);
  int beg = node.textbeg();
  int num = node.textnum();
  return DblNumMat(dim(), num, false, _textpos.data()+beg*dim());
}
DblNumVec FMM3d::textval(int gni)
{
  Let3d::Node& node=_let->node(gni);
  int beg = node.textbeg();
  int num = node.textnum();
  return DblNumVec(tdof()*num, false, _textval.data()+beg*tdof());
}
DblNumVec FMM3d::tdneden(int gni)
{
  Let3d::Node& node=_let->node(gni);
  int idx = node.tndeidx();
  return DblNumVec(datasize(DE), false, _tdneden.data()+idx*datasize(DE));
}
DblNumVec FMM3d::tdncval(int gni)
{
  Let3d::Node& node=_let->node(gni);
  int idx = node.tndeidx();
  return DblNumVec(datasize(DC), false, _tdncval.data()+idx*datasize(DC));
}


