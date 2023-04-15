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

using std::cerr;
using std::endl;

using std::istringstream;

// ---------------------------------------------------------------------- 
int FMM3d::setup(map<string,string>& opts)
{
  //-----------------------------------------------------
  map<string,string>::iterator mi;
  mi = opts.find("-" + prefix() + "np"); iA(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>_np; }
  
  //-----------------------------------------------------
  iA(_srcpos!=NULL && _srcnor!=NULL && _trgpos!=NULL);
  iA((*_srcpos).m()==dim() && (*_trgpos).m()==dim());  //nothing to do
  //1. _let
  _let = new Let3d(prefix()+"let3d_");
  _let->srcpos()=_srcpos;  _let->trgpos()=_trgpos;  _let->ctr()=_ctr;  _let->rootlvl()=_rootlvl; //cerr<<" rootlvl "<<_rootlvl<<endl;
  iC( _let->setup(opts) );
  //2. decide _eq_mm and _mul_mm, and get matmgnt based on that
  switch(_knl.kt()) {
	 //laplace kernels
  case KNL_LAP_S_U: _knl_mm = Kernel3d(KNL_LAP_S_U, _knl.coefs()); break;
  case KNL_LAP_D_U: _knl_mm = Kernel3d(KNL_LAP_S_U, _knl.coefs()); break;
	 //stokes kernels
  case KNL_STK_S_U: _knl_mm = Kernel3d(KNL_STK_F_U, _knl.coefs()); break;
  case KNL_STK_S_P: _knl_mm = Kernel3d(KNL_LAP_S_U, vector<double>()); break;
  case KNL_STK_D_U: _knl_mm = Kernel3d(KNL_STK_F_U, _knl.coefs()); break;
  case KNL_STK_D_P: _knl_mm = Kernel3d(KNL_LAP_S_U, vector<double>()); break;
	 //navier kernels
  case KNL_NAV_S_U: _knl_mm = Kernel3d(KNL_NAV_S_U, _knl.coefs()); break;
  case KNL_NAV_D_U: _knl_mm = Kernel3d(KNL_NAV_S_U, _knl.coefs()); break;
  default: iA(0);
  }
  _mul_mm = 1; //for the time being
  
  _matmgnt  = MatMgnt3d::getmmptr(_knl_mm, _np);
  //3. self setup
  iC( srcdata() );
  iC( trgdata() );
  //-----------------------------------------------------
  return (0);
}
// ---------------------------------------------------------------------- 
int FMM3d::srcdata()
{
  //1. create vecs
  int sndecnt = _let->sndecnt();
  int sextcnt = _let->sextcnt();
  _sextpos.resize(dim(), sextcnt);
  _sextnor.resize(dim(), sextcnt);
  _sextden.resize(sextcnt * sdof());
  _supeden.resize(sndecnt * datasize(UE));
  _supcval.resize(sndecnt * datasize(UC));
  
  //2. gather the position using the pos scatter
  vector<int> ordvec;  iC( _let->upOrderCollect(ordvec) );
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if(_let->tag(gni) & LET_SRCNODE) { //contributor
		if(_let->terminal(gni)==true) { //terminal cbtr
		  DblNumMat sextpos(this->sextpos(gni));
		  DblNumMat sextnor(this->sextnor(gni));
		  vector<int>& curvis = _let->node(gni).sownvis();
		  for(int k=0; k<curvis.size(); k++) {
			 int poff = curvis[k];
			 for(int d=0; d<dim(); d++) {
				sextpos(d,k) = (*_srcpos)(d,poff);//pos
				sextnor(d,k) = (*_srcnor)(d,poff);//nor
			 }
		  }
		}
	 }
  }
  return (0);
}
// ---------------------------------------------------------------------- 
int FMM3d::trgdata()
{
  //1. create vecs
  int tndecnt = _let->tndecnt();
  int textcnt = _let->textcnt();
  _textpos.resize(dim(), textcnt);
  _textval.resize(textcnt * tdof());
  _tdneden.resize(tndecnt * datasize(DE));
  _tdncval.resize(tndecnt * datasize(DC));
  
  //2. gather data from _trgpos
  vector<int> ordvec; iC( _let->upOrderCollect(ordvec) );
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if(_let->tag(gni) & LET_TRGNODE) {
		if(_let->terminal(gni)==true) {
		  DblNumMat textpos(this->textpos(gni));
		  vector<int>& curvis = _let->node(gni).townvis();
		  for(int k=0; k<curvis.size(); k++) {
			 int poff = curvis[k];
			 for(int d=0; d<dim(); d++)
				textpos(d,k) = (*_trgpos)(d, poff);
		  }
		}
	 }
  }
  
  //3. allocate ENExt
  _ndevec.resize( _let->ndevec().size() );
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if(_let->tag(gni) & LET_TRGNODE) {
		//V
		Let3d::Node& gg = _let->node(gni);
		_ndevec[gni].Vinnum() = gg.Vnodes().size();
		for(vector<int>::iterator vi=gg.Vnodes().begin(); vi!=gg.Vnodes().end(); vi++) {
		  _ndevec[*vi].Votnum() ++;
		}
	 }
  }
  
  return (0);
}
