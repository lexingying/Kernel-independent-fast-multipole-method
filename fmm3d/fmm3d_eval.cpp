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

// ---------------------------------------------------------------------- 
int FMM3d::eval(const DblNumVec& srcden, DblNumVec& trgval)
{
  //-----------------------------------
  iA(srcden.m()==sdof()*(*_srcpos).n());  iA(trgval.m()==tdof()*(*_trgpos).n());
  
  //cerr<<"fmm src and trg numbers "<<pglbnum(_srcpos)<<" "<<pglbnum(_trgpos)<<endl;
  int sdof = this->sdof();
  int tdof = this->tdof();
  
  //1. zero out vecs
  setvalue(trgval, 0.0);
  
  setvalue(_sextden, 0.0);
  setvalue(_supeden, 0.0);
  setvalue(_supcval, 0.0);
  
  setvalue(_textval, 0.0);
  setvalue(_tdneden, 0.0);
  setvalue(_tdncval, 0.0);
  //clock_t ck0, ck1;
  //CLOCKING;
  //ck0 = clock();
  vector<int> ordvec; iC( _let->upOrderCollect(ordvec) ); //BOTTOM UP

  //2. for cbtr, load extden
  //ck0 = clock();
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if(_let->tag(gni) & LET_SRCNODE) {
		if(_let->terminal(gni)==true) {
		  DblNumVec sextden(this->sextden(gni));
		  vector<int>& curvis = _let->node(gni).sownvis();
		  for(int k=0; k<curvis.size(); k++) {
			 int poff = curvis[k];
			 for(int d=0; d<sdof; d++) {
				sextden(k*sdof+d) = srcden(poff*sdof+d);
			 }
		  }
		}
	 }
  }
  //ck1 = clock();  cout<<"load  "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;  
  
  //3. up computation
  //ck0 = clock();
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if(_let->tag(gni) & LET_SRCNODE) {		//GNTra gnt = _let->gni2gnt(gni);
		if(_let->depth(gni)>=2) {
		  DblNumVec supcvalgni(supcval(gni));
		  DblNumVec supedengni(supeden(gni));
		  if(_let->terminal(gni)==true) {
			 //S2M
			 iC( SE2UC_dgemv(sextpos(gni), sextnor(gni), _let->center(gni), _let->radius(gni), sextden(gni), supcvalgni) );
		  } else {
			 //M2M
			 for(int a=0; a<2; a++) for(int b=0; b<2; b++) for(int c=0; c<2; c++) {
				Index3 idx(a,b,c);
				int chi = _let->child(gni, idx);
				if(_let->tag(chi) & LET_SRCNODE) {
				  iC( _matmgnt->UE2UC_dgemv(_let->depth(chi)+_rootlvl, idx, supeden(chi), supcvalgni) );
				}
			 }
		  }
		  //M2M
		  iC( _matmgnt->UC2UE_dgemv(_let->depth(gni)+_rootlvl, supcvalgni, supedengni) );
		}
	 }
  }
  //ck1 = clock();  cout<<"up  "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;
  
  ordvec.clear();  iC( _let->dnOrderCollect(ordvec) );
  //U
  //ck0 = clock();
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if(_let->tag(gni) & LET_TRGNODE) { //evaluator
		if( _let->terminal(gni)==true ) { //terminal		  //DblNumVec etextval(this->etextval(gni));
		  Let3d::Node& gg = _let->node(gni);
		  DblNumVec textvalgni(textval(gni));
		  DblNumMat textposgni(textpos(gni));
		  for(vector<int>::iterator vi=gg.Unodes().begin(); vi!=gg.Unodes().end(); vi++) {
			 //S2T
			 iC( SE2TC_dgemv(sextpos(*vi), sextnor(*vi), textposgni, sextden(*vi), textvalgni) );
		  }
		}
	 }
  }
  //ck1 = clock();  cout<<"u list "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;
  //V
  //ck0 = clock();
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if( _let->tag(gni) & LET_TRGNODE) { //evaluator		//GNTra gnt = _let->gni2gnt(gni);
		Point3 gnictr(_let->center(gni));
		double D = 2.0 * _let->radius(gni);
		DblNumVec tdncval(this->tdncval(gni));
		Let3d::Node& gg = _let->node(gni);
		for(vector<int>::iterator vi=gg.Vnodes().begin(); vi!=gg.Vnodes().end(); vi++) {
		  Point3 victr(_let->center(*vi));
		  Index3 idx;		  for(int d=0; d<dim(); d++)			 idx(d) = int(round( (victr[d]-gnictr[d])/D ));
		  Node& srcptr = node(*vi);
		  Node& trgptr = node(gni);
		  if(srcptr.Votcnt()==0) {
			 srcptr.effden().resize( _matmgnt->effdatasize(UE) );			 setvalue(srcptr.effden(), 0.0);//1. resize effden
			 iC( _matmgnt->plnden2effden(_let->depth(gni)+_rootlvl, supeden(*vi),  srcptr.effden()) );			 //2. transform from upeden to effden
		  }
		  if(trgptr.Vincnt()==0) {
			 trgptr.effval().resize( _matmgnt->effdatasize(DC) );			 setvalue(trgptr.effval(), 0.0);			 //1. resize effval
		  }
		  //M2L		  //int tmplvl = _let->depth(gni)+_rootlvl;
		  iC( _matmgnt->UE2DC_dgemv(_let->depth(gni)+_rootlvl, idx, srcptr.effden(), trgptr.effval()) );
		  
		  srcptr.Votcnt()++;
		  trgptr.Vincnt()++;
		  if(srcptr.Votcnt()==srcptr.Votnum()) {
			 srcptr.effden().resize(0);			 //1. resize effden to 0
			 srcptr.Votcnt()=0;
		  }
		  if(trgptr.Vincnt()==trgptr.Vinnum()) {
			 iC( _matmgnt->effval2plnval(_let->depth(gni)+_rootlvl, trgptr.effval(), tdncval) );			 //1. transform from effval to dncval
			 trgptr.effval().resize(0); //2. resize effval to 0
			 trgptr.Vincnt()=0;
		  }
		}
	 }
  }
  //ck1 = clock();  cout<<"v list "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;
  //W
  //ck0 = clock();
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if( _let->tag(gni) & LET_TRGNODE ) {
		if( _let->terminal(gni)==true ) {
		  DblNumVec textval_gni(this->textval(gni));
		  Let3d::Node& gg = _let->node(gni);
		  for(vector<int>::iterator vi=gg.Wnodes().begin(); vi!=gg.Wnodes().end(); vi++) {
			 if(_let->terminal(*vi) && _let->node(*vi).sextnum()*sdof<_matmgnt->plndatasize(UE)) { //use ext instead
				//S2T
				iC( SE2TC_dgemv(sextpos(*vi), sextnor(*vi), textpos(gni), sextden(*vi), textval_gni) );
			 } else {
				//M2T
				int vni = *vi;				//GNTra vnt = _let->gni2gnt(*vi);
				iC( UE2TC_dgemv(_let->center(vni), _let->radius(vni), textpos(gni), supeden(*vi), textval_gni) );
			 }
		  }
		}
	 }
  }
  //ck1 = clock();  cout<<"w list "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;
  //X
  //ck0 = clock();
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if( _let->tag(gni) & LET_TRGNODE) {		//GNTra gnt = _let->gni2gnt(gni);
		Let3d::Node& gg = _let->node(gni);
		DblNumVec textval_gni(textval(gni));
		DblNumVec tdncval_gni(tdncval(gni));
		for(vector<int>::iterator vi=gg.Xnodes().begin(); vi!=gg.Xnodes().end(); vi++) {
		  if(_let->terminal(gni) && _let->node(gni).textnum()*tdof<_matmgnt->plndatasize(DC)) { //use ext instead
			 iC( SE2TC_dgemv(sextpos(*vi), sextnor(*vi), textpos(gni), sextden(*vi), textval_gni) );
		  } else {
			 //S2L
			 iC( SE2DC_dgemv(sextpos(*vi), sextnor(*vi), _let->center(gni), _let->radius(gni), sextden(*vi), tdncval_gni) );
		  }
		}
	 }
  }
  //ck1 = clock();  cout<<"x list "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;
  //7. combine
  //ck0 = clock();
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if( _let->tag(gni) & LET_TRGNODE ) { //evaluator		//GNTra gnt = _let->gni2gnt(gni);
		if(_let->depth(gni)>=3) {
		  int pargni = _let->parent(gni);		  //GNTra pargnt = _let->parent(gnt);
		  Index3 chdidx( _let->path(gni)-2 * _let->path(pargni) );
		  //L2L
		  DblNumVec tdncval_gni(tdncval(gni));
		  iC( _matmgnt->DE2DC_dgemv(_let->depth(pargni)+_rootlvl, chdidx, tdneden(pargni), tdncval_gni) );
		}
		if(_let->depth(gni)>=2) {
		  //L2L
		  DblNumVec tdneden_gni(tdneden(gni));
		  iC( _matmgnt->DC2DE_dgemv(_let->depth(gni)+_rootlvl, tdncval(gni), tdneden_gni) );
		}
		if(_let->terminal(gni)) {
		  //L2T
		  DblNumVec textval_gni(textval(gni));
		  iC( DE2TC_dgemv(_let->center(gni), _let->radius(gni), textpos(gni), tdneden(gni), textval_gni) );
		}
	 }
  }
  //ck1 = clock();  cout<<"combine  "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;
  
  //8. save tdtextval
  //ck0 = clock();
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if( _let->tag(gni) & LET_TRGNODE ) {
		if( _let->terminal(gni)==true ) {
		  DblNumVec textval(this->textval(gni));
		  vector<int>& curvis = _let->node(gni).townvis();
		  for(int k=0; k<curvis.size(); k++) {
			 int poff = curvis[k];
			 for(int d=0; d<tdof; d++) {
				trgval(poff*tdof+d) = textval(k*tdof+d);
			 }
		  }
		}
	 }
  }
  //ck1 = clock();  cout<<"save  "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;
  //CLOCKING
  //ck1 = clock();  cout<<"save  "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;
  
  //----------------
  //if(usechk()==1) {	 iC( check(srcden,trgval) );  }
  return (0);
}



