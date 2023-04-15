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
#include "let3d.hpp"

using std::min;
using std::max;
using std::set;
using std::queue;
using std::ofstream;
using std::cerr;
using std::endl;
using std::istringstream;

//-----------------------------------------------------------------------
Let3d::Let3d(const string& p):
  ComObject(p), _srcpos(NULL), _trgpos(NULL), _ctr(0.0), _rootlvl(0),
  _ptsmax(150), _maxlevel(10)
{
}

Let3d::~Let3d()
{
}

// ---------------------------------------------------------------------- 
int Let3d::setup(map<string,string>& opts)
{
  //------------
  map<string,string>::iterator mi;
  mi = opts.find("-" + prefix() + "ptsmax"); iA(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>_ptsmax; }
  mi = opts.find("-" + prefix() + "maxlevel"); iA(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>_maxlevel; }
  //------------
  iC( srcdata() );
  iC( trgdata() );  //iC( print("tree") );
  return 0;
}

// ---------------------------------------------------------------------- 
int Let3d::srcdata()
{
  //-----------------------------------------
  //gdata
  DblNumMat& pos = *(_srcpos);  iA( pos.m()==dim() );
  
  vector<Node>& ndevec = this->ndevec(); ndevec.clear();
  vector< vector<int> > visvec;  //vector<int> lsmvec; //local sum
  vector<int> lsmvec; //glb 
  
  ndevec.push_back( Node(-1,-1, Index3(0,0,0), 0) );
  visvec.push_back( vector<int>() );
  vector<int>& curvis = visvec[0];
  Point3 bbmin(ctr()-Point3(rad()));
  Point3 bbmax(ctr()+Point3(rad()));
  for(int k=0; k<pos.n(); k++) {
	 Point3 tmp(pos.clmdata(k));	 //posarr+(k-plclstart)*dim());
	 iA(tmp>=bbmin && tmp<=bbmax);//LEXING: IMPORANT
	 curvis.push_back(k);
  }
  lsmvec.push_back( curvis.size() );
  
  int level = 0;
  int arrbeg = 0;
  int arrend = 1;
  while(arrbeg<arrend) {
	 //1.
	 int arrcnt = arrend;
	 for(int k=arrbeg; k<arrend; k++) {
		//---
		if( lsmvec[k]>ptsmax() && level<maxlevel()-1 ) {
		  ndevec[k].chd() = arrcnt;
		  arrcnt = arrcnt + pow2(dim());
		  //children's ess		  //for(int ord=0; ord<8; ord++) {
		  for(int a=0; a<2; a++) for(int b=0; b<2; b++) for(int c=0; c<2; c++) {
			 ndevec.push_back( Node(k,-1, 2*ndevec[k].path()+Index3(a,b,c), ndevec[k].depth()+1) ); //par, chd
			 visvec.push_back( vector<int>() );
			 lsmvec.push_back( 0 );
		  }
		  //children's vis
		  Point3 curctr( center(k) ); //get center of current node
		  for(vector<int>::iterator pi=visvec[k].begin(); pi!=visvec[k].end(); pi++) {
			 Point3 tmp(pos.clmdata(*pi));
			 Index3 idx;			 for(int j=0; j<dim(); j++)				idx(j) = (tmp(j) >= curctr(j));
			 int chdgni = child(k, idx);
			 visvec[chdgni].push_back(*pi);
		  }
		  visvec[k].clear(); //VERY IMPORTANT
		  //children's lsm		  //for(IdxIter ii=ran.begin(); ii!=ran.end(); ++ii) {		  //for(int ord=0; ord<8; ord++) {
		  for(int a=0; a<2; a++) for(int b=0; b<2; b++) for(int c=0; c<2; c++) {
			 int chdgni = child( k, Index3(a,b,c) );
			 lsmvec[chdgni] = visvec[chdgni].size();
		  }
		}
	 }
	 level++;
	 arrbeg = arrend;
	 arrend = arrcnt;
  }
  _level = level; //SET LEVEL

  vector<int> ordvec;  iC( dnOrderCollect(ordvec) );
  //set other parts of essvec
  int cnt = 0;
  int sum = 0;
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if(lsmvec[gni]>0) {
		ndevec[gni].tag() = ndevec[gni].tag() | LET_SRCNODE;
		ndevec[gni].sndeidx() = cnt;
		cnt++;
		if(ndevec[gni].chd()==-1) {
		  ndevec[gni].sextbeg() = sum;
		  ndevec[gni].sextnum() = lsmvec[gni];
		  sum += lsmvec[gni];
		  ndevec[gni].sownvis() = visvec[gni];
		}
	 }
  }
  _sndecnt = cnt;  _sextcnt = sum; //SET S cnts
  return 0;
}

// ---------------------------------------------------------------------- 
int Let3d::trgdata()
{
  //-----------------------------------------
  //edata
  DblNumMat& pos = *(_trgpos);  iA( pos.m()==dim() );
  
  vector<Node>& ndevec = this->ndevec();
  vector< vector<int> > visvec; visvec.resize(ndevec.size() );
  vector<int> lsmvec;           lsmvec.resize(ndevec.size(), 0);
  
  vector<int>& curvis = visvec[0];
  Point3 bbmin(ctr()-Point3(rad()));
  Point3 bbmax(ctr()+Point3(rad()));  //cerr<<" bbx "<<bbmin<<" "<<bbmax<<endl;
  for(int k=0; k<pos.n(); k++) {
	 Point3 tmp(pos.clmdata(k));
	 iA(tmp>=bbmin && tmp<=bbmax);	 //LEXING: IMPORTANT
	 curvis.push_back(k);
  }
  lsmvec[0] = curvis.size();
  
  vector<int> ordvec;  iC( dnOrderCollect(ordvec) );
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 Node& curnd = ndevec[gni];
	 vector<int>& curvis = visvec[gni];	 //int&         curlsm = lsmvec[gni];
	 if(curnd.chd()!=-1) { //not terminal
		//children's vis
		Point3 curctr( center(gni) );
		for(vector<int>::iterator pi=curvis.begin(); pi!=curvis.end(); pi++) {
		  Point3 tmp(pos.clmdata(*pi));
		  Index3 idx;		  for(int j=0; j<dim(); j++)			 idx(j) = (tmp(j)>=curctr(j));
		  int chdgni = child(gni, idx);
		  vector<int>& chdvis = visvec[chdgni];
		  chdvis.push_back(*pi);
		}
		curvis.clear(); //VERY IMPORTANT
		//children's lsm		//		for(int ord=0; ord<8; ord++) {
		for(int a=0; a<2; a++) for(int b=0; b<2; b++) for(int c=0; c<2; c++) {
		  int chdgni = child(gni, Index3(a,b,c));
		  lsmvec[chdgni] = visvec[chdgni].size();
		}
	 }
  }
  //set EVTR
  int cnt = 0;
  int sum = 0;
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if(lsmvec[gni]>0) { //evtr node
		ndevec[gni].tag() = ndevec[gni].tag() | LET_TRGNODE;
		ndevec[gni].tndeidx() = cnt;
		cnt ++;
		if(ndevec[gni].chd()==-1) { //terminal
		  ndevec[gni].textbeg() = sum;
		  ndevec[gni].textnum() = lsmvec[gni];
		  sum += lsmvec[gni];
		  ndevec[gni].townvis() = visvec[gni];
		}
	 }
  }
  _tndecnt = cnt;  _textcnt = sum;
  
  //set USER
  for(int i=0; i<ordvec.size(); i++) {
	 int gni = ordvec[i];
	 if(ndevec[gni].tag() & LET_TRGNODE) { //a evtr		//if(gptrvec[gni]==NULL) {  //calculate the U,V,W,X lists
		iC( calgnext(gni) );
	 }
  }
  //cerr<<usndecnt<<" "<<usextcnt<<endl;
  return 0;
}

// ---------------------------------------------------------------------- 
int Let3d::calgnext(int gni)
{
  vector<Node>&  ndevec = this->ndevec();
  
  set<int> Uset, Vset, Wset, Xset;
  int curgni = gni;
  if(root(curgni)==false) {
	 //GNTra curgnt = gni2gnt(curgni);
	 //GNTra pargnt = parent( curgnt);	 //int pargni = parent(curgni); //int   pargni = gnt2gni(pargnt);
	 int pargni = parent(curgni);
	 
	 Index3 minidx(0);
	 Index3 maxidx(pow2(depth(curgni)));
	 //IdxRan ran(-2,4);	 for(IdxIter ii=ran.begin(); ii!=ran.end(); ++ii) {
	 for(int i=-2; i<4; i++)		for(int j=-2; j<4; j++)		  for(int k=-2; k<4; k++) {
		Index3 trypath( 2*path(pargni) + Index3(i,j,k) );
		if(trypath >= minidx &&
			trypath <  maxidx &&
			trypath != path(curgni)) {		  //GNTra wntgnt(curgnt.depth(), path);
		  int resgni = findgnt(depth(curgni), trypath);		  //int resgni = resgnt.gni();
		  bool adj = adjacent(resgni, curgni);
		  if( depth(resgni)<depth(curgni) ) 
			 if(adj)
				if(terminal(curgni))
				  Uset.insert(resgni);
				else
				  ;
			 else
				Xset.insert(resgni);
		  if( depth(resgni)==depth(curgni) )
			 if(!adj) {
				Index3 bb(path(resgni)-path(curgni));
				assert( bb.linfty()<=3 );
				Vset.insert(resgni);
			 }
			 else
				if(terminal(curgni)) {
				  queue<int> rest;
				  rest.push(resgni);
				  while(rest.empty()==false) {
					 int fntgni = rest.front(); rest.pop();					 //int fntgni = fntgnt.gni();
					 if(adjacent(fntgni, curgni)==false)
						Wset.insert( fntgni );
					 else
						if(terminal(fntgni))
						  Uset.insert(fntgni);
						else { //IdxRan fran(0,2); for(IdxIter fii=fran.begin(); fii!=fran.end(); ++fii) //for(int ord=0; ord<8; ord++)
						  for(int a=0; a<2; a++) for(int b=0; b<2; b++) for(int c=0; c<2; c++)
							 rest.push( child(fntgni, Index3(a,b,c)) );
						}
				  }
				}
		}
	 }
  }
  if(terminal(curgni))
	 Uset.insert(curgni);
  
  for(set<int>::iterator si=Uset.begin(); si!=Uset.end(); si++)
	 if(ndevec[*si].tag() & LET_SRCNODE)		ndevec[gni].Unodes().push_back(*si);
  for(set<int>::iterator si=Vset.begin(); si!=Vset.end(); si++)
	 if(ndevec[*si].tag() & LET_SRCNODE)		ndevec[gni].Vnodes().push_back(*si);
  for(set<int>::iterator si=Wset.begin(); si!=Wset.end(); si++)
	 if(ndevec[*si].tag() & LET_SRCNODE)		ndevec[gni].Wnodes().push_back(*si);
  for(set<int>::iterator si=Xset.begin(); si!=Xset.end(); si++)
	 if(ndevec[*si].tag() & LET_SRCNODE)		ndevec[gni].Xnodes().push_back(*si);
  
  return (0);
}
// ---------------------------------------------------------------------- 
int Let3d::dnOrderCollect(vector<int>& ordvec)
{
  ordvec.clear();
  for(int i=0; i<ndevec().size(); i++)
	 ordvec.push_back(i);
  iA(ordvec.size()==ndevec().size());
  return 0;
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "Let3d::upOrderCollect"
int Let3d::upOrderCollect(vector<int>& ordvec)
{
  ordvec.clear();
  for(int i=ndevec().size()-1; i>=0; i--)
	 ordvec.push_back(i);
  iA(ordvec.size()==ndevec().size());
  return 0;
}
// ---------------------------------------------------------------------- 
int Let3d::child(int gni, const Index3& idx)
{
  assert(idx>=Index3(0) && idx<Index3(2));
  return node(gni).chd() + (idx(0)*4+idx(1)*2+idx(2));
}
Point3 Let3d::center(int gni) //center of a node
{
  Point3 ll( ctr() - Point3(rad()) );
  int tmp = pow2(depth(gni));
  Index3 pass(path(gni));
  Point3 res;
  for(int d=0; d<dim(); d++) {
	 res(d) = ll(d) + (2*rad()) * (pass(d)+0.5) / double(tmp);
  }
  return res;
}
double Let3d::radius(int gni) //radius of a node
{
  return rad()/double(pow2(depth(gni)));
}
// ---------------------------------------------------------------------- 
int Let3d::findgnt(int wntdepth, const Index3& wntpath)
{
  int cur = 0;  //cerr<<"GOOD "<<path(cur)<<"     ";
  Index3 leftpath(wntpath);
  while(depth(cur)<wntdepth && terminal(cur)==false) {
	 int dif = wntdepth-depth(cur);
	 int tmp = pow2(dif-1);
	 Index3 choice( leftpath/tmp );
	 leftpath -= choice*tmp;
	 cur = child(cur, choice);	 //cur = child(cur, IdxIter(0,2,true,choice) );	 //cerr<<path(cur)<<"["<<choice<<" "<<tmp<<"]"<<"     ";
  }  //cerr<<endl;
  return cur;
}
// ---------------------------------------------------------------------- 
bool Let3d::adjacent(int me, int yo)
{
  int md = max(depth(me),depth(yo));
  Index3 one(1);
  Index3 mectr(  (2*path(me)+one) * pow2(md - depth(me))  );
  Index3 yoctr(  (2*path(yo)+one) * pow2(md - depth(yo))  );
  int merad = pow2(md - depth(me));
  int yorad = pow2(md - depth(yo));
  Index3 dif( abs(mectr-yoctr) );
  int rad  = merad + yorad;
  return
	 ( dif <= Index3(rad) ) && //not too far
	 ( dif.linfty() == rad ); //at least one edge touch
}

