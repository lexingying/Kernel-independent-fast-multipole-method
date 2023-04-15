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
#ifndef _LET3D_HPP_
#define _LET3D_HPP_

#include "common/vec3t.hpp"
#include "common/nummat.hpp"
#include "comobject.hpp"

using std::vector;

enum {
  LET_SRCNODE = 1,
  LET_TRGNODE = 2
};

//---------------------------------------------------------------------------
// unique identifier: a set of points and its bounding box 
class Let3d: public ComObject
{
public:
  //---------------------------------------
  class Node {
  protected:
	 int _par, _chd;
	 Index3 _path;
	 int _depth;
	 int _tag; //s empty, t empty ...
	 
	 int _sndeidx, _sextbeg, _sextnum; //index for the source nodes
	 vector<int> _sownvis;
	 
	 int _tndeidx, _textbeg, _textnum; //index for the target nodes
	 vector<int> _townvis;
	 
	 vector<int> _Unodes, _Vnodes, _Wnodes, _Xnodes;
	 
  public:
	 Node(int p, int c, Index3 t, int d):
		_par(p), _chd(c), _path(t), _depth(d), _tag(false),
		_sndeidx(0), _sextbeg(0), _sextnum(0),
		_tndeidx(0), _textbeg(0), _textnum(0) {;}
	 int& par() { return _par; }
	 int& chd() { return _chd; }
	 Index3& path() { return _path; }
	 int& depth() { return _depth; }
	 int& tag() { return _tag; }
	 
	 int& sndeidx() { return _sndeidx; }	 int& sextbeg() { return _sextbeg; }	 int& sextnum() { return _sextnum; }
	 vector<int>& sownvis() { return _sownvis; }
	 
	 int& tndeidx() { return _tndeidx; }	 int& textbeg() { return _textbeg; }	 int& textnum() { return _textnum; }
	 vector<int>& townvis() { return _townvis; }
	 
	 vector<int>& Unodes() { return _Unodes; }
	 vector<int>& Vnodes() { return _Vnodes; }
	 vector<int>& Wnodes() { return _Wnodes; }
	 vector<int>& Xnodes() { return _Xnodes; }
  };
  //----------------------------------------------
protected:
  //PARAMS(REQ)
  DblNumMat* _srcpos;
  DblNumMat* _trgpos;
  Point3 _ctr;
  int _rootlvl;
  //PARAMS(OPT)
  int _ptsmax;
  int _maxlevel; //AT MOST 16
  //COMPONENTS
  vector<Node> _ndevec;
  int _level; //total number of levels
  int _sndecnt, _sextcnt;
  int _tndecnt, _textcnt;
public:
  Let3d(const string& p);
  ~Let3d();
  //MEMBER ACCESS
  DblNumMat*& srcpos() { return _srcpos; }
  DblNumMat*& trgpos() { return _trgpos; }
  Point3& ctr() { return _ctr; }
  int& rootlvl() { return _rootlvl; }
  double rad() { return pow(2.0, -_rootlvl); }
  int& ptsmax() { return _ptsmax; }
  int& maxlevel() { return _maxlevel; }
  //SETUP AND USE
  int setup(map<string,string>& opts);
  int print(const char*);
  //access
  vector<Node>& ndevec() { return _ndevec; }
  int sndecnt() { return _sndecnt; }
  int sextcnt() { return _sextcnt; }
  int tndecnt() { return _tndecnt; }
  int textcnt() { return _textcnt; }
  //construction
  int srcdata();
  int trgdata();
  //---------------------------------------------
  //LOCAL
  int calgnext(int gni); //build U,V,W,X lists
  int dnOrderCollect(vector<int>&); //top down ordering of the nodes
  int upOrderCollect(vector<int>&); //bottom up ordering of the nodes
  //node access
  Node& node(int gni) { return _ndevec[gni]; }
  //tree traversal and informational retrieval
  //GNTra   gni2gnt(int gni);  //int     gnt2gni(GNTra gnt) { return gnt.gni(); }
  bool    root(int gni)     { return node(gni).par()==-1; } //no par
  bool    terminal(int gni) { return node(gni).chd()==-1; } //no chd
  int     parent(int gni)   { assert(node(gni).par()!=-1); return node(gni).par(); }
  int     child( int gni, const Index3& idx); //LEXING { assert(node(gni).chd()!=-1); return node(gni).chd()+chdi2o(idx); }
  Index3  path(int gni)     { return node(gni).path(); }
  int     depth(int gni)    { return node(gni).depth(); }  //LEXING: int     child( int gni, int ord) { return node(gni).chd() + ord; }
  int     tag(int gni)      { return node(gni).tag(); }
  
  Point3  center(int gni);  //LEXING: { GNTra gnt=gni2gnt(gni); return center(gnt); }
  double  radius(int gni);  //LEXING: { GNTra gnt=gni2gnt(gni); return radius(gnt); }
  
  int findgnt(int depth, const Index3& path);  //LEXING: GNTra   findgnt(int depth, const Index3& path); //
  bool adjacent(int, int);                //LEXING: bool    adjacent(GNTra, GNTra);
  
  int dim() const { return 3; }  //int chdi2o(Index3 idx); //child index to order  //Index3 chdo2i(int ord); //child order to index
};
  
//extra stuff
//int  plclnum(Vec pos) { int tmp; VecGetLocalSize(pos, &tmp); return tmp/dim(); }
//int  pglbnum(Vec pos) { int tmp; VecGetSize(     pos, &tmp); return tmp/dim(); }
//void plclran(Vec pos, int& beg, int& end) { VecGetOwnershipRange(pos, &beg, &end); beg=beg/dim(); end=end/dim(); }

//bool    root(GNTra gnt)     { return root(gnt.gni()); }
//bool    terminal(GNTra gnt) { return terminal(gnt.gni()); }
//GNTra   parent(GNTra);
//GNTra   child( GNTra, const Index3& ii);
//Point3  center(GNTra);
//double  radius(GNTra);




#endif
