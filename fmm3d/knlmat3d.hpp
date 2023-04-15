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
#ifndef _KNLMAT3D_HPP_
#define _KNLMAT3D_HPP_

#include "common/vec3t.hpp"
#include "comobject.hpp"
#include "kernel3d.hpp"

//----------------------------------------------------------------------------------
class KnlMat3d: public ComObject
{
protected:
  //---PARAMS (REQ)
  DblNumMat* _srcpos; //src position
  DblNumMat* _srcnor; //src normal
  DblNumMat* _trgpos; //trg position
  Kernel3d _knl;
public:
  KnlMat3d(const string& p):  ComObject(p), _srcpos(NULL), _srcnor(NULL), _trgpos(NULL) {;}
  virtual ~KnlMat3d() { }
  //MEMBER ACESS
  DblNumMat*& srcpos() { return _srcpos; }
  DblNumMat*& srcnor() { return _srcnor; }
  DblNumMat*& trgpos() { return _trgpos; }
  Kernel3d& knl()    { return _knl; }
  //SETUP and USE
  virtual int setup(map<string, string>&)=0;
  virtual int eval(const DblNumVec& srcden, DblNumVec& trgval) = 0;
  //OTHER ACCESS
  int dim() { return 3; }
  int sdof() { return _knl.sdof(); }
  int tdof() { return _knl.tdof(); }
};

#endif


