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
/*! \file */
#ifndef _NUMVEC_HPP_
#define _NUMVEC_HPP_

#include "commoninc.hpp"

using std::ostream;
using std::ios_base;
using std::endl;

template <class F>
class NumVec
{
public:
  int  _m;
  bool _owndata;
  F* _data;
public:
  NumVec(int m=0): _m(m), _owndata(true)  {
	 if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
  }
  NumVec(int m, bool owndata, F* data): _m(m), _owndata(owndata) {
	 if(_owndata) {
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
		if(_m>0) memcpy( _data, data, _m*sizeof(F) );
	 } else {
		_data = data;
	 }
  }
  NumVec(const NumVec& C): _m(C._m), _owndata(C._owndata)  {
	 if(_owndata) {
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
		if(_m>0) memcpy( _data, C._data, _m*sizeof(F) );
	 } else {
		_data = C._data;
	 }
  }
  ~NumVec() {
	 if(_owndata) {
		if(_m>0) { delete[] _data; _data = NULL; }
	 }
  }
  NumVec& operator=(const NumVec& C)  {
	 if(_owndata) { 
		if(_m>0) { delete[] _data; _data = NULL; }
	 }
	 _m = C._m; _owndata=C._owndata;
	 if(_owndata) {
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
		if(_m>0) memcpy( _data, C._data, _m*sizeof(F) );
	 } else {
		_data =C._data;
	 }
	 return *this;
  }
  void resize(int m)  {
	 assert(_owndata==true);
	 if(m !=_m) {
		if(_m>0) { delete[] _data; _data = NULL; }
		_m = m;
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL);  memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
	 }
  }
  const F& operator()(int i) const  {	 assert(i>=0 && i<_m);
	 return _data[i]; 
  }
  F& operator()(int i)  {	 assert(i>=0 && i<_m);
	 return _data[i]; 
  }
  
  F* data() const { return _data; }
  int m () const { return _m; }
};

template <class F> inline ostream& operator<<( ostream& os, const NumVec<F>& vec)
{
  os<<vec.m()<<endl;
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=0; i<vec.m(); i++)	 os<<" "<<vec(i);
  os<<endl;
  return os;
}
template <class F> inline void setvalue(NumVec<F>& V, F val)
{
  for(int i=0; i<V.m(); i++)
	 V(i) = val;
}
template <class F> inline void clear(NumVec<F>& V)
{
  memset(V.data(), 0, V.m()*sizeof(F));
}

typedef NumVec<bool>   BolNumVec;
typedef NumVec<int>    IntNumVec;
typedef NumVec<double> DblNumVec; //typedef NumVec<double> SclNumVec;


#endif


