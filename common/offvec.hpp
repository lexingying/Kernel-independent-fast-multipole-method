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
#ifndef _OFFVEC_HPP_
#define _OFFVEC_HPP_

#include "commoninc.hpp"

template <class F>
class OffVec
{
public:
  int  _m;
  int  _s;
  bool _owndata;
  F* _data;
public:
  OffVec(int m=0, int s=0): _m(m), _s(s), _owndata(true)  {
	 if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
  }
  OffVec(int m, int s, bool owndata, F* data): _m(m), _s(s), _owndata(owndata) {
	 if(_owndata) {
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
		if(_m>0) memcpy( _data, data, _m*sizeof(F) );
	 } else {
		_data = data;
	 }
  }
  OffVec(const OffVec& C): _m(C._m), _s(C._s), _owndata(C._owndata)  {
	 if(_owndata) {
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
		if(_m>0) memcpy( _data, C._data, _m*sizeof(F) );
	 } else {
		_data = C._data;
	 }
  }
  ~OffVec() {
	 if(_owndata) {
		if(_m>0) { delete[] _data; _data = NULL; }
	 }
  }
  OffVec& operator=(const OffVec& C)  {
	 if(_owndata) { 
		if(_m>0) { delete[] _data; _data = NULL; }
	 }
	 _m=C._m; _s=C._s; _owndata=C._owndata;
	 if(_owndata) {
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
		if(_m>0) memcpy( _data, C._data, _m*sizeof(F) );
	 } else {
		_data =C._data;
	 }
	 return *this;
  }
  void resize(int m, int s)  {
	 assert(_owndata==true);
	 if(m !=_m) {
		if(_m>0) { delete[] _data; _data = NULL; }
		_m = m;		_s = s;
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL);  memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
	 }
  }
  const F& operator()(int i) const  {
	 assert(i>=_s && i<_m+_s);
	 return _data[i-_s]; 
  }
  F& operator()(int i)  {
	 assert(i>=_s && i<_m+_s);
	 return _data[i-_s]; 
  }
  
  F* data() const { return _data; }
  int m() const { return _m; }
  int s() const { return _s; }
};

template <class F> inline ostream& operator<<( ostream& os, const OffVec<F>& vec)
{
  os<<vec.m()<<" "<<vec.s()<<endl;
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=vec.s(); i<vec.m()+vec.s(); i++)	 os<<" "<<vec(i);
  os<<endl;
  return os;
}
template <class F> inline void setvalue(OffVec<F>& vec, F val)
{
  for(int i=vec.s(); i<vec.s()+vec.m(); i++)	 vec(i) = val;
}
template <class F> inline void clear(OffVec<F>& vec)
{
  memset(vec.data(), 0, vec.m()*sizeof(F));
}

typedef OffVec<bool>   BolOffVec;
typedef OffVec<int>    IntOffVec;
typedef OffVec<double> DblOffVec;


#endif


