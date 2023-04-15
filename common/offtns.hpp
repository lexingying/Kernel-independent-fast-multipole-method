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
#ifndef _OFFTNS_HPP_
#define _OFFTNS_HPP_

#include "offmat.hpp"

template <class F>
class OffTns
{
public:
  int _m, _n, _p;
  int _s, _t, _u;
  bool _owndata;
  F* _data;
public:
  OffTns(int m=0, int n=0, int p=0, int s=0, int t=0, int u=0): _m(m), _n(n), _p(p), _s(s), _t(t), _u(u), _owndata(true) {
	 if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
  }
  OffTns(int m, int n, int p, int s, int t, int u, bool owndata, F* data): _m(m), _n(n), _p(p), _s(s), _t(t), _u(u), _owndata(owndata) {
	 if(_owndata) {
		if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
		if(_m>0 && _n>0 && _p>0) { for(int i=0; i<_m*_n*_p; i++) _data[i] = data[i]; }
	 } else {
		_data = data;
	 }
  }
  OffTns(const OffTns& C): _m(C._m), _n(C._n), _p(C._p), _s(C._s), _t(C._t), _u(C._u), _owndata(C._owndata) {
	 if(_owndata) {
		if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
		if(_m>0 && _n>0 && _p>0) { for(int i=0; i<_m*_n*_p; i++) _data[i] = C._data[i]; }
	 } else {
		_data = C._data;
	 }
  }
  ~OffTns() { 
	 if(_owndata) { 
		if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
	 }
  }
  OffTns& operator=(const OffTns& C) {
	 if(_owndata) { 
		if(_m>0 && _n>0 && _p>0) { delete[] _data; _data = NULL; } 
	 }
	 _m = C._m; _n=C._n; _p=C._p;	_s = C._s; _t=C._t; _u=C._u; _owndata=C._owndata;
	 if(_owndata) {
		if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
		if(_m>0 && _n>0 && _p>0) { for(int i=0; i<m*n*p; i++) _data[i] = C._data[i]; }
	 } else {
		_data = C._data;
	 }
	 return *this;
  }
  void resize(int m, int n, int p, int s, int t, int u)  {
	 assert( _owndata==true );
	 if(_m!=m || _n!=n || _p!=p) {
		if(_m>0 && _n>0 && _p>0) { delete[] _data; _data = NULL; } 
		_m = m; _n = n; _p=p; _s=s; _t=t; _u=u;
		if(_m>0 && _n>0 && _p>0) { _data = new F[_m*_n*_p]; assert( _data!=NULL ); } else _data=NULL;
	 }
  }
  const F& operator()(int i, int j, int k) const  {
	 assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t && k>=_u && k<_p+_u);
	 return _data[(i-_s) + (j-_t)*_m + (k-_u)*_m*_n];
  }
  F& operator()(int i, int j, int k)  {
	 assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t && k>=_u && k<_p+_u);
	 return _data[(i-_s) + (j-_t)*_m + (k-_u)*_m*_n];
  }
  
  F* data() const { return _data; }
  int m() const { return _m; }
  int n() const { return _n; }
  int p() const { return _p; }
  int s() const { return _s; }
  int t() const { return _t; }
  int u() const { return _u; }
};

template <class F> inline ostream& operator<<( ostream& os, const OffTns<F>& tns)
{
  os<<tns.m()<<" "<<tns.n()<<" "<<tns.p()<<" "<<tns.s()<<" "<<tns.t()<<" "<<tns.u()<<endl;
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=tns.s(); i<tns.s()+tns.m(); i++) {
	 for(int j=tns.t(); j<tns.t()+tns.n(); j++) {
		for(int k=tns.u(); k<tns.u()+tns.p(); k++) {
		  os<<" "<<tns(i,j,k);
		}
		os<<endl;
	 }
	 os<<endl;
  }
  return os;
}
template <class F> inline void setvalue(OffTns<F>& M, F val)
{
  for(int i=0; i<M.m(); i++)
	 for(int j=0; j<M.n(); j++)
		for(int k=0; k<M.p(); k++)
		  M(i,j,k) = val;
  return;
}

typedef OffTns<bool>   BolOffTns;
typedef OffTns<int>    IntOffTns;
typedef OffTns<double> DblOffTns;
typedef OffTns<cpx>    CpxOffTns;

#endif




