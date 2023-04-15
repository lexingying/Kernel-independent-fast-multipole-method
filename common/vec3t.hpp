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
#ifndef  _VEC3T_HPP_
#define  _VEC3T_HPP_

#include "commoninc.hpp"

using std::istream;
using std::ostream;
using std::min;
using std::max;
using std::abs;

///Common VECtor Template
template <class F>
class Vec3T {
private:
  F _v[3];
public:
  enum{ X=0, Y=1, Z=2 };
  //------------CONSTRUCTOR AND DESTRUCTOR 
  Vec3T()              { _v[0]=F(0);    _v[1]=F(0);    _v[2]=F(0); }
  Vec3T(F f)           { _v[0]=f;       _v[1]=f;       _v[2]=f;}
  Vec3T(const F* f)    { _v[0]=f[0];    _v[1]=f[1];    _v[2]=f[2]; }
  Vec3T(F a,F b,F c)   { _v[0]=a;       _v[1]=b;       _v[2]=c; }
  Vec3T(const Vec3T& c){ _v[0]=c._v[0]; _v[1]=c._v[1]; _v[2]=c._v[2]; }
  ~Vec3T() {}
  //------------POINTER and ACCESS
  operator F*()             { return &_v[0]; }
  operator const F*() const { return &_v[0]; }
  F* array()                { return &_v[0]; }  //access array
  F& operator()(int i)             { assert(i<3); return _v[i]; }
  const F& operator()(int i) const { assert(i<3); return _v[i]; }
  F& operator[](int i)             { assert(i<3); return _v[i]; }
  const F& operator[](int i) const { assert(i<3); return _v[i]; }
  F& x()             { return _v[0];}
  F& y()             { return _v[1];}
  F& z()             { return _v[2];}
  const F& x() const { return _v[0];}
  const F& y() const { return _v[1];}
  const F& z() const { return _v[2];}
  //------------ASSIGN
  Vec3T& operator= ( const Vec3T& c ) { _v[0] =c._v[0]; _v[1] =c._v[1]; _v[2] =c._v[2]; return *this; }
  Vec3T& operator+=( const Vec3T& c ) { _v[0]+=c._v[0]; _v[1]+=c._v[1]; _v[2]+=c._v[2]; return *this; }
  Vec3T& operator-=( const Vec3T& c ) { _v[0]-=c._v[0]; _v[1]-=c._v[1]; _v[2]-=c._v[2]; return *this; }
  Vec3T& operator*=( const F& s )     { _v[0]*=s;       _v[1]*=s;       _v[2]*=s;       return *this; }
  Vec3T& operator/=( const F& s )     { _v[0]/=s;       _v[1]/=s;       _v[2]/=s;       return *this; }
  //-----------LENGTH...
  F l1( void )     const  { F sum=F(0); for(int i=0; i<3; i++) sum=sum+abs(_v[i]); return sum; }
  F linfty( void ) const  { F cur=F(0); for(int i=0; i<3; i++) cur=max(cur,abs(_v[i])); return cur; }
  F l2( void )     const  { F sum=F(0); for(int i=0; i<3; i++) sum=sum+_v[i]*_v[i]; return sqrt(sum); }
  F length( void ) const  { return l2(); }
  Vec3T dir( void )    const  { F a=l2(); return (*this)/a; }
};

//-----------BOOLEAN OPS
template <class F> inline bool operator==(const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)   res = res && (a(i)==b(i));  return res;
}
template <class F> inline bool operator!=(const Vec3T<F>& a, const Vec3T<F>& b) {
  return !(a==b);
}
template <class F> inline bool operator> (const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)   res = res && (a(i)> b(i));  return res; 
}
template <class F> inline bool operator< (const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)   res = res && (a(i)< b(i));  return res; 
}
template <class F> inline bool operator>=(const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)	res = res && (a(i)>=b(i));  return res; 
}
template <class F> inline bool operator<=(const Vec3T<F>& a, const Vec3T<F>& b) {
  bool res = true;  for(int i=0; i<3; i++)   res = res && (a(i)<=b(i));  return res; 
}

//-----------NUMERICAL OPS
template <class F> inline Vec3T<F> operator- (const Vec3T<F>& a) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = -a[i]; return r;
}
template <class F> inline Vec3T<F> operator+ (const Vec3T<F>& a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]+b[i]; return r; 
}
template <class F> inline Vec3T<F> operator- (const Vec3T<F>& a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]-b[i]; return r;
}
template <class F> inline Vec3T<F> operator* (F scl, const Vec3T<F>& a) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec3T<F> operator* (const Vec3T<F>& a, F scl) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec3T<F> operator/ (const Vec3T<F>& a, F scl) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]/scl;  return r;
}

template <class F> inline F operator* (const Vec3T<F>& a, const Vec3T<F>& b) {
  F sum=F(0); for(int i=0; i<3; i++) sum=sum+a(i)*b(i); return sum;
}
template <class F> inline F dot       (const Vec3T<F>& a, const Vec3T<F>& b) {
  return a*b;
}
template <class F> inline Vec3T<F> operator^ (const Vec3T<F>& a, const Vec3T<F>& b) {
  return Vec3T<F>(a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0)); 
}
template <class F> inline Vec3T<F> cross     (const Vec3T<F>& a, const Vec3T<F>& b) { 
  return a^b; 
}

//-------------ew OPS
template <class F> inline Vec3T<F> min(const Vec3T<F>& a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = min(a[i], b[i]); return r;
}
template <class F> inline Vec3T<F> max(const Vec3T<F>& a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = max(a[i], b[i]); return r;
}
template <class F> inline Vec3T<F> abs(const Vec3T<F>& a) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = abs(a[i]); return r;
}
template <class F> inline Vec3T<F> ewmul(const Vec3T<F>&a, const Vec3T<F>& b) {
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]*b[i]; return r;
}
template <class F> inline Vec3T<F> ewdiv(const Vec3T<F>&a, const Vec3T<F>& b) { 
  Vec3T<F> r;  for(int i=0; i<3; i++) r[i] = a[i]/b[i]; return r;
}
//---------------INOUT
template <class F> istream& operator>>(istream& is, Vec3T<F>& a) {
  for(int i=0; i<3; i++) is>>a[i]; return is;
}
template <class F> ostream& operator<<(ostream& os, const Vec3T<F>& a) { 
  for(int i=0; i<3; i++) os<<a[i]<<" "; return os;
}

//---------------------------------------------------------
/// MOST COMMONLY USED
typedef Vec3T<double> Point3;
typedef Vec3T<int>    Index3;



#endif
