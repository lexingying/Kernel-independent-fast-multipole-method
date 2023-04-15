#ifndef _MEXAUX_HPP_
#define _MEXAUX_HPP_

#include "mex.h"
#include "matrix.h"

#include "fmm3d.hpp"

using std::cerr;
using std::endl;

//integer
inline void mex2cpp(const mxArray*& md, int& cd);
inline void cpp2mex(const int& cd, mxArray*& md);
//double
inline void mex2cpp(const mxArray*& md, double& cd);
inline void cpp2mex(const double& cd, mxArray*& md);
//dblnummat
inline void mex2cpp(const mxArray*& md, DblNumMat& cd);
inline void cpp2mex(const DblNumMat& cd, mxArray*& md);
//dblnumvec
inline void mex2cpp(const mxArray*& md, DblNumVec& cd);
inline void cpp2mex(const DblNumVec& cd, mxArray*& md);
//point3
inline void mex2cpp(const mxArray*& md, Point3& cd);
inline void cpp2mex(const Point3& cd, mxArray*& md);
//Kernel3d
inline void mex2cpp(const mxArray*& md, Kernel3d& cd);
inline void cpp2mex(const Kernel3d& cd, mxArray*& md);
//vector<double>
inline void mex2cpp(const mxArray*& md, vector<double>& cd);
inline void cpp2mex(const vector<double>& cd, mxArray*& md);
//vector<T>
template <class T> inline void mex2cpp(const mxArray*& md, vector<T>& cd);
template <class T> inline void cpp2mex(const vector<T>& cd, mxArray*& md);

//----------------------int
inline void mex2cpp(const mxArray*& md, int& cd)
{
  cd = int(mxGetScalar(md));
  return;
}
inline void cpp2mex(const int& cd, mxArray*& md)
{
  md = mxCreateDoubleScalar(cd);
  return;
}
//----------------------double
inline void mex2cpp(const mxArray*& md, double& cd)
{
  cd = mxGetScalar(md);
  return;
}
inline void cpp2mex(const double& cd, mxArray*& md)
{
  md = mxCreateDoubleScalar(cd);
  return;
}
//----------------------DblNumMat
inline void mex2cpp(const mxArray*& md, DblNumMat& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);
  double* xr = mxGetPr(md);
  cd.resize(m,n);
  int cnt = 0;
  for(int j=0; j<n; j++)
	 for(int i=0; i<m; i++) {
		cd(i,j) = xr[cnt];
		cnt++;
	 }
  return;
}
inline void cpp2mex(const DblNumMat& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  md = mxCreateDoubleMatrix(m, n, mxREAL);
  double* xr = mxGetPr(md);
  int cnt = 0;
  for(int j=0; j<n; j++)
	 for(int i=0; i<m; i++) {
		xr[cnt] = cd(i,j);
		cnt++;
	 }
  return;
}
//----------------------DblNumMat
inline void mex2cpp(const mxArray*& md, DblNumVec& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md); iA(n==1);
  double* xr = mxGetPr(md);
  cd.resize(m);
  for(int i=0; i<m; i++)
	 cd(i) = xr[i];
  return;
}
inline void cpp2mex(const DblNumVec& cd, mxArray*& md)
{
  int m = cd.m();
  md = mxCreateDoubleMatrix(m, 1, mxREAL);
  double* xr = mxGetPr(md);
  for(int i=0; i<m; i++)
	 xr[i] = cd(i);
  return;
}
//----------------------Point3
inline void mex2cpp(const mxArray*& md, Point3& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);  iA(m==3 && n==1);
  double* xr = mxGetPr(md);
  cd = Point3(xr[0], xr[1], xr[2]);
  return;
}
inline void cpp2mex(const Point3& cd, mxArray*& md)
{
  md = mxCreateDoubleMatrix(3, 1, mxREAL);
  double* xr = mxGetPr(md);
  xr[0] = cd(0);  xr[1] = cd(1);  xr[2] = cd(2);
  return;
}
//----------------------Kernel3d
inline void mex2cpp(const mxArray*& md, Kernel3d& cd)
{
  int m = mxGetM(md); iA(m==2);
  int n = mxGetN(md); iA(n==1);  //cerr<<m<<" "<<n<<endl;
  { const mxArray* tt = mxGetCell(md, 0);  mex2cpp(tt, cd.kt()); }
  { const mxArray* tt = mxGetCell(md, 1);  mex2cpp(tt, cd.coefs()); }
  return;
}
inline void cpp2mex(const Kernel3d& cd, mxArray*& md)
{
  int m = 2;
  int n = 1;
  md = mxCreateCellMatrix(m,n);
  mxArray* ss;
  cpp2mex(cd.kt(), ss);  mxSetCell(md, 0, ss);
  cpp2mex(cd.coefs(), ss);  mxSetCell(md, 1, ss);
  return;
}
//----------------------vector<int>
inline void mex2cpp(const mxArray*& md, vector<double>& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md); iA(n==1);  //cerr<<m<<"   "<<n<<endl;
  double* xr = mxGetPr(md);
  cd.resize(m*n);
  for(int i=0; i<m*n; i++)
	 cd[i] = xr[i];
  return;
}
inline void cpp2mex(const vector<double>& cd, mxArray*& md)
{
  int m = cd.size();
  int n = 1;
  md = mxCreateDoubleMatrix(m, n, mxREAL);
  double* xr = mxGetPr(md);
  for(int i=0; i<m*n; i++)
	 xr[i] = cd[i];
  return;
}
//----------------------vector<...>
template <class T> inline void mex2cpp(const mxArray*& md, vector<T>& cd)
{
  int m = mxGetM(md); 
  int n = mxGetN(md);  assert(n==1);
  cd.resize(m*n);
  for(int ci=0; ci<m*n; ci++) {
	 const mxArray*tt = mxGetCell(md, ci);
	 mex2cpp(tt, cd[ci]);
  }
  return;
}
template <class T> inline void cpp2mex(const vector<T>& cd, mxArray*& md)
{
  int n = cd.size();
  md = mxCreateCellMatrix(n, 1);
  for(int ci=0; ci<n; ci++) {
	 mxArray* ss;	 cpp2mex(cd[ci], ss);
	 mxSetCell(md, ci, ss);
  }
  return;
}

#endif
