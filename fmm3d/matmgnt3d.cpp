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
#include "matmgnt3d.hpp"
#include "common/vecmatop.hpp"

using std::cerr;
using std::endl;

// ---------------------------------------------------------------------- 
double MatMgnt3d::_wsbuf[16384];

// ---------------------------------------------------------------------- 
MatMgnt3d::MatMgnt3d(): _np(6), _forplan(NULL), _invplan(NULL)
{
}
// ---------------------------------------------------------------------- 
MatMgnt3d::~MatMgnt3d()
{
  //-------------------------------------------------------
  if(_forplan!=NULL) {	 rfftwnd_destroy_plan(_forplan); _forplan=NULL;  }
  if(_invplan!=NULL) { 	 rfftwnd_destroy_plan(_invplan); _invplan=NULL;  }  //cerr<<"matmgnt3d free"<<endl;
}
double MatMgnt3d::alt()
{
  return pow(0.1, _np+1);
}
// ---------------------------------------------------------------------- 
// ---------------------------------------------------------------------- 
int MatMgnt3d::setup()
{
  //--------------------------------------------------------
  _hom = _knl.hom();
  if(_hom==true) {
	 _knl.homsdeg(_degvec); iA(_degvec.size()==sdof());
  }
  iC( splposcal(_np,   1.0, _splpos[UE]) );
  iC( splposcal(_np+2, 3.0, _splpos[UC]) );
  iC( splposcal(_np,   3.0, _splpos[DE]) );
  iC( splposcal(_np,   1.0, _splpos[DC]) );
  
  iC( regposcal(_np,   1.0, _regpos    ) ); //only one regpos
  
  _forplan = rfftw3d_create_plan(2*_np,2*_np,2*_np, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
  _invplan = rfftw3d_create_plan(2*_np,2*_np,2*_np, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
  //--------------------------------------------------------
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::report()
{
  cerr<<"matrix map size"<<endl;
  cerr<<_uc2ue.size()<<" "<<_ue2uc.size()<<" "<<_dc2de.size()<<" "<<_de2dc.size()<<" "<<_ue2dc.size()<<endl;
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::plndatasize(int tp)
{
  if(tp==UE || tp==DE)
	 return _splpos[tp].n()*sdof();
  else
	 return _splpos[tp].n()*tdof();
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::effdatasize(int tp)
{
  //int et = eq().et();
  int effnum = (2*_np+2)*(2*_np)*(2*_np);
  if(tp==UE || tp==DE)
	 return effnum*sdof();
  else
	 return effnum*tdof();
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::UC2UE_dgemv(int l, const DblNumVec& chk, DblNumVec& den)
{
  DblNumMat& _UC2UE = (_hom==true) ? _uc2ue[0] : _uc2ue[l];
  double R = (_hom==true) ? 1 : 1.0/pow(2.0,l);
  //---------compute matrix
  if(_UC2UE.m()==0) {	 //cerr<<"UC2UE compute"<<endl;
	 //set matrix
	 DblNumMat ud2c(plndatasize(UC), plndatasize(UE));
	 DblNumMat chkpos(dim(),splpos(UC).n());	 clear(chkpos);	 iC( daxpy(R, splpos(UC), chkpos) ); //scale
	 DblNumMat denpos(dim(),splpos(UE).n());	 clear(denpos);	 iC( daxpy(R, splpos(UE), denpos) ); //scale
	 
	 iC( _knl.kernel(denpos, denpos, chkpos, ud2c) );
	 _UC2UE.resize(plndatasize(UE), plndatasize(UC)); //_memused[0] += plnnum(UE)*dof()*plnnum(UC)*dof()*sizeof(double);
	 iC( pinv(ud2c, alt(), _UC2UE) );
  }
  //---------matvec
  if(_hom==true) {
	 //matvec
	 int sdof = this->sdof();
	 DblNumVec tmpden(sdof*splpos(UE).n(), false, _wsbuf);	 clear(tmpden);
	 iC( dgemv(1.0, _UC2UE, chk, 1.0, tmpden) );
	 //scale
	 vector<double> sclvec(sdof);	 for(int s=0; s<sdof; s++)		sclvec[s] = pow(2.0, - l*_degvec[s]);
	 int cnt = 0;
	 for(int i=0; i<splpos(UE).n(); i++)
		for(int s=0; s<sdof; s++) {
		  den(cnt) = den(cnt) + tmpden(cnt) * sclvec[s];		  cnt++;
		}
  } else {
	 iC( dgemv(1.0, _UC2UE, chk, 1.0, den) );
  }
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::UE2UC_dgemv(int l, Index3 idx, const DblNumVec& den, DblNumVec& chk)
{
  OffTns<DblNumMat>& _UE2UC = (_hom==true) ? _ue2uc[0] : _ue2uc[l];
  double R = (_hom==true) ? 1 : 1.0/pow(2.0, l);  //OffTns<DblNumMat>& _UE2UC = _ue2uc[l];  double R       = 1.0/pow(2.0, l);
  if(_UE2UC.m()==0)	 _UE2UC.resize(2,2,2,0,0,0);
  DblNumMat& _UE2UCii = _UE2UC(idx(0), idx(1), idx(2));
  //---------compute matrix
  if(_UE2UCii.m()==0) {	 //cerr<<"UE2UC compute"<<endl;
	 _UE2UCii.resize(plndatasize(UC), plndatasize(UE)); //_memused[1] += plnnum(UC)*dof()*plnnum(UE)*dof()*sizeof(double);
	 DblNumMat chkpos(dim(),splpos(UC).n());	 clear(chkpos);	 iC( daxpy(2.0*R, splpos(UC), chkpos) ); //scale
	 DblNumMat denpos(dim(),splpos(UE).n());	 clear(denpos);	 iC( daxpy(R, splpos(UE), denpos) ); //scale
	 for(int i=0; i<dim(); i++) for(int j=0; j<splpos(UE).n(); j++)	denpos(i,j) = denpos(i,j) + (2*idx(i)-1)*R;//shift
	 
	 iC( _knl.kernel(denpos, denpos, chkpos, _UE2UCii) );
  }
  //---------matvec
  if(_hom==true) {
	 int sdof = this->sdof();
	 DblNumVec tmpden(sdof*splpos(UE).n(), false, _wsbuf);	 clear(tmpden);
	 vector<double> sclvec(sdof);	 for(int s=0; s<sdof; s++)		sclvec[s] = pow(2.0, l*_degvec[s]);
	 int cnt = 0;
	 for(int i=0; i<splpos(UE).n(); i++)
		for(int s=0; s<sdof; s++) {
		  tmpden(cnt) = den(cnt) * sclvec[s];		  cnt++;
		}
	 iC( dgemv(1.0, _UE2UCii, tmpden, 1.0, chk) );
  } else {
	 iC( dgemv(1.0, _UE2UCii, den, 1.0, chk) );
  }
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::DC2DE_dgemv(int l, const DblNumVec& chk, DblNumVec& den)
{
  DblNumMat& _DC2DE = (_hom==true) ? _dc2de[0]: _dc2de[l];
  double R = (_hom==true) ? 1 : 1.0/pow(2.0,l);
  //---------compute matrix
  if(_DC2DE.m()==0) {	 //cerr<<"DC2DE compute"<<endl;
	 DblNumMat dd2c(plndatasize(DC), plndatasize(DE));
	 DblNumMat chkpos(dim(),splpos(DC).n());		clear(chkpos);	 iC( daxpy(R, splpos(DC), chkpos) ); //scale
	 DblNumMat denpos(dim(),splpos(DE).n());		clear(denpos);	 iC( daxpy(R, splpos(DE), denpos) ); //scale
	 
	 iC( _knl.kernel(denpos, denpos, chkpos, dd2c) );//matrix
	 _DC2DE.resize(plndatasize(DE), plndatasize(DC)); //_memused[2] += plndnenum()*dof()*plndncnum()*dof()*sizeof(double);
	 iC( pinv(dd2c, alt(), _DC2DE) );
  }
  //---------matvec
  if(_hom==true) {
	 int sdof = this->sdof();
	 DblNumVec tmpden(sdof*splpos(DE).n(), false, _wsbuf);	 clear(tmpden);
	 iC( dgemv(1.0, _DC2DE, chk, 1.0, tmpden) );
	 //scale
	 vector<double> sclvec(sdof);	 for(int s=0; s<sdof; s++)		sclvec[s] = pow(2.0, - l*_degvec[s]);
	 int cnt = 0;
	 for(int i=0; i<splpos(DE).n(); i++)
		for(int s=0; s<sdof; s++) {
		  den(cnt) = den(cnt) + tmpden(cnt) * sclvec[s];		  cnt++;
		}
  } else {
	 iC( dgemv(1.0, _DC2DE, chk, 1.0, den) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::DE2DC_dgemv(int l, Index3 idx, const DblNumVec& den, DblNumVec& chk)
{
  OffTns<DblNumMat>& _DE2DC = (_hom==true) ? _de2dc[0] : _de2dc[l];
  double R = (_hom==true) ? 1 : 1.0/pow(2.0, l);  //OffTns<DblNumMat>& _DE2DC = _de2dc[l];  double R       = 1.0/pow(2.0, l);
  if(_DE2DC.m()==0)	 _DE2DC.resize(2,2,2,0,0,0);
  DblNumMat& _DE2DCii = _DE2DC(idx[0], idx[1], idx[2]);
  
  //---------compute matrix
  if(_DE2DCii.m()==0) {	 //cerr<<"DE2DC compute"<<endl;
	 _DE2DCii.resize(plndatasize(DC), plndatasize(DE)); //_memused[3] += plndncnum()*dof()*plndnenum()*dof()*sizeof(double);
	 DblNumMat denpos(dim(),splpos(DE).n());		  clear(denpos);	 iC( daxpy(R, splpos(DE), denpos) ); //scale
	 DblNumMat chkpos(dim(),splpos(DC).n());		  clear(chkpos);	 iC( daxpy(0.5*R, splpos(DC), chkpos) ); //scale
	 for(int i=0; i<dim(); i++) for(int j=0; j<splpos(DC).n(); j++) chkpos(i,j) = chkpos(i,j) + (double(idx(i))-0.5)*R;
	 
	 iC( _knl.kernel(denpos, denpos, chkpos, _DE2DCii) );
  }
  //---------matvec
  if(_hom==true) {
	 int sdof = this->sdof();
	 DblNumVec tmpden(sdof*splpos(DE).n(), false, _wsbuf);	 clear(tmpden);
	 vector<double> sclvec(sdof);	 for(int s=0; s<sdof; s++)		sclvec[s] = pow(2.0, l*_degvec[s]);
	 int cnt = 0;
	 for(int i=0; i<splpos(DE).n(); i++)
		for(int s=0; s<sdof; s++) {
		  tmpden(cnt) = den(cnt) * sclvec[s];		  cnt++;
		}
	 iC( dgemv(1.0, _DE2DCii, tmpden, 1.0, chk) );
  } else {
	 iC( dgemv(1.0, _DE2DCii, den, 1.0, chk) );
  }
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::plnden2effden(int l, const DblNumVec& plnden, DblNumVec& effden)
{
  DblNumVec regden(regpos().n()*sdof()); clear(regden);
  if(_hom==true) {
	 int sdof = this->sdof();
	 DblNumVec tmpden(sdof*splpos(UE).n(), false, _wsbuf);	 clear(tmpden);
	 vector<double> sclvec(sdof);	 for(int s=0; s<sdof; s++)		sclvec[s] = pow(2.0, l*_degvec[s]);
	 int cnt = 0;
	 for(int i=0; i<splpos(UE).n(); i++)
		for(int s=0; s<sdof; s++) {
		  tmpden(cnt) = plnden(cnt) * sclvec[s];		  cnt++;
		}
	 iC( splden2regden(tmpden, regden) );
  } else {
	 iC( splden2regden(plnden, regden) );
  }
  rfftwnd_real_to_complex(_forplan, sdof(), regden.data(), sdof(), 1, (fftw_complex*)(effden.data()), sdof(), 1);
  
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::splden2regden(const DblNumVec& splden, DblNumVec& regden)
{
  int np = _np;
  int rgnum = 2*np;
  int sdof = this->sdof();
  int count=0;
  //the order of iterating is the same as SampleGrid
  for(int i=0; i<np; i++)
	 for(int j=0; j<np; j++)
		for(int k=0; k<np; k++) {
		  if(i==0 || i==np-1 || j==0 || j==np-1 || k==0 || k==np-1) {
			 //the position is fortran style
			 int rgoff = (k+np/2)*rgnum*rgnum + (j+np/2)*rgnum + (i+np/2);
			 for(int f=0; f<sdof; f++) {
				regden(sdof*rgoff + f) += splden(sdof*count + f);
			 }
			 count++;
		  }
		}  //iC( PetscLogFlops(np*np*np*dof) );
  return 0;
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::effval2plnval(int level, const DblNumVec& effval, DblNumVec& plnval)
{
  DblNumVec regval(regpos().n()*tdof());
  rfftwnd_complex_to_real(_invplan, tdof(), (fftw_complex*)(effval.data()), tdof(), 1, regval.data(), tdof(), 1);
  iC( regval2splval(regval, plnval) );
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::regval2splval(const DblNumVec& regval, DblNumVec& splval)
{
  int np = _np;
  int rgnum = 2*np;
  int tdof = this->tdof();
  int count=0;
  //the order of iterating is the same as SampleGrid
  for(int i=0; i<np; i++)
	 for(int j=0; j<np; j++)
		for(int k=0; k<np; k++) {
		  if(i==0 || i==np-1 || j==0 || j==np-1 || k==0 || k==np-1) {
			 //the position is fortran style
			 int rgoff = (k+np/2)*rgnum*rgnum + (j+np/2)*rgnum + (i+np/2);
			 for(int f=0; f<tdof; f++) {
				splval(tdof*count + f) += regval(tdof*rgoff + f);
			 }
			 count++;
		  }
		}
  return 0;
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::UE2DC_dgemv(int l, Index3 idx, const DblNumVec& effden, DblNumVec& effval)
{
  OffTns<DblNumMat>& _UE2DC = (_hom==true) ? _ue2dc[0] : _ue2dc[l];
  double R = (_hom==true) ? 1.0 : 1.0/pow(2.0, l); //OffTns< DblNumMat >& _UE2DC = _ue2dc[l];  double R       = 1.0/pow(2.0, l);
  if(_UE2DC.m()==0)	 _UE2DC.resize(7,7,7,-3,-3,-3);
  DblNumMat& _UE2DCii = _UE2DC(idx[0], idx[1], idx[2]);
  
  int effnum = (2*_np+2)*(2*_np)*(2*_np);
  int sdof = this->sdof();
  int tdof = this->tdof();
  //---------compute matrix
  if(_UE2DCii.m()==0) { //compute it if necessary
	 //-----------------------	 //cerr<<"UE2DC(FFT) compute"<<endl;	 //COMPUTE FFT	 //Index3 idx = ii.idx();	 
	 iA( idx.linfty()>1 );
	 DblNumMat denpos(dim(),1);	 for(int i=0; i<dim(); i++)		denpos(i,0) = double(idx(i))*2.0*R; //shift
	 DblNumMat chkpos(dim(),regpos().n());	 clear(chkpos);	 iC( daxpy(R, regpos(), chkpos) );
	 DblNumMat tt(regpos().n()*tdof, sdof);
	 iC( _knl.kernel(denpos, denpos, chkpos, tt) );
	 // move data to tmp
	 DblNumMat tmp(tdof,regpos().n()*sdof);
	 for(int k=0; k<regpos().n();k++) {
		for(int i=0; i<tdof; i++)
		  for(int j=0; j<sdof; j++) {
			 tmp(i,j+k*sdof) = tt(i+k*tdof,j);
		  }
	 }
	 _UE2DCii.resize(tdof*sdof, effnum); //_memused[4] += dof*dof*effnum(UE)*sizeof(double);
	 //forward FFT from tmp to _UE2DCii;
	 rfftwnd_real_to_complex(_forplan, sdof*tdof, tmp.data(), sdof*tdof, 1, (fftw_complex*)(_UE2DCii.data()), sdof*tdof, 1);
  }
  //---------matvec
  /*
  //scale
  DblNumVec tmpden(effden.m(), false, _wsbuf);  clear(tmpden);
  fftw_complex* denptr = NULL;// = (fftw_complex*)(tmpden.data());
  if(_hom==true) {
	 vector<double> sclvec(sdof);	 for(int s=0; s<sdof; s++)		sclvec[s] = pow(2.0, l*_degvec[s]);
	 int cnt = 0;
	 for(int i=0; i<effden.m()/sdof; i++)
		for(int s=0; s<sdof; s++) {
		  tmpden(cnt) = effden(cnt) * sclvec[s];		cnt++;
		}
	 denptr = (fftw_complex*)(tmpden.data());
  } else {
	 denptr = (fftw_complex*)(effden.data());
  }
  */
  //fft mult
  double nrmfc = 1.0/double(regpos().n());
  fftw_complex* matptr = (fftw_complex*)(_UE2DCii.data());
  fftw_complex* denptr = (fftw_complex*)(effden.data());
  fftw_complex* chkptr = (fftw_complex*)(effval.data());
  int matstp = sdof*tdof;
  int denstp = sdof;
  int chkstp = tdof;
  
  double newalpha = nrmfc;
  for(int i=0; i<tdof; i++)
	 for(int j=0; j<sdof; j++) {
		int matoff = j*tdof + i;
		int denoff = j;
		int chkoff = i;
		iC( cptwvv(effnum/2, newalpha, matptr+matoff, matstp, denptr+denoff, denstp, chkptr+chkoff, chkstp) );
	 }
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::locpos(int tp, Point3 ctr, double rad, DblNumMat& pos)
{
  const DblNumMat& bas = splpos(tp);
  pos.resize(dim(), bas.n());
  for(int i=0; i<dim(); i++)
	 for(int j=0; j<pos.n(); j++)
		pos(i,j) = ctr(i) + rad * bas(i,j);
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::splposcal(int np, double R, DblNumMat& pos)
{
  int n = np*np*np - (np-2)*(np-2)*(np-2);
  pos.resize(dim(),n);
  double step = 2.0/(np-1);
  double init = -1.0;
  int count = 0;
  for(int i=0; i<np; i++)
	 for(int j=0; j<np; j++)
		for(int k=0; k<np; k++) {
		  if(i==0 || i==np-1 || j==0 || j==np-1 || k==0 || k==np-1) {
			 double x = init + i*step;
			 double y = init + j*step;
			 double z = init + k*step;
			 pos(0,count) = R*x;
			 pos(1,count) = R*y;
			 pos(2,count) = R*z;
			 count++;
		  }
		}
  iA(count==n);
  return 0;
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::regposcal(int np, double R, DblNumMat& pos)
{
  int n = 2*np*2*np*2*np;
  pos.resize(dim(), n);
  double step = 2.0/(np-1);
  int count = 0;
  for(int k=0; k<2*np; k++)
	 for(int j=0; j<2*np; j++)
		for(int i=0; i<2*np; i++) {
		  int gi = (i<np) ? i : i-2*np;
		  int gj = (j<np) ? j : j-2*np;
		  int gk = (k<np) ? k : k-2*np;
		  pos(0, count) = R * gi*step;
		  pos(1, count) = R * gj*step;
		  pos(2, count) = R * gk*step;
		  count ++;
		}
  iA(count==n);
  return 0;
}

// ---------------------------------------------------------------------- 
inline int MatMgnt3d::cptwvv(int n, double alpha, fftw_complex* x, int incx, fftw_complex* y, int incy, fftw_complex* z, int incz)
{
  for(int i=0; i<n; i++) {
	 (*z).re += alpha * ( (*x).re * (*y).re - (*x).im * (*y).im);
	 (*z).im += alpha * ( (*x).re * (*y).im + (*x).im * (*y).re);
	 x = x + incx;
	 y = y + incy;
	 z = z + incz;
  }  //iC( PetscLogFlops( 10*n ) );
  return 0;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
vector<MatMgnt3d> MatMgnt3d::_mmvec;

MatMgnt3d* MatMgnt3d::getmmptr(Kernel3d knl, int np)
{
  for(int i=0; i<_mmvec.size(); i++)
	 if(_mmvec[i].knl()==knl && _mmvec[i].np()==np)
		return &(_mmvec[i]);
  
  _mmvec.push_back( MatMgnt3d() );
  int last = _mmvec.size()-1;
  MatMgnt3d* tmp = &(_mmvec[last]); //get the last one
  tmp->knl() = knl;  tmp->np() = np;
  tmp->setup();
  return tmp;
}
