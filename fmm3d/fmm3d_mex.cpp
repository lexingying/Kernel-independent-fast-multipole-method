#include "mex.h"
#include "matrix.h"

#include "fmm3d.hpp"

#include "mexaux.hpp"

using namespace std;

//digital curvelet transform
extern void _main();

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if(nrhs!=8)
	 mexErrMsgTxt("8 inputs required");
  if(nlhs!=1)
	 mexErrMsgTxt("1 outputs required");
  
  DblNumMat srcpos; mex2cpp(prhs[0], srcpos);
  DblNumMat srcnor; mex2cpp(prhs[1], srcnor);
  DblNumMat trgpos; mex2cpp(prhs[2], trgpos);
  
  Point3 ctr; mex2cpp(prhs[3], ctr); //cerr<<ctr<<endl;
  int rootlvl; mex2cpp(prhs[4], rootlvl); //cerr<<rootlvl<<endl;
  Kernel3d knl; mex2cpp(prhs[5], knl); //cerr<<knl.kt()<<endl;
  int np; mex2cpp(prhs[6], np); iA(np==4 || np==6 || np==8); //cerr<<np<<endl;
  
  DblNumVec srcden; mex2cpp(prhs[7], srcden);
  
  char tmp[100];
  map<string, string> opts;
  if(       np==4) {
	 sprintf(tmp, "%d", np); opts["-fmm3d_np"] = string(tmp);
	 sprintf(tmp, "%d", 150);opts["-fmm3d_let3d_ptsmax"] = string(tmp);
	 sprintf(tmp, "%d", 10); opts["-fmm3d_let3d_maxlevel"] = string(tmp);
  } else if(np==6) {
	 sprintf(tmp, "%d", np); opts["-fmm3d_np"] = string(tmp);
	 sprintf(tmp, "%d", 150);opts["-fmm3d_let3d_ptsmax"] = string(tmp);
	 sprintf(tmp, "%d", 10); opts["-fmm3d_let3d_maxlevel"] = string(tmp);
  } else if(np==8) {
	 sprintf(tmp, "%d", np); opts["-fmm3d_np"] = string(tmp);
	 sprintf(tmp, "%d", 150);opts["-fmm3d_let3d_ptsmax"] = string(tmp);
	 sprintf(tmp, "%d", 10); opts["-fmm3d_let3d_maxlevel"] = string(tmp);
  }
  
  vector<MatMgnt3d> mmvec;
  FMM3d* fmm = new FMM3d("fmm3d_");
  fmm->srcpos()=&srcpos;  fmm->srcnor()=&srcpos;
  fmm->trgpos()=&trgpos;
  fmm->ctr() = Point3(0,0,0);  fmm->rootlvl() = 0;
  fmm->knl() = knl;
  fmm->mmvptr() = &mmvec;
  iC( fmm->setup(opts) );
  
  iA( srcden.m()==knl.sdof()*srcpos.n()); //check size
  DblNumVec trgval(knl.tdof()*trgpos.n());
  iC( fmm->eval(srcden, trgval) );
  iC( fmm->check(srcden, trgval, 20) );
  
  cpp2mex(trgval, plhs[0]);
  
  delete fmm;
  
  return;
}
