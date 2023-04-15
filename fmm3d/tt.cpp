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

#include "fmm3d.hpp"

using namespace std;

int optionsCreate(const char* optfile, map<string,string>& options)
{
  options.clear();
  ifstream fin(optfile);
  if(fin.good()==false) {
	 cerr<<"wrong option file"<<endl;	 exit(1);
  }
  string name;  fin>>name;
  while(fin.good()) {
	 char cont[100];	 fin.getline(cont, 99);
	 options[name] = string(cont);
	 fin>>name;
  }
  fin.close();
  return 0;
}

int main(int argc, char** argv)
{
  srand48( (long)time(NULL) );
  iA(argc==2);
  map<string,string> opts;
  optionsCreate(argv[1], opts);
  
  //1. allocate random data
  map<string,string>::iterator mi;
  mi = opts.find("-numsrc"); assert(mi!=opts.end());
  int numsrc;  { istringstream ss((*mi).second);  ss>>numsrc; }
  mi = opts.find("-numtrg"); assert(mi!=opts.end());
  int numtrg;  { istringstream ss((*mi).second);  ss>>numtrg; }
  mi = opts.find("-kt"); assert(mi!=opts.end());
  int kt;  { istringstream ss((*mi).second);  ss>>kt; } //cerr<<kt<<endl;
  
  vector<double> tmp(2);	 tmp[0] = 1;	 tmp[1] = 0.25; //coefs in the kernel
  Kernel3d knl(kt, tmp);
  
  DblNumMat srcpos(3, numsrc);
  for(int i=0; i<numsrc; i++) {
	 srcpos(0,i) = (2.0*drand48()-1.0);
	 srcpos(1,i) = (2.0*drand48()-1.0);
	 srcpos(2,i) = (2.0*drand48()-1.0);
  }
  DblNumMat trgpos(3, numtrg);
  for(int i=0; i<numtrg; i++) {
	 trgpos(0,i) = (2.0*drand48()-1.0);
	 trgpos(1,i) = (2.0*drand48()-1.0);
	 trgpos(2,i) = (2.0*drand48()-1.0);
  }
  
  int sdof = knl.sdof();
  int tdof = knl.tdof();
  DblNumVec srcden(sdof * numsrc);
  for(int i=0; i<numsrc; i++) {
	 for(int d=0; d<sdof; d++)
		srcden(d + i*sdof) = drand48(); //(2.0*drand48()-1.0);
  }
  DblNumVec trgval(tdof * numtrg);
  
  //2. allocate fmm 
  clock_t ck0, ck1;
  vector<MatMgnt3d> mmvec;
  
  FMM3d* fmm = new FMM3d("fmm3d_");
  fmm->srcpos()=&srcpos;  fmm->srcnor()=&srcpos;
  fmm->trgpos()=&trgpos;
  fmm->ctr() = Point3(0,0,0); // CENTER OF THE TOPLEVEL BOX
  fmm->rootlvl() = 0;         // 2^(-rootlvl) is the RADIUS OF THE TOPLEVEL BOX
  fmm->knl() = knl;
  
  ck0 = clock();
  iC( fmm->setup(opts) );
  ck1 = clock();  cout<<"fmm setup used "<<double(ck1-ck0)/CLOCKS_PER_SEC<<"secs "<<endl;
  
  //3. run fmm
  for(int i=0; i<3; i++) {
	 ck0 = clock();
	 iC( fmm->eval(srcden, trgval) );
	 ck1 = clock();  cout<<"fmm eval used "<<double(ck1-ck0)/CLOCKS_PER_SEC<<"secs "<<endl;
  }
  
  //4. check
  ck0 = clock();
  iC( fmm->check(srcden, trgval, 20) );
  ck1 = clock();  cout<<"fmm check used "<<double(ck1-ck0)/CLOCKS_PER_SEC<<"sec "<<endl;
  
  delete fmm;
  
  return 0;
}
