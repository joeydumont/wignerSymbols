#ifndef FORTRANLINKAGE_H
#define FORTRANLINKAGE_H

extern "C"
{
	extern void drc3jj_wrap(double,double,double,double,double*,double*,double*,int,int*);
	extern void drc6j_wrap(double,double,double,double,double*,double*,double*,int,int*);
}

#endif // FORTRANLINKAGE_H