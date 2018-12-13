/**
 * A best exemplar finder.  Scans over the entire image (using a
 * sliding window) and finds the exemplar which minimizes the sum
 * squared error (SSE) over the to-be-filled pixels in the target
 * patch. 
 *
 * @author Sooraj Bhat
 */
#include "mex.h"
#include <limits.h>

// (m,n) : size of Ip; (mm,nn):size of the entire image; img, Ip, toFill: mask for Ip, to be filled; sourceRegion for 
void bestexemplarhelper(const int mm, const int nn, const int m, const int n, 
			const double *img, const double *Ip, 
                        const double *depth, const double *Dp,
			const mxLogical *toFill, const mxLogical *sourceRegion,
			double *best) 
{
  register int i,j,ii,jj,ii2,jj2,M,N,I,J,ndx,ndx2,mn=m*n,mmnn=mm*nn;
  double patchErr=0.0,err=0.0,bestErr=1000000000.0;
  const double alpha = 1.0;
  /* foreach patch */
  N=nn-n+1;  M=mm-m+1;    // Search window boundary.
  for (j=1; j<=N; ++j) {  // For each column..
    J=j+n-1;
    for (i=1; i<=M; ++i) {
      I=i+m-1;
      /*** Calculate patch error ***/
      /* foreach pixel in the current patch */
      for (jj=j,jj2=1; jj<=J; ++jj,++jj2) {
	for (ii=i,ii2=1; ii<=I; ++ii,++ii2) {
	  ndx=ii-1+mm*(jj-1);                        // Notice the array arrangemeent is Fortran style.
	  if (!sourceRegion[ndx])
	    goto skipPatch;
	  ndx2=ii2-1+m*(jj2-1);
	  if (!toFill[ndx2]) {
	    err=img[ndx      ] - Ip[ndx2    ]; patchErr += err*err;
	    err=depth[ndx      ] - depth[ndx2    ]; patchErr += alpha*err*err;
	    err=img[ndx+=mmnn] - Ip[ndx2+=mn]; patchErr += err*err;
	    err=depth[ndx] - depth[ndx2]; patchErr += alpha*err*err;
	    err=img[ndx+=mmnn] - Ip[ndx2+=mn]; patchErr += err*err;
	  }
	}
      }
      /*** Update ***/
      if (patchErr < bestErr) {
	bestErr = patchErr; 
	best[0] = i; best[1] = I;
	best[2] = j; best[3] = J;
      }
      /*** Reset ***/
    skipPatch:
      patchErr = 0.0; 
    }
  }
}

/* best = bestexemplarhelper(mm,nn,m,n,img,Ip,toFill,sourceRegion); */
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{
  int mm,nn,m,n;
  double *img,*Ip,*best;
  double *depth, *Dp;
  mxLogical *toFill,*sourceRegion;

  /* Extract the inputs */
  mm = (int)mxGetScalar(prhs[0]);    // (mm, nn) size of image.
  nn = (int)mxGetScalar(prhs[1]);    // 
  m  = (int)mxGetScalar(prhs[2]);    // (m, n) size of patch.
  n  = (int)mxGetScalar(prhs[3]);    //
  img = mxGetPr(prhs[4]);            // Image data.
  Ip  = mxGetPr(prhs[5]);            // Image patch.
  depth = mxGetPr(prhs[6]);          // Depth Gradient data.
  Dp  = mxGetPr(prhs[7]);            // Depth Gradient patch.
  toFill = mxGetLogicals(prhs[8]);
  sourceRegion = mxGetLogicals(prhs[9]);
  
  /* Setup the output */
  plhs[0] = mxCreateDoubleMatrix(4,1,mxREAL);
  best = mxGetPr(plhs[0]);
  best[0]=best[1]=best[2]=best[3]=0.0;

  /* Do the actual work */
  bestexemplarhelper(mm,nn,m,n,img,Ip,depth, Dp,toFill,sourceRegion,best);
}
