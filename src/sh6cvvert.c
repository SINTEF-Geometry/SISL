#include "sisl-copyright.h"
#define SH6CVVERT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6cvvert(SISLCurve *pc1, SISLCurve *pc2, double *cpar1, double *cpar2)
#else
void sh6cvvert(pc1,pc2,cpar1,cpar2)
   SISLCurve *pc1;
   SISLCurve *pc2;
   double *cpar1;
   double *cpar2;
#endif
/*
*********************************************************************
* 
* PURPOSE    : Estimate the parameter values the closest vertices of
*              two B-spline curves
* 
* 
* 
* INPUT      : pc1      - Pointer to first curve.
*              pc1      - Pointer to second curve.
*
*
* OUTPUT     : cpar1    - Parameter value of closest vertex of the 1. curve.
*              cpar2    - Parameter value of closest vertex of the 2. curve.
* 
* 
* METHOD     : Find the vertices of the two curves that are closest 
*              to each other. Regard these vertices as Schoenberg
*              points, and compute the parameter values.
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      : s6dist - Distance between two points.           
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 09.00.
*
*********************************************************************
*/
{
  int ki,kj, kh;           /* Counters.   */
  int kdim = pc1->idim; /* Dimension of geometry space.      */
  int kminc1;           /* Number of closest vertex of 1. curve. */
  int kminc2;           /* Number of closest vertex of 2. curve. */
  int kn1 = pc1->in;    /* Number of coefficients of 1. curve. */
  int kn2 = pc2->in;    /* Number of coefficients of 2. curve. */
  int kk1 = pc1->ik;    /* Order of 1. curve. */
  int kk2 = pc2->ik;   /* Order of 2. curve. */
  double tdist;           /* Distance.   */
  double tmin = HUGE;     /* Minimum distance.  */
  double tpar;            /* Used to compute parameter values.   */
  double *s1,*s2;         /* Pointers into arrays.   */
  
  /* Find position of closest vertices. */
  
  for (s1=pc1->ecoef, ki=0; ki<kn1; s1+=kdim, ki++)
    for (s2=pc2->ecoef, kj=0; kj<kn2; s2+=kdim, kj++)
      {
	for (tdist=0.0, kh=kdim-1; kh>=0; kh--)
	  tdist += (s2[kh]-s1[kh])*(s2[kh]-s1[kh]);
	//	tdist = s6dist(s1,s2,kdim);
	if (tdist < tmin)
	  {
	    tmin = tdist;
	    kminc1 = ki;
	    kminc2 = kj;
	   }
	}
  
  /* Estimate parameter values of vertices.  */
  
  for (ki=kminc1+1, s1=pc1->et+ki, tpar=0.0; 
   ki<kminc1+kk1; tpar+=(*s1), s1++, ki++);
  *cpar1 = tpar/(double)(kk1-1);
  
  for (ki=kminc2+1, s1=pc2->et+ki, tpar=0.0; 
   ki<kminc2+kk2; tpar+=(*s1), s1++, ki++);
  *cpar2 = tpar/(double)(kk2-1);

  return;
}
