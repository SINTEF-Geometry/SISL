
/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

#define SH1831

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     sh1831(SISLCurve *pc1, SISLCurve *pc2, int isign, double epoint[], 
	    double enorm[], double aepsge, int *jstat)
#else
void sh1831(pc1, pc2, isign, epoint, enorm, aepsge, jstat)
     SISLCurve *pc1;
     SISLCurve *pc2;
     int isign;
     double epoint[];
     double enorm[];
     double aepsge;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Interception test between two curves. Test if it is
*              possible to split the curves by a given plane.
*
*
*
* INPUT      : pc1    - First curve.
*              pc2    - Second curve.
*              epoint - Point in plane.
*              enorm  - Plane normal.
*              aepsge - Geometric tolerance.
*
*
*
* OUTPUT     : jstat  - status messages  
*                              = 1      : Possibility of inner intersections.
*                              = 0      : No possibility of inner intersections.
*                              < 0      : error
*
*
* METHOD     : 
*
*
* REFERENCES :
*
*-
* CALLS      : s6scpr - Scalar product between two vectors.
*              s6diff - Difference vector between two vectors.
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 94-02.
*
*********************************************************************
*/
{
  int kpos = 0;          /* Position of error.               */
  int ki;
  int kdim;              /* Dimension of space.              */
  int kbez1, kbez2;      /* Indicates if the curves are of type Bezier. */
  int ksignprev = 0;     /* Sign of distance between previous curve and plane.*/
  int ksign1 = 0;        /* Sign of distance between curve and plane.*/
  int ksign2 = 0;        /* Sign of distance between curve and plane.*/
  double tdist;          /* Distance between coefficient and plane.     */
  double *s1;            /* Pointer to coefficient of curve. */
  double sdiff[3];       /* Difference vector.               */
  
  /* Test dimension of geometry space. */
  
  kdim = pc1->idim;
  if (kdim != 2 && kdim != 3) goto err105;
  if (kdim != pc2->idim) goto err106;
  
  /* Test if the curves are Bezier curves. */
  
  kbez1 = (pc1->ik == pc1->in) ? 1 : 0;
  kbez2 = (pc2->ik == pc2->in) ? 1 : 0;
  
  /* For each curve, compute the distance between the coefficients of the
     curve and the given plane.    */
  
  for (s1=pc1->ecoef, ki=0; ki<pc1->in; ki++, s1+=kdim)
  {
     s6diff(epoint, s1, kdim, sdiff);
     tdist = s6scpr(sdiff, enorm, kdim);
     
     if (fabs(tdist) <= aepsge && !kbez1 && !(ki==0 || ki==pc1->in-1)) break;
     ksign2 = (DEQUAL(tdist,DZERO)) ? 0 : ((tdist > 0) ? 1 : -1);
     if (ksign1*ksign2 < 0) break;
     ksign1 = ksign2;
  }
  if (ki < pc1->in)
  {
     *jstat = 1;
     goto out;
  }

  ksignprev = isign*ksign1;
  ksign1 = 0;
  for (s1=pc2->ecoef, ki=0; ki<pc2->in; ki++, s1+=kdim)
  {
     s6diff(epoint, s1, kdim, sdiff);
     tdist = s6scpr(sdiff, enorm, kdim);
     
     if (fabs(tdist) <= aepsge && !kbez2 && !(ki==0 || ki==pc2->in-1)) break;
     ksign2 = (DEQUAL(tdist,DZERO)) ? 0 : ((tdist > 0) ? 1 : -1);
     if (ksign1*ksign2 < 0) break;
     if (ksignprev*ksign1 > 0) break;
     ksign1 = ksign2;
  }
  if (ki < pc2->in)
  {
     *jstat = 1;
     goto out;
  }
  
  goto out;
  
  /* Error in input. Dimension not equal to 2 or 3.  */
  
 err105: *jstat = -105;
  s6err("sh1831",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimensions conflicting.  */
  
 err106: *jstat = -106;
  s6err("sh1831",*jstat,kpos);
  goto out;
  
 out:
  
  return;
}
