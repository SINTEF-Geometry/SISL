/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: s1920.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1920

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1920(SISLCurve *pc1,double edir[],int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1920(pc1,edir,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
     SISLCurve    *pc1;
     double   edir[];
     int      idim;
     double   aepsco;
     double   aepsge;
     int      *jpt;
     double   **gpar;
     int      *jcrv;
     SISLIntcurve ***wcurve;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find the extremal points/intervals of the curve pc1 in 
*              the direction edir.
*
*
*
* INPUT      : pc1    - Pointer to the curve.
*              edir   - The direction in which the extremal point(s)
*                       and/or interval(s) are to be calculated. If
*                       idim=1 a positive value indicates the maximum
*                       of the B-spline function and a negative value
*                       the minimum. If the dimension is greater that
*                       1 the array contains the coordinates of the
*                       direction vector.
*              idim   - Dimension of the space in which the vector edir
*                       lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single extremal points.
*              gpar   - Array containing the parameter values of the
*                       single extremal points in the parameter
*                       interval of the curve. The points lie continuous. 
*                       Extremal curves are stored in wcurve.
*              jcrv   - Number of extremal curves.
*              wcurve - Array containing descriptions of the extremal
*                       curves. The curves are only described by points
*                       in the parameter interval. The curve-pointers points
*                       to nothing. (See descrjption of Intcurve
*                       in intcurve.dcl).
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The scalar-product of the coefficients of the curve with
*              the vector edir is calculated giving a curve with dimension
*              1. Then the maxima of this curve is calculated.
*              
*
*
* REFERENCES :
*
*- 
* CALLS      : s1880 - Put extremal points/intervals into output format.
*              s1161 - Find maxima of one-dimensional curve.
*              s6scpr - Scalar-product between two vectors.
*              newCurve   - Create new curve.                        
*              newObject  - Create new object.
*              newIntpt   - Create new extremal point.
*              freeObject - Free the space occupied by a given object.
*              freeIntdat - Free space occupied by an extremal list.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-10.
* REVISED BY : Mike Floater, SI, 1991-04 for a rational curve.
*
*********************************************************************
*/                                                               
{                                                                     
  int ikind;               /* Type of curve pc1 is                     */
  int kstat = 0;           /* Local status variable.                     */
  int kpos = 0;            /* Position of error.                         */
  int kn;                  /* Number of vertices of curve.               */
  int kk;                  /* Order of curve.                            */
  double tmax;             /* Estimate of maximal value of 1-dim. curve. */
  double *st;              /* Pointer to knotvector of curve.            */
  double *scoef;           /* Pointer to vertices of curve.              */
  double *sc = NULL;       /* Pointer to vertices of curve in maxima
			      calculation.                               */
  double *spar = NULL;     /* Values of extremal points in the parameter 
			      area of the second object. Empty in this case.*/
  double *s1,*s2,*sstop;   /* Pointers used to traverse double-arrays.   */
  SISLIntdat *qintdat = NULL;  /* Maximum results.                     */
  SISLCurve *qc = NULL;        /* Pointer to curve in maxima calculation.    */
  SISLObject *qo1 = NULL;      /* Pointer to object in maxima calculation. */
  
  *jpt = 0;
  *jcrv = 0;
  
  
  /* Check dimension.  */
  
  if (pc1 -> idim != idim) goto err106;
  
  /* Describe curve with local variables.  */
  
  kn = pc1 -> in;
  kk = pc1 -> ik;
  st = pc1 -> et;
  ikind = pc1 -> ikind;
  
  if(ikind == 2 || ikind == 4)
  {
      scoef = pc1 -> rcoef;
      /* Allocate space for coeffecients of new curve.  */
      
      if ((sc = newarray(2*kn,double)) == NULL) goto err101; 
      
      /* Compute scalar-product of curve-vertices and direction vector. */
      /* Copy over weights. */
      
      for (s1=scoef,s2=sc,sstop=s2+2*kn; s2<sstop; s1+=idim+1,s2+=2)
      {
          *s2 = s6scpr(s1,edir,idim);
          *(s2+1) = *(s1+idim);
      }
  }
  else
  {
      scoef = pc1 -> ecoef;
      /* Allocate space for coeffecients of new curve.  */
      
      if ((sc = newarray(kn,double)) == NULL) goto err101; 
      
      /* Compute scalar-product of curve-vertices and direction vector. */
      
      for (s1=scoef,s2=sc,sstop=s2+kn; s2<sstop; s1+=idim,s2++)
      {
          *s2 = s6scpr(s1,edir,idim);
      }
  }

  
  /* Create new curve.  */
  
  qc = newCurve(kn,kk,st,sc,pc1->ikind,1,1);
  if (qc == NULL) goto err101;
  
  /* Create new object and connect curve to object.  */
  
  qo1 = newObject(SISLCURVE);
  if (qo1 == NULL) goto err101;
  qo1 -> c1 = qc;
  
  tmax = -(double)HUGE;
  
  /* Find maxima. */
  s1161(qo1,&tmax,aepsge,&qintdat,&kstat);
  if (kstat < 0) goto error;
  
  if (qintdat)
    {
      
      /* Express maximal points/intervals on output format.  */
      s1880(1,0,&qintdat->ipoint,qintdat->vpoint,&qintdat->ilist,qintdat->vlist
	    ,jpt,gpar,&spar,jcrv,wcurve,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* Extremal points/intervals found.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1920",*jstat,kpos);
  goto out;
  
  /* Dimensions conflicting.  */
  
 err106: *jstat = -106;
  s6err("s1920",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1920",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free allocated space.  */
  
  if (sc) freearray(sc);
  if (spar) freearray(spar);
  if (qintdat) freeIntdat(qintdat);
  if (qo1) freeObject(qo1);
  
  
  return;
}                                               

