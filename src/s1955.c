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
 * $Id: s1955.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1955

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1955(SISLCurve *pc1,SISLCurve *pc2,double aepsco,double aepsge,
	   int *jpt,double **gpar1,double **gpar2,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1955(pc1,pc2,aepsco,aepsge,jpt,gpar1,gpar2,jcrv,wcurve,jstat)
     SISLCurve    *pc1;
     SISLCurve    *pc2;
     double   aepsco;
     double   aepsge;
     int      *jpt;
     double   **gpar1;
     double   **gpar2;
     int      *jcrv;
     SISLIntcurve ***wcurve;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find the closest point between the two curve pc1 and pc2.
*
*
*
* INPUT      : pc1    - Pointer to the first curve in the closest point
*                       problem.
*              pc2    - Pointer to the second curve in the closest point
*                       problem.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single closest points.
*              gpar1  - Array containing the parameter values of the
*                       single closest points in the parameter intervals
*                       of the first curve. The points lie continuous. 
*                       Closest curves are stored in wcurve.
*              gpar2  - Array containing the parameter values of the
*                       single closest points in the parameter intervals
*                       of the second curve. The points lie continuous. 
*              jcrv   - Number of closest curves.
*              wcurve - Array containing descriptions of the closest
*                       curves. The curves are only described by points
*                       in the parameter area. The curve-pointers points
*                       to nothing. (See descrjption of Intcurve
*                       in intcurve.dcl).
*                       If the curves given as input are degnenerate an
*                       intersection point can be returned as an intersection
*                       curve. Use s1327 to decide if an intersection curve
*                       is a point on one of the curves.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Generate the difference surface between the two curves,
*              and find those points of this surface closest to origo.
*              
*
*
* REFERENCES :
*
*- 
* CALLS      : s1954    - Find closest points between a surface and a point.
*              s1956    - Create difference surface between two curves.
*              freeSurf - Free space occupied by a surface.
*              newIntcurve - Create new closest curve.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  int kstat = 0;               /* Local status variable.                 */
  int kpos = 0;                /* Position of error.                     */
  int ki,kj;                   /* Counters.                              */
  int kdim;                    /* Dimension of the space in which the 
				  curves lie.                            */
  double *sorigo = NULL;       /* Array representing origo.              */
  double *spar1 = NULL;        /* Parameter-values in the parameter
				  interval of the first curve of the
				  closest curve.                         */
  double *spar2 = NULL;        /* Parameter-values in the parameter
				  interval of the second curve of the
				  closest curve.                         */
  double *ssing = NULL;        /* Single closest points between 
				  difference surface and origo.          */
  SISLSurf *qsdiff = NULL;         /* Difference surface.                    */
  SISLIntcurve *qcurve;            /* Pointer to closest curve.              */
  int k1=1,k2=2,k4=4;          /* Constants                              */
  
  /* Test input.  */
  
  kdim = pc1 -> idim;
  if (kdim != pc2 -> idim) goto err106;
  
  /* Allocate space for the point origo.  */
  
  if ((sorigo = new0array(kdim,double)) == NULL) goto err101;
  
  /* Generate the difference surface between the curves.  */
  
  s1956(pc1,pc2,&qsdiff,&kstat);
  if (kstat < 0) goto error;
  
  if (kstat > 0)
    {
      
      /* The curves define a closest interval. Allocate space
	 for points defining the interval.                      */
      
      if ((spar1 = newarray(k2,double)) == NULL) goto err101;
      if ((spar2 = newarray(k2,double)) == NULL) goto err101;
      
      if (kstat == 1)
	{
	  
	  /* The curves are equal.  */
	  
	  spar1[0] = *(pc1->et + pc1->ik - 1);
	  spar1[1] = *(pc1->et + pc1->in);
	  spar2[0] = *(pc2->et + pc2->ik - 1);
	  spar2[1] = *(pc2->et + pc2->in);
	  
	}
      else if (kstat == 2)
	{
	  
	  /* The curves are equal, but oppositely directed.  */
	  
	  spar1[0] = *(pc1->et + pc1->ik - 1);
	  spar1[1] = *(pc1->et + pc1->in);
	  spar2[0] = *(pc2->et + pc2->in);
	  spar2[1] = *(pc2->et + pc2->ik - 1);
	}
      
      *wcurve = NULL;
      if ((*wcurve = newarray(k1,SISLIntcurve*)) == NULL) goto err101;
      *jcrv = 1;
      
      **wcurve = NULL;
      if ((**wcurve = newIntcurve(k2,k1,k1,spar1,spar2,k4)) == NULL) goto err101;
      
    }
  else
    {
      
      /* Find the closest points between the difference function and origo.*/
      
      s1954(qsdiff,sorigo,kdim,aepsco,aepsge,jpt,&ssing,jcrv,wcurve,&kstat);
      if (kstat < 0) goto error;
      
      if (*jpt > 0)
	{
	  
	  /* Allocate space for output of single closest points.  */
	  
	  *gpar1 = *gpar2 = NULL;
	  if ((*gpar1 = newarray(*jpt,double)) == NULL) goto err101;
	  if ((*gpar2 = newarray(*jpt,double)) == NULL) goto err101;
	  
	  /* Copy single closest points to output arrays.  */
	  
	  for (ki=0; ki<(*jpt); ki++)
	    {
	      *(*gpar1 + ki) = ssing[2*ki];
	      *(*gpar2 + ki) = ssing[2*ki+1];
	    }
	}
      
      for (kj=0; kj<(*jcrv); kj++)
	{
	  
	  /* Correct parameter-arrays of closest curve.  */
	  
	  qcurve = *(*wcurve + kj);
	  qcurve->ipar1 = qcurve->ipar2 = 1;
	  
	  /* Allocate space for new parameter arrays.  */
	  
	  spar1 = spar2 = NULL;
	  if ((spar1 = newarray(qcurve->ipoint,double)) == NULL) goto err101;
	  if ((spar2 = newarray(qcurve->ipoint,double)) == NULL) goto err101;
	  
	  /* Copy parameter-values to new arrays.  */
	  
	  for (ki=0; ki<qcurve->ipoint; ki++)
	    {
	      spar1[ki] = *(qcurve->epar1+2*ki);
	      spar2[ki] = *(qcurve->epar1+2*ki+1);
	    }
	  
	  freearray(qcurve->epar1);
	  qcurve -> epar1 = spar1;
	  qcurve -> epar2 = spar2;
	}
    }
  
  /* Closest points found.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1955",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimensions conflicing.  */
  
 err106: *jstat = -106;
  s6err("s1955",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1955",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free space occupied local variables.  */
  
  if (ssing != NULL) freearray(ssing);
  if (sorigo != NULL) freearray(sorigo);
  if (qsdiff != NULL) freeSurf(qsdiff);
  
  return;
}
