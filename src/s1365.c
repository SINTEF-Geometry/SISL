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
 * $Id: s1365.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1365

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1365(SISLSurf *ps, double aoffset, double aepsge, double amax,
	   int idim, SISLSurf **rs, int *jstat)
#else
void s1365(ps,aoffset,aepsge,amax,idim,rs,jstat)
     SISLSurf   *ps;
     double aoffset;
     double aepsge;
     double amax;
     int    idim;
     SISLSurf   **rs;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To create a B-spline approximating the offset surface of
*              a B-spline surface.
*
*
*
* INPUT      : ps     - The input B-spline surface.
*              aoffset- The offset distance.
*                       If idim=2 a positive signe on this value put the
*                       offset on the side of the positive normal vector,
*                       and a negative sign puts the offset on the sign
*                       of the negative normal vector.
*                       If idim=3 the offset is determined by the cross
*                       product of the tangent vector and the anorm vector.
*                       The offset distance is multiplied by this vector.
*              aepsge - Maximal deviation allowed between true offset surface
*                       and the approximated offset surface.
*              amax   - Maximal stepping length. Is negleceted if amax<=aepsge
*                       If amax==0 then a maximal step length of the longest
*                       box side is used.
*              idim   - The dimension of the space (2 or 3).
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              rs     - Pointer the approximated offset surface
*
* METHOD     : The 4 edge curves of the surface are extracted. Offset curves
*              of these 4 edge curves are approximated and a common
*              basis for the two pairs of opposite offset curves is calculated.
*              Vertices are recomputed. Finally data-reduction is performed.
*
* EXAMPLE OF USE:
*
*
* REFERENCES :
*
*
*-
* CALLS      : s1435     - Pick a given edge-curve of a B-spline surface.
*              s1360     - Approximate the offset curve of a B-spline curve.
*              s1933     - Put a set of curves on a common basis.
*              s1366     - Create a B-spline surface approximating the offset
*                          surface of a B-spline surface
*              s1345     - Data reduction of a B-spline surface. 
*
* WRITTEN BY : Per Evensen,  SI, 89-3.
* REWISED BY : Per Evensen,  SI, 90-9; Corrected start/end parameter values of
*                                      common curves.
* REWISED BY : Vibeke Skytt, SINTEF, 9403. Introduced data reduction and freed
*                                          scratch used by local curves.
*
*********************************************************************
*/
{
  SISLCurve *pc[4];
  SISLCurve *rc[4],*rc13[2],*rc24[2];
  int kstat = 0;     /* Local status variable.                           */
  int kpos = 0;      /* Position of error.                               */
  int kdim;          /* Dimension of the space in which the surface lies.*/
  int knbcrv = 2;    /* Number of curves in set. */
  int kopen = 1;     /* Flag telling that the resulting surface should
			be open in both parameter directions. */
  int kn13,kord13;   /* Number of vertices and order of edge curves along  
			1. parameter direction. */
  int kn24,kord24;   /* Number of vertices and order of edge curves along  
			2. parameter direction. */
  int nder[4];       /* Number of edges along each surface edge.         */
  double sp[4];      /* Parameter value of edge in constatnt direction.  */
  double toffset = (double)0.0; /* Local offset value for extraction of edge-
				   curves. */
  double snorm[3];   /* Local normal vector for extraction of edge-curves. */
  double tstart1,tend1; /* Endpoints of parameter interval in first 
			   direction.                                     */
  double tstart2,tend2; /* Endpoints of parameter interval in second 
			   direction.                                     */
  double *sknot13=NULL;/* Pointer to common knot-vector of edge curves along
			  1. parameter direction. */
  double *scoef13=NULL;/* Pointer to vertices expressed in same basis of edge
			  curves along 1. parameter direction. */
  double *sknot24=NULL;/* Pointer to common knot-vector of edge curves along
			  2. parameter direction. */
  double *scoef24=NULL;/* Pointer to vertices expressed in same basis of edge
			  curves along 2. parameter direction. */
  int kj,kk;            /* Loop controller. */
  double tepsge = (double)0.5*aepsge;  /* Local tolerance.             */
  SISLSurf *qs = NULL;                 /* Intermediate offset surface. */
  double seps[3];
  double sedgeeps[12];
  int ledgefix[4];
  int kopt;
  int ktmax;
  double smaxerr[3];
  
  /* Initialization of variables */
  kdim = ps -> idim;
  for (kk=0; kk<4; kk++) 
  {
     nder[kk] = 1;
     pc[kk] = NULL;
     rc[kk] = NULL;
  }
  for (kk=0; kk<3; kk++) snorm[kk] = DNULL;
  
  /* Set parameters for the data reduction. */
  
  for (kk=0; kk<idim; kk++)
     seps[kk] = pow(tepsge, (double)1/(double)idim);
  for (kj=0; kj<4; kj++)
  {
     ledgefix[kj] = 0;
     for (kk=0; kk<kdim; kk++) sedgeeps[kj*4+kk] = DNULL;
  }
  kopt = 3;
  ktmax = 20;
  
  /* Fetch the 4 edge-curves of surface */
  
  for (kk=0; kk<4; kk++)
    {
      s1435(ps,kk,&pc[kk],&sp[kk],&kstat);
      if (kstat < 0) goto error;
    }
  
  /* Create a B-spline curve approximating the offset curve of the 4 edges */
  
  for (kk=0; kk<4; kk++)
    {
      s1360(pc[kk],toffset,tepsge,snorm,amax,kdim,&rc[kk],&kstat);
      if (kstat<0) goto error;
    }
  
  /* Rearrange the pointers to the 4 edge curves. */
  
  rc13[0] = rc[0];
  rc13[1] = rc[2];
  rc24[0] = rc[1];
  rc24[1] = rc[3];
  
  /* Fetch endpoints of parameter intervals.  */
  
  tstart1 = *(ps->et1 + ps->ik1 - 1);
  tend1 = *(ps->et1 + ps->in1);
  tstart2 = *(ps->et2 + ps->ik2 - 1);
  tend2 = *(ps->et2 + ps->in2);
  
  /* Put the edge curves along 1. parameter direction into common basis. */

  s1933(knbcrv,rc13,tstart1,tend1,&sknot13,&kn13,&kord13,&kstat);
  if (kstat<0) goto error;
  
  /* Put the edge curves along 2. parameter direction into common basis. */
  
  s1933(knbcrv,rc24,tstart2,tend2,&sknot24,&kn24,&kord24,&kstat);
  if (kstat<0) goto error;
  
  /* Create a B-spline surface approximating the offset surface of a B-spline
     surface. */
  
  s1366(ps,aoffset,tepsge,amax,idim,sknot13,kn13,kord13,
	sknot24,kn24,kord24,&qs,&kstat);
  if (kstat<0) goto error;
  
  /* Try to reduce the amount of data necessary to represent the
     surface. */
  
  s1345(qs, seps, ledgefix, sedgeeps, REL_COMP_RES, kopt, ktmax, rs, 
	smaxerr, &kstat);
  if (kstat < 0) goto error;
  
  /* Surface approximated. */
  
  *jstat = 0;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1365",*jstat,kpos);
  goto out;
  
  out: 
     if (qs != NULL) freeSurf(qs);
     for (kk=0; kk<4; kk++)
     {
	if (pc[kk] != NULL) freeCurve(pc[kk]);
	if (rc[kk] != NULL) freeCurve(rc[kk]);
     }
     
     return;
}
