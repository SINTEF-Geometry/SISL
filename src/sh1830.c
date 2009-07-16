/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: sh1830.c,v 1.2 2001-03-19 15:59:05 afr Exp $
 *
 */


#define SH1830

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
  sh1830(SISLObject *po1,SISLObject *po2,double aepsge,int *jstat)
#else
void sh1830(po1,po2,aepsge,jstat)
     SISLObject  *po1;
     SISLObject  *po2;
     double aepsge;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Perform improved box-test in curve/surface intersection.
*              The test is based on the main tangent of the curve and
*              the main normal of the surface.
*
*
*
* INPUT      : po1    - Pointer to object involved in intersection problem.
*              po2    - Pointer to object involved in intersection problem.
*
*
*
* OUTPUT     : jstat  - status messages  
*                              = 2      : Only edge-intersections possible.
*                              = 1      : Possibility of intersection.
*                              = 0      : No possibility of intersection.
*                              < 0      : error
*
*
* METHOD     : Let the main tangent of the curve be the vector from the
*              first to the last vertex, and let the main normal of the
*              surface be the cross-product of the main diagonals. Perform
*              rotated box-tests with respect to this two vectors.
*
*
* REFERENCES :
*
*-
* CALLS      : sh1834  - Perform rotated box-test.
*              s6diff - Difference vector between two vectors.
*              s6crss - Cross-product of two vectors.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* REWISED BY : Vibeke Skytt, SI, 91-01.
*
*********************************************************************
*/
{
  int kstat = 0;         /* Local status variable.           */
  int kpos = 0;          /* Position of error.               */
  int kdim;              /* Dimension of space.              */
  int knc;               /* Number of vertices of curve.     */
  int kn1,kn2;           /* Number of vertices of surface.   */
  double *scurve;        /* Vertices of curve.               */
  double *ssurf;         /* Vertices of surface.             */
  double *stan = SISL_NULL;   /* Main tangent of curve.           */
  double *sdiag1 = SISL_NULL; /* First main diagonal of surface.  */
  double *sdiag2 = SISL_NULL; /* Second main diagonal of surface. */
  double *snorm = SISL_NULL;  /* Main normal of surface.          */
  SISLCurve *qcurve;     /* Pointer to curve.                */
  SISLSurf *qsurf;       /* Pointer to surface.              */
  
  /* Test input.  */
  
  if (!((po1->iobj == SISLSURFACE && po2->iobj == SISLCURVE) ||
	(po1->iobj == SISLCURVE && po2->iobj == SISLSURFACE))) 
     goto err121;
  
  /* Set pointers to objects.  */
  
  if (po1->iobj == SISLSURFACE)
  {
     qsurf = po1->s1;  qcurve = po2->c1;
  }
  else
  {
     qsurf = po2->s1;  qcurve = po1->c1;
  }
  
  /* Test dimension.  */
  
  kdim = qsurf -> idim;
  if (kdim != 3) goto err104;
  if (kdim != qcurve -> idim) goto err106;
  
  /* Allocate space for local arrays.  */
  
  if ((stan = newarray(kdim,double)) == SISL_NULL) goto err101;
  if ((sdiag1 = newarray(kdim,double)) == SISL_NULL) goto err101;
  if ((sdiag2 = newarray(kdim,double)) == SISL_NULL) goto err101;
  if ((snorm = newarray(kdim,double)) == SISL_NULL) goto err101;
  
  /* Describe curve with local parameters.  */
  
  knc = qcurve->in;
  scurve = qcurve->ecoef;
  
  /* Describe surface with local parameters.  */
  
  kn1 = qsurf->in1;
  kn2 = qsurf->in2;
  ssurf = qsurf->ecoef;
  
  /* Fetch main tangent of curve.  */
  
  s6diff(scurve+(knc-1)*kdim,scurve,kdim,stan);
  
  /* Fetch main diagonals of surface.  */
  
  s6diff(ssurf+(kn1*kn2-1)*kdim,ssurf,kdim,sdiag1);
  
  s6diff(ssurf+kn1*(kn2-1)*kdim,ssurf+(kn1-1)*kdim,kdim,sdiag2);
  
  /* Compute main normal of surface.  */
  
  s6crss(sdiag1,sdiag2,snorm);
  
  /* Perform rotated box-test.  */
  
  sh1834(po1,po2,aepsge,kdim,stan,snorm,&kstat);
  if (kstat < 0) goto error;
  
  if (kstat == 1)
    {
      kstat = 0;
      
      sh1834(po1,po2,aepsge,kdim,snorm,stan,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* Improved box-test performed.  */
  
  *jstat = kstat;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("sh1830",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 3.  */
  
 err104: *jstat = -104;
  s6err("sh1830",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimensions conflicting.  */
  
 err106: *jstat = -106;
  s6err("sh1830",*jstat,kpos);
  goto out;
  
  /* Error in kind of object.  */
  
  err121: *jstat = -121;
  s6err("s1930",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("sh1830",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free space occupied by local arrays.  */
  
  if (stan != SISL_NULL) freearray(stan);
  if (sdiag1 != SISL_NULL) freearray(sdiag1);
  if (sdiag2 != SISL_NULL) freearray(sdiag2);
  if (snorm != SISL_NULL) freearray(snorm);
  
  return;
}
