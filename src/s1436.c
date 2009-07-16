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
 * $Id: s1436.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1436

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1436(SISLSurf *ps1,double apar,SISLCurve **rcurve,int *jstat)
#else
void s1436(ps1,apar,rcurve,jstat)
     SISLSurf   *ps1;
     double apar;
     SISLCurve  **rcurve;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make constant parameter curve in the surface. The 
*              constant parameter value used is apar and is in the 
*              second parameter direction.
*
*
*
* INPUT      : ps1    - Surface.
*              apar   - Parameter value to use whe picking out constant
*                       parameter curve in second parameter direction.
*
*
*
* OUTPUT     : rcurve - Constant parameter curve.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The surface is treated as a (ps1->idim*ps1->in1)-
*              dimensional curve with ps1->in2 vertices. The value
*              of this curve at apar is calculated. This value is
*              then viewed as the vertices of a (ps1->idim)-dimensional
*              curve. The curve with these vertices and the knot-vector
*              ps1->in2 is the curve we are looking for.
*
*
* REFERENCES :
*
*-
* CALLS      : s1221     - Evaluate curve in given parameter value.
*              newCurve  - Create and initialize new curve.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* REWISED BY : Per Evensen,  SI, 89-3; Prepared for rational description.
* REWISED BY : Per Evensen,  SI, 90-9; Corrected arguments in last call to newCurve.
* MODIFIED BY : Mike Floater,  SI, 91-01; Use rcoef instead of ecoef if rational.
*
*********************************************************************
*/                                     
{
  int kstat = 0;     /* Local status variable.                        */
  int kpos = 0;      /* Position of error.                            */
  int kder = 0;      /* Number of derivatives of curve to evalutate.  */
  int kleft = 0;     /* Parameter used in evaluation of curve.        */
  int kind = 0;      /* Kind of curve
                         = 1 : Polynomial B-spline curve.
                         = 2 : Rational B-spline curve.
                         = 3 : Polynomial Bezier curve.
                         = 4 : Rational Bezier curve.                 */
  int kdim;
  double *scoef;     /* Pointer to vertices.                          */
  double *scurve = SISL_NULL; /* Vertices of constant parameter curve.     */
  SISLCurve *qc = SISL_NULL;  /* Pointer to intermediate curve.         */

  /* Prepare for rational description. */

  kdim = ps1->idim;
  kind = ps1->ikind;
  if(ps1->ikind == 2 || ps1->ikind == 4)
  {
      scoef = ps1->rcoef;
      kdim = kdim+1;
  }
  else
  {
      scoef = ps1->ecoef;
  }

  /* Create the curve describing the surface as a curve.  */
  
  if ((qc = newCurve(ps1->in2,ps1->ik2,ps1->et2,scoef,1,
		     ps1->in1*kdim,0)) == SISL_NULL) goto err101;
  
  /* Allocate space for value of curve.  */
  
  if ((scurve = newarray(ps1->in1*kdim,double)) == SISL_NULL) goto err101;
  
  /* Evaluate this curve at given parameter value.  */
  
  s1221(qc,kder,apar,&kleft,scurve,&kstat);
  if (kstat < 0) goto error;
  
  /* Create constant parameter curve.  */
  
  *rcurve = newCurve(ps1->in1,ps1->ik1,ps1->et1,scurve,kind,ps1->idim,1);
  if (*rcurve == SISL_NULL) goto err101;
  
  /* Set periodicity flag.      */
	
  (*rcurve)->cuopen = ps1->cuopen_1;
  
  /* Curve picked.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1436",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1436",*jstat,kpos);
  goto out;
  
 out: 
  
  /* Free space occupied by local arrays.  */
  
  if (scurve != SISL_NULL) freearray(scurve);
  if (qc != SISL_NULL) freeCurve(qc);
  
  return;
}
