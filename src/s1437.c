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
 * $Id: s1437.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1437

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1437(SISLSurf *ps1,double apar,SISLCurve **rcurve,int *jstat)
#else
void s1437(ps1,apar,rcurve,jstat)
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
*              first parameter direction.
*
*
*
* INPUT      : ps1    - Surface.
*              apar   - Parameter value to use whe picking out constant
*                       parameter curve in first parameter direction.
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
* METHOD     : First the parameter directions of the surface is 
*              interchanged. Then the surface is viewed as a
*              (ps1->idim*ps1->in2)-dimensional curve with
*              knot vector ps1->et1, and evaluated in apar.
*              The result of the evalutation is when viewed as
*              a curve with ps1->in2 vertices and knot vector
*              ps1->et2, the curve we are looking for.
*
*
* REFERENCES :
*
*-
* CALLS      : s1221    - Evaluate curve in given parameter value.
*              s6chpar  - Change parameter direction of vertices of
*                         surface.
*              newCurve - Create and initialize new curve.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* REWISED BY : Per Evensen,  SI, 89-3; Prepared for rational description.
* REWISED BY : Per Evensen,  SI, 90-9; Corrected arguments in last call to newCurve.
*
*********************************************************************
*/                                     
{
  int kstat = 0;     /* Local status variable.                           */
  int kpos = 0;      /* Position of error.                               */
  int kind = 0;      /* Kind of curve
                         = 1 : Polynomial B-spline curve.
                         = 2 : Rational B-spline curve.
                         = 3 : Polynomial Bezier curve.
                         = 4 : Rational Bezier curve.                 */
  int kdim;          /* Dimension of the space in which the surface lies.*/
  int kder = 0;      /* Number of derivatives of curve to evaluate.      */
  int kleft = 0;     /* Parameter used in evalutation of curve.          */
  double *ecoef = SISL_NULL;  /* Pointer to vertices                          */
  double *scoef = SISL_NULL;  /* Vertices of surface with changed parameter
			    directions.                                  */
  double *scurve = SISL_NULL; /* Vertices of constant parameter curve.        */
  SISLCurve *qc = SISL_NULL;  /* Intermediate curve.                       */
  
  /* Get dimension of space.  */
  
  kdim = ps1 -> idim;
  kind = ps1->ikind;
                       
  /* Prepare for rational description. */

  if(ps1->ikind == 2 || ps1->ikind == 4)
  {
      ecoef = ps1->rcoef;
      kdim = kdim+1;
  }
  else
  {
      ecoef = ps1->ecoef;
  }

  
  /* Allocate space for coefficients of constant parameter curve
     and of surface with changed parameter direction.                */
  
  if ((scurve = newarray(kdim*ps1->in2,double)) == SISL_NULL) goto err101;
  if ((scoef = newarray(kdim*ps1->in1*ps1->in2,double)) == SISL_NULL) goto err101;
  
  /* Change parameter directions of surface.  */
  
  s6chpar(ecoef,ps1->in1,ps1->in2,kdim,scoef);
  
  /* Create curve to evaluate.  */
  
  qc = newCurve(ps1->in1,ps1->ik1,ps1->et1,scoef,1,kdim*ps1->in2,0);
  if (qc == SISL_NULL) goto err101;
  
  /* Evaluate this curve in given parameter value.  */
  
  s1221(qc,kder,apar,&kleft,scurve,&kstat);
  if (kstat < 0) goto error;
  
  /* Create constant parameter curve.  */
  
  *rcurve = newCurve(ps1->in2,ps1->ik2,ps1->et2,scurve,kind,ps1->idim,1);
  if (*rcurve == SISL_NULL) goto err101;
  
  /* Set periodicity flag.      */
		       
  (*rcurve)->cuopen = ps1->cuopen_2;		       
  
  /* Curve picked.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1437",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1437",*jstat,kpos);
  goto out;
  
 out: 
  
  /* Free space occupied by local arrays.  */
  
  if (scoef != SISL_NULL) freearray(scoef);
  if (scurve != SISL_NULL) freearray(scurve);
  if (qc != SISL_NULL) freeCurve(qc);
  
  return;
}
