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
 * $Id: sh1463.c,v 1.2 2001-03-19 15:59:04 afr Exp $
 *
 */


#define SH1463

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
typedef void (*fshapeProc)(double [],double [],int,int,int *);
#else
typedef void (*fshapeProc)();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  sh1463(fshapeProc fshape,
	 SISLCurve *vboundc[],int icurv,double etwist[],
	 double etang[],double eder[],int *jstat)
#else	 
void sh1463(fshape,vboundc,icurv,etwist,etang,eder,jstat)
     fshapeProc    fshape;
     double        etwist[],etang[],eder[];
     SISLCurve     *vboundc[];
     int           icurv,*jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given a four sided vertex region, evaluate the first
*              blending surface in the corner lying in the middle of the
*              vertex region. Compute the tangent vectors in the middle 
*              vertex along the inner boundaries of the region.
*
*
*
* INPUT      : fshape  - Application driven routine that gives the user an
*                        ability to change the middle point of the region
*                        (the vertex at which the blending surfaces meet),
*                        and the tangent vectors in the middle point along
*                        the curves which divedes the region. 
*              vboundc - Position and cross-tangent curves around the vertex
*                        region. For each edge of the region position and cross-
*                        tangent curves are given. The curves follow each other
*                        around the region and are oriented counter-clock-wise.
*                        The dimension of the array is 8.
*              icurv   - Number of sides. icurv = 4.
*              etwist  - Twist-vectors of the corners of the vertex region. The
*                        first element of the array is the twist in the corner
*                        before the first edge, etc. The dimension of the array
*                        is 4*kdim.
*                       
*
* OUTPUT     : etang   - Tangent vectors at the midpoint of the vertex region.
*                        The dimension is icurv*idim.
*              eder    - Value, first and second derivative of the first blending
*                        surface in the corner at the midpoint. The sequence is the
*                        following : Value, 1. derivative in 1. parameter direction,
*                        1. derivative in the 2. parameter direction, 2. derivative
*                        in the 1. parameter direction, mixed derivative and 2.
*                        derivative in the 2. parameter direction. Dimension 6*idim.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Perform rectangular blending of the whole region. Evaluate
*              the resulting surface in the midpoint.
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1401 - Rectangular blending. 
*             s1421        - Surface evaluation.   
*             s1706        - Turn orientation of curve. 
*             freeCurve    - Free space occupied by a curve.
*             newCurve     - Create new curve.     
*
* WRITTEN BY : Vibeke Skytt, SI, 05.90.
*
*********************************************************************
*/
{
  int kstat = 0;     /* Status variable.  */
  int kder = 2;      /* Number of derivatives of surface to evaluate. */
  int ki;            /* Counters.         */
  int kdim = 3;      /* Dimension of geometry space. */
  int klfs = 0;      /* Parameter used in surface evaluation. */
  int klft = 0;      /* Parameter used in surface evaluation. */
  double spar[2];    /* Parameter value of midpoint of rectanular patch. */
  double sder[18];   /* Value and derivatives of patch in the midpoint.  */
  double snorm[3];   /* Normal of patch in the midpoint.      */
  SISLSurf *qsurf = SISL_NULL;  /* Rectangular blending patch.  */
  SISLCurve *qc[8];        /* Copy of edge curves.         */
  SISLCurve *qpt;          /* Pointer to curve.            */

  /* Make a copy of the edge curves.  */

  for (ki=0; ki<8; ki++)
    {
      qpt = vboundc[ki];
      
      /* Test dimension of curves. */

      if (qpt->idim != kdim) goto err104;
      
      qc[ki] = newCurve(qpt->in,qpt->ik,qpt->et,qpt->ecoef,qpt->ikind,kdim,1);
      if (qc[ki] == SISL_NULL) goto err101;
    }
    
  /* Turn the orientation of the curves at the 3. and 4. edge. */

  s1706(qc[4]);
  s1706(qc[5]);  
  s1706(qc[6]);  
  s1706(qc[7]);  

  /* Compute Coon's patch.  */

  s1401(qc,etwist,&qsurf,&kstat);
  if (kstat < 0) goto error;
   
  /* Evaluate the surface in the middle.  */

  spar[0] = (double)0.5*(*(qsurf->et1+qsurf->ik1-1) + *(qsurf->et1+qsurf->in1));
  spar[1] = (double)0.5*(*(qsurf->et2+qsurf->ik2-1) + *(qsurf->et2+qsurf->in2));
  
  s1421(qsurf,kder,spar,&klfs,&klft,sder,snorm,&kstat);
  if (kstat < 0) goto error;
  
  /* Compute tangent vectors at the midpoint of the vertex region as 
     half the tangent vectors at the rectangular patch.  */

  for (ki=0; ki<kdim; ki++)
    {
      etang[ki] = (double)0.5*sder[2*kdim+ki];
      etang[kdim+ki] = (double)0.5*sder[kdim+ki];
      etang[2*kdim+ki] = -(double)0.5*sder[2*kdim+ki];
      etang[3*kdim+ki] = -(double)0.5*sder[kdim+ki];
    }
  
  /* Application driven routine to alter the midpoint and tangents in the
     midpoint.  */

  fshape(sder,etang,kdim,icurv,&kstat);
  if (kstat < 0) goto error;
  
  /* Copy value and 1. derivatives of first patch.  */

  memcopy(eder,sder,kdim,DOUBLE);
  memcopy(eder+kdim,etang+3*kdim,kdim,DOUBLE);
  memcopy(eder+2*kdim,etang,kdim,DOUBLE);
  
  /* Compute 2. derivatives.  */

  for (ki=0; ki<kdim; ki++)
    {
      eder[3*kdim+ki] = (double)0.1*sder[5*kdim+ki];
      eder[4*kdim+ki] = (double)0.1*sder[4*kdim+ki];
      eder[5*kdim+ki] = (double)0.1*sder[3*kdim+ki];
    }
  
  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation. */

  err101 :
  *jstat = -101;
  goto out;
 
  /* Error in input. Dimension not equal to 3.  */

  err104 :
    *jstat = -104;
  goto out;
   
  /* Error in a lower level function.  */

 error:
  *jstat = kstat;
  goto out;
  
  out :

    /* Free space occupied by local curves and surface. */

    for (ki=0; ki<8; ki++)
      if (qc[ki] != SISL_NULL) freeCurve(qc[ki]);
  if (qsurf != SISL_NULL) freeSurf(qsurf);
  
    return;
}
