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
 * $Id: sh1464.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH1464

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
typedef void (*fshapeProc)(double [],double [],int,int,int *);
#else
typedef void (*fshapeProc)();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  sh1464(fshapeProc fshape,SISLCurve *vboundc[],int icurv,
	 double etwist[],double etang[],double eder[],int *jstat)
#else      
void sh1464(fshape,vboundc,icurv,etwist,etang,eder,jstat)
     fshapeProc    fshape;
     double        etwist[],etang[],eder[];
     SISLCurve     *vboundc[];
     int           icurv,*jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given a five sided vertex region, evaluate the first
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
*                        The dimension of the array is 10.
*              icurv   - Number of sides. icurv = 5.
*              etwist  - Twist-vectors of the corners of the vertex region. The
*                        first element of the array is the twist in the corner
*                        before the first edge, etc. The dimension of the array
*                        is icurve*kdim.
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
* METHOD     : Evaluate the Gregory Charrot function in the midpoint of the
*              vertex region. Compute the wanted derivatives.
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : sh1467  - Evaluate Gregory Charrot function over 
*			 5-sided region.                       
*
* WRITTEN BY : Vibeke Skytt, SI, 03.90.
*
*********************************************************************
*/
{
  int kstat = 0;         /* Status variable.    */
  int kder = 2;          /* Number of derviatives to evaluate.  */
  int ki;                /* Counter.  */
  int kdim = 3;          /* Dimension of geometry space.     */
  double tlambda = (double)1.0/sqrt((double)5.0);           /* Constant.  */
  double tl1 = (double)2.0*tlambda*tan(PI/(double)5.0);     /* Constant.  */
  double tl2 = sin((double)0.3*PI);                         /* Constant.  */
  double tconst1 = tl1/(double)2.0 - tlambda;               /* Constant.  */
  double tconst2 = tl2 - tlambda;                           /* Constant.  */
  double sbar[5];        /* Barycentric coordinates of the blending function. */
  double sder[18];       /* Value and derivatives of blending function.       */

  /* Set up the barycentric coordinates of the midpoint of the region. */

  sbar[0] = sbar[1] = sbar[2] = sbar[3] = sbar[4] = tlambda;
  
  /* Evaluate the Gregory Charrot function at the midpoint. */

  sh1467(vboundc,etwist,kder,sbar,sder,&kstat);
  if (kstat < 0) goto error;
   
  /* Compute tangent vectors.  */

  for (ki=0; ki<kdim; ki++)
    {
      etang[ki] = sder[kdim+ki]*tconst1 + sder[2*kdim+ki]*tconst2;
      etang[kdim+ki] = -sder[kdim+ki]*tlambda + sder[2*kdim+ki]*tconst1;
      etang[2*kdim+ki] = sder[kdim+ki]*tconst1 - sder[2*kdim+ki]*tlambda;
      etang[3*kdim+ki] = sder[kdim+ki]*tconst2 + sder[2*kdim+ki]*tconst1;
      etang[4*kdim+ki] = sder[kdim+ki]*tconst2 + sder[2*kdim+ki]*tconst2;
    }
  
  /* Application driven routine to alter the midpoint and tangents in the
     midpoint.  */

  fshape(sder,etang,kdim,icurv,&kstat);
  if (kstat < 0) goto error;
  
  /* Copy value and 1. derivatives of first patch.  */

  memcopy(eder,sder,kdim,DOUBLE);
  memcopy(eder+kdim,etang+4*kdim,kdim,DOUBLE);
  memcopy(eder+2*kdim,etang,kdim,DOUBLE);
  
  /* Compute 2. derivatives.  */

  for (ki=0; ki<kdim; ki++)
    {
      eder[3*kdim+ki] = sder[3*kdim+ki]*tconst2*tconst2 
	+ (double)2.0*sder[4*kdim+ki]*tconst2*tconst2 
	  + sder[5*kdim+ki]*tconst2*tconst2;
      eder[4*kdim+ki] = sder[3*kdim+ki]*tconst1*tconst2 
	+ sder[4*kdim+ki]*tconst2*(tconst1+tconst2)
	  + sder[5*kdim+ki]*tconst2*tconst2;
      eder[5*kdim+ki] = sder[3*kdim+ki]*tconst1*tconst1 
	- (double)2.0*sder[4*kdim+ki]*tconst1*tconst2 
	  + sder[5*kdim+ki]*tconst2*tconst2;
    }
  
  *jstat = 0;
  goto out;
  
  /* Error in a lower level function.  */

 error:
  *jstat = kstat;
  goto out;
  
  out :
    return;
}
