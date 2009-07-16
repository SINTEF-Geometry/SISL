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
 * $Id: sh1460.c,v 1.3 2001-03-19 15:59:03 afr Exp $
 *
 */


#define SH1460

#include "sislP.h"


typedef void (*fshapeProc)(
#if defined(SISLNEEDPROTOTYPES)
                           double [],
  			   double [],
                           int,
                           int,
                           int *
#endif
);

typedef  void (*fevalmidProc)(
#if defined(SISLNEEDPROTOTYPES)
                              fshapeProc,
                              SISLCurve *[],
                              int,
                              double [],
                              double [],
                              double [],
                              int *
#endif
);

#if defined(SISLNEEDPROTOTYPES)
void sh1460(fshapeProc fshape,SISLCurve *vboundc[],int icurv,
	    SISLSurf ***wsurf,int *jstat)
#else
void sh1460(fshape,vboundc,icurv,wsurf,jstat)
     fshapeProc fshape;
     SISLCurve  *vboundc[];
     int        icurv;
     SISLSurf   ***wsurf;
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE   : To create a first derivative geometry continuous blend over a
*             a 3-, 4-, 5- or 6-sided region in space. The boundary of the
*             region are B-spline curves and the cross boundary derivatives
*             are given as B-spline curves.
*
*
* INPUT     : fshape    - Application driven routine that gives the user an
*                         ability to change the middle point of the region
*                         (the vertex at which the blending surfaces meet),
*                         and the tangent vectors in the middle point along
*                         the curves which divide the region.
*             vboundc   - Pointers to the boundary curves. All curves must
*                         have be parameterized counter clockwise around
*                         the patch.
*                         For each edge a position curve and a cross derivative
*                         curve is given. The sequence is the following :
*                         Pointer to position curve of first edge, pointer to
*                         cross derivative curve of first edge, pointer
*                         to position
*                         curve of second edge etc. Dimension is 2*icurv.
*             icurv     - (3, 4, 5 or 6), Number of edges of the region.
*
*
* OUTPUT    : wsurf     - wsurf[0:icurve-1] are pointers to
*                         the blending surfaces
*             jstat     - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* Use        : void shape();
*              SISLCurve uboundc[KCURVE];
*              SISLSurf  usurf[KSURF];
*              int jstat,icurve;
*
*              icurve = number of edges of the region (i.e. 3, 4, 5, or 6).
*              uboundc[0] = curve one;
*              .
*              .
*
*              sh1460(shape,uboundc,icurve,usurf,&jstat);
*
*
* Method     : Hahn's method. Split the region in n 4-sided regions. Define
*              position and cross tangent conditions along the inner edges
*              such that G1-continuity is satisfied. Define the blending
*              patches over the 4-sided regions as Coons patches. The order
*              of the blending surfaces is equal to 7 in one parameter
*              direction and 5 in the other.
*
*
* References : Joerg Hahn : "Filling Polygonal Holes with Rectangular Patches"
*              Theory and Practice of Geometric Modeling, Blackburn Oct 1988.
*
*-
* Calls      : sh1461 - Routine performing Hahn's method.
*              sh1462 - Evaluate midpoint of 3-sided region.
*              sh1463 - Evaluate midpoint of 4-sided region.
*              sh1464 - Evaluate midpoint of 5-sided region.
*              sh1465 - Evaluate midpoint of 6-sided region.
*              s6err - Error messages.
*
* Written by : Vibeke Skytt, SI, May 90.
*
*********************************************************************
*/
{
  int kstat = 0;         /* Status variable   */
  int kpos = 0;          /* Position of error */

  /* Define pointer to midpoint eveluator of vertex region. */

  fevalmidProc fevalmid;
/*
 #if defined(SISLNEEDPROTOTYPES)
   void (*fevalmid)(fshapeProc,SISLCurve *[],int,double [],double [],
 		   double [],int *);
 #else
   void (*fevalmid)();
 #endif
 */
  /* Allocate scratch for output surfaces.  */

  *wsurf = SISL_NULL;
  if ((*wsurf = newarray(icurv,SISLSurf*)) == SISL_NULL) goto err101;

  /* Perform Hahn's method. Pass the routine evaluating the midpoint of the
     region as a parameter to the routine performing Hahn's method.         */

  if (icurv == 3)
    {
      /* 3-sided region.  First set pointer to midpoint evaluator. */

      fevalmid = sh1462;
      sh1461(fshape,fevalmid,vboundc,icurv,*wsurf,&kstat);
      if (kstat < 0) goto error;
    }
  else if (icurv == 4)
    {
      /* 4-sided region.  First set pointer to midpoint evaluator. */

       fevalmid = sh1463 ;
       sh1461(fshape,fevalmid,vboundc,icurv,*wsurf,&kstat);
       if (kstat < 0) goto error;
    }
  else if (icurv == 5)
    {
      /* 5-sided region.  First set pointer to midpoint evaluator. */

      fevalmid = sh1464;
      sh1461(fshape,fevalmid,vboundc,icurv,*wsurf,&kstat);
      if (kstat < 0) goto error;
    }
  else if (icurv == 6)
    {
      /* 6-sided region.  First set pointer to midpoint evaluator. */

      fevalmid = sh1465;
      sh1461(fshape,fevalmid,vboundc,icurv,*wsurf,&kstat);
      if (kstat < 0) goto error;
    }
  else
    {
      /* Illegal number of edges, not implemented.  */

      goto err105;
    }

  /* Blending performed.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

  err101:
    *jstat = -101;
    s6err("sh1460",*jstat,kpos);
    goto out;

  /* Error in input. Wrong number of edges.  */

  err105 :
    *jstat = -105;
    s6err("sh1460",*jstat,kpos);
    goto out;

  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
    s6err("sh1460",*jstat,kpos);
    goto out;

  out :
    return;
}
