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
 * $Id: s1953.c,v 1.3 2001-03-19 15:58:57 afr Exp $
 *
 */


#define S1953

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1953(SISLCurve *pcurve,double epoint[],int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1953(pcurve,epoint,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
     SISLCurve    *pcurve;
     double   epoint[];
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
* PURPOSE    : Find the closest point between the curve pcurve and the
*              point epoint.
*
*
*
* INPUT      : pcurve - Pointer to the curve in the closest point problem.
*              epoint - The point in the closest point problem.
*              idim   - Dimension of the space in which epoint lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single closest points.
*              gpar   - Array containing the parameter values of the
*                       single closest points in the parameter interval
*                       of the curve. The points lie continuous.
*                       Closest curves are stored in wcurve.
*              jcrv   - Number of closest curves.
*              wcurve - Array containing descriptions of the closest
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
* METHOD     : Put the equation of the curve into the equation of a
*              circel/sphere with center epoint and radius zero. The
*              result is a curve of dimension one. Find the minimum
*              points of this curve.
*
*
*
* REFERENCES :
*
*-
* CALLS      : s1321,s1370,s1920,freeCurve.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 1994-08.  Updated to handle
*              input curve with periodic basis correctly.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                    */
  int kpos = 0;             /* Position of error.                        */
  int kdim = 1;             /* Dimension of curve in extremal problem.   */
  double tradius = 0;       /* Radius of circel/sphere describing point. */
  double tdir = -1;         /* Direction of extremal value.              */
  double sarray[16];        /* Matrix describing circel/sphere.          */
  SISLCurve *qc = SISL_NULL;     /* Curve of which to find extremal points.   */
  SISLCurve *qkreg = SISL_NULL;  /* Input curve with ensured k-regular basis. */
  int ratflag = 0;          /* Flag to indicate if curve is rational.    */
  int ki;                   /* Counter.                                  */


  /* Ensure input has k-regular knot basis. */

  if ( pcurve -> cuopen == SISL_CRV_PERIODIC )
  {
    /* Cyclic (periodic) curve */

    make_cv_kreg(pcurve, &qkreg, &kstat);
    if ( kstat < 0 )  goto error;
  }
  else
    qkreg = pcurve;


  /* Test input.  */

  if ( idim != 2 && idim != 3 )  goto err105;
  if ( qkreg -> idim != idim )  goto err106;

  if ( qkreg -> ikind == 2 || qkreg -> ikind == 4 )  ratflag = 1;

  /* Make a matrix of dimension (idim+1)*(idim+1) to describe
     the circel/shpere.                                        */

  s1321(epoint, tradius, idim, kdim, sarray, &kstat);
  if ( kstat < 0 )  goto error;

  /* Put curve into circel/sphere equation.  */

  s1370(qkreg, sarray, idim, kdim, ratflag, &qc, &kstat);
  if ( kstat < 0 )  goto error;

  /* Find minimum points of the new curve.  */

  s1920(qc, &tdir, kdim, aepsco, aepsge, jpt, gpar, jcrv, wcurve, &kstat);
  if ( kstat < 0 )  goto error;

  /* Handle periodicity (remove extraneous points) */
  if ( *jpt > 1  &&  idim > 1  &&  pcurve -> cuopen == SISL_CRV_PERIODIC )
  {
    for ( ki=0;  ki < (*jpt);  ki++ )
    {
      if ( (*gpar)[ki] == pcurve -> et[pcurve->in] )
      {
	(*jpt)--;
	(*gpar)[ki] = (*gpar)[*jpt];
	ki--;
      }
    }
  }

  /* Closest points/intervals found.  */

  *jstat = 0;
  goto out;


  /* Error in input. Dimension not equal to 2 or 3.  */

 err105:
  *jstat = -105;
  s6err("s1953",*jstat,kpos);
  goto out;

  /* Dimensions conflicting.  */

 err106:
  *jstat = -106;
  s6err("s1953",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

 error:
  *jstat = kstat;
  s6err("s1953",*jstat,kpos);
  goto out;

 out:

  /* Free allocated space.  */

  if ( qkreg  &&  qkreg != pcurve )  freeCurve(qkreg);
  if (qc)  freeCurve(qc);

  return;
}
