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
 * $Id: s1239.c,v 1.3 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1239

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/
#if defined(SISLNEEDPROTOTYPES)
static void s1239_s9sort(double [],int [],int);
#else
static void s1239_s9sort();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
s1239(SISLCurve *pcpar,int ipar,double apar,SISLCurve *pcurve,double aepsco,
	   double aepsge,SISLCurve **vpartc,int imax,int *jpartc,int *jstat)
#else
void s1239(pcpar,ipar,apar,pcurve,aepsco,aepsge,vpartc,imax,jpartc,jstat)
     SISLCurve  *pcpar;
     int    ipar;
     double apar;
     SISLCurve  *pcurve;
     double aepsco;
     double aepsge;
     SISLCurve  **vpartc;
     int    imax;
     int    *jpartc;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Fetch the parts of the constant parameter curve pcpar
*              that lies within the closed curve pcurve. Used in
*              drawing constant parameter lines of a surface limited
*              by a closed B-spline curve.
*
*
*
* INPUT      : pcpar  - Constant parameter curve in a surface.
*              ipar   - Constant parameter direction of the surface
*                       where the curve lies. ipar=0 or ipar=1.
*              apar   - Constant parameeter value of the curve
*                       corresponding to the surface.
*              pcurve - The closed curve that the wanted curve-
*                       segments lie within.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*              imax   - Maximum number of curve-segments that can be
*                       stored.
*
*
*
* OUTPUT     : vpartc - Array containing curve-segments lying within
*                       pcurve.
*              jpartc - Number of curve-segments in wpartc.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s1221,s1227,s1710,s1712,s1850,freeCurve,s1239_s9sort.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* Changed by : Per OEyvind Hvidsten, SINTEF, 94-11.
*              Added free(spt), thus removing a memory leak.
*
*********************************************************************
*/
{
  int kstat = 0;      /* Local status variable.                            */
  int kpos = 0;       /* Position of error.                                */
  int ki,kj;          /* Counters.                                         */
  int kpt;            /* Number of intersection points found.              */
  int kcrv;           /* Number of intersection curves found.              */
  int kleft = 0;      /* Parameter used in curve evaluation.               */
  int kder = 1;       /* Number of derivatives of curve to evaluate.       */
  int kpar = 0;       /* Number of times the closed curve crosses the
			 parameter interval of constant parameter curve.   */
  int *linter = SISL_NULL; /* Indicates kind of intersection point.
			 = 0 : The point does not belong to an interval.
			 = 1 : Belongs to an interval that touches line.
			 = 2 : Belongs to an interval that crosses line.   */
  int kinter;         /* Code for one element of linter.                   */
  int kpartc = 0;     /* Counter of curves to return.                      */
  double tzero = (double)0.0;   /* Zero.                                   */
  double tstart,tend; /* End-parameters of constant parameter curve.       */
  double tclstart,tclend; /* End-parameters of closed curve.               */
  double t1,t2;       /* Help variables.                                   */
  double tchange;     /* Help varaible.                                    */
  double spoint[2];   /* SISLPoint on straight line describing the parameter
			 interval of the constant parameter curve.         */
  double snorm[2];    /* Normal to straight line.                          */
  double *spt = SISL_NULL; /* Intersection points between curve and line.       */
  double sder1[4];    /* Position and derivative of curve.                 */
  double sder2[4];    /* Position and derivative of curve.                 */
  double *spar = SISL_NULL;      /* Points where pcurve crosses the parameter
			       interval of pcpar.                          */
  SISLIntcurve **ucrv;   /* Intersection curves between pcurve and line.*/
  SISLCurve *qc1 = SISL_NULL; /* First part of subdivided curve.             */
  SISLCurve *qc2 = SISL_NULL; /* Second part of subdivided curve.            */

  /* Test input.  */

  if (pcurve -> idim != 2) goto err108;

  /* Fetch ends of parameter interval of closed curve. */

  tclstart = *(pcurve->et + pcurve->ik - 1);
  tclend   = *(pcurve->et + pcurve->in);

  /* Describe the parameter interval of the constant parameter curve
     as a straight line. First fetch endpoints of interval.           */

  tstart = *(pcpar->et + pcpar->ik - 1);
  tend = *(pcpar->et + pcpar->in);

  if (ipar == 0)
    {
      spoint[0] = apar;  spoint[1] = tstart;
      snorm[0] = tend-tstart;   snorm[1] = tzero;
    }
  else
    {
      spoint[0] = tstart;    spoint[1] = apar;
      snorm[0] = tzero;       snorm[1] = tend-tstart;
    }

  /* Find all intersection between this straight line and the curve pcurve.*/

  s1850(pcurve,spoint,snorm,2,aepsco,aepsge,&kpt,&spt,&kcrv,&ucrv,&kstat);
  if (kstat < 0) goto error;

  /* Make sure that only one end-point of the closed curve is represented.*/

  for (kj=0,ki=0; ki<kpt; ki++)
    {
      if (s6equal(spt[ki],tclstart,(double)1.0) ||
	  s6equal(spt[ki],tclend,(double)1.0))
	{
	  if (kj == 0) kj++;
	  else
	    {
	      for (kj=ki+1; kj<kpt; kj++) spt[kj-1] = spt[kj];
	      kpt--;
	    }
	}
    }

  /* Allocate space for curve interval information.  */

  if (kpt+2*kcrv > 0)
    {
      if ((spar = newarray(kpt+2*kcrv,double)) == SISL_NULL) goto err101;
      if ((linter = new0array(kpt+2*kcrv,int)) == SISL_NULL) goto err101;
    }

  /* Discuss intersection intervals.  */

  for (ki=0; ki<kcrv; ki++)
    {

      /* Fetch parameters of endpoints of interval.  */

      t1 = *((*(ucrv+ki))->epar1);
      t2 = *((*(ucrv+ki))->epar1 + (*(ucrv+ki))->ipoint - 1);
      if (t2 < t1)
	{
	  tchange = t1;  t1 = t2;  t2 = tchange;
	}

      /* Evaluate curve in endpoints.  */

      s1227(pcurve,kder,t1,&kleft,sder1,&kstat);
      if (kstat < 0) goto error;

      s1221(pcurve,kder,t2,&kleft,sder2,&kstat);
      if (kstat < 0) goto error;

      /* Test if the curve crosses the constant parameter line during
	 the interval.                                                 */

      if (sder1[2+ipar]*sder2[2+ipar] > tzero)
	kinter = 2;    /* The curve crosses the straight line.  */
      else kinter = 1;  /* The curve only touch the straight line.  */

      /* Get endpoints of the part of the intersection interval that lies
	 inside the parameter interval of the constant parameter curve.   */

      t1 = MAX(MIN(sder1[1-ipar],sder2[1-ipar]),tstart);
      t2 = MIN(MAX(sder1[1-ipar],sder2[1-ipar]),tend);

      /* Remember endpoints of interval.  */

      linter[kpar] = kinter;   spar[kpar++] = t1;
      if (kpar == imax) goto err190;
      linter[kpar] = kinter;   spar[kpar++] = t2;
      if (kpar == imax) goto err190;
    }

  /* Discuss intersection points.  */

  for (ki=0; ki<kpt; ki++)
    {

      /* Compute position and right derivative in point.  */

      s1221(pcurve,kder,spt[ki],&kleft,sder1,&kstat);
      if (kstat < 0) goto error;

      if (!(s6equal(spt[ki],tclstart,(double)1.0) ||
	    s6equal(spt[ki],tclend,(double)1.0)))
	{

	  /* The point is not an end-point. Compute position and left
	     derivative in point.                                      */

	  s1227(pcurve,kder,spt[ki],&kleft,sder2,&kstat);
	  if (kstat < 0) goto error;
	}
      else
	{

	  /* Evaluate curve in the other end-point.  */

	  t1 = (s6equal(spt[ki],tclstart,(double)1.0)) ? tclend : tclstart;
	  s1221(pcurve,kder,t1,&kleft,sder2,&kstat);
	  if (kstat < 0) goto error;
	}

      if (DEQUAL(sder1[2]*sder2[3] - sder1[3]*sder2[2],(double)0.0))
	{

	  /* The derivative of pcurve is continuous in the point. */

	  if (DEQUAL(sder1[2+ipar],(double)0.0))
	    {

	      /* Possible touching point. Subdivide the curve in the point. */

	      s1710(pcurve,spt[ki],&qc1,&qc2,&kstat);
	      if (kstat < 0) goto error;

	      /* Test if the curve crosses the constant parameter line.  */

	      t1 = apar;  kj = 0;
	      while (qc1 && DEQUAL(t1,apar))
		{
		  kj++;
		  t1 = *(qc1->ecoef+2*(qc1->in-kj-1)+ipar);
		}

	      t2 = apar;  kj = 0;
	      while (qc2 && DEQUAL(t2,apar))
		{
		  kj++;
		  t2 = *(qc2->ecoef+2*kj+ipar);
		}

	      /* Free curve-parts generated at subdivision.  */

	      if (qc1 != SISL_NULL) freeCurve(qc1);  qc1 = SISL_NULL;
	      if (qc2 != SISL_NULL) freeCurve(qc2);  qc2 = SISL_NULL;
	    }

	  if (DNEQUAL(sder1[2+ipar],(double)0.0) ||
	      (t1-apar)*(t2-apar) < (double)0.0)
	    {

	      /* The curve crosses the parameter line.  */

	      t1 = sder1[1-ipar];
	      if (t1 < tstart) t1 = tstart;
	      else if (t1 > tend) t1 = tend;

	      /* Remember point.  */

	      spar[kpar++] = t1;
	      if (kpar == imax) goto err190;
	    }
	}
      else
	{

	  /* Derivative discontinuous in point. Test if the curve crosses
	     the constant parameter line.                                 */

	  if (sder1[2+ipar]*sder2[2+ipar] > (double)0.0)
	    {

	      /* The curve crosses the line. */

	      t1 = sder1[1-ipar];
	      if (t1 < tstart) t1 = tstart;
	      else if (t1 > tend) t1 = tend;

	      /* Remember point.  */

	      spar[kpar++] = t1;
	      if (kpar == imax) goto err190;
	    }
	}
    }

  /* Sort the arrays spar and linter.  */

  s1239_s9sort(spar,linter,kpar);

  /* Pick the curve-segments that is to be returned from this routine. */

  ki = 0;
  while (ki < kpar-1)
    {

      /* Find index of end of curve-segment.  */

      kj = ki + (linter[ki] == 0 && linter[ki+1] > 0) + 1;
      while (kj < kpar-1 && ((ki < kj-1 && linter[kj] == 1) ||
			     (ki == kj-1 && linter[kj] == 2))) kj += ((linter[kj+1] > 0) + 1);
      kj = MIN(kj,kpar-1);

      if (DNEQUAL(spar[kj],spar[ki]))
	{

	  /* Pick the part of the constant parameter curve lying between
	     spar[ki] and spar[kj].                                       */

	  s1712(pcpar,spar[ki],spar[kj],vpartc+kpartc,&kstat);
	  if (kstat < 0) goto error;

	  /* Prepare for next curve segment.  */

	  ki = kj + 1;
	  kpartc++;
	}
      else ki = kj + 1;
    }

  /* Wanted curve-segments picked.  */

  *jpartc = kpartc;
  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

 err101: *jstat = -101;
  s6err("s1239",*jstat,kpos);
  goto out;

  /* Error in input. Dimension not equal to 2.  */

 err108: *jstat = -108;
  s6err("s1239",*jstat,kpos);
  goto out;

  /* Error in storing curves. Array too small.  */

 err190: *jstat = -190;
  s6err("s1239",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("s1239",*jstat,kpos);
  goto out;

 out:

  /* Free space occupied by local arrays.  */

  if (spt != SISL_NULL) free(spt);
  if (spar != SISL_NULL) freearray(spar);
  if (linter != SISL_NULL) freearray(linter);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1239_s9sort(double epar[],int nint[],int ipar)
#else
static void s1239_s9sort(epar,nint,ipar)
     double epar[];
     int    nint[];
     int    ipar;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Sort the array epar and change nint correspondingly.
*
*
*
* INPUT      : ipar   - Number of elements in the arrays epar and nint.
*
*
* INPUT/OUTPUT : epar  - Array to be sorted.
*                nint  - Array where the elements corresponds to the
*                        elements of epar.
*
*********************************************************************
*/
{
  int ki,kj;      /* Counters.                                          */
  int kchange;    /* Help variable used to change two elements of nint. */
  double tchange; /* Help variable used to change two elements of epar. */

  for (ki=0; ki<ipar; ki++)
    for (kj=ki+1; kj<ipar; kj++)
      if (epar[kj] < epar[ki])
	{
	  tchange = epar[ki];  epar[ki] = epar[kj];  epar[kj] = tchange;
	  kchange = nint[ki];  nint[ki] = nint[kj];  nint[kj] = kchange;
	}

  return;
}
