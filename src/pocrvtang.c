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
 * $Id: pocrvtang.c,v 1.2 2001-03-19 15:58:40 afr Exp $
 *
 */


#define PO_CRV_TANG

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   po_crv_tang(SISLCurve *pcurve, double point[], double ang_tol,
	    double guess_par, double *iter_par, int *jstat)
#else
      void po_crv_tang(pcurve, point, ang_tol, guess_par, iter_par, jstat)
	 SISLCurve   *pcurve;
	 double       point[];
	 double       ang_tol;
	 double       guess_par; 
	 double      *iter_par;
	 int         *jstat;
#endif
/********************************************************************
*                                                                   
* PURPOSE    : Find parameter value on a B-Spline curve for
*              the construction of a tangent from a point to the
*              curve (in 2D).
*
*
*
* INPUT      : pcurve    - Pointer to the curve.
*              point     - The point.
*              ang_tol   - The anular tolerance (in radians).
*              guess_par - Start parameter value on the B-spl crv.
*
*
* OUTPUT     : iter_par  - The found parameter value on the B-spl crv.
*              jstat     - status messages  
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Convert the problem to a onedimention zero problem
*              by using s1893.
*              
*
*
* REFERENCES :
*
*- 
* CALLS      : 
*
* WRITTEN BY : Ulf J. Krystad, SI, 92-03.
*
*********************************************************************
*/                                                               
{                                                                     
  int kstat   = 0;             /* Local status variable.                       */
  int kpos    = 0;             /* Position of error.                           */
  int kdim    = 2;             /* Legal dimension.                             */
  int der0    = 0;             /* Derivative indicator.                        */  
  int der1    = 1;             /* Derivative indicator.                        */  
  int kleft   = 0;             /* Knot navigator.                              */  
  int narr    = 1;             /* No of parallell problems (s1893)             */  
  double zero = 0;             /* Value to iteration.                          */ 
  double iter_ang;             /* The error angle                              */
  double diff[2];              /* Difference vector.                           */
  double sder[4];              /* Result from s1221.                           */
  double trans_arr[9];         /* Matrix describe the curve product (s1893)    */
  SISLCurve *testcurve = SISL_NULL; /* The curve to iterate on.                     */
  SISLPoint *p1 = newPoint (&zero, 1, 0);
  /* ------------------------------------------------------------------------- */
  
  /* Test input. ---------------------------------------------- */
  if (!p1) goto err101;
  if (pcurve->idim != kdim) goto err105;
  
  /* Set up matrix */
  trans_arr[0] = DZERO;
  trans_arr[1] = -1;
  trans_arr[2] = DZERO;
  trans_arr[3] = 1;
  trans_arr[4] = DZERO;
  trans_arr[5] = DZERO;
  trans_arr[6] = -point[1];
  trans_arr[7] =  point[0];
  trans_arr[8] = 1;

  
  /* Convert the problem. ------------------------------------- */
  s1893(pcurve, trans_arr,kdim+1,narr=1,der0=0,der1=1, &testcurve, &kstat);
  if (kstat < 0) goto error;
  
  
  /* Iteration. ---------------------------------------------- */
  s1771(p1, testcurve, REL_COMP_RES, testcurve->et[testcurve->ik - 1],
	testcurve->et[testcurve->in], guess_par, iter_par, &kstat); 
  if (kstat < 0) goto error;
  
  /* Analyse the result */
  s1221(pcurve, der1=1, *iter_par, &kleft, sder, &kstat);
  if (kstat < 0) goto error;

  s6diff(sder,point,kdim,diff);
  iter_ang = s6ang(sder+2, diff, kdim);
  
  if (iter_ang < ang_tol)   *jstat = 1;
  else   *jstat = 2;

  goto out;
  

  /* EXITS. ------------------------------------------------ */

  /* Error in space allocation.  */
 err101: *jstat = -101;
  s6err("po_crv_tang",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 2.  */
 err105: *jstat = -105;
  s6err("po_crv_tang",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  error : *jstat = kstat;
  s6err("po_crv_tang",*jstat,kpos);
  goto out;
  
 out:
  /* Free allocated space.  */
  if (testcurve) freeCurve(testcurve);
  if (p1)        freePoint(p1);
 
  return;
}                                               
                                           
                       
