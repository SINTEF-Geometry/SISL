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
 * $Id: s1931unit.c,v 1.2 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1931UNIT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1931unit (int inbcrv, SISLCurve ** vpcrv, double **gknot2,
       double **gcoef2, int *jn2, int *jord2, int *jstat)
#else
void
   s1931unit (inbcrv, vpcrv, gknot2, gcoef2, jn2, jord2, jstat)
     int inbcrv;
     SISLCurve **vpcrv;
     double **gknot2;
     double **gcoef2;
     int *jn2;
     int *jord2;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : Given a set of curve, put these on a common basis
*              using a unit knot vector.
*              (Common knot-vector of length (jn2+jord2).
*              The vertices are recomputed according to this new basis.
*
* INPUT      : inbcrv - Number of curves in the curve-set.
*              vpcrv  - Array (length inbcrv) of pointers to the
*                       curves in the curve-set.
*
* OUTPUT     : gknot2 - Common knot-vector (new basis) for the curves.
*                       (jn2+jord2).
*              gcoef2 - The vertices of the inbcrv curves
*                       expressed in the new basis. Stored in sequence
*                       first curve, second curve,...
*              jn2    - The no. of vertices in each of the inbcrv curves.
*              jord2  - The order of the new representation of the curves.
*              jstat  - Output status:
*                        < 0: Error.
*                        = 0: Ok.
*                        > 0: Warning.
*
* NOTE	     : The maximal order of the B-spline basis in the lofting
*	       direction is no longer used. In earlier version (s1337)
*	       the maximal order iord1 was not found, nor calculated.
*
* CALLS      : s1349,s1933,s1932,s6err.
*
* WRITTEN BY : Christophe R. Birkeland, SI, 1991-07
*
*********************************************************************
*/
{
   int ki;                      /* Counter.                               */
  double tstart = 0;           	/* Start parameter-value of B-spline basis
			         * to be made 			 	  */
  double tstop = 1;          	/* Stop parameter-value of B-spline basis
			         * to be made			  	  */
  int kstat = 0;		/* Status variable.                         */
  int kpos = 0;			/* Position of error.                       */
  
  SISLCurve **tmp_vpcrv = SISL_NULL;  /* Temporary array for curve pointers   */
  SISLCurve *pcrv=SISL_NULL;

  *jstat = 0;

  /* s1349 make all the curves k-regular, copy curves not to destroy cyclic bases */
  
  tmp_vpcrv = new0array(inbcrv, SISLCurve*);
  if (tmp_vpcrv == SISL_NULL) goto err101;
    
  /* Copy all curves */
  for (ki=0 ; ki<inbcrv ; ki++)
    { 
       pcrv = SISL_NULL;
       pcrv = newCurve(vpcrv[ki]->in,vpcrv[ki]->ik,vpcrv[ki]->et,vpcrv[ki]->ecoef,
		       vpcrv[ki]->ikind,vpcrv[ki]->idim,1);
       if (pcrv==SISL_NULL) goto err101;
       tmp_vpcrv[ki] = pcrv;	       
    }  

  /* Be sure that all curves have got an open description. */

  s1349 (inbcrv, tmp_vpcrv, &kstat);
  if (kstat < 0) goto error;

  /* Find common basis for all B-spline curves. */

  s1933 (inbcrv, tmp_vpcrv, tstart, tstop, gknot2, jn2, jord2, &kstat);
  if (kstat < 0) goto error;

  /* Express the curves in the already found basis. */

  s1932 (inbcrv, tmp_vpcrv, tstart, tstop, *gknot2, *jn2, *jord2, gcoef2, &kstat);
  if (kstat < 0) goto error;

  goto out;

  err101:
    *jstat = -101;
    s6err("s1931unit",*jstat,kpos);
    goto out;

  /* Error in lower level routine.  */

  error:
    *jstat = kstat;
    s6err ("s1931unit", *jstat, kpos);
    goto out;

  out:
    /* Release allocated curve pointer array and curves */

    if (tmp_vpcrv != SISL_NULL)
    {
      
      for (ki=0 ; ki<inbcrv ; ki++)
	{
	  if (tmp_vpcrv[ki] != SISL_NULL) freeCurve(tmp_vpcrv[ki]); 
	}
      if (tmp_vpcrv != SISL_NULL) freearray(tmp_vpcrv);      
    }
 
    return;
}
