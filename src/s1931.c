/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: s1931.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1931

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1931 (int inbcrv, SISLCurve ** vpcrv, double **gknot2,
       double **gcoef2, int *jn2, int *jord2, int *jstat)
#else
void
s1931 (inbcrv, vpcrv, gknot2, gcoef2, jn2, jord2, jstat)
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
* PURPOSE    : Given a set of curve, put these on a common basis.
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
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              Error message corrected
*
*********************************************************************
*/
{
  int ki;                       /* Counter.                               */
  double tstart;            	/* Start parameter-value of B-spline basis
			         * to be made 			 	  */
  double tstop;          	/* Stop parameter-value of B-spline basis
			         * to be made			  	  */
  int kstat = 0;		/* Status variable.                         */
  int kpos = 0;			/* Position of error.                       */
  SISLCurve **tmp_vpcrv = NULL;  /* Temporary array for curve pointers   */
  SISLCurve *pcrv=NULL;

  *jstat = 0;

  /* VSK, 10.92. Set parameter interval in the curve direction as
     the medium of the parameter interval of the input curves. */
  
  for (tstart=DNULL, tstop=DNULL, ki=0; ki<inbcrv; ki++)
  {
     tstart += *(vpcrv[ki]->et + vpcrv[ki]->ik - 1);
     tstop += *(vpcrv[ki]->et + vpcrv[ki]->in);
  }
  tstart /= (double)inbcrv;
  tstop /= (double)inbcrv;
  
  /* s1349 make all the curves k-regular, copy curves not to destroy cyclic bases */
  
  tmp_vpcrv = new0array(inbcrv, SISLCurve*);
  if (tmp_vpcrv == NULL) goto err101;
  
  /* Copy all curves */
  for (ki=0 ; ki<inbcrv ; ki++)
    { 
       pcrv = NULL;
       pcrv = newCurve(vpcrv[ki]->in,vpcrv[ki]->ik,vpcrv[ki]->et,vpcrv[ki]->ecoef,
		       vpcrv[ki]->ikind,vpcrv[ki]->idim,1);
       if (pcrv==NULL) goto err101;
       tmp_vpcrv[ki] = pcrv;	       
    }  

  /* Be sure that all curves have got an open description. */

  s1349 (inbcrv, tmp_vpcrv, &kstat);
  if (kstat < 0)  goto error;

  /* Find common basis for all B-spline curves. */

  s1933 (inbcrv, tmp_vpcrv, tstart, tstop, gknot2, jn2, jord2, &kstat);
  if (kstat < 0) goto error;

  /* Express the curves in the already found basis. */

  s1932 (inbcrv, tmp_vpcrv, tstart, tstop, *gknot2, *jn2, *jord2, gcoef2, &kstat);
  if (kstat < 0) goto error;

  goto out;

  err101:
    *jstat = -101;
    s6err("s1931",*jstat,kpos); 
    goto out;
  
  /* Error in lower level routine.  */

  error:
    *jstat = kstat;
    s6err ("s1931", *jstat, kpos);
    goto out;

  out:
    /* Release allocated curve pointer array and curves */

    if (tmp_vpcrv != NULL)
    {      
      for (ki=0 ; ki<inbcrv ; ki++)
	{
	  if (tmp_vpcrv[ki] != NULL) freeCurve(tmp_vpcrv[ki]); 
	}
      if (tmp_vpcrv != NULL) freearray(tmp_vpcrv); 
    }  
  return;
}
