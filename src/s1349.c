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
 * $Id: s1349.c,v 1.3 2001-03-19 15:58:46 afr Exp $
 *
 */


#define S1349

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1349(int inbcrv,SISLCurve *vpcrv[],int *jstat)
#else
void s1349(inbcrv,vpcrv,jstat)
     int   inbcrv;
     SISLCurve *vpcrv[];
     int   *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To convert a B-spline curve with the first ik knots not
*              equal to et(ik) and the last ik knots not equal to
*              et(in+1) to a representation with the first ik knots
*              equal to et(ik), and the last ik knots equal to et(in+1).
*
* INPUT      : inbcrv - No. of curves in the curve-set.
*              vpcrv  - Array (length inbcrv) of pointers to curves
*                       to the curves in the curve-set.
*
* OUTPUT     : jstat  - status messages:
*                        > 0      : warning
*                        = 0      : ok
*                        < 0      : error
*              vpcrv  - Array (length inbcrv) to curves in the curve-set..
*                       (Changed description.)
*
* METHOD     : For each curve, the multiplicity of the end vertices is
*              counted, and those curves not having the right multiplicity
*              are converted to a representation with the right multiplicity.
*-
* CALLS      : s1712,s6err
*
* WRITTEN BY : A. M. Ytrehus SI, Oslo, Norway.  Sep.  1988
* Revised by : Tor Dokken, SI, Oslo, Norway. 26. Feb. 1989
* Revised by : Paal Fugelli, SINTEF, Oslo 02/08-1994. Fixed memory leak.
*
*********************************************************************
*/
{
  SISLCurve **wp = SISL_NULL;      /* Local pointer to current curve.           */
  SISLCurve *qc2 = SISL_NULL;      /* Pointer to new curve-object.              */
  int kvert;              /* No. of vertices in current curve.           */
  int kord;               /* Order of current curve.                     */
  double *sknot = SISL_NULL;   /* Pointer to knot-vector of current curve.    */
  int kk,kr;              /* Loop controllers.                           */
  double *sp1,*sp2;       /* Pointers to sknot.                          */
  int kmul1,kmul2;        /*                                             */
  double tval1,tval2;     /*                                             */
  int kstat = 0;          /* Status variable.                            */
  int kpos = 0;           /* Position of error.                          */

  wp = vpcrv;
  for (kk=0; kk<inbcrv; kk++)
    {
      /* Make local pointers to description of curve */
      sknot = (*wp) -> et;
      kvert = (*wp) -> in;
      kord = (*wp) -> ik;

      /* Count multiplicity of start-knot. */
      kmul1 = 0;
      sp1 = sknot + kord - 1;
      tval1 = *sp1;
      for (kr=0; kr<kord; kr++)
	{
          if (*sp1 == tval1) kmul1++;
          sp1--;
	}

      /* Count multiplicity of end-knot. */
      kmul2 = 0;
      sp2 = sknot + kvert;
      tval2 = *sp2;
      for (kr=0; kr<kord; kr++)
	{
          if (*sp2 == tval2) kmul2++;
          sp2++;
	}

      /* If the multiplicity of both end-knots equals kord, the curve is ok. */
      if (kmul1 != kord || kmul2 != kord)
	{
          /* Both ends do not have multiplicity of order kord.
             Create a new curve-object.                        */
          s1712((*wp),tval1,tval2,&qc2,&kstat);
          if (kstat<0) goto error;
	  if ((*wp)) freeCurve(*wp);  /* PFU 02/08-1994 */
          *wp = qc2;
          qc2 = SISL_NULL;
	}
      wp++;
    }

  *jstat = 0;
  goto out;

/* Error in lower level routine. */
error: *jstat = kstat;
       s6err("s1349",*jstat,kpos);
       goto out;
out:
return;
}
