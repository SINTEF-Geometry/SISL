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
 * $Id: s1932.c,v 1.2 2001-03-19 15:58:56 afr Exp $
 *
 */


#define S1932

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1932 (int inbcrv, SISLCurve ** crvarr, double start, double stop, double *et,
       int in, int iordr, double **iright, int *jstat)
#else
void
s1932 (inbcrv, crvarr, start, stop, et, in, iordr, iright, jstat)
     int inbcrv;
     SISLCurve **crvarr;
     double start;
     double stop;
     double *et;
     int in;
     int iordr;
     double **iright;
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To express a set of B-spline curves described in different
*		B-spline bases using a basis reflecting the continuity in
*		all the B-spline bases used in the description of the
*		different curves.
*
*
* INPUT      :	inbcrv	- Number of curves in the curve-set.
*		crvarr	- Array (length inbcrv) of pointers to the
*			  curves in the curve-set.
*		start	- Start parameter-value of B-spline basis to be made.
*		stop	- Stop parameter-value of B-spline basis to be made.
*		et	- The knot vector of the basis in which the curves
*			  are to be expressed.
*		in	- The number of B-spline bases in which the curves
*			  are to be expressed.
*		iordr	- The order of the basis in which the curves are to
*			  be expressed.
*
* OUTPUT     :  iright	- Array containing the right hand side of the equation
*			  system. Sequence first curve, second curve etc...
*               jstat   - Output status:
*                          < 0: Error.
*                          = 0: Ok.
*                          > 0: Warning.
*
* METHOD     : 	The description of the curves are lifted to the order given
*		as input. The lifted knot vectors are then mapped into the
*		parameter range given by astart and astop. Then the curves
*		are expressed in the refined basis given et. This basis must
*		be calculated to be refinement of the lifted and mapped knot
*		vectors of the input curves.
*
* REFERENCES :  Fortran revised version:
*               T.Dokken, SI, 1989-02
*
* CALLS      :  s1750,s1934,s1936,s6err.
*
*
* WRITTEN BY :  Christophe R. Birkeland, SI, 1991-07
*
*********************************************************************
*/
{
  int ki, kj, kp;		/* Loop control parameters 		*/
  int kpos = 0;			/* Error position indicator		*/
  int kordr;			/* Highest curve-order used in
				 * description				*/
  int idim;			/* Dimension of the space where the
				 * curves lie				*/
  int kstat;			/* Status variable from lower level
				 * routines				*/
  SISLCurve *crv = SISL_NULL;	/* SISLCurve, sent to s1936		*/
  double *kdumcf = SISL_NULL;	/* Contains curve-coefficients sent over
				 * from s1936				*/


  *jstat = 0;

  /* Initailzation of variables */

  idim = crvarr[0]->idim;
  kordr = 0;

  for (ki = 0; ki < inbcrv; ki++)
    if ((crvarr[ki]->in <crvarr[ki]->ik) ||(crvarr[ki]->ik < 1))
      goto err112;


  /* Find highest order used in description */

  for (ki = 0; ki < inbcrv; ki++)
    kordr = MAX (kordr, crvarr[ki]->ik);

  if (kordr > iordr)
    goto err151;


  /* Lift the order of the description of the curve, then refine the
     description of the curve, and copy the result into the right hand
     side of the equation ssystem. */

  /* Allocate array kdumcf and output array iright */

  kdumcf = newarray (idim * in, DOUBLE);
  if (kdumcf == SISL_NULL)
    goto err101;
  *iright = newarray (in *idim * inbcrv, DOUBLE);
  if (*iright == SISL_NULL)
    goto err101;

  kp = 0;

  for (ki = 0; ki < inbcrv; ki++)
    {
      /* Increase the order of the basis */

      s1750 (crvarr[ki], iordr, &crv, &kstat);
      if (kstat < 0)
	goto error;


      /* Normalize the lifted basis */

      s1934 (crv->et, crv->in, crv->ik, start, stop, &kstat);
      if (kstat < 0)
	goto error;


      /* Express the curve using the refined basis */

      s1936 (crv, et, in, kdumcf, &kstat);
      if (kstat < 0)
	goto error;

      if (crv != SISL_NULL)
	freeCurve (crv);
      crv = SISL_NULL;	


      /* Copy coefficients into right hand side of equation system */

      for (kj = 0; kj < (in *idim); kj++)
	(*iright)[kj + kp] = kdumcf[kj];
      kp += in *idim;
    }
  goto out;

  /* Memory error */

  err101:
    *jstat = -101;
    s6err ("s1932", *jstat, kpos);
    goto out;

  /* Error in curve description */ 
  
  err112:
     *jstat = -112;
     s6err("s1932",*jstat,kpos);
     goto out;
 
  /* Error in input */

  err151:  
    *jstat = -151;
    s6err ("s1932", *jstat, kpos);
    goto out;

  /* Error in lower level routines */

  error:
    *jstat = kstat;
    s6err ("s1932", *jstat, kpos);
    goto out;

  out:
   /* Free scratch occupied by local array.  */
   
   if (kdumcf != SISL_NULL) freearray(kdumcf);

   return;
}
