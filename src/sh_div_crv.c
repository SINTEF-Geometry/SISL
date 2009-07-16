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
 * $Id: sh_div_crv.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH_DIV_CRV

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    sh_div_crv (SISLCurve * pc, int which_end, double aepsge, SISLCurve ** rcnew, int *jstat)
#else
void 
   sh_div_crv (pc, which_end, aepsge, rcnew, jstat)
     SISLCurve *pc;
     int which_end;
     double aepsge;
     SISLCurve **rcnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE     :To factorize a bezier curve J(x) over the interval [a,b]
*              into 
*              oldcurve = (x-a)/(b-a)*newcurve 
*              when which_end eq 0   (requires C0 = 0)
*                 and
*              oldcurve = (b-x)/(b-a)*newcurve 
*              when which_end eq 1   (requires Cn = 0).
*
*
*
*
* INPUT      : pc           - Oldcurve to factorize.
*              which_end     - Branch parameter for zero point.
*              aepsge       - Geometry tolerance.
*
*
*
* OUTPUT     : rcnew      -The new curve.
*              jstat     - status messages
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
* CALLS      :
*              
*
* WRITTEN BY : Ulf J. Krystad, SI, 92-12.
* MODIFIED BY :
*
**********************************************************************/
{
  int kpos = 0;			/* Position of error.               */
  int ki,kj;                    /* Loop control                     */
  int kn,kk,kdim;               /* Attributes of inut curve         */			/* Position of error.               */
  double a,b;                   /* Bezier interval                  */
  double *et_new = SISL_NULL;        /* New knot array                   */
  double *ecoef_new = SISL_NULL;     /* New coefficient array            */
  SISLCurve *qc = SISL_NULL;		/* Pointer to new curve-object.     */


  /* Check that we have a curve. */
  if (!pc)
    goto err150;

  /* Minimum order allowed is 3. */
  if (pc->ik < 3)
     goto err151;
  
  /* The curve has to be of bezier type. */
  if (pc->in != pc->ik)
     goto err152;

  kn = pc->in;
  kk = pc->ik;
  a  = pc->et[kk-1];
  b  = pc->et[kn];
  kdim = pc->idim;

    /* Test if the corresponding coeficient is zero. */
/*  if (which_end == 0)
  {
     for (ki=0; ki < kdim; ki++)
	if (fabs(pc->ecoef[ki]) > aepsge)
	   goto err153;
  }
  else
  {
     for (ki=(kn-1)*kdim; ki < kn*kdim; ki++)
	if (fabs(pc->ecoef[ki]) > aepsge)
	   goto err153;
  }
  */
  
  /* create knot array. __________________________________________*/
  if ((et_new= newarray(kn+kk-2,DOUBLE)) == SISL_NULL) goto err101;

  for (ki=0; ki < kk-1; ki++)
  et_new[ki] = a;

    for (; ki < kn+kk-2; ki++)
  et_new[ki] = b;

  /* create coeficient array. _________________________________ */
  if ((ecoef_new= newarray(kdim*(kn-1),DOUBLE)) == SISL_NULL) goto err101;

  if (which_end)
     for (ki=0; ki < kn-1; ki++)
	for (kj=0; kj < kdim; kj++)
	   ecoef_new[ki*kdim +kj] = pc->ecoef[ki*kdim +kj]*(kn-1)/(kn-1-ki);
  else
     for (ki=0; ki < kn-1; ki++)
	for (kj=0; kj < kdim; kj++)
	   ecoef_new[ki*kdim +kj] = pc->ecoef[(ki+1)*kdim + kj]*(kn-1)/(ki+1);
  
  
  /* Create factor curve */
  if ((qc = newCurve (kn-1, kk-1, et_new, ecoef_new, pc->ikind, kdim, 2))
      == SISL_NULL) goto err101;

  *rcnew = qc;
  *jstat = 0;
  goto out;

/* ERROR EXITS ___________________________________________ */

/* Error. No input curve.  */
err150:
  *jstat = -150;
  s6err ("sh_div_crv", *jstat, kpos);
  goto out;


/* Error. order less than 3.  */
err151:
  *jstat = -151;
  s6err ("sh_div_crv", *jstat, kpos);
  goto out;

/* Error. Not a bezier curve.  */
err152:
  *jstat = -152;
  s6err ("sh_div_crv", *jstat, kpos);
  goto out;


/* Error in allocation.*/

err101:
  if (et_new) freearray(et_new);
  if (ecoef_new) freearray(ecoef_new);
  *jstat = -101;
  s6err ("sh_div_crv", *jstat, kpos);
  goto out;

out:
;
}
