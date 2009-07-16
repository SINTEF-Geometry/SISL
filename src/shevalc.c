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
 * $Id: shevalc.c,v 1.3 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SHEVALC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
shevalc(SISLCurve *pc1,int ider,double ax,double aepsge,int *ileft,
	     double eder[],int *jstat)
#else
void shevalc(pc1,ider,ax,aepsge,ileft,eder,jstat)
     SISLCurve *pc1;
     int ider;
     double ax;
     double aepsge;
     int *ileft;
     double eder[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To compute the value and ider first derivatives of the
*              B-spline curve pointed to by pc1, at the point with
*              parameter value ax. Use the filtered coefficients of the
*              curve.
*
*
*
* INPUT      : pc1    - Pointer to the curve for which position
*                       and derivatives are to be computed.
*              ider   - The number of derivatives to compute.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*              ax     - The parameter value at which to compute
*                       position and derivatives.
*              aepsge - Geometry resolution.
*
*
*
* INPUT/OUTPUT : ileft - Pointer to the interval in the knot vector
*                        where ax is located. If et is the knot vector,
*                        the relation
*
*                          et[ileft] <= ax < et[ileft+1]
*
*                        should hold. (If ax == et[in] then ileft should
*                        be in-1. Here in is the number of B-spline
*                        coefficients.)
*                        If ileft does not have the right value upon
*                        entry to the routine, its value will be changed
*                        to the value satisfying the above condition.
*
*
*
* OUTPUT     : eder   - Double array of dimension [(ider+1)*idim]
*                       containing the position and derivative vectors.
*                       (idim is the number of components of each B-spline
*                       coefficient, i.e. the dimension of the Euclidean
*                       space in which the curve lies.)
*                       These vectors are stored in the following order:
*                       First the idim components of the position vector,
*                       then the idim components of the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
*                       (The C declaration of eder as a two dimensional array
*                       would therefore be eder[ider+1,idim].)
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     :
*
* REFERENCES :
*
*-
* CALLS      : s1221    - Evaluate curve.
*              s1991    - Make the direction cone of a curve.
*              newCurve - Create new curve object.
*              freeCurve - Free scratch occupied by curve object.
*
* WRITTEN BY :  Vibeke Skytt, SI, 04.91.
* CORRECTED BY: UJK, SI, 06.91
* Modified by : Paal Fugelli, SINTEF, Oslo, Norway, 09.94. Replaced
*               code using 'pdir->esmooth' and s1991() with 'qc = pc1;'
*               and modified the free'ing of 'qc'.
*               Added check for 1D rationals.
*********************************************************************
*/
{
  int kstat=0;        /* Local status variable.                          */
  int kdim = pc1->idim;  /* Dimension of geometry space.                 */
  double *scoef=SISL_NULL;    /* Array storing filtered coefficients.         */
  double *s1,*s2,*s3,*s4; /* Pointers into coefficient arrays.           */
  SISLCurve *qc = SISL_NULL;   /* Curve to evaluate.                          */

  /* Make sure that the filtered coefficients of the curve exist.  */

  if (kdim == 1)
  {

    /*
     * PFU 09-94.
     * There should never be a rational 1D curve here according to UJK, but
     * I (PFU) added a test just in case...
     * A rational curve would have caused a memory usage error in newCurve
     * when trying to divide out the weights from the coefs.
     * This could result in a core dump (division by zero) since the data
     * would be "garbage".
     *
     * If future changes requires this to handle rational 1D curves, this
     * must be updated to use rcoef when input is rational.
     *
     */

    if ( pc1->ikind == 2 || pc1->ikind == 4 )
      goto err151;

     /* Create filtered coefficients. */

     if ((scoef = newarray(pc1->in,DOUBLE)) == SISL_NULL) goto err101;

     for (s1=pc1->ecoef, s2=scoef, s3=s1+pc1->in; s1<s3; s1=s4)
     {
	*s2 = *s1;
	for (s2++, s4=s1+1; s4<s3; s4++, s2++)
	{
	   if (fabs((*s4)-(*s1)) < aepsge) *s2 = *s1;
	   else break;
	}
     }

     /* Create curve object.  */

     if ((qc = newCurve(pc1->in,pc1->ik,pc1->et,scoef,pc1->ikind,
			kdim,0)) == SISL_NULL) goto err101;
  }
  else
    qc = pc1;

  /*
   * This previously used AN ANACRONISM ('pdir->esmooth') - taken out
   * (Confirmed by VSK).
   */



  /* Evaluate curve.  */

  s1221(qc,ider,ax,ileft,eder,&kstat);
  if (kstat < 0) goto error;

  /* UJK Let's have a normal exit possibility !*/
  *jstat = 0;
  goto out;


  /* Error in input (1D rationals is not handled) */
 err151:
  *jstat = -151;
  goto out;

  /* Error in scratch allocation.  */
 err101:
  *jstat = -101;
  goto out;

  /* Error in lower level routine.  */

 error:
  *jstat = kstat;
  goto out;

out:
   /* Free scratch occupied by local objects. */

   if (scoef != SISL_NULL) freearray(scoef);
   if (qc != SISL_NULL && qc != pc1 ) freeCurve(qc);

   return;
}
