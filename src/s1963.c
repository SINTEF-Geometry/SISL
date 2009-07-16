/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

#define S1963

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1963(SISLCurve *pc,double eeps[],int ilend,int irend,int iopen, 
	   int itmax, SISLCurve **rc,int *jstat)
#else
void s1963(pc,eeps,ilend,irend,iopen,itmax,rc,jstat)
     SISLCurve  *pc;
     double eeps[];
     int    ilend;
     int    irend;
     int    iopen;
     int    itmax;
     SISLCurve  **rc;
     int    *jstat;
#endif
/*
***********************************************************
*
* Purpose: To approximate the input spline-curve by a cubic spline-
*          curve with error less than eeps in each of the kdim components.
*
* Input :
*         Pc     - Pointer to curve.
*         Eeps   - Array (length kdim) giving the desired accuracy of
*                  the spline-approximation in each component.
*         Ilend  - The no. of derivatives that are not allowed to change
*                  at the left end of the curve.
*                  The (0 - (ilend-1)) derivatives will be kept fixed.
*                  If ilend <0, this routine will set it to 0.
*                  If ilend<kord, this routine will set it to kord.
*         Irend  - The no. of derivatives that are not allowed to change
*                  at the right end of the curve.
*                  The (0 - (irend-1)) derivatives will be kept fixed.
*                  If irend <0, this routine will set it to 0.
*                  If irend<kord, this routine will set it to kord.
*         iopen  - Open/closed parameter.
*                        =  1 : Produce open curve.
*                        =  0 : Produce closed, non-periodic curve if possible.
*                        = -1 : Produce closed, periodic curve if possible.
*         Itmax  - Max. no. of iterations.
*
* Output:
*         Jstat  - Output status:
*                   < 0 : Error.
*                   = 0 : Ok.
*                   > 0 : Warning:
*         Rc     - Pointer to curve.
*-
* Method: First sampling points on the input-curve are determined such
*         that the cubic hermite spline interpolant to these points will
*         deviate with less than half the tolerance from the input-spline,
*         and this hermite spline interpolant is run through data
*         reduction to reduce the no. of knots.
*
* Calls: s1379, s1940, s1355, s6err
*
* Written by: C.R. Birkeland, Si, April 1993.
* Changed and renamed by : Vibeke Skytt, SINTEF Oslo, 02.95. 
*                          Introduced periodicity.
**********************************************************
*/
{
  int i,j;                   /* Loop control variables           */
  int index = 0;
  int idim = pc->idim;       /* Space dimension                  */
  int stat = 0;              /* Error control parameters         */
  int kpos = 0;
  int km;
  int leftknot = 0;
  double *error1 = SISL_NULL;
  double *error2 = SISL_NULL;
  double *epar = SISL_NULL;
  double *derive = SISL_NULL;
  double *kp = SISL_NULL;
  double *kder = SISL_NULL;
  SISLCurve *ocurve = SISL_NULL;
  
  
  /* Check input-curve. */

  if (!pc) goto err150;
  s1707(pc, &stat);
  if (stat<0) goto error;
  
  /* Check input */

  if (pc->in < pc->ik || pc->ik < 1 || idim < 1) goto err103;
  if (pc->ik < 5)
    { 
      /* Increase order of curve to 4 (cubic) */

      s1750(pc, 4, &ocurve, &stat);
      if (stat < 0) goto error;

      /* Curve is now a cubic spline 
       * Call reduction routine      */

      if( (error2 = newarray( idim, DOUBLE )) == SISL_NULL) goto err101;
      s1940( ocurve, eeps, ilend, irend, iopen, itmax, rc, 
	    error2, &stat);
      if (stat < 0) goto error;

      /* Success !  Go out ! */

      goto out;
    }

  /* Set local tolerance */

  if( (error1 = newarray(idim, DOUBLE)) == SISL_NULL) goto err101;

  for (i=0; i<idim; i++)
    error1[i] = 0.5*eeps[i];

  /* Determine the sample points */

  s1355( pc, error1, &epar, &km, &stat );
  if (stat < 0) goto error;

  /* Compute positions and derivatives of input spline
   * at the given sampling points                      */

  derive = newarray( idim * 2, DOUBLE );
  kp     = newarray( idim * km, DOUBLE );
  kder   = newarray( idim * km, DOUBLE );
  if (derive == SISL_NULL || kp == SISL_NULL || kder == SISL_NULL) goto err101;

  for(i=0; i<km; i++)
    {
      if(i<km-1 && epar[i] == epar[i+1])
	s1227( pc, 1, epar[i], &leftknot, derive, &stat );
      else
	s1221( pc, 1, epar[i], &leftknot, derive, &stat );

      for(j=0; j<idim; j++, index++)
	{
	  kp[index] = derive[j];
	  kder[index] = derive[j+idim];
	}
    }
  if (stat < 0) goto error;

  /* Compute Hermite interpolant */

  s1379( kp, kder, epar, km, idim, &ocurve, &stat );
  if (stat < 0) goto error;

  /* Compute datareduction on the cubic hermite interpolant */

  if( (error2 = newarray( idim, DOUBLE )) == SISL_NULL) goto err101;
  s1940( ocurve, error1, ilend, irend, iopen, itmax, rc, 
	error2, &stat);
  if (stat < 0) goto error;

  /* Success ! */
  
  *jstat = 0;
  goto out;
  
  /* Allocation error. */

 err101: 
   *jstat = -101;
   s6err("s1963",*jstat,kpos);
   goto out;
  
  /* Error in input */

  err103: 
    *jstat = -103; 
    s6err("s1963",*jstat,kpos);
    goto out;
  
  /* Empty curve. */

  err150: 
    *jstat = -150;
    s6err("s1963",*jstat,kpos);
    goto out;
  
  /* Error in lower level routine. */

  error: 
    *jstat = stat;
    s6err("s1963",*jstat,kpos);
    goto out;

  /* Exit */

  out:
    /* Free allocated arrays */

    if( error1 != SISL_NULL) freearray(error1);
    if( error2 != SISL_NULL) freearray(error2);
    if( epar   != SISL_NULL) freearray(epar);
    if( derive != SISL_NULL) freearray(derive);
    if( kp     != SISL_NULL) freearray(kp);
    if( kder   != SISL_NULL) freearray(kder);

    /* Free local SISL-curves */

    if( ocurve != SISL_NULL) freeCurve(ocurve);

    return;
}
