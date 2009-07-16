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
 * $Id: s1380.c,v 1.2 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1380

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1380(double ep[],double ev[],int im,int idim,int ipar,
	   SISLCurve **rcurve,int *jstat)
#else
void s1380(ep,ev,im,idim,ipar,rcurve,jstat)
     double ep[];
     double ev[];
     int    im;
     int    idim;
     int    ipar;
     SISLCurve  **rcurve;
     int    *jstat;
#endif
/*
************************************************************************
*
* Purpose: To compute the cubic Hermit interpolant to the data given
*          by the points ep and the derivatives ev. 
*          The curve is represented as a B-spline curve.
*
* Input:
*          ep     - Array containing the point in sequence
*                   (x,y,..,x,y,..)
*          ev     - Array containing the derivatives in sequence
*                   (x,y,..,x,y,..)
*          im     - Number of point and derivatives
*          idim   - The dimension of the space the points and derivatives
*                   lie in
*          ipar   - Type of parametrization
*                    1 - Parametrization using cordlength between point
*                  !=1 - Uniform parametrization
* Output:
*          rcurve - Pointer to the curve produced
*          jstat  - Status variable
*                    < 0 - Error.
* Method:
*     First the parmaterization is calculated and then the interpolation
*     is performed using s1379.
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI 1988-11
*
*********************************************************************
*/                                                               
{
  int ki;             /* Loop variables                              */
  int kpek1,kpek2;    /* Pointers into point array                   */
  int kstat;          /* Status variable                             */
  int kpos=0;         /* Position of error                           */
  double *spar=SISL_NULL;  /* Pointer to parametrization array            */
  
  
  
  /* Check input */        
  
  if (im < 2)   goto err181;
  if (idim < 1) goto err102;
  
  /* Allocate array for parametrization */
  
  spar = newarray(im,DOUBLE);
  if (spar == SISL_NULL) goto err101;
  
  spar[0] = (double)0.0;
  
  if (ipar == 1)                 
    {
      /*  Cord length parametrization */
      
      kpek1 = 0;
      for (ki=1 ; ki<im ; ki++)
	{
	  kpek2 = kpek1 + idim;
	  spar[ki] = spar[ki-1] + s6dist(&ep[kpek2],&ep[kpek1],idim);
	  kpek1 = kpek2;
	}
    }
  else
    {
      /*  Uniform parametrization */
      for (ki=0;ki<im;ki++)
	spar[ki] = ki;
    }
  
  /* Calculate Hermite interpolant */
  
  s1379(ep,ev,spar,im,idim,rcurve,&kstat);
  if (kstat<0) goto error;
  
  /* Calculation completed */
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation. Return zero. */
  
  
  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1380",*jstat,kpos);
  goto out;
  
  
  /* Dimension less than 1*/
 err102: *jstat = -102;
  s6err("s1380",*jstat,kpos);
  goto out;
  
  /* Too few interpolation conditions */
  
 err181: *jstat = -181;
  s6err("s1380",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine */
  
 error:  *jstat = kstat;
  s6err("s1380",*jstat,kpos);
  goto out;
  
  
 out:
  if (spar != SISL_NULL) freearray(spar);
  
  return;
}
