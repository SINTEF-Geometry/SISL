/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright                                                             */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s1225.c,v 1.3 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1225

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1225(SISLCurve *curve,int der,double parvalue,int *leftknot,
      double derive[],double curvature[], double *radius_of_curvature,
      int *jstat)
#else
void s1225(curve,der,parvalue,leftknot,derive,curvature,
	   radius_of_curvature,jstat)
     SISLCurve  *curve;
     int    der;
     double parvalue;
     int    *leftknot;
     double derive[];
     double curvature[];
     double *radius_of_curvature;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate position, first derivative, curvature and radius of
*              curvature of a curve at a given parameter value, from the
*              left hand side.
*
*
*
* INPUT      : curve    - Pointer to the curve for which position
*                       and derivatives are to be computed.
*              der      - The number of derivatives to compute.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*              parvalue - The parameter value at which to compute
*                       position and derivatives.
*
*
*
* INPUT/OUTPUT : leftknot - Pointer to the interval in the knot vector
*                        where ax is located. If et is the knot vector,
*                        the relation
*
*                          et[ileft] < parvalue <= et[ileft+1]
*
*                        should hold. (If parvalue == et[ik-1] then ileft
*                        should be ik-1. Here in is the number of B-spline
*                        coefficients.)
*                        If ileft does not have the right value upon
*                        entry to the routine, its value will be changed
*                        to the value satisfying the above condition.
*
*
*
* OUTPUT     : derive   - Double array of dimension [(ider+1)*idim]
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
*              curvature - Array of dimension idim
*              radius_of_curvature -  A radius of curvature =-1, indicates
*                       that the radius of curvature is infinit.
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : The derivatives are evaluated from the left hand
*              side by  s1227
*              The curvature and the radius of curvature are evaluated
*              by s1307
*
* REFERENCES :
*
*-
* CALLS      : s1227, s1307
*
* WRITTEN BY : Cathrine Tegnander
* MODIFIED BY :
* REVISED BY :
*
*********************************************************************
*/
{
  int ider=2;                 /* minimum number of derivatives needed to
                                 calculated the curvature */
  int kdim = curve -> idim;   /* copy curve attribute to local parameter */
  int kstat=0;                /* local status variable */
  int kpos = 0;               /* local error position  */
  int iknot = 0;              /* local version of leftknot */
  double *iderive = SISL_NULL;     /* pointer to array used to store the position
                                 and the first and second derivatives  */
  double *egeo = SISL_NULL;        /* pointer to store curvature and radius of
                                 curvature */

  iderive= newarray(3*kdim,DOUBLE);
  if (iderive == SISL_NULL)
    goto err101;
  egeo= newarray(3*kdim+1,DOUBLE);
  if (egeo == SISL_NULL)
    goto err101;



  /* Evaluate the derivatives */
    if (der<2)
      {
	s1227(curve,ider,parvalue,&iknot,iderive,&kstat);
	if (kstat<0) goto error;

	/* Copy position */
	memcopy(derive,iderive,(der+1)*kdim,DOUBLE);
      }
  else
    {
      s1227(curve,der,parvalue,&iknot,derive,&kstat);
      if (kstat<0) goto error;

      	/* Copy position */
	memcopy(iderive,derive,3*kdim,DOUBLE);
    }
    *leftknot = iknot;

  /* Evaluate the curvature and the radius_of_curvature */
  s1307(iderive,kdim,egeo,&kstat);
  if (kstat<0) goto error;

  /* Copy position */
  memcopy(curvature,egeo+kdim*2,kdim,DOUBLE);

  *radius_of_curvature = egeo[kdim*3];

  /* Free memory */
 freearray(iderive);
 freearray(egeo);

 /* Successful computations.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation */
 err101:
  *jstat = -101;
  s6err("s1226",*jstat,kpos);


  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err("S1227",*jstat,kpos);
  goto out;

 out: return;

}
