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
 * $Id:
 *
 */

#define S2551

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2551( SISLCurve *curve,
       double     parvalue,
       int       *leftknot,
       double     derive[],
       double    *curvature,
       int       *jstat )
#else
void s2551( curve, parvalue, leftknot, derive, curvature, jstat )
     SISLCurve  *curve;
     double      parvalue;
     int        *leftknot;
     double      derive[];
     double     *curvature;
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the curvature of a curve at a given parameter value,
*              from the right hand side.
*
*
*
* INPUT      : curve    - Pointer to the curve.
*              parvalue - The parameter value at which to compute
*                         curvature.
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
*              curvature - The curvature value computed
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : The derivatives are evaluated from the right hand
*              side by  s1221
*              The curvature are evaluated by s1307
*
* REFERENCES :
*
*-
* CALLS      : s1221, s1307
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int kdim = curve -> idim;   /* copy curve attribute to local parameter  */
  int kstat = 0;              /* local status variable                    */
  int kpos = 0;               /* local error position                     */
  double *egeo = SISL_NULL;        /* pointer to store pos, tangent vector,
				 curvature vector and radius of curvature */



  /* Allocate local arrays */

  egeo = newarray( 3*kdim+1, DOUBLE );

  if ( egeo == SISL_NULL )      goto err101;


  /* Evaluate the derivatives */

  s1221( curve, 2, parvalue, leftknot, derive, &kstat );

  if ( kstat < 0 ) goto error;


  /* Evaluate the curvature vector and the radius_of_curvature */

  s1307( derive, kdim, egeo, &kstat );

  if ( kstat < 0 ) goto error;


  /* Evaluate curvature */

  *curvature = s6length( egeo+2*kdim, kdim, &kstat );




 /* Successful computations.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation */
 err101:
  *jstat = -101;
  s6err( "s2551",*jstat,kpos);

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err( "s2551", *jstat, kpos );
  goto out;


 out:
  /* Free local arrays */

 if ( egeo != SISL_NULL ) freearray( egeo );

 return;

}
