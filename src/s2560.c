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

#define S2560

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2560( SISLCurve *curve,
       double     parvalue,
       int       *leftknot,
       double     derive[],
       double     p[],
       double     t[],
       double     n[],
       double     b[],
       int       *jstat )
#else
void s2560( curve, parvalue, leftknot, derive, p, t, n, b, jstat )
     SISLCurve  *curve;
     double      parvalue;
     int        *leftknot;
     double      derive[];
     double      p[];
     double      t[];
     double      n[];
     double      b[];
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the Frenet Frame (t,n,b) of a curve
*              at a given parameter value, from the right hand side.
*
*
*
* INPUT      : curve    - Pointer to the curve.
*              parvalue - The parameter value at which to compute
*                         the Frenet Frame.
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
* OUTPUT     : derive   - Double array of dimension [3*idim]
*                       containing the position and derivative vectors.
*                       (idim is the number of components of each B-spline
*                       coefficient, i.e. the dimension of the Euclidean
*                       space in which the curve lies.)
*                       These vectors are stored in the following order:
*                       First the idim components of the position vector,
*                       then the idim components of the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
*
*         (t,n,b) - The Frenet Frame (in 3D) computed. Each of the vectors
*                   (t,n,b) are of dim. 3, and the data are
*                   stored like this: tx(parvalue), ty(parvalue), tz(parvalue).
*
*         p     - 3D curve posistions at parvalue.
*
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : The derivatives are evaluated from the right hand
*              side by  s1221(), and the Frenet Frame are evaluated by s2561()
*
* REFERENCES :
*
*-
* CALLS      : s1221(), s2561()
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





  /* Evaluate the derivatives */

  s1221( curve, 2, parvalue, leftknot, derive, &kstat );

  if ( kstat < 0 ) goto error;


  /* Evaluate the Frenet Frame based on the derivatives in "derive" */

  s2561( derive, kdim, p, t, n, b, &kstat );

  if ( kstat < 0 ) goto error;



 /* Successful computations.  */

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err( "s2560", *jstat, kpos );
  goto out;


 out:
 return;

}
