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

#define S2554

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2554( SISLCurve *curve,
       double     parvalue,
       int       *leftknot,
       double     derive[],
       double    *torsion,
       int       *jstat )
#else
void s2554( curve, parvalue, leftknot, derive, torsion, jstat )
     SISLCurve  *curve;
     double      parvalue;
     int        *leftknot;
     double      derive[];
     double     *torsion;
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the torsion of a curve at a given parameter value,
*              from the right hand side.
*
*
*
* INPUT      : curve    - Pointer to the curve.
*              parvalue - The parameter value at which to compute
*                         torsion.
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
* OUTPUT     : derive   - Double array of dimension [4*idim]
*                       containing the position and derivative vectors.
*                       (idim is the number of components of each B-spline
*                       coefficient, i.e. the dimension of the Euclidean
*                       space in which the curve lies.)
*                       These vectors are stored in the following order:
*                       First the idim components of the position vector,
*                       then the idim components of the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
*              torsion - The torsion value computed
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : The derivatives are evaluated from the right hand
*              side by  s1221(), and the torsion are evaluated by s2555()
*
* REFERENCES :
*
*-
* CALLS      : s1221(), s2555()
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


  /* Check input */

  if ( kdim != 3 )
  {
    *torsion = 0.0;
    goto out;
  }



  /* Evaluate the derivatives */

  s1221( curve, 3, parvalue, leftknot, derive, &kstat );

  if ( kstat < 0 ) goto error;


  /* Evaluate the torsion based on the derivatives in "derive" */

  s2555( derive, torsion, &kstat );

  if ( kstat < 0 ) goto error;



 /* Successful computations.  */

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err( "s2554", *jstat, kpos );
  goto out;


 out:
 return;

}
