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


#define S2555

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2555( double     derive[],
       double    *torsion,
       int       *jstat )
#else
void s2555( derive, torsion, jstat )
     double      derive[];
     double     *torsion;
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the torsion of a curve based on the
*              derivatives in "derive". We assume that the
*              curve lies in IR^3.
*
*
* INPUT      :
*            derive   - Double array of dimension [4*3]
*                       containing the position and derivative vectors.
*                       These vectors are stored in the following order:
*                       First the position vector, then the tangent vector,
*                       then the 3 components of the second derivative
*                       vector, and so on.
*
*
* OUTPUT     : torsion - The torsion value computed
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : See formula in the book of Do Carmo
*              (Differential Geometry of Curves and Surfaces).
*
* REFERENCES : Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*
*
* CALLS       : s6length(), s6scpr()
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int kstat = 0;        /* Local status variable       */
  double crpr[3], a;    /* Temp. variable.             */


  /* Evaluate torsion */

  s6crss( derive + 3, derive + 6, crpr );

  a = s6length( crpr, 3, &kstat );

  if ( a != 0.0 )
    {
      *torsion = ( s6scpr( derive + 9, crpr, 3 ) ) / ( a*a ) ;
    }
  else
    {
      *torsion = 0.0;

      goto war002;
    }



 /* Successful computations.  */

  *jstat = 0;
  goto out;


  /* Torsion undefined, since the curvature = 0.0.  */
war002:
  *jstat = 2;
  goto out;


 out:

 return;

}
