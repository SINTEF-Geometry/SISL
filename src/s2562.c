/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of SISL.
 *
 * SISL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * SISL is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with SISL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using SISL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the SISL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "sisl-copyright.h"

/*
 *
 * $Id:
 *
 */

#define S2562

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2562( SISLCurve *curve,
       double     ax[],
       int        num_ax,
       int        val_flag,
       double     p[],
       double     t[],
       double     n[],
       double     b[],
       double     val[],
       int       *jstat )
#else
void s2562( curve, ax, num_ax, val_flag, p, t, n, b, val, jstat )
     SISLCurve  *curve;
     double      ax[];
     int         num_ax;
     int         val_flag;
     double      p[];
     double      t[];
     double      n[];
     double      b[];
     double     val[];
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the 3D position, the Frenet Frame (t,n,b) and
*              geometric property (curvature, torsion or variation of
*              curvature) of a curve at given parameter values
*              ax[0],...,ax[num_ax-1].
*              These data are needed to produce spike plots (using the
*              Frenet Frame and the geometric property) and circular
*              tube plots (using circular in the normal plane (t,b),
*              where the radius is equal to the geometric property times
*              a scaling factor for visual effects).
*
*
* INPUT      : curve    - Pointer to the curve.
*              ax       - The parameter values
*          num_ax       - No. of parameter values
*          val_flag     - Compute geometric property
*                                = 1  : curvature
*                                = 2  : torsion
*                                = 3  : variation of curvature
*
*
*
* OUTPUT     :
*         (t,n,b) - The Frenet Frame (in 3D) computed. Each of the arrays
*                   (t,n,b) are of dim. 3*num_ax, and the data are
*                   stored like this: tx(ax[0]), ty(ax[0]), tz(ax[0]),
*                   ...,tx(ax[num_ax-1]), ty(ax[num_ax-1]), tz(ax[num_ax-1]).
*
*         p     - 3D curve positions at ax[0],...,ax[num_ax-1]
*
*         val   -  Geometric property (curvature, torsion or
*                  variation of curvature) of a curve at given parameter
*                  values ax[0],...,ax[num_ax-1].
*
*        jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
** METHOD     : See formulas in given in the reference.
*
*
* REFERENCES :   Westgaard G., Heimann J.:
*                Analysis and Visualization of Geometric Curve Properties,
*                in M. D{\ae}hlen, T. Lyche and L. L. Schumaker (Eds.)
*                Mathematical Methods in CAGD IV, 1998.
*
*
* CALLS      : s1221(), s2561(), s1307(), s6length(), s2555(), s2558().
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int i, m;                 /* Index variable                      */
  int kstat = 0;            /* Local status variable               */
  int kpos = 0;             /* Local error position                */
  int leftknot = 0;         /* Knot index                          */
  int kdim = curve -> idim; /* Dim. of the space in which the
			       curve lies                          */
  double *derive = SISL_NULL;    /* Derivatives                         */
  double *egeo = SISL_NULL;      /* Array to store pos, tangent vector,
                               curvature vector and radius of
			       curvature, used by s1307()          */




  /* Allocate local arrays */

  derive = newarray( 4*kdim, DOUBLE );

  if ( derive == SISL_NULL )      goto err101;

  egeo = newarray( 3*kdim+1, DOUBLE );

  if ( egeo == SISL_NULL )        goto err101;


  /* Evaluate in all ax[i] positions. */

  m = 0;
  for ( i = 0; i < num_ax; i++ )
  {

    /* Evaluate the derivatives */

    s1221( curve, 3, ax[i], &leftknot, derive, &kstat );

    if ( kstat < 0 ) goto error;


    /* Evaluate 3D position (p) and the 3D Frenet Frame (t,n,b) */

    s2561( derive, kdim, p + m, t + m, n + m, b + m, &kstat );

    if ( kstat < 0 ) goto error;

    m = m + 3;


    /* Evaluate the geometric property */

    if ( val_flag == 1 ) /* Curvature */
      {
	s1307( derive, kdim, egeo, &kstat );

	if ( kstat < 0 ) goto error;

	val[ i ] = s6length( egeo+2*kdim, kdim, &kstat );
      }
    else if ( val_flag == 2 ) /* Torsion */
      {
	s2555( derive, val + i, &kstat );

	if ( kstat < 0 ) goto error;
      }
    else if ( val_flag == 3 )  /* Variation of curvature */
      {
	s2558( derive, kdim, val + i, &kstat );

	if ( kstat < 0 ) goto error;
      }
  }




 /* Successful computations.  */
  *jstat = 0;
  goto out;

  /* Error in scratch allocation */
 err101:
  *jstat = -101;
  s6err( "s2562", *jstat, kpos );


  /* Error in lower level routine.  */
 error:  *jstat = kstat;
  s6err( "s2562", *jstat, kpos );
  goto out;


 out:
  /* Free local arrays */
 if ( derive != SISL_NULL ) freearray( derive );
 if ( egeo   != SISL_NULL ) freearray( egeo   );

 return;

}
