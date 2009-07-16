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
