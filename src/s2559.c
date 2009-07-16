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

#define S2559

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2559( SISLCurve *curve,
       double     ax[],
       int        num_ax,
       double     p[],
       double     t[],
       double     n[],
       double     b[],
       int       *jstat )
#else
void s2559( curve, ax, num_ax, p, t, n, b, jstat )
     SISLCurve  *curve;
     double      ax[];
     int         num_ax;
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
*              at given parameter values ax[ 0 ],...,ax[ num_ax - 1 ].
*
*
*
* INPUT      : curve    - Pointer to the curve.
*              ax       - The parameter values
*          num_ax       - No. of parameter values
*
*
*
*
* OUTPUT     :
*         (t,n,b) - The Frenet Frame (in 3D) computed. Each of the arrays
*                   (t,n,b) are of dim. 3*num_ax, and the data are
*                   stored like this: tx(ax[0]), ty(ax[0]), tz(ax[0]),
*                   ...,tx(ax[num_ax-1]), ty(ax[num_ax-1]), tz(ax[num_ax-1]).
*
*         p     - 3D curve posistions at  ax[ 0 ],...,ax[ num_ax - 1 ]
*
*        jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : Call the SISL rutine s2560() num_ax times
*
* REFERENCES :
*
*
* CALLS      : s2560()
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int i, m;                   /* Index variable                      */
  int kstat = 0;              /* local status variable               */
  int kpos = 0;               /* local error position                */
  int leftknot = 0;           /* Knot index                          */
  double *derive = SISL_NULL;      /* Derivatives                         */



  /* Allocate local arrays */

  derive = newarray( 3*curve->idim, DOUBLE );

  if ( derive == SISL_NULL )      goto err101;


  /* Evaluate the Frenet Frame in all ax[i] positions. */

  m = 0;
  for ( i = 0; i < num_ax; i++ )
  {
    s2560( curve, ax[i], &leftknot, derive,
	   p + m, t + m, n + m, b + m, &kstat );

    m = m + 3;

    if ( kstat < 0 ) goto error;
  }




 /* Successful computations.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation */
 err101:
  *jstat = -101;
  s6err( "s2559", *jstat, kpos );

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err( "s2559", *jstat, kpos );
  goto out;


 out:
  /* Free local arrays */

 if ( derive != SISL_NULL ) freearray( derive );

 return;

}
