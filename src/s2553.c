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

#define S2553

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2553( SISLCurve *curve,
       double     ax[],
       int        num_ax,
       double     torsion[],
       int       *jstat )
#else
void s2553( curve, ax, num_ax, torsion, jstat )
     SISLCurve  *curve;
     double      ax[];
     int         num_ax;
     double      torsion[];
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the torsion of a curve at given parameter values
*              ax[ 0 ],...,ax[ num_ax - 1 ].
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
*             torsion - The "num_ax" torsion values computed
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : Call the SISL rutine s2554() num_ax times
*
* REFERENCES :
*
*
* CALLS      : s2554()
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int i;                      /* Index variable                      */
  int kstat = 0;              /* local status variable               */
  int kpos = 0;               /* local error position                */
  int leftknot = 0;           /* Knot index                          */
  double *derive = SISL_NULL;      /* Derivatives                         */



  /* Allocate local arrays */

  derive = newarray( 4*curve->idim, DOUBLE );

  if ( derive == SISL_NULL )      goto err101;


  /* Evaluate the torsion in all ax[i] positions. */

  for ( i = 0; i < num_ax; i++ )
  {
    s2554( curve, ax[i], &leftknot, derive, torsion + i, &kstat );

    if ( kstat < 0 ) goto error;
  }




 /* Successful computations.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation */
 err101:
  *jstat = -101;
  s6err( "s2553",*jstat,kpos);

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err( "s2553", *jstat, kpos );
  goto out;


 out:
  /* Free local arrays */

 if ( derive != SISL_NULL ) freearray( derive );

 return;

}
