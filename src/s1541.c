/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id:
 *
 */
#define S1541

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1541   ( SISLCurve  *pc1,
	  int         npol,
	  double      ebder[],
	  int         ileft[],
	  double      eder[],
	  int        *jstat )

#else
void
s1541 ( pc1, npol, ebder, ileft, eder, jstat )
     SISLCurve*  pc1;
     int         npol;
     double      ebder[];
     int         ileft[];
     double      eder[];
     int        *jstat;

#endif
/*
*********************************************************************
*
* PURPOSE:	Given a (polynomial) spline curve pc1 and
*		preevaluated basis functions (using s1540()) on an
*               npol polyline, calculate the 3D positions on
*		that polyline (eder).
*
* INPUT:	pc1  	-  the spline curve
*		npol    -  the number of polyline points
*              ebder   -  array containing the basis values
*                           B(ax[0   ],i0-k+1),...,B(ax[0   ],i0)
*                           B(ax[1   ],i1-k+1),...,B(ax[1   ],i1)
*                            :                :
*                           B(ax[m1-1],im1-1-k+1),...,B(ax[m1-1],im1-1)
*
*         	ileft   -  ileft[i] <= ti < ileft[i] + 1
*			   (exception for ti == in : ileft[ti] = n-1)
*
* OUTPUT:	eder	-  contains the 3D polyline points
*              jstat    - Status messages
*
*                         = 0 : Ok.
*                         < 0 : Error.
*
*
* METHOD:
*
* WRITTEN BY:	Geir Westgaard, SINTEF, Oslo, November 1999
*
*********************************************************************
*/
{
   int m = 0, m1 = 0;
   int my, my1;
   int np, i;
   int ik;
   double bas;
   double cx, cy, cz;
   double* ecoef = SISL_NULL;


   /* Check the input. */

   if ( pc1->idim != 3 ) goto err104;


   /* Set input to local variables. */

   ik    = pc1 -> ik;
   ecoef = pc1 -> ecoef;


   for ( np = 0; np < npol; np++ )
   {
     my  = ileft[ np ] - ik;

     cx = cy = cz = 0.0;

     for ( i = 0; i < ik; i++ )
     {
       my++;
       my1 = 3*my;
       bas = ebder[ m1++ ];

       cx += ecoef[ my1   ]*bas;
       cy += ecoef[ my1+1 ]*bas;
       cz += ecoef[ my1+2 ]*bas;
     }
     eder[ m++ ] = cx;
     eder[ m++ ] = cy;
     eder[ m++ ] = cz;
   }


  /* Successful computations.  */

   *jstat = 0;
   goto out;

/* Error in input, crv->idim != 3 */
 err104: *jstat = -104;
         s6err( "s1541", *jstat, 0 );
         goto out;

out: return;

}
