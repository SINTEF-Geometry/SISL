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
* REVISED BY:   Vibeke Skytt, SINTEF, Dec. 2006. Allow dimension different from 3
*                                                and rational curves
*
*********************************************************************
*/
{
   int m = 0, m1 = 0;
   int my, my1;
   int np, i, j;
   int ik;
   double bas;
   double* ecoef = SISL_NULL;
   int kdim = pc1->idim;
   double scratch[4];
   double *xyz = NULL;
   int krat = (pc1->ikind == 2 || pc1->ikind == 4);


   /* Check the input. */

   //if ( pc1->idim != 3 ) goto err104;

   if (krat)
       kdim++;

   if (kdim > 4)
   {
       if ((xyz = newarray(kdim, DOUBLE)) == SISL_NULL)
	   goto err101;
   }
   else xyz = scratch;


   /* Set input to local variables. */

   ik    = pc1 -> ik;
   ecoef = (krat) ? pc1->rcoef : pc1 -> ecoef;


   for ( np = 0; np < npol; np++ )
   {
     my  = ileft[ np ] - ik;
     
     for (j=0; j<kdim; j++)
	 xyz[j] = 0.0;

     for ( i = 0; i < ik; i++ )
     {
       my++;
       my1 = kdim*my;
       bas = ebder[ m1++ ];

       for (j=0; j<kdim; j++)
	   xyz[j] += ecoef[my1+j]*bas;
     }

     if (krat)
     {
	 for (j=0; j<pc1->idim; j++)
	     xyz[j] /= xyz[pc1->idim];
     }

     for (j=0; j<pc1->idim; j++)
	 eder[m++] = xyz[j];
   }


  /* Successful computations.  */

   *jstat = 0;
   goto out;

/* Error in input, crv->idim != 3 */
 err104: *jstat = -104;
         s6err( "s1541", *jstat, 0 );
         goto out;

 err101: *jstat = -101;
         s6err( "s1541", *jstat, 0 );
         goto out;

out: 
	 if (xyz != SISL_NULL && xyz != scratch)
	     freearray(xyz);
	 return;

}
