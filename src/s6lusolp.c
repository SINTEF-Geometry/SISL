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
 * $Id: s6lusolp.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6LUSOLP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6lusolp(double ea[],double eb[],int nl[],int im,int *jstat)
#else
void s6lusolp(ea,eb,nl,im,jstat)
     double ea[];
     double eb[];
     int    nl[];
     int    im;
     int    *jstat;
#endif
/*
************************************************************************
*
***********************************************************************
*
*   PURPOSE : Solve the equationsystem LU=eb by forth- and back-
*             substitution. L and U are stored in ea.
*
*
*   INPUT   : ea   - The factorized coeffecient-matrix.
*             im   - The number of equations.
*             nl   - Ordering of lines in the matrix.
*
*
*   INPUT/OUTPUT : eb - the right side of the equationsystem and
*                       the found values for the unknowns.
*
*   OUTPUT  : jstat  - Status variable.
*                        < 0 : Error
*                        = 0 : ok
*                        > 0 : Warning
*
*                                                                       
*   METHOD  : Solve on the equation-system LU=eb.
*
*
*   REFERENCES : Cheney & Kincaid : Numerical Mathematics and
*                                   Computing.
*
*-
*   CALLS      :
*
*   WRITTEN BY : Vibeke Skytt, SI, 86-10.
*
************************************************************************
*/
{
  int kpos = 0;      /* Position of error.                             */
  int ki,kj;         /* Counters.                                      */
  double *sx = SISL_NULL; /* Array used to keep solution of equation system
			internally.                                    */
  double tdiv;       /* Dividend in expression.                        */
  
  /* Allocate space for local array.  */
  
  if ((sx = newarray(im,double)) == SISL_NULL) goto err101;
  
  for (ki=0; ki<im-1; ki++)
    {
      /*  Gauss on right side of equation  */
      
      for (kj=ki+1; kj<im; kj++)      
	eb[nl[kj]] -= eb[nl[ki]]*ea[ki+nl[kj]*im];
    }
  
  tdiv = ea[im-1+nl[im-1]*im];
  if (DEQUAL(tdiv,DZERO)) goto warn1;
  sx[im-1] = eb[nl[im-1]]/tdiv;
  
  for (ki=im-2; ki>=0; ki--)
    {
      /*  Backwards substitution.   */
      
      for (kj=ki+1; kj<im; kj++)
	eb[nl[ki]] -= sx[kj]*ea[kj+nl[ki]*im];
      
      tdiv = ea[ki+nl[ki]*im];
      if (DEQUAL(tdiv,DZERO)) goto warn1;
      sx[ki] = eb[nl[ki]]/tdiv;
    }   
  for (ki=0; ki<im; ki++) eb[ki] = sx[ki];
  
  /* Equation system solved.  */
  
  *jstat = 0; 
  goto out;

/* Singular equation system.  */

warn1 : *jstat = 1;
        goto out;

/* Error in space allocation.  */

err101: *jstat = -101;
        s6err("s6lusolp",*jstat,kpos);
        goto out;

out:

/* Free space occupied by local array.  */

if (sx != SISL_NULL) freearray(sx);

return;
}
