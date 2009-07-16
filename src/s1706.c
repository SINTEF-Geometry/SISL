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
 * $Id: s1706.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1706

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1706(SISLCurve *pc)
#else
void s1706(pc)
     SISLCurve *pc;
#endif
/*
*******************************************************************
*
*********************************************************************
*
* PURPOSE    : Turn the direction of a curve.
*              The start of the new parameter is the same as the start
*              of the old parameter.
*              This rutine turn the direction of the orginal curve.
*              If you want a copy with a turned direction, just
*              make a copy and turn the direction of the copy.
*
*
*
* INPUT      : pc      -The curve.
*
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* REVISED BY : Johannes Kaasa, SI, Sep 1991 (Introduced NURBS).
*
********************************************************************/
{
  int  kk=pc->ik;             /* Order of the input curve.             */
  int  kn=pc->in;             /* Number of vertices in the input curve.*/
  int  kdim=pc->idim;         /* Dimensjon of the space in whice curve
				 lies.                                 */
  register double *s1,*s2;
  register double *s3; 	       /* Pointers used in loop.               */
  register double t1,t2;       /* Help variables.                      */
  
  /* Now curve to turn. */
  
  if (!pc) goto out;
  
  /* Here we are turning the knot vector such that the first
     element have the same value as the old first element. */
  
  for (s1=pc->et,s2=s1+kk+kn-1,t1=(*s1)+(*s2); s1<=s2; s1++,s2--)
    {
      t2 = *s1;
      *s1 = t1 - *s2;
      *s2 = t1 - t2;
    }
  
  /* Here we just turn the vertices. */
  
  for (s1=pc->ecoef,s2=s1+kdim*(kn-1); s1<s2; s2-=2*kdim)
    for (s3=s1+kdim; s1<s3; s1++,s2++)
      {
	t1 = *s1;
	*s1 = *s2;
	*s2 = t1;
      }

  /* If necessary turn rational vertices. */

  if (pc->ikind == 2 || pc->ikind == 4)
    {
      kdim++;
      for (s1=pc->rcoef,s2=s1+kdim*(kn-1); s1<s2; s2-=2*kdim)
        for (s3=s1+kdim; s1<s3; s1++,s2++)
          {
	    t1 = *s1;
            *s1 = *s2;
	    *s2 = t1;
          }
    }
  
 out:
  return;
}

