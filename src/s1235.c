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
 * $Id: s1235.c,v 1.3 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1235

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1235(double et[],int in,int ik,int *jnbreak,double **gbreak,int *jstat)
#else
void s1235(et,in,ik,jnbreak,gbreak,jstat)
     double et[];
     int    in;
     int    ik;
     int    *jnbreak;
     double **gbreak;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find break points in a knot vector. The first and last
*              parameter values are break values.
*
*
*
* INPUT      : et     - Knot vector to find break points in.
*              in     - Number of vertices of the curve corresponding
*                       to the knot vector.
*              ik     - Order of the curve corresponding to et.
*
*
*
* OUTPUT     : jnbreak - Number of break points found.
*              gbreak  - Array containing parameter values of break points.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The knot vector has a break point at a knot if the 
*              multiplicity of the knot is ik-1.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* REVISED BY : Vibeke Skytt, SINTEF, 9801. Correction in loop counter
*                                          for the periodic case.
*
*********************************************************************
*/
{
  int kpos = 0;   /* Position of error.                    */
  int kj;         /* Counter.                              */
  int kbreak;     /* Current number of break points found. */
  int kmult;      /* Multiplisity of current knot.         */
  double tprev;   /* Value of previous knot.               */
  double *sbreak; /* Pointer into break point array.       */
  double *st;     /* Pointer into knot vector.             */
  
  /* Allocate space for an array that is at least as great as the
     number of break points.                                       */
  
  *gbreak = SISL_NULL;
  if ((*gbreak = newarray(in+2,double)) == SISL_NULL) goto err101;
  
  /* Set local pointer to and counter of break points.  */
  
  sbreak = *gbreak;
  kbreak = 0;
  
  /* Find break point in start of parameter interval and internal breaks. */
  
  tprev = et[ik-1];
  kmult = ik - 1;
  for (st=et+ik,kj=ik; kj<in; st++,kj++)
    {
      
      if (*st == tprev) kmult++;
      else
	{
	  if (kmult >= ik-1)
	    {
	      
	      /* New break point found.  */
	      
	      *(sbreak++) = tprev;
	      kbreak++;
	    }
	  tprev = *st;
	  kmult = 1;
	}
    }
  
  /* Find break point in end of interval.  */
  
  if (et[in] != tprev && kmult >= ik-1)
    {
      
      /* Remember last internal break point.  */
      
      *(sbreak++) = tprev;
      kbreak++;
    }
  *(sbreak++) = et[in];
  kbreak++;
  
  /* Reduce break point array to correct size.  */
  
  if (kbreak < in+2)
    if ((*gbreak = increasearray(*gbreak,kbreak,double)) == SISL_NULL) goto err101;
  
  /* Break points found.  */
  
  *jnbreak = kbreak;
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: 
  *jstat = -101;
  s6err("s1235",*jstat,kpos);
  goto out;
  
 out: 
  return;
}
