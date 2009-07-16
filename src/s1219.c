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
 * $Id: s1219.c,v 1.2 2000-05-23 08:44:05 vsk Exp $
 *
 */


#define S1219

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1219(double *et,int ik,int in,int *ileft,double ax,int *jstat)
#else
void s1219(et,ik,in,ileft,ax,jstat)
     double *et;
     int    ik;
     int    in;
     int    *ileft;
     double ax;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To localize the point ax in the array et.
*              The output ileft should satisfy the relations
*                          
*                    et[ileft] <= ax < et[ileft+1].
* 
*              There are two exceptions to this. If (ax >= et[in])
*              then ileft should be in-1 (this corresponds to extending
*              the polynomial piece between et[in-1] and et[in] to the
*              right of the natural parameter interval.
*              Similarly, if (ax < et[ik-1]) then ileft should still be
*              ik-1.
*
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              in     - The dimension of the spline space associated with
*                       the knot vector et.
*              ax     - The point at which the B-spline values and derivatives
*                       are to be computed.
*
*                
*
* INPUT/OUTPUT : ileft - Pointer to the interval in the knot vector
*                       where ax is located, check the relations above.
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : The aim is to do as little work as possible in the cases
*              where ileft has the right or almost the right value.
*              First of all we make sure that ileft has a legal value
*              (a value in the range ik-1 to in-1). Then we check
*              if the current value is OK.
*              If it is not we check that ax is in the interior of et
*              or if the right value is obtained by either increasing
*              or decreasing ileft by 1. If the right value still has
*              not been found we do a binary search.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Knut Moerken, University of Oslo, August 1988.
*
*********************************************************************
*/                                     
{
  int kpos=0;         /* The position of the error.                      */
  int kleft;          /* Local version of ileft to avoid the pointer.    */
  
  /* Check the input. */
  
  if (ik < 1) goto err110;
  
  if (in < ik) goto err111;
  
  if (et[ik-1] == et[ik] || et[in-1] == et[in]) goto err112;
  
  /* Make sure that kleft is in the legal range. */
  
  kleft = min(max(ik-1,*ileft),in-1);
  
  /* Check if the current value of kleft is acceptable. */
  
  if (et[kleft] <= ax && ax < et[kleft+1]) ;
  
  /* Check if ax is outside (et[ik-1],et[in]). */
  
  else if (ax >= et[in-1])
    kleft = in - 1;
  else if (ax <= et[ik-1])
    kleft = ik - 1;
  
  /* Check if it is sufficient to increase or decrease kleft by one. */
  
  else if (et[kleft+1] <= ax && ax < et[kleft+2])
    kleft += 1;
  else if (kleft > 0 && et[kleft-1] <= ax && ax < et[kleft])
    kleft -= 1;
  
  /* Last resort - a binary search. */
  else
    {
      
      /* kmin and kmax gives the upper and lower limits on the possible values
	 of kleft.                                                       */
      
      int kmin,kmax;
      
      kmin = ik - 1; kmax = in - 1;
      kleft = (kmin+kmax)/2;
      
      while (ax < et[kleft] || et[kleft+1] <= ax)
	{
	  if (ax < et[kleft])
	    kmax = kleft;
	  else
	    kmin = kleft;
	  
	  kleft = (kmin+kmax)/2;
	}
    }
  
  *ileft = kleft;
  
  /* Successful computations.  */
  
  *jstat = 0;
  goto out;
  
  /* Polynomial order less than 1. */
 err110: *jstat = -110;
  s6err("s1219",*jstat,kpos);
  goto out;
  
  /* Fewer B-splines than the order. */
 err111: *jstat = -111;
  s6err("s1219",*jstat,kpos);
  goto out;
  
  /* Error in knot vector.
     (The first or last interval of the knot vector is empty.) */
 err112: *jstat = -112;
  s6err("s1219",*jstat,kpos);
  goto out;
  
 out: return;
}
