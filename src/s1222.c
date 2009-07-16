/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1995 by                                                     */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s1222.c,v 1.5 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1222

#include "sislP.h"
 


#if defined(SISLNEEDPROTOTYPES)
void 
   s1222(double et[], int ik, int in, int ibase, 
	   double ax, int ider, double ebder[], int *jstat)
#else
void s1222(et, ik, in, ibase, ax, ider, ebder, jstat)
     double et[];
     int    ik;
     int    in;
     int    ibase;
     double ax;
     int    ider;
     double ebder[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the value and ider first derivatives of the
*              B-spline base function starting at et[ibase], associated
*              with the knot vector et at the point ax.
*              REMARK: THE OUTPUT ARRAY MUST HAVE DIMENSION IK*(IDER + 1),
*                      BECAUSE OF INTERNAL USE !!!!!!!!!!!!!!!!!!!!!!
*
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              in     - The dimension of the spline space associated with
*                       the knot vector et.
*              ibase  - The B-spline base function to evaluate (starting at
*                       et[ibase]).
*              ax     - The point at which the B-spline value and derivatives
*                       are to be computed.
*              ider   - The number of derivatives to be computed.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*
*                
*
* OUTPUT     : ebder  - Double array of containing value of the B-splines and 
*                       its derivatives at the point ax. These numbers are 
*                       stored in the following order:
*                       Position, 1. derivative, ... .
*                       REMARK: THIS ARRAY MUST HAVE DIMENSION IK*(IDER + 1),
*                               BECAUSE OF INTERNAL USE !!!!!!!!!!!!!!!!!!!!!!
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : This is based on s1220, bit with all unneccessary calculations
*              removed (we calculate only a band of the recursion triangle).
*
* CALLS      : 
*
* WRITTEN BY : Johannes Kaasa, SINTEF, October 1995.
*
*********************************************************************
*/                                     
{
   int kpos=0;          /* The position of the error.                      */
   int ki, kj, kl;      /* Index in for loop.                              */
   int degree;          /* Degree of the basis function, ik - 1.           */
   int numb_places;     /* Number of places, position and derivatives.     */
   int offset;          /* Basis function offset, kleft - ibase.           */
   int pre_deriv;       /* Columns before derivatives are calculated.      */
   int startA;          /* Start value in for loops.                       */
   int stopA;           /* End value in for loops.                         */
   int indexA;          /* Array index.                                    */
   int indexB;          /* Array index.                                    */
   int kder;            /* Local version of ider. All derivatives of order
			   higher than ik-1 are zero so
			   kder = min(ik-1,ider).                          */
   int kleft;           /* Pointer to the interval in the knot vector
                           where ax is located.                            */
   int idx1;            /* Array index.                                    */
   int lknot, rknot;    /* Support limits in the knot array.               */
   double width;        /* Width of knot support.                          */ 
   double td1, td2;     /* Denominators in the weight ratios.              */
   double tw1, tw2;     /* The weight ratios.                              */
   double ts1, ts2;     /* Variants of td1 and td2 for derivatives.        */

   *jstat = 0;
  
   /* Initiation. */
  
   degree = ik - 1;
  
   /* Check the knot vector. */
  
   if (in < ik || ik < 1 || ibase < 0 || ibase > (in - 1))
      goto err112;
   
   /* Check number of derivatives. */
   
   if (ider < 0)
      goto err178;
  
   /* Find the position of the evaluation point. */
   
   if ((ax < et[ibase] && et[ibase] > et[ik - 1])
       || (ax > et[ibase + ik] && et[ibase + ik] < et[in]))
   {
      stopA = ik*(ider + 1);
      for (ki = 0; ki < stopA; ki++)
	 ebder[ki] = 0.;
      goto out;
   }
   else
   {
      for (kleft = max(ibase, ik - 1); kleft < (ibase + ik); kleft++)
      {
	 if (ax < et[kleft + 1] || kleft == (in - 1))
	    break;
      }
   }
   offset = kleft - ibase;
   
   kder = min(degree, ider);
   numb_places = kder + 1;
   pre_deriv = degree - kder;
   
   /* Start the order iteration for initial positions. */
   
   ebder[0] = 1.;
   
   for (ki = 1; ki < ik; ki++)
   {
      
      startA = min(offset + 1, ki);
      stopA  = max(0, (offset - ik + ki));
      
      /* Compute the first basis function. */
      
      lknot = kleft - startA + 1;
      rknot = lknot + ki;
      width = et[rknot] - et[lknot];
      if (width < REL_PAR_RES)
	 goto err112;
      td2 = (double)1./width;
      tw2 = (et[rknot] - ax)*td2;
      ts2 = ki*td2;
      
      indexB = numb_places*startA;
      indexA = indexB - numb_places;
      
      if (ki < (offset + 1))
      {
	 
	 /* Inside area of calculation. */
	 
	 ebder[indexB] = tw2*ebder[indexA];
	 
	 /* If necessary calculate derivatives. */
	 
	 if (ki > pre_deriv && kder > 0)
	 {
	    for (kl = 1; kl < (ki - pre_deriv + 1); kl++)
	       ebder[indexB + kl] = - ts2*ebder[indexA + kl - 1];
	    
	 }
	 
      }
      
      lknot++;
      rknot++;
      indexB = indexA;
      indexA -= numb_places;
     
      /* Go through the interior basis functions. */
      
      for (kj = (startA - 1); kj > stopA; kj--)
      {
	 
	 width = et[rknot] - et[lknot];
	 if (width < REL_PAR_RES)
	    goto err112;
	 td1 = td2;
	 td2 = (double)1./width;
	 tw1 = (double)1. - tw2;
	 tw2 = (et[rknot] - ax)*td2;
	 
	 /* If necessary calculate derivatives. */
	 
	 if (ki > pre_deriv && kder > 0)
	 {
	    ts1 = ts2;
	    ts2 = ki*td2;
	    for (kl = (ki - pre_deriv); kl > 0; kl--)
	       ebder[indexB + kl] = ts1*ebder[indexB + kl - 1] 
		  - ts2*ebder[indexA + kl - 1];
	    
	 }
	 
	 /* Position. */
	 
	 ebder[indexB] = tw1*ebder[indexB] + tw2*ebder[indexA];
	 
	 indexB = indexA;
	 indexA -= numb_places;
	 lknot++;
	 rknot++;
      }
      
      /* Compute the last basis function. */
      
      if (ki < (ik - offset))
      {
	 
	 /* Inside area of calculation. */
	 
	 tw1 = (double)1. - tw2;
	 
	 /* If necessary calculate derivatives. */
	 
	 if (ki > pre_deriv && kder > 0)
	 {
	    ts1 = ts2;
	    for (kl = (ki - pre_deriv); kl > 0; kl--)
	       ebder[indexB + kl] = ts1*ebder[indexB + kl - 1]; 
	    
	 }
	 
	 /* Position. */
	 
	 ebder[indexB] = tw1*ebder[indexB];
	 
      }
      
   }
   
   /* Copy out the result. */
   
   idx1 = offset*numb_places;
   for (ki = 0; ki < numb_places; ki++)
      ebder[ki] = ebder[idx1 + ki];
   
   stopA = ik*(ider + 1);
   for ( ; ki < stopA; ki++)
      ebder[ki] = 0.;
   
   /* Successful computations.  */
   
   goto out;
   
   /* Error in knot vector (something wrong in kleft calculation, 
      we must have et[kleft] < et[kleft + 1]). */
   
   err112: 
      *jstat = -112;
   s6err("s1222",*jstat,kpos);
   goto out;
   
   /* Illegal derivative requested. */
   
   err178: 
      *jstat = -178;
   s6err("s1222",*jstat,kpos);
   goto out;


   out: 
      return;
}
