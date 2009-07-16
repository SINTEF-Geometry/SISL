/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

#define SFNDINTVL

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s6fndintvl(double *et,int ik,int in,int *ileft,
	      double ax1,double ax2,int mu_max,int *jstat)
#else
void s6fndintvl(et,ik,in,ileft,ax1,ax2,mu_max,jstat)
     double *et;
     int    ik;
     int    in;
     int    *ileft;
     double ax1;
     double ax2;
     int mu_max;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To examine if two parameter values are separated by a knot
*              of multiplicity more than mu_max including boarders.
*                          
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              in     - The dimension of the spline space associated with
*                       the knot vector et.
*              ax1    - First parameter value
*              ax2    - Second parameter value
*              mu_max - Maximum allowed multiplicity.
*
*                
*
* INPUT/OUTPUT : ileft - Pointer to the interval in the knot vector
*                       where the first separating knot is located.
*                *jstat - Status:
*                         = 0 : No separating knot found.
*                         = 1 : Separating knot found.
*     			  < 0 : Error.
*
*
* METHOD     :
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SINTEF Oslo, 18.07.93.
*
*********************************************************************
*/                                     
{
  int kpos=0;           /* The position of the error.                      */
  int kstat=0;          /* Local status                                    */
  int kleft_1=*ileft;   /* Local version of ileft to avoid the pointer.    */
  int kleft_2=*ileft;   /* Local version of ileft to avoid the pointer.    */
  int mu = 0;           /* Knot mltiplicity.                               */
  double tmp;
  double tval;
  /* _____________________________________________________________________ */
  *jstat = 0;

  /* Sort position */
  if (ax1 > ax2)
  {
     tmp = ax1;
     ax1 = ax2;
     ax2 = tmp;
  }
  
  
  /* Find knot navigators */
  s1219(et,ik,in,&kleft_1,ax1,&kstat);
  if (kstat < 0) goto error;
  
  tval = et[kleft_1+1];
  while (tval < ax2 && tval < et[in])
  {
     mu = s6knotmult(et,ik,in,&kleft_2,tval, &kstat);
     if (mu > mu_max)
     {
	*jstat = 1;
	*ileft = kleft_2;
	break;
     }
     tval = et[kleft_2 +1];
  }
  
  /* Successful computations.  */
  goto out;
  
  
  /* Error */
 error: *jstat = kstat;
  s6err("s6fndintvl",*jstat,kpos);
  goto out;
  
 out: return;
}
