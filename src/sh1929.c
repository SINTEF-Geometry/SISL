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
 * $Id: sh1929.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */

#define SH1929

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      sh1929(double etau[],int in,int ik,int imu,double et[],int im,
	     int ij,double eah[],int *jmuprm,int *jnu,int *jstat)
#else
void sh1929(etau,in,ik,imu,et,im,ij,eah,jmuprm,jnu,jstat)
   double etau[];
   int in;
   int ik;
   int imu;
   double et[];
   int im;
   int ij;
   double eah[];
   int *jmuprm;
   int *jnu;
   int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : To compute the nonzero discrete B-splines at ij of order
*              ik, given the knot vector etau and the refined knot vector
*              et. Equivalently to compute the nonzero elements of row ij
*              of the B-spline refinement matrix from the knot vector etau to
*              the knot vector et.
* 
* INPUT      : etau   - Real array of length (in+ik) containing the original 
*                       knot vector.
*              in     - The dimension of the spline space corresponding
*                       to etau.
*	       ik     - The order of the spline space.
*              imu    - Integer with the property: 
*                          etau[imu] <= et[ij] < etau[imu+1]
*              et     - Real array of length (im+ik) containing the refined
*                       knot vector.
*              im     - The dimension of the spline space corresponding to et.
*              ij     - The point at which the value of the discrete B-splines
*                       is to be evaluated, or equivalently, the row of 
*                       the B-spline refinement matrix that is to be computed.
*                       The integer ij is the index of an element of the
*                       knot vector with 1 <= ij <= im. In addition it is
*                       assumed that there is at least one nonzero discrete
*                       B-spline at ij.
*
* 
* OUTPUT     : eah    - Real array of length (ik) containing the value of the
*                       (at most) ik nonzero discrete B-splines at ij, with 
*                       the last nonzero value in eah[ik-1].
*              jmuprm - Integer giving the position of the last nonzero 
*                       discrete B-spline at ij in the etau knot vector.
*              jnu    - Integer giving the number of nonzero discrete
*                       B-splines at ij reduced by one, i.e.
*                       #nonzero B-splines=inu+1. This may not be true near
*                       the beginning and end of etau. In general the nonzero 
*                       discrete B-splines are the ones with index i
*                       satisfying MAX(*jmuprm-*jnu,1) <= i <= MIN(*jmuprm,in).
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*             
* 
* METHOD     : 
*
*
* REFERENCES : Lyche, T. and Moerken, K., Making the Oslo algorithm more
*              efficient, SIAM Journal of numerical analyis, June 1986.
*              
*
* USE        :
*
*-
* CALLS      :   
*
* WRITTEN BY : Vibeke Skytt, SI, 05.92, on the basis of a routine
*              written by Tom Lyche and Knut Moerken, 12.85.
*
*********************************************************************
*/
{ 
   int ki,kp;
   int kih,kil,kiu,kkh;
   int k1 = ik - 1;
   double tbeta,tbeta1;
   double tj;
   double td1,td2;
   double *ssi = SISL_NULL;
   
   /* Allocate scratch for a local array of length (ik-1) to hold the
      new knots among et[ij+1],...,et[ij+ik-1].  */
   
   if ((ssi = newarray(ik-1,DOUBLE)) == SISL_NULL) goto err101;

   /* Determine jmuprm.  */
   
   for (ki=ij+1, *jmuprm=imu; DEQUAL(et[ki],etau[*jmuprm]) && ki<ij+ik;
        ki++, (*jmuprm)--);
   kih = (*jmuprm) + 1;
   *jnu = 0;
   
   /* Determine the new knots (together with jnu) and store them in ssi.  */

   for (kp=1; kp<=k1; kp++)
     {
	if (DEQUAL(et[ij+kp],etau[kih])) kih++;
	else
	{
	   ssi[*jnu] = et[ij+kp];
	   (*jnu)++;
	}
     }
   
   /* Compute the eah-array.
      The discrete B-splines are built up one order at a time starting with
      splines of order one. With *jnu knots the number of nonzero discrete
      B-splines of order ik at ij will be *jnu+1.   */
   
   eah[k1] = (double)1.0;
   for (kp=0; kp<(*jnu); kp++)
     {
	tbeta1 = (double)0.0;
	tj = ssi[kp];
	kih = kp + ik - (*jnu);
	if (kp >= (*jmuprm))
	   tbeta1 = (tj-etau[0])*eah[ik-(*jmuprm)-1]/(etau[kih]-etau[0]);
	kil = (*jmuprm) - kp;
	kiu = in + (*jnu) - kp;
	if (kil < 1 ) kil = 1;
	if (kiu > (*jmuprm)) kiu = (*jmuprm);
	for (ki=kil; ki<=kiu; ki++)
	  {
	     td1 = tj - etau[ki];
	     td2 = etau[ki+kih] - tj;
	     tbeta = eah[ki+ik-(*jmuprm)-1]/(td1+td2);
	     eah[ki+ik-(*jmuprm)-2] = td2*tbeta + tbeta1;
	     tbeta1 = td1*tbeta;
	  }
	kkh = kiu + ik - (*jmuprm) - 1;
	eah[kkh] = tbeta1;
	if (kiu < (*jmuprm)) 
	   eah[kkh] = tbeta1 + 
	   (etau[in+ik-1]-tj)*eah[kkh+1]/(etau[in+ik-1]-etau[kiu+1]);
     }
   
   /* Discrete B-splines computed.  */

   *jstat = 0;
   goto out;
   
   /* Error in space allocation.  */
   
   err101: *jstat = -101;
   goto out;
   
   out:
      /* Free scratch used for local array.  */
      
      if (ssi != SISL_NULL) freearray(ssi);
	  
      return;
}
   
   
