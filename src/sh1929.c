/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of SISL.
 *
 * SISL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * SISL is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with SISL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using SISL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the SISL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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
   
   
