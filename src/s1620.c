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
 * $Id: s1620.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */

#define S1620

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1620(double epoint[],int inbpnt1, int inbpnt2, int ipar,
	   int iopen1, int iopen2, int ik1, int ik2, int idim,
	   SISLSurf **rs,int *jstat)
#else
void s1620(epoint,inbpnt1,inbpnt2,ipar,
	   iopen1,iopen2,ik1,ik2,idim,rs,jstat)
     double epoint[];
     int inbpnt1;
     int inbpnt2;
     int ipar;
     int iopen1;
     int iopen2;
     int ik1;
     int ik2;
     int idim;
     SISLSurf **rs;
     int *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate a B-spline surface using the input points as
*              control vertices. The parametrization is calculated
*              according to ipar.
*
* INPUT      : epoint   - The array containing the points to be used as
*                         controlling vertices of the B-spline curve.
*              inbpnt1  - The number of points in first par. direction
*              inbpnt2  - The number of points in second par. direction
*              ipar     - Flag showing the desired parametrization to be 
*                         used:
*                          = 1: Mean accumulated cord-length parameterization.
*                          = 2: Uniform parametrization.
*              iopen1   - Open/close condition (Open=1,Close=0,Periodic=-1)
*              iopen2   - Open/close condition (Open=1,Close=0,Periodic=-1)
*              ik1      - The order of the surface in first direction.
*              ik2      - The order of the surface in second direction.
*              idim     - The dimension of the space
*
* OUTPUT     : rs       - Pointer to the surface.
*              jstat    - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : First the parametrization of the surface is calculated.
*              If more than ik adjacent vertices are equal then the
*              superfluous vertices are removed. Then the knots are
*              calculated.
*
* REFERENCES : See 'equivalent routine' for curves: s1620
*
*                                                 
* CALLS      : s1902,s1750,newSurf,s6err
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF OSLO, June 1993.
* REWISED BY : Vibeke Skytt, SINTEF, 9403. Change in input parameters
*                                          iopen1 and iopen2.
*
*********************************************************************
*/
{
  int kstat;          /* Status variable                                 */
  int kn1, kn2;       /* The number of B-splines in each directions, 
			 i.e., the dimension of the spline space 
			 associated with the knot vector.                */
  int kk1, kk2;       /* The polynomial orders of the surafce.           */
  int kpos=0;         /* Position of error                               */
  int j;              /* Counter for loop control                        */
  int kopen1, kopen2; /* Local open/closed parameter. Closed,
			 non-periodic is treated as an open curve.       */
  double *par1=SISL_NULL;  /* Pointer to parameterization array in first
		       * direction. */
  double *par2=SISL_NULL;  /* Pointer to parameterization array in second
		       * direction. */
  double *scoef=SISL_NULL; /* Pointer to vertex array                         */
  double *knot1=SISL_NULL; /* Pointer to knot vector in first par. direction  */
  double *knot2=SISL_NULL; /* Pointer to knot vector in second par. direction */
  SISLSurf *qs=SISL_NULL;
  
  /* Set local open/closed parameter. */
  
  kopen1 = (iopen1 == SISL_CRV_PERIODIC) ? 0 : 1;
  kopen2 = (iopen2 == SISL_CRV_PERIODIC) ? 0 : 1;
  
  /* Control input     */
  
  kk1 = MIN(inbpnt1,ik1);
  kk2 = MIN(inbpnt2,ik2);
  if (ik1 < 2 || ik2 < 2) goto err109;
  if  (iopen1 != SISL_CRV_OPEN && iopen1 != SISL_CRV_CLOSED &&
       iopen1 != SISL_CRV_PERIODIC) goto err109;
  if  (iopen2 != SISL_CRV_OPEN && iopen2 != SISL_CRV_CLOSED &&
       iopen2 != SISL_CRV_PERIODIC) goto err109;
  
  /* Generate parametrizations */
  
  s1528(idim, inbpnt1, inbpnt2, epoint, ipar, iopen1, iopen2,
	&par1, &par2, &kstat);
  if(kstat < 0) goto error;

  /* Find knot vector in first parameter direction */

  s1902(par1,inbpnt1+(iopen1==SISL_CRV_CLOSED),kk1,kopen1,&knot1,&kstat);
  if (kstat < 0 || knot1 == SISL_NULL) goto error;

  /* Find knot vector in second parameter direction */

  s1902(par2,inbpnt2+(iopen2==SISL_CRV_CLOSED),kk2,kopen2,&knot2,&kstat);
  if (kstat < 0 || knot2 == SISL_NULL) goto error;
  
  /* Allocate space for vertice array   */
  
  scoef = newarray((inbpnt1+kk1-1)*(inbpnt2+kk2-1)*idim,DOUBLE);
  if (scoef == SISL_NULL) goto err101;
  
  /* Check if closed surface in first direction. If closed, add the 
     (ik-1) first points to the vertice array for each j=1..inbpnt2. */
  
  if (iopen1 == SISL_CRV_PERIODIC)
    {
      kn1 = inbpnt1 + kk1 - 1;
      for (j=0; j<inbpnt2; j++)
	{
	  memcopy (&scoef[j*kn1*idim],&epoint[j*inbpnt1*idim],
		   inbpnt1*idim,DOUBLE);  
	  memcopy(&scoef[(j*kn1+inbpnt1)*idim],&epoint[j*inbpnt1*idim],
		  (kk1-1)*idim,DOUBLE);
	}
    }
  else if (iopen1 == SISL_CRV_CLOSED)
  {
      kn1 = inbpnt1 + 1;
      for (j=0; j<inbpnt2; j++)
	{
	  memcopy (&scoef[j*kn1*idim],&epoint[j*inbpnt1*idim],
		   inbpnt1*idim,DOUBLE);  
	  memcopy(&scoef[(j*kn1+inbpnt1)*idim],&epoint[j*inbpnt1*idim],
		  idim,DOUBLE);
	}
  }
  else
    {
      /* Surface is open in first direction */

      kn1 = inbpnt1;
      memcopy (scoef,epoint,inbpnt1*inbpnt2*idim,DOUBLE);
    }

  /* Check if closed surface in second direction. If closed, add the 
     (kk2-1)*kn1 first points to the vertice array. */
  
  if (iopen2 == SISL_CRV_PERIODIC)
    {
      kn2 = inbpnt2 + kk2 - 1;
      memcopy(&scoef[inbpnt2*kn1*idim],scoef,
	      (kk2-1)*kn1*idim,DOUBLE);
    }
  else if (iopen2 == SISL_CRV_CLOSED)
  {
     kn2 = inbpnt2 + 1;
     memcopy(&scoef[inbpnt2*kn1*idim],scoef,
	     kn1*idim,DOUBLE);
  }
  else
     
    /* Surface is open in first direction */
    
    kn2 = inbpnt2;
  

  /* Make surface */
  
  qs = newSurf(kn1, kn2, kk1, kk2, knot1, knot2, scoef, 1, idim, 1);
  if (!qs) goto err101;                

  /* Set peridicity parameters. */
  
  qs->cuopen_1 = iopen1;
  qs->cuopen_2 = iopen2;
  
  /* Increase the orders if the orders were lowered when controlling input */
  
  if (kk1 < ik1 || kk2 < ik2)
    {
      s1387(qs, ik1, ik2, &qs, &kstat);
      if (kstat< 0) goto error;
    }
  
  if (qs != SISL_NULL) *rs = qs;
  *jstat = 0;
  goto out;
  
  /* Error in memory allocation */
  
  err101: 
    *jstat = -101;
    s6err("s1620",*jstat,kpos);
    goto out;
  
  /* Error in input, order less than 2 */
  
  err109: 
    *jstat = -109;
    s6err("s1620",*jstat,kpos);
    goto out;
  
  /* Error in lower level function */  
  
  error:  
    *jstat = kstat;
    s6err("s1620",*jstat,kpos); 
    goto out;
  
  out:
    if (knot1 != SISL_NULL) freearray(knot1);
    if (knot2 != SISL_NULL) freearray(knot2);
    if (par1  != SISL_NULL) freearray(par1);
    if (par2  != SISL_NULL) freearray(par2);
    if (scoef != SISL_NULL) freearray(scoef);
    return;
}    


