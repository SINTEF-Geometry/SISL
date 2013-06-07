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

#define S1528

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1528(int idim, int m1, int m2, double points[], int ipar, 
	   int iopen1, int iopen2, double **par1, double **par2, int *jstat)
#else
void s1528(idim,m1,m2,points,ipar,iopen1,iopen2,par1,par2,jstat)
     int idim;
     int m1;
     int m2;
     double points[];
     int ipar;
     int iopen1;
     int iopen2;
     double **par1;
     double **par2;
     int *jstat;
#endif
/*
************************************************************************
*
* PURPOSE: To compute a suitable parametrization for a given set of points
*          according to descriptor ipar.
*
* INPUT:
*          idim   - The space dimension.
*          m1     - The number of points in first par. direction
*          m2     - The number of points in second par. direction
*          points - Array of dimension idim*m1*m2 containing the
*                   points.
*          ipar   - Flag showing the desired parametrization to be used:
*                   = 1: Mean accumulated cord-length parameterization.
*                   = 2: Uniform parametrization.
*          iopen1 - Flag telling if the surface should be open in 1. par. dir.:
*                     1 : The curve should be open.
*                     0 : The curve should be closed.
*                    -1 : The curve should be closed and periodic.
*          iopen2 - Flag telling if the curve should be open in 2. par. dir.
*
* OUTPUT:
*          par1   - Array containing parametrization for first direction.
*                   The dimension is m1+(iopen1!=1).
*          par2   - Array containing parametrization for second direction.
*                   The dimension is m2+(iopen2!=1).
*          jstat  - Status variable
*                    < 0 - Memory allocation error.
*
* METHOD:
*
* REFERENCES :
*
* CALLS:
*
* WRITTEN BY: Christophe Rene Birkeland, SINTEF, May 1993.
* REVISED BY: Vibeke Skytt, SINTEF, 0394. Introduced open/closed parameter.
*
*********************************************************************
*/
{
  int i,j;             /* Loop variables                              */
  int kpek1, kpek2, kpek3; /* Loop variables                          */
  int kpr1, kpr2;      /* Type of parametrization flag                */
  int kpos=0;          /* Position of error                           */
  int idimm1;          /*  =  idim * m1                               */
  double tdist;        /* Help parameter.                             */
  double *local_par1=SISL_NULL; /* Parametrization arrays used in routine  */
  double *local_par2=SISL_NULL;


  /* Allocate arrays for parametrizations */

  local_par1 = newarray(m1+(iopen1 != SISL_CRV_OPEN), DOUBLE); 
  local_par2 = newarray(m2+(iopen2 != SISL_CRV_OPEN), DOUBLE); 
  if(local_par1 == SISL_NULL || local_par2 == SISL_NULL) goto err101; 
  local_par1[0] = (double) 0.0;
  local_par2[0] = (double) 0.0;

  /* Compute parametrization */

  kpr1 = ipar; kpr2 = ipar;

  if (ipar == 1)    
    {
      /* Chord length parametrization in
       * first direction */

      kpek1 = 0;    
      idimm1 = idim*m1; 
      for (i=1; i<m1; i++)
	{    
	  local_par1[i] = local_par1[i-1];      
	  kpek2 = kpek1 + idim;
	  tdist = (double)0;
	  for (j=0, kpek3=0 ; j<m2; j++, kpek3+=idimm1) 
	    {
	      tdist += s6dist(&points[kpek2 + kpek3],
				      &points[kpek1 + kpek3], idim);          
	    }
	  local_par1[i] += tdist/(double)m2;
	  kpek1 = kpek2;
	} 
      
      if (iopen1 != SISL_CRV_OPEN)
      {
	 local_par1[m1] = local_par1[m1-1];      
	 kpek2 = 0;
	 tdist = (double)0;
	 for (j=0, kpek3=0 ; j<m2; j++, kpek3+=idimm1) 
	 {
	    tdist += s6dist(&points[kpek2 + kpek3],
				    &points[kpek1 + kpek3], idim);          
	 }
	 local_par1[m1] += tdist/(double)m2;
      } 
      
      if (local_par1[m1-1] == (double)0.) kpr1 = 2;
      
      /* Chord length parametrization in second direction */
      
      kpek1 = 0; 
      for (j=1; j<m2; j++)
	{    
	  local_par2[j] = local_par2[j-1];
	  kpek2 = kpek1 + idimm1; 
	  tdist = (double)0;
	  for (i=0, kpek3=0 ;i<m1 ; i++, kpek3+=idim)      
	    {                
	      tdist += s6dist(&points[kpek2 + kpek3],
				      &points[kpek1 + kpek3], idim);
	    }
	  local_par2[j] += tdist/(double)m1;
	  kpek1 = kpek2;     
	}
      
      if (iopen2 != SISL_CRV_OPEN)
      {
	 local_par2[m2] = local_par2[m2-1];      
	 kpek2 = 0;
	 tdist = (double)0;
	 for (i=0, kpek3=0 ; i<m1; i++, kpek3+=idim) 
	 {
	    tdist += s6dist(&points[kpek2 + kpek3],
				    &points[kpek1 + kpek3], idim);          
	 }
	 local_par2[m2] += tdist/(double)m1;
      } 
      
      if(local_par2[m2-1] == (double)0.) kpr2 = 2; 
    }
  if (kpr1 == 2)    
    {
      /*Uniform parametrization in first direction */
      
      for (i=1; i<m1+(iopen1!=SISL_CRV_OPEN); i++)
	local_par1[i] = i; 
    }
  if (kpr2 == 2)
    {
      /* Uniform parametrization in second direction */
      
      for (i=1; i<m2+(iopen2!=SISL_CRV_OPEN); i++) 
	local_par2[i] = i;
    }


  /* Success */

  *par1 = local_par1;
  *par2 = local_par2;
  *jstat = 0;
  goto out;


  /* Error in space allocation */

  err101: *jstat = -101;
    s6err("s1531",*jstat,kpos);
    goto out;

  out:
    return;
}
