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
 * $Id: sh1839.c,v 1.2 2001-03-19 15:59:06 afr Exp $
 *
 */


#define SH1839

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void sh1839(SISLObject *po1,SISLObject *po2,double aepsge,int *jstat)
#else
void sh1839(po1,po2,aepsge,jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Perform improved box-test in intersection involving
*              a surface.
*
*
*
* INPUT      : po1     - Pointer to the first object in the intersection.
*                        This object is expected to be a surface.
*              po2     - Pointer to the second object in the intersection.
*              aepsge  - Geometry resolution.
*
*
*
* OUTPUT     : *jstat  - status messages  
*                            = 2      : Only edge-intersections possible.
*                            = 1      : Internal intersections possible.
*                            = 0      : No intersections possible.
*                            < 0      : error
*
*
* METHOD     : Pick the diagonal of the surface and the tangents in
*              the corners of the surface if these might be different 
*              from the diagonals. Rotate coordinate system according
*              to the vectors picked and perform a box-test.
*
*
* REFERENCES :
*
*-
* CALLS      : sh1834 - Perform rotated box-test.
*
* WRITTEN BY : Tor Dokken, SI, OSlo, Norway.
* REWISED BY : Vibeke Skytt, SI, 91-01.
*
*********************************************************************
*/
{
  int kstat = 0;           /* Local status variable.                        */
  int kpos = 0;            /* Position of error.                            */
  int ki,kj;               /* Counter.                                      */
  int kdim;                /* Dimension of geometry space.                  */
  int kn1,kn2;             /* Number of vertices in each parameter
			      direction of surface.                         */
  int kk1,kk2;             /* Order in each parameter direction of surface. */
  int kvec;                /* Number of direction vectors to calculate.     */
  int klap;                /* Indcates whether SISLbox of surface overlap point.*/
  double tang1,tang2;      /* Angles between direction vectors.             */
  double *scoef;           /* Vertices of surface.                          */
  register double *s1,*s2,
  *s3,*s4;                 /* Pointers used to traverse arrays.             */
  double *sdir = SISL_NULL;     /* Array containing direction vectors.           */
  
  /* Test kind of first object.  */
  
  if (po1->iobj != SISLSURFACE) goto err122;
    
  /* Copy surface to local parameters.  */
  
  kdim = po1 -> s1 -> idim;
  kn1 = po1 -> s1 -> in1;
  kn2 = po1 -> s1 -> in2;
  kk1 = po1 -> s1 -> ik1;
  kk2 = po1 -> s1 -> ik2;
  scoef = po1 -> s1 -> ecoef;
  
  /* Find number of rotations to make.  */
  
  if (kk1 > 2 || kk2 > 2) kvec = 10; else kvec = 2;
  
  /* Allocate space for vectors with which the x-axis is to be parallell. */
  
  sdir = newarray(kvec*kdim,double);
  if (sdir == SISL_NULL) goto err101;
  
   if (kvec == 2)
   {
  /* Make diagonal from lower left to upper right corner of patch.  
     s1 points to the array which contains the results, s3 points to the
     lower left corner and s4 to the upper right corner.                 */
  
  for (s1=sdir,s2=s1+kdim,s3=scoef,s4=scoef+kdim*(kn1*kn2-1); s1<s2;
       s1++,s3++,s4++)
    *s1 = *s4 - *s3;
  
  /* Make diagonal from upper left to lower right corner of patch. s1
     points to the array which contains the results, s3 points to the
     upper left and s4 to the lower right corner.                       */
  
  for (s1=sdir+kdim,s2=s1+kdim,s3=scoef+kdim*kn1*(kn2-1),
       s4=scoef+kdim*(kn1-1); s1<s2; s1++,s3++,s4++)
    *s1 = *s4 - *s3;
   }
   
  if (kvec > 2)
    {
      
      /* The surface is not linear in both parameter directions. Make
	 horizontal and vertical tangent in lower left corner. s1 points
	 to the array which contain the results and s3 to the corner.    */
      
      for (s1=sdir+2*kdim,s2=s1+kdim,s3=scoef; s1<s2; s1++,s3++)
	{       
	  *s1 = *(s3+kdim) - *s3;
	  *(s1+kdim) = *(s3+kdim*kn1) - *s3;
	}
      
      /* Make the horizontal and vertical tangent in lower right corner. s1
	 points to the array which contain the results and s3 to the corner.*/
      
      for (s1=sdir+4*kdim,s2=s1+kdim,s3=scoef+kdim*(kn1-1); s1<s2; s1++,s3++)
	{
	  *s1 = *(s3-kdim) - *s3;
	  *(s1+kdim) = *(s3+kdim*kn1) - *s3;
	}
      
      /* Make the horizontal and vertical tangent in upper left corner. s1
	 points to the result array and s3 to the corner.                  */
      
      for (s1=sdir+6*kdim,s2=s1+kdim,s3=scoef+kdim*kn1*(kn2-1);s1<s2;s1++,s3++)
	{
	  *s1 = *(s3+kdim) - *s3;
	  *(s1+kdim) = *(s3-kdim*kn1) - *s3;
	}
      
      /* Make the horizontal and vertical tangent in upper right corner. 
	 s1 points to the result array and s3 to the corner.             */
      
      for (s1=sdir+8*kdim,s2=s1+kdim,s3=scoef+kdim*(kn1*kn2-1);s1<s2;s1++,s3++)
	{
	  *s1 = *(s3-kdim) - *s3;
	  *(s1+kdim) = *(s3-kn1*kdim) - *s3;
	}
    }
  
  /* Rotate coordinate system according to the vectors found and perform
     box-test. First use the diagonal vectors.                           */
  
  klap = 1;
  if (kvec == 2)
  {
  sh1834(po1,po2,aepsge,kdim,sdir,sdir+kdim,&kstat);
  if (kstat < 0) goto error;
  klap = kstat;
  
  if (klap == 1)
    {
       sh1834(po1,po2,aepsge,kdim,sdir+kdim,sdir,&kstat);
      if (kstat < 0) goto error;
      klap = kstat;
    }
  }
  
  /* If the box-tests performed till now show overlap and the surface
     is non-linear in at least one direction rotate the geometry according
     to the tangent information gathered.                                 */
  /* First remove superfluous rotation directions.                        */
  
  for (ki=4; ki<kvec; )
  {
     for (kj=2; kj<4; kj+=2)
     {
	/* Test if the found vectors are aproximately equal.  */
	
	tang1 = s6ang(sdir+ki*kdim,sdir+kj*kdim,kdim);
	tang2 = s6ang(sdir+(ki+1)*kdim,sdir+(kj+1)*kdim,kdim);
	
	if (tang1 < ANGULAR_TOLERANCE && tang2 < ANGULAR_TOLERANCE) break;
     }
     
     if (kj < 4)
     {
	/* Remove set of rotation vectors.  */
	
	if (ki+2 < kvec)
	  memmove(sdir+ki*kdim, sdir+(ki+2)*kdim, (kvec-ki-2)*kdim*sizeof(double));
	//memcopy(sdir+ki*kdim, sdir+(ki+2)*kdim, (kvec-ki-2)*kdim, DOUBLE);
	kvec -= 2;
     }
     else ki+=2;
  }
  
  ki = 2;
  while (ki<kvec && klap == 1)
    {
       sh1834(po1,po2,aepsge,kdim,sdir+ki*kdim,sdir+(ki+1)*kdim,&kstat);
      if (kstat < 0) goto error;
      klap = kstat;
      
      if (klap && 
	  fabs(s6ang(sdir+ki*kdim,sdir+(ki+1)*kdim,kdim)-PIHALF) 
	  > ANGULAR_TOLERANCE)
      {
	 /* VSK, 01/93. Use the other partial derivative as x-axis in 
	    the rotation. */
	 
	 sh1834(po1,po2,aepsge,kdim,sdir+(ki+1)*kdim,sdir+ki*kdim,&kstat);
	 if (kstat < 0) goto error;
	 klap = kstat;
      }
      ki += 2;
    }
  
  /* Improved boxtest performed.  */
  
  *jstat = klap;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("sh1839",*jstat,kpos);
  goto out;
  
  /* Error in input. Unexpected object found.  */
  
 err122: *jstat = -122;
  s6err("sh1839",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("sh1839",*jstat,kpos);
  goto out;
  
 out: 
  
  /* Free allocated space.  */
  
  if (sdir != SISL_NULL) freearray(sdir);
  
  return;
}
