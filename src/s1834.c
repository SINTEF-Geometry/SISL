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
 * $Id: s1834.c,v 1.3 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1834

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void
s1834_s9mat2d(double [],double []);
static void
s1834_s9mat3d(double [],double [],double []);
#else
static void s1834_s9mat2d();
static void s1834_s9mat3d();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
s1834(double ecoef1[],int in1,double ecoef2[],int in2,int idim,
	   double edir1[],double edir2[],int *jstat)
#else
void s1834(ecoef1,in1,ecoef2,in2,idim,edir1,edir2,jstat)
     double ecoef1[];
     int    in1;
     double ecoef2[];
     int    in2;
     int    idim;
     double edir1[];
     double edir2[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Perform a box-test to check if two objects overlap in 
*              the rotated coordinate system where edir1 defines the 
*              x-axis and (edir1 x edir2) defines the z-axis if
*              idim = 3.
*
*
*
* INPUT      : ecoef1 - Coefficients of the first object.
*              in1    - Number of coefficients in the first object.
*              ecoef2 - Coefficients of the second object.
*              in2    - Number of coefficients in the second object.
*              idim   - The dimension of the space in which the objects
*                       lie. idim = 2 or idim = 3.
*              edir1  - First direction vector. Defines x-axis in rotated
*                       coordinate system.
*              edir2  - Second direction vector.
*
*
*
* OUTPUT     : jstat  - status messages  
*                                 = 2      : Boundaries just touch.
*                                 = 1      : Rotated SISLbox overlaps point.
*                                 = 0      : No overlap.
*                                 < 0      : error
*
*
* METHOD     : The coordinate system is rotated such that if idim = 2,
*              the x-axis of the new coordinate system is parallell to 
*              the vector edir1. If idim = 3, the cross-product of edir1
*              and edir2 is rotated to be parallell to the z-axis and
*              edir1 rotated to be parallell to the x-axis. The objects
*              are moved into this rotated coordinate system and a 
*              box-test is performed.
*
*
* REFERENCES :
*
*-
* CALLS      : s6norm - Normalize vector and compute length of original
*                       vector.
*              s6scpr - Compute scalar product of two vectors.
*              s6crss - Make cross product of two 3D vectors. 
*              s1834_s9mat2d - Set up rotation matrix for 2 dimensional space.
*              s1834_s9mat3d - Set up rotation matrix for 3 dimensional space.
*                              
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
*
*********************************************************************
*/
{                                   
  int kpos = 0;    /* Position of error.                         */
  int kopen;       /* Indicates if boxes overlap as open sets.   */
  int kclose;      /* Indicates if boxes oberlap as closed sets. */
  double tleng;    /* Scalar product between two vectors.        */
  double *smat;    /* Rotation matrix.                           */
  double *smin1,*smax1; /* Extremal values of first rotated box. */
  double *smin2,*smax2; /* Extremal values of second rotated box. */
  double *s1,*s2,*s3,*s4;  /* Pointers used to traverse arrays.  */
  
  smat = smin1 = smin2 = smax1 = smax2 = SISL_NULL;
  
  /* Test input.  */
  
  if (idim != 2 && idim != 3) goto err105;
  
  /* Allocate space for local parameters.  */
  
  smin1 = newarray(idim,double);
  smin2 = newarray(idim,double);
  smax1 = newarray(idim,double);
  smax2 = newarray(idim,double);
  smat = new0array(idim*idim,double);
  if (smin1 == SISL_NULL || smin2 == SISL_NULL || smax1 == SISL_NULL ||
      smax2 == SISL_NULL || smat == SISL_NULL) goto err101;
  
  /* Initialize min and max vectors.  */
  
  for (s2=smin1+idim; smin1<s2; smin1++,smin2++,smax1++,smax2++)
    {
      *smin1 = *smin2 = (double)9999999999.0;
      *smax1 = *smax2 = (double)-9999999999.0;
    }
  smin1 -= idim;
  smin2 -= idim;
  smax1 -= idim;
  smax2 -= idim;
  
  /* Find the rotation matrix.  */
  
  if (idim == 2)
    
    /* After normalization edir1[0] will contain the cosine of the 
       rotation angle and edir1[1] will contain the sine.           */
    
    s1834_s9mat2d(smat,edir1);
  else
    
    /* Set up the rotation matrix when idim = 3. (edir1 x edir2) is
       rotated to be parallell to the z-axis and edir1 to be parallell
       to the x-axis.                                                   */
    
    s1834_s9mat3d(smat,edir1,edir2);
  
  /* The objects is moved into the new coordinate system by rotating
     them using the rotation matrix.                                 */
  
  /* Rotate first object and make a SISLbox around the object. */
  
  for (s1=smat,s3=smat+idim*idim; s1<s3; s1+=idim,smin1++,smax1++)
    for (s2=ecoef1,s4=s2+idim*in1; s2<s4; s2+=idim)
      {
	tleng = s6scpr(s1,s2,idim);
	*smin1 = MIN(*smin1,tleng);
	*smax1 = MAX(*smax1,tleng);
      }
  smin1 -= idim;
  smax1 -= idim;
  
  /* Rotate second object and make a SISLbox around the object. */
  
  for (s1=smat,s3=smat+idim*idim; s1<s3; s1+=idim,smin2++,smax2++)
    for (s2=ecoef2,s4=s2+idim*in2; s2<s4; s2+=idim)
      {
	tleng = s6scpr(s1,s2,idim);
	*smin2 = MIN(*smin2,tleng);
	*smax2 = MAX(*smax2,tleng);
      }
  smin2 -= idim;
  smax2 -= idim;
  
  /* Check if the boxes overlap.  */
  
  kopen = 0;
  kclose = 1;
  for (s2=smin1+idim; smin1<s2; smin1++,smax1++,smin2++,smax2++)
    {
      if (DEQUAL(MIN(*smin1,*smin2),MAX(*smax1,*smax2)))
	continue;
      
      if (*smin1 > *smax2 || *smin2 > *smax1)
	
	/* The objects do not intersect as closed sets.  */
	
	kclose = 0;
      else if (DEQUAL(*smin1,*smax2) || DEQUAL(*smin2,*smax1))
	
	/* The objects do not intersect as open sets.  */
	
	kopen = 1;
    }
  smin1 -= idim;
  smax1 -= idim;
  smin2 -= idim;
  smax2 -= idim;
  
  if (kopen == 0 && kclose == 1)
    
    /* The objects intersect as open sets.  */
    
    *jstat = 1;
  else if (kopen == 1)
    
    /* The objects intersect as closed sets and not as open.  */
    
    *jstat = 2;                                            
  else
    
    /* The objects do not intersect.  */
    
    *jstat = 0;
  
  /* Box-test permformed.  */
  
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1834",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 2 or 3.  */
  
 err105: *jstat = -105;
  s6err("s1834",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free space occupied by local arrays.  */
  
  if (smin1 != SISL_NULL) freearray(smin1);
  if (smin2 != SISL_NULL) freearray(smin2);
  if (smax1 != SISL_NULL) freearray(smax1);
  if (smax2 != SISL_NULL) freearray(smax2);
  if (smat != SISL_NULL) free0array(smat);
  
  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1834_s9mat2d(double emat[],double edir[])
#else
static void s1834_s9mat2d(emat,edir)
     double emat[];
     double edir[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Set up rotation matrix in two dimensions when the x-axis
*              is supposed to be rotated to be parallell with edir.
*
* INPUT      : edir  - Direction of rotated x-axis.
*
* OUTPUT     : emat  - Rotation matrix. emat is supposed to be
*                      initialized to zero before this routine is entered.
*
*********************************************************************
*/
{          
  int kstat = 0;   /* Local status variable.              */
  double tlength;  /* Length of vector edir.              */
  double sdir[2];  /* Normalized vertion of vector edir.  */
  
  tlength = s6norm(edir,2,sdir,&kstat);
  if (kstat == 0)
    
    /* Length of edir equal to zero. Let the rotation matrix be
       the identity matrix.                                      */
    
    emat[0] = emat[3] = (double)1.0;
  else
    {          
      
      /* Make rotation matrix.  */
      
      emat[0] = sdir[0];
      emat[1] = -sdir[1];
      emat[2] = sdir[1];
      emat[3] = sdir[0];
    }
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1834_s9mat3d(double emat[],double edir1[],double edir2[])
#else
static void s1834_s9mat3d(emat,edir1,edir2)
     double emat[];
     double edir1[];
     double edir2[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Set up rotation matrix in three dimensions when edir1
*              is supposed to be rotated to be parallell with the x-axis
*              and (edir1 x edir2) to be parallell with the z-axis.
*
* INPUT      : edir1 - First direction vector.
*              edir2 - Second direction vector.
*
* OUTPUT     : emat  - Rotation matrix. The matrix is supposed to be
*                      initialized to zero before this routine is enterd.
*
*********************************************************************
*/
{          
  int kstat = 0;    /* Local status variable.                         */
  double snorm[3];  /* Cross-product of edir1 and edir2.              */
  double sdir[3];   /* Normalized vertion of edir1.                   */
  double *s1;       /* Pointer into emat array.                       */
  double tleng1,tleng2; /* Length of snorm and edir1 respectively.    */
  double ta1,ta2,ta3,tb1,tb2,tb3,td1,td2,tl1,tl2,tl3; /* Help variables. */
  
  /* Calculate cross-product of edir1 and edir2.  */
  
  s6crss(edir1,edir2,snorm);
  
  /* Normalize snorm.  */
  
  tleng1 = s6norm(snorm,3,snorm,&kstat);
  
  /* Normalize edir1.  */
  
  tleng2 = s6norm(edir1,3,sdir,&kstat);
  
  /* Initialize help variables.  */
  
  ta1 = snorm[0];
  ta2 = snorm[1];
  ta3 = snorm[2];
  tl1 = sqrt(ta2*ta2+ta3*ta3);
  
  /* Set up rotation matrix.  */
  
  if ((DEQUAL(tleng1,DZERO) || DEQUAL(tl1,DZERO)) && DEQUAL(tleng2,DZERO))
    
    /* The rotation matrix is the identity matrix.  */
    
    emat[0] = emat[4] = emat[8] = (double)1.0;
  else if (DEQUAL(tleng1,DZERO) || DEQUAL(tl1,DZERO))
    {
      
      /* The rotation matrix is supposed to rotate edir1 to be parallell
	 to the x-axis.                                                   */
      
      tb1 = sdir[0];
      tb2 = sdir[1];
      tb3 = sdir[2];
      tl3 = sqrt(tb1*tb1+tb2*tb2);
      
      if (DEQUAL(tl3,DZERO)) emat[0] = emat[4] = emat[8] = (double)1.0;
      else
	{
	  s1      = emat;    
	  *(s1++) = tb1;
	  *(s1++) = tb2;
	  *(s1++) = tb3;
	  *(s1++) = -tb2/tl3;
	  *(s1++) = tb1/tl3;  
	  *(s1++) = DZERO;
	  *(s1++) = -tb1*tb3/tl3;
	  *(s1++) = -tb2*tb3/tl3;  
	  *(s1++) = tl3;
	}
    }
  else
    {
      td1 = edir1[0]/tl1;
      td2 = (ta3*edir1[1] - ta2*edir1[2])/tl1;
      tl2 = sqrt(td1*td1+td2*td2);
      
      if (DEQUAL(tl2,DZERO))
	{
	  
	  /* The normal snorm is rotated to be parallell to the z-axis. */
	  
	  s1      = emat;
	  *(s1++) = tl1;
	  *(s1++) = -ta1*ta2/tl1;
	  *(s1++) = -ta1*ta3/tl1;
	  *(s1++) = DZERO;
	  *(s1++) = ta3/tl1;
	  *(s1++) = -ta2/tl1;
	  *(s1++) = ta1;
	  *(s1++) = ta2;
	  *(s1++) = ta3;
	}
      else
	{
	  
	  /* The normal is rotated to be parallell to the z-axis and edir1
	     to be parallell to the x-axis.                                 */
	  
	  s1      = emat;
	  *(s1++) = td1*tl1/tl2;
	  *(s1++) = (-ta1*ta2*td1 + ta3*td2)/(tl1*tl2);
	  *(s1++) = (-ta1*ta3*td1 - ta2*td2)/(tl1*tl2);
	  *(s1++) = -td2*tl1/tl2;
	  *(s1++) = (ta1*ta2*td2 + ta3*td1)/(tl1*tl2);
	  *(s1++) = (ta1*ta3*td2 - ta2*td1)/(tl1*tl2);
	  *(s1++) = ta1;
	  *(s1++) = ta2;
	  *(s1++) = ta3;
	}
    }
}
