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
 * $Id: sh1790.c,v 1.2 2001-03-19 15:59:05 afr Exp $
 *
 */


#define SH1790

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void sh1790(SISLObject *po1,SISLObject *po2,int itype,
	    double aepsge,int *jstat)
#else
void sh1790(po1,po2,itype,aepsge,jstat)
     SISLObject *po1;
     SISLObject *po2;
     int    itype;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Perform a box-test on the two object given by
*              po1 and po2 and check if the boxes overlap.
*
*
*
* INPUT      : po1    - First object.
*              po2    - Second object.
*              itype  - Kind of box to make.
*                       = 0 : Do not expand box.
*                       = 1 : Make a totally expanded box.
*                       = 2 : Make a box expanded in the inner of the
*                             object, and reduced along the edges/endpoints.
*                       If itype>=10, it is interpreted as itype-10. This
*                       code is present to reduce the number of sides in
*                       the box to be made.
*	       aepsge - Geometry resolution.
*                                                                     
*
* OUTPUT     : jstat  - status messages  
*                            = 5      : Danger of shadow area in point
*                                       intersection.
*			     = 4      : One box is collapsed into a point.
*			     = 3      : Both boxes inside geometry resolution.
*                            = 2      : Overlap as closed sets only.
*                            = 1      : Overlap as open sets.
*                            = 0      : No overlap.
*                            < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-03.
*              Arne Laksaa, SI, 89-07.
* REWISED BY : Vibeke Skytt, SI, 91-01. Tolerance dependant boxes.
* CORRECTED BY: UJK and ALA, SI, 91-08. Removed edge test
* REWISED BY : VSK, SI, 92-10. Transform tolerance to other object when a point
*                              is included. Special check on "shadow area",
*                              which may arise in point-object intersection
*                              where the point lies outside the reduced box
*                              of the other object, but still within a 
*                              distance less than the tolerance from the
*                              object itself (in particular curve). This
*                              is possible near the endpoints/edges of the
*                              other object.
*********************************************************************
*/                                     
{
  int kstat = 0;        /* Local status error.                        */
  int kpos = 0;         /* Position of error.                         */
  int ktype = itype%10; /* Type of box to make.                       */
  int kant;             /* Number of sides in all boxes.              */
  int kdim1,kdim2;      /* Dimension of space.	    	   	      */
  int ki,kj=0;          /* Counters.                                  */
  int kpttest=0;        /* Indicates wether one object is a point and
			   the geometry has dimension greater than 1.
			   In that case the boundaries of the box is
			   placed further out in the min- and max-arrays. */
  int kshadow = 0;      /* Indicates if there is a risk of a shadow areas */
  int kbez = 0;         /* Indicates if any object is a bezier object.    */
  double tdist;         /* Distance between boxes.                    */
  double teps1;         /* To be used in 1D.                          */
  double teps2;         /* Two times aepsge if expanded box.          */
  double tepsge = aepsge; /* Tolerance to be used locally. In 1D 
			     tepsge=2*aepsge since the point is not
			     expended.                                */
  double t1,t2,t3,t4;   /* Help variables.                            */
  double t01,t02,t03,t04;   /* Help variables.                        */
  double *tmin1,*tmax1; /* Smallest and larges value of the vertices of
			   first object in each box in all dimension. */
  double *tmin2,*tmax2; /* Smallest and larges value of the vertices of
                           second object in each box in all dimension.*/
  
  /* Set local tolerance.  */
  
  teps2 = (ktype == 0) ? aepsge : (double)2.0*aepsge;
  
  /* Check if one object is a point. */
  
  if (po1->iobj == SISLPOINT) kdim1 = po1->p1->idim;
  else if (po1->iobj == SISLCURVE) kdim1 = po1->c1->idim;
  else if (po1->iobj == SISLSURFACE) kdim1 = po1->s1->idim;
  
  if (kdim1 == 1)
  {
     if (ktype != 0) teps2 = (double)3.0*aepsge;
     teps1 = (ktype == 0) ? aepsge : aepsge+aepsge;
     
     if (po1->iobj == SISLPOINT || po2->iobj == SISLPOINT)
	tepsge += aepsge;
  }
  else
  {
     teps1 = DZERO;
     if (po1->iobj == SISLPOINT || po2->iobj == SISLPOINT)
	kpttest = 1;
  }

	       
  /* Compute the box of the first object.  */ 
  
  sh1992(po1,itype,tepsge,&kstat);
  if (kstat < 0) goto error;
  
  /* VSK. 06.93. Adjust microbox tolerance in bezier case. */
 
  if (ktype != 0 && kstat == 1) kbez = 1;

  /* Set pointers to box boundaries.  */
  
  if (po1->iobj == SISLPOINT)
  {
     tmin1 = po1->p1->pbox->e2min[ktype];
     tmax1 = po1->p1->pbox->e2max[ktype];
  }
  else if (po1->iobj == SISLCURVE)
  {
     tmin1 = po1->c1->pbox->e2min[ktype];
     tmax1 = po1->c1->pbox->e2max[ktype];
  }       
  else if (po1->iobj == SISLSURFACE)
  {
     tmin1 = po1->s1->pbox->e2min[ktype];
     tmax1 = po1->s1->pbox->e2max[ktype];
  }         
  else goto err121;
  
  /* Make requested box.  */
		 
  sh1992(po2,itype,tepsge,&kstat);
  if (kstat < 0) goto error;

  /* VSK. 06.93. Adjust microbox tolerance in bezier case. */

  if (ktype != 0 && kstat == 1) kbez += 1;
  if (kdim1 == 1 && kbez > 0) teps2 -= (double)2.0*aepsge;
  else if (kbez > 0) teps2 -= aepsge;
  if (kbez) teps1 = aepsge;

  /* Set pointers to box boundaries.  */
  
  if (po2->iobj == SISLPOINT)
  {
     tmin2 = po2->p1->pbox->e2min[ktype];
     tmax2 = po2->p1->pbox->e2max[ktype];
     kdim2 = po2->p1->idim;
  }  
  else if (po2->iobj == SISLCURVE)
  {
     tmin2 = po2->c1->pbox->e2min[ktype];
     tmax2 = po2->c1->pbox->e2max[ktype];
     kdim2 = po2->c1->idim;
  }       
  else if (po2->iobj == SISLSURFACE)
  {
     tmin2 = po2->s1->pbox->e2min[ktype];
     tmax2 = po2->s1->pbox->e2max[ktype];
     kdim2 = po2->s1->idim;
  }       
  else goto err121;
  
  /* Check dimension. */
  
  if (kdim1 != kdim2 ) goto err106;
  else
    if (kdim1 < 1 )   goto err105;
  
  
  /* Compute total number of SISLbox edges. */
  
  if (itype < 10 && kdim1 == 3) kant = 9;
  else
    if (itype < 10 && kdim1 == 2) kant = 4;
    else           kant = kdim1;
  
  
  /* For each dimension in all boxes perform box-test. First check
     that we point on the right place in the min- and max-arrays.  */	
  
  if (kpttest)
  {
     tmin1 += kant;
     tmin2 += kant;
     tmax1 += kant;
     tmax2 += kant;
  }
  
  for (ki=0; ki<kant; ki++,tmin1++,tmax1++,tmin2++,tmax2++)
    {
      /* Sorting: t1-t2 The SISLbox with largest max value.
	 t3-t4 The other box. */
      
      if (*tmax1 > *tmax2)
	{
	  t1 = *tmax1;   t01 = *(tmax1 - kant*kpttest);
	  t2 = *tmin1;   t02 = *(tmin1 - kant*kpttest);
	  t3 = *tmax2;   t03 = *(tmax2 - kant*kpttest);
	  t4 = *tmin2;   t04 = *(tmin2 - kant*kpttest);
	}
      else
	{
	  t1 = *tmax2;   t01 = *(tmax2 - kant*kpttest);
	  t2 = *tmin2;   t02 = *(tmin2 - kant*kpttest);
	  t3 = *tmax1;   t03 = *(tmax1 - kant*kpttest);
	  t4 = *tmin1;   t04 = *(tmin1 - kant*kpttest);
	}
      tdist = t2 - t3;

      /* Use non-expanded boxes when testing if minibox is possible. */
      
      if (t01 - min(t02,t04) <= teps2)
	kj++;                 /* Minibox is possible. */
      else if (kdim1 == 1 && t01-t04 <= teps1 && t03-t02 <= teps1)
	kj++;                 /* Minibox is possible. */
      else if (t3 < t2 && (tdist > teps2 || !kpttest))
	{
	  *jstat = 0;           /* No overlap. */
	  goto out;
	}
      else if (t3 < t2)
	 kshadow = 1;           /* Danger of shadow area.  */

      /* UJK and ALA, 91-08: The tolerance is now
	 incorporated in the box. 
	 else if (kdim1 != 1 && t3 - t2 < tepsge && t3 - t4 > teps2 &&
	 t1 - t2 > teps2)
	 {
	 *jstat = 2;           Only edge touching possible. 
	 goto degenerate;
	 }
	 else if ( kdim1 == 1 && (t1 - t4 < tepsge || t3 - t2 < tepsge))
	 {
	 *jstat = 2;            Only edge touching possible. 
	 goto out;
	 } */
	
      /* else possible overlap. */
    }
  
  if (kj == kant)
    *jstat = 3;                   /* Minibox found. */
  else if (kshadow)
     *jstat = 5;                  /* Danger of shadow area.  */
  else
    *jstat = 1;                   /* Overlap.  */
  
  
  /* Box-test performed. */
  
  /* degenerate: */
  /* Test if one of the objects has collapsed. Use non-expanded box.  */
  if (kdim1 != 1 && po1->iobj > SISLPOINT)
    {
      if (po1->iobj == SISLCURVE)
	{
	  tmin1 = po1->c1->pbox->e2min[ktype];
	  tmax1 = po1->c1->pbox->e2max[ktype];
	}
      
      else if (po1->iobj == SISLSURFACE)
	{
	  tmin1 = po1->s1->pbox->e2min[ktype];
	  tmax1 = po1->s1->pbox->e2max[ktype];
	}
      
      /* Make sure that we are correct place in the min- and max-arrays. */
      
      if (kpttest)
      {
	 tmin1 += kant;
	 tmax1 += kant;
      }
      
      for(ki=0;ki<kdim1;ki++,tmin1++,tmax1++)
	if (fabs(tmax1[0]-tmin1[0]) > tepsge) break;
      
      if (ki == kdim1+1)
	{
	  *jstat = 4;
	  goto out;
	}
    }
  
  if (kdim1 != 1 && po2->iobj > SISLPOINT)
    {
      if (po2->iobj == SISLCURVE)
	{
	  tmin1 = po2->c1->pbox->e2min[ktype];
	  tmax1 = po2->c1->pbox->e2max[ktype];
	}
      
      else if (po2->iobj == SISLSURFACE)
	{
	  tmin1 = po2->s1->pbox->e2min[ktype];
	  tmax1 = po2->s1->pbox->e2max[ktype];
	}
      
      /* Make sure that we are correct place in the min- and max-arrays. */
      
      if (kpttest)
      {
	 tmin1 += kant;
	 tmax1 += kant;
      }
      
      for(ki=0;ki<kdim1;ki++,tmin1++,tmax1++)
	if (fabs(tmax1[0]-tmin1[0]) > tepsge) break;
      
      if (ki == kdim1+1)
	{
	  *jstat = 4;
	  goto out;
	}
    }
  goto out;
  
  /* Dimensions conflicting. */
  
 err106: *jstat = -106;
  s6err("sh1790",*jstat,kpos);
  goto out;
  
  /* Dimensions less than one. */
  
 err105: *jstat = -105;
  s6err("sh1790",*jstat,kpos);
  goto out;
  
  /* Kind of object does not exist. */
  
 err121: *jstat = -121;
  s6err("sh1790",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine. */
  
 error:  *jstat = kstat;
  s6err("sh1790",*jstat,kpos);
  goto out;
  
 out:	return;
}





