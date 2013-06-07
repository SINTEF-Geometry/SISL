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
 * $Id: s1790.c,v 1.2 2001-03-19 15:58:54 afr Exp $
 *
 */


#define S1790

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1790(SISLObject *po1,SISLObject *po2,double aepsge,int *jstat)
#else
void s1790(po1,po2,aepsge,jstat)
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
* PURPOSE    : Perform a box-test on the two object given by
*              po1 and po2 and check if the boxes overlap.
*
*
*
* INPUT      : po1    - First object.
*              po2    - Second object.
*	       aepsge - Geometry resolution.
*                                                                     
*
* OUTPUT     : jstat  - status messages  
*			     = 4      : One SISLbox is collapsed into a point.
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
*
*********************************************************************
*/                                     
{
  int kstat = 0;        /* Local status error.                        */
  int kpos = 0;         /* Position of error.                         */
  int kant;             /* Number of vertices in all boxes.           */
  int kbez = 0;         /* Flag to mark bezier curve or patch.         */
  int kdim;	      /* Dimension of space.			    */
  int ki,kj=0;          /* Counters.                                  */
  double t1,t2,t3,t4;   /* Help variables.                            */
  double *tmin1,*tmax1; /* Smallest and larges value of the vertices of
			   first object in each SISLbox in all dimension. */
  double *tmin2,*tmax2; /* Smallest and larges value of the vertices of
			   second object in each SISLbox in all dimension.*/
  
  /* Check kind of first object. */
  
  if (po1->iobj == SISLPOINT)
    {
      /* Fetch dimention of the object. */
      
      kdim = po1->p1->idim;
      
      
      /* Check if the SISLbox have been computed. */
      
      if (po1->p1->pbox == SISL_NULL)
	{
	  /* If not compute a box. */
	  
	  s1992(po1,&kstat);
	  if (kstat<0) goto error;
	}
      
      /* Fetch the SISLbox boarder. */
      
      tmax1 = po1->p1->pbox->emax;
      tmin1 = po1->p1->pbox->emin;
    }
  else
    if (po1->iobj == SISLCURVE)
      {
	/* Fetch dimention of the object. */
	
	kdim = po1->c1->idim;
	
	/* Check if we have a bezier curve. */
	
	if (po1->c1->in == po1->c1->ik) kbez = 1;
	
	/* Check if the SISLbox have been computed. */
	
	if (po1->c1->pbox == SISL_NULL)
	  {
	    /* If not compute a box. */
	    
	    s1992(po1,&kstat);
	    if (kstat<0) goto error;
	  }
	
	/* Fetch the SISLbox boarder. */
	
	tmax1 = po1->c1->pbox->emax;
	tmin1 = po1->c1->pbox->emin;
      }
    else
      if (po1->iobj == SISLSURFACE)
	{
	  /* Fetch dimention of the object. */
	  
	  kdim = po1->s1->idim;
	  
	  /* Check if we have a bezier patch. */
	  
	  if (po1->s1->in1 == po1->s1->ik1 &&
	      po1->s1->in2 == po1->s1->ik2)    kbez = 1;
	  
	  /* Check if the SISLbox have been computed. */
	  
	  if (po1->s1->pbox == SISL_NULL)
	    {
	      /* If not compute a box. */
	      
	      s1992(po1,&kstat);
	      if (kstat<0) goto error;
	    }
	  
	  /* Fetch the SISLbox boarder. */
	  
	  tmax1 = po1->s1->pbox->emax;
	  tmin1 = po1->s1->pbox->emin;
	}
      else  goto err121;
  
  
  /* Check kind of second object. */
  
  if (po2->iobj == SISLPOINT)
    {
      /* Fetch dimention of the object. */
      
      ki = po2->p1->idim;
      
      
      /* Check if the SISLbox have been computed. */
      
      if (po2->p1->pbox == SISL_NULL)
	{
	  /* If not compute a box. */
	  
	  s1992(po2,&kstat);
	  if (kstat<0) goto error;
	}
      
      /* Fetch the SISLbox boarder. */
      
      tmax2 = po2->p1->pbox->emax;
      tmin2 = po2->p1->pbox->emin;
    }
  else
    if (po2->iobj == SISLCURVE)
      {
	/* Fetch dimention of the object. */
	
	ki = po2->c1->idim;
	
	/* Check if we have a bezier curve. */
	
	if (po2->c1->in == po2->c1->ik) kbez = 1;
	
	/* Check if the SISLbox have been computed. */
	
	if (po2->c1->pbox == SISL_NULL)
	  {
	    /* If not compute a box. */
	    
	    s1992(po2,&kstat);
	    if (kstat<0) goto error;
	  }
	
	/* Fetch the SISLbox boarder. */
	
	tmax2 = po2->c1->pbox->emax;
	tmin2 = po2->c1->pbox->emin;
      }
    else
      if (po2->iobj == SISLSURFACE)
	{
	  /* Fetch dimention of the object. */
	  
	  ki = po2->s1->idim;
	  
	  /* Check if we have a bezier patch. */
	  
	  if (po2->s1->in1 == po2->s1->ik1 &&
	      po2->s1->in2 == po2->s1->ik2)    kbez = 1;
	  
	  /* Check if the SISLbox have been computed. */
	  
	  if (po2->s1->pbox == SISL_NULL)
	    {
	      /* If not compute a box. */
	      
	      s1992(po2,&kstat);
	      if (kstat<0) goto error;
	    }
	  
	  /* Fetch the SISLbox boarder. */
	  
	  tmax2 = po2->s1->pbox->emax;
	  tmin2 = po2->s1->pbox->emin;
	}
      else goto err121;
  
  /* Check dimension. */
  
  if (ki != kdim ) goto err106;
  else
    if (kdim < 1 )   goto err105;
  
  
  /* Compute total number of SISLbox edges. */
  
  if (kdim == 3) kant = 12;
  else
    if (kdim == 2) kant = 4;
    else           kant = kdim;
  
  
  /* For each dimension in all boxes perform box-test.  */	
  
  for (ki=0; ki<kant; ki++,tmin1++,tmax1++,tmin2++,tmax2++)
    {
      /* Sorting: t1-t2 The SISLbox with largest max value.
	 t3-t4 The other box. */
      
      if (*tmax1 > *tmax2)
	{
	  t1 = *tmax1;
	  t2 = *tmin1;
	  t3 = *tmax2;
	  t4 = *tmin2;
	}
      else
	{
	  t1 = *tmax2;
	  t2 = *tmin2;
	  t3 = *tmax1;
	  t4 = *tmin1;
	}
      
      /* SISLPoint intersection in 3D must have a tolerance on the box. */
      if ((po1->iobj == SISLPOINT || po2->iobj == SISLPOINT) && kdim != 1)
	{
	  t1 += 0.1 *  aepsge;
	  t2 -= 0.1 *  aepsge;
	}
      
      if (t1 - min(t2,t4) <= aepsge)
	kj++;                 /* Minibox is possible. */
      else if (t3 < t2)
	{
	  *jstat = 0;           /* No overlap. */
	  goto out;
	}
      else if (kdim != 1 && t3 - t2 <= aepsge && t3 - t4 > aepsge &&
	       t1 - t2 > aepsge && kbez)
	{
	  *jstat = 2;           /* Only edge touching possible.*/
	  goto degenerate;
	}
      else if ( kdim == 1 && (t1 - t4 <= aepsge || t3 - t2 <= aepsge) &&
	       kbez)
	{
	  *jstat = 2;           /* Only edge touching possible.*/
	  goto out;
	}
      /* else possible overlap. */
    }
  
  if (kj == kant)
    *jstat = 3;                   /* Minibox found. */
  else
    *jstat = 1;                   /* Overlap.  */
  
  
  /* Box-test performed. */
  
 degenerate:
  /* Test if one of the objects has collapsed. */
  if (kdim != 1 && po1->iobj > SISLPOINT)
    {
      if (po1->iobj == SISLCURVE)
	{
	  tmin1 = po1->c1->pbox->emin;
	  tmax1 = po1->c1->pbox->emax;
	}
      
      else if (po1->iobj == SISLSURFACE)
	{
	  tmin1 = po1->s1->pbox->emin;
	  tmax1 = po1->s1->pbox->emax;
	}
      
      for(ki=0;ki<kdim;ki++,tmin1++,tmax1++)
	if (DNEQUAL(*tmin1,*tmax1)) break;
      
      if (ki == kdim)
	{
	  *jstat = 4;
	  goto out;
	}
    }
  
  if (kdim != 1 && po2->iobj > SISLPOINT)
    {
      if (po2->iobj == SISLCURVE)
	{
	  tmin1 = po2->c1->pbox->emin;
	  tmax1 = po2->c1->pbox->emax;
	}
      
      else if (po2->iobj == SISLSURFACE)
	{
	  tmin1 = po2->s1->pbox->emin;
	  tmax1 = po2->s1->pbox->emax;
	}
      
      for(ki=0;ki<kdim;ki++,tmin1++,tmax1++)
	if (DNEQUAL(*tmin1,*tmax1)) break;
      
      if (ki == kdim)
	{
	  *jstat = 4;
	  goto out;
	}
    }
  goto out;
  
  /* Dimensions conflicting. */
  
 err106: *jstat = -106;
  s6err("s1790",*jstat,kpos);
  goto out;
  
  /* Dimensions less than one. */
  
 err105: *jstat = -105;
  s6err("s1790",*jstat,kpos);
  goto out;
  
  /* Kind of object does not exist. */
  
 err121: *jstat = -121;
  s6err("s1790",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine. */
  
 error:  *jstat = kstat;
  s6err("s1790",*jstat,kpos);
  goto out;
  
 out:	return;
}





