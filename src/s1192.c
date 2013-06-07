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
 * $Id: s1192.c,v 1.2 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1192

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/
#if defined(SISLNEEDPROTOTYPES)
static void s1192_s9mbox(double [],int,int,double,double *,double *,
			 int *,int *);
#else
static void s1192_s9mbox();
#endif
  
#if defined(SISLNEEDPROTOTYPES)
void 
s1192(SISLObject *po,double aepsge,int *jstat)
#else
void s1192(po,aepsge,jstat)
     SISLObject *po;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make a SISLbox on the control-polygons given by
*              ecoef to the object. If dimension is 2 then one
*	       SISLbox rotated 45 degree is also made. If dimension
*	       is 3 then one SISLbox roteted 45 degree around each
*	       main axes, in all 4 boxes is made.
*
*
*
* INPUT      : po     - SISLObject to treat.
*              aepsge - Geometry resolution.
*
* OUTPUT     : jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                                                     
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
* WRITTEN BY : Ulf J. krystad, SI, 89-06.
*
*********************************************************************
*/                                     
{
  int kpos = 0;                        /* Position of error.   */
  
  
  if (po -> iobj == SISLPOINT)
    {
      if (po->p1->idim != 1) goto err105;
      
      if (po->p1->pbox == SISL_NULL)
	{
	  if ((po->p1->pbox = newbox(po->p1->idim))==SISL_NULL)
	    goto err101;
	  
	  s1192_s9mbox(po->p1->ecoef,1,1,aepsge,
		 po->p1->pbox->emax,po->p1->pbox->emin,
		 &po->p1->pbox->imax,&po->p1->pbox->imin);
	  
	  
	}
    }
  else
    if (po -> iobj == SISLCURVE)
      {
	if (po->c1->idim != 1) goto err105;
	if (po->c1->pbox == SISL_NULL)
	  {
	    if ((po->c1->pbox = newbox(po->c1->idim))==SISL_NULL)
	      goto err101;
	    
	    s1192_s9mbox(po->c1->ecoef,po->c1->in,1,aepsge,
		   po->c1->pbox->emax,po->c1->pbox->emin,
		   &po->c1->pbox->imax,&po->c1->pbox->imin);
	    
	  }
      }
    else
      if (po -> iobj == SISLSURFACE)
	{
	  if (po->s1->idim != 1) goto err105;
	  if (po->s1->pbox == SISL_NULL)
	    {
	      if ((po->s1->pbox = newbox(po->s1->idim))==SISL_NULL)
		goto err101;
	      
	      s1192_s9mbox(po->s1->ecoef,po->s1->in1,po->s1->in2,aepsge,
		     po->s1->pbox->emax,po->s1->pbox->emin,
		     &po->s1->pbox->imax,&po->s1->pbox->imin);
	    }
	}
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1192",*jstat,kpos);
  goto out;
  
  /* Dimension not equal one.  */
  
 err105: *jstat = -105;
  s6err("s1192",*jstat,kpos);
  goto out;
  
 out: ;
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1192_s9mbox(double ecoef[], int in1,int in2,double aepsge,
		   double *cmax, double *cmin,int *jmax,int *jmin)
#else
static void s1192_s9mbox(ecoef,in1,in2,aepsge,cmax,cmin,jmax,jmin)
     double ecoef[];
     int    in1;
     int    in2;
     double aepsge;
     double *cmax;
     double *cmin;
     int    *jmin;
     int    *jmax;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make a SISLbox on the control-polygons given by
*              ecoef to the object, one-dimensional case.
*
*
* INPUT      : ecoef  - Control-polygon.
*              in1    - Number of vertices in control-polygon 1. direction.
*              in2    - Number of vertices in control-polygon 2. direction.
*              aepsge - The geometry resolution.
*                                                                     
*
* OUTPUT     : cmax   - Maximum value of the box.
*	       cmin   - Minimum value of the box.
*	       jmax   - The ecoef index of the maximum value.
*	       jmin   - The ecoef index of the minimum value.
*                       NB : If jmin,jmax has the value of a corner, it is 
*                       guaranteed that no inner vertice is close(within aepsge)
*                       to the max (min) value.

*
* METHOD     : We treat the corners and the inner points separately and get
*              the result by comparing max(min) corners to ditto inners.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SI, 89-06.
*
*********************************************************************
*/                                     
{
  int ki,kj,li[4];         /* Counters.  */
  int icorn;               /* Number of corners in object.  */
  int kmin, kmax;          /* Index for max and min corner value.  */
  double tmin, tmax;       /* Max and min corner value.  */
  
  /* Compute the indexes of the (up to four) corners. */
  li[0] = 0;
  li[1] = in1 -1;
  li[2] = in1*(in2 - 1);
  li[3] = in1*in2 - 1;
  
  /* Set number of corners. 
     for point, curve, surface. */
  if(in1 == 1)
    {
      if(in2 == 1) 
	icorn = 0;
      else
	icorn = 2;
    }
  else
    icorn = 4;
  
  /* Now find the max and min corner. */
  tmax = tmin = ecoef[li[0]];
  kmin = kmax = 0;
  
  for (ki = 1; ki < icorn; ki++)
    {
      if (ecoef[li[ki]] > tmax)
	{
	  tmax = ecoef[li[ki]];
	  kmax = li[ki];
	}
      
      if (ecoef[li[ki]] < tmin)
	{
	  tmin = ecoef[li[ki]];
	  kmin = li[ki];
	}
    }
  
  /* Now find the max and min for the inner of the object. */
  *cmax = tmax - (double)1000.0;
  *jmax = -1;
  *cmin = tmin + (double)1000.0;
  *jmin = -1;
  
  for (ki = 0; ki < icorn - 1; ki++)
    for (kj = li[ki] + 1; kj < li[ki + 1]; kj++)
      {
	if (ecoef[kj] > *cmax)
	  {
	    *cmax = ecoef[kj];
	    *jmax = kj;
	  }
	
	if (ecoef[kj] < *cmin)
	  {
	    *cmin = ecoef[kj];
	    *jmin = kj;
	  }
      }
  
  
  /* At last compare the corner values against the interior ones */
  
  if (tmax > *cmax + aepsge)
    { 
      *cmax = tmax;
      *jmax = kmax;
    }
  
  if (tmin < *cmin - aepsge)
    { 
      *cmin = tmin;
      *jmin = kmin;
    }
}
