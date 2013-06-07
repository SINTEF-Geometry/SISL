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
 * $Id: s2535.c,v 1.3 2006-05-02 15:06:04 sbr Exp $
 *
 */

#define S2535

#include "sislP.h"



#if defined(SISLNEEDPROTOTYPES)
void
   s2535(SISLSurf *surf, int u_continuity, int v_continuity, int *u_surfnumb, 
	 int *v_surfnumb, SISLSurf ***patches, int *stat)
#else
void s2535(surf, u_continuity, v_continuity, u_surfnumb, v_surfnumb, patches,
	   stat)
     SISLSurf *surf;
     int      u_continuity;
     int      v_continuity;
     int      *u_surfnumb;
     int      *v_surfnumb;
     SISLSurf ***patches;
     int      *stat;
#endif
/*
***************************************************************************
*
***************************************************************************
* PURPOSE : To split a surface in order to meet given continuity requirements.
*
*
* INPUT   : surf         - The original surface.
*           u_continuity - Desired continuity of the surface in u direction:
*                          = 0 : Positional continuity,
*                          = 1 : Tangential continuity,
*                          and so on.
*                          SISL only accepts surfaces of continuity 0 or larger.
*                          If the surface is to be intersected with another,
*                          the continuity must be 1 or larger to find all the
*                          intersection curves.
*           v_continuity - Desired continuity of the surface in v direction:
*                          = 0 : Positional continuity,
*                          = 1 : Tangential continuity,
*                          and so on.
*                          SISL only accepts surfaces of continuity 0 or larger.
*                          If the surface is to be intersected with another,
*                          the continuity must be 1 or larger to find all the
*                          intersection curves.
*
*
* OUTPUT  : u_surfnumb   - Number of surface patches in u direction.
*           v_surfnumb   - Number of surface patches in v direction.
*           patches      - Array of patches with the given continuity. The
*                          array index runs fastest in the u direction.
*                          If u_surfnumb and v_surfnumb are both 1, no split
*                          is necessary and patches is empty.
*           stat         - Status message.
*                          > 0      : Warning.
*                          = 0      : Ok.
*                          < 0      : Error.
*
*
* METHOD  :
*
*
* CALLS   : s1711()
*
* WRITTEN BY : Johannes Kaasa, SINTEF, Oslo, Norway, Aug. 1995.
*
*********************************************************************
*/
{
   int ki, kj;        /* Array indexes.                                   */
   int mult;          /* Knot multiplicity.                               */
   int u_mult;        /* Knot multiplicity giving a split in u direction. */
   int v_mult;        /* Knot multiplicity giving a split in v direction. */
   double *u_splitpar = SISL_NULL; /* Split parameters in u direction.         */
   double *v_splitpar = SISL_NULL; /* Split parameters in v direction.         */
   SISLSurf *u_surf;  /* Surface pointer used in the surface splitting.   */
   SISLSurf *u_lsurf; /* Surface pointer used in the surface splitting.   */
   SISLSurf *u_rsurf; /* Surface pointer used in the surface splitting.   */
   SISLSurf *v_surf;  /* Surface pointer used in the surface splitting.   */
   SISLSurf *v_lsurf; /* Surface pointer used in the surface splitting.   */
   SISLSurf *v_rsurf; /* Surface pointer used in the surface splitting.   */
   
   /* Check input. */
   
   if (surf == SISL_NULL || u_continuity < 0 || v_continuity < 0) 
      goto err150;
   
   /* Initiation and allocation of split parameters. */
   
   u_mult = max(surf->ik1 - u_continuity, 1);
   v_mult = max(surf->ik2 - v_continuity, 1);
   
   if ((u_splitpar = newarray(((int) floor((double) surf->in1/u_mult)) - 1, DOUBLE)) 
       == SISL_NULL) goto err101;
   if ((v_splitpar = newarray(((int) floor((double) surf->in2/v_mult)) - 1, DOUBLE)) 
       == SISL_NULL) goto err101;
   
   *u_surfnumb = 0;
   *v_surfnumb = 0;
   
   /* Generate the split points. */
   
   ki = surf->ik1;
   while(ki < surf->in1)
   {
      mult = 1;
      kj = ki;
      while (DEQUAL(surf->et1[kj], surf->et1[ki]))
      {
	 mult++;
	 kj++;
      }
      if (mult > u_mult)
      {
	 u_splitpar[*u_surfnumb] = surf->et1[ki];
	 *u_surfnumb += 1;
      }
      ki = kj;
   }
   
   ki = surf->ik2;
   while(ki < surf->in2)
   {
      mult = 1;
      kj = ki;
      while (DEQUAL(surf->et2[kj], surf->et2[ki]))
      {
	 mult++;
	 kj++;
      }
      if (mult > v_mult)
      {
	 v_splitpar[*v_surfnumb] = surf->et2[ki];
	 *v_surfnumb += 1;
      }
      ki = kj;
   }   


   *u_surfnumb += 1;
   *v_surfnumb += 1;
   
   if (*u_surfnumb == 1 && *v_surfnumb == 1)
      goto out;
   
   /* Allocate the output array. */
   
   if ((*patches = newarray((*u_surfnumb)*(*v_surfnumb), SISLSurf*)) 
       == SISL_NULL) goto err101;
   
   /* Split the surfaces. */
   
   v_surf = surf;
   for (ki = 0; ki < (*v_surfnumb - 1); ki++)
   {
      
      /* Split in v direction. */
      
      s1711(v_surf, 2, v_splitpar[ki], &v_lsurf, &v_rsurf, stat);
      if (*stat < 0) goto error;
      
      if (v_surf != surf && v_surf != SISL_NULL) freeSurf(v_surf);
      v_surf = v_rsurf;
      
      /* Split in u direction. */
      
      u_surf = v_lsurf;
      for (kj = 0; kj < (*u_surfnumb - 1); kj++)
      {
	 s1711(u_surf, 1, u_splitpar[kj], &u_lsurf, &u_rsurf, stat);
	 if (*stat < 0) goto error;
	 
	 if (u_surf != SISL_NULL) freeSurf(u_surf);
         u_surf = u_rsurf;
	 
	 (*patches)[ki*(*u_surfnumb) + kj] = u_lsurf;
      }
      
      /* Take care of the last column. */
      
      (*patches)[ki*(*u_surfnumb) + (*u_surfnumb - 1)] = u_surf;
   }
   
   /* Split the last row. */
   
   u_surf = v_surf;
   for (kj = 0; kj < (*u_surfnumb - 1); kj++)
   {
      s1711(u_surf, 1, u_splitpar[kj], &u_lsurf, &u_rsurf, stat);
      if (*stat < 0) goto error;
      
      if (u_surf != surf && u_surf != SISL_NULL) freeSurf(u_surf);
      u_surf = u_rsurf;
      
      (*patches)[(*v_surfnumb - 1)*(*u_surfnumb) + kj] = u_lsurf;
   }
   
   /* Take care of the last column. */
   
   (*patches)[(*v_surfnumb)*(*u_surfnumb) - 1] = u_surf;

   goto out;



   /* ---------------------- ERROR EXITS ------------------------------- */

   /* Error in space allocation */
   
 err101: 
   *stat = -101;
   s6err("s2535",*stat,0);
   goto out;

   /* Error in input. */
   
 err150:
   *stat = -150;
   s6err("s2535", *stat, 0);
   goto out;
   
   /* Error in lower level routine. */
   
 error:
   s6err("s2535", *stat, 0);
   goto out;

   /* ---------------------- NORMAL EXIT ------------------------------- */

 out:
   if (u_splitpar != SISL_NULL) freearray(u_splitpar);
   if (v_splitpar != SISL_NULL) freearray(v_splitpar);
   return;

}
