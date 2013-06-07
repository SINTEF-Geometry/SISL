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
 * $Id: s2536.c,v 1.3 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S2536

#include "sislP.h"



#if defined(SISLNEEDPROTOTYPES)
void
   s2536(SISLSurf *surf, int u_continuity, int v_continuity, 
      int *u_surfnumb, int *v_surfnumb, SISLSurf ***mehlum_surf, int *stat)
#else
void s2536(surf, u_continuity, v_continuity, u_surfnumb, 
	   v_surfnumb, mehlum_surf, stat)
     SISLSurf *surf;
     int      u_continuity;
     int      v_continuity;
     int      *u_surfnumb;
     int      *v_surfnumb;
     SISLSurf ***mehlum_surf;
     int      *stat;
#endif
/*
***************************************************************************
*
***************************************************************************
* PURPOSE : To interpolate the Mehlum curvature of a B-SPLINE or NURBS 
*           surface, as a NURBS surface. 
*           The desired continuity of the Mehlum-curvature surface is given as 
*           input, this may lead to a patchwork of output surfaces. 
*           THE INTERPOLATION RESULTS IN A HIGH ORDER SURFACE (IF THE ORIGINAL 
*           SURFACE IS A B-SPLINE OF ORDER K, THE RESULT IS OF ORDER 12K - 17, 
*           IN THE NURBS CASE THE ORDER IS 48K - 53). TO AVOID UNSTABILITY
*           BECAUSE OF THIS, A MAX. ORDER IS APPLIED. THIS MAY LEAD TO AN
*           APPROXIMATION INSTEAD OF AN INTERPOLATION.
*
*
* INPUT   : surf         - The original surface.
*           u_continuity - Desired continuity of the Mehlum-curvature surfaces
*                          in u direction:
*                          = 0 : Positional continuity,
*                          = 1 : Tangential continuity,
*                          and so on.
*                          SISL only accepts surfaces of continuity 0 or larger.
*                          If the surface is to be intersected with another,
*                          the continuity must be 1 or larger to find all the
*                          intersection curves.
*           v_continuity - Desired continuity of the Mehlum-curvature surfaces
*                          in v direction:
*                          = 0 : Positional continuity,
*                          = 1 : Tangential continuity,
*                          and so on.
*                          SISL only accepts surfaces of continuity 0 or larger.
*                          If the surface is to be intersected with another,
*                          the continuity must be 1 or larger to find all the
*                          intersection curves.
*
*
* OUTPUT  : u_surfnumb   - Number of Mehlum-curvature surface patches
*                          in u direction.
*           v_surfnumb   - Number of Mehlum-curvature surface patches
*                          in v direction.
*           mehlum_surf  - The Mehlum curvature interpolation surfaces.
*                          This will be a pointer to an array of length
*                          u_surfnum*v_surfnumb of SISLSurf pointers,
*                          where the indexing is running fastest in the
*                          u direction.
*           stat         - Status message.
*                          > 0      : Warning.
*                          = 2      : The surface is degenerated.
*                          = 0      : Ok.
*                          < 0      : Error.
*
*
* METHOD  :
*
*
* CALLS   : s2535(), s2534()
*
* WRITTEN BY : Johannes Kaasa, SINTEF, Oslo, Norway, Aug. 1995.
*
*********************************************************************
*/
{
   int ki;                     /* Index in for loops.                    */
   int u_multinc;              /* Increased knot multiplicity.           */
   int v_multinc;              /* Increased knot multiplicity.           */
   int newik1;                 /* New order in u direction.              */
   int newik2;                 /* New order in v direction.              */
   int eval_dim;               /* Evaluation dimention.                  */
   int max_order = 20;         /* Max. order of the curvature surface.   */
   SISLSurf *temp = SISL_NULL;      /* Temporary surface.                     */
   SISLSurf *regular = SISL_NULL;   /* k-regular surface.                     */
   SISLSurf **org_surf = SISL_NULL; /* Array of pointers to original patches. */
   
   
   /* Check input. */
   
   if (surf == SISL_NULL || u_continuity < 0 || v_continuity < 0)
      goto err150;
   
   /* Curvature continuity decrease. */
   
   u_multinc = 2;
   v_multinc = 2;
   
   /* Make sure we have a k-regular surface. */

   if (surf->cuopen_1 == SISL_SURF_PERIODIC || 
       surf->cuopen_2 == SISL_SURF_PERIODIC )
   {
      make_sf_kreg(surf, &regular, stat);
      if (*stat < 0) goto error;
   }
   else
      regular = surf;
   
   /* Split the surface to meet the continuity requirements. */
   
   s2535(regular, (u_continuity + u_multinc), (v_continuity + v_multinc), 
	 u_surfnumb, v_surfnumb, &org_surf, stat);
   if (*stat < 0) goto error;
   
   /* Allocate output array. */
   
   if ((*mehlum_surf = newarray((*u_surfnumb)*(*v_surfnumb), SISLSurf*)) 
       == SISL_NULL) goto err101;
   for (ki = 0; ki < (*u_surfnumb)*(*v_surfnumb); ki++)
      (*mehlum_surf)[ki] = SISL_NULL;
   
   /* Calculate curvature order. */
   
   if (regular->ikind == 1 || regular->ikind == 3)
   {
      newik1 = 12*regular->ik1 - 17;
      newik2 = 12*regular->ik2 - 17;
   }
   else
   {
      newik1 = 48*regular->ik1 - 53;
      newik2 = 48*regular->ik2 - 53;
   }
   newik1 = min(newik1, max_order);
   newik2 = min(newik2, max_order);
   
   eval_dim = 2;
   
   /* Interpolate. */
   
   if (*u_surfnumb == 1 && *v_surfnumb == 1)
   {
      s2534(regular, u_multinc, v_multinc, newik1, newik2, s2515, eval_dim,
	    &temp, stat);
      if (*stat < 0) goto error;
      if (*stat == 2) goto war002;
      
      (*mehlum_surf)[0] = newSurf(temp->in1, temp->in2, temp->ik1, temp->ik2,
				 temp->et1, temp->et2, temp->ecoef, 2, 1, 1);
      
      if (temp != SISL_NULL)
      {
	 freeSurf(temp);
	 temp = SISL_NULL;
      }
   }
   else
   {
      for (ki = 0; ki < (*u_surfnumb)*(*v_surfnumb); ki++)
      {
	 s2534(org_surf[ki], u_multinc, v_multinc, newik1, newik2, s2515, 
	       eval_dim, &temp, stat);
         if (*stat < 0) goto error;
	 if (*stat == 2) goto war002;
      
         (*mehlum_surf)[ki] = newSurf(temp->in1, temp->in2, temp->ik1, 
				     temp->ik2, temp->et1, temp->et2, 
				     temp->ecoef, 2, 1, 1);
	 
	 if (temp != SISL_NULL)
         {
	    freeSurf(temp);
	    temp = SISL_NULL;
         }
      }
   }


   goto out;



   /* ---------------------- ERROR EXITS ------------------------------- */
   
   /* The surface is degenerated at (u,v) */
   
 war002:
   if (*mehlum_surf != SISL_NULL)
   {
      for (ki = 0; ki < ((*u_surfnumb)*(*v_surfnumb)); ki++)
	 if ((*mehlum_surf)[ki] != SISL_NULL) freeSurf((*mehlum_surf)[ki]);
      freearray(*mehlum_surf);
   }
   *u_surfnumb = 0;
   *v_surfnumb = 0;
   goto out;
  
  /* Error in space allocation */
   
 err101: 
   *stat = -101;
   s6err("s2536", *stat, 0);
   goto out;

   /* Error in input. */
   
 err150:
   *stat = -150;
   s6err("s2536", *stat, 0);
   goto out;

   /* Error in lower level routine. */
   
 error:
   s6err("s2536", *stat, 0);
   goto out;

   /* ---------------------- NORMAL EXIT ------------------------------- */

 out:
   if (regular != surf ) freeSurf(regular);
   if (org_surf != SISL_NULL)
   {
      for (ki = 0; ki < ((*u_surfnumb)*(*v_surfnumb)); ki++)
	 if (org_surf[ki] != SISL_NULL) freeSurf(org_surf[ki]);
      freearray(org_surf);
   }
   
   return;

}
