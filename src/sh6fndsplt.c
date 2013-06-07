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
 * $Id: sh6fndsplt.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6FINDSPLIT

#include "sislP.h"

/* 
 extern int nmbcall;
extern int nmb0;
extern int nmb1;
extern int nmb2;
extern int nmb3;
extern int nmb4;
extern int nmbsuccess; 
*/

#if defined(SISLNEEDPROTOTYPES)
void
sh6findsplit (SISLSurf *ps1, SISLSurf *ps2, double aepsge, int *jstat)
#else
void
sh6findsplit (ps1, ps2, aepsge, jstat)
   SISLSurf *ps1;
   SISLSurf *ps2;
   double aepsge;
   int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Try to intercept intersection between two B-spline
*              surfaces by finding a splitting geometry between the
*              two surfaces. This geometry is a plane, a sphere, a
*              cylinder or a torus.
*
*
* INPUT      : ps1      - First surface.
*              ps2      - Second surface.
*              aepsge   - Geometry tolerance. 
*
*
* OUTPUT     : jstat    - status messages
*                                = 1   : Intersection still possible.
*                                = 0   : Ok. No intersection is possible.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : 
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SI, 93-05.
*
*********************************************************************
*/
{
   int kstat = 0;   /* Local status variable.  */
   int kdim = ps1->idim;  /* Dimension of space. */
   int ratflag = 0; /* Indicates if surface is rational. */
   double tepsge;   /* Local tolerance.        */
   double tdist;    /* Large radius of torus.  */
   double trad;     /* Radius of sphere, cylinder or torus. */
   double scentre[3];  /* Centre of splitting geometry.     */
   double saxis[3];    /* Axis of splitting geometry.       */
   double simpli[16];   /* Array containing torus info.      */
   double splitgeom[16];         /* Matrix description of a sphere
				    or cylinder.                   */
   SISLSurf *qs1 = SISL_NULL; /* 1D surface.                     */
   SISLSurf *qs2 = SISL_NULL; /* 1D surface.                     */
      
   /* Still overlap. Try to find splitting geometry object. */
   
   sh6splitgeom(ps1, ps2, aepsge, scentre, saxis, &tdist,
		&trad, &kstat);
   if (kstat < 0) goto error;
   
   /* 
   if (kstat == 0) nmb0++;
   else if (kstat == 1) nmb1++;
   else if (kstat == 2) nmb2++; 
   else if (kstat == 3) nmb3++; 
   else if (kstat == 4) nmb4++; 
 */
   
   /* If kstat = 0 is returned, no splitting geometry is found,
      and no further interception is to be tried.  */
   
   if (kstat > 0)
   {
      if (kstat == 1)
      {
	 /* The splitting geometry is a plane. Set the two surfaces
	    into the plane equation.  */
	 
	 s1329 (ps1, scentre, saxis, kdim, &qs1, &kstat);
	 if (kstat < 0)
	    goto error;
	 s1329 (ps2, scentre, saxis, kdim, &qs2, &kstat);
	 if (kstat < 0)
	    goto error;
	 
	 
	 /* Set local tolerance.  */
	 
	 tepsge = aepsge;
      }
      else if (kstat == 2 || kstat == 3)
      {
	 if (kstat == 2)
	 {
	    /* The splitting geometry object is a sphere.  
	       Make a matrix of dimension (idim+1)x(idim+1) describing a hyper
	       sphere as an implicit function.      	      */
	    
	    s1321(scentre,trad,kdim,1,splitgeom,&kstat);
	    if (kstat < 0) goto error;
	    
	 }
	 else if (kstat == 3)
	 {
	    /* The splitting geometry object is a cylinder.
	       Make a matrix of dimension (idim+1)x(idim+1) describing a 
	       cylinder as an implicit function.           */
	    
	    s1322(scentre,saxis,trad,kdim,1,splitgeom,&kstat);
	    if (kstat < 0) goto error;
	 }
	 /* 
	 * Put the description of the surfaces into the implicit
	 * equation for the sphere or cylinder.
	 * ----------------------------------------------------------
	 */
	 
	 ratflag = (ps1->ikind == 2 || ps1->ikind == 4) ? 1 : 0;
	 s1320(ps1,splitgeom,1,ratflag,&qs1,&kstat);
	 if (kstat < 0) goto error;
	 
	 ratflag = (ps2->ikind == 2 || ps2->ikind == 4) ? 1 : 0;
	 s1320(ps2,splitgeom,1,ratflag,&qs2,&kstat);
	 if (kstat < 0) goto error;
	 
	 /* Set up local tolerance. */
	 
	 tepsge = (double)2.0*trad*aepsge;
      }
      else if (kstat == 4)
      {
	 /* Set surfaces into torus equation. */
	 
	 /* 
	 * Put the information concerning the torus in the following sequence
	 * into simpli: Center, normal, big radius, small radius 
	 * ------------------------------------------------------------------
	 */
	 
	 memcopy(simpli,scentre,3,DOUBLE);
	 memcopy(simpli+3,saxis,3,DOUBLE);
	 simpli[6] = tdist;
	 simpli[7] = trad;
	 
	 /* 
	 * Put surfaces into torus equation 
	 * -------------------------------
	 */ 
	 
	 s1378(ps1,simpli,1001,kdim,&qs1,&kstat);
	 if (kstat<0) goto error;
	 
	 s1378(ps2,simpli,1001,kdim,&qs2,&kstat);
	 if (kstat<0) goto error;
	 
	 /* Set up local tolerance. */
	 
	 tepsge = (double)8.0*aepsge*trad*tdist*tdist;
      }	 
	 
      
      /* Make box of first 1D surface. */
      
      sh1992su(qs1,2,tepsge,&kstat);
      if (kstat < 0) goto error;
      
      /* Make box of second 1D surface. */
      
      sh1992su(qs2,2,tepsge,&kstat);
      if (kstat < 0) goto error;
      
      /* Check if the boxes overlap.  */
      
      if (qs1->pbox->e2min[2][0] > qs2->pbox->e2max[2][0] ||
	  qs1->pbox->e2max[2][0] < qs2->pbox->e2min[2][0])
      {
	 /* 
	 nmbsuccess++; 
	 */
	 
	 /* No intersection is possible.  */
	 
	 *jstat = 0;
      }
      else *jstat = 1;  /* Mark possibility of intersection.  */
   }
   else *jstat = 1;  /* Mark possibility of intersection.  */

   goto out;
   
   error : *jstat = kstat;
   goto out;
   
   out:
      if (qs1) freeSurf(qs1);
      if (qs2) freeSurf(qs2);
      
      return;
}
