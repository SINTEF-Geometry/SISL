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
 * $Id: s1796.c,v 1.3 2001-03-19 15:58:54 afr Exp $
 *
 */
#define S1796

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1796(SISLCurve *pc1,SISLCurve *pc2,double aepsge,double aang,int *jstat)
#else
void s1796(pc1,pc2,aepsge,aang,jstat)
     SISLCurve  *pc1;
     SISLCurve  *pc2;
     double aepsge;
     double aang;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To make an extra simple case test by using the
*	       projection of tangents into a plan span by
*	       the center axes of the direction cones.
*
*
*
* INPUT      : pc1,pc2  - The B-spline curves.
*              aang     - The angel beetween the center axes
*                         in the direction cones.
*              aepsge -   Geometry resolution.
*
*
*
* OUTPUT     : jstat  - status messages  
*                                         > 1      : Simpel case
*                                         = 0      : No simpel case
*                                         < 0      : error
*
*
* METHOD     : We are making a plan span by the center axes af the
*              cones surrounding the orientating surfaces
*	       on the unit sphere. Then we project each tangent
*	       to this plan and compute the angel beetween these
*	       projection and the senter of the cone. If the new
*              projected cones do not overlap we have simel case.
*             
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-07.
*
*********************************************************************
*/
{
  int kstat = 0;     /* Local status variable.                           */
  int kpos = 0;      /* Position of the error.                           */
  int turned = 0;    /* Use as mark if dir of curve2 is turned.		 */
  int kn;            /* Number of vertices of curve.                     */
  int kdim;	     /* Dimension of the space in which the objects lie. */
  int kin;           /* The index to the vertice to treat.               */
  int k1,k2;         /* Control variables in loop.                       */
  double *t=SISL_NULL;    /* Tangent at each coeficient.                      */
  double tlen;       /* The length of a vector.                          */
  double *scen1;     /* The orginal basis vector to the projection plan. */
  double *scen2;     /* The computed basis vector to the projection plan.*/
  double tang;	     /* An angle between two vectors.		         */
  double tang1=DZERO;/* An angle between two vectors.			 */
  double tang2=DZERO;/* An angle between two vectors.			 */
  double t1,t2;      /* Help variables.				         */
  
  
  /* Initialate space dimentions. */
  
  kdim = pc1 -> idim;
  
  
  /* Allocate local used array. */
  
  if ((t = newarray(2*kdim,double)) == SISL_NULL) goto err101;
  
  /* We have to turn the direction into the smallest angel. */
  
  if (aang > PIHALF)
  {
    aang = PI - aang;
    turned = 1;
  }
  
  /* The first basis vector. */
  
  scen1 = pc1->pdir->ecoef;
  
  /* We must orthonormalize the second basis vector. */
  
  scen2 = t + kdim;
  tlen = s6scpr(scen1,pc2->pdir->ecoef,kdim);
  for (k1=0; k1 < kdim; k1++)
    scen2[k1] = pc2->pdir->ecoef[k1] - tlen*scen1[k1];
  tlen = s6length(scen2,kdim,&kstat);
  for (k1=0; k1 < kdim; k1++)
    scen2[k1] /= tlen;
  
  if (turned)  
     for (k1=0; k1 < kdim; k1++)    scen2[k1] = -scen2[k1];
  
  
  /* Here we are treating each patch in the control polygon separately.*/
  
  for (kn=pc1->in,k2=0,kin=0; kin < kn-1; kin++)
    {
      
      /* Here we make an aproximative tangents to the curve
	 using the control polygon. The tangents are also normalized
	 by deviding with its own length. */
      
      for (tlen=DZERO,k1=0; k1 < kdim; k1++,k2++)
	{
	  t[k1] = pc1->pdir->esmooth[k2+kdim] - pc1->pdir->esmooth[k2];
	  tlen += t[k1]*t[k1];
	}
      
      tlen = sqrt(tlen);
      
      if (tlen > aepsge)
	for (k1=0; k1 < kdim; k1++) t[k1] /= tlen;
      else
	for (k1=0; k1 < kdim; k1++) t[k1] = scen1[k1];
      
      t2 = scen2[0]*t[0];
      for (k1=1; k1<kdim; k1++)
	t2 += scen2[k1]*t[k1];
      
      if (t2 <= DZERO) continue;
      
      t1 = scen1[0]*t[0];
      for (k1=1; k1<kdim; k1++)
	t1 += scen1[k1]*t[k1];
      
      tang = t1/sqrt(t1*t1 + t2*t2);
      
      if (tang >= DZERO) tang = min((double)1,tang);
      else               tang = max((double)-1,tang);
      
      tang = acos(tang);
      
      tang1 = max(tang1,tang);
    }
  
  /* The first basis vector. */
  
  scen1 = pc2->pdir->ecoef;
  
  /* We must orthonormalize the second basis vector. */
  
  scen2 = t + kdim;
  tlen = s6scpr(scen1,pc1->pdir->ecoef,kdim);
  for (k1=0; k1 < kdim; k1++)
    scen2[k1] = pc1->pdir->ecoef[k1] - tlen*scen1[k1];
  tlen = s6length(scen2,kdim,&kstat);
  for (k1=0; k1 < kdim; k1++)
    scen2[k1] /= tlen;
  
  if (turned)  
     for (k1=0; k1 < kdim; k1++)    scen2[k1] = -scen2[k1];
  
  /* Here we are treating each patch in the control polygon separately.*/
  
  for (kn =pc2->in,k2=0,kin=0; kin < kn-1; kin++)
    {
      
      /* Here we make an aproximative tangents to the curve
	 using the control polygon. The tangents are also normalized
	 by deviding with its own length. */
      
      for (tlen=DZERO,k1=0; k1 < kdim; k1++,k2++)
	{
	  t[k1] = pc2->pdir->esmooth[k2+kdim] - pc2->pdir->esmooth[k2];
	  tlen += t[k1]*t[k1];
	}
      
      tlen = sqrt(tlen);
      
      if (tlen > aepsge)
	for (k1=0; k1 < kdim; k1++) t[k1] /= tlen;
      else
	for (k1=0; k1 < kdim; k1++) t[k1] = scen1[k1];
      
      
      t2 = scen2[0]*t[0];
      for (k1=1; k1<kdim;k1++)
	t2 += scen2[k1]*t[k1];
      
      if (t2 <= DZERO) continue;
      
      t1 = scen1[0]*t[0];
      for (k1=1; k1<kdim; k1++)
	t1 += scen1[k1]*t[k1];
      
      tang = t1/sqrt(t1*t1 + t2*t2);
      
      if (tang >= DZERO) tang = min((double)1,tang);
      else               tang = max((double)-1,tang);
      
      tang = acos(tang);
      
      tang2 = max(tang2,tang);
    }
  
  
  
  /* Performing a simple case check. */
  
  if (tang1 + tang2 <= aang)
    *jstat = 1;       /* A simpel case.*/
  else
    *jstat = 0;
  
  goto out;
  
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1796",*jstat,kpos);
  goto out;
    
  
  /* Free local used memory. */
  
 out:    if (t != SISL_NULL) freearray(t);
  
}
