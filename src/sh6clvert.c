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
 * $Id: sh6clvert.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */
#define SH6CLOSEVERT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6closevert(SISLCurve *pcurve,SISLSurf *psurf,double *cpar1,
	     double epar2[])
#else
void sh6closevert(pcurve,psurf,cpar1,epar2)
   SISLCurve *pcurve;
   SISLSurf *psurf;
   double *cpar1;
   double epar2[];
#endif
/*
*********************************************************************
*
* PURPOSE    : Estimate the parameter values the closest vertices of
*              a B-spline curve and a B-spline surface.
*
*
*
* INPUT      : pcurve   - Pointer to curve.
*              psurf    - Pointer to surface.
*
*
* OUTPUT     : cpar1    - Parameter value of closest vertex of the curve.
*              epar2    - Parameter value of closest vertex of the surface.
*
*
* METHOD     : Find the vertices of the curve and the surface that are
*              closest to each other. Regard these vertices as Schoenberg
*              points, and compute the parameter values.
*
*
* REFERENCES :
*
*
* USE        :
*
*-
* CALLS      : s6dist - Distance between two points.
*
*
* WRITTEN BY : Vibeke Skytt, SI, 06.92.
*
*********************************************************************
*/
{
  int ki,kj,kl;           /* Counters.   */
  int kdim = pcurve->idim; /* Dimension of geometry space.      */
  int kminc;              /* Number of closest vertex of curve. */
  int kmins1,kmins2;      /* Numbers of closest vertex of surface. */
  int kn = pcurve->in;    /* Number of coefficients of curve. */
  int kn1 = psurf->in1;   /* Number of coefficients of surface, 1. par. dir. */
  int kn2 = psurf->in2;   /* Number of coefficients of surface, 2. par. dir. */
  int kk = pcurve->ik;    /* Order of curve. */
  int kk1 = psurf->ik1;   /* Order of surface, 1. par. dir. */
  int kk2 = psurf->ik2;   /* Order of surface, 2. par. dir. */
  double tdist;           /* Distance.   */
  double tmin = HUGE;     /* Minimum distance.  */
  double tpar;            /* Used to compute parameter values.   */
  double *s1,*s2;         /* Pointers into arrays.   */

  /* Find position of closest vertices. */

  for (s1=pcurve->ecoef, ki=0; ki<kn; s1+=kdim, ki++)
    for (s2=psurf->ecoef, kj=0; kj<kn1; kj++)
      for (kl=0; kl<kn2; s2+=kdim, kl++)
	{
	   tdist = s6dist(s1,s2,kdim);
	   if (tdist < tmin)
	   {
	      tmin = tdist;
	      kminc = ki;
	      kmins1 = kj;
	      kmins2 = kl;
	   }
	}

  /* Estimate parameter values of vertices.  */

  for (ki=kminc+1, s1=pcurve->et+ki, tpar=DZERO;
   ki<kminc+kk; tpar+=(*s1), s1++, ki++);
  *cpar1 = tpar/(double)(kk-1);

  for (ki=kmins1+1, s1=psurf->et1+ki, tpar=DZERO;
   ki<kmins1+kk1; tpar+=(*s1), s1++, ki++);
  epar2[0] = tpar/(double)(kk1-1);

  for (ki=kmins2+1, s1=psurf->et2+ki, tpar=DZERO;
   ki<kmins2+kk2; tpar+=(*s1), s1++, ki++);
  epar2[1] = tpar/(double)(kk2-1);


  goto out;

 out:
  return;
}
