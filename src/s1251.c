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
 * $Id: s1251.c,v 1.3 2001-03-19 15:58:43 afr Exp $
 *
 */

#define S1251

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1251(SISLCurve *pcurve,double aepsco,double *clength,int *jstat)
#else
void s1251(pcurve,aepsco,clength,jstat)
     SISLCurve  *pcurve;
     double aepsco;
     double *clength;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Calculate the length of a B-spline curve. The length
*              calculated will not deviate more than aepsco from the
*              real length of the curve.
*
*
*
* INPUT      : pcurve - Pointer to curve.
*              aepsco - Computer resolution.
*
*
*
* OUTPUT     : clength - The length of the curve.
*              jstat   - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6dist,s1251,s1710,freeCurve.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* CHANGED BY : Vibeke Skytt, SINTEF Oslo, 95-05. Introduced stop criterion
*                            on absolute lenght of the curve, and subdivided
*                            initially into Bezier curves.
* CORRECTED BY : Ulf J Krystad, SINTEF Oslo, 96-03. Corrected no of bez seg.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int ki;                   /* Counter.                                    */
  int kn;                   /* Number of vertices of curve.                */
  int kk;                   /* Order of curve.                             */
  int kdim;                 /* Dimension of the space in which the curve lies.*/
  int knbez;                /* Number of Bezier segments of curve.         */
  double tpol = (double)0.0;/* Length of control polygon.                  */
  double tdist =(double)0.0;/* Distance between first and last vertex.     */
  double tdum;              /* Help variable.                              */
  double tmid;              /* Midpoint of parameter interval of curve.    */
  double tlength1,tlength2; /* Length of sub-curves.                       */
  double *s1;               /* Pointer used to traverse coefficient array. */
  SISLCurve *qc1 = SISL_NULL;        /* First sub-curve.                            */
  SISLCurve *qc2 = SISL_NULL;        /* Second sub-curve.                           */

  /* Copy properties of curve to local parameters. */

  kn = pcurve -> in;
  kk = pcurve -> ik;
  kdim = pcurve -> idim;

  /* Calculate length of control polygon. */

  for (ki=1,s1=pcurve->ecoef+kdim; ki<kn; ki++,s1+=kdim)
    tpol += s6dist(s1-kdim,s1,kdim);

  /* Calculate distance from first to last vertex.  */

  tdist = s6dist(pcurve->ecoef,pcurve->ecoef+(kn-1)*kdim,kdim);

  /* Test if the length of the curve can be approximated by the
     distance from the first to the last vertex.                */

  if (DEQUAL(tpol + tdist,(double)0.0)) tdum = (double)0.0;
  else tdum = (tpol - tdist)/(tpol + tdist);

  if (tdum < aepsco)

    /* The length of the curve is found.  */

    *clength = tdist;
  else if ((tdist <= REL_COMP_RES && tpol <= (double)10*REL_COMP_RES) ||
	   tpol-tdist <= REL_COMP_RES)

     /* Do not subdivide any further. */

     *clength = (double)0.5*(tdist+tpol);
  else if (pcurve->ik == pcurve->in)
    {

      /* Subdivide the curve at the midpoint. */

      tmid = ((double)0.5)*(*(pcurve->et+pcurve->ik-1) + *(pcurve->et+kn));

      s1710(pcurve,tmid,&qc1,&qc2,&kstat);
      if (kstat < 0) goto error;

      /* Compute length of first sub-curve. */

      s1251(qc1,aepsco,&tlength1,&kstat);
      if (kstat < 0) goto error;

      /* Compute length of second sub-curve.  */

      s1251(qc2,aepsco,&tlength2,&kstat);
      if (kstat < 0) goto error;

      *clength = tlength1 + tlength2;
    }
  else
  {
     /* Make a sequence of Bezier curves. */

     s1730(pcurve, &qc1, &kstat);
     if (kstat < 0) goto error;

     /* Pick next Bezier curve and compute the length of the curve. */

     /* UJK, 960329, ref BEOrd29367 */
     /* knbez = qc1->in/kk - 1;     */

     knbez = qc1->in/kk;

     tlength1 = DZERO;
     for (ki=0; ki<knbez; ki++)
     {
	/* Represent the current segment as a curve. */

	if (qc1->ikind == 1 || qc1->ikind == 3)
	   qc2 = newCurve(kk, kk, qc1->et+ki*kk, qc1->ecoef+ki*kk*kdim,
			  qc1->ikind, kdim, 0);
	else
	   qc2 = newCurve(kk, kk, qc1->et+ki*kk, qc1->rcoef+ki*kk*(kdim+1),
			  qc1->ikind, kdim, 0);

	if (qc2 == SISL_NULL) goto err101;

      /* Compute length of the current sub-curve.  */

      s1251(qc2,aepsco,&tlength2,&kstat);
      if (kstat < 0) goto error;

      tlength1 += tlength2;

      if (qc2 != SISL_NULL) freeCurve(qc2);
      qc2 = SISL_NULL;
     }

     *clength = tlength1;
  }

  /* Length of curve found.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation. */

  err101: *jstat = -101;
  s6err("s1251",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("s1251",*jstat,kpos);
  goto out;

 out:

  /* Free space occupied by sub-curves.  */

  if (qc1 != SISL_NULL) freeCurve(qc1);
  if (qc2 != SISL_NULL) freeCurve(qc2);

  return;
}
