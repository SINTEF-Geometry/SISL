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
 * $Id: sh1779_at.c,v 1.2 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH1779_AT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh1779_at (SISLObject * po1, SISLObject * po2, SISLIntpt * pintpt,
	   int *jstat)
#else
void
sh1779_at (po1, po2, pintpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt *pintpt;
     int *jstat;
#endif


/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology AT in 3-dimensional curve-surface
*              intersection.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pintpt   - Current intersection point.
*
*
* OUTPUT     : pintpt   - Intersection point after updating pre-topology.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : sh6gettop  -  Get topology of intersection point.
*              sh6settop  -  Set topology of intersection point.
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, 06.91.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int kpar1, kpar2;		/* Index of parameter value of object.     */
  int kn;			/* Number of vertices of curve.            */
  int kk;			/* Order of curve.                         */
  int lleft[2];			/* Array storing pre-topology information. */
  int lright[2];		/* Array storing pre-topology information. */
  int *ll1, *ll2, *lr1, *lr2;	/* Pointers into pre-topology arrays.   */
  double tref;			/* Referance value in equality test.       */
  double *st;			/* Knot vector of curve.                   */
  double *sptpar = pintpt->epar;/* Pointer to parameter values of int.pt.  */
  SISLCurve *qc;		/* Pointer to the curve.                   */
  SISLSurf *qs;			/* Pointer to the surface.                 */
  double sf_low_lim[2];
  double sf_high_lim[2];
  /* ---------------------------------------------------------------------- */
  /* Don't make pretop for help points ! */
  if (sh6ishelp (pintpt))
    {
      *jstat = 0;
      goto out;
    }

  /* Set pointers into the arrays storing pre-topology information. */
  if (po1->iobj == SISLCURVE)
    {
      qc = po1->c1;
      qs = po2->s1;
      kpar1 = 0;
      kpar2 = 1;
      ll1 = lleft;
      lr1 = lright;
      ll2 = lleft + 1;
      lr2 = lright + 1;
    }
  else
    {
      qc = po2->c1;
      qs = po1->s1;

      kpar1 = 2;
      kpar2 = 0;
      ll1 = lleft + 1;
      lr1 = lright + 1;
      ll2 = lleft;
      lr2 = lright;
    }

  kk = qc->ik;
  kn = qc->in;
  st = qc->et;
  tref = st[kn] - st[kk - 1];

  sf_low_lim[0] = qs->et1[qs->ik1 - 1] + REL_COMP_RES;
  sf_low_lim[1] = qs->et2[qs->ik2 - 1] + REL_COMP_RES;
  sf_high_lim[0] = qs->et1[qs->in1] - REL_COMP_RES;
  sf_high_lim[1] = qs->et2[qs->in2] - REL_COMP_RES;

  sh6gettop (pintpt, -1, lleft, lright, lleft + 1, lright + 1, &kstat);
  if (kstat < 0)
    goto error;
  /* Check endpoint of curve. */

  if (DEQUAL (sptpar[kpar1] + tref, st[kk - 1] + tref))
    *ll1 = SI_AT;
  if (DEQUAL (sptpar[kpar1] + tref, st[kn] + tref))
    *lr1 = SI_AT;

  /* Update pre-topology of intersection point.  */
  sh6settop (pintpt, -1, lleft[0], lright[0], lleft[1], lright[1], &kstat);
  if (kstat < 0)
    goto error;

  *jstat = 0;
  goto out;

  /* Error lower level routine.  */
error:*jstat = kstat;
  goto out;

out:
  return;
}
