/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl-copyright.h"

/*
 *
 * $Id: shmkhlppts.c,v 1.4 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SHMKHLPPTS

#include "sislP.h"
#if defined(SISLNEEDPROTOTYPES)
void
shmkhlppts (SISLObject * po1, SISLObject * po2, double aepsge,
	SISLIntdat ** rintdat, SISLEdge * vedge[], int *jnewpt, int *jstat)
#else
void
shmkhlppts (po1, po2, aepsge, rintdat, vedge, jnewpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **rintdat;
     SISLEdge *vedge[];
     int *jnewpt;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Make help points (and pretopology) for main points.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              aepsge   - Geometry resolution.
*              vedge    - SISLEdge intersection objects to the two
*                         objects in intersection problem.
*
*
* INPUT/OUTPUT : rintdat  - Intersection data of object-object intersection.
*
* OUTPUT       : jnewpt   - Number of new intersection points.
*                jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : sh1781          - Find pre-topology of 1D curve-point.
*              sh1780          - Find pre-topology of curve-curve.
*              sh1779          - Find pre-topology of 3D curve-surface.
*              sh1786          - Find pre-topology of 2D point-curve.
*              sh1787          - Find pre-topology of 2D point-surface.
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, 09.91
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 09-94. Fixed over-running
*              of 'up' array.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int knum = 0;			/* Number of intpt on edges.               */
  int ki;			/* Counter.                                */
  int kdim;			/* Dimension of geometry space.            */
  int knewpt = 0;		/* Number of new intersection points.      */
  int kobj;			/* Number of obj, used in s6idint          */
  int index1, index2;		/* Dummy in this context                   */
  SISLIntpt **up = SISL_NULL;	/* Array of poiners to intersection point. */
  /*  SISLIntpt *lup[3];*/		/* Array of poiners to intersection point. */
  SISLIntpt *qptint = SISL_NULL;	/* Pointer to internal intersection point. */
  SISLIntpt *qpt = SISL_NULL;	/* Pointer to intersection point.          */
  /* --------------------------------------------------------------------- */

  /* Init */
  *jstat = 0;
  *jnewpt = 0;

  /* Test if an intersection data structure exist.  */
  if (*rintdat == SISL_NULL)
    goto out;


  /* Fetch dimension of geometry space. */
  if (po1->iobj == SISLPOINT)

    kdim = po1->p1->idim;
  else if (po1->iobj == SISLCURVE)
    kdim = po1->c1->idim;
  else
    kdim = po1->s1->idim;

  /* Treat only cases:
     crv vs pt 1D
     crv vs crv
     crv vs sf
     crv vs pt 2D
     sf vs pt 2D
     */

  if (!(((po1->iobj == SISLCURVE && po2->iobj >= SISLCURVE) ||
	 (po2->iobj == SISLCURVE && po1->iobj >= SISLCURVE)) ||
	(kdim == 1 && (po1->iobj + po2->iobj) == (SISLPOINT + SISLCURVE)) ||
	(kdim == 2 && (po1->iobj + po2->iobj) >= (SISLPOINT + SISLCURVE))))
    goto out;

  /* Compute number of intersection points on edges, 0 1 or 2. */
  if (vedge[0] == SISL_NULL)
    knum = 0;
  else
    knum = vedge[0]->ipoint;

  if (vedge[1] != SISL_NULL)
    knum += vedge[1]->ipoint;


  if (knum > 0)
    {
      sh6edgpoint (vedge, &up, &knum, &kstat);
      if (kstat < 0)
	goto error;
    }

  if (knum == 2)
    {
      /* when two edge points, check if they are connected */
      sh6getlist (up[0], up[1], &index1, &index2, &kstat);
      if (kstat == 0)
	knum = 0;
    }

  if (knum == 0) /* BOH & ALA Added: 200993 */
  {
    /* Task performed.  */

    *jstat = 0;
    goto out;
  }

  /* Copy pointer of edge points into local pointer array */
  /*for (ki = 0; ki < knum; ki++)

    lup[ki] = up[ki]; */

  /* Get the internal point if any */
  if (po1->iobj == SISLPOINT)
    kobj = 2;
  else
    kobj = 1;

    s6idint (po1, po2, *rintdat, &qptint, kobj);
    if (qptint)
    {
       qpt = qptint;
       ki=-1;
    }
    else
    {
       ki = 0;
       qpt = up[0];
    }

  for (; ki < knum; ki++ )
    {

      if (ki >= 0) qpt = up[ki];

      /* Browse on the dimension of geometry space and the type of
         the input objects.     */

      if (kdim == 1 && ((po1->iobj == SISLCURVE && po2->iobj == SISLPOINT)
		     || (po2->iobj == SISLCURVE && po1->iobj == SISLPOINT)))
	{
	  /* Compute pre-topology in one-dimensional curve-level value
             intersection.            */

	  sh1781 (po1, po2, aepsge, rintdat, qpt, &knewpt, &kstat);
	  if (kstat < 0)
	    goto error;
	  *jnewpt += knewpt;
	}
      else if (po1->iobj == SISLCURVE && po2->iobj == SISLCURVE)
	{
	  /* curve-curve intersection.  */
	  sh1780 (po1, po2, aepsge, rintdat, qpt, &knewpt, &kstat);
	  if (kstat < 0)
	    goto error;
	  *jnewpt += knewpt;
	}
      else if (kdim == 2 &&
	       ((po1->iobj == SISLCURVE && po2->iobj == SISLPOINT)
		|| (po2->iobj == SISLCURVE && po1->iobj == SISLPOINT)))
	{
	  /* 2 dimensional point-curve intersection.  */

	  sh1786 (po1, po2, aepsge, rintdat, qpt, &knewpt, &kstat);
	  if (kstat < 0)
	    goto error;
	  *jnewpt += knewpt;
	}
      else if (kdim == 2 &&
	       ((po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT)
		|| (po2->iobj == SISLSURFACE && po1->iobj == SISLPOINT)))
	{
	  /* 2 dimensional point-surface intersection.  */

	  sh1787 (po1, po2, aepsge, rintdat, qpt, &knewpt, &kstat);
	  if (kstat < 0)
	    goto error;
	  *jnewpt += knewpt;
	}
      else if (kdim == 3 &&
	       ((po1->iobj == SISLCURVE && po2->iobj == SISLSURFACE) ||
		(po1->iobj == SISLSURFACE && po2->iobj == SISLCURVE)))
	{
	  /* Surface-curve intersection in 3-dimensional geometry space. */

	  sh1779 (po1, po2, aepsge, rintdat, qpt, &knewpt, &kstat);
	  if (kstat < 0)
	    goto error;
	  *jnewpt += knewpt;

	}
    }

  /* Task performed.  */

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
  if (up != SISL_NULL)
    freearray (up);

  return;
}
