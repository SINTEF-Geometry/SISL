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
 * $Id: sh6trmlist.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH6TRIMLIST

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6trimlist (SISLIntpt * pt, SISLIntpt *** ptlist, int *no_of_points,
	     int *no_alloc)
#else
void
sh6trimlist (pt, ptlist, no_of_points, no_alloc)
     SISLIntpt *pt;
     SISLIntpt ***ptlist;
     int *no_of_points;
     int *no_alloc;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To find the maximal trim boundary (of coincidence)
*              containg the given point pt.
*
*
* INPUT      : pt            - Pointer to int point to be examined
*
*
* INPUT/OUTP:  ptlist        - Pointer to an array containing
*                              pointers to all intersection points
*                              that are contigous trim neighbours.
*               no_of_points - Number of points in ptlist array
*               no_alloc     - Allocation size of ptlist array
* OUTPUT     :
*              jstat     - Error flag.
*                         jstat =  0  => OK.
*                         jstat = -1  => Data structure inconsistent.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway. October 91.
*
*********************************************************************
*/
{
  int clean_up = FALSE;		/* Clean up on top level */
  int incr = 20;		/* Allocation size       */
  int ki;			/* Loop control          */
  /* --------------------------------------------------- */


  /* Check if point is a TRIM point */
  if (pt->iinter != SI_TRIM)
    goto out;

  /* Check if point is treated */
  if (pt->marker == -90)
    goto out;

  /* Mark point as treated */
  pt->marker = -90;


  if (*no_alloc <= *no_of_points)
    {
      if (*no_alloc == 0)
	{
	  clean_up = TRUE;
	  (*no_alloc) += incr;
	  *ptlist = newarray (*no_alloc, SISLIntpt *);
	  if (*ptlist == SISL_NULL)
	    goto out;
	}
      else
	{
	  clean_up = FALSE;
	  (*no_alloc) += incr;
	  *ptlist = increasearray (*ptlist, *no_alloc, SISLIntpt *);
	  if (*ptlist == SISL_NULL)
	    goto out;
	}
    }

  /* Fill in */
  (*ptlist)[*no_of_points] = pt;
  (*no_of_points)++;

  /* Treat all neighbours */
  for (ki = 0; ki < pt->no_of_curves; ki++)
    sh6trimlist (pt->pnext[ki], ptlist, no_of_points, no_alloc);


/* Must unmark the points in array if no_alloc == 0 */
  if (clean_up)
    for (ki = 0; ki < (*no_of_points); ki++)
      (*ptlist)[ki]->marker = 0;

  goto out;


out:
  return;
}
