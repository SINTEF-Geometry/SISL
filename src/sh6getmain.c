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
 * $Id: sh6getmain.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETMAIN

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
SISLIntpt *
sh6getmain (SISLIntpt * pt)
#else
SISLIntpt *
sh6getmain (pt)
     SISLIntpt *pt;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given a help point, find the unique main point it is
*              linked to. If there is no such main point, return
*              SISL_NULL.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
*
* OUTPUT     :
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
* CHANGED BY: Ulf J. Krystad, SI, Oslo, Norway. September 91.
*********************************************************************
*/
{
  int ki;			/* Loop control */
  int kstat;			/* Local status */
  int more = TRUE;		/* Loop control */
  SISLIntpt *mainpt = SISL_NULL;
  SISLIntpt *pt1 = SISL_NULL;
  SISLIntpt *pt2 = SISL_NULL;
  SISLIntpt *prev = SISL_NULL;
  SISLIntpt *pcurr = SISL_NULL;
  SISLIntpt *pnext = SISL_NULL;
  /* ------------------------------------------------------------- */


  if (!sh6ishelp (pt))
    goto out;

  for (ki = 0; ki < pt->no_of_curves; ki++)
    {
      if (sh6ismain (pt1 = sh6getnext (pt, ki)))
	{
	  mainpt = pt1;
	  break;
	}
    }

  if (!mainpt)
    {
      /* No close neighbour is main, check along list
         if not meeting point. */
      sh6getnhbrs (pt, &pt1, &pt2, &kstat);
      if (kstat == 1)
	{
	  /* Terminator, go towards other end */
	  prev = pt;
	  pcurr = pt1;
	  more = TRUE;

	  while ((!mainpt) && more)
	    {
	      sh6getother (pcurr, prev, &pnext, &kstat);
	      if (kstat < 0)
		goto error;

	      if (pnext && (pnext != pt))
		{
		  if (sh6ismain (pnext))
		    mainpt = pnext;
		  else
		    {
		      prev = pcurr;
		      pcurr = pnext;
		      pnext = SISL_NULL;
		    }
		}
	      else
		more = FALSE;

	    }
	}

      else if (kstat == 0)
	{
	  /* Two neighbours, search both directions */
	  for (ki = 0, prev = pt, pcurr = pt1, more = TRUE; (!mainpt) && (ki < 2);
	       ki++, prev = pt, pcurr = pt2, more = TRUE)

	    while ((!mainpt) && more)
	      {
		sh6getother (pcurr, prev, &pnext, &kstat);
		if (kstat < 0)
		  goto error;

		if (pnext && (pnext != pt))
		  {
		    if (sh6ismain (pnext))
		      mainpt = pnext;
		    else
		      {
			prev = pcurr;
			pcurr = pnext;
			pnext = SISL_NULL;
		      }
		  }
		else
		  more = FALSE;

	      }
	}
    }

  goto out;

  /* ------------------------------------------------------------- */
error:mainpt = SISL_NULL;
  s6err ("sh6getmain", kstat, 0);
  goto out;



out:
  return mainpt;
}

