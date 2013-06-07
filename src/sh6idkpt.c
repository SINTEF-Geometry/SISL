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
 * $Id: sh6idkpt.c,v 1.2 2001-03-19 15:59:08 afr Exp $
 *
 */


#define S6IDKPT


#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
sh6idkpt (SISLIntdat ** pintdat, SISLIntpt ** pintpt, int join, int *jstat)
#else
void
sh6idkpt (pintdat, pintpt, join, jstat)
     SISLIntdat **pintdat;
     SISLIntpt **pintpt;
     int join;
     int *jstat;
#endif


/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To remove an intersection point pintpt from pintdat.
*              pintpt is removed from all lists which it lies in if any.
*              If pintpt has exactly two neighbours, they are joined
*              together if the option join is selected.
*              After disconnection is done, pintpt is killed. If pintdat
*              is empty pintdat is killed and set to SISL_NULL.
*
*
*
* INPUT/OUTPUT:pintpt   - Pointer to a pointer to new intersection point.
*              pintdat  - Pointer to a pointer to intersection data.
*              join     - Flag for whether the lists are repaired.
*			   --ALA-- and kill all help-points connected
*			  to this point if this point is a main point.
*
*
* OUTPUT  :    jstat    - status messages
*                               = 2      : Pintpt is not in pintdat.
*                               = 1      : Pintpt is SISL_NULL
*                               = 0      : OK!
*                               < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : sh6err      - Gives error message.
*              freeIntpt  - free instant of intpt structure.
*
* WRITTEN BY : Ulf J. Krystad, 06.91.
*
*********************************************************************
*/
{
  int ki;			/* Counters.    */
  int knum;
  int kstat = 0;
  SISLIntpt *pnhbr_1 = SISL_NULL;	/* First neighbour  */
  SISLIntpt *pnhbr_2 = SISL_NULL;	/* Second neighbour */
  SISLIntpt *help_pt = SISL_NULL;	/* help point */
  int crv_dir_1 = 0;
  int crv_dir_2 = 0;
  int index1 = 0;
  int index2 = 0;
  int dummy;
  /* ------------------------------------------------*/
  
  *jstat = 0;
  
  if ((*pintpt) == SISL_NULL)
  {
     *jstat = 1;
     goto out;
  }
  
  if (join)
  {
     /* ALA-- We first remove all help point if this point is a main point. */
     if (sh6ismain(*pintpt))
	for (ki = 0; ki < (*pintpt)->no_of_curves; ki++)
	{
	   if (sh6ishelp(help_pt = sh6getnext(*pintpt, ki)))
	   {
	      sh6idkpt (pintdat, &help_pt, 1, &kstat);
	      if (kstat < 0)
		 goto error;
	   }
	}
     
     /* Remember the two neighbours */
     sh6getnhbrs (*pintpt, &pnhbr_1, &pnhbr_2, &kstat);
     if (kstat < 0)
	goto error;
     
     
     if (pnhbr_1 && pnhbr_2)
     {
	/* Two neighbours, remember crv_dir */
	sh6getlist (*pintpt, pnhbr_1, &dummy, &index1, &kstat);
	if (kstat < 0)
	   goto error;		/* Error. */
	if (kstat == 1)
	   goto err1;		/* pt1 and pt2 are not linked. */
	
	sh6getlist (*pintpt, pnhbr_2, &dummy, &index2, &kstat);
	if (kstat < 0)
	   goto error;		/* Error. */
	if (kstat == 1)
	   goto err1;		/* pt1 and pt2 are not linked. */
	
	crv_dir_1 = pnhbr_1->curve_dir[index1];
	crv_dir_2 = pnhbr_2->curve_dir[index2];
     }
  }

  
  for (; (*pintpt)->no_of_curves;)
  {
     /* Disconnect all */
     sh6disconnect (*pintpt, (*pintpt)->pnext[0], &kstat);
     if (kstat < 0)
	goto error;
  }
  
  /* Connect the two neighbours */
  if (pnhbr_1 && pnhbr_2)
  {
     sh6connect (pnhbr_1, pnhbr_2, &kstat);
     if (kstat < 0)
	goto error;
     
     /* UJK, MESZ 930617: Don't bother with curve_dir when 
	the points already were connected. */
     if (kstat != 1)
     {
	sh6getlist (pnhbr_1, pnhbr_2, &index1, &index2, &kstat);
	if (kstat < 0)
	   goto error;		/* Error. */
	if (kstat == 1)
	   goto err1;		/* pt1 and pt2 are not linked. */
	
	pnhbr_1->curve_dir[index1] = crv_dir_1;
	pnhbr_2->curve_dir[index2] = crv_dir_2;
     }
  }
  
  if ((*pintdat) == SISL_NULL)
  {
     freeIntpt (*pintpt);
     (*pintpt) = SISL_NULL;
     
     *jstat = 1;
     goto out;
  }
  
  
  /* Find pintpt in pintdat. */
  
  for (knum = -1, ki = 0; ki < (*pintdat)->ipoint; ki++)
  {
     if ((*pintdat)->vpoint[ki] == (*pintpt))
     {
	knum = ki;
	break;
     }
  }
  
  
  if (knum == -1)
     *jstat = 1;
  else
  {
     (*pintdat)->vpoint[knum] = (*pintdat)->vpoint[(*pintdat)->ipoint - 1];
     ((*pintdat)->ipoint)--;
     (*pintdat)->vpoint[(*pintdat)->ipoint] = SISL_NULL;
     
     
     
     if ((*pintdat)->ipoint == 0)
     {
	freeIntdat (*pintdat);
	(*pintdat) = SISL_NULL;
     }
  }
  
  freeIntpt (*pintpt);
  (*pintpt) = SISL_NULL;
  goto out;
  
  
err1:
  *jstat = -1;
  goto out;
  
error:
  *jstat = kstat;
  goto out;

out:;
}
