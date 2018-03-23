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

#define SH6GETSEGDIV

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
int sh6getsegdiv(SISLSurf *ps, int idiv, double epar[], int *seg_flag)
#else
  int sh6getsegdiv(ps, idiv, epar, seg_flag)
     SISLSurf *ps;
     int idiv;
     double epar[];
     int *seg_flag;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : 
*
*
*
* INPUT      : ps     - Surface to fetch segmentation information
*              idiv   - Parameter directions where a segmentation parameter
*                       is wanted
*                                                                     
*
* OUTPUT     : epar   - Segmentation parameters
*              return value - Indicates found segmentation parameters
*			     = 0      : No parameters found
*			     = 1      : Segmentation in 1. par. dir. found
*                            = 2      : Segmentation in 2. par. dir. found
*                            = 3      : Segmentation in both par. dir. found
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 2018-02
*
*********************************************************************
*/                                     
{
  int ki;
  int found = 0;
  double tmid;
  double tmindist;
  int prio;
  double curr_par;
  int curr_prio;
  double par;
  double tstart, tend;

  *seg_flag = 0;

  if ((idiv == 1 || idiv == 3) && ps->seg1)
    {
      tstart = ps->et1[ps->ik1-1];
      tend = ps->et1[ps->in1];
      tmid = 0.5*(tstart + tend);
      tmindist = HUGE;
      prio = 0;
      for (ki=0; ki<ps->seg1->num_seg; ++ki)
	{
	  curr_par = ps->seg1->seg_val[ki];
	  curr_prio = ps->seg1->seg_type[ki];
	  if (curr_prio > prio ||
	      (curr_prio == prio && fabs(curr_par-tmid) < tmindist))
	    {
	      par = curr_par;
	      prio = curr_prio;
	      tmindist = fabs(curr_par-tmid);
	    }
	}
      if (DNEQUAL(par, tstart) && par > tstart && DNEQUAL(par, tend) &&
	  par < tend)
	{
	  epar[0] = par;
	  found += 1;
	  if (prio == TANGENTIAL_BELT_LEFT || prio == TANGENTIAL_BELT_RIGHT)
	    (*seg_flag) += prio;
	}
    }

  if ((idiv == 2 || idiv == 3) && ps->seg2)
    {
      tstart = ps->et2[ps->ik2-1];
      tend = ps->et2[ps->in2];
      tmid = 0.5*(tstart + tend);
      tmindist = HUGE;
      prio = 0;
      for (ki=0; ki<ps->seg2->num_seg; ++ki)
	{
	  curr_par = ps->seg2->seg_val[ki];
	  curr_prio = ps->seg2->seg_type[ki];
	  if (curr_prio > prio ||
	      (curr_prio == prio && fabs(curr_par-tmid) < tmindist))
	    {
	      par = curr_par;
	      prio = curr_prio;
	      tmindist = fabs(curr_par-tmid);
	    }
	}
      if (DNEQUAL(par, tstart) && par > tstart && DNEQUAL(par, tend) &&
	  par < tend)
	{
	  epar[1] = par;
	  found += 2;
	  if (prio == TANGENTIAL_BELT_LEFT || prio == TANGENTIAL_BELT_RIGHT)
	    {
	      if ((*seg_flag) == 0)
		(*seg_flag) += prio;
	      else
		(*seg_flag) = 0;  /* Tangential belt segmentation is
				     expected in one parameter dir. only */
	    }
	}
    }

  return found;
}
