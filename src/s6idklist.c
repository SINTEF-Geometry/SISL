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
 * $Id: s6idklist.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDKLIST

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idklist(SISLIntdat **pintdat,SISLIntlist *pintlist,int *jstat)
#else
void s6idklist(pintdat,pintlist,jstat)
     SISLIntdat  **pintdat;
     SISLIntlist *pintlist;
     int     *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To remove an intersection list including all intersection points
*              in the list. The mother pintdat is updated.
*              If pintdat is empty, pintdat is killed and set to SISL_NULL.
*
*
*
* INPUT/OUTPUT:pintlist - Pointer to a list.
*              pintdat  - Pointer to a pointer to intersection data.
*
*
* OUTPUT     :jstat    - status messages  
*                               = 1      : Pintlist is not in pintdat.
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
* CALLS      : s6err      - Gives error message.
*              s6idkpt    - Kills a point.
*              freeIntpt  - free instant of intpt structure.
*
* WRITTEN BY : Ulf J. Krystad, SI, 08.89.
*
*********************************************************************
*/                                     
{
  SISLIntpt *qkillpt,*qnext,*qdum1,*qdum2;
  
  int ki,knum,kstat;  
  
  *jstat = 0;
  
  /* We have to be sure that we have an intdat structure. */
  
  if ((*pintdat) == SISL_NULL)
    goto out;
  
  if (pintlist == SISL_NULL)
    {
      *jstat = 1;
      goto out;
    }
  
  /* Now we have to find the index in the vlist array in pintdat. */
  
  
  for (ki=0,knum = -1; ki < (*pintdat)->ilist; ki++)
    if ((*pintdat)->vlist[ki] == pintlist)
      {
	knum = ki;
	break;
      }
  
  if (knum == -1)
    /* Not in the pintdat list. */
    *jstat = 1;
  else
    {
      pintlist->plast->pcurve = SISL_NULL;
      
      /* Kill all points in the list. */
      for (ki=0,qkillpt=pintlist->pfirst,qnext=qkillpt->pcurve;
	   qnext!=SISL_NULL;
	   qkillpt=qnext,qnext=qnext->pcurve)
	{
	  s6idkpt(pintdat,&qkillpt,&qdum1,&qdum2,&kstat);
	  if (kstat < 0) goto error;
	}
      s6idkpt(pintdat,&qkillpt,&qdum1,&qdum2,&kstat);
      if (kstat < 0) goto error;
      
      /* Update pintdat. */
      if ((*pintdat) != SISL_NULL)
	{
	  (*pintdat)->vlist[knum] = (*pintdat)->vlist[(*pintdat)->ilist-1];
	  ((*pintdat)->ilist)--;
	  (*pintdat)->vlist[(*pintdat)->ilist] = SISL_NULL;
	}
      freeIntlist(pintlist);
    }
  
  goto out;  
  
  error : *jstat = kstat;
  s6err("s6idklist",*jstat,0);
  goto out;                       
  
  out: ;
}
