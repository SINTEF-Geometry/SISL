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
 * $Id: s6idcon.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDCON

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void
s6idcon_s9turn(SISLIntpt *);
static void
s6idcon_s9endturn(SISLIntdat *,SISLIntpt *);
#else
static void s6idcon_s9turn();
static void s6idcon_s9endturn();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
s6idcon(SISLIntdat **pintdat,SISLIntpt **pintpt1,SISLIntpt **pintpt2,int *jstat)
#else
void s6idcon(pintdat,pintpt1,pintpt2,jstat)
     SISLIntdat **pintdat;
     SISLIntpt  **pintpt1;
     SISLIntpt  **pintpt2;
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To connect two intersection points in pintdat into a list.
*              If pintdat is SISL_NULL a new pintdat is also made.
*              If  one of pintpt is close to an other intersection point
*              the object pintpt is pointing to is freed, and
*              pintpt is set to point to the already inserted point.
*
*
*
* INPUT      : pintpt1  - Pointer to a pointer to new intersection point.
*              pintpt2  - Pointer to a pointer to new intersection point.
*              pintdat  - Pointer to a pointer to intersection date.
*
*
* OUTPUT     : jstat  - status messages  
*                               = 3      : Only "one" junction point.
*                               = 2      : Only one point.
*                               = 1      : Already connected.
*                               = 0      : Connection done.
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
*              s6idnpt    - Insert a new intpt structure.
*              copyIntpt  - Copy an intpt structure.
*              newIntdat  - Create new intdat structure.
*
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  int kstat;                /* Local status variable.                     */
/*guen  int kpos;   */               /* Position of error.                         */
/*guen changed into:*/
  int kpos=0;                 /* Position of error.                         */

  int kfirst1,kfirst2;      /* To mark if the point is first in the list. */
  int ki1,ki2;              /* Counters                                   */
  SISLIntpt *qpt1,*qpt2;
  
  
  /* First we have to be sure that pintdat contain the two points. */
  
  s6idnpt(pintdat,pintpt1,1,&kstat);
  if (kstat < 0) goto error;
  
  s6idnpt(pintdat,pintpt2,1,&kstat);
  if (kstat < 0) goto error;
  
  
  qpt1 = *pintpt1;
  qpt2 = *pintpt2;
  
  
  /* Then we have to be sure that we do not have the same points as
     copies, junction points. */
  
  if (qpt1->iinter == 2 || qpt2->iinter == 2)
    {
      if (qpt1->iinter == 2 && qpt2->iinter == 2)
	{
	  for (ki1=0; ki1 < qpt1->ipar; ki1++)
	    if (qpt1->epar[ki1] != qpt2->epar[ki1]) break;
	  
	  if (ki1 == qpt1->ipar)
	    {
	      *jstat = 3;
	      goto out;
	    }
	}
      
      if (qpt1->iinter == 2)
	{
	  for (ki1=0; ki1 < (*pintdat)->ipoint; ki1++)
	    {
	      for (ki2=0; ki2 < qpt1->ipar; ki2++)
		if (qpt1->epar[ki2] != (*pintdat)->vpoint[ki1]->epar[ki2])
		  break;
	      
	      if (ki2 == qpt1->ipar)
		{
		  /* UJK && ALA 19.09.90 qpt1 changed to qpt2. */
		  
		  if (qpt2->pcurve == (*pintdat)->vpoint[ki1] || 
		      (*pintdat)->vpoint[ki1]->pcurve == qpt2)
		    {
		      /* The points are already connected. */
		      *jstat = 1;
		      goto out;
		    }
		}
	    }
	}
      
      if (qpt2->iinter == 2)
	{
	  for (ki1=0; ki1 < (*pintdat)->ipoint; ki1++)
	    {
	      for (ki2=0; ki2 < qpt2->ipar; ki2++)
		if (qpt2->epar[ki2] != (*pintdat)->vpoint[ki1]->epar[ki2])
		  break;
	      
	      if (ki2 == qpt2->ipar)
		{
		  /* UJK && ALA 19.09.90 qpt2 changed to qpt1. */
		  if (qpt1->pcurve == (*pintdat)->vpoint[ki1] || 
		      (*pintdat)->vpoint[ki1]->pcurve == qpt1)
		    {
		      /* The points are already connected. */
		      *jstat = 1;
		      goto out;
		    }
		}
	    }
	}
    }
  
  
  
  if (qpt1 == qpt2)
    /* There is only one point. */
    *jstat = 2;
  if (qpt1->pcurve == qpt2 || qpt2->pcurve == qpt1)
    /* The points are already connected. */
    *jstat = 1;
  else
    {
      /* We have to be sure that if one of the points is in the end of 
	 a list than this point is the first point. */
      
      if (qpt1->pcurve != SISL_NULL && qpt2->pcurve == SISL_NULL)
        {
	  SISLIntpt *pt;
	  
	  pt = qpt1;
	  qpt1 = qpt2;
	  qpt2 = pt;
        }
      
      /* Computing the index of the point pointing to the first point.    */
      
      for (ki1=0; ki1<(*pintdat)->ipoint; ki1++)
        if ((*pintdat)->vpoint[ki1]->pcurve == qpt1)
	  break;
      
      if ( ki1 < (*pintdat)->ipoint)
        kfirst1 = 0;
      else
        kfirst1 = 1;
      
      /* Computing the index of the point pointing to the sescond point.  */
      
      for (ki2=0; ki2<(*pintdat)->ipoint; ki2++)
        if ((*pintdat)->vpoint[ki2]->pcurve == qpt2)
	  break;
      
      if ( ki2 < (*pintdat)->ipoint)
        kfirst2 = 0;
      else
        kfirst2 = 1;
      
      /* If the first point is not at end, than we have to
	 reorganize the first list.  */
      
      if (qpt1->pcurve != SISL_NULL)
        {
	  if (kfirst1)
	    s6idcon_s9turn(qpt1);                  /* First point is at start. */
	  else                               /* First point is internal. */
	    {
	      /* We have a junction point. We therfor make a copy of
		 this point, and set this copy to the first point. */
	      
	      qpt1->iinter = 2;
	      
	      if((qpt1 = copyIntpt(qpt1)) == SISL_NULL) goto err101;
	      
	      s6idnpt(pintdat,&qpt1,0,&kstat);
	      if (kstat < 0) goto error;
	    }
        }
      
      
      if (kfirst2)                             /*Second point is at start.*/
        qpt1->pcurve = qpt2;
      else if (qpt2->pcurve == SISL_NULL)     /* Second point is at end. */
        {
	  s6idcon_s9endturn(*pintdat,qpt2);
	  qpt1->pcurve = qpt2;
        }
      else                          /* Second point is an internal point. */
        {
	  /* We have a junction point. We therfor make a copy of
	     this point, and set the first point  to point to this copy. */
	  
	  qpt2->iinter = 2;
	  
	  if((qpt2 = copyIntpt(qpt2)) == SISL_NULL) goto err101;
	  
	  s6idnpt(pintdat,&qpt2,0,&kstat);
	  if (kstat < 0) goto error;
	  
	  qpt1->pcurve = qpt2;
        }
      *jstat = 0;
    }
  
  goto out;
  

/* Error in space allocation.  */

err101: *jstat = -101;
        s6err("s6idcon",*jstat,kpos);
        goto out;

/* Error in sub function.  */

error:  *jstat = kstat;
        s6err("s6idcon",*jstat,kpos);
        goto out;

 out: ;
}

#if defined(SISLNEEDPROTOTYPES)
static void
s6idcon_s9turn(SISLIntpt *pt)
#else
static void s6idcon_s9turn(pt)
     SISLIntpt *pt;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To turn a list of intersection points where we have
*              a pointer at the start of the list.
*
*
* INPUT      : pt      - Pointer to the first intersection point in list.
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
*
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  register SISLIntpt *pt1,*pt2;/* Help pointer to traverse lists.*/
  
  pt1 = pt->pcurve;
  pt2 = pt1->pcurve;
  pt->pcurve = SISL_NULL;  
  pt1->pcurve = pt;
  
  while (pt2 != SISL_NULL)
    {
      pt  = pt1;
      pt1 = pt2;
      pt2 = pt2->pcurve;
      pt1->pcurve = pt;
    }
}

#if defined(SISLNEEDPROTOTYPES)
static void
s6idcon_s9endturn(SISLIntdat *pintdat,SISLIntpt *pt)
#else
static void s6idcon_s9endturn(pintdat,pt)
     SISLIntdat *pintdat;
     SISLIntpt  *pt;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To turn a list of intersection points where we
*              just have the last element in the list.
*
*
* INPUT      : pt      - Pointer to the last intersection point in list.
*
*
* OUTPUT     : pintdat - Intersection dates where the list to be
*                        turned is.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  register int ki;
  
  while(1)
    {
      for (ki=0; ki < pintdat->ipoint; ki++)
        if (pintdat->vpoint[ki]->pcurve == pt)
	  break;
      
      if (ki < pintdat->ipoint)
	pt = pintdat->vpoint[ki];
      else
	break;
    }
  
  s6idcon_s9turn(pt);
}
