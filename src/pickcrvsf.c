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
 * $Id: pickcrvsf.c,v 1.2 2001-03-19 15:58:40 afr Exp $
 *
 */


#define PICK_CRV_SF

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   pick_crv_sf(SISLObject *po1, SISLObject *po2,int ipar,
	       SISLIntpt *pt1,SISLIntpt *pt2,SISLCurve **rcrv,
	       int *jstat)
#else
void pick_crv_sf(po1,po2,ipar,pt1,pt2,rcrv,jstat)
SISLObject *po1;
SISLObject *po2;
int ipar;
	       SISLIntpt *pt1;
	       SISLIntpt *pt2;
	       SISLCurve **rcrv;
	       int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Pick a curve along a given parameter direction
*              in a surface.
*
*
* INPUT      : po1  - first object in intersection.
*              po2  - second object in intersection.
*              pt1  - first intersection point.
*              pt2  - second intersection point.
*              ipar - index of constant parameter value
*
*
*
* OUTPUT     : rcrv - SISL curve.
*              jstat  - status messages  
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
* CALLS      : s1436 - Pick curve with constant second parameter.
*              s1437 - Pick curve with constant first parameter.
*
* WRITTEN BY : UJK, SI, 91-09.
*
*********************************************************************
*/                                     
{
  int kstat = 0;        /* Local status parameter.                        */
  int kpos = 0;         /* Position of error.                             */
  int index=0;          /* Index of other par dir in surf.                */
  int first_const;      /* Flag, const first direction or not             */
  double tpar;          /* Parameter value of curve in constant parameter
			   direction.                                     */
  SISLSurf *ps1=SISL_NULL;   /* Pointer to surf to pick crv from               */
  SISLCurve *pick_crv=SISL_NULL;/* Picked curve before trimming.               */
  /* -------------------------------------------------------------------- */
  if (ipar < 0 || ipar >= po1->iobj + po2->iobj) goto errinp;

  if (ipar >= po1->iobj)
  {
     /* pick from second object (must be a sf) */
     if (po2->iobj != SISLSURFACE) goto errinp;
     ps1 = po2->s1;
     index = (ipar == po1->iobj) ? po1->iobj + 1 : po1->iobj;
  }
  else
  {
     /* pick from first object (must be a sf) */
     if (po1->iobj != SISLSURFACE) goto errinp;
     ps1 = po1->s1;
     index = (ipar == 0) ? 1 : 0;
  }
  
  if (ipar < index) first_const = TRUE;
  else first_const = FALSE;
  tpar = pt1->epar[ipar];
  
  
  if (first_const == FALSE)
    {
       /* Pick curve with constant second parameter.  */
       s1436(ps1,tpar,&pick_crv,&kstat);
       if (kstat < 0) goto error;
    }
  else 
    {
       /* Pick curve with constant first parameter.  */
       s1437(ps1,tpar,&pick_crv,&kstat);
       if (kstat < 0) goto error;
    }
  
  /* SISLCurve picked, now trim it.  */
  if (DEQUAL(pt1->epar[index], pick_crv->et[pick_crv->ik-1]) &&
      DEQUAL(pt2->epar[index], pick_crv->et[pick_crv->in]))
    {
       /* Return the whole curve */
       (*rcrv)  = pick_crv;
       pick_crv = SISL_NULL;
    }
  
  
  else if(DEQUAL(pt1->epar[index], pick_crv->et[pick_crv->in]) &&
	  DEQUAL(pt2->epar[index], pick_crv->et[pick_crv->ik-1]))
    {
       /* Return the whole curve, but turn it first */
       /* Return the whole curve */
       (*rcrv)  = pick_crv;
       pick_crv = SISL_NULL;
       s1706(*rcrv); 
    } 
  else
  {
     /* Return a part of the curve */ 
     double amin = min(pt1->epar[index], pt2->epar[index]);
     double amax = max(pt1->epar[index], pt2->epar[index]);
     
     if (pick_crv->cuopen == SISL_CRV_PERIODIC)
	s1713(pick_crv,amin,amax,rcrv,&kstat);
     else
	s1712(pick_crv,amin,amax,rcrv,&kstat);
     
     if (kstat < 0) goto error;
     
     if (pt1->epar[index] > pt2->epar[index]) s1706(*rcrv); 

    }
  
  *jstat = 0;
  goto out;
  
  /* Error in input.  */
  errinp : *jstat = -1;
  s6err("pick_crv_sf",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */
  error : *jstat = kstat;
  s6err("pick_crv_sf",*jstat,kpos);
  goto out;
  
 out: if (pick_crv) freeCurve(pick_crv);
}
