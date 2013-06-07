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
 * $Id: s6idget.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDGET

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idget(SISLObject *po1,SISLObject *po2,int ipar,double apar,SISLIntdat *pintdat,
	     SISLIntdat **rintdat,int *jstat)
#else
void s6idget(po1,po2,ipar,apar,pintdat,rintdat,jstat)
     SISLObject *po1;
     SISLObject *po2;
     int    ipar;
     double apar;
     SISLIntdat *pintdat;
     SISLIntdat **rintdat;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To insert all points in one intdat with one less number
*              of parameters into rintdat. New copies is made with
*              the missing parameter. All listes is also to
*              be uppdated.
*
*
*
* INPUT      : po1     - First object in the intersection.
*              po2     - Second object in the intersection.
*              ipar    - Number of the parameter that is missing in rintdat.
*              apar    - Parameter value of the missing parameter.
*              pintdat - Pointer to intersection data on the mother problem.
*
*
* OUTPUT     : rintdat  - Pointer to a pointer to intersection data.
*              jstat    - status messages  
*                               = 0      : Inserting done.
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
* WRITTEN BY : Arne Laksaa, 06.89.
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              Error message corrected
*
*********************************************************************
*/                                     
{
  int kstat;                    /* Local status variable.                 */
  int kpos=0;                   /* Position of error.                     */
  int ki,kj,kn;
  double tstart[4];
  double tend[4];
  double spar[4];               /* Storing uppdated parametervalues.      */
  SISLIntpt *qpt;
  
  
  if (pintdat == SISL_NULL)
    goto out;
  
  if (po1->iobj == SISLCURVE)
    {
      tstart[0] = po1->c1->et[po1->c1->ik - 1];
      tend[0]   = po1->c1->et[po1->c1->in];
    }
  else if (po1->iobj == SISLSURFACE)
    {
      tstart[0] = po1->s1->et1[po1->s1->ik1 - 1];
      tend[0]   = po1->s1->et1[po1->s1->in1];
      tstart[1] = po1->s1->et2[po1->s1->ik2 - 1];
      tend[1]   = po1->s1->et2[po1->s1->in2];
    }
  
  if (po2->iobj == SISLCURVE)
    {
      tstart[po1->iobj]   = po2->c1->et[po2->c1->ik - 1];
      tend[po1->iobj]     = po2->c1->et[po2->c1->in];
    }
  else if (po2->iobj == SISLSURFACE)
    {
      tstart[po1->iobj]   = po2->s1->et1[po2->s1->ik1 - 1];
      tend[po1->iobj]     = po2->s1->et1[po2->s1->in1];
      tstart[po1->iobj+1] = po2->s1->et2[po2->s1->ik2 - 1];
      tend[po1->iobj+1]   = po2->s1->et2[po2->s1->in2];
    }
  
  tstart[ipar] = tend[ipar] = apar;
  
  
  /* Uppdate the array. */
  
  for (ki=0; ki<pintdat->ipoint; ki++)
    {
      for (kj=0; kj<pintdat->vpoint[ki]->ipar; kj++)
	if((DNEQUAL(pintdat->vpoint[ki]->epar[kj],tstart[kj]) &&
	    pintdat->vpoint[ki]->epar[kj] < tstart[kj]) ||
	   (DNEQUAL(pintdat->vpoint[ki]->epar[kj],tend[kj]) &&
	    pintdat->vpoint[ki]->epar[kj] > tend[kj]))
	  break;
      
      if (kj == pintdat->vpoint[ki]->ipar)
	{
	  for(kn=0; kn<ipar; kn++) spar[kn] = pintdat->vpoint[ki]->epar[kn];
	  for(; kn<pintdat->vpoint[ki]->ipar-1; kn++)
	    spar[kn] = pintdat->vpoint[ki]->epar[kn+1];
	  
	  /* Than we can insert all new intersection points in rintdat. */
	  
	  qpt = newIntpt(pintdat->vpoint[ki]->ipar-1,spar,DZERO);
	  if (qpt == SISL_NULL) goto err101;
	  
	  s6idnpt(rintdat,&qpt,1,&kstat);
	  if (kstat < 0) goto error;
	}
    }
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */

  err101: 
    *jstat = -101;
    s6err("s6idget",*jstat,kpos);
    goto out;

  /* Error in sub function.  */

  error: 
    *jstat = kstat;
    s6err("s6idget",*jstat,kpos);
    goto out;

  out: ;
}
