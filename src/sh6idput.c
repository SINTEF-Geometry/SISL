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
 * $Id: sh6idput.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6IDPUT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6idput (SISLObject * po1, SISLObject * po2,
	  SISLIntdat ** rintdat, SISLIntdat * pintdat,
      int inr, double apar, SISLIntpt *** outintpt, int *npoint, int *jstat)
#else
void
sh6idput (po1, po2, rintdat, pintdat, inr, apar, outintpt, npoint, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat **rintdat;
     SISLIntdat *pintdat;
     int inr;
     double apar;
     SISLIntpt ***outintpt;
     int *npoint;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To insert all points in one intdat with one less number
*              of parameters into rintdat. New copies is made with
*              the missing parameter. Connections are kept.
*              The geometry data of the object not enhanced, is kept.
*
*
*
* INPUT      : po1      - First object in the intersection.
*              po2      - Second object in the intersection.
*              pintdat  - Pointer to intersection data with one less
*                         parameter than rintdat.
*              inr      - Number of the parameter that is missing in pintdat.
*              apar     - Parameter value of the missing parameter.
*
*
* OUTPUT     : rintdat  - Pointer to a pointer to intersection data.
*              outintpt - Help array containing pointers to the enhanced
*                         points in rintdat.
*              npoints -  Size of outinpt.
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
* CALLS      : s6err        - Gives error message.
*              sh6idnpt     - Insert a new intpt structure.
*              sh6idfcross  - Fetch cross intersections.
*              sh6idrmcross - Remove cross intersections.
*
* WRITTEN BY : Ulf J. Krystad, SI, 06.91.
* CHANGED BY : Vibeke Skytt, SI, 12.92.  Remove cross intersections in
*                                        closed and degenerate coincidence.
*********************************************************************
*/
{
  int kstat;			/* Local status variable.               */
  int kpos = 0;			/* Position of error.                   */
  int ki, kj;			/* Counters                             */
  int keep_first;		/* Flag, which object is not enhanced   */
  int kant;			/* Number of parameters in new points.  */
  int ind1, ind2;		/* Indexes (not used)                   */
  int no;			/* No. of doubles to copy into geo_aux  */
  double *scoef = SISL_NULL;		/* Pointer to array copying into geo_aux*/
  double *spar = SISL_NULL;		/* Storing uppdated parametervalues.    */
  SISLIntpt **uintpt = SISL_NULL;	/* Help array while getting connections */
  int iinter;
  double *nullp = SISL_NULL;
  /* VSK. Remove cross intersections. ----------------------------  */
  int kcross = 1;     /* Indicates existence of cross intersections. */
  int kncross = 0;    /* Number of cross intersections.              */
  int kpt;            /* Index in uintpt.                            */
  SISLIntpt *ucross[4];  /* Cross intersections.                     */

  *npoint = 0;

  /* Find out which object the parameter belongs to */
  if (inr < po1->iobj)
    keep_first = 0;
  else
    keep_first = 1;

  /* Do we have an intdat structure? */
  if (pintdat == SISL_NULL)
    {
      *jstat = 0;
      goto out;
    }

  /* Computing number of new parameter direction. */
  kant = pintdat->vpoint[0]->ipar + 1;


  if (inr < 0 || inr >= kant)
    goto err191;

  *npoint = pintdat->ipoint;

  /* Allocate an array for intersection points. */
  if ((uintpt = newarray (pintdat->ipoint, SISLIntpt *)) == SISL_NULL)
    goto err101;

  /* Allocate an array for parametervalues. */
  if ((spar = newarray (kant, double)) == SISL_NULL)
    goto err101;


  /* Enhance all intersection points. */
  for (ki = 0; ki < pintdat->ipoint; ki++)
    {
      /* Insert the missing parameter value. */

      for (kj = 0; kj < inr; kj++)
	spar[kj] = pintdat->vpoint[ki]->epar[kj];
      spar[kj] = apar;
      for (kj++; kj < kant; kj++)
	spar[kj] = pintdat->vpoint[ki]->epar[kj - 1];

      iinter = pintdat->vpoint[ki]->iinter;

      uintpt[ki] = hp_newIntpt (kant, spar, pintdat->vpoint[ki]->adist,
				iinter,
				pintdat->vpoint[ki]->left_obj_1[0],
				pintdat->vpoint[ki]->right_obj_1[0],
				pintdat->vpoint[ki]->left_obj_2[0],
				pintdat->vpoint[ki]->right_obj_2[0],
			     (keep_first ? pintdat->vpoint[ki]->size_1 : 0),
			     (keep_first ? 0 : pintdat->vpoint[ki]->size_2),
		     (keep_first ? pintdat->vpoint[ki]->geo_data_1 : nullp),
		    (keep_first ? nullp : pintdat->vpoint[ki]->geo_data_2));

      if (uintpt[ki] == SISL_NULL)
	goto err101;

      /* Store info from lower level object */
      if (keep_first)
	{
	  /* Store second object geometry in geo_aux */
	  no = pintdat->vpoint[ki]->size_2;
	  scoef = pintdat->vpoint[ki]->geo_data_2;
	}
      else
	{
	  /* Store first object geometry in geo_aux */
	  no = pintdat->vpoint[ki]->size_1;
	  scoef = pintdat->vpoint[ki]->geo_data_1;
	}
      /*      if (no > 0)
        	memcopy (uintpt[ki]->geo_aux, scoef, (no < 6) ? no : 6, DOUBLE); */



    }



  /* Insert all new intersection points in rintdat. */
  for (ki = 0; ki < pintdat->ipoint; ki++)
    {
      sh6idnpt (rintdat, &uintpt[ki], 1, &kstat);
      if (kstat < 0)
	goto error;
    }

  /* Transform the connections. */
  for (ki = 0; ki < pintdat->ipoint; ki++)
    {
      for (kj = ki + 1; kj < pintdat->ipoint; kj++)
	{
	  sh6getlist (pintdat->vpoint[ki], pintdat->vpoint[kj],
		      &ind1, &ind2, &kstat);
	  if (kstat < 0)
	    goto error;
	  if (kstat == 0)
	    {
	      sh6idcon (rintdat, &uintpt[ki], &uintpt[kj], &kstat);
	      if (kstat < 0)
		goto error;
	    }
	}

      if (sh6ismain (pintdat->vpoint[ki]) &&
	  sh6nmbmain (pintdat->vpoint[ki], &kstat))
	{
	  sh6tomain (uintpt[ki], &kstat);
	  if (kstat < 0)
	    goto error;
	}
    }
  
  if (po1->iobj > SISLPOINT && po2->iobj > SISLPOINT)
  {
     /* There is a possibility for cross intersections. Check the
	intersection data.  */
     
     kpt = 0;
     while (kpt < (*npoint))
     {
	kncross = 0;
	ucross[kncross] = uintpt[kpt];
	kncross = 1;
	
	/* Fetch cross intersections.  */
	
	sh6idfcross(*rintdat,ucross,&kncross,po1->iobj,po2->iobj,&kstat);
	kcross = kstat;
	
	if (kcross)
	{
	   /* Remove cross intersections. */
	   
	   sh6idrmcross(po1, po2, rintdat, ucross, kncross, uintpt, 
			*npoint, &kstat);
	   if (kstat < 0) goto error;
	   
	   if (kstat)
	   {
	      /* Points have been removed. Update uintpt.  */
	      
	      for (kj=0; kj<*npoint; kj++)
		 if (uintpt[kj] == SISL_NULL)
		 {
		    uintpt[kj] = uintpt[(*npoint)-1];
		    kj--;
		    (*npoint)--;
		 }
	   }
	   else kpt++;
	}
	else kpt++;
     }
  }
     
  *jstat = 0;
  goto out;


/* Error in inserted parameter number.  */

err191:*jstat = -191;
  s6err ("sh6idput", *jstat, kpos);
  goto out;


/* Error in space allocation.  */

err101:*jstat = -101;
  s6err ("sh6idput", *jstat, kpos);
  goto out;

/* Error in sub function.  */

error:*jstat = kstat;
  s6err ("sh6idput", *jstat, kpos);
  goto out;

out:*outintpt = uintpt;
  if (spar != SISL_NULL)
    freearray (spar);
}
