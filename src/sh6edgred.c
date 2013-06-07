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
 * $Id: sh6edgred.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6EDGRED

/* commented out - guen Fri Oct 25 12:59:43 MEZ 1991
 include <stdio.h>
   commented out - guen Fri Oct 25 12:59:43 MEZ 1991 */

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
sh6edgred (SISLObject * po1, SISLObject * po2,
	   SISLIntdat * pintdat, int *jstat)
#else
void
sh6edgred (po1, po2, pintdat, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat *pintdat;
     int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Reduse main int point to help point if they are not legal.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pintdat  - Intersection data structure.
*
*
* OUTPUT     : jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      :
*
* REFERENCES :
*
* WRITTEN BY : Arne Laksaa, SI, Oslo, sep. -91.
* CORRECTED BY: UJK
*********************************************************************
*/
{
  int kstat, gstat, i, ki;
  int change = FALSE;
  int change_2 = FALSE;
  int num = 0;
  SISLIntpt *pt1 = SISL_NULL;
  SISLIntpt *pt2 = SISL_NULL;
  SISLIntpt *pcurr = SISL_NULL;

  if (pintdat != SISL_NULL)
    {
      do
	{
	  change_2 = FALSE;
	  /* If trim point is internal and one neighbours, change to help
	     point, if two neighbours unite till one of them */
	  do
	    {
	      change = FALSE;
	      for (i = 0; i < pintdat->ipoint; i++)
		{
		  pcurr = pintdat->vpoint[i];
		  if (pcurr->iinter == SI_TRIM)
		    {
		      sh6isinside (po1, po2, pcurr, &kstat);
		      if (kstat < 0)
			goto error;
		      if (kstat == 1)
			{
			  num = sh6nmbmain (pcurr, &kstat);
			  if (kstat < 0)
			    goto error;
			  if (num == 1)
			    {
			      sh6tohelp (pcurr, &kstat);
			      change = TRUE;
			    }
			  else if (num == 2)
			    {
			      sh6getnhbrs (pcurr, &pt1, &pt2, &gstat);
			      if (kstat < 0)
				goto error;
			      if (pt1->iinter == SI_TRIM &&
				  pt2->iinter == SI_TRIM)
				{
				  sh6idunite (&pintdat, &pt1, &pcurr,
					      DZERO, &kstat);
				  if (kstat < 0)
				    goto error;
				  change = TRUE;
				}
			    }
			}
		    }
		}
	  } while (change);


	  /* For a trim point on the edge with only one trim
             neighbour on an edge, unit till the other edge
             neighbour and change status of neighbour*/
	  do
	    {
	      change = FALSE;
	      for (i = 0; i < pintdat->ipoint; i++)
		{
		  pt1 = pt2 = SISL_NULL;
		  pcurr = pintdat->vpoint[i];
		  if (pcurr->iinter == SI_TRIM)
		    {
		      sh6isinside (po1, po2, pcurr, &kstat);
		      if (kstat < 0)
			goto error;
		      if (kstat == 2)
			{
			  for (ki = 0; ki < pcurr->no_of_curves; ki++)
			    {
			      pt1 = pcurr->pnext[ki];
			      if (pt1->iinter == SI_TRIM)
				{
				  sh6comedg (po1, po2, pcurr, pt1, &kstat);
				  if (kstat < 0)
				    goto error;
				    if (kstat)
				    {
				       if (pt2)
					  {
					     pt2 = SISL_NULL;
					     break;
					  }
					  else
					  pt2 = pt1;
				    }
				}
			    }
			  if (pt2)
			    {
			       /* sh6idunite (&pintdat, &pt2, &pcurr,
				              DZERO, &kstat);  */
			       /* UJK, 12.08.93  */
			      /* sh6idkpt (&pintdat, &pcurr, 1, &kstat);
			     sh6disconnect(pcurr,pt2,&kstat); */
			     pcurr->iinter = SI_SING;

			      /*------------------- */
			      /* If no trim neighbours on common
			         edge, remove trim status. */
			      pcurr = pt2;
			      kstat = 0;

			      for (ki = 0; ki < pcurr->no_of_curves; ki++)
				{
				  pt1 = pcurr->pnext[ki];
				  if (pt1->iinter == SI_TRIM)
				    {
				      sh6comedg (po1, po2, pcurr, pt1, &kstat);
				      if (kstat < 0)
					goto error;
				      if (kstat)
					break;

				    }

				}
			      /* -------------------- */
			      if (!kstat)
				pcurr->iinter = SI_SING;
			      change = TRUE;
			      change_2 = TRUE;
			    }

			}

		    }
		}
	  } while (change);
      } while (change_2);


      /* Reduce internal stuff */
      sh6red (po1, po2, pintdat, &kstat);
      if (kstat < 0)
	goto error;


      /* General edge treatment */

      /* UJK, aug 93, spesial branch for crv/crv */
      if (po1->iobj == SISLCURVE &&
	  po2->iobj == SISLCURVE )
      {
	 do
	 {
	    change = 0;
	    for (i = 0; i < pintdat->ipoint; i++)
	    {
	       if (sh6ismain (pintdat->vpoint[i]))
	       {
		  sh6getnhbrs (pintdat->vpoint[i], &pt1, &pt2, &gstat);
		  if (gstat == 1)
		  {
		     double parval;
		     SISLCurve *pcu=SISL_NULL;
		     if (pintdat->vpoint[i]->epar[0] == pt1->epar[0])
		     {
			parval = pintdat->vpoint[i]->epar[1];
			pcu    = po2->c1;
		     }
		     else if (pintdat->vpoint[i]->epar[1] == pt1->epar[1])
		     {
			parval = pintdat->vpoint[i]->epar[0];
			pcu    = po1->c1;
		     }
		     
		     if (pcu &&
			 parval > pcu->et[pcu->ik-1] &&
			 parval < pcu->et[pcu->in] )
			
			
		     {
			sh6tohelp (pintdat->vpoint[i], &kstat);
			if (kstat < 0)
			   goto error;
			change = 1;
		     }
		  }
	       }
	    }
	 } while (change);
      }
      else 
      { 
	 do
	 {
	    change = 0;
	    for (i = 0; i < pintdat->ipoint; i++)
	    {
	       if (sh6ismain (pintdat->vpoint[i]))
	       {
		  sh6isinside (po1, po2, pintdat->vpoint[i], &kstat);
		  if (kstat < 0)
		     goto error;
		  
		  /* ALA and VSK. Test if the point lies on edge in 
		     one or two objects.         */
		  if (kstat == 2 || kstat == 5)
		  {
		     sh6getnhbrs (pintdat->vpoint[i], &pt1, &pt2, &gstat);
		     if (gstat == 1)
		     {
			sh6comedg (po1, po2, pintdat->vpoint[i], pt1, &gstat);
			
			/* ALA and VSK. Test if the points lie on the same
			   edge in both objects if it lies on an edge in
			   both objects.                  */
			if ((kstat == 2 && gstat > 0) ||
			    (kstat == 5 && gstat == 3))
			{
			   sh6tohelp (pintdat->vpoint[i], &kstat);
			   if (kstat < 0)
			      goto error;
			   change = 1;
			}
		     }
		  }
	       }
	    }
	 } while (change);
      }
    }



  *jstat = 0;
  goto out;

  /* Error lower level routine.  */

error:(*jstat) = kstat;
  s6err ("sh6edgred", *jstat, 0);
  goto out;

out:
  return;
}

