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
 * $Id: s1905.c,v 1.3 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1905

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1905 (double econd1[], int ntype1[], int inpt1, int ik, int idim,
       int iopen, double **gcond2, int **mtype2, int *jnpt2, int *jstat)
#else
void
s1905 (econd1, ntype1, inpt1, ik, idim, iopen, gcond2, mtype2, jnpt2, jstat)
     double econd1[];
     int ntype1[];
     int inpt1;
     int ik;
     int idim;
     int iopen;
     double **gcond2;
     int **mtype2;
     int *jnpt2;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : Check legality of interpolation conditions, and adjust
*              to legal conditions if possible.
*
* INPUT      : econd1 - Array of interpolation conditions. Dimension
*                       is inpt1*idim.
*              ntype1 - Array containing kind of condition. Dimension
*                       is inpt1.
*                       =  0 : A point is given.
*                       =  d : The d'th derivatative condition to the
*                              previous point is given.
*                       = -d : The d'th derivatative condition to the
*                              next point is given.
*              inpt1  - Number of original interpolation conditions.
*              ik     - Order of interpolating curve.
*              idim   - Dimension of geometry space.
*              iopen - Indicates if the curve is to be open, closed or
*                       periodic.
*
* OUTPUT     : gcond2 - Adjusted interpolation conditions.
*              mtype2 - Type of adusted conditions. See description of
*                       ntype1.
*              jnpt2  - Number of adusted interpolation conditions.
*                       The dimension is inpt1*idim.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 91-04.
*
*********************************************************************
*/
{
   int kstat = 0;
  int kpos = 0;
  int kmaxpt = inpt1 + ik*(iopen != SISL_CRV_OPEN); /* Maximum number of 
						    adjusted conditions.  */
  int knpt = 0;			/* Current number of adjusted conditions. */
  int ki, kj, kl;		/* Counters.                              */
  int ktype;			/* Kind of interpolation condition.       */
  int kneg;			/* Indicates negative type indicator.     */
  int kder;			/* Order of differentiation.              */
  int *lder = SISL_NULL;		/* Kind of derivative.                    */
  double *sdum = SISL_NULL;		/* Help array.                            */
  double tdist;                 /* Distance between first and last point. */
  double tref;                  /* Referance value.                       */

  *jstat = 0;


  /* Allocate scratch for output arrays. Make sure that the arrays
     are large enough.  */

  if ((*gcond2 = newarray (kmaxpt * idim, DOUBLE)) == SISL_NULL)
    goto err101;
  if ((*mtype2 = newarray (kmaxpt, INT)) == SISL_NULL)
    goto err101;

  /* Allocate scratch for local arrays.  */

  if ((lder = new0array (ik, INT)) == SISL_NULL)
    goto err101;
  if ((sdum = newarray (idim, DOUBLE)) == SISL_NULL)
    goto err101;


  /* Find first positional condition.  */

  for (ki = 0; ki < inpt1; ki++)
    if (ntype1[ki] == 0)
      break;

  lder[0] = 1;
  (*mtype2)[0] = 0;
  memcopy (*gcond2, econd1 + ki * idim, idim, DOUBLE);
  knpt++;
  ki++;

  /* Move any derivative conditions to the first point after the position. */

  for (kj = ki - 2; kj >= 0 && ntype1[kj] < 0; kj--)
    {
      ktype = abs (ntype1[kj]);
      if (ktype >= ik)
	continue;		/* Not a legal derivative condition. */
      if (lder[ktype])
	continue;		/* Derivative condition already given. */
      lder[ktype] = 1;
      (*mtype2)[knpt] = ktype;
      memcopy ((*gcond2) + knpt * idim, econd1 + kj * idim, idim, DOUBLE);
      knpt++;
    }

  /* Copy the remaining derivative conditions of the first point. */

  for (; ki < inpt1 && ntype1[ki] > 0; ki++)
    {
      ktype = ntype1[ki];
      if (ktype >= ik)
	continue;		/* Not a legal derivative condition. */
      if (lder[ktype])
	continue;		/* Derivative condition already given. */
      lder[ktype] = 1;
      (*mtype2)[knpt] = ktype;
      memcopy ((*gcond2) + knpt * idim, econd1 + ki * idim, idim, DOUBLE);
      knpt++;
    }

  /* Traverse the remaining interpolation conditions and copy legal
     conditions. */

  for (; ki < inpt1; ki = kj)
    {
      /* Initiate array indicating occupied derivatives to zero. */

      for (kj = 0; kj < ik; kj++)
	lder[kj] = 0;

      /* Copy all conditions corresponding to current position. */

      kneg = 1;
      for (kj = ki; kj < inpt1 && (kneg || ntype1[kj] > 0); kj++)
	{
	  ktype = abs (ntype1[kj]);
	  if (ktype == 0)
	    kneg = 0;		/* Position condition reached. */
	  if (ktype >= ik)
	    continue;		/* Not a legal derivative condition. */
	  if (lder[ktype])
	    continue;		/* Derivative condition already given. */
	  lder[ktype] = 1;
	  (*mtype2)[knpt] = ntype1[kj];
	  memcopy ((*gcond2) + knpt * idim, econd1 + kj * idim, idim, DOUBLE);
	  knpt++;
	}
    }

  if (iopen != SISL_CRV_OPEN)
  {
     /* Closed curve requested. Let the first positional interpolation
	condition also be the last condition. First fetch derivative
	conditions.     */
     
     /* Test first if the first and last interpolation points is equal
	already.   */
     
     for (kj=ki-1; kj<=0; kj--)
        if (ntype1[kj] == 0) break;
     for (kl=0;  kl<inpt1; kl++)
        if (ntype1[kl] == 0) break;
     tdist = s6dist(econd1+kl*idim,econd1+kj*idim,idim);
     tref = MAX(s6length(econd1+kl*idim,idim,&kstat),
		s6length(econd1+kj*idim,idim,&kstat));
     
     if (DNEQUAL(tdist+tref,tref))
     {
	/* Initiate array indicating occupied derivatives to zero. */
	
	for (kj = 0; kj < ik; kj++)
	   lder[kj] = 0;
	
	/* Fetch derivative conditions to next point. */
	
	for (kj = ki; kj < inpt1 && ntype1[kj] < 0; kj++)
	{
	   ktype = abs (ntype1[kj]);
	   if (ktype >= ik)
	      continue;		/* Not a legal derivative condition. */
	   if (lder[ktype])
	      continue;		/* Derivative condition already given. */
	   lder[ktype] = 1;
	   (*mtype2)[knpt] = ntype1[kj];
	   memcopy ((*gcond2) + knpt * idim, econd1 + kj * idim, idim, DOUBLE);
	   knpt++;
	}
	
	/* Fetch derivative conditions prior to first point. */
	
	for (kj = 0; kj < inpt1 && ntype1[kj] > 0; kj++)
	{
	   ktype = ntype1[kj];
	   if (ktype >= ik)
	      continue;		/* Not a legal derivative condition. */
	   if (lder[ktype])
	      continue;		/* Derivative condition already given. */
	   lder[ktype] = 1;
	   (*mtype2)[knpt] = -ktype;
	   memcopy ((*gcond2) + knpt * idim, econd1 + kj * idim, idim, DOUBLE);
	   knpt++;
	}
	
	/* Find first interpolation point and copy it. */
	
	for (; ntype1[kj] != 0; kj++) ;
	
	(*mtype2)[knpt] = 0;
	memcopy ((*gcond2) + knpt * idim, econd1 + kj * idim, idim, DOUBLE);
	knpt++;
     }
  }

  /* Make sure that the last interpolation conditions are
     of decreasing order of interpolation. */

  kneg = 1;
  for (kder = 0, ki = knpt - 1; ki >= 0 && (kneg || (*mtype2)[ki] < 0); ki--, kder--)
    {
      if ((*mtype2)[ki] != kder)
	{
	  for (kj = ki - 1; kj >= 0 && (kneg || (*mtype2)[kj] < 0) &&
	       (*mtype2)[kj] != kder; kj--) ;
	  if ((*mtype2)[kj] == kder)
	    {
	      /* Interchange interpolation conditions. */

	      memcopy (sdum, (*gcond2) + kj * idim, idim, DOUBLE);
	      memcopy ((*gcond2) + kj * idim, (*gcond2) + ki * idim, idim, DOUBLE);
	      memcopy ((*gcond2) + ki * idim, sdum, idim, DOUBLE);
	      (*mtype2)[kj] = (kder == 0) ? -(*mtype2)[ki] : (*mtype2)[ki];
	      (*mtype2)[ki] = kder;
	      if (kder == 0)
		{
		  for (kj++; kj<ki; kj++)
		    (*mtype2)[kj] *= -1;
		}
	    }
	}
      kneg = 0;
    }

  /* Adjust size of output arrays.  */

  *gcond2 = increasearray (*gcond2, knpt * idim, DOUBLE);
  if (*gcond2 == SISL_NULL)
    goto err101;
  *mtype2 = increasearray (*mtype2, knpt, INT);
  if (*mtype2 == SISL_NULL)
    goto err101;

  /* Conditions adjusted.  */

  *jnpt2 = knpt;
  goto out;


  /* Error in scratch allocation.  */

err101:
  *jstat = -101;
  s6err ("s1905", *jstat, kpos);
  goto out;

out:
  /* Free scratch occupied by local array. */

  if (lder != SISL_NULL)
    freearray (lder);
  if (sdum != SISL_NULL)
    freearray (sdum);

  return;
}
