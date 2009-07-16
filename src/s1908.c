/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s1908.c,v 1.2 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1908

#include "sislP.h"
#define MAX_SIZE  30

#if defined(SISLNEEDPROTOTYPES)
void
s1908 (double econd1[], int ntype1[], double epar[], int inpt1, int ik, int idim,
       int iopen, double **gcond2, int **mtype2, double *mpar[],
       int *jnpt2, int *jstat)
#else
void
s1908 (econd1, ntype1, epar, inpt1, ik, idim, iopen,
       gcond2, mtype2, mpar, jnpt2, jstat)
     double econd1[];
     int ntype1[];
     double epar[];
     int inpt1;
     int ik;
     int idim;
     int iopen;
     double *gcond2[];
     int *mtype2[];
     double *mpar[];
     int *jnpt2;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : Check legality of interpolation conditions, and adjust
*              to legal conditions if possible. The legal parameter values
*	       are also sent back.
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
*	       epar   - The parameter values. In the open case the lenght
*                       of this array is inpt1. In the closed or periodic
*                       case the length is inpt1+1. The last entry contains
*                       the parametrization of the repeted start point.
*                       (if the endpoint is equal to the startpoint of
*                        the interpolation the lenght of the array should
*                        be equal to inpt1 also in the closed case).
*              inpt1  - Number of original interpolation conditions.
*              ik     - Order of interpolating curve.
*              idim   - Dimension of geometry space.
*              iopen - Indicates if the curve is to be open, closed or
*                       periodic.
*
*
* OUTPUT     : gcond2 - Adjusted interpolation conditions.
*              mtype2 - Type of adusted conditions. See description of
*                       ntype1.
*	       mpar   - The adjusted parameter values.
*              jnpt2  - Number of adusted interpolation conditions.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 91-04.
* REVISED BY : Trond Vidar Stensby, SI, 91-07
*
*********************************************************************
*/
{
   int kstat = 0;
  int kmaxpt = inpt1+ik*(iopen!=SISL_CRV_OPEN);  /* Maximum number of 
						 adjusted conditions. */
  int knpt = 0;			/* Current number of adjusted conditions. */
  int ki, kj, kl;		/* Counters.                              */
  int ktype;			/* Kind of interpolation condition.       */
  int kneg;			/* Indicates negative type indicator.     */
  int kder;			/* Order of differentiation.              */
  int lderarray[MAX_SIZE];
  int alloc_needed=FALSE;
  int *lder = SISL_NULL;		/* Kind of derivative.                    */
  double *sdum = SISL_NULL;	        /* Help array.                            */
  double tdist;                 /* Distance between first and last point. */
  double tref;                  /* Referance value.                       */

  if ((sdum = newarray (idim, DOUBLE)) == SISL_NULL)
     goto err101;

  /* Allocate scratch for output arrays. Make sure that the arrays
     are large enough.  */

  *gcond2 = SISL_NULL;
  if ((*gcond2 = newarray (kmaxpt * idim, DOUBLE)) == SISL_NULL)
    goto err101;
  *mtype2 = SISL_NULL;
  if ((*mtype2 = newarray (kmaxpt, INT)) == SISL_NULL)
    goto err101;
  *mpar = SISL_NULL;
  if ((*mpar = newarray (kmaxpt, DOUBLE)) == SISL_NULL)
    goto err101;

  /* Allocate scratch for local arrays.  */

  if (ik > MAX_SIZE)
    {
       if ((lder = new0array (ik, INT)) == SISL_NULL)
         goto err101;
       alloc_needed = TRUE;
    }
  else
    {
       lder = lderarray;
       memzero(lder,MAX_SIZE,INT);
    }

  
  /* Find first positional condition.  */

  for (ki = 0; ki < inpt1; ki++)
    if (ntype1[ki] == 0)
      break;

  lder[0] = 1;
  (*mtype2)[0] = 0;
  (*mpar)[0] = epar[ki];
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
      (*mpar)[knpt] = epar[kj];
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
      (*mpar)[knpt] = epar[ki];
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
	  (*mpar)[knpt] = epar[kj];
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
	   (*mpar)[knpt] = epar[kj];
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
	   (*mpar)[knpt] = epar[inpt1];
	   memcopy ((*gcond2) + knpt * idim, econd1 + kj * idim, idim, DOUBLE);
	   knpt++;
	}
	
	/* Find first interpolation point and copy it. */
	
	for (kj = 0; ntype1[kj] != 0; kj++) ;
	for (ki = inpt1 - 1; ntype1[ki] != 0; ki--) ;
	
	(*mtype2)[knpt] = 0;
	
	/* UJK & CBI      
	   Replacing     (*mpar)[knpt] = epar[ki] + epar[inpt1]; 
	   with                                                   */
	
	(*mpar)[knpt] = epar[inpt1]; 
	
	memcopy ((*gcond2) + knpt * idim, econd1 + kj * idim, idim, DOUBLE);
	knpt++;
     }
  }

  /* Make sure that the last interpolation conditions are
     of decreasing order of interpolation. */

  kneg = 1;
  for (kder = 0, ki = knpt - 1; ki >= 0 && (kneg ||
					    (*mtype2)[ki] < 0); ki--, kder--)
    {
      if ((*mtype2)[ki] != kder)
	{
	  for (kj = ki - 1; kj >= 0 && (kneg || (*mtype2)[kj] < 0) &&
	       (*mtype2)[kj] != kder; kj--) ;
	  if ((*mtype2)[kj] == kder)
	    {
	      /* Interchange interpolation conditions. */

	      memcopy (sdum, (*gcond2) + kj * idim, idim, DOUBLE);
	      memcopy ((*gcond2) + kj * idim, (*gcond2) + ki * idim,
		       idim, DOUBLE);
	      memcopy ((*gcond2) + ki * idim, sdum, idim, DOUBLE);
	      (*mtype2)[kj] = (kder == 0) ? -(*mtype2)[ki] : (*mtype2)[ki];
	      (*mtype2)[ki] = kder;
	    }
	}
      kneg = 0;
    }

/* UJK & CBI: The following additional parameter value is not used

  if (iopen != SISL_CRV_OPEN)
    {
      (*mpar)[knpt] = (*mpar)[knpt - 1];
      kpar = knpt + 1;
    }
  else
    kpar = knpt;

*/

  /* Conditions adjusted.  */

  *jnpt2 = knpt;
  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

  err101:
    *jstat = -101;
    s6err ("s1908", *jstat, 0);
    goto out;

  out:
    /* Free scratch occupied by local array. */

    if (alloc_needed) freearray (lder);
    if(sdum != SISL_NULL)  freearray(sdum);
    return;
}
