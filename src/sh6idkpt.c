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
