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
 * $Id: sh6topohlp.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH6GETTOPHLP

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6gettophlp (SISLIntpt * pt, int pretop[4], int case_2d, int *jstat)
#else
void
sh6gettophlp (pt, pretop, case_2d, jstat)
     SISLIntpt *pt;
     int pretop[4];
     int case_2d;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Get pre-topology information, by traversing all help
*              points connected in a network.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
* INPUT/OUTPUT:
*              pretop   - pre-topology data.
*              case_2d  - flage 2d surf point.
* OUTPUT     : pt       - Pointer to the Intpt.
*              jstat    - Error flag.
*                         jstat =  0  => OK.
*                         jstat = -2  => Error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : UJK, SI, Oslo, Norway. October 91.
*
*********************************************************************
*/
{
  int loc_top[4];
  int ki;

  *jstat = 0;

/* Check pt. */

  if (pt == SISL_NULL)
    goto err2;
/* Only help points are treated */
  if (sh6ishelp (pt) && pt->marker == 0)
    {
      /* To avoid infinite loops : */
      pt->marker = -10;

      sh6gettop (pt, 0, loc_top, loc_top + 1, loc_top + 2, loc_top + 3, jstat);
      if (*jstat < 0)
	goto out;

      if (case_2d)
      {
	 /* Spesial treatment 2D surf point */
	 for (ki=0; ki<4; ki++)
	    if (loc_top[ki] == SI_IN) pretop[ki] = SI_IN;
	    else if (loc_top[ki] == SI_OUT && pretop[ki] != SI_IN)
	       pretop[ki] = SI_OUT;
      }
      else
      {
	 /* Overrule ? */
	 for (ki = 0; ki < 4; ki++)
	    if ((pretop[ki] == SI_UNDEF ||
		 pretop[ki] == SI_ON) &&
		loc_top[ki] != SI_UNDEF &&
		loc_top[ki] != SI_ON)
	       pretop[ki] = loc_top[ki];
      }
      
      for (ki = 0; ki < pt->no_of_curves; ki++)
	sh6gettophlp (pt->pnext[ki],  pretop, case_2d, jstat);

      /* Data is set. */

    }


  goto out;


err2:
  /* Error in input. pt is SISL_NULL. */

  *jstat = -2;
  s6err ("sh6gettophlp", *jstat, 0);
  goto out;


out:
  return;
}
