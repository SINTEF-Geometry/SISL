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
 * $Id: sh6getmain.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETMAIN

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
SISLIntpt *
sh6getmain (SISLIntpt * pt)
#else
SISLIntpt *
sh6getmain (pt)
     SISLIntpt *pt;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given a help point, find the unique main point it is
*              linked to. If there is no such main point, return
*              SISL_NULL.
*
*
* INPUT      : pt       - Pointer to the Intpt.
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
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
* CHANGED BY: Ulf J. Krystad, SI, Oslo, Norway. September 91.
*********************************************************************
*/
{
  int ki;			/* Loop control */
  int kstat;			/* Local status */
  int more = TRUE;		/* Loop control */
  SISLIntpt *mainpt = SISL_NULL;
  SISLIntpt *pt1 = SISL_NULL;
  SISLIntpt *pt2 = SISL_NULL;
  SISLIntpt *prev = SISL_NULL;
  SISLIntpt *pcurr = SISL_NULL;
  SISLIntpt *pnext = SISL_NULL;
  /* ------------------------------------------------------------- */


  if (!sh6ishelp (pt))
    goto out;

  for (ki = 0; ki < pt->no_of_curves; ki++)
    {
      if (sh6ismain (pt1 = sh6getnext (pt, ki)))
	{
	  mainpt = pt1;
	  break;
	}
    }

  if (!mainpt)
    {
      /* No close neighbour is main, check along list
         if not meeting point. */
      sh6getnhbrs (pt, &pt1, &pt2, &kstat);
      if (kstat == 1)
	{
	  /* Terminator, go towards other end */
	  prev = pt;
	  pcurr = pt1;
	  more = TRUE;

	  while ((!mainpt) && more)
	    {
	      sh6getother (pcurr, prev, &pnext, &kstat);
	      if (kstat < 0)
		goto error;

	      if (pnext && (pnext != pt))
		{
		  if (sh6ismain (pnext))
		    mainpt = pnext;
		  else
		    {
		      prev = pcurr;
		      pcurr = pnext;
		      pnext = SISL_NULL;
		    }
		}
	      else
		more = FALSE;

	    }
	}

      else if (kstat == 0)
	{
	  /* Two neighbours, search both directions */
	  for (ki = 0, prev = pt, pcurr = pt1, more = TRUE; (!mainpt) && (ki < 2);
	       ki++, prev = pt, pcurr = pt2, more = TRUE)

	    while ((!mainpt) && more)
	      {
		sh6getother (pcurr, prev, &pnext, &kstat);
		if (kstat < 0)
		  goto error;

		if (pnext && (pnext != pt))
		  {
		    if (sh6ismain (pnext))
		      mainpt = pnext;
		    else
		      {
			prev = pcurr;
			pcurr = pnext;
			pnext = SISL_NULL;
		      }
		  }
		else
		  more = FALSE;

	      }
	}
    }

  goto out;

  /* ------------------------------------------------------------- */
error:mainpt = SISL_NULL;
  s6err ("sh6getmain", kstat, 0);
  goto out;



out:
  return mainpt;
}

