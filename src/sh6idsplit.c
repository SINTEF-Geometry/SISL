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
 * $Id: sh6idsplit.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define S6IDSPLIT


#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    sh6idsplit (SISLIntdat ** pintdat, SISLIntpt * psource, int *jstat)
#else
void
   sh6idsplit (pintdat, psource, jstat)
     SISLIntdat **pintdat;
     SISLIntpt *psource;
     int *jstat;
#endif


/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To split an intersection point into several instanses so
*              that the each of the instances has exactly one (main) neighbour.
*              
*
*
* INPUT:       psource  - Pointer to an intersection point.
* 
* 
* INPUT/OUTPUT:pintdat  - Pointer to a pointer to intersection data.
*
*
* OUTPUT  :    jstat    - status messages
*                               = 0      : OK
*                               < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*   
*
* WRITTEN BY : Ulf J. Krystad, 04.92.
*
*********************************************************************
*/
{
  int ki;			/* Counters.                         */
  int no_main;			/* No of neighbours (main points)    */
  int test= FALSE;              /* No equality testing when inserted
				   in pintdat                        */
  int kstat = 0;                /* Local status.                     */
  SISLIntpt *pneighb = SISL_NULL;	/* Current neighbour                 */
  SISLIntpt *pshadow = SISL_NULL;	/* Current copy of source point      */
  /* ------------------------------------------------*/
  
  *jstat = 0;
  
  if (psource == SISL_NULL)
    {
       *jstat = 1;
       goto out;
    }
  
  /* Get number of neighbours */
  no_main = sh6nmbmain (psource, &kstat);
  if (kstat < 0)
    goto error;
  
  for (ki=psource->no_of_curves - 1; no_main > 1; ki--)
    {
       pneighb = sh6getnext(psource, ki);
       if (!pneighb) goto error;
       if (sh6ismain(pneighb))
	 {
	    pshadow = hp_copyIntpt(psource);
	    sh6idnpt(pintdat, &pshadow, test=FALSE, &kstat);
	    if (kstat < 0) goto error;
	    
	    sh6insertpt(psource, pneighb, pshadow, &kstat);
	    if (kstat < 0) goto error;
	    
	    sh6disconnect(psource, pshadow, &kstat);
	    if (kstat < 0) goto error;
	    no_main--;
	 }
    }
  goto out;
  
  
error:
  *jstat = kstat;
  goto out;

out:;
}
