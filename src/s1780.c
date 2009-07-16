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
 * $Id: s1780.c,v 1.2 2001-03-19 15:58:53 afr Exp $
 *
 */


#define S1780

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1780(SISLCurve *pc1,SISLCurve *pc2,SISLIntpt *vipt[],int *jstat)
#else
void s1780(pc1,pc2,vipt,jstat)
     SISLCurve *pc1;
     SISLCurve *pc2;
     SISLIntpt *vipt[];
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Test if two curves originating from a curve/curve
*              intersection problem, coincide beetween two
*              intersection points. The part of curves beetween
*              the intersection points must not have internal
*              knots.
*
*
* INPUT      : pc1     - First curve in intersection problem.
*              pc2     - Second curve in intersection problem.
*              vipt[]  - Array of structures of intersection points.
*
*
*
* OUTPUT     : jstat  - status messages  
*                                         = 1      : Coincide.
*                                         = 0      : No coincide.
*                                         < 0      : error
*
*
* METHOD     : Let k be the maximum order of the two Bezier-curves.
*              Test if the curves or derivatives of the curves are 
*              identical in k points. Then the curves coincide.
*
*
* REFERENCES :
*
*-
* CALLS      : s1221  - Evaluate curve in given parameter value.
*              s6par  - Test if two vectors are parallel.
*
* WRITTEN BY : Arne Laksaa, SI, 89-07.
*
*********************************************************************
*/
{ 
  int kstat = 0;     /* Local status variable.                            */
  int kpos = 0;      /* Position of error;                                */
  int ki,kj;         /* Counters.                                         */
  int kdim;          /* Dimension of the space in which the curves lie.   */
  int kord;          /* Number of points in which the curves have to be
		        equal if the curves coincide.                     */
  int kder;          /* Number of derivatives of curve to evaluate.       */
  int klefs = 0;     /* Parameter used when evaluating first curve.       */
  int kleft = 0;     /* Parameter used when evaluating second curve.      */
  double *sder1 = SISL_NULL; /* Value of first curve in given parameter value. */
  double *sder2 = SISL_NULL; /* Value of second curve in given parameter value.*/
  double tpar1,tpar2;   /* Parameter values in intersection points.       */
  double tang;          /* An angel beetween to vectors.                  */


  /* Initiate to no coincide. */

  *jstat = 0;

  kdim = pc1 -> idim;
  if (kdim != pc2->idim) goto err106;

  /* Test if Bezier curves beetween intersection points.  */

  tpar1 = min(vipt[0]->epar[0],vipt[1]->epar[0]);
  tpar2 = max(vipt[0]->epar[0],vipt[1]->epar[0]);

  for (ki=0;  pc1->et[ki] <= tpar1; ki++);
  ki--;

  for (kj=0;  pc1->et[kj] < tpar2; kj++);

  if (kj-ki > 1)
    goto out;

  tpar1 = min(vipt[0]->epar[1],vipt[1]->epar[1]);
  tpar2 = max(vipt[0]->epar[1],vipt[1]->epar[1]);

  for (ki=0;  pc2->et[ki] <= tpar1; ki++);
  ki--;

  for (kj=0;  pc2->et[kj] < tpar2; kj++);

  if (kj-ki >1)
    goto out;


  /* Find number of points/derivatives in which the curves have to 
     be equal.                                                     */

  kord = MAX(pc1->ik,pc2->ik);
  kder = kord/2;
  kder = MAX(kder,kord-kder);

  /* Allocate space for local arrays.  */

  if ((sder1 = newarray((2*(kder+1)*kdim),double)) == SISL_NULL) goto err101;
  sder2 = sder1 + (kder+1)*kdim;



  /* Test if curves coincide beetween intersection points.  */

  if (kder > 1)
    {

      /* Test if the first kder derivatives of the curves are equal in
	 the first intersection points of the curves.   */

      /* Evaluate curves in first endpoint.  */

      s1221(pc1,kder,vipt[0]->epar[0],&kleft,sder1,&kstat);
      if (kstat < 0) goto error;

      s1221(pc2,kder,vipt[0]->epar[1],&klefs,sder2,&kstat);
      if (kstat < 0) goto error;

      for (ki=1; ki<kder; ki++)
	{

	  /* Test if the ki'th derivatives of the curves are parallel. */

	  tang = s6ang(sder1+(ki*kdim),sder2+(ki*kdim),kdim);

	  tang = min(tang,fabs(PI - tang));

	  if (tang > ANGULAR_TOLERANCE) goto out;
	}
    }

  kder = kord - kder;
  if (kder > 1)
    {

      /* Test if the first kder derivatives of the curves are equal in 
	 the endpoints of the curves.                                   */

      /* Evaluate curves in second endpoint. */

      s1221(pc1,kder,vipt[1]->epar[0],&kleft,sder1,&kstat);
      if (kstat < 0) goto error;

      s1221(pc2,kder,vipt[1]->epar[1],&klefs,sder2,&kstat);
      if (kstat < 0) goto error;

      for (ki=1; ki<kder; ki++)
	{

	/* Test if the ki'th derivatives of the curves are parallel. */
	  tang = s6ang(sder1+(ki*kdim),sder2+(ki*kdim),kdim);

	  tang = min(tang,fabs(PI - tang));

	  if (tang > ANGULAR_TOLERANCE) goto out;
	}
    }

  /* Test performed.  */

  *jstat = 1;
  goto out;

  /* Error in space allocation.  */

err101: *jstat = -101;
        s6err("s1780",*jstat,kpos);
        goto out;
  
  /* Error in input. Dimensions of curves conflicting. */

err106: *jstat = -106;
        s6err("s1780",*jstat,kpos);
        goto out;

  /* Error in lower level routine. */

error : *jstat = kstat;
        s6err("s1780",*jstat,kpos);
        goto out;

 out:

  /* Free space occupied by local arrays.  */

  if (sder1 != SISL_NULL) freearray(sder1);
}
