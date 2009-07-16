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
 * $Id: sh6connect.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6CONNECT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6connect (SISLIntpt * pt1, SISLIntpt * pt2, int *jstat)
#else
void
sh6connect (pt1, pt2, jstat)
     SISLIntpt *pt1;
     SISLIntpt *pt2;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Connect the two points.
*              If they are already connected give message.
*
*
* INPUT      : pt1      - Pointer to first Intpt.
*              pt2      - Pointer to second Intpt.
*              jstat    - Error flag.
*                        jstat =  0  => Successful
*                        jstat =  1  => Points already connected.
*                        jstat = -1  => Illegal to connect.
*                        jstat = -2  => Error in data structure.
*                        jstat = -3  => Error in subfunction.
*                        jstat = -4  => Selfconnecting not legal.
*
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. May 91.
* CORRECTED BY : Ulf J. Krystad, SI, Oslo, Norway. July 91.
*********************************************************************
*/
{
  int kstat;			/* error flag. */
  int index1, index2;		/* dummy indices.           */
  int num;			/* Number of main point pinters.  */

  *jstat = 0;
  
  if (pt1 == pt2)
    goto err4;

  /* Check if pt1 and pt2 are already connected. */

  sh6getlist (pt1, pt2, &index1, &index2, &kstat);
  if (kstat < 0)
    goto err3;
  if (kstat < 1)		/* Already connected. */
    {
      *jstat = 1;
      goto out;
    }

  /* Check that we can connect pt1. There are restrictions if it
     it a help point.  */

  if (sh6ishelp (pt1))		/* pt1 is a help point */
    {
      /* UJK, this is NO invariant */
      /*if (pt1->no_of_curves > 2)
         goto err2;
         if (pt1->no_of_curves == 2)
         goto err1; */

      if (sh6ismain (pt2))	/* pt2 is a main point. */
	{
	  num = sh6nmbmain (pt1, &kstat);

	  /* UJK, If invar does not hold, MAKE it hold */
	  /* if (num > 1)
	    goto err2;
	  if (num == 1)
	    goto err1; */
	  if (num >= 1)
	    sh6tomain (pt1, &kstat);
	  if (kstat < 0)
	    goto err2;

	  /* pt1 cannot be connected to two main points. */
	}
    }

  /* Check that we can connect pt2. There are restrictions if it
     it a help point.  */

  if (sh6ishelp (pt2))		/* pt2 is a help point */
    {
      /* UJK, this is NO invariant */
      /*if (pt2->no_of_curves > 2)
         goto err2;
         if (pt2->no_of_curves == 2)
         goto err1; */

      if (sh6ismain (pt1))	/* pt1 is a main point. */
	{
	  num = sh6nmbmain (pt2, &kstat);

	  /* UJK, If invar does not hold, MAKE it hold */
	  /*if (num > 1)
	    goto err2;
	  if (num == 1)
	    goto err1; */
	  if (num >= 1)
	    sh6tomain (pt2, &kstat);
	  if (kstat < 0)
	    goto err2;

	  /* pt2 cannot be connected to two main points. */
	}
    }

  /* Now make the connection. */


  /* Point pt1 to pt2. */

  /* Check if we need to reallocate the pnext and curve_dir arrays. */

  if (pt1->no_of_curves > pt1->no_of_curves_alloc)
    goto err2;
  if (pt1->no_of_curves == pt1->no_of_curves_alloc)
    {
      pt1->no_of_curves_alloc += 4;
      pt1->pnext = increasearray (pt1->pnext,
				  pt1->no_of_curves_alloc, SISLIntpt *);
      pt1->curve_dir = increasearray (pt1->curve_dir,
				      pt1->no_of_curves_alloc, int);
      /* UJK, Must have size of pretop arrays increased */
      pt1->left_obj_1 = increasearray (pt1->left_obj_1,
				       pt1->no_of_curves_alloc, int);
      pt1->left_obj_2 = increasearray (pt1->left_obj_2,
				       pt1->no_of_curves_alloc, int);
      pt1->right_obj_1 = increasearray (pt1->right_obj_1,
					pt1->no_of_curves_alloc, int);
      pt1->right_obj_2 = increasearray (pt1->right_obj_2,
					pt1->no_of_curves_alloc, int);
    }

  /* Set new pointer to new position in array. */
  /* Set new curve direction to 0 for now. */

  pt1->pnext[pt1->no_of_curves] = pt2;
  pt1->curve_dir[pt1->no_of_curves] = 0;

  /* Increment no_of_curves. */

  pt1->no_of_curves++;


  /* Point pt2 to pt1. */

  /* Check if we need to reallocate the pnext and curve_dir arrays. */

  if (pt2->no_of_curves > pt2->no_of_curves_alloc)
    goto err2;
  if (pt2->no_of_curves == pt2->no_of_curves_alloc)
    {
      pt2->no_of_curves_alloc += 4;
      /* UJK, pt1->pnext chaged to pt2->pnext */
      pt2->pnext = increasearray (pt2->pnext,
				  pt2->no_of_curves_alloc, SISLIntpt *);
      pt2->curve_dir = increasearray (pt2->curve_dir,
				      pt2->no_of_curves_alloc, int);
      /* UJK, Must have size of pretop arrays increased */
      pt2->left_obj_1 = increasearray (pt2->left_obj_1,
				       pt2->no_of_curves_alloc, int);
      pt2->left_obj_2 = increasearray (pt2->left_obj_2,
				       pt2->no_of_curves_alloc, int);
      pt2->right_obj_1 = increasearray (pt2->right_obj_1,
					pt2->no_of_curves_alloc, int);
      pt2->right_obj_2 = increasearray (pt2->right_obj_2,
					pt2->no_of_curves_alloc, int);
    }

  /* Set new pointer to new position in array. */
  /* Set new curve direction to 0 for now. */

  pt2->pnext[pt2->no_of_curves] = pt1;
  pt2->curve_dir[pt2->no_of_curves] = 0;

  /* Increment no_of_curves. */

  pt2->no_of_curves++;



  goto out;

  /* Illegal to connect. */
  /*err1:

  *jstat = -1;
  s6err ("sh6connect", *jstat, 0);
  goto out; */

  /* Error in data structure. */
err2:

  *jstat = -2;
  s6err ("sh6connect", *jstat, 0);
  goto out;

  /* Error in subfunction. */
err3:

  *jstat = -3;
  s6err ("sh6connect", *jstat, 0);
  goto out;

err4:
  /* Selfconnecting not legal */
  *jstat = -4;
  s6err ("sh6connect", *jstat, 0);
  goto out;


out:
  return;
}
