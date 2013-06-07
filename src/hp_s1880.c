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
 * $Id: hp_s1880.c,v 1.5 2003-01-10 12:53:36 vsk Exp $
 *
 */


#define HP_S1880

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    hp_s1880(SISLObject * po1, SISLObject * po2,
	     int ideg,
	     int ipar1, int ipar2, SISLIntdat *pintdat,
	     int *jpar, double **gpar1, double **gpar2, int **pretop,
	     int *jcrv, SISLIntcurve *** wcrv, int *jsurf, SISLIntsurf *** wsurf,
	     int *jstat)
#else
void
   hp_s1880(po1, po2, ideg, ipar1, ipar2, pintdat, jpar, gpar1, gpar2, pretop, jcrv,
       wcrv, jsurf, wsurf, jstat)
     SISLObject *po1;
     SISLObject *po2;
     int ideg;
     int ipar1;
     int ipar2;
     SISLIntdat *pintdat;
     int *jpar;
     double **gpar1;
     double **gpar2;
     int **pretop;
     int *jcrv;
     SISLIntcurve ***wcrv;
     int *jsurf;
     SISLIntsurf ***wsurf;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Transform intersection points and curves from internal
*              format in the recursive part of intersection routines
*              to output format.
*
*
*
* INPUT      : po1    - Pointer first object.
*              po2    - Pointer second object.
*              ideg   - Type of implicit geometry.
*              ipar1  - Number of parameter directions of first object.
*              ipar2  - Number of parameter directions of second object.
*              pintdat - SISLIntdat object.
*
*
* OUTPUT     : jpar   - Number of single intersection points.
*              gpar1  - Parameter values of the single intersection points
*                       in the parameter area of the first object.
*              gpar2  - Parameter values of the single intersection points
*                       in the parameter area of the second object.
*              pretop - Array of pretopology information for the points.
*              jcrv   - Number of intersection curves.
*              wcrv   - Array containing description of intersection curves.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : newIntcurve - Create a new instance of Intcurve.
*              freeIntpt   - Free space occupied by intersection point.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
* REWRITTEN BY : Ulf J. Krystad, SI, 91-06
* REWRITTEN BY : Arne Laksaa & Michael Floater, SI, 91-07
*                 Use pintdat directly as input. Fetch pretop info.
* REWRITTEN BY : UJK, SI, 93-01
*                Exact curve treatment.
*********************************************************************
*/
{
  int kstat;			/* Local status                                */
  int kpos = 0;			/* Position of error.                          */
  int ki, kj, kk;		/* Counters.                                   */
  int kpoint;			/* Number of points in an intersection list.   */
  int ktype;			/* Kind of intersection curve. (See SISLIntcurve). */
  int kpt;			/* Used to find number of single intersection points.*/
  int index;			/* Array index for next point at start        */
  int sing_1, sing_2;		/* Sing point flag at start and end           */
  double *spar1, *spar2;	/* Values of points belonging to an intersection
				     curve in the parameter area of the objects
				     involved in the intersection.                 */
  double *stpar1, *stpar2, *stpar3;	/* Pointers used to travers arrays
				           containing parameter values.        */
  int *top1;	                /* Pointers used to travers pretop. */
  SISLIntcurve **ucrv;		/* Pointer used to traverse *wcrv array.     */
  SISLIntsurf **usurf;		/* Pointer used to traverse *wsurf array.    */
  SISLIntpt *qpt;		/* Pointer to an intersection point.         */
  SISLIntpt *qprev;		/* Pointer to an intersection point.         */
  SISLIntpt *qnext;		/* Pointer to an intersection point.         */
  SISLIntpt *qpfirst;		/* Pointer to first intersection point.      */
  SISLIntpt *qplast;		/* Pointer to last intersection point.       */

  int jpt = pintdat->ipoint;
  SISLIntpt **vpoint = pintdat->vpoint;
  int jlist = pintdat->ilist;
  SISLIntlist **vlist = pintdat->vlist;
  SISLObject *qo2 = SISL_NULL;
  int kdir,kdir1=-1,kdir2=-1;
  int exact=FALSE, exact_treat=FALSE;
  int log_test = 0;
  double dummy;
  /* ------------------------------------------------------------------ */

  for (ki = 1; ki < 5; ki++)
     log_test |= 1 << ki;

  if (po1->iobj == SISLSURFACE && ideg != 0)
  {
     /* Surf vs Implicit geometry */
     qo2 = newObject (SISLPOINT);
     exact_treat = TRUE;
  }
  else  if (po1->iobj == SISLSURFACE &&
	    po2->iobj == SISLSURFACE &&
	    ideg ==0)
  {
     /* Surf vs Implicit geometry */
     qo2 = po2;
     exact_treat = TRUE;
  }


  /* Initiate output.  */

  *gpar1 = *gpar2 = SISL_NULL;
  *wcrv = SISL_NULL;
  *wsurf = SISL_NULL;

  *jcrv = 0;
  *jsurf = 0;

  /* Allocate space for intersection curve and surface array.  */

  *wcrv = newarray (jlist, SISLIntcurve *);
  if (jlist > 0 && *wcrv == SISL_NULL)
    goto err101;
  *wsurf = newarray (jlist, SISLIntsurf *);
  if (jlist > 0 && *wcrv == SISL_NULL)
    goto err101;

  /* Transfer curve-information from vlist array to wcrv and wsurf arrais. */

  ucrv  = *wcrv;
  usurf = *wsurf;

  for (kpt = ki = 0; ki < jlist; ki++)
  {
     qpfirst = qpt = (*vlist)->pfirst;
     qplast = (*vlist)->plast;
     index = (*vlist)->ind_first;
     kpoint = (*vlist)->inumb;
     if (kpoint == 0)
	goto err137;

     if (qpfirst->iinter == SI_TRIM && qpfirst == qplast)
     {
	/* Create new intersection surf.  */

	*usurf = newIntsurf(*vlist);
	if (*usurf == SISL_NULL)
	   goto err101;

	/* Copy pretopology
	   memcopy((*usurf)->pretop,(*vlist)->pretop,4,int); */

	kpt += kpoint-1;
	usurf++;
	(*jsurf)++;
     }
     else
     {
	if (qpfirst->iinter == SI_SING ||
	    (sh6nmbmain (qpfirst,&kstat)) > 2)
	   sing_1 = TRUE;
	else
	   sing_1 = FALSE;

	if (qplast->iinter == SI_SING ||
	    (sh6nmbmain (qplast,&kstat)) > 2)
	   sing_2 = TRUE;
	else
	   sing_2 = FALSE;


	/* Allocate space for arrays containing parameter values of points
	   in intersection curves.                                          */

	spar1 = newarray (ipar1 * kpoint, double);
	spar2 = newarray (ipar2 * kpoint, double);
	if ((ipar1 > 0 && spar1 == SISL_NULL) ||
	    (ipar2 > 0 && spar2 == SISL_NULL))
	   goto err101;

	/* Collect parameter values of the points in this intersection list
	   and distribute values to the objects in the intersection.         */

	kj = 0;
	stpar1 = spar1;
	stpar2 = spar2;
	while (qpt != SISL_NULL && kj < kpoint)
	{
	   stpar3 = qpt->epar;
	   for (kk = 0; kk < ipar1; kk++)
	      *(stpar1++) = *(stpar3++);
	   for (kk = 0; kk < ipar2; kk++)
	      *(stpar2++) = *(stpar3++);

	   /* Reduce no of single points */
	   if (qpt->marker != -99)
	   {
	      kpt++;

	      /* Flag point */
	      qpt->marker = -99;
	   }
	   if (qpt == qpfirst)
	   {
	      qprev = qpt;
	      qpt = qpt->pnext[index];
	   }
	   else
	   {
	      sh6getother (qpt, qprev, &qnext, &kstat);
	      qprev = qpt;
	      qpt = qnext;
	   }
	   kj++;
	}

	/* Find type of intersection curve.  */

	if (sing_1 && sing_2)
	   /* Both ends junction */
	   ktype = 7;
	else if (qpfirst == qplast)
	   /* Closed curve, not singular */
	   ktype = 2;
	else if (sing_1)
	   /* Junction at start */
	   ktype = 5;
	else if (sing_2)
	   /* Junction at end */
	   ktype = 6;
	else
	   /* Open and clean */
	   ktype = 4;


	exact = FALSE;
	/* UJK, March 1995, when curve type is 9 and the curve is
	   an iso-line in both surfaces, we want to return both
	   ppar1 and ppar2. The logic is simple:
	   kdir1 > -1 (ie eq 0 or 1) means constant in 1. surf.
	   kdir2 > -1 (ie eq 2 or 3) means constant in 2. surf. 
	   The object space curve pgeom is picked from the 1. surf
	   if kdir1 is set and from 2. surf if only kdir2 is set. */
	
	/* UJK, January 1993, if exact curve mark it with type 9. */
	kdir1 = kdir2 = -1;
	if (exact_treat &&
	    kj == 2 &&
	    (qpfirst->curve_dir[(*vlist)->ind_first] & log_test))
	{
	   /* Constant parameter curve */
	   for (kdir = 0; kdir < qpfirst->ipar; kdir++)
	      if (qpfirst->curve_dir[(*vlist)->ind_first] &
		  (1 << (kdir + 1)))
	      {
		 exact = TRUE;
		 ktype = 9;
		 if (kdir >= po1->iobj) kdir2 = kdir;
		 else                   kdir1 = kdir;
	      }
	}
	
	if (kdir1>-1) kdir = kdir1;
	else          kdir = kdir2;
	
	/* Create new intersection curve.  */
	*ucrv = newIntcurve (kj, ipar1, ipar2, spar1, spar2, ktype);
	if (*ucrv == SISL_NULL)
	   goto err101;
	
	/* Copy pretopology */
	memcopy((*ucrv)->pretop,(*vlist)->pretop,4,int);
	
	
	/* UJK, January 1993, if exact curve mark it with type 9. */
	if (exact)
	{
	   
	   pick_crv_sf (po1, qo2, kdir, qpfirst, 
			qplast, &(*ucrv)->pgeom, &kstat);
	   if (kstat < 0)
	      goto error;
	   
	   /* UJK, Pick 2D line  */
	   
	   if (kdir2 >= po1->iobj)
	   {
	      s1602(&(qpfirst->epar[po1->iobj]),
		    &(qplast->epar[po1->iobj]),
		    2,
		    2,
		    (*ucrv)->pgeom->et[(*ucrv)->pgeom->ik - 1],
		    &dummy,
		    &(*ucrv)->ppar2,
		    &kstat);
	      if (kstat < 0) goto error;
	   }
	   
	   if (kdir1 >= 0)
	   {
	      s1602(qpfirst->epar,
		    qplast->epar,
		    2,
		    2,
		    (*ucrv)->pgeom->et[(*ucrv)->pgeom->ik - 1],
		    &dummy,
		    &(*ucrv)->ppar1,
		    &kstat);
	      
	      if (kstat < 0) goto error;
	   }
	   
	}
	
	
	ucrv++;
	(*jcrv)++;
     }
     vlist++;
  }

  /* Find number of single intersection points.  */

  kpt = jpt - kpt;
  if (kpt < 0) goto err137;

  /* Create arrays to keep parameter values of intersection points.  */

  *gpar1 = newarray (ipar1 * kpt, double);
  *gpar2 = newarray (ipar2 * kpt, double);
  *pretop = newarray (4 * kpt, int);
  if ((ipar1 * kpt > 0 && *gpar1 == SISL_NULL)
      || (ipar2 * kpt > 0 && *gpar2 == SISL_NULL)
      || (4 * kpt > 0 && *pretop == SISL_NULL))
    goto err101;

  /* Copy parameters of single intersection points into output-arrays. */

  kj = 0;
  stpar1 = *gpar1;
  stpar2 = *gpar2;
  top1 = *pretop;
  for (ki = 0; ki < jpt; ki++)
    {
      qpt = *vpoint;
      if (qpt != SISL_NULL)
	{
	  if (sh6ismain(qpt) && qpt->marker != -99)
	    {
	      kj++;
	      stpar3 = qpt->epar;
	      for (kk = 0; kk < ipar1; kk++)
		*(stpar1++) = *(stpar3++);
	      for (kk = 0; kk < ipar2; kk++)
		*(stpar2++) = *(stpar3++);
              *(top1++) = qpt->left_obj_1[0];
              *(top1++) = qpt->right_obj_1[0];
              *(top1++) = qpt->left_obj_2[0];
              *(top1++) = qpt->right_obj_2[0];
	    }
	}

      vpoint++;
    }

  *jpar = kj;

  /* Adjust output arrays to correct length.  */

  if ((*jcrv) < jlist)
  {
     if ((*jcrv) > 0)
     {
        if (((*wcrv) = increasearray (*wcrv, *jcrv, SISLIntcurve *)) == SISL_NULL)
           goto err101;
     }
     else
     {
        if (*wcrv != SISL_NULL)
	freearray (*wcrv);
        *wcrv = SISL_NULL;
     }
  }
  if ((*jsurf) < jlist)
  {
     if ((*jsurf) > 0)
     {
        if (((*wsurf) = increasearray (*wsurf, *jsurf, SISLIntsurf *)) == SISL_NULL)
           goto err101;
     }
     else
     {
        if (*wsurf != SISL_NULL)
	freearray (*wsurf);
        *wsurf = SISL_NULL;
     }
  }
  if (kj * ipar1 > 0)
    {
      if ((*gpar1 = increasearray (*gpar1, kj * ipar1, double)) == SISL_NULL)
	goto err101;
    }
  else
    {
      if (*gpar1 != SISL_NULL)
	freearray (*gpar1);
      *gpar1 = SISL_NULL;
    }
  if (kj * ipar2 > 0)
    {
      if ((*gpar2 = increasearray (*gpar2, kj * ipar2, double)) == SISL_NULL)
	goto err101;
    }
  else
    {
      if (*gpar2 != SISL_NULL)
	freearray (*gpar2);
      *gpar2 = SISL_NULL;
    }
  if (kj  > 0)
    {
      if ((*pretop= increasearray (*pretop, kj * 4, int)) == SISL_NULL)
	goto err101;
    }
  else
    {
      if (*pretop != SISL_NULL)
	freearray (*pretop);
      *pretop = SISL_NULL;
    }

  /* Intersections copied to output format.  */

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

err101:*jstat = -101;
  s6err ("hp_s1880", *jstat, kpos);
  goto out;

  /* Error in data-strucuture. Expected intersection point not found. */

err137:*jstat = -137;
  s6err ("hp_s1880", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */
error:
  *jstat = kstat;
  s6err ("hp_s1880", *jstat, kpos);
  goto out;


out:
   if (po1->iobj == SISLSURFACE && ideg != 0)
      freeObject (qo2);
return;
}
