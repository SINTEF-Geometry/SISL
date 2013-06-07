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
 * $Id: sh1782.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH1782

#include "sislP.h"



#if defined(SISLNEEDPROTOTYPES)
static void
sh1782_s9sf_pt (SISLObject *, SISLObject *, double, SISLIntdat **,
		SISLIntpt **, int, int, int *);
static void
sh1782_s9sf_cu (SISLObject *, SISLObject *, double, SISLIntdat **,
		SISLIntpt **, int, int, int *);
static void
sh1782_s9sf_sf (SISLObject *, SISLObject *, double, SISLIntdat **,
		SISLIntpt **, int, int, int *);
#else
static void sh1782_s9sf_pt ();
static void sh1782_s9sf_cu ();
static void sh1782_s9sf_sf ();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
sh1782 (SISLObject * po1, SISLObject * po2, double aepsge,
	SISLIntdat * pintdat, int ipar, double apar,
	SISLIntdat ** rintdat, int *jnewpt, int *jstat)
#else
void
sh1782 (po1, po2, aepsge, pintdat, ipar, apar, rintdat, jnewpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat *pintdat;
     int ipar;
     double apar;
     SISLIntdat **rintdat;
     int *jnewpt;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Insert all points in one intersection data structure
*              into another intersection data structure with one more
*              parameter direction.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              aepsge   - Geometry resolution.
*              pintdat  - Intersection data structure of a problem in a
*                         reduced parameter space.
*              ipar     - Missing parameter direction in the reduced space.
*              apar     - Value of the missing parameter for all intersection
*                         points.
*
*
* OUTPUT     : rintdat  - Intersection data of object-object intersection.
*              jnewpt   - Number of new intersection points.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : sh1782_s9sf_pt  - Set edge info 1D and 2D surface-point.
*              sh1782_s9sf_cu  - Set edge info 3D surface-curve.
*              sh1782_s9sf_sf  - Set edge info 3D surface-surface.
*              sh6idput        - Put point into higher dimensional parameter
*                                space.
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, 06.91
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                     */
  int kdim;			/* Dimension of geometry space.         */
  int knpoint;			/* Number of int.pt. in array.          */
  SISLIntpt **uintpt = SISL_NULL;	/* Array storing intersection points.   */

  *jnewpt = 0;

  /* Test if an intersection data structure exist.  */

  if (pintdat == SISL_NULL)
    {
      *jstat = 0;
      goto out;
    }

  /* Fetch dimension of geometry space. */

  if (po1->iobj == SISLPOINT)
    kdim = po1->p1->idim;
  else if (po1->iobj == SISLCURVE)
    kdim = po1->c1->idim;
  else
    kdim = po1->s1->idim;

  /* Insert the intersection points from pintdat into the structure
     *rintdat belonging to the current problem.  */

  sh6idput (po1, po2, rintdat, pintdat, ipar, apar, &uintpt, &knpoint, &kstat);
  if (kstat < 0)
    goto error;

  if (knpoint == 0)
    {
      *jstat = 0;
      goto out;
    }


  /* Browse on the dimension of geometry space and the type of
     the input objects.     */

  if (kdim <= 2 && ((po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT)
		    || (po2->iobj == SISLSURFACE && po1->iobj == SISLPOINT)))
    {
      /* Compute pre-topology in one-dimensional surface-level value
         intersection.           */

      sh1782_s9sf_pt (po1, po2, aepsge, rintdat, uintpt, knpoint,
		      ipar, &kstat);
      if (kstat < 0)
	goto error;

    }

  else if (kdim == 3 && 
	 ((po1->iobj == SISLSURFACE && po2->iobj == SISLCURVE) ||
	  (po1->iobj == SISLCURVE && po2->iobj == SISLSURFACE )))
    {
      /* Surface-surface intersection in 3-dimensional geometry space. */

      sh1782_s9sf_cu (po1, po2, aepsge, rintdat, uintpt, knpoint,
		      ipar, &kstat);
      if (kstat < 0)
	goto error;
    }

  else if (kdim == 3 && po1->iobj == SISLSURFACE
	   && po2->iobj == SISLSURFACE)
    {
      /* Surface-surface intersection in 3-dimensional geometry space. */

      sh1782_s9sf_sf (po1, po2, aepsge, rintdat, uintpt, knpoint,
		      ipar, &kstat);
      if (kstat < 0)
	goto error;
    }

  /* Task performed.  */

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

error:*jstat = kstat;
  goto out;

out:if (uintpt != SISL_NULL)
    freearray (uintpt);
}




#if defined(SISLNEEDPROTOTYPES)
static void
sh1782_s9sf_pt (SISLObject * po1, SISLObject * po2, double aepsge,
		SISLIntdat ** rintdat, SISLIntpt ** uintpt, int kpoint,
		int ipar, int *jstat)
#else
static void
sh1782_s9sf_pt (po1, po2, aepsge, rintdat, uintpt, kpoint, ipar, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **rintdat;
     SISLIntpt **uintpt;
     int kpoint;
     int ipar;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology data in 1-dimensional surface-point
*              intersection.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              aepsge   - Geometry resolution.
*              rintdat  - Intersection data structure.
*              uintpt   - Intersection points.
*              knpoint  - No of int points
*              ipar     - Parameter direction of missing parameter in the
*                         lower dimensional parameter space.
*
* OUTPUT     : uintpt   - Updated intersection point.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
*
* METHOD     :
*
* CALLS      : sh6getgeom -  Get geometry of intersection point.
*              sh6gettop  -  Get topology of intersection point.
*              sh6settop  -  Set topology of intersection point.
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway. 07.91.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int kdim;			/* Dimension of geometry space.            */
  double tdum, tsign;		/* Parameters used to browse between cases.*/
  double *sder;			/* Surface geometry values.                */
  double *sdum;			/* Normal of surface(dummy)                */
  int ind;			/* Index for curve's par val.              */
  int ki, kj, klow, khigh;	/* Help indexes into uintpt                */
  int index1, index2;		/* Dummy in this context                   */
  int kant;			/* Edge indicator                          */
  int *edge_f, *edge_l;		/* Pointers to atribute edge in intpt      */
  SISLSurf *ps;			/* The 1D Surface                          */
  /* --------------------------------------------------------------------- */

  /* Test input */
  kdim = (po1->iobj == SISLSURFACE) ? po1->s1->idim : po2->s1->idim;
  if (kdim != 1 && kdim !=2)
    goto err106;

  if (ipar == 0)
    ind = 1;
  else
    ind = 0;

  if (po1->iobj == SISLSURFACE)
    tsign = (double) 1.0;
  else
    tsign = -(double) 1.0;

  /* Loop for all points */
  for (ki = 0; ki < kpoint; ki++)
    for (kj = ki + 1; kj < kpoint; kj++)
      /* Connected ? */
      {
	sh6getlist (uintpt[ki], uintpt[kj], &index1, &index2, &kstat);
	if (kstat < 0)
	  goto error;
	if (!kstat)
	  {
	    /* They are connected, sort on edge par val*/
	    if (uintpt[ki]->epar[ind] <
		uintpt[kj]->epar[ind])
	      {
		klow = ki;
		khigh = kj;
	      }
	    else
	      {
		klow = kj;
		khigh = ki;
	      }

	    /* Fetch geometry of surface */
	    sh6getgeom ((po1->iobj == SISLSURFACE) ? po1 : po2,
			(po1->iobj == SISLSURFACE) ? 1 : 2,
			uintpt[klow], &sder, &sdum, aepsge, &kstat);

	    if (kstat < 0)
	      goto error;

	    if (ipar == 0)
	      tdum = tsign * sder[1];
	    else
	      tdum = -tsign * sder[2];

	    if (tdum > DZERO)
	      sh6setdir (uintpt[klow], uintpt[khigh], &kstat);
	    else if (tdum < 0)
	      sh6setdir (uintpt[khigh], uintpt[klow], &kstat);

	      sh6setcnsdir (uintpt[klow], uintpt[khigh], ipar, &kstat);

	    /* Mark edge intersection curve */
	    if (sh6ismain (uintpt[ki]) && sh6ismain (uintpt[kj]) &&
		po1->o1 == po1 &&
		po2->o1 == po2)
	      {

		ps = po1->s1;
		edge_f = &(uintpt[ki]->edge_1);
		edge_l = &(uintpt[kj]->edge_1);

		if (po2->iobj == SISLSURFACE)
		  {
		    edge_f = &(uintpt[ki]->edge_2);
		    edge_l = &(uintpt[kj]->edge_2);
		    ps = po2->s1;
		  }

		kant = 0;
		if (ipar == 0)
		  {
		    if (DEQUAL (uintpt[ki]->epar[ipar],
				ps->et1[ps->ik1 - 1]))
		      kant = -1;
		    else if (DEQUAL (uintpt[ki]->epar[ipar],
				     ps->et1[ps->in1]))
		      kant = 1;
		  }
		else
		  {
		    if (DEQUAL (uintpt[ki]->epar[ipar], ps->et2[ps->ik2 - 1]))
		      kant = 1;
		    else if (DEQUAL (uintpt[ki]->epar[ipar],
				     ps->et2[ps->in2]))
		      kant = -1;
		  }

		if (kant)
		  {
		    if ((double) kant * tdum > DZERO)
		      *edge_f = *edge_l = SI_RIGHT;
		    else if ((double) kant * tdum < DZERO)
		      *edge_f = *edge_l = SI_LEFT;

		  }

	      }

	  }
      }

  *jstat = 0;
  goto out;

  /* Error in input. Incorrect dimension.  */

err106:*jstat = -106;
  goto out;

  /* Error lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
  return;
}


#if defined(SISLNEEDPROTOTYPES)
static void
sh1782_s9sf_cu (SISLObject * po1, SISLObject * po2, double aepsge,
		SISLIntdat ** rintdat, SISLIntpt ** uintpt, int kpoint,
		int ipar, int *jstat)
#else
static void
sh1782_s9sf_cu (po1, po2, aepsge, rintdat, uintpt, kpoint, ipar, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **rintdat;
     SISLIntpt **uintpt;
     int kpoint;
     int ipar;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology data in 3-dimensional surface-curve
*              intersection.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              aepsge   - Geometry resolution.
*              rintdat  - Intersection data structure.
*              uintpt   - Intersection points.
*              kpoint   - No of int points
*              ipar     - Parameter direction of missing parameter in the
*                         lower dimensional parameter space.
*
* OUTPUT     : uintpt   - Updated intersection point.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : s6crss    -  Cross product between two vectors.
*              s6scpr    -  Scalar product between two vectors.
*              sh6getgeo -  Fetch geometry information of intersection point.
*
* REFERENCES :
*
* WRITTEN BY : ALA & VSK okt. 1992.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int kdim;			/* Dimension of geometry space.            */
  double tdum;			/* Parameter used to browse on the cases.  */
  double *val_s1;		/* Surface geometry                        */
  double *val_s2;		/* Surface geometry                        */
  double *snorm1;		/* Normal vector of surface.               */
  double *snorm2;		/* Normal vector of surface.               */
  double *sdir;	        	/* Direction vector of intersection curve. */
  double *stang;		/* Tangent vector of edge curve.           */
  int ind;			/* Index for curve's par val.              */
  int ki, kj, klow, khigh;	/* Help indexes into uintpt                */
  int index1, index2;		/* Dummy in this context                   */
  int kcurve;			/* Which object is the curve?		   */
  /* --------------------------------------------------------------------- */
  
  *jstat = 0;

  /* Test input */
  
  if (po1->iobj == SISLCURVE && po2->iobj == SISLSURFACE)
  {
     kcurve = 1;
     kdim = po1->c1->idim;
     if (kdim != 3)
	goto err104;
     if (kdim != po2->s1->idim)
	goto err106;
  }
  else if (po1->iobj == SISLSURFACE && po2->iobj == SISLCURVE)
  {
     kcurve = 2;
     kdim = po1->s1->idim;
     if (kdim != 3)
	goto err104;
     if (kdim != po2->c1->idim)
	goto err106;
  }
  else
    goto err104;
  
  /* Get index for curve par val  */
  if (kcurve == 1 && ipar == 1)
    ind = 2;
  else  if ((kcurve == 1 && ipar == 2) ||
	    (kcurve == 2 && ipar == 0))
    ind = 1;
  else  if (kcurve == 2 && ipar == 1)
    ind = 0;
  else
     goto out;
  
  /* Loop for all points */
  for (ki = 0; ki < kpoint; ki++)
    for (kj = ki + 1; kj < kpoint; kj++)
      /* Connected ? */
      {
	sh6getlist (uintpt[ki], uintpt[kj], &index1, &index2, &kstat);
	if (kstat < 0)
	  goto error;
	if (!kstat)
	  {
	    /* They are connected, sort on edge par val*/
	    if (uintpt[ki]->epar[ind] <
		uintpt[kj]->epar[ind])
	      {
		klow = ki;
		khigh = kj;
	      }
	    else
	      {
		klow = kj;
		khigh = ki;
	      }

	    /* Get geometry first object */
	    sh6getgeom (po1, 1, uintpt[klow], &val_s1, &snorm1, aepsge, &kstat);
	    if (kstat < 0)
	      goto error;

	    /* Get geometry second object */
	    sh6getgeom (po2, 2, uintpt[klow], &val_s2, &snorm2, aepsge, &kstat);
	    if (kstat < 0)
	      goto error;

	    /* Fetch tangent of edge curve.  */
	    if (ipar == 0)
	      stang = val_s1 + 6;
	    else if (ipar == 1 && kcurve == 2)
	      stang = val_s1 + 3;
	    else if (ipar == 1 && kcurve == 1)
	      stang = val_s2 + 6;
	    else
	      stang = val_s2 + 3;

	    /* The direction of the intersection curve is set along
	       the direction of the curve in the intersection.  */
	    if (kcurve == 1)
	       sdir = val_s1 + 3;
	    else
	       sdir = val_s2 + 3;

	    tdum = s6scpr (sdir, stang, kdim);

	    if (tdum > 0)
	      sh6setdir (uintpt[klow], uintpt[khigh], &kstat);
	    else if (tdum < 0)
	      sh6setdir (uintpt[khigh], uintpt[klow], &kstat);

	      sh6setcnsdir (uintpt[klow], uintpt[khigh], ipar, &kstat);
	  }
      }
  *jstat = 0;
  goto out;


  /* Error in input. Dimension not equal to 3. */

err104:*jstat = -104;
  goto out;

  /* Error in input. Conflicting dimensions.  */

err106:*jstat = -106;
  goto out;

  /* Error lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
  return;
}



#if defined(SISLNEEDPROTOTYPES)
static void
sh1782_s9sf_sf (SISLObject * po1, SISLObject * po2, double aepsge,
		SISLIntdat ** rintdat, SISLIntpt ** uintpt, int kpoint,
		int ipar, int *jstat)
#else
static void
sh1782_s9sf_sf (po1, po2, aepsge, rintdat, uintpt, kpoint, ipar, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **rintdat;
     SISLIntpt **uintpt;
     int kpoint;
     int ipar;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology data in 3-dimensional surface-surface
*              intersection.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              aepsge   - Geometry resolution.
*              rintdat  - Intersection data structure.
*              uintpt   - Intersection points.
*              knpoint  - No of int points
*              ipar     - Parameter direction of missing parameter in the
*                         lower dimensional parameter space.
*
* OUTPUT     : uintpt   - Updated intersection point.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : s6crss    -  Cross product between two vectors.
*              s6scpr    -  Scalar product between two vectors.
*              sh6getgeo -  Fetch geometry information of intersection point.
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway. 07.91.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int kdim;			/* Dimension of geometry space.            */
  double tdum;			/* Parameter used to browse on the cases.  */
  double *val_s1;		/* Surface geometry                        */
  double *val_s2;		/* Surface geometry                        */
  double *snorm1;		/* Normal vector of surface.               */
  double *snorm2;		/* Normal vector of surface.               */
  double sdir[3];		/* Direction vector of intersection curve. */
  double *stang;		/* Tangent vector of edge curve.           */
  int ind;			/* Index for curve's par val.              */
  int ki, kj, klow, khigh;	/* Help indexes into uintpt                */
  int index1, index2;		/* Dummy in this context                   */
  int *edge_f, *edge_l;		/* Pointers to atribute edge in intpt      */
  int kant;			/* Edge indicator                          */
  SISLSurf *ps;			/* The 3D Surface enhanced                 */
  /* --------------------------------------------------------------------- */

  /* Test input */
  kdim = po1->s1->idim;
  if (kdim != 3)
    goto err104;
  if (kdim != po2->s1->idim)
    goto err106;

  /* Get index for curve par val  */
  if (ipar == 0)
    ind = 1;
  else if (ipar == 1)
    ind = 0;
  else if (ipar == 2)
    ind = 3;
  else
    ind = 2;

  /* Loop for all points */
  for (ki = 0; ki < kpoint; ki++)
    for (kj = ki + 1; kj < kpoint; kj++)
      /* Connected ? */
      {
	sh6getlist (uintpt[ki], uintpt[kj], &index1, &index2, &kstat);
	if (kstat < 0)
	  goto error;
	if (!kstat)
	  {
	    /* They are connected, sort on edge par val*/
	    if (uintpt[ki]->epar[ind] <
		uintpt[kj]->epar[ind])
	      {
		klow = ki;
		khigh = kj;
	      }
	    else
	      {
		klow = kj;
		khigh = ki;
	      }




	    /* Get geometry first object */
	    sh6getgeom (po1, 1, uintpt[klow], &val_s1, &snorm1, aepsge, &kstat);
	    if (kstat < 0)
	      goto error;

	    /* Get geometry second object */
	    sh6getgeom (po2, 2, uintpt[klow], &val_s2, &snorm2, aepsge, &kstat);
	    if (kstat < 0)
	      goto error;

	    /* Get geometry lower level object */
	    /* Fetch tangent of edge curve.  */
	    if (ipar == 0)
	      stang = val_s1 + 6;
	    else if (ipar == 1)
	      stang = val_s1 + 3;
	    else if (ipar == 2)
	      stang = val_s2 + 6;
	    else
	      stang = val_s2 + 3;

	    /* Compute direction of intersection curve.  */
	    s6crss (snorm1, snorm2, sdir);

	    tdum = s6scpr (sdir, stang, kdim);

	    if (tdum > 0)
	      sh6setdir (uintpt[klow], uintpt[khigh], &kstat);
	    else if (tdum < 0)
	      sh6setdir (uintpt[khigh], uintpt[klow], &kstat);

	      sh6setcnsdir (uintpt[klow], uintpt[khigh], ipar, &kstat);
	    /* Mark edge intersection curve */
	    if (sh6ismain (uintpt[ki]) && sh6ismain (uintpt[kj]) &&
		po1->o1 == po1 &&
		po2->o1 == po2)
	      {

		if (ipar < 2)
		  /* First surf */
		  {
		    edge_f = &(uintpt[ki]->edge_1);
		    edge_l = &(uintpt[kj]->edge_1);
		    ps = po1->s1;
		  }
		else
		  {
		    edge_f = &(uintpt[ki]->edge_2);
		    edge_l = &(uintpt[kj]->edge_2);
		    ps = po2->s1;
		  }

		kant = 0;
		if (ipar == 0 || ipar == 2)
		  {
		    if (DEQUAL (uintpt[ki]->epar[ipar],
				ps->et1[ps->ik1 - 1]))
		      kant = -1;
		    else if (DEQUAL (uintpt[ki]->epar[ipar],
				     ps->et1[ps->in1]))
		      kant = 1;
		  }
		else
		  {
		    if (DEQUAL (uintpt[ki]->epar[ipar], ps->et2[ps->ik2 - 1]))
		      kant = 1;
		    else if (DEQUAL (uintpt[ki]->epar[ipar],
				     ps->et2[ps->in2]))
		      kant = -1;
		  }

		if (kant)
		  {
		    if ((double) kant * tdum > DZERO)
		      *edge_f = *edge_l = SI_RIGHT;
		    else if ((double) kant * tdum < DZERO)
		      *edge_f = *edge_l = SI_LEFT;

		  }





	      }
	  }
      }
  *jstat = 0;
  goto out;


  /* Error in input. Dimension not equal to 3. */

err104:*jstat = -104;
  goto out;

  /* Error in input. Conflicting dimensions.  */

err106:*jstat = -106;
  goto out;

  /* Error lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
  return;
}
