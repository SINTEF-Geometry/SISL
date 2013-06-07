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
 * $Id: sh1780.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH1780

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
sh1780 (SISLObject * po1, SISLObject * po2, double aepsge,
	SISLIntdat ** rintdat, SISLIntpt * pintpt, int *jnewpt, int *jstat)
#else
void
sh1780 (po1, po2, aepsge, rintdat, pintpt, jnewpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **rintdat;
     SISLIntpt *pintpt;
     int *jnewpt;
     int *jstat;
#endif

 /* UPDATE: (ujk) : Must tune the test on when to use local
  info (=subtask no ?),
  As in sh1779 it is necessary to determine when a helppoint
  is close enough to a mainpoint to be neglected. */
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Generate help points and set pre-topology data in
*              curve-curve intersection.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              aepsge   - Geometry resolution.
*              rintdat  - Intersection data structure.
*              pintpt   - Current intersection point.
*
*
* OUTPUT     : pintpt   - Intersection point with updated pre-topology info.
*              jnewpt   - Number of new int.pt. created.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : shevalc    -  ? (s1221 used )Evaluate curve using filtered coefficients.
*              s1421      -  Evaluate surface using filtered coefficients.
*              sh1783     -  March along two curves as long as they coincide.
*              s6ang      -  Angle between two vectors.
*              s6length   -  Length of vector.
*              s6dist     -  Distance between two points.
*              sh6idcon   -  Connect intersection points.
*              sh6genhbrs -  Get neighbours.
*              sh6getgeo  -  Get geometric info from intersection point.
*              sh6settop  -  Set topology of intersection point.
*              hp_newIntpt   -  Create new intersection point.
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SI, Oslo, Norway. 04.91.
*              Ulf J. Krystad, SI, 06.91
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int ki;			/* Counters.                               */
  int kleft1 = 0;               /* Parameters to the evaluator.            */
  int kdim;			/* Dimension of geometry space.            */
  int kpos = 0;			/* Current position in output array.       */
  int kdir1, kdir2;		/* Directions in which to march the curves.*/
  int kk1, kk2;			/* Orders of the two curves.               */
  int kn1, kn2;			/* Number of vertices in the curves.       */
  int lleft[2];			/* Array storing pre-topology information. */
  int lright[2];		/* Array storing pre-topology information. */
  double tref;			/* Reference value in equality test.       */
  double *st1, *st2;		/* Pointers to knot vectors of curves.     */
  double sder[6];		/* Result of curve evaluation.             */
  double stang1[3];		/* Tangent vector of curve.                */
  double stang2[3];		/* Tangent vector of level value.          */
  double slast[3];		/* Last parameter value of coincidence.    */
  double snext[3];		/* First parameter value outside interval
			           of coincidence.                         */
  double *ret_val;		/* Pointer to geo data from sh6getgeom     */
  double *ret_norm;		/* Pointer to geo data from sh6getgeom     */
  double *sptpar = pintpt->epar;/* Parameter array of int.pt.        */
  SISLIntpt *uintpt[2];		/* Pointer to new intersection points.     */
  double *nullp = SISL_NULL;

  /* Don't make pretop for help points ! */
  if (sh6ishelp (pintpt))
    {
      *jstat = 0;
      goto out;
    }

  /* Test dimension of geometry space.  */

  kdim = po1->c1->idim;
  if (kdim > 3)
    goto err108;
  if (kdim != po2->c1->idim)
    goto err106;

  /* Express the curve by local parameters.  */

  kn1 = po1->c1->in;
  kk1 = po1->c1->ik;
  st1 = po1->c1->et;
  kn2 = po2->c1->in;
  kk2 = po2->c1->ik;
  st2 = po2->c1->et;
  tref = MAX (st1[kn1] - st1[kk1 - 1], st2[kn2] - st2[kk2 - 1]);

  /* Fetch already existing topology. */
  sh6gettop (pintpt, -1, lleft, lright, lleft + 1, lright + 1, &kstat);

  /* Fetch geometry information, first curve.  */
  sh6getgeom (po1, 1, pintpt, &ret_val, &ret_norm, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  /* Local copy of curve tangent */
  memcopy (stang1, ret_val + kdim, kdim, DOUBLE);

  /* Fetch geometry information,second curve.  */
  sh6getgeom (po2, 2, pintpt, &ret_val, &ret_norm, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  /* Local copy of curve tangent */
  memcopy (stang2, ret_val + kdim, kdim, DOUBLE);

  /* Compute the angle between the tangent vectors of the curves
     in the current intersection point, and check if marching is
     necessary to compute the pre-topology information.  */

  /* UPDATE (ujk) : tune */
  if (s6ang (stang1, stang2, kdim) <= ANGULAR_TOLERANCE)
    {
      /* Perform marching in positive direction of the first curve.  */

      kdir1 = 1;
      kdir2 = (s6scpr (stang1, stang2, kdim) >= DZERO) ? 1 : -1;

      /* Check if the intersection point is situated at the endpoint
	 of a curve.             */

      if (DEQUAL (sptpar[0] + tref, st1[kn1] + tref) ||
	  (kdir2 == 1 && DEQUAL (sptpar[1] + tref, st2[kn2] + tref)) ||
	  (kdir2 == -1 && DEQUAL (sptpar[1] + tref, st2[kk2 - 1] + tref)))
	{
	}
      else
	{
	  /* Perform marching.  */

	  sh1783 (po1->c1, po2->c1, aepsge, sptpar, kdir1, kdir2, slast,
		  snext, &kstat);
	  if (kstat < 0)
	    goto error;

	  if (kstat > 0)
	    {
	      /* An intersection interval is found. */
	      /* Set pre-topology */

	      lright[0] = SI_ON;
	      if (kdir2 == 1)
		lright[1] = SI_ON;
	      else
		lleft[1] = SI_ON;
	    }
	  else
	    {
	      /* Create help point. First fetch geometry information. */

	      s1221 (po1->c1, 0, slast[0], &kleft1, sder, &kstat);
	      if (kstat < 0)
		goto error;

	      s1221 (po1->c1, 0, snext[0], &kleft1, sder + kdim, &kstat);
	      if (kstat < 0)
		goto error;
	      s6diff (sder + kdim, sder, kdim, stang1);

	      s1221 (po2->c1, 0, slast[1], &kleft1, sder, &kstat);
	      if (kstat < 0)
		goto error;

	      s1221 (po2->c1, 0, snext[1], &kleft1, sder + kdim, &kstat);
	      if (kstat < 0)
		goto error;
	      s6diff (sder + kdim, sder, kdim, stang2);

	      /* Discuss directions of vectors and set up pre-topology
	         information in one direction of the curves.             */

	      if ((stang1[0] * stang2[1] - stang1[1] * stang2[0]) * (double) kdir2
		  < DZERO)
		lright[0] = SI_OUT;
	      else
		lright[0] = SI_IN;

	      if (kdir2 == 1)
		lright[1] = (lright[0] == SI_IN) ? SI_OUT : SI_IN;
	      else
		lleft[1] = (lright[0] == SI_OUT) ? SI_OUT : SI_IN;

	      /* UPDATE (ujk) : tune */
	      if (s6dist (sptpar, slast, 2) > (double) 0.05 * tref)
		{
		  /* Create help point. Set pre-topology data as SI_UNDEF. */

		  uintpt[kpos] = SISL_NULL;
		  if ((uintpt[kpos] = hp_newIntpt (2, slast, DZERO, -SI_ORD,
				     SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
					       0, 0, nullp, nullp)) == SISL_NULL)
		    goto err101;

		  kpos++;
		}
	    }
	}

      /* Perform marching in negative direction of the first curve.  */

      kdir1 = -1;
      kdir2 = -kdir2;

      /* Check if the intersection point is situated at the endpoint
	 of a curve.             */

      if (DEQUAL (sptpar[0] + tref, st1[kk1 - 1] + tref) ||
	  (kdir2 == 1 && DEQUAL (sptpar[1] + tref, st2[kn2] + tref)) ||
	  (kdir2 == -1 && DEQUAL (sptpar[1] + tref, st2[kk2 - 1] + tref)))
	{
	}
      else
	{
	  /* Perform marching.  */

	  sh1783 (po1->c1, po2->c1, aepsge, sptpar, kdir1, kdir2, slast,
		  snext, &kstat);
	  if (kstat < 0)
	    goto error;

	  if (kstat > 0)
	    {
	      /* An intersection interval is found. Set pre-topology. */

	      lleft[0] = SI_ON;
	      if (kdir2 == 1)
		lright[1] = SI_ON;
	      else
		lleft[1] = SI_ON;
	    }
	  else
	    {
	      /* Create help point. First fetch geometry information. */

	      s1221 (po1->c1, 0, slast[0], &kleft1, sder, &kstat);
	      if (kstat < 0)
		goto error;

	      s1221 (po1->c1, 0, snext[0], &kleft1, sder + kdim, &kstat);
	      if (kstat < 0)
		goto error;
	      s6diff (sder + kdim, sder, kdim, stang1);

	      s1221 (po2->c1, 0, slast[1], &kleft1, sder, &kstat);
	      if (kstat < 0)
		goto error;

	      s1221 (po2->c1, 0, snext[1], &kleft1, sder + kdim, &kstat);
	      if (kstat < 0)
		goto error;
	      s6diff (sder + kdim, sder, kdim, stang2);

	      /* Discuss directions of vectors and set up pre-topology
	         information in one direction of the curves.             */

	      if ((stang1[0] * stang2[1] - stang1[1] * stang2[0]) * (double) kdir2
		  < DZERO)
		lleft[0] = SI_OUT;
	      else
		lleft[0] = SI_IN;

	      if (kdir2 == -1)
		lleft[1] = (lleft[0] == SI_IN) ? SI_OUT : SI_IN;
	      else
		lright[1] = (lleft[0] == SI_OUT) ? SI_OUT : SI_IN;

	      /* UPDATE (ujk) : tune */
	      if (s6dist (sptpar, slast, 2) > (double) 0.05 * tref)
		{
		  /* Create help point. Set pre-topology data as SI_UNDEF. */

		  uintpt[kpos] = SISL_NULL;
		  if ((uintpt[kpos] = hp_newIntpt (2, slast, DZERO, -SI_ORD,
				     SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
					       0, 0, nullp, nullp)) == SISL_NULL)
		    goto err101;

		  kpos++;
		}
	    }
	}
    }
  else
    {
      /* The pretopology may be computed using local information. */

      if (stang1[0] * stang2[1] - stang1[1] * stang2[0] < DZERO)
	{
	  lleft[0] = SI_IN;
	  lright[0] = SI_OUT;
	  lleft[1] = SI_OUT;
	  lright[1] = SI_IN;
	}
      else
	{
	  lleft[0] = SI_OUT;
	  lright[0] = SI_IN;
	  lleft[1] = SI_IN;
	  lright[1] = SI_OUT;
	}


    }

  /* Update pre-topology of intersection point.  */
  /* UPDATE (ujk), index = -1 ?? */
  sh6settop (pintpt, -1, lleft[0], lright[0], lleft[1], lright[1], &kstat);

  /* Join intersection points, and set pretopology of help points.  */

  for (ki = 0; ki < kpos; ki++)
    {
      sh6idnpt (rintdat, &uintpt[ki], 1, &kstat);
      if (kstat < 0)
	goto error;

      if (sh6ishelp (uintpt[ki]) && uintpt[ki]->no_of_curves == 0)
	{
	  sh6settop (uintpt[ki], -1, *(pintpt->left_obj_1), *(pintpt->right_obj_1),
		     *(pintpt->left_obj_2), *(pintpt->right_obj_2), &kstat);

	  /* UPDATE (ujk) : Transfer pintpt to main point ?? */
	  /* Mark that an intersection interval is found.  */
	  sh6idcon (rintdat, &uintpt[ki], &pintpt, &kstat);
	  if (kstat < 0)
	    goto error;
	}
    }

  /* Pre-topology information computed. */

  *jnewpt = kpos;
  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

err101:*jstat = -101;
  goto out;

  /* Error in input. Conflicting dimensions.  */

err106:*jstat = -106;
  goto out;

  /* Error in input. Dimension not equal to 2. */

err108:*jstat = -108;
  goto out;

  /* Error lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
  return;
}

