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
 * $Id: sh1779.c,v 1.2 2001-03-19 15:59:05 afr Exp $
 *
 */


#define SH1779

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh1779 (SISLObject * po1, SISLObject * po2, double aepsge,
	SISLIntdat ** rintdat, SISLIntpt * pintpt,
	int *jnewpt, int *jstat)
#else
void
sh1779 (po1, po2, aepsge, rintdat, pintpt, jnewpt, jstat)
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
    As in sh1780 it is necessary to determine when a helppoint
    is close enough to a mainpoint to be neglected.
    sh1784 does not tell if it march outside the surface,
    this might be a problem. */

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology data in 3-dimensional curve-surface
*              intersection. Also find help points if necessary.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              aepsge   - Geometry resolution.
*              rintdat  - Intersection data structure.
*              pintpt   - Current intersection point.
*
*
* OUTPUT     : pintpt   - Intersection point after updating pre-topology.
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
*              sh1784     -  March along curve as long as it coincide
*                            with a surface.
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
*              Ulf J. Krystad, SI, 06.91.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int ki;			/* Counters.                               */
  int kleft1 = 0, kleft2 = 0;	/* Parameters to the evaluator.            */
  int kdim;			/* Dimension of geometry space.            */
  int kpos = 0;			/* Current position in int.pt. array.      */
  int kpar1, kpar2;		/* Index of parameter value of object.     */
  int kn;			/* Number of vertices of curve.            */
  int kk;			/* Order of curve.                         */
  int kmarch = 0;		/* Indicates if marching is necessary.     */
  int lleft[2];			/* Array storing pre-topology information. */
  int lright[2];		/* Array storing pre-topology information. */
  int *ll1, *ll2, *lr1, *lr2;	/* Pointers into pre-topology arrays.   */
  double tref;			/* Referance value in equality test.       */
  double *st;			/* Knot vector of curve.                   */
  double sder[9];		/* Result of curve evaluation.             */
  double stang[3];		/* Tangent vector of curve.                */
  double snorm[3];		/* Normal vector of surface.               */
  double slast[3];		/* Last parameter value of coincidence.    */
  double snext[3];		/* First parameter value outside interval
			           of coincidence.                         */
  double *ret_val;		/* Pointer to geo data from sh6getgeom     */
  double *ret_norm;		/* Pointer to geo data from sh6getgeom     */
  double *sptpar = pintpt->epar;/* Pointer to parameter values of int.pt.  */
  SISLCurve *qc;		/* Pointer to the curve.                   */
  SISLSurf *qs;			/* Pointer to the surface.                 */
  SISLIntpt *uintpt[2];		/* Array containing new intersection points. */
  SISLIntpt *qpt1, *qpt2;	/* Intersection points in list.            */
  double *nullp = SISL_NULL;
  double sf_low_lim[2];
  double sf_high_lim[2];

  /* Don't make pretop for help points ! */
  if (sh6ishelp (pintpt))
    {
      *jstat = 0;
      goto out;
    }

  /* Set pointers into the arrays storing pre-topology information. */

  if (po1->iobj == SISLCURVE)
    {
      qc = po1->c1;
      qs = po2->s1;

      kpar1 = 0;
      kpar2 = 1;
      ll1 = lleft;
      lr1 = lright;
      ll2 = lleft + 1;
      lr2 = lright + 1;
    }
  else
    {
      qc = po2->c1;
      qs = po1->s1;

      kpar1 = 2;
      kpar2 = 0;
      ll1 = lleft + 1;
      lr1 = lright + 1;
      ll2 = lleft;
      lr2 = lright;
    }

  /* Get pre-topology of intersection point.  */
  sh6gettop (pintpt, -1, lleft, lright, lleft + 1, lright + 1, &kstat);

  /* Describe curve partly by local parameters. */

  kdim = qc->idim;
  kk = qc->ik;
  kn = qc->in;
  st = qc->et;
  tref = st[kn] - st[kk - 1];

  sf_low_lim[0] = qs->et1[qs->ik1 - 1] + REL_COMP_RES;
  sf_low_lim[1] = qs->et2[qs->ik2 - 1] + REL_COMP_RES;
  sf_high_lim[0] = qs->et1[qs->in1] - REL_COMP_RES;
  sf_high_lim[1] = qs->et2[qs->in2] - REL_COMP_RES;

  /* Fetch geometry information, curve.  */
  sh6getgeom ((po1->iobj == SISLCURVE) ? po1 : po2,
	      (po1->iobj == SISLCURVE) ? 1 : 2,
	      pintpt, &ret_val, &ret_norm, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  /* Local copy of curve tangent */
  memcopy (stang, ret_val + kdim, kdim, DOUBLE);

  /* Fetch geometry information, surface.  */
  sh6getgeom ((po1->iobj == SISLSURFACE) ? po1 : po2,
	      (po1->iobj == SISLSURFACE) ? 1 : 2,
	      pintpt, &ret_val, &ret_norm, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  /* Local copy of surface normal */
  memcopy (snorm, ret_norm, kdim, DOUBLE);


  /* (ALA) Test if local information may be used to compute pre-topology. */
  s6length (snorm, kdim, &kstat);
  s6length (snorm, kdim, &ki);

  if (!kstat || !ki || fabs (PIHALF - s6ang (snorm, stang, kdim)) < 0.05)
    {
      /* Check if the intersection point lies at the start point of
         the curve. */

      if (DEQUAL (sptpar[kpar1] + tref, st[kn] + tref))
	;
      else
	{
	  /* Check if the intersection point is member of a list
             in this parameter direction of the curve. */

	  qpt1 = qpt2 = SISL_NULL;
	  kmarch = 1;

	  /* UPDATE (ujk) : only one list ? */
	  sh6getnhbrs (pintpt, &qpt1, &qpt2, &kstat);
	  if (kstat < 0)
	    goto error;

	  kmarch = 0;
	  if (qpt1 != SISL_NULL && qpt1->epar[kpar1] > sptpar[kpar1])
	    *lr1 = SI_ON;
	  else if (qpt2 != SISL_NULL && qpt2->epar[kpar1] > sptpar[kpar1])
	    *lr1 = SI_ON;
	  else
	    kmarch = 1;
	}

      if (kmarch)
	{
	  /* Perform marching to compute pre-topology. March first in the
             positive direction of the curve.  */

	  sh1784 (qc, qs, aepsge, sptpar, (kpar1 == 0), 1, slast, snext, &kstat);
	  if (kstat < 0)
	    goto error;

	  if (kstat == 1)
	    {
	      /* The endpoint of the curve is reached. */
	      ;
	    }
	  else if (kstat == 2)
	    ;
	  else
	    {

	      if (slast[kpar2] > sf_high_lim[0] ||
		  slast[kpar2 + 1] > sf_high_lim[1] ||
		  slast[kpar2] < sf_low_lim[0] ||
		  slast[kpar2 + 1] < sf_low_lim[1])
		;
	      else
		{


		  /* Create help point. First fetch geometry information. */

		  s1221 (qc, 0, slast[kpar1], &kleft1, sder, &kstat);
		  if (kstat < 0)
		    goto error;

		  s1221 (qc, 0, snext[kpar1], &kleft1, sder + kdim, &kstat);
		  if (kstat < 0)
		    goto error;
		  s6diff (sder + kdim, sder, kdim, stang);

		  s1421 (qs, 1, slast + kpar2, &kleft1, &kleft2, sder, snorm, &kstat);
		  if (kstat < 0)
		    goto error;

		  /* Discuss tangent- and normal vector, and set up pre-topology
	             in one direction of the curve. 		   */

		  if (s6scpr (snorm, stang, kdim) > DZERO)
		    *lr1 = SI_OUT;
		  else
		    *lr1 = SI_IN;

		  /* UPDATE (ujk) : Tuning on distance */
		  if (s6dist (sptpar, slast, 3) > (double) 0.05 * tref)
		    {
		      /* Create help point. Set pre-topology data as undefined. */
		      /* UPDATE (ujk) : If calculated values is stored, kder must
		         be 1 for curve and 2 for surface (sh6getgeom). Should
		         shevalc be used in stead of s1221 ? */

		      uintpt[kpos] = SISL_NULL;
		      if ((uintpt[kpos] = hp_newIntpt (3, slast, DZERO, -SI_ORD,
					      lleft[0], lright[0], lleft[1],
				    lright[1], 0, 0, nullp, nullp)) == SISL_NULL)
			goto err101;

		      kpos++;
		    }
		}
	    }
	}

      /* Check if the intersection point lies at the end point of
         the curve. */

      kmarch = 0;
      if (DEQUAL (sptpar[kpar1] + tref, st[kk - 1] + tref))
	;
      else
	{
	  /* Check if the intersection point is member of a list
             in this parameter direction of the curve. */

	  qpt1 = qpt2 = SISL_NULL;
	  kmarch = 1;

	  /* UPDATE (ujk) : only one list ? */
	  /* UPDATE (ujk) : only one list ? */
	  sh6getnhbrs (pintpt, &qpt1, &qpt2, &kstat);
	  if (kstat < 0)
	    goto error;

	  kmarch = 0;
	  if (qpt1 != SISL_NULL && qpt1->epar[kpar1] < sptpar[kpar1])
	    *ll1 = SI_ON;
	  else if (qpt2 != SISL_NULL && qpt2->epar[kpar1] < sptpar[kpar1])
	    *ll1 = SI_ON;
	  else
	    kmarch = 1;

	}

      if (kmarch)
	{
	  /* March in the negative direction of the curve. */

	  sh1784 (qc, qs, aepsge, sptpar, (kpar1 == 0), -1, slast, snext, &kstat);
	  if (kstat < 0)
	    goto error;

	  if (kstat == 1)
	    {
	      /* The endpoint of the curve is reached. */
	      ;
	    }
	  else if (kstat == 2)
	    ;
	  else
	    {
	      if (slast[kpar2] > sf_high_lim[0] ||
		  slast[kpar2 + 1] > sf_high_lim[1] ||
		  slast[kpar2] < sf_low_lim[0] ||
		  slast[kpar2 + 1] < sf_low_lim[1])
		;
	      else
		{

		  /* Create help point. First fetch geometry information. */

		  s1221 (qc, 0, slast[kpar1], &kleft1, sder, &kstat);
		  if (kstat < 0)
		    goto error;

		  s1221 (qc, 0, snext[kpar1], &kleft1, sder + kdim, &kstat);
		  if (kstat < 0)
		    goto error;
		  s6diff (sder + kdim, sder, kdim, stang);

		  s1421 (qs, 1, slast + kpar2, &kleft1, &kleft2, sder, snorm, &kstat);
		  if (kstat < 0)
		    goto error;

		  /* Discuss tangent- and normal vector, and set up pre-topology
	             in one direction of the curve. 		   */

		  if (s6scpr (snorm, stang, kdim) > DZERO)
		    *ll1 = SI_OUT;
		  else
		    *ll1 = SI_IN;

		  /* UPDATE (ujk) : Tuning on distance */
		  if (s6dist (sptpar, slast, 3) > (double) 0.05 * tref)
		    {
		      /* Create help point. Set pre-topology data as undefined. */
		      /* UPDATE (ujk) : If calculated values is stored, kder must
		         be 1 for curve and 2 for surface (sh6getgeom). Should
		         shevalc be used in stead of s1221 ? */

		      uintpt[kpos] = SISL_NULL;
		      if ((uintpt[kpos] = hp_newIntpt (3, slast, DZERO, -SI_ORD,
					      lleft[0], lright[0], lleft[1],
				    lright[1], 0, 0, nullp, nullp)) == SISL_NULL)
			goto err101;

		      kpos++;
		    }
		}
	    }
	}
    }
  else
    {
      /* Pre-topology data of the curve may be computed from
         local information. */

      if (s6scpr (snorm, stang, kdim) > DZERO)
	{
	  *ll1 = SI_IN;
	  *lr1 = SI_OUT;
	}
      else
	{
	  *ll1 = SI_OUT;
	  *lr1 = SI_IN;
	}

    }

  /* Update pre-topology of intersection point.  */
  /* UPDATE (ujk), index = -1 ?? */
  sh6settop (pintpt, -1, lleft[0], lright[0], lleft[1], lright[1], &kstat);

  /* Join intersection points, and set pretopology of help points.  */

  for (ki = 0; ki < kpos; ki++)
    {
      /* Help point ? */
      if (sh6ishelp (uintpt[ki]))
	sh6settop (uintpt[ki], -1, *(pintpt->left_obj_1), *(pintpt->right_obj_1),
		   *(pintpt->left_obj_2), *(pintpt->right_obj_2), &kstat);

      sh6idcon (rintdat, &uintpt[ki], &pintpt, &kstat);
      if (kstat < 0)
	goto error;
    }

  /* Pre-topology information computed. */

  *jnewpt = kpos;
  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

err101:*jstat = -101;
  goto out;


  /* Error lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
  return;
}
