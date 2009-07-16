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
 * $Id: sh1787.c,v 1.4 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH1787

#include "sislP.h"



#if defined(SISLNEEDPROTOTYPES)
void
    sh1787 (SISLObject * po1, SISLObject * po2, double aepsge,
	SISLIntdat ** rintdat, SISLIntpt * pintpt, int *jnewpt,
	int *jstat)
#else
void
   sh1787 (po1, po2, aepsge, rintdat, pintpt, jnewpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **rintdat;
     SISLIntpt *pintpt;
     int *jnewpt;
     int *jstat;
#endif
/*
******************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology data in n-dimensional surface-point
*              intersection.
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
* CALLS      : shevalc  -  Evaluate surface using filtered coefficients.
*              s6ang    -  Angle between two vectors.
*              s6idcon  -  Connect two intersection points.
*              hp_newIntpt -  Create new intersection point.
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. 09.91.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int kdim;			/* Dimension of geometry space.              */
  int kn1;			/* Nmb vertices of surface in 1st direc.     */
  int kn2;			/* Nmb vertices of surface in 1st direc.     */
  int kk1;			/* Order of surface in 1st direction.        */
  int kk2;			/* Order of surface in 1st direction.        */
  int kpos = 0;			/* Current position in int.pt. array.        */
  int lleft[2];			/* Array storing pre-topology information.   */
  int lright[2];		/* Array storing pre-topology information.   */
  int *ll1, *ll2, *lr1, *lr2;	/* Pointers into pre-topology arrays.        */
  double tpoint[3];		/* Value of point to intersect.              */
  double sder[21];		/* Result of surface evaluation.             */
  double *st1;			/* First knot vector of surface.             */
  double *st2;			/* Second knot vector of surface.            */
  double tref1;			/* Referance value in equality test.         */
  double tref2;			/* Referance value in equality test.         */
  SISLSurf *qs;		        /* Pointer to current surface.               */
  double *ret_val;		/* Pointer to geo data from sh6getgeom       */
  double *ret_norm;		/* Pointer to geo data from sh6getgeom       */
  int i;                        /* Loop variable.                            */
  double cross;                 /* utang x vtang.                            */
  double in_out[2];             /* To be used in touchy situations           */
  /* ----------------------------------------------------------------------  */

  /* Don't make pretop for help points ! */
  /* Oh, yes ?, 2D is some nice case ! */
  /* if (sh6ishelp (pintpt))
     {
     *jstat = 0;
     goto out;
     }
     */

  /* Set pointers into the arrays storing pre-topology information. */
  if (po1->iobj == SISLSURFACE)
    {
      ll1 = lleft;
      lr1 = lright;
      ll2 = lleft + 1;
      lr2 = lright + 1;
    }
  else
    {
      ll1 = lleft + 1;
      lr1 = lright + 1;
      ll2 = lleft;
      lr2 = lright;
    }

  /* Get pre-topology information. */
  sh6gettop (pintpt, -1, lleft, lright, lleft + 1, lright + 1, &kstat);
  if (kstat < 0)
    goto error;

  /* Test dimension of geometry space. */
  if (po1->iobj == SISLSURFACE)
      qs = po1->s1;
  else
      qs = po2->s1;

  kdim = qs->idim;
  if (kdim != 2)
    goto err106;

  /* Store surface information in local parameters. */

  kn1 = qs->in1;
  kn2 = qs->in2;
  kk1 = qs->ik1;
  kk2 = qs->ik2;
  st1 = qs->et1;
  st2 = qs->et2;
  tref1 = st1[kn1] - st1[kk1 - 1];
  tref2 = st2[kn2] - st2[kk2 - 1];

  /* Fetch geometry information, point.  */
  sh6getgeom ((po1->iobj == SISLPOINT) ? po1 : po2,
	      (po1->iobj == SISLPOINT) ? 1 : 2,
	      pintpt, &ret_val, &ret_norm, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  for(i=0; i<kdim; i++)
      tpoint[i]=ret_val[i];


  /* Fetch geometry information, surface.  */
  sh6getgeom ((po1->iobj == SISLSURFACE) ? po1 : po2,
	      (po1->iobj == SISLSURFACE) ? 1 : 2,
	      pintpt, &ret_val, &ret_norm, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  for(i=0; i<kdim*3; i++)
       sder[i]=ret_val[i];

/* Set normal vector from the 2D tangent vectors. */

  cross = sder[kdim]*sder[2*kdim+1] + sder[kdim+1]*sder[2*kdim];

  /*  Could improve this test. */
  if (fabs(cross) > ANGULAR_TOLERANCE)
    {
      /* Compute pre-topology using local information.  */

      if (cross > 0)
	{
	  *ll1 = SI_UNDEF;
	  *lr1 = SI_UNDEF;
	  *ll2 = SI_IN;
	  *lr2 = SI_OUT;
	}
      else
	{
	  *ll1 = SI_UNDEF;
	  *lr1 = SI_UNDEF;
	  *ll2 = SI_OUT;
	  *lr2 = SI_IN;
	}

    }
  else if (qs->pdir && qs->pdir->ecoef &&
	   (DNEQUAL(qs->pdir->ecoef[0],DZERO) ||
	    DNEQUAL(qs->pdir->ecoef[1],DZERO)))
    {
       /* March to find help points.
	   Not implemented yet. */
       /* UJK, I'm not sure, but something like this should work :
	  Remeber we are in a simple case situation !
	  */
       in_out[0] =  (double)1.0;
       in_out[1] = -(double)1.0;

       if (s6scpr(qs->pdir->ecoef,in_out,kdim) > 0)
	  {
	     *ll1 = SI_UNDEF;
	     *lr1 = SI_UNDEF;
	     if (*ll2 == SI_UNDEF &&
		 *lr2 == SI_UNDEF)
	     {
		*ll2 = SI_IN;
		*lr2 = SI_OUT;
	     }
	     else if (!((*ll2 == SI_IN && *lr2 == SI_IN) ||
		      (*ll2 == SI_OUT && *lr2 == SI_OUT)))
		{
		   if (*ll2 != SI_IN) *ll2 = SI_IN;
		}
	  }
	  else
	  {
	     *ll1 = SI_UNDEF;
	     *lr1 = SI_UNDEF;
	     if (*ll2 == SI_UNDEF &&
		 *lr2 == SI_UNDEF)
	     {
		*ll2 = SI_OUT;
		*lr2 = SI_IN;
	     }
	     else if (!((*ll2 == SI_IN && *lr2 == SI_IN) ||
		      (*ll2 == SI_OUT && *lr2 == SI_OUT)))
		{
		   if (*lr2 != SI_IN) *lr2 = SI_IN;
		}
	  }

    }

  /* Update pretopology of intersection point.  */

  sh6settop (pintpt, -1, lleft[0], lright[0], lleft[1], lright[1], &kstat);
  if (kstat < 0)
    goto error;

  /* Pre-topology information computed. */

  *jnewpt = kpos;
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
