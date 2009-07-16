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
 * $Id: sh6idget.c,v 1.2 2001-03-19 15:59:08 afr Exp $
 *
 */


#define SH6IDGET

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6idget (SISLObject * po1, SISLObject * po2, int ipar, double apar, SISLIntdat * pintdat,
	  SISLIntdat ** rintdat, double aepsge, int *jstat)
#else
void
sh6idget (po1, po2, ipar, apar, pintdat, rintdat, aepsge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     int ipar;
     double apar;
     SISLIntdat *pintdat;
     SISLIntdat **rintdat;
     double aepsge;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To pick out points from one intdat(pintdat) and copy them over
*              to another intdat(rintdat) which has one parameter
*              dimension lower. The parameter index to be removed is
*              ipar. Only those points in pintdat which has the value
*              apar (ie epar[ipar]==apar) are copied.
*              Help points are transformed to main points and
*              the geometry data of the object not reduced is copied.
*              Pretopology for the first connection is copied.
*
*
*
*
* INPUT      : po1     - First object in the intersection.
*              po2     - Second object in the intersection.
*              ipar    - Number of the parameter that is missing in rintdat.
*              apar    - Parameter value of the missing parameter.
*              pintdat - Pointer to intersection data on the mother problem.
*              aepsge  - Geometry tolerance
*
* OUTPUT     : rintdat  - Pointer to a pointer to intersection data.
*              jstat    - status messages
*                               = 0      : Inserting done.
*                               < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6err      - Gives error message.
*              sh6idnpt   - Insert a new intpt structure.
*
* WRITTEN BY : UJK, 06.91.
*
*********************************************************************
*/
{
  int kstat;			/* Local status variable.                 */
  int kpos = 0;			/* Position of error.                     */
  int ki, kj, kn, kl;
  int keep_first;		/* Flag, which object is not reduced      */
  double tstart[4];
  double tend[4];
  double spar[4];		/* Storing uppdated parametervalues.      */
  double tlow, thigh;
  double help_arr[4];
  double thelp;
  SISLIntpt *qpt = SISL_NULL, *pinter=SISL_NULL;
  double *nullp = SISL_NULL;
  int found = FALSE;
  int ind_div, ind_other;
  SISLObject *qo_div = SISL_NULL, *qo_other = SISL_NULL;
  int kleftt = 0, klefts = 0;
  double point[3];
  int log_ind;

  /* Find out which object the parameter belongs to */
 if (ipar < po1->iobj)
 {
   if(ipar == 1)	log_ind = 0;
   else 		log_ind=1;
   qo_div = po1;
   qo_other = po2;
   ind_div = 0;
   ind_other = po1->iobj;
      keep_first = 0;
 }
 else
 {
   if(ipar == po1->iobj)	log_ind = ipar +1;
   else 			log_ind=ipar-1;
   qo_div = po2;
   qo_other = po1;
   ind_div = po1->iobj;
   ind_other = 0;
   keep_first = 1;
    }

  if (pintdat == SISL_NULL)
    goto out;

  /* ----------------------------------------- */

  for (ki = 0; ki < pintdat->ipoint; ki++)
  {
    sh6isinside (po1, po2, pintdat->vpoint[ki], &kstat);
    if (kstat < 0)
      goto error;

    if (kstat)
    {
      for (kj = 0; kj < (pintdat->vpoint[ki])->no_of_curves;kj++)
      {
	qpt = sh6getnext (pintdat->vpoint[ki], kj);
	sh6isinside (po1, po2, qpt, &kstat);
	if (kstat < 0)
	  goto error;

	/* For surface, check on curve_dir */
	if (kstat &&
	    (qo_div->iobj == SISLCURVE ||
	     (qo_div->iobj == SISLSURFACE &&
	      (pintdat->vpoint[ki]->curve_dir[kj] & (1 << (log_ind + 1))))))
	{
	  /* curve: */
	  tlow = pintdat->vpoint[ki]->epar[ipar];
	  thigh = qpt->epar[ipar];
	  if (thigh < tlow)
	  {
	    thelp = thigh;
	    thigh = tlow;
	    tlow = thelp;
	  }

	  if (apar > tlow && apar < thigh)
	  {
	    found = TRUE;
	    break;
	  }
	}
      }
    }
    if (found)
    {

      for (kl=0; kl<qpt->ipar; kl++)
	help_arr[kl] = (double)0.5*(pintdat->vpoint[ki]->epar[kl]
				    + qpt->epar[kl]);
      help_arr[ipar] = apar;

      /* Prepare for iteration. */
      if (qo_div->iobj == SISLCURVE)
      {
	kleftt=0;
	s1221 (qo_div->c1, 0, apar, &kleftt, point, &kstat);
	if (kstat < 0)
	  goto error;
      }
      else
      {
	kleftt=0;
	klefts=0;
	s1421 (qo_div->s1, 0, help_arr + ind_div,
	       &kleftt, &klefts, point, nullp, &kstat);
	if (kstat < 0)
	  goto error;
      }


      sh6ptobj (point, qo_other, aepsge, help_arr + ind_other,
		    help_arr + ind_other, &kstat);
      if (kstat == 1)
      {
	/* Point found, insert */

	pinter = hp_newIntpt (qpt->ipar, help_arr, DZERO, qpt->iinter,
			      SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
			      0, 0, nullp, nullp);
	if (pinter == SISL_NULL)
	  goto err101;

	sh6insert (&pintdat, pintdat->vpoint[ki], qpt, &pinter, &kstat);
	if (kstat < 0)
	  goto error;
      }
    }
    if (qo_div->iobj == SISLCURVE)
      break;
    else
      found = FALSE;
  }



  /* ----------------------------------------- */



  if (po1->iobj == SISLCURVE)
  {
    tstart[0] = po1->c1->et[po1->c1->ik - 1];
    tend[0] = po1->c1->et[po1->c1->in];
  }
  else if (po1->iobj == SISLSURFACE)
  {
    tstart[0] = po1->s1->et1[po1->s1->ik1 - 1];
    tend[0] = po1->s1->et1[po1->s1->in1];
    tstart[1] = po1->s1->et2[po1->s1->ik2 - 1];
    tend[1] = po1->s1->et2[po1->s1->in2];
  }

  if (po2->iobj == SISLCURVE)
  {
    tstart[po1->iobj] = po2->c1->et[po2->c1->ik - 1];
    tend[po1->iobj] = po2->c1->et[po2->c1->in];
  }
  else if (po2->iobj == SISLSURFACE)
  {
    tstart[po1->iobj] = po2->s1->et1[po2->s1->ik1 - 1];
    tend[po1->iobj] = po2->s1->et1[po2->s1->in1];
    tstart[po1->iobj + 1] = po2->s1->et2[po2->s1->ik2 - 1];
    tend[po1->iobj + 1] = po2->s1->et2[po2->s1->in2];
  }

  /* Fix pick values for reduced paramater. */
  tstart[ipar] = tend[ipar] = apar;


  /* Uppdate the array. */

  for (ki = 0; ki < pintdat->ipoint; ki++)
  {
    for (kj = 0; kj < pintdat->vpoint[ki]->ipar; kj++)
      if ((DNEQUAL (pintdat->vpoint[ki]->epar[kj], tstart[kj]) &&
	   pintdat->vpoint[ki]->epar[kj] < tstart[kj]) ||
	  (DNEQUAL (pintdat->vpoint[ki]->epar[kj], tend[kj]) &&
	   pintdat->vpoint[ki]->epar[kj] > tend[kj]))
	break;

    if (kj == pintdat->vpoint[ki]->ipar)
    {
      for (kn = 0; kn < ipar; kn++)
	spar[kn] = pintdat->vpoint[ki]->epar[kn];
      for (; kn < pintdat->vpoint[ki]->ipar - 1; kn++)
	spar[kn] = pintdat->vpoint[ki]->epar[kn + 1];

      /* Point accepted, insert into rintdat. */
      /* VSK. Let the point be a normal main point.  */

      qpt = hp_newIntpt (pintdat->vpoint[ki]->ipar - 1, spar,
			 pintdat->vpoint[ki]->adist,1,
			 pintdat->vpoint[ki]->left_obj_1[0],
			 pintdat->vpoint[ki]->right_obj_1[0],
			 pintdat->vpoint[ki]->left_obj_2[0],
			 pintdat->vpoint[ki]->right_obj_2[0],
			 (keep_first ? pintdat->vpoint[ki]->size_1 : 0),
			 (keep_first ? 0 : pintdat->vpoint[ki]->size_2),
			 (keep_first ? pintdat->vpoint[ki]->geo_data_1 : nullp),
			 (keep_first ? nullp : pintdat->vpoint[ki]->geo_data_2));

      if (qpt == SISL_NULL)
	goto err101;

      sh6idnpt (rintdat, &qpt, 1, &kstat);
      if (kstat < 0)
	goto error;
    }
  }

  *jstat = 0;
  goto out;


/* Error in space allocation.  */

err101:*jstat = -101;
  s6err ("sh6idget", *jstat, kpos);
  goto out;

/* Error in sub function.  */

error:*jstat = kstat;
  s6err ("sh6idget", *jstat, kpos);
  goto out;

out:;
}
