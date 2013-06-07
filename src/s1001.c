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
 * $Id: s1001.c,v 1.7 2001-03-19 15:58:40 afr Exp $
 *
 */


#define S1001

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
  s1001 (SISLSurf * ps, double min1, double min2,
		double max1, double max2,
		SISLSurf ** rsnew, int *jstat)
#else
void
s1001 (ps, min1, min2, max1, max2, rsnew, jstat)
     SISLSurf *ps;
     double min1;
     double min2;
     double max1;
     double max2;
     SISLSurf **rsnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : To pick a part of a surface.
*              The surface produced will always be k-regular, i.e. with
*              k-tupple end knots.
*
*
*
* INPUT      : ps	- Surface to pick a part of.
*	       min1     - Min value 1. parameter direction.
*	       min2     - Min value 2. parameter direction.
*	       max1     - Max value 1. parameter direction.
*	       max2     - Max value 1. parameter direction.
*
*
*
* OUTPUT     : rsnew	- The new, refined surface.
*              jstat	- status messages
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
*
* WRITTEN BY : Ulf J. Krystad, SI, 04.92.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 94-08. Added error propagation.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 94-11. Looked at handling of
*              periodics - since this routine is used to convert a periodic into
*              a k-regular surface, it must always return a k-regular surface,
*              even if the whole periodic parameter range (ik-1 to in) is picked.
*              To do this correctly the 'cuopen' flags must indicate closed after
*              the picking of the whole periodic range.
*
**********************************************************************/
{
  int kstat;			/* Local status variable.	       */
  int kpos = 0;			/* Position of error.		       */
  int kleft1 = 0;		/* Knot navigator.		       */
  int kleft2 = 0;		/* Knot navigator.		       */
  int kleft3 = 0;		/* Knot navigator.		       */
  int kleft4 = 0;		/* Knot navigator.		       */
  int kdim = ps->idim;		/* Dimension of geometry space.        */
  int kkind = ps->ikind;	/* Kind of surface.                    */
  int kn1;			/* Number of vertices in 1. par. dir.  */
  int kn2;			/* Number of vertices in 2. par. dir.  */
  int cuopen_1, cuopen_2;	/* Open flags for the new surface.     */
  int change_1,change_2;	/* Flag, need to change surf in dir ?  */
  int wholeperi1 = FALSE;       /* Flag, pick whole peri. param. range */
  int wholeperi2 = FALSE;       /* Flag, pick whole peri. param. range */
  double *st1=SISL_NULL;		/* Knot vector in 1. par. dir.         */
  double *st2=SISL_NULL;		/* Knot vector in 2. par. dir.         */
  double *scoef1 = SISL_NULL;	/* Coefficients of input curve to
			           refinement in 1. par. dir.          */
  double *scoef2 = SISL_NULL;	/* Coefficients of refined surface.    */
  double *scoef  = SISL_NULL;	/* Coefficients of refined surface.    */
  SISLCurve *qc1 = SISL_NULL;	/* Input curve to pick curve.          */
  SISLCurve *qc2 = SISL_NULL;	/* Output curve from pick curve.       */
  SISLCurve *qc3 = SISL_NULL;	/* Output curve from pick curve.       */
  double *oldcoef;           	/* Pointer to vertices of old surf.    */
  /* ----------------------------------------------------------------- */

  if(kkind == 2 || kkind == 4)
  {
     oldcoef = ps->rcoef;
     kdim++;
  }
  else
  {
     oldcoef = ps->ecoef;
  }

  kleft1=ps->ik1-1;
  kleft2=ps->in1;
  kleft3=ps->ik2-1;
  kleft4=ps->in2;
  change_1 = change_2 = TRUE;

  if ( min1 == ps->et1[ps->ik1 -1]  &&  max1 == ps->et1[ps->in1] )
  {
    if ( s6knotmult(ps->et1,ps->ik1,ps->in1,
		    &kleft1,ps->et1[ps->ik1-1],&kstat) == ps->ik1 &&
	 s6knotmult(ps->et1,ps->ik1,ps->in1,
		    &kleft2,ps->et1[ps->in1],&kstat) == ps->ik1 )
      change_1 = FALSE;
    else
      wholeperi1 = ( ps->cuopen_1 == SISL_SURF_PERIODIC );
  }

  if ( min2 == ps->et2[ps->ik2 -1]  &&  max2 == ps->et2[ps->in2] )
  {
    if ( s6knotmult(ps->et2,ps->ik2,ps->in2,
		    &kleft3,ps->et2[ps->ik2-1],&kstat) == ps->ik2 &&
	 s6knotmult(ps->et2,ps->ik2,ps->in2,
		    &kleft4,ps->et2[ps->in2],&kstat) == ps->ik2 )
      change_2 = FALSE;
    else
      wholeperi2 = ( ps->cuopen_2 == SISL_SURF_PERIODIC );
  }

  if (change_1)
    {
       /* Treat the first parameter direction of the
	  surface. First express the surface as a curve.  */
       if ((scoef1 = newarray (kdim * ps->in1 * ps->in2, double)) == SISL_NULL)
	 goto err101;

       /* Change parameter directions of surface.  */
       s6chpar (oldcoef, ps->in1, ps->in2, kdim, scoef1);

       /* Create curve.  */
       qc1 = newCurve (ps->in1, ps->ik1, ps->et1, scoef1, 1, kdim * ps->in2, 0);
       if (qc1 == SISL_NULL)
	 goto err101;
       qc1->cuopen = ps->cuopen_1;

       /* Pick part of curve */
       s1713 (qc1, min1, max1, &qc2, &kstat);
       if (kstat < 0)
	 goto error;

       /* Change parameter directions of the coefficient array of
	  the refined curve.     */

       if ((scoef2 = newarray (qc2->in *ps->in2 * kdim, DOUBLE)) == SISL_NULL)
	 goto err101;
       s6chpar (qc2->ecoef, ps->in2, qc2->in, kdim, scoef2);

       /* Set local parameters of refined surface. */

       kn1 = qc2->in;
       kn2 = ps->in2;
       st1 = qc2->et;
       st2 = ps->et2;
       if ( wholeperi1 )
	 cuopen_1 = SISL_SURF_CLOSED;
       else
	 cuopen_1 = qc2->cuopen;

       /* Free curve used as input to s1713. */
       if (qc1)
	 freeCurve (qc1);
       qc1 = SISL_NULL;
    }

  else
    {
       /* Set local parameters of input surface. */

       kn1 = ps -> in1;
       kn2 = ps -> in2;
       st1 = ps -> et1;
       st2 = ps -> et2;
       scoef2   = oldcoef;
       cuopen_1 = ps->cuopen_1;
    }

  if (change_2)
    {
       /* Treat the first parameter direction of the
	  surface. First express the surface as a curve.  */

       if ((qc1 = newCurve (kn2, ps->ik2, st2, scoef2, 1, kn1 * kdim, 0))
	   == SISL_NULL)
	 goto err101;
       qc1->cuopen = ps->cuopen_2;

       /* Pick part of curve */
       s1713 (qc1, min2, max2, &qc3, &kstat);
       if (kstat < 0)
	 goto error;


       /*	Set local parameters of the refined surface. */
       kn2 = qc3->in;
       st2 = qc3->et;
       scoef = qc3->ecoef;
       if ( wholeperi2 )
	 cuopen_2 = SISL_SURF_CLOSED;
       else
	 cuopen_2 = qc3->cuopen;

       /* Free curve used as input to s1713. */
       if (qc1)
	 freeCurve (qc1);
       qc1 = SISL_NULL;
    }
  else
    {
       scoef = scoef2;
       cuopen_2 = ps->cuopen_2;
    }

  /* Express result as a surface.  */
  if ((*rsnew = newSurf (kn1, kn2, ps->ik1, ps->ik2, st1, st2,
			 scoef, kkind, ps->idim, 1)) == SISL_NULL)
    goto err101;


  (*rsnew)->cuopen_1 = cuopen_1;
  (*rsnew)->cuopen_2 = cuopen_2;

  /* Task done  */

  *jstat = 0;
  goto out;

  /* ---------------------- ERROR EXITS ------------------------------- */
  /* Error in scratch allocation.  */

err101:
  *jstat = -101;
  s6err ("s1001", *jstat, kpos);
  goto out;

  /* Error in lower level routine.  */

error:
  *jstat = kstat;
  s6err ("s1001", *jstat, kpos);
  goto out;

  out:
     /* Free scratch occupied by local arrays and objects.  */

     if (change_1)
       {
	  if (scoef1) freearray (scoef1);
	  if (scoef2) freearray (scoef2);
	  scoef1 = SISL_NULL;
	  scoef2 = SISL_NULL;
       }

     if (qc1) freeCurve (qc1);
     if (qc2) freeCurve (qc2);
     if (qc3) freeCurve (qc3);
}
