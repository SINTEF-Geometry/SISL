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
 * $Id: sh6evalint.c,v 1.3 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6EVALINT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
sh6evalint (SISLObject * ob1, SISLObject * ob2, double eimpli[], int ideg,
	    SISLIntpt * pt, double aepsge,
	    double *curve_3d[], double *curve_2d_1[], double *curve_2d_2[],
	    int *jstat)
#else
void
sh6evalint (ob1, ob2, eimpli, ideg,
	    pt, aepsge,
	    curve_3d, curve_2d_1, curve_2d_2,
	    jstat)
     SISLObject *ob1;
     SISLObject *ob2;
     double eimpli[];
     int ideg;
     SISLIntpt *pt;
     double aepsge;
     double *curve_3d[];
     double *curve_2d_1[];
     double *curve_2d_2[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt, get the tangent and curvature of intersection
*
*
* INPUT      : ob1	- Pointer to geometric object(first surf).
*	       ob2	- Pointer to geometric object(second surf).
*              eimpli   - Array containing descr. of implicit surf
*	       ideg     - Type of impl surf.
*	       pt       - Pointer to the Intpt.
*              aepsge   - Absolute tolerance
*
*
* OUTPUT     : curve_3d	   - Interscetion data, object space (as in s1304, s1306)
*              curve_2d_1  - Interscetion data, param space (as in s1304, s1306)
*              curve_2d_2  - Interscetion data, param space (as in s1304)
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway. July 91.
* REVICED BY : UJK, Changed to deal with SISL_Curves also.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Sept. 1994. Avoid array bound
*              over-run in memcopy.
*********************************************************************
*/
{
  int dim;			/* Geometric dimension. */
  int kstat;
  int kpos = 1;			/* Position indicator ofr errors          */
  int left1 = 0, left2 = 0;	/* Knot navigators in s1421               */
  int kder = 2;			/* Numb of derivatives                    */
  int silhouett;		/* Flag silhouett case                    */
  int ki;			/* Variable used in loop                  */
  int ksize;			/* Size of output from s1421 or getgeom   */
  double *geom1 = SISL_NULL;		/* Output values from s1421 or getgeom    */
  double con_tang[3];		/* Constant tangent.                      */
  double *norm1 = SISL_NULL;		/* Output values from s1421 or getgeom    */
  double *geom2 = SISL_NULL;		/* Output values from s1421 or getgeom    */
  double *norm2 = SISL_NULL;		/* Output values from s1421 or getgeom    */
  double normimpl[3];		/* Normal of impl surf                    */
  double right_dir[3];		/* Right direction of 3D intersect. curve */
  double dot;			/* Scalar product */
  double dummy[6];
  double ang;
  double min_hp_ang = 0.00000000001;
  *jstat = 0;
  con_tang[0] = (double) 1.0;
  con_tang[1] = DZERO;
  con_tang[2] = DZERO;


  if (ob1->iobj != SISLSURFACE && ob1->iobj != SISLCURVE)
    goto errinp;
  if (!pt)
    goto errinp;
  if (ideg < 0)
    goto errinp;

  if (ob1->iobj == SISLSURFACE )
    dim = ob1->s1->idim;
  else
    dim = ob1->c1->idim;

  if (dim > 3 || dim < 1)
    goto errinp;

  *curve_3d = pt->geo_track_3d;
  *curve_2d_1 = pt->geo_track_2d_1;
  *curve_2d_2 = pt->geo_track_2d_2;

  if (pt->evaluated)
    goto out;

  if (ideg == 0)
    {
       /* No implicit geometry involved */
       kpos = 1;
       if (ob2->iobj != SISLSURFACE && ob2->iobj != SISLCURVE)
	 goto errinp;

       if (ob2->iobj == SISLCURVE)
	 {
	    /* At least the second object is a spline curve,
	       use this one */
	    kpos = 2;
	    if (ob2->c1->idim > 3) goto errinp;

	    /* Get geometry of first surface */
	    sh6getgeom (ob1, 1, pt, &geom1, &norm1, aepsge, &kstat);
	    if (kstat < 0)
	      goto error;

	    /* Get geometry of objects */
	    sh6getgeom (ob2, 2, pt, &geom2, &norm2, aepsge, &kstat);
	    if (kstat < 0)
	      goto error;

	    /* The number of elements to copy is given by pt->size_<obnr>
	       and we have obnr=2  (PFU 05/09-94) */
	    memcopy(*curve_3d,geom2,pt->size_2,double);

	 }
       else
	 {


	    /* Two 3d surfaces */
	    kpos = 3;
	    if (ob2->iobj != SISLSURFACE)
	      goto errinp;
	    if ((dim = ob2->s1->idim) != 3)
	      goto errinp;

	    /* Get geometry of first surface */
	    sh6getgeom (ob1, 1, pt, &geom1, &norm1, aepsge, &kstat);
	    if (kstat < 0)
	      goto error;

	    /* Get geometry of second surface */
	    sh6getgeom (ob2, 2, pt, &geom2, &norm2, aepsge, &kstat);
	    if (kstat < 0)
	      goto error;

	    /* Get normal direction */
	    s6crss (norm1, norm2, right_dir);
	    if (kstat < 0)
	      goto error;

	    /* Compute angle. */
	    ang = s6ang(norm1, norm2,3);
	    if (ang < min_hp_ang)
	    {
	       /* The point is a singular meeting point.*/
	       if (pt->iinter == SI_ORD) pt->iinter = SI_SING;
	    }

	    /* Get tangent and curvature */
	    s1304 (geom1, geom2, pt->epar, pt->epar + 2,
		   *curve_3d, *curve_2d_1, *curve_2d_2, &kstat);
	    if (kstat < 0)
	      goto error;

	    if ((dot = s6scpr (right_dir, *curve_3d + 3, 3)) < DZERO)
	      {
		 /* Change direction for tangent */
		 for (ki = 0; ki < 3; ki++)
		   (*curve_3d)[ki + 3] *= -(double) 1;
		 for (ki = 0; ki < 2; ki++)
		   {
		      (*curve_2d_1)[ki + 2] *= -(double) 1;
		      (*curve_2d_2)[ki + 2] *= -(double) 1;
		   }
	      }
	 }

       pt->evaluated = TRUE;
    }
  else
    {
       /* Implicit cases */
       if (ideg == 2000)
	 {
	    /* Here we treat the cases
	       spline surf vs implicit analytic curve
	       spline curve vs implicit analytic curve
	       spline curve vs implicit analytic surf
	       in all these cases only 3D posisition is necessary */

	    /* Clean up from 1D or 2D result */
	    if (pt->geo_data_1)
	      freearray (pt->geo_data_1);
	    if (pt->geo_data_2)
	      freearray (pt->geo_data_2);
	    pt->geo_data_1 = SISL_NULL;
	    pt->size_1 = 0;
	    pt->geo_data_2 = SISL_NULL;
	    pt->size_2 = 0;

	    /* Get the right values are computed */
	    sh6getgeom (ob1, 1, pt, &geom1, &norm1, aepsge, &kstat);
	    if (kstat < 0)
	      goto error;

     	    memcopy(*curve_3d,geom1,dim,double);
     	    memcopy((*curve_3d)+dim,con_tang,dim,double);

	 }
       else
	 {
	    if (ideg == 1003 || ideg == 1004 || ideg == 1005)
	      {
		 /* Silhouette cases, B-spline surface */
		 kpos = 3;
		 ksize = 33;
		 silhouett = TRUE;
		 kder = 3;

	      }
	    else
	      {
		 /* Analytic surf vs B-spline surface */

		 kpos = 4;
		 ksize = 21;
		 silhouett = FALSE;
		 kder = 2;
	      }

	    if (pt->size_1 != ksize)
	      {
		 /* Clean up from 1D result */
		 if (pt->geo_data_1)
		   freearray (pt->geo_data_1);
		 if (pt->geo_data_2)
		   freearray (pt->geo_data_2);
		 pt->geo_data_1 = SISL_NULL;
		 pt->size_1 = 0;
		 pt->geo_data_2 = SISL_NULL;
		 pt->size_2 = 0;


		 if ((pt->geo_data_1 = newarray (ksize, DOUBLE))
		     == SISL_NULL)
		   goto err101;
		 pt->size_1 = ksize;
		 geom1 = pt->geo_data_1;
		 norm1 = pt->geo_data_1 + ksize - 3;

		 s1422 (ob1->s1, kder, pt->iside_1, pt->iside_2,
			pt->epar, &left1, &left2, geom1,
			norm1, &kstat);
		 if (kstat < 0)
		   goto error;
	      }
	    else
	      {
		 /* The right values are computed */
		 sh6getgeom (ob1, 1, pt, &geom1, &norm1, aepsge, &kstat);
		 if (kstat < 0)
		   goto error;

	      }


	    /* Get normal of implicit surface */
	    s1331 (geom1, eimpli, ideg, kder = -1, dummy, normimpl, &kstat);
	    if (kstat < 0)
	      goto error;

	    /* Get the right direction of the intersection curve */
	    if (silhouett)
	      {
		 ang = 1.5; /* Not used */
		 memcopy (right_dir, normimpl, 3, DOUBLE);
		 for (ki=0;ki<3;ki++) right_dir[ki] *= -(double)1.0;
	      }
	    else
	    {
	       /* Compute angle. */
	       ang = s6ang(norm1, normimpl,3);
	       s6crss (norm1, normimpl, right_dir);
	    }

	    /* Get tangent and curvature to the real intersection. */
	    s1306 (geom1, pt->epar,
		   eimpli, ideg, *curve_3d, *curve_2d_1, &kstat);
	    if (kstat < 0)
	      goto error;
	    if (kstat == 2)
	    {
	       /* The point is a singular meeting point.*/
	       if (pt->iinter == SI_ORD) pt->iinter = SI_SING;
	    }
	    else if (kstat == 10)
	    {
	       /* The point is a singular non-meeting point.
		  Tangent found, but sign might be wrong. */
	       if (pt->iinter == SI_ORD || pt->iinter == SI_SING )
		  pt->iinter = SI_TOUCH;
	    }
	    else if (ang < min_hp_ang)
	    {
	       /* The point is a singular meeting point.*/
	       if (pt->iinter == SI_ORD) pt->iinter = SI_SING;
	    }
	    else
	    if ((dot = s6scpr (right_dir, *curve_3d + 3, 3)) < DZERO)
	      {
		 /* Change direction for tangent */
		 for (ki = 0; ki < 3; ki++)
		   (*curve_3d)[ki + 3] *= -(double) 1;
		 for (ki = 0; ki < 2; ki++)
		   (*curve_2d_1)[ki + 2] *= -(double) 1;

	      }
	 }


       pt->evaluated = TRUE;

    }

  *jstat = 0;
  goto out;

  /* ---------- ERROR EXITS --------------------------- */
  /* Error in alloc  */
  err101:
     *jstat = -101;
  s6err ("shevalint", *jstat, kpos);
  goto out;

  /* Error in lower level */
  error:
     *jstat = kstat;
  s6err ("shevalint", *jstat, kpos);
  goto out;

  /* Error in input */
  errinp:
     *jstat = -200;
  s6err ("shevalint", *jstat, kpos);
  goto out;


  out:;
}
