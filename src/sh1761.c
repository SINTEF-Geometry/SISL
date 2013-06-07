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
 * $Id: sh1761.c,v 1.7 2001-03-19 15:59:04 afr Exp $
 *
 */


#define SH1761

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh1761 (SISLObject * po1, SISLObject * po2, double aepsge, SISLIntdat ** pintdat, int *jstat)
#else
void
sh1761 (po1, po2, aepsge, pintdat, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **pintdat;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*          NOTE : Comments for further developments/tasks starts
*                 with <comment-sign> UPDATE :
*
*
* PURPOSE    : Find all intersections between two objects (of type
*              point, curve or surface). In this rouine the outer
*              edges/endpoints of the objects are treated.
*
*
*
* INPUT      : po1    - First object in the intersection.
*              po2    - Second object in the intersection.
*              aepsge - Geometry resolution.
*              *jstat    - Flag
*                          = 202 : Complicated point-surface intersection
*                                  in 3D. Perform extra interception test.
*
*
*
*
*
*
* OUTPUT     : intdat - Intersection dates found beetween the two
*                       objects.
*              jstat  - status messages
*                                       = 1      : Intersection found
*                                       = 0      : no intersection
*                                       < 0      : error
*
*
* METHOD     : This function is computing point/point intersection
*              otherwise it is computing edge/endpoint intersections
*              by recurson on one end object and the other object,
*              and futher calling a rutine for computing intersections
*              in the inner of an object.
*
*
* REFERENCES :
*
*-
* CALLS      : s6err      - Gives error message.
*              s6dist     - Compute the distance beetween two point.
*              s1435      - Pick edge curve from a surface.
*              s1438      - Pick endpoint from a curve.
*              sh1762      - Find the intersections in the inner of the objects.
*              sh1790      - Perform BOX test.
*              sh6idnpt    - Put a new intpt to given intdat.
*              sh6idput    - Put contence of one intdat in an other intdat.
*              sh6idlis    - Compute list elements from given intdat.
*              sh6idalledg - Uppdate edges from given intdat.
*              freeEdge   - Free space occupied by given edge.
*              freePoint  - Free space occupied by given point.
*              freeCurve  - Free space occupied by given curve.
*              freeObject - Free space occupied by given object.
*              freeIntdat - Free space occupied by given intdat.
*              newEdge    - Create new edge-structure.
*              hp_newIntpt   - Create new intersection point-structure.
*              newObject  - Create new object.
*
* WRITTEN BY : Arne Laksaa, 05.89.
*              UJK, 06.91 newi
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov.1994. Fixed memory
*              leaks from 'qt', 'qedge', and 'qintdat' - occured on error exits.
*********************************************************************
*/
{
  int kstat = 0;		/* Local status variable.                 */
  int kpos = 0;			/* Position of error.                     */
  int ktotal = 1;		/* Make totally expanded box.             */
  int kxintercept = (*jstat == 202);  /* Extra interception               */
  double tpar;			/* Help variable used for parameter value
				   and geometric distance.                */
  SISLObject *po1_kreg=SISL_NULL;    /* Pointer to first object converted to
				   k-regular basis.                       */
  SISLObject *po2_kreg=SISL_NULL;    /* Pointer to second object converted to
				   k-regular basis.                       */
  SISLIntpt *qt = SISL_NULL;         /* Temporary intersection point. */
  SISLEdge *qedge[2];	        /* Edges for use in s1862().      */
  SISLIntdat *qintdat = SISL_NULL;	/* Intdat for use in recurson.    */

  double *nullp = SISL_NULL;
  int idummy;

  qedge[0] = qedge[1] = SISL_NULL;  /* PFU - to fix memory leak */


  /* Ensure K-regularity on B-spline basis ________________*/
  if (po1->iobj == SISLCURVE)
    {
       if (po1->c1->cuopen == SISL_CRV_PERIODIC)
	 {
	    if ((po1_kreg = newObject (SISLCURVE)) == SISL_NULL)
	      goto err101;
	    make_cv_kreg(po1->c1, &po1_kreg->c1, &kstat);
	    if (kstat < 0) goto error;
	 }
       else po1_kreg = po1;

    }
  else if (po1->iobj == SISLSURFACE)
    {
       if (po1->s1->cuopen_1 == SISL_CRV_PERIODIC ||
	   po1->s1->cuopen_2 == SISL_CRV_PERIODIC)
	 {
	    if ((po1_kreg = newObject (SISLSURFACE)) == SISL_NULL)
	      goto err101;
	    make_sf_kreg(po1->s1, &po1_kreg->s1, &kstat);
	    if (kstat < 0) goto error;
	 }
       else po1_kreg = po1;
    }
  else po1_kreg = po1;



  if (po2->iobj == SISLCURVE)
    {
       if (po2->c1->cuopen == SISL_CRV_PERIODIC)
	 {
	    if ((po2_kreg = newObject (SISLCURVE)) == SISL_NULL)
	      goto err101;
	    make_cv_kreg(po2->c1, &po2_kreg->c1, &kstat);
	    if (kstat < 0) goto error;
	 }
       else po2_kreg = po2;
    }
  else if (po2->iobj == SISLSURFACE)
    {
       if (po2->s1->cuopen_1 == SISL_CRV_PERIODIC ||
	   po2->s1->cuopen_2 == SISL_CRV_PERIODIC)
	 {
	    if ((po2_kreg = newObject (SISLSURFACE)) == SISL_NULL)
	      goto err101;
	    make_sf_kreg(po2->s1, &po2_kreg->s1, &kstat);
	    if (kstat < 0) goto error;
	 }
       else po2_kreg = po2;
    }

  else po2_kreg = po2;

  /* End of ensure K-regularity on B-spline basis ________________*/


  if (po1_kreg->iobj == SISLPOINT && po2_kreg->iobj == SISLPOINT)
    {
      /* Control the dimension of the two points. */

      if (po1_kreg->p1->idim != po2_kreg->p1->idim)
	goto err106;

      /* Computing the distanse beetween the points. */

      tpar = s6dist (po1_kreg->p1->ecoef, po2_kreg->p1->ecoef, po1_kreg->p1->idim);

      if (tpar <= aepsge)
	{
	  /* PFU, moved to top to fix memory leak: SISLIntpt *qt = SISL_NULL; */

	  *jstat = 1;		/* Mark intersection found. */
	  /* UJK , newi */
	  /* Making intersection point. */

	  qt = hp_newIntpt (0, &tpar, DZERO, SI_ORD,
			    SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
			    0, 0, nullp, nullp);
	  if (qt == SISL_NULL)
	    goto err101;

	  /* Uppdating pintdat. */

	  sh6idnpt (pintdat, &qt, 1, &kstat);
	  if (kstat < 0)
	    goto error;

	  qt = SISL_NULL;      /* PFU - to fix memory leak */
	}
      else
	*jstat = 0;
    }
  else
    {
      /* Test if intersection is possible (perform box-test).  */

      sh1790 (po1_kreg, po2_kreg, ktotal, aepsge, &kstat);
      if (kstat < 0)
	goto error;

      /* We may have four different values on kstat.
	 kstat = 0 : Intersection not possible.
	 kstat = 1 : The two boxes overlap.
	 kstat = 2 : The two "bezier" boxes is just touching.
	 kstat = 3 : The two boxes is both inside a microbox of aepsge. */


      if (kstat)
	{
	  int ki, kj;		/* Counters.                      */
	  int kedg = 0;		/* Number of parameter direction. */
	  SISLObject *qo1 = SISL_NULL;  /* Help pointer.                  */
	  SISLObject *qo2 = SISL_NULL;  /* Help pointer.                  */

	  /* PFU, moved to top to fix memory leak:
	   *  SISLEdge *qedge[2];
	   *  SISLIntdat *qintdat = SISL_NULL;
	   */

	  qintdat = SISL_NULL; /* PFU to fix memory leak. */
	  qedge[0] = qedge[1] = SISL_NULL;
	  *jstat = 0;

	  for (kj = 0, qo1 = po1_kreg, qo2 = po2_kreg; kj < 2; kj++, qo1 = po2_kreg, qo2 = po1_kreg)
	    if (qo1->iobj == SISLPOINT)
	      qedge[kj] = SISL_NULL;	/* Not necessary to compute edge intersection.*/
	    else if (qo1->iobj == SISLCURVE)
	      {
		if ((qedge[kj] = newEdge (2)) == SISL_NULL)
		  goto err101;

		for (ki = 0; ki < 2; ki++)
		  {
		    if (qo1->edg[ki] == SISL_NULL)
		      {
			if ((qo1->edg[ki] = newObject (SISLPOINT)) == SISL_NULL)
			  goto err101;

			/* Pick out end point from a curve. */

			s1438 (qo1->c1, ki, &(qo1->edg[ki]->p1), &tpar, &kstat);
			if (kstat < 0)
			  goto error;
		      }
		    else
		      tpar = (ki == 0 ? qo1->c1->et[qo1->c1->ik - 1] :
			      qo1->c1->et[qo1->c1->in]);

		    /* Recursiv computing of end intersection. */

		    sh1761 (kj == 0 ? qo1->edg[ki] : qo2, kj == 0 ? qo2 : qo1->edg[ki],
			    aepsge, &qintdat, &kstat);
		    if (kstat < 0)
		      goto error;


		    if (kstat)
		      {
			*jstat = 1;	/* Mark intersection found. */

			/* Put intersection found on edges into pintdat. */


			/* Get parameter border values of qo2. */

			/* UJK , newi */
			sh1782 (po1_kreg, po2_kreg, aepsge, qintdat, kedg, tpar,
				pintdat, &idummy, &kstat);
			if (kstat < 0)
			  goto error;


			/* Uppdate edge structure. */
			s6idedg (po1_kreg, po2_kreg, kj + 1, 1, tpar, *pintdat,
			  &qedge[kj]->prpt[ki], &qedge[kj]->ipoint, &kstat);
			if (kstat < 0)
			  goto error;
		      }

		    if (qintdat != SISL_NULL)
		      freeIntdat (qintdat);
		    qintdat = SISL_NULL;
		  }
		kedg++;
	      }

	    else if (qo1->iobj == SISLSURFACE)
	      {
		int kpar;	/* Parameter direction.            */

		if ((qedge[kj] = newEdge (4)) == SISL_NULL)
		  goto err101;

		for (ki = 0; ki < 4; ki++)
		  {
		    if (qo1->edg[ki] == SISL_NULL)
		      {
			if ((qo1->edg[ki] = newObject (SISLCURVE)) == SISL_NULL)
			  goto err101;

			/* Pick out edge curve from a surface. */

			s1435 (qo1->s1, ki, &(qo1->edg[ki]->c1), &tpar, &kstat);
			if (kstat < 0)
			  goto error;
		      }
		    else
		      tpar = (ki == 0 ? qo1->s1->et2[qo1->s1->ik2 - 1] :
			      (ki == 1 ? qo1->s1->et1[qo1->s1->in1] :
			       (ki == 2 ? qo1->s1->et2[qo1->s1->in2] :
				qo1->s1->et1[qo1->s1->ik1 - 1])));

		    /* Recursiv computing of edge intersection. */

		    if (kj == 0)
		      sh1761 (qo1->edg[ki], qo2, aepsge, &qintdat, &kstat);
		    else
		      sh1761 (qo2, qo1->edg[ki], aepsge, &qintdat, &kstat);
		    if (kstat < 0)
		      goto error;


		    if (kstat)
		      {
			*jstat = 1;	/* Mark intersection found. */

			/* Compute Parameter direction of edge. */

			kpar = ((ki == 0 || ki == 2) ? 2 : 1);


			/* Put intersection found on edges into pintdat. */

			/* UJK , newi */
			sh1782 (po1_kreg, po2_kreg, aepsge, qintdat, kpar + kedg - 1, tpar,
				pintdat, &idummy, &kstat);
			if (kstat < 0)
			  goto error;


			/* Uppdate edge structure. */
			s6idedg (po1_kreg, po2_kreg, kj + 1, kpar, tpar, *pintdat,
			  &qedge[kj]->prpt[ki], &qedge[kj]->ipoint, &kstat);
			if (kstat < 0)
			  goto error;
		      }

		    if (qintdat != SISL_NULL)
		      freeIntdat (qintdat);
		    qintdat = SISL_NULL;
		  }
		kedg += 2;
	      }

	    else
	      goto err121;

	  /* Before we enter internal intersection and subdivision we
	     initiate pointers to top level objects. */

	  if (po1_kreg->o1 == SISL_NULL)
	    po1_kreg->o1 = po1_kreg;
	  if (po2_kreg->o1 == SISL_NULL)
	    po2_kreg->o1 = po2_kreg;

	  /* Find the intersections in the inner of the object.  */

	  /* NEWI (ujk) Must treat helppoint on edges */
	  if (qedge[0] != SISL_NULL)
	    freeEdge (qedge[0]);
	  qedge[0] = SISL_NULL;
	  if (qedge[1] != SISL_NULL)
	    freeEdge (qedge[1]);
	  qedge[1] = SISL_NULL;

	  if (po1_kreg->iobj == SISLPOINT)
	    qedge[0] = SISL_NULL;
	  else if ((qedge[0] = newEdge (2 * po1_kreg->iobj)) == SISL_NULL)
	    goto err101;

	  if (po2_kreg->iobj == SISLPOINT)
	    qedge[1] = SISL_NULL;
	  else if ((qedge[1] = newEdge (2 * po2_kreg->iobj)) == SISL_NULL)
	    goto err101;

	  sh6idalledg (po1_kreg, po2_kreg, *pintdat, qedge, &kstat);
	  if (kstat < 0)
	    goto error;

	  kstat = (kxintercept) ? 202 : 0;
	  sh1762 (po1_kreg, po2_kreg, aepsge, pintdat, qedge, &kstat);
	  if (kstat < 0)
	    goto error;
	  else if (kstat)
	    *jstat = 1;


	  /* Free the edges used in s1762. */

	  if (qedge[0] != SISL_NULL)
	    freeEdge (qedge[0]);
	  qedge[0] = SISL_NULL;
	  if (qedge[1] != SISL_NULL)
	    freeEdge (qedge[1]);
	  qedge[1] = SISL_NULL;


	  /* UJK, edge reduction rules */
	     sh6edgred (po1_kreg, po2_kreg, (*pintdat), &kstat);

	  /* Organize the list in pintdat. */
	    sh6idlis (po1_kreg, po2_kreg, pintdat, aepsge, &kstat);
	  if (kstat < 0)
	    goto error;

          /* Convert any degenerate intersection curves (Intlists)
	     to intersection points (Intpts).
	     This may be necessary when
	     either of the intersection objects po1_kreg and po2_kreg
	     are parametrically degenerate. */
/*
	  sh6degen(po1_kreg,po2_kreg,pintdat,aepsge,&kstat);
	  if (kstat < 0) goto error;
*/
	  if (qintdat)  freeIntdat(qintdat);  /* PFU                  */
	  qintdat = SISL_NULL;                     /* - to fix memory leak */
	}
      else
	*jstat = 0;
    }

  goto out;

  /* Error in space allocation.  */

err101:*jstat = -101;
  s6err ("sh1761", *jstat, kpos);
  goto out;

  /* Error. Dimensions conflicting.  */

err106:*jstat = -106;
  s6err ("sh1761", *jstat, kpos);
  goto out;

  /* Error. Kind of object does not exist.  */

err121:*jstat = -121;
  s6err ("sh1761", *jstat, kpos);
  goto out;

  /* Error in lower order routine.  */

error:*jstat = kstat;
  s6err ("sh1761", *jstat, kpos);
  goto out;

  out:
     /* Kill kreg temporary geometry */
     if (po1_kreg && po1_kreg != po1)
       {
	  freeObject(po1_kreg);
	  po1_kreg = SISL_NULL;
       }

     if (po2_kreg && po2_kreg != po2)
       {
	  freeObject(po2_kreg);
	  po2_kreg = SISL_NULL;

       }

     if (qt)  freeIntpt(qt);             /* PFU                   */
     if (qintdat)  freeIntdat(qintdat);  /* - to fix memory leak. */
     if (qedge[0])  freeEdge(qedge[0]);
     if (qedge[1])  freeEdge(qedge[1]);

}
