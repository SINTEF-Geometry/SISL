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
 * $Id: sh1762.c,v 1.18 2007-08-06 13:09:12 vsk Exp $
 *
 */


#define SH1762

#include "sislP.h"

/*
extern int nmbcall;
extern int nmb0;
extern int nmb1;
extern int nmb2;
extern int nmbsuccess;
*/

/*
* Level counter.
* --------------
*/

static int xc = 0;
static int xmax = 0;

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void sh1762_s9mic (SISLObject *, SISLObject *, SISLIntdat **, SISLEdge **[], int *);
static void sh1762_s9num (SISLObject *, SISLObject *, int *, int *);
static void sh1762_s9div (SISLObject *, SISLObject *, double, int, int, SISLObject *[], SISLEdge *[], SISLIntdat **, int *);
static void sh1762_s9subdivpt (SISLObject *, SISLObject *, double, int, int, SISLEdge *[], SISLIntdat **, int *, SISLIntpt **, double[], int *);
static void sh1762_s9update (SISLObject *, SISLObject *, double, SISLIntdat **, SISLEdge **[], int *);
static void sh1762_s9con (SISLObject *, SISLObject *, double, SISLIntdat **, SISLEdge *[], int *);
static void sh1762_s9intercept (SISLObject *, SISLObject *, double, int, SISLIntpt *[], int *);
static void sh1762_s9coincide (SISLObject *, SISLObject *, double, int, SISLIntpt *[], int *);
static void sh1762_s9toucharea (SISLObject *, SISLObject *, double, int, SISLIntpt *[], int *);
/*static void sh1762_s9edgpoint (SISLEdge *[], SISLIntpt ***, int *, int *); */
static void sh1762_s9edgsscon (SISLEdge *[], SISLSurf *, SISLSurf *, SISLIntdat *, int, double, int *);
static void sh1762_s9edgpscon (SISLEdge *, double, SISLSurf *, int, SISLIntdat *, double, int *);
static void sh1762_s9simple (SISLObject *, SISLObject *, SISLEdge *[], int *);
/* static void sh1762_s9reex (SISLObject *, SISLObject *, SISLEdge *[], double, SISLIntdat *, int *); */
static void sh1762_s9ptiter (SISLObject *, SISLObject *, double, SISLIntdat **, SISLEdge *[], int *);
static int sh1762_is_taboo(SISLSurf *, SISLSurf *, SISLIntpt *, int, int *);
static double sh1762_sflength(SISLSurf *, int, int *);
#else
static void sh1762_s9mic ();
static void sh1762_s9num ();
static void sh1762_s9div ();
static void sh1762_s9subdivpt ();
static void sh1762_s9update ();
static void sh1762_s9con ();
static void sh1762_s9intercept ();
static void sh1762_s9coincide ();
static void sh1762_s9toucharea ();
/*static void sh1762_s9edgpoint ();*/
static void sh1762_s9edgsscon ();
static void sh1762_s9edgpscon ();
static void sh1762_s9simple ();
/* static void sh1762_s9reex (); */
static void sh1762_s9ptiter ();
static int sh1762_is_taboo();
static double sh1762_sflength();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
sh1762 (SISLObject * po1, SISLObject * po2, double aepsge,
	SISLIntdat ** pintdat, SISLEdge * vedge[], int *jstat)
#else
void
sh1762 (po1, po2, aepsge, pintdat, vedge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **pintdat;
     SISLEdge *vedge[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*          NOTE : Comments for further developments/tasks starts
*                 with UPDATE :
*
*
* PURPOSE    : SISLObject - object intersection. Treat the inner of the
*              object.
*
*
*
* INPUT      : po1       - Pointer to first object
*              po2       - Pointer to second object
*              aepsge    - Geometry resolution.
*              vedge[2]  - Pointers to structure of edge-intersections.
*              *jstat    - Flag
*                          = 202 : Complicated point-surface intersection
*                                  in 3D. Perform extra interception test.
*
*
* INPUT/OUTPUT : pintdat - Pointer to intersection data.
*
*
*
* OUTPUT     : jstat     - status messages
*                                         = 1      : intersection found
*                                         = 0      : no intersection
*                                         < 0      : error
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6err      - Treating error situation.
*              s1741      - Simple Case test for intersections.
*              sh1790     - Box test.
*              s1791      - Test if possible to subdivide
*              s1792      - Computing midpoint of parameter intervall.
*              s1770      - Curve/curve iteration.
*              s1771      - Point/curve iteration.
*              s1772      - Curve/surface iteration.
*              s1773      - Point/surface iteration.
*              s1231      - Subdivide curve.
*              s1711      - Subdivide surface.
*              sh1761     - Object/object intersection.
*              s1435      - Pick an edge curve from a surface.
*              s1438      - Pick an end point from a curve.
*              sh6idnpt    - New intpoint in intdat.
*              s6idkpt    - Kill intpoint in intdat.
*              s6idcpt    - Copy intpoint in intdat.
*              sh6idcon    - Connect intpoint in intdat.
*              s6idput    - put intpoint from one intdat to another
*                           intdat with one paramete more.
*              s6idput    - get intpoint from one intdat to another
*                           intdat with one paramete less.
*              s6idedg    - Get intpoint on edges from intdat.
*              s6idint    - Get internal intpoint from intdat.
*              newPoint   - Create new point structure.
*              newObject  - Create new object structure.
*              hp_newIntpt   - Create new intpt structure.
*              newEdge    - Create new edge structure.
*              freeObject - Free space occupied by a given object.
*              freeCurve  - Free space occupied by a given curve.
*              freeEdge   - Free space occupied by a given edge.
*              freeIntdat - Free space occupied by a given intdat.
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*              UJK        , SI, 89-04.
*              Arne Laksaa, SI, 89-06.
*              UJK newi         91-06.
*********************************************************************
*/
{
  int kpos = 0;			/* Position of error.                 */
  int kstat = 0;		/* Local error status.                */
  int kdiv1 = 0;		/* Parameter direction of subdivsion. */
  int kdiv2 = 0;		/* Parameter direction of subdivsion. */
  int ki, ki1, ki2;		/* Counters.                          */
  int at_bottom=TRUE;           /* Flag, true on bottom level of recur*/
  int knewpt=0;                 /* No of points made in prtop part    */
  int kexpand = 2;		/* Expand box in the inner of object. */
  int kxintercept = (*jstat == 202);  /* Extra interception           */
  /* int knum;  */                   /* Number of intersection points at edges. */
  SISLObject *uob1[4];		/* Pointers to subdivided object.     */
  SISLObject *uob2[4];		/* Pointer to object to subdivide.    */

  int debug_flag=0;

  /*  FOR DEBUGGING define debug_flag as an extern variable, i.e.:
   *
   *                    extern int debug_flag;
   */

    if (debug_flag)
    {
       if ((po1->iobj == SISLSURFACE && po1->s1->idim == 1) ||
           (po1->iobj == SISLSURFACE  && po2->iobj == SISLSURFACE))
	   {
    	           /*	if (po1->s1->et1[0] >= 3.3 &&
		        po1->s1->et1[po1->s1->in1] <= 3.6 &&
		        po1->s1->et2[0] >= 0.7 &&
		        po1->s1->et2[po1->s1->in2] <= 0.9)
		        {
		   */
	   int knum;
	   int ipar = 2;
	   int kj, ki;
	   SISLIntpt **up = SISL_NULL;  /* Array of poiners to intersection point.*/

	   sh6edgpoint (vedge, &up, &knum, &kstat);
	   if (kstat < 0)
	      goto error;
	   if (debug_flag == 1)
	   {
	      printf("\n___________________________________________________");

	      printf("\n par val(1) :%#10.10g %#10.10g %#10.10g %#10.10g ",
		     po1->s1->et1[0],
		     po1->s1->et1[po1->s1->in1],
		     po1->s1->et2[0],
		     po1->s1->et2[po1->s1->in2]);
	      if (po2->iobj == SISLSURFACE)
	      {
		 ipar = 4;
		 printf("\n par val(2) :%#10.10g %#10.10g %#10.10g %#10.10g ",
			po2->s1->et1[0],
			po2->s1->et1[po2->s1->in1],
			po2->s1->et2[0],
			po2->s1->et2[po2->s1->in2]);
	      }
	      printf("\n No of pts: %d",knum);
	      for (ki = 0; ki < knum; ki++)
	      {
		 printf("\n point %d :",ki);
		 for (kj = 0; kj < ipar; kj++)
		    printf(" %#10.10g", up[ki]->epar[kj]);
	      }
	   }
	   else              /* if (debug_flag == 2) */
	   {
	      printf("fg: black \n");
	      printf("lin: \n%#10.10g %#10.10g \n",
		     po1->s1->et1[0],
		     po1->s1->et2[0]);

	      printf("%#10.10g %#10.10g \n",
		     po1->s1->et1[0],
		     po1->s1->et2[po1->s1->in2]);

	      printf("%#10.10g %#10.10g \n",
		     po1->s1->et1[po1->s1->in1],
		     po1->s1->et2[po1->s1->in2]);

	      printf("%#10.10g %#10.10g \n",
		     po1->s1->et1[po1->s1->in1],
		     po1->s1->et2[0]);

	      printf("%#10.10g %#10.10g \n",
		     po1->s1->et1[0],
		     po1->s1->et2[0]);

	   }

	   if (up) freearray(up);

	                     /*	   }  */
     }
  }

  xc++;
  xmax = MAX (xmax, xc);
  /*  printf("Max : %d \n",xc); */


  for (ki = 0; ki < 4; ki++)
    uob1[ki] = uob2[ki] = SISL_NULL;

  /* Initiate to no intersection. */

  *jstat = 0;

  /* Test if intersection is possible (perform box-test).  */

  /*  box_nmb++;
  time_before = clock();  */

  sh1790 (po1, po2, kexpand, aepsge, &kstat);

  /*  time_used = clock() - time_before;
  box_time += time_used; */
  if (kstat < 0)
    goto error;

  /*  printf("Box test. Status = %d \n",kstat); */

  /* We may have tree different values on kstat.
     kstat = 1 : The two boxes overlapp.
     kstat = 2 : The two "bezier" boxes is just touching.
     kstat = 3 : The two boxes is both inside a microbox of aepsge.
     kstat = 4 : One of the objects is degenerated to one 3D point.
     kstat = 5 : Danger of shadow area in point object intersection,
                 dimension > 1.   */

  if (kstat == 5)
  {
     /* VSK, 92-10.
	Either make sure that there is no overlap, or find the intersection.
	The situation that there is an intersection point in point-object
	intersection when dim > 1 where the usual box test fails to
	recognize the possibility may arise near the endpoints/edges of
	the other object. */

     sh1762_s9ptiter(po1, po2, aepsge, pintdat, vedge, &kstat);
     if (kstat < 0) goto error;

     /* kstat = 0 : No overlap.
	kstat = 1 : The boxes overlap, and the intersection is found. */

     if (kstat == 1) *jstat = 1;
  }

  else if (kstat == 4)

    goto out;

  else if (kstat == 3)
    {
      /* Microbox found.*/

      sh1762_s9mic (po1, po2, pintdat, &vedge, &kstat);
      if (kstat < 0)
	goto error;
      else
	*jstat = kstat;		/* Possible uppdating intersection found. */
    }
  else if (kstat == 1)
    {
      /* Simple Case test (more than one intersection possible?)  */

      /* UJK, div until bezier, due to problems in silhouettes */
       /* Must be opened again for silhouettes NO/YES?/NO!/...
	  ???????????????????????????????????
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       if ((po1->iobj == SISLSURFACE && po1->s1->idim == 1 &&
	    (po1->s1->ik1 != po1->s1->in1 || po1->s1->ik2 != po1->s1->in2)) ||
	   (po2->iobj == SISLSURFACE && po2->s1->idim == 1 &&
	    (po2->s1->ik1 != po2->s1->in1 || po2->s1->ik2 != po2->s1->in2)))
	  kstat = 0;
       else
       { */

	  s1741 (po1, po2, aepsge, &kstat);
	  if (kstat < 0)
	    goto error;
	  //else if (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE &&
	  //	   vedge[0]->ipoint + vedge[1]->ipoint > 0 && !kstat)
	  //  sh1762_s9simple (po1, po2, vedge, &kstat);
	  if (kstat < 0)
	    goto error;
	  /* } */
      /* We may have two different values on kstat.
	 kstat = 0 : No simple case.
	 kstat = 1 : Simple case (surfaces possible simple case). */

      if (kstat == 1)
      {
	 /* Possible simple Case, update intersection list. */

	 sh1762_s9update (po1, po2, aepsge, pintdat, &vedge, &kstat);
	 if (kstat < 0)
	    goto error;

	 /* We may have two different values on kstat.
	    kstat = 0 : No simple case, more than two edge intersection.
	    kstat = 1 : Intersection found. */

	 if (kstat == 1)
	    *jstat = 1;		/*Updating found intersection. */
      }

      /* UJK,20.01.93, Don't skip s9con when not success in s9update.
	 removed else.*/
      if (kstat ==0)
      {
	 /* UJK, 17.12.92, for a 1D surface of bezier type
	    there may be a posibility of dividing out edge
	    curve intersections */
	 if (po1->iobj == SISLSURFACE && po1->s1->idim ==1)
	 {
	    sh_1d_div(po1, po2, aepsge, pintdat, vedge, &kstat);
	    if (kstat < 0)
	       goto error;
	    if (kstat == 1)
	       *jstat = 1;		/*Updating found intersection. */
	 }
	 else if (po2->iobj == SISLSURFACE && po2->s1->idim == 1)
	 {
	    sh_1d_div(po2, po1, aepsge, pintdat, vedge, &kstat);
	    if (kstat < 0)
	       goto error;
	    if (kstat == 1)
	       *jstat = 1;		/*Updating found intersection. */
	 }

	 else
	 {

	    /* Check for interval intersection. */

	   kstat = (kxintercept) ? 202 : 0;
	    sh1762_s9con (po1, po2, aepsge, pintdat, vedge, &kstat);
	    if (kstat < 0)
	       goto error;

	    /*  printf("sh1762_s9con. Status = %d \n",kstat); */

	    /* We may have two different values on kstat.
	       kstat = 0 : No intervall intersection.
	       kstat = 1 : Intervall intersection found.
	       kstat = 2 : Intersection not possible */

	    if (kstat == 1)
	       *jstat = 1;		/*Updating found intersection. */
	 }
      }


      if (kstat == 0)
	{
	  /* Find number of possible subdivision directions.
	     kdiv1 and kdiv2 may have 4 difference values :
	     kdiv = 0 : Subdivision not possible.
	     kdiv = 1 : Subdivision in first parameter direction.
	     kdiv = 2 : Subdivision in second parameter direction.
	     kdiv = 3 : Subdivision in both parameter directions. */

	  sh1762_s9num (po1, po2, &kdiv1, &kstat);
	  if (kstat < 0)
	    goto error;

	  sh1762_s9num (po2, po1, &kdiv2, &kstat);
	  if (kstat < 0)
	    goto error;

	  if (kdiv1 + kdiv2 == 0)
	    {
	      /* There is two almost plane parallel objects, and
		 there is nothing at the edges (otherwise the
		 intersections should be found by s9con). Then the
		 only possibility is that there is no intersection. */
	       /* VSK, 11-92. Since partial coincidence is not
		  implemented, there might be intersections on the
		  edges. Check this.   This should not be necessary
		  any more.

	        Check if there are intersection points on edges.

	       if (vedge[0] == SISL_NULL)
		  knum = 0;
	       else
		  knum = vedge[0]->ipoint;

	       if (vedge[1] != SISL_NULL)
		  knum += vedge[1]->ipoint;



	       if (knum > 0)
	       {
		   Do something that makes the routine terminate
		     until partial coincidence is implemented.

		  sh1762_s9mic(po1, po2, pintdat, &vedge, &kstat);
		  if (kstat < 0) goto error;

		  *jstat = kstat;
	       }
	       else
	       {
		  *jstat = 0;
		  goto out;
	       } */

	       *jstat = 0;
	       goto out;
	    }
	  else
	    {
	      SISLEdge *uedge[2];	/* Array of pointers to edges
					      to use in subproblems.    */


	      /* We do not have simple case and it is possible to
		 subdivide. We therefor subdivide and update the
		 edge intersection and then do a recurcive call
		 to treat the sub problems. Curves are subdivided
		 into two, surfaces into four. We can therefor get
		 up to sexteen recurcive calls.*/



	      /***** Treating objects on sub problems. *****/

	      if (kdiv1 > 0)	/* New objects for subdivision of po1. */
		{
		  for (ki = 0; ki < (kdiv1 < 3 ? 2 : 4); ki++)
		    {
		      if ((uob1[ki] = newObject (po1->iobj)) == SISL_NULL)
			goto err101;

		      /* Initiate o1 pointer to point to top level object. */

		      uob1[ki]->o1 = po1->o1;
		    }

		  /* Subdivide the po1 object. */

		  sh1762_s9div (po1, po2, aepsge, 1, kdiv1, uob1, vedge, pintdat, &kstat);
		  if (kstat < 0)
		    goto error;
		  else if (kstat == 1)
		    *jstat = 1;
		}


	      if (kdiv2 > 0)	/* New objects for subdivision of po2. */
		{
		  for (ki = 0; ki < (kdiv2 < 3 ? 2 : 4); ki++)
		    {
		      if ((uob2[ki] = newObject (po2->iobj)) == SISL_NULL)
			goto err101;

		      /* Initiate o1 pointer to point to top level object. */

		      uob2[ki]->o1 = po2->o1;
		    }

		  /* Subdivide the po2 object. */

		  sh1762_s9div (po1, po2, aepsge, 2, kdiv2, uob2, vedge, pintdat, &kstat);
		  if (kstat < 0)
		    goto error;
		  else if (kstat == 1)
		    *jstat = 1;
		}


	      /***** Recursion. *****/

	      if (kdiv1 == 0)	/* Only second object subdivided. */
		for (ki = 0; ki < (kdiv2 < 3 ? 2 : 4); ki++)
		  {
		    /***** Treating edges on sub problems. *****/

		    /* Making new edge object to sub problems. */

		    if (po1->iobj == SISLPOINT)
		      uedge[0] = SISL_NULL;
		    else if ((uedge[0] = newEdge (vedge[0]->iedge)) == SISL_NULL)
		      goto err101;
		    if ((uedge[1] = newEdge (vedge[1]->iedge)) == SISL_NULL)
		      goto err101;

		    /* Update edge intersection on sub problems. */

		    sh6idalledg (po1, uob2[ki], *pintdat, uedge, &kstat);
		    if (kstat < 0)
		      goto error;

		    at_bottom = FALSE;
		    kstat = (kxintercept) ? 202 : 0;
		    sh1762 (po1, uob2[ki], aepsge, pintdat, uedge, &kstat);
		    if (kstat < 0)
		      goto error;
		    else
		      *jstat = *jstat || kstat;

		    if (uedge[0] != SISL_NULL)
		      freeEdge (uedge[0]);
		    if (uedge[1] != SISL_NULL)
		      freeEdge (uedge[1]);
		  }
	      else if (kdiv2 == 0)	/* Only first object subdivided.   */
		for (ki = 0; ki < (kdiv1 < 3 ? 2 : 4); ki++)
		  {
		    /***** Treating edges on sub problems. *****/

		    /* Making new edge object to sub problems. */

		    if ((uedge[0] = newEdge (vedge[0]->iedge)) == SISL_NULL)
		      goto err101;
		    if (po2->iobj == SISLPOINT)
		      uedge[1] = SISL_NULL;
		    else if ((uedge[1] = newEdge (vedge[1]->iedge)) == SISL_NULL)
		      goto err101;

		    /* Update edge intersection on sub problems. */

		    sh6idalledg (uob1[ki], po2, *pintdat, uedge, &kstat);
		    if (kstat < 0)
		      goto error;

		    at_bottom = FALSE;
		    kstat = (kxintercept) ? 202 : 0;
		    sh1762 (uob1[ki], po2, aepsge, pintdat, uedge, &kstat);
		    if (kstat < 0)
		      goto error;
		    else
		      *jstat = *jstat || kstat;

		    if (uedge[0] != SISL_NULL)
		      freeEdge (uedge[0]);
		    if (uedge[1] != SISL_NULL)
		      freeEdge (uedge[1]);
		  }
	      else		/* Both objects subdivided.        */
		for (ki1 = 0; ki1 < (kdiv1 < 3 ? 2 : 4); ki1++)
		  for (ki2 = 0; ki2 < (kdiv2 < 3 ? 2 : 4); ki2++)
		    {
		      /***** Treating edges on sub problems. *****/

		      /* Making new edge object to sub problems. */

		      if ((uedge[0] = newEdge (vedge[0]->iedge)) == SISL_NULL)
			goto err101;
		      if ((uedge[1] = newEdge (vedge[1]->iedge)) == SISL_NULL)
			goto err101;

		      /* Update edge intersection on sub problems. */

		      sh6idalledg (uob1[ki1], uob2[ki2], *pintdat, uedge, &kstat);
		      if (kstat < 0)
			goto error;


		      at_bottom = FALSE;
		      kstat = (kxintercept) ? 202 : 0;
		      sh1762 (uob1[ki1], uob2[ki2], aepsge, pintdat, uedge, &kstat);
		      if (kstat < 0)
			goto error;
		      else
			*jstat = *jstat || kstat;

		      if (uedge[0] != SISL_NULL)
			freeEdge (uedge[0]);
		      if (uedge[1] != SISL_NULL)
			freeEdge (uedge[1]);
		    }
	    }
	}
    }


  /* Must update vedge before going into reex */
  /* if (vedge[0] != SISL_NULL)
  {
     knedge1 = vedge[0]->iedge;
     freeEdge (vedge[0]);
     if ((vedge[0] = newEdge (knedge1)) == SISL_NULL)
        goto err101;
  }
  if (vedge[1] != SISL_NULL)
  {
     knedge2 = vedge[1]->iedge;
     freeEdge (vedge[1]);
     if ((vedge[1] = newEdge (knedge2)) == SISL_NULL)
	goto err101;
  }*/

  /* Making new edge object to sub problems. */


  /* sh6idalledg (po1, po2, *pintdat, vedge, &kstat);
  if (kstat < 0)
     goto error; */


  /* UPDATE (ujk): s9reex must be changed, interface = ?
     Now it connects points on edge when they are
     connected to an internal point ?*/
  /* Now changed! ALA and MSF.  */

  /* UJK, VSK, ALA, 09.02.93, don't need it any longer !? */
  /* sh1762_s9reex (po1, po2, vedge, aepsge, *pintdat, &kstat);
     if (kstat < 0)
     goto error; */

  /* VSK, 10.92. Set status if reex takes action.  */

  /* *jstat = MAX(*jstat,kstat);

  if (debug_flag && kstat)
     printf("\n Output reex: %d \n",kstat); */

  /* Reduction rules */

  sh6red (po1, po2, (*pintdat), &kstat);

  /* Make help points and pretopology at bottom */

  if (at_bottom)
    shmkhlppts (po1, po2, aepsge, pintdat, vedge, &knewpt, &kstat);

  /* UJK, aug.92, If we make help points, status must be set !,
     are there other updating statuses that we've missed ? */
  if (knewpt) *jstat = 1;

  /* Intersections in the inner found.  */

  goto out;

  /* Error in space allocation.         */

err101:*jstat = -101;
  s6err ("sh1762", *jstat, kpos);
  goto out;

  /* Error in lower level routine.      */

error:*jstat = kstat;
  s6err ("sh1762", *jstat, kpos);
  goto out;

  /* Free the space that is  allocated. */

out:
  for (ki = 0; ki < 4; ki++)
    {
      if (uob1[ki] != SISL_NULL)
	freeObject (uob1[ki]);
      if (uob2[ki] != SISL_NULL)
	freeObject (uob2[ki]);
    }
  xc--;

}

#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9mic (SISLObject * po1, SISLObject * po2, SISLIntdat ** rintdat,
	      SISLEdge ** vedge[], int *jstat)
#else
static void
sh1762_s9mic (po1, po2, rintdat, vedge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat **rintdat;
     SISLEdge **vedge[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE     : Treat intersection when it is not possible to
*               subdivide any futher, and it is not simple case.
*
*
*
* INPUT      : (*vedge)[2] - SISLEdge intersection objects to the two
*                         objects in intersection problem.
*
*
*
* OUTPUT     : rintdat - intersection data.
*              jstat   - status messages
*                               = 1     : Intersection found.
*                               = 0     : Intersection not found.
*                               < 0     : error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*
*********************************************************************
*/
{
  int kpos = 0;			/* Position of error.                      */
  int kstat = 0;		/* Local error status.                     */
  int knum = 0;			/* Number of intpt on edges.               */
  /*int klist1, klist2;	*/	/* List index in iintpt.                   */
  int ind1, ind2;		/* Help index in up array.                 */
  double *spar = SISL_NULL;		/* Array to store parameter values.        */
  SISLIntpt **up = SISL_NULL;	/* Array of poiners to intersection point. */
  double *nullp = SISL_NULL;
  double tepsge = 0.0000001;    /* Tolerance used in merging of points.    */

  /* Initiate to no new intersection point. */

  *jstat = 0;

  /* Compute number of intersection points on edges. */

  if ((*vedge)[0] == SISL_NULL)
    knum = 0;
  else
    knum = (*vedge)[0]->ipoint;

  if ((*vedge)[1] != SISL_NULL)
    knum += (*vedge)[1]->ipoint;


  if (knum > 0)
    {
      /* sh1762_s9edgpoint ((*vedge), &up, &knum, &kstat); */
      sh6edgpoint ((*vedge), &up, &knum, &kstat);
      if (kstat < 0)
	goto error;
    }


  if (knum > 1)
    {
      int kturn, ki;

      /* We have more than one intersection point on the edges,
	 we therefor have to treat these problem. */

      if ((po1->iobj == SISLPOINT && po1->p1->idim <= 2) ||
	  (po2->iobj == SISLPOINT && po2->p1->idim <= 2) ||
	  (po1->iobj == SISLCURVE && po2->iobj == SISLPOINT && knum == 2) ||
	  (po1->iobj == SISLPOINT && po2->iobj == SISLCURVE && knum == 2))
	{
	  SISLObject *qo1, *qo2;

	  /* In dimension one and two this function is not
	     a degenenerate treatment function, it is a coincidence
	     function. */

	  if (po1->iobj == SISLPOINT)
	    {
	      qo1 = po1;
	      qo2 = po2;
	      kturn = 0;
	    }
	  else
	    {
	      qo2 = po1;
	      qo1 = po2;
	      kturn = 1;
	    }

	  if (qo2->iobj == SISLSURFACE)
	    {

	       /* Trim area found */
	       for (ki=0; ki<(*rintdat)->ipoint; ki++)
	       {
		  sh6isinside(po1,po2,(*rintdat)->vpoint[ki],&kstat);
		  if (kstat < 0) goto error;
		  if (kstat)
		  {
		     sh6tomain((*rintdat)->vpoint[ki], &kstat);
		     if (kstat < 0) goto error;
		     (*rintdat)->vpoint[ki]->iinter = SI_TRIM;
		  }
	       }

	       /* UJK 18.09.90 Must set intersection found status */
	       *jstat = 1;
	       goto out;

	    }
	  else if (qo2->iobj == SISLCURVE && knum == 2)
	    {
	      double tres;

	      tres = (qo2->c1->et[qo2->c1->in] -
		      qo2->c1->et[qo2->c1->ik - 1]) /
		(qo2->o1->c1->et[qo2->o1->c1->in] -
		 qo2->o1->c1->et[qo2->o1->c1->ik - 1]);

	      if (tres > REL_PAR_RES)
		{
		  /* UJK newi :Main points, curve point 1+2D, connect */
		  sh6idcon (rintdat, up, up + 1, &kstat);
		  if (kstat < 0)
		    goto error;
		  /* Sort points */
		  ind1 = 0;
		  ind2 = 1;
		  if (up[0]->epar[0] > up[1]->epar[0])
		    {
		      ind1 = 1;
		      ind2 = 0;
		    }
		  sh6setdir (up[ind1], up[ind2], &kstat);
		  if (kstat < 0)
		    goto error;



		  /* Set pretopology */
		  /* No, it's there already */

		  /*		  ind1 = 1;
	          ind2 = 0;
	          if (up[0]->epar[0] < up[1]->epar[0])
	        {
	          ind1 = 0;
	          ind2 = 1;
	          }

	          sh6getlist (up[ind1], up[ind2], &klist1, &klist2, &kstat);
	          if (kstat != 0)
	        {
	          kstat = -1;
	          goto error;
	          }
	          if (kturn)
	        {
	          sh6settop (up[ind1], -1,
	          SI_AT, SI_ON, SI_UNDEF, SI_ON, &kstat);
	          if (kstat < 0)
	          goto error;
	          sh6settop (up[ind2], -1,
	          SI_ON, SI_AT, SI_ON, SI_UNDEF, &kstat);
	          if (kstat < 0)
	          goto error;
	          }
	          else
	        {
	          sh6settop (up[ind1], -1,
	          SI_UNDEF, SI_ON, SI_AT, SI_ON, &kstat);
	          if (kstat < 0)
	          goto error;
	          sh6settop (up[ind2], -1,
	          SI_ON, SI_UNDEF, SI_ON, SI_AT, &kstat);
	          if (kstat < 0)
	          goto error;
	          }

	          */

		  /* UJK 18.09.90 Must set intersection found status */
		  *jstat = 1;
		  goto out;
		}
	    }
	}


	/* VSK to treat degenerated curves.  */

	if (po1->iobj == SISLCURVE && po2->iobj == SISLCURVE && knum >= 2)
	{
	   /* The two curves is within a microbox. The intersection will
	      be represented with two points that are connectd. Merge
	      the rest of the points into one of the two remaining.   */

	   for (ki=1; ki<knum-1; ki++)
	     {
		sh6idnewunite(po1,po2,rintdat,&up[0],&up[ki],DZERO,
			      tepsge,&kstat);
		if (kstat < 0) goto error;
             }

	   sh6connect(up[0],up[knum-1],&kstat);
	   if (kstat < 0) goto error;

	   /* Update edge structure.  */

      if ((*vedge)[0] != SISL_NULL)
	{
	  ki = (*vedge)[0]->iedge;
	  freeEdge ((*vedge)[0]);
	  (*vedge)[0] = SISL_NULL;
	  if (((*vedge)[0] = newEdge (ki)) == SISL_NULL)
	    goto err101;
	}
      if ((*vedge)[1] != SISL_NULL)
	{
	  ki = (*vedge)[1]->iedge;
	  freeEdge ((*vedge)[1]);
	  (*vedge)[1] = SISL_NULL;
	  if (((*vedge)[1] = newEdge (ki)) == SISL_NULL)
	    goto err101;
	}

          sh6idalledg (po1, po2, *rintdat, *vedge, &kstat);
          if (kstat < 0)
            goto error;

           *jstat = 1;
	   goto out;
	}

      /* We have more than one intersection point on the edges.
	 We therefor kill these points and
	 try to find a new intersection point. */



      for (ki = 1; ki < knum; ki++)
	{
	  /* UJK newi, unite the points : */
	   sh6idnewunite (po1, po2, rintdat, &up[0], &up[ki], (double) 0.5,
			  tepsge, &kstat);
	  if (kstat < 0)
	    goto error;
	}

      if ((*vedge)[0] != SISL_NULL)
	{
	  ki = (*vedge)[0]->iedge;
	  freeEdge ((*vedge)[0]);
	  if (((*vedge)[0] = newEdge (ki)) == SISL_NULL)
	    goto err101;
	}
      if ((*vedge)[1] != SISL_NULL)
	{
	  ki = (*vedge)[1]->iedge;
	  freeEdge ((*vedge)[1]);
	  if (((*vedge)[1] = newEdge (ki)) == SISL_NULL)
	    goto err101;
	}
      /* UJK newi, one point kept : */
      knum = 1;
    }



  if (knum == 0)
    {
      int kpar = 0;
      SISLIntpt *qt;


      /* There is no intersection points on the edges.
	 We therfore make one new intersection point with parameter
	 values in senter of each object. */


      /* Number of parameter values of object 1. */

      if (po1->iobj == SISLCURVE)
	kpar = 1;
      else if (po1->iobj == SISLSURFACE)
	kpar = 2;
      else
	kpar = 0;

      /* Number of parameter values of object 2. */

      if (po2->iobj == SISLCURVE)
	kpar++;
      else if (po2->iobj == SISLSURFACE)
	kpar += 2;


      /* Allocate array to store midpoint parameter values. */

      if ((spar = newarray (kpar, double)) == SISL_NULL)
	goto err101;


      /* Compute midpoint parameter values. */

      if (po1->iobj == SISLCURVE)
	{
	  spar[0] = (po1->c1->et[po1->c1->ik - 1] +
		     po1->c1->et[po1->c1->in]) * (double) 0.5;
	  kpar = 1;
	}
      else if (po1->iobj == SISLSURFACE)
	{
	  spar[0] = (po1->s1->et1[po1->s1->ik1 - 1] +
		     po1->s1->et1[po1->s1->in1]) * (double) 0.5;
	  spar[1] = (po1->s1->et2[po1->s1->ik2 - 1] +
		     po1->s1->et2[po1->s1->in2]) * (double) 0.5;
	  kpar = 2;
	}
      else
	kpar = 0;

      if (po2->iobj == SISLCURVE)
	{
	  spar[kpar] = (po2->c1->et[po2->c1->ik - 1] +
			po2->c1->et[po2->c1->in]) * (double) 0.5;
	  kpar++;
	}
      else if (po2->iobj == SISLSURFACE)
	{
	  spar[kpar] = (po2->s1->et1[po2->s1->ik1 - 1] +
			po2->s1->et1[po2->s1->in1]) * (double) 0.5;
	  spar[kpar + 1] = (po2->s1->et2[po2->s1->ik2 - 1] +
			    po2->s1->et2[po2->s1->in2]) * (double) 0.5;
	  kpar += 2;
	}

      *jstat = 1;		/* Mark intersection found. */


      /* Making intersection point. */
      /* UJK newi */
      /* UPDATE: ? Be aware of this situation, can it occur ? */

      qt = hp_newIntpt (kpar, spar, DZERO, SI_ORD,
			SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
			0, 0, nullp, nullp);
      if (qt == SISL_NULL)
	goto err101;

      /* Uppdating pintdat. */

      sh6idnpt (rintdat, &qt, 1, &kstat);
      if (kstat < 0)
	goto error;
    }

  goto out;

/* Error in space allocation.         */

err101:*jstat = -101;
  s6err ("sh1762_s9mic", *jstat, kpos);
  goto out;

/* Error in lower level routine.      */

error:*jstat = kstat;
  s6err ("sh1762_s9mic", *jstat, kpos);
  goto out;

out:if (spar != SISL_NULL)
    freearray (spar);
  if (up != SISL_NULL)
    freearray (up);
}

#if defined(SISLNEEDPROTOTYPES)
static double sh1762_sflength(SISLSurf *psurf, int idir, int *jstat)
#else
static double sh1762_sflength(psurf, idir, jstat)
    SISLSurf *psurf;
    int idir;
    int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Estimate the surface length in a given parameter direction.
*
*
*
* INPUT      : psurf  - The surface.
*              idir   - The parameter direction.
*
*
* OUTPUT     : return value - Estimated surface length.
*              jstat  - status messages
*                         = 0     : no error.
*                         < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 99-05.
*
*********************************************************************
*/
{
  int kstat = 0;
  int kleft1 = 0, kleft2 = 0;
  int ki;
  int kdim = psurf->idim;
  double spar[2];  /* Parameter value in which to evaluate. */
  double sint[2];  /* Interval between parameter values.    */
  double sder[12]; /* Points on the surface.                */
  int kneval;      /* Number of points to evaluate.         */
  double tlength = 0.0;  /* Estimated length of surface.    */

  kneval = (idir == 1) ? psurf->ik1 : psurf->ik2;
  kneval = max(2, min(kneval, 4));

  /* Set first parameter in which to evaluate. */
  if (idir == 1)
    {
      spar[0] = psurf->et1[psurf->ik1-1];
      spar[1] = (double)0.5*(psurf->et2[psurf->ik2-1]+psurf->et2[psurf->in2]);

      sint[0] = (psurf->et1[psurf->in1] - spar[0])/(double)(kneval-1);
      sint[1] = 0.0;
    }
  else
    {
      spar[0] = (double)0.5*(psurf->et1[psurf->ik1-1]+psurf->et1[psurf->in1]);
      spar[1] = psurf->et2[psurf->ik2-1];

      sint[0] = 0.0;
      sint[1] = (psurf->et2[psurf->in2] - spar[1])/(double)(kneval-1);
    }

  /* Evaluate points. */

  for (ki=0; ki<kneval; ki++, spar[0]+=sint[0], spar[1]+=sint[1])
    {
      s1424(psurf, 0, 0, spar, &kleft1, &kleft2, sder+ki*kdim, &kstat);
      if (kstat < 0)
	goto error;
    }

  /*  Compute the distance between the points. */

  for (tlength=0.0, ki=1; ki<kneval; ki++)
    tlength += s6dist(sder+(ki-1)*kdim, sder+ki*kdim, kdim);

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */
  error:
  *jstat = kstat;
  s6err ("sh1762_sflength", *jstat, 0);
  goto out;

  out:
  return tlength;
}
  
#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9num (SISLObject * po, SISLObject * poref, int *jdiv, int *jstat)
#else
static void
sh1762_s9num (po, poref, jdiv, jstat)
     SISLObject *po;
     SISLObject *poref;
     int *jdiv;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find number of possible subdivisions of an object.
*
*
*
* INPUT      : po     - SISLObject to subdevide.
*              poref  - The other object in intersection.
*
*
* OUTPUT     : jdiv   - Possible subdivisions of object.
*                         = 0     : No subdivision.
*                         = 1     : Subdivision in first parameter direction.
*                         = 2     : Subdivision in second parameter direction.
*                         = 3     : Subdivision in both parameter direction.
*              jstat  - status messages
*                         = 0     : no error.
*                         < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*
*********************************************************************
*/
{
  int kstat = 0;
  int kgtpi1=0, kgtpi2=0;
  double tang1=DZERO, tang2=DZERO;
  int not_case_2d;
  int kbez1=1, kbez2=1;

  /* Init. */

  *jdiv = 0;

  if (po->iobj < SISLPOINT || po->iobj > SISLSURFACE)
    goto err121;
  if (poref->iobj < SISLPOINT || poref->iobj > SISLSURFACE)
    goto err121;

  if (po->iobj == SISLPOINT)
    goto out;

  kgtpi1 = 10;
  tang1 = HUGE;

  kgtpi2 = 0;
  tang2 = (double) 0.0;  /* VSK. 030394. Changed tang1 into tang2. */

  /* Get attributes from object to divide. */
  if (po->iobj == SISLCURVE)
    {
      if (po->c1->pdir != SISL_NULL)
	{
	  kgtpi1 = po->c1->pdir->igtpi;
	  tang1 = po->c1->pdir->aang;
	}
      kbez1 = (po->c1->ik == po->c1->in);
    }
  else
    {
      if (po->s1->pdir != SISL_NULL)
	{
	  kgtpi1 = po->s1->pdir->igtpi;
	  tang1 = po->s1->pdir->aang;
	}
      kbez1 = (po->s1->ik1 == po->s1->in1 && po->s1->ik2 == po->s1->in2);
    }

  /* Get attributes from referance object. */
  if (poref->iobj == SISLCURVE)
    {
      if (poref->c1->pdir != SISL_NULL)
	{
	  kgtpi2 = poref->c1->pdir->igtpi;
	  tang2 = poref->c1->pdir->aang;
	}
      kbez2 = (poref->c1->ik == poref->c1->in);
    }
  else if (poref->iobj == SISLSURFACE)
    {
      if (poref->s1->pdir != SISL_NULL)
	{
	  kgtpi2 = poref->s1->pdir->igtpi;
	  tang2 = poref->s1->pdir->aang;
	}
      kbez2 = (poref->s1->ik1 == poref->s1->in1 &&
	       poref->s1->ik2 == poref->s1->in2);
    }

    if (poref->iobj == SISLPOINT && poref->p1->idim == 2)
       not_case_2d = FALSE;
    else
       not_case_2d = TRUE;


    /* Test for number of division directions.     */
  /*---------------------------------------------*/
  /* If linear, we do not subdivide.             */
  if (kgtpi1 == 0 && tang1 <= ANGULAR_TOLERANCE/10.0 && not_case_2d)
    *jdiv = 0;

  else if (po->iobj == SISLCURVE && poref->iobj == SISLSURFACE)
    /* Subdivide curve. */
    {
      if (s1791 (po->c1->et, po->c1->ik, po->c1->in))
	*jdiv = 1;

      else
	*jdiv = 0;

    }

  else if (kgtpi1 == 0 && tang1 < SIMPLECASE / (double) 2.0 && kbez1 == 1 &&
	   (kgtpi2 != 0 || tang2 > tang1 * (double) 2.0))
    *jdiv = 0; 

  else if (po->iobj == SISLCURVE)
    {
      if (s1791 (po->c1->et, po->c1->ik, po->c1->in))
	*jdiv = 1;

      else
	*jdiv = 0;
    }
  else if (po->iobj == SISLSURFACE)
    {
	double tsfp1, tsfp2, tref;
	tref = 5.0;

	tsfp1 = sh1762_sflength(po->s1, 1, &kstat);
	if (kstat < 0)
	  goto error;

	tsfp2 = sh1762_sflength(po->s1, 2, &kstat);
	if (kstat < 0)
	  goto error;

	if (s1791 (po->s1->et1, po->s1->ik1, po->s1->in1)  &&
	  !(po->s1->ik1 == 2 && tsfp1 < tref*tsfp2))
	*jdiv = 1;

      else
	*jdiv = 0;

	if (s1791 (po->s1->et2, po->s1->ik2, po->s1->in2) &&
	  !(po->s1->ik2 == 2 && tsfp2 < tref*tsfp1))
	*jdiv += 2;

    }
  goto out;


  /* Error in lower level routine. */
  error:
  *jstat = kstat;
  s6err ("sh1762_s9num", *jstat, 0);
  goto out;

  /* Error. Kind of object does not exist.  */
err121:
  *jstat = -121;
  s6err ("sh1762_s9num", *jstat, 0);

out:;
}

#if defined(SISLNEEDPROTOTYPES)
static int 
sh1762_is_taboo(SISLSurf *psurf1, SISLSurf *psurf2, SISLIntpt *pintpt, 
		int idir, int *jstat)
#else
static int
sh1762_is_taboo(psurf1, psurf2 ,pintpt, idir, jstat)
    SISLSurf *psurf1;
    SISLSurf *psurf2;
    SISLIntpt *pintpt;
    int idir;
    int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if the intersection curve passing through
*	       the point is always parallel to an iso-curve. 
*
*
*
* INPUT      : psurf1   - 1. surface in intersection problem.
*              psurf2   - 2. surface in intersection problem or SISL_NULL.
*              pintpt   - Intersection point.
*              idir     - Parameter direction in surface.
*
*
* OUTPUT     : return value : 1 = is taboo point, 0 = no taboo point
*              jstat    - Status messages
*                          = 0     : OK
*                          < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 99-05.
*
*********************************************************************
*/
{
   static double parallel    = 0.01;
   static double fuzzy_angle = 1e-4;
   static double tol = (double) 1000000.0 * REL_COMP_RES;

   int kstat = 0;
   int is_taboo = 0;
   double derivs1[9], derivs2[9], norm[3], nor1[3], nor2[3], angle;
   double abs_tang1[2], abs_tang2[2];
   double tmax;
   int ilfs = 0, ilft = 0;

   if (psurf1->idim == 2)
     return 0;

   /* Test input. */

   if (psurf2 && (psurf1->idim != psurf2->idim || psurf1->idim != 3))
     goto err104;

   if (!psurf2 && psurf1->idim != 1)
     goto err105;

   if (psurf2)
     {
       /* Evaluate the intersection point in both surfaces. */

       s1421(psurf1, 1, &pintpt->epar[0], &ilfs, &ilft, derivs1, norm, &kstat);
       if (kstat < 0)
	 goto error;

       s1421(psurf2, 1, &pintpt->epar[2], &ilfs, &ilft, derivs2, norm, jstat);
       if (kstat < 0)
	 goto error;

       s6crss(derivs2+3, derivs2+6, nor2);
       s6crss(derivs1+3, derivs1+6, nor1);

       /* If we have a singularity, we don't declare it as taboo. */

       angle = s6ang(nor1, nor2, 3);

       abs_tang1[0] = fabs(s6scpr(derivs1+6, nor2, 3));
       abs_tang1[1] = fabs(s6scpr(derivs1+3, nor2, 3));

       abs_tang2[0] = fabs(s6scpr(nor1, derivs2+6, 3));
       abs_tang2[1] = fabs(s6scpr(nor1, derivs2+3, 3));

       if (angle < fuzzy_angle)
	 is_taboo = 0;
       else if (idir == 1 && abs_tang1[0] < parallel*abs_tang1[1])
	 is_taboo = 1;
       else if (idir == 2 && abs_tang1[1] < parallel*abs_tang1[0])
	 is_taboo = 1;
       else 
	 is_taboo = 0;
     }
   else 
     {
       /* Evaluate the intersection point. */

       s1421(psurf1, 1, &pintpt->epar[0], &ilfs, &ilft, derivs1, norm, &kstat);
       if (kstat < 0)
	 goto error;

       /* If we have a singularity, we don't declare it as taboo. */

       tmax = sqrt(derivs1[1]*derivs1[1] + derivs1[2]*derivs1[2]);
       if (tmax < tol)
	  /* The length of the surface normal is less than the 
	     given tolerance*/
	is_taboo = 0;

       else if (idir == 1 && fabs(derivs1[2]) < parallel*tmax)
	 is_taboo = 1;
       else if (idir == 2 && fabs(derivs1[1]) < parallel*tmax)
	 is_taboo = 1;
       else 
	 is_taboo = 0;
     }

   *jstat = 0;
   goto out;

   /* Error in lower order routine. */
  error:
  *jstat = kstat;
  s6err ("sh1762_is_taboo", *jstat, 0);
   goto out;

  /* Error. Dimension not equal to 3.  */
err104:
  *jstat = -104;
  s6err ("sh1762_is_taboo", *jstat, 0);
   goto out;

  /* Error. Conflicting dimensions.  */
err105:
  *jstat = -105;
  s6err ("sh1762_is_taboo", *jstat, 0);
   goto out;

out:
   return is_taboo;
}



#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9subdivpt (SISLObject * po1, SISLObject * po2, double aepsge,
		   int iobj, int idiv, SISLEdge * vedge[], SISLIntdat ** pintdat,
		   int *fixflag, SISLIntpt ** rpt, double epar[], int *jstat)
#else
static void
sh1762_s9subdivpt (po1, po2, aepsge, iobj, idiv, vedge, pintdat, fixflag, rpt, epar, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     int iobj;
     int idiv;
     SISLEdge *vedge[];
     SISLIntdat **pintdat;
     int *fixflag;
     SISLIntpt **rpt;
     double epar[];
     int *jstat;
#endif
 /* UPDATE ujk: Set up strategy for finding subdiv. point.
    Use of shsing(singular iteration) pluss segmentation. */
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find an appropriate subdivision point of an object when
*              the other object in the subdivision and eventual
*              intersections found on the edges, are known.
*
*
*
* INPUT      : po1      - First object in intersection.
*              po2      - Second object in intersection.
*              aepsge   - Geometry resolution.
*              iobj     - Number of object to divide.
*              vedge[]  - Intersection on edges.
*              pintdat  - Intersection data.
*
*
* OUTPUT     : fixflag  - Indicates if the subdivision point is fixed.
*              rpt      - An existing internal intersection point used as
*                         subdividing point. If no such point is found,
*                         *rpt == SISL_NULL.
*              epar     - Parameter values of subdivision point. The number of
*                         elements used is equal to the number of parameter
*                         directions of the object to be subdivide. The dimension
*                         is equal to 2.
*              jstat    - Status messages
*                          = 1     : Intersection found
*                          = 0     : no intersection found.
*                          < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
* REWISED BY : Vibeke Skytt, SI, 90-10.
* REWRITTEN BY : Vibeke Skytt, SINTEF, 94-02.
*
*********************************************************************
*/
{
   int kstat = 0;
   int kpos = 0;
   int kpar;            /* First index of subdivision point in the
			   parameter value of in intersection point.        */
   int kfound;          /* Indicates if in intersection point / extremal
			   point is to be used.                             */
   int kf1=0, kf2=0;    /* Indicates if an internal intersection point is
			   legal in the parameter directions of a surface.  */
   double tdel;		/* Parameter used to measure closeness to an edge.  */
   double tdel1, tdel2;	/* Parameters used to measure closeness to an edge. */
   double tstart, tend; /* Endparameters of curve.                          */
   double tstart2, tend2; /* Endparameters of second curve.                 */
   double tpar;         /* Parameter value of subdivision point. */
   double tpar2;        /* Parameter value of point from iteration. */
   double sstart[2], send[2];  /* Endparameters of surface.      */
   double spar[2];      /* Parameter value of subdivision point. */
   double spar2[2];     /* Parameter value of subdivision point. */
   double sparsave[2];  /* Parameter value of subdivision point. */
   SISLObject *qo1;	/* Pointer to the object that is to be subdivided. */
   SISLObject *qo2;	/* Pointer to the other object.          */
   SISLIntpt *qpt = SISL_NULL;  /* An internal intersection point.    */

   /* Set pointer to subdivision object. */

   qo1 = (iobj == 1 ? po1 : po2);
   qo2 = (iobj == 1 ? po2 : po1);
   kpar = (iobj == 1 ? 0 : po1->iobj);

  *jstat = 0;

  /* Branch on subdivision object. */

  if (qo1->iobj == SISLCURVE)
  {
     /* Find a proper subdivision value of the curve.  First set when a point
	is to close to an edge to be used as a subdivision point. */

     tdel = (double) 0.01 *(qo1->c1->et[qo1->c1->in] -
			    qo1->c1->et[qo1->c1->ik - 1]);

     /* Try to find an internal intersection point. */

     s6idint (po1, po2, *pintdat, &qpt, iobj);
     if (!(3*qo1->c1->ik > qo1->c1->in) &&
	 /* if (qo1->c1->ik != qo1->c1->in && */
	 (sh6ismain (qpt)) && sh6nmbhelp (qpt,&kstat) == 0)
	qpt = SISL_NULL;

      if (qpt != SISL_NULL)
	{
	  /* Internal intersection point found. */

	  tpar = qpt->epar[kpar];

	  if (tpar < (qo1->c1->et[qo1->c1->ik - 1] + tdel) ||
	      tpar > (qo1->c1->et[qo1->c1->in] -tdel))
	       qpt = SISL_NULL;  /* Do not use the point as a subdivision point. */
	}

      if (qpt == SISL_NULL &&
	  vedge[iobj - 1]->ipoint == 0 && qo1->c1->ik == qo1->c1->in)
      {
	 /* No internal intersection is found. The curve is of Bezier type,
	    and there is no intersection on the endpoints of the curve.
	    Then we try to iterate in order to find an intersection or
	    closest point to use as a subdivision point. Branch on the
	    various kind of other objects involved in the intersection. */

	 tstart = qo1->c1->et[qo1->c1->ik - 1];
	 tend = qo1->c1->et[qo1->c1->in];
	 tpar = (tstart + tend) * (double) 0.5;
	 kfound = 1;

	 if (qo2->iobj == SISLPOINT)
	 {
	    /* ALA & UJK start 31/10/90. */
	    if (qo2->p1->idim == 1)
	       s1172 (qo1->o1->c1, tstart, tend,
		      tpar, &tpar, &kstat);
	    else
	    {
	       kstat = 1;   /* Use quick iteration. */
	       s1771 (qo2->o1->p1, qo1->o1->c1, aepsge, tstart, tend,
		      tpar, &tpar, &kstat);
	    }
	    if (kstat < 0)
	       goto error;
	 }

	 else if (qo2->iobj == SISLCURVE)
	 {
	    tstart2 = qo2->c1->et[qo2->c1->ik - 1];
	    tend2 = qo2->c1->et[qo2->c1->in];
	    tpar2 = (tstart + tend) * (double) 0.5;
	    tdel2 = (double)0.01*(tend - tstart);

	    s1770 (qo1->o1->c1, qo2->o1->c1, aepsge, tstart, tstart2, tend,
		   tend2, tpar, tpar2, &tpar, &tpar2, &kstat);
	    if (kstat < 0)
	       goto error;

	    /* Test the subdivision point towards the endpoint of the
	       second curve. */

	    if (tpar2 < tstart2+tdel2 || tpar2 > tend2-tdel2)
	       kfound = 0;
	 }

	 else if (qo2->iobj == SISLSURFACE)
	 {
	    sstart[0] = qo2->s1->et1[qo2->s1->ik1 - 1];
	    sstart[1] = qo2->s1->et2[qo2->s1->ik2 - 1];

	    send[0] = qo2->s1->et1[qo2->s1->in1];
	    send[1] = qo2->s1->et2[qo2->s1->in2];

	    spar[0] = (sstart[0] + send[0]) * (double) 0.5;
	    spar[1] = (sstart[1] + send[1]) * (double) 0.5;

	    tdel1 = (double)0.01* (send[0] - sstart[0]);
	    tdel2 = (double)0.01* (send[1] - sstart[1]);

	    kstat = 1;    /* Use quick iteration. */
	    s1772 (qo1->o1->c1, qo2->o1->s1, aepsge, tstart, sstart, tend,
		   send, tpar, spar, &tpar, spar, &kstat);
	    if (kstat < 0)
	       goto error;

	    /* Test the subdivision point towards the edges of the surface. */

	    if (spar[0] < sstart[0]+tdel1 || spar[0] > send[0]-tdel1 ||
		spar[1] < sstart[1]+tdel2 || spar[1] > send[1]-tdel2)
	       kfound = 0;
	 }

	 /* Test the subdivision point towards the edges of the subdivision
	    curve. */

	 if (!kfound ||
	     tpar < tstart+tdel || tpar > tend-tdel)

	    /* Use the midpoint of the curve as subdivision point. */

	    tpar = s1792 (qo1->c1->et, qo1->c1->ik, qo1->c1->in);
      }
      else if (qpt == SISL_NULL)
	 /* Use the midpoint as a subdivision point. */

	 tpar = s1792 (qo1->c1->et, qo1->c1->ik, qo1->c1->in);

      /* Set output variables  */

      epar[0] = tpar;
      *rpt = qpt;
  }
  else if (qo1->iobj == SISLSURFACE)
  {
     /* Find a subdivision point of the surface. Branch on the other
	object involved in the intersection. First set the endparameters
	of the surface and when a point is to close to an edge
	to be used as a subdivision point. */

     sstart[0] = qo1->s1->et1[qo1->s1->ik1 - 1];
     sstart[1] = qo1->s1->et2[qo1->s1->ik2 - 1];

     send[0] = qo1->s1->et1[qo1->s1->in1];
     send[1] = qo1->s1->et2[qo1->s1->in2];

     tdel1 = (double) 0.01 *(send[0] - sstart[0]);
     tdel2 = (double) 0.01 *(send[1] - sstart[1]);

     /* In the Bezier case, search for an internal intersection point. */

     if (qo1->s1->ik1 == qo1->s1->in1 && qo1->s1->ik2 == qo1->s1->in2)
	s6idint (po1, po2, *pintdat, &qpt, iobj);
     if (qpt != SISL_NULL)
     {
	/* Internal intersection point found. */
	sparsave[0] = spar[0] = qpt->epar[kpar];
	sparsave[1] = spar[1] = qpt->epar[kpar + 1];
	kf1 = kf2 = 1;

	/* Test the point towards the edges of the surface. */

	if (spar[0] < sstart[0] + tdel1 || spar[0] > send[0] - tdel1)
	{
	   kf1--;
	   qpt = SISL_NULL;
	}
	if (spar[1] < sstart[1] + tdel2 || spar[1] > send[1] - tdel2)
	{
	   kf2--;
	   qpt = SISL_NULL;
	}
     }

     kfound = 0;   /* If no iteration is tryed, use the midpoint. */
     if ((!qpt) && qo2->iobj != SISLSURFACE &&
	 !(qo2->iobj == SISLPOINT && qo2->p1->idim == 1) &&
	 qo1->s1->ik1 == qo1->s1->in1 && qo1->s1->ik2 == qo1->s1->in2)
     {
	/* No internal intersection is found. The second object is not a
	   surface, and the subdivision surface is of Bezier type.
	   Prepare for iteration. */

	spar[0] = (sstart[0] + send[0]) * (double) 0.5;
	spar[1] = (sstart[1] + send[1]) * (double) 0.5;
	kfound = 3;

	if (qo2->iobj == SISLPOINT)
	{
	   s1773 (qo2->o1->p1, qo1->o1->s1, aepsge, sstart, send, spar,
		  spar, &kstat);
	   if (kstat < 0)
	      goto error;
	}

	else if (qo2->iobj == SISLCURVE)
	{
	   tstart = qo2->c1->et[qo2->c1->ik - 1];
	   tend = qo2->c1->et[qo2->c1->in];
	   tpar = (tstart + tend) * (double) 0.5;
	   tdel = (double)0.01*(tend - tstart);

	   kstat = 1;
	   s1772 (qo2->o1->c1, qo1->o1->s1, aepsge, tstart, sstart,
		  tend, send, tpar, spar, &tpar, spar, &kstat);
	   if (kstat < 0)
	      goto error;

	   /* Control the edges of the curve. */

	   if (tpar < tstart+tdel || tpar > tend-tdel)
	      kfound = 0;
	}

	/* Test the edges of the surface to be subdivided. */

	   if (spar[0] < sstart[0]+tdel1 || spar[0] > send[0]-tdel1)
	      kfound--;
	   if (spar[1] < sstart[1]+tdel2 || spar[1] > send[1]-tdel2)
	      kfound -= 2;
	}

     if ((!qpt) && (!(kfound==3) && qo2->iobj != SISLSURFACE &&
		    !(qo2->iobj == SISLPOINT && qo2->p1->idim == 1)))
	 {
	    /* Use the midpoint of the surface as a subdivision point. */

	    if (kfound != 1)
	       spar[0] = s1792 (qo1->s1->et1, qo1->s1->ik1, qo1->s1->in1);
	    if (kfound != 2)
	       spar[1] = s1792 (qo1->s1->et2, qo1->s1->ik2, qo1->s1->in2);

	    /* Test if this subdivision point is too close to an existing
	       inner intersection point. */

	    if (kf1 && fabs(spar[0]-sparsave[0]) < tdel1)
	       spar[0] = sparsave[0];
	    if (kf2 && fabs(spar[1]-sparsave[1]) < tdel2)
	       spar[1] = sparsave[1];
	 }

     if ((!qpt) && (qo2->iobj == SISLSURFACE ||
		    (qo2->iobj == SISLPOINT &&
		     (qo2->p1->idim == 1 || qo2->p1->idim == 2))))
     {
	SISLPtedge *qptedg;	/* Pointer used to traverse int. points on edges. */
	SISLIntpt *pt1 = SISL_NULL;  /* Intersection point on edge. */
	SISLIntpt *pt2 = SISL_NULL;  /* Intersection point on edge. */
	SISLIntpt *ptsing1 = SISL_NULL; /* Singular intersection point on edge. */
	SISLIntpt *ptsing2 = SISL_NULL; /* Singular intersection point on edge. */
	SISLIntpt *pcurr;          /* Current intersection point.          */
	int kj;                    /* Counter.                             */
	double tmean[2];           /* Middle parameter of the surface.     */
	double tpar1=HUGE, tpar2=HUGE;  /* Used for comparisement with
					   intersection point.             */
	int ktype1=-10, ktype2=-10;     /* As previous.                    */

	/* There is a surface-surface intersection or an intersection
	   between a surface and a point in 1D. In both cases intersection
	   curves are the expected output. Start by logging the intersection
	   points at the edges. */
	/* If the surface is almost a Bezier surface, make it Bezier. */

	s9simple_knot(qo1->s1, idiv, spar, fixflag, &kstat);
	if ( kstat < 0 ) goto error;

	memcopy(sparsave, spar, 2, DOUBLE);
	if (((*fixflag) == 1 || (*fixflag) == 3) &&
	    (spar[0] < sstart[0]+tdel1 || spar[0] > send[0]-tdel1))
	   *fixflag -= 1;
	if (((*fixflag) == 2 || (*fixflag) == 3) &&
	    (spar[1] < sstart[1]+tdel2 || spar[1] > send[1]-tdel2))
	   *fixflag -= 2;

	if ( *fixflag < 3 )
	{
	   /* In at least one parameter direction there is a freedom
	      of the subdivision point.                               */

	   /* Set the middle parameter.  */

	   tmean[0] = s1792 (qo1->s1->et1, qo1->s1->ik1, qo1->s1->in1);
	   tmean[1] = s1792 (qo1->s1->et2, qo1->s1->ik2, qo1->s1->in2);

	   if (!(*fixflag == 1) && vedge[iobj - 1]->ipoint > 0)
	   {
	      /* Search for intersection points on the edges in the
		 first paramter direction, i.e. edge 1 and 3. Find the
		 intersection point closest to the middle parameter value
		 and distinguish between ordinary intersection points and
		 singular or almost singular (touchy) points. */

	      /* Loop for edges no 1 and 3*/
	      for (kj = 0; kj < 3; kj += 2)
		 /* Loop for all points on edge*/
		 for (qptedg = vedge[iobj - 1]->prpt[kj]; qptedg != SISL_NULL;
	       qptedg = qptedg->pnext)
		 {
		    pcurr = qptedg->ppt;

		    /* Test if the point is too close to an edge. */

		    if (pcurr->epar[kpar] < sstart[0]+tdel1 ||
			pcurr->epar[kpar] > send[0]-tdel1) continue;

		    /* Check if the intersection curve passing through
		       the point is always parallel to an iso-curve. */
		    
		    if (sh1762_is_taboo(qo1->s1, 
					(qo2->iobj == SISLSURFACE) ? 
					qo2->s1 : SISL_NULL, 
					pcurr, 1, &kstat))
		      continue;

		    if (kstat < 0)
		      goto error;

		    if (pcurr->iinter == SI_SING)
		    {
		       /* Test if the singular/near singular point is the one
			  closest to the middle point.  */

		       if (!ptsing1 || fabs(pcurr->epar[kpar]-tmean[0]) <
			   fabs(ptsing1->epar[kpar]-tmean[0]))
			  ptsing1 = pcurr;
		    }
		    else
		    {
		       /* Test if the intersection point is the one closest
			  to the middle. */

		       if (!pt1 || fabs(pcurr->epar[kpar]-tmean[0]) <
			   fabs(pt1->epar[kpar]-tmean[0]))
			  pt1 = pcurr;
		    }
		 }
	   }

	   if (!(*fixflag == 2) && vedge[iobj - 1]->ipoint > 0)
	   {
	      /* Search for intersection points on the edges in the
		 second paramter direction, i.e. edge 2 and 4. Find the
		 intersection point closest to the middle parameter value
		 and distinguish between ordinary intersection points and
		 singular or almost singular (touchy) points. */

	      /* Loop for edges no 2 and 4*/
	      for (kj = 1; kj < 4; kj += 2)
		 /* Loop for all points on edge*/
		 for (qptedg = vedge[iobj - 1]->prpt[kj]; qptedg != SISL_NULL;
	       qptedg = qptedg->pnext)
		 {
		    pcurr = qptedg->ppt;

		    /* Test if the point is too close to an edge. */

		    if (pcurr->epar[kpar+1] < sstart[1]+tdel2 ||
			pcurr->epar[kpar+1] > send[1]-tdel2) continue;

		    /* Check if the intersection curve passing through
		       the point is always parallel to an iso-curve. */
		    
		    if (sh1762_is_taboo(qo1->s1,  
					(qo2->iobj == SISLSURFACE) ? 
					qo2->s1 : SISL_NULL, 
					pcurr, 2, &kstat))
		      continue;

		    if (kstat < 0)
		      goto error;

		    if (pcurr->iinter == SI_SING)
		    {
		       /* Test if the singular/near singular point is the one
			  closest to the middle point.  */

		       if (!ptsing2 || fabs(pcurr->epar[kpar+1]-tmean[1]) <
			   fabs(ptsing2->epar[kpar+1]-tmean[1]))
			  ptsing2 = pcurr;
		    }
		    else
		    {
		       /* Test if the intersection point is the one closest
			  to the middle. */

		       if (!pt2 || fabs(pcurr->epar[kpar+1]-tmean[1]) <
			   fabs(pt2->epar[kpar+1]-tmean[1]))
			  pt2 = pcurr;
		    }
		 }
	   }

	   if (qo1->s1->idim == 1)
	   {
	      /* One-dimensional case. Iterate to find an extremal point. */

	      if (!(*fixflag == 1) && ptsing1)
	      {
		 /* Set startpoint to iteration. */

		 spar[0] = ptsing1->epar[kpar];
		 spar[1] = ptsing1->epar[kpar+1];
	      }
	      else if (!(*fixflag == 2) && ptsing2)
	      {
		 spar[0] = ptsing2->epar[kpar];
		 spar[1] = ptsing2->epar[kpar+1];
	      }
	      else
	      {
		 /* No (almost) singular intersection point is found
		    at the edge. */

		 spar[0] = (double)0.5*(sstart[0] + send[0]);
		 spar[1] = (double)0.5*(sstart[1] + send[1]);
	      }

	      /* Perform iteration. */

	      kfound = 0;
	      s1174 (qo1->o1->s1, sstart, send, spar, spar, &kstat);
	      if (kstat < 0)
		goto error;
	      if (kstat == 1)
		{
		   /* An extremal point is found. Test if it is too close
		      to an edge. */

		   kfound = 3;
		   if (spar[0] < sstart[0]+tdel1 || spar[0] > send[0]-tdel1)
		      kfound--;
		   if (spar[1] < sstart[1]+tdel2 || spar[1] > send[1]-tdel2)
		      kfound -= 2;
		}

	      if (*fixflag == 0 && ptsing2 && ptsing1)
	      {
		 /* Try a second iteration for an extremal point
		    in order to find a subdivision parameter in the
		    second parameter direction.  */

		 spar2[0] = ptsing2->epar[kpar];
		 spar2[1] = ptsing2->epar[kpar+1];

		 s1174(qo1->o1->s1, sstart, send, spar2, spar2, &kstat);
		 if (kstat < 0)
		    goto error;
		 if (kstat == 1)
		 {
		    /* An extremal point is found. Test it against the edges. */

		    if (!(spar2[1] < sstart[1]+tdel2 || spar2[1] > send[1]-tdel2))
		    {
		       spar[1] = spar2[1];
		       if (kfound < 2) kfound += 2;
		    }
		 }
	      }

	      /* Set intermediate subdivision point. */

	      if (*fixflag == 1)
		 spar[0] = sparsave[0];
	      else if (kfound == 1 || kfound == 3)
		 (*fixflag)++;
	      else if (ptsing1)
	      {
		 spar[0] = ptsing1->epar[kpar];
		 (*fixflag)++;
	      }
	      else if (pt1)
	      {
		 spar[0] = pt1->epar[kpar];
		 (*fixflag)++;
	      }
	      else
		 spar[0] = tmean[0];

	      if (*fixflag == 2)
		 spar[1] = sparsave[1];
	      else if (kfound == 2 || kfound == 3)
		 (*fixflag) += 2;
	      else if (ptsing2)
	      {
		 spar[1] = ptsing2->epar[kpar+1];
		 (*fixflag) += 2;
	      }
	      else if (pt2)
	      {
		 spar[1] = pt2->epar[kpar+1];
		 (*fixflag) += 2;
	      }
	      else
		 spar[1] = tmean[1];
	   }
	   else
	   {
	      /* Surface-surface intersection. Set intermediate
		 subdivision point. */

	      if (*fixflag == 1)
		 spar[0] = sparsave[0];
	      else if (ptsing1 && pt1)
	      {
		 if (fabs(ptsing1->epar[kpar]-tmean[0]) <=
		     fabs(pt1->epar[kpar]-tmean[0]))
		    spar[0] = ptsing1->epar[kpar];
		 else
		    spar[0] = pt1->epar[kpar];
		 (*fixflag)++;
	      }
	      else if (ptsing1)
	      {
		 spar[0] = ptsing1->epar[kpar];
		 (*fixflag)++;
	      }
	      else if (pt1)
	      {
		 spar[0] = pt1->epar[kpar];
		 (*fixflag)++;
	      }
	      else spar[0] = tmean[0];

	      if (*fixflag == 2)
		 spar[1] = sparsave[1];
	      else if (ptsing2 && pt2)
	      {
		 if (fabs(ptsing2->epar[kpar+1]-tmean[1]) <=
		     fabs(pt2->epar[kpar+1]-tmean[1]))
		    spar[1] = ptsing2->epar[kpar+1];
		 else
		    spar[1] = pt2->epar[kpar+1];
		 (*fixflag) += 2;
	      }
	      else if (ptsing2)
	      {
		 spar[1] = ptsing2->epar[kpar+1];
		 (*fixflag) += 2;
	      }
	      else if (pt2)
	      {
		 spar[1] = pt2->epar[kpar+1];
		 (*fixflag) += 2;
	      }
	      else spar[1] = tmean[1];

	   }
	}

	/* Test if the found subdivision value lies very close to an
	   existing intersection point. In that case move the subdivision
	   point to the intersection point. The two parameter directions
	   are treated separately.  */

	if ((*pintdat) && (*pintdat)->ipoint > 0)
	   for (kj=0; kj<(*pintdat)->ipoint; kj++)
	   {
	      pcurr = (*pintdat)->vpoint[kj];

	      if ((*fixflag)==1 || (*fixflag)==3)
	      {
		 if (fabs(spar[0]-pcurr->epar[kpar]) < (double)0.001*tdel1)
		 {
		    if (fabs(spar[0]-pcurr->epar[kpar]) < fabs(tpar1-spar[0]) &&
			ktype1 <= pcurr->iinter)
		    {
		       tpar1 = pcurr->epar[kpar];
		       ktype1 = pcurr->iinter;
		    }
		 }
	      }
	      else
	      {
		 if (fabs(spar[0]-pcurr->epar[kpar]) < (double)0.1*tdel1 &&
		     pcurr->epar[kpar] >= sstart[0]+tdel1 &&
		     pcurr->epar[kpar] <= send[0]-tdel1)
		 {
		    if (fabs(spar[0]-pcurr->epar[kpar]) < fabs(tpar1-spar[0]) &&
			ktype1 <= pcurr->iinter)
		    {
		       tpar1 = pcurr->epar[kpar];
		       ktype1 = pcurr->iinter;
		    }
		 }
	      }

	      if ((*fixflag)==2 || (*fixflag)==3)
	      {
		 if (fabs(spar[1]-pcurr->epar[kpar+1]) < (double)0.001*tdel2)
		 {
		    if (fabs(spar[1]-pcurr->epar[kpar+1]) < fabs(tpar2-spar[1]) &&
			ktype2 <= pcurr->iinter)
		    {
		       tpar2 = pcurr->epar[kpar+1];
		       ktype2 = pcurr->iinter;
		    }
		 }
	      }
	      else
	      {
		 if (fabs(spar[1]-pcurr->epar[kpar+1]) < (double)0.1*tdel2 &&
		     pcurr->epar[kpar+1] >= sstart[1]+tdel2 &&
		     pcurr->epar[kpar+1] <= send[1]-tdel2)
		 {
		    if (fabs(spar[1]-pcurr->epar[kpar+1]) < fabs(tpar2-spar[1]) &&
			ktype2 <= pcurr->iinter)
		    {
		       tpar2 = pcurr->epar[kpar+1];
		       ktype2 = pcurr->iinter;
		    }
		 }
	      }
	   }

	if (ktype1 > -10 && tpar1 > sstart[0]+tdel1 && tpar1 < send[0]-tdel1)
	{
	   if (!((*fixflag == 1 || *fixflag == 3) && ktype1 < 0))
	      spar[0] = tpar1;
	}
	if (ktype2 > -10 && tpar2 > sstart[1]+tdel2 && tpar2 < send[1]-tdel2)
	{
	   if (!((*fixflag == 2 || *fixflag == 3) && ktype2 < 0))
	      spar[1] = tpar2;
	}
     }

      /* Set output variables.  */

     epar[0] = spar[0];
     epar[1] = spar[1];
     *rpt = qpt;
     *fixflag = ((*fixflag) >=2) ? 1 : 0;
  }
  else goto err122;  /* Unexpected kind of object. */

  goto out;

/* Error. Unexpected kind of object.  */

err122:*jstat = -122;
  s6err ("sh1762_s9subdivpt", *jstat, kpos);
  goto out;

/* Error in lower level routine.  */

error:*jstat = kstat;
  s6err ("sh1762_s9subdivpt", *jstat, kpos);
  goto out;

out:
   return;

}

#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9div (SISLObject * po1, SISLObject * po2, double aepsge,
	      int iobj, int idiv, SISLObject * wob[], SISLEdge * vedge[],
	      SISLIntdat ** pintdat, int *jstat)

#else
static void
sh1762_s9div (po1, po2, aepsge, iobj, idiv, wob, vedge, pintdat, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     int iobj;
     int idiv;
     SISLObject *wob[];
     SISLEdge *vedge[];
     SISLIntdat **pintdat;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Uppdating intersection data if possible.
*
*
*
* INPUT      : po1      - First object in intersection.
*              po2      - Second object in intersection.
*              aepsge   - Geometry resolution.
*              iobj     - Number of object to divide.
*              idiv     - Subdivision direction.
*                          = 0     : No subdivision.
*		           = 1     : Subdivision in first parameter
*                                                          direction.
*		           = 2     : Subdivision in second parameter
*                                                          direction.
*		           = 3     : Subdivision in first and second
*                                                parameter direction.
*              vedge[]  - Intersection on edges.
*
*
* OUTPUT     : pintdat  - Intersection data.
*              wob[]    - Pointers to the subdivided objects.
*              jstat    - Status messages
*                          = 1     : Intersection found
*                          = 0     : no intersection found.
*                          < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*
*********************************************************************
*/
{
  int kpos = 0;			/* Position of error.      */
  int kstat = 0;		/* Local status variable.           */
  int ki, kn;			/* Counters.                        */
  int kpar;			/* First parameter direction corresponding the the object
			           that is to be subdivided.        */
  int knum;			/* Total number of points in the data structures of the
			           original problem, and the problem of reduced dimension. */
  double tdel;			/* Parameter used to measure closeness between an
			           intersection point and the subdivision point.    */
  double spar[2];		/* Parameter values of subdividing point. If a curve is to
			           be subdivided, only the first element is used.           */
  SISLCurve *qcrv = SISL_NULL;       /* Mother curve of subdivision curve.               */
  SISLObject *qso = SISL_NULL;	/* Subdivision object, it is a point if a curve is subdivided
			           and a curve if a surface is subdivided.          */
  SISLObject *qmotherobj = SISL_NULL; /* Mother object of subdivision object.            */
  SISLObject *qo1 = SISL_NULL;	/* Pointer to the object that is to be subdivided.  */
  SISLObject *qo2 = SISL_NULL;	/* Pointer to the other object.                     */
  SISLObject *qs1 = SISL_NULL, *qs2 = SISL_NULL;	/* Subobjects to be used in the case where a surface
				       is to be subdivided in both parameter direction to
				       store intermediate subsurfaces.               */
  SISLIntpt *qpt = SISL_NULL;	/* An eventual found inner intersection used as subdivision
			           point. If the subdivision point is found in another way,
			           qpt = SISL_NULL.                                      */
  SISLIntdat *qintdat = SISL_NULL;	/* Data structure of intersection problem with lower dim. */
  SISLIntpt *qp;		/* Closest intersection point to the subdiv point */
  SISLEdge *uedge[2];		/* Edge intersections of subproblem.                 */
  int fixflag = 0;		/* UJK 31.10.90 */
  int idummy;

  /* Fetch subdivision point of object.  */

  sh1762_s9subdivpt (po1, po2, aepsge, iobj, idiv, vedge, pintdat, &fixflag, &qpt, spar, &kstat);
  if (kstat < 0)
    goto error;

  qo1 = (iobj == 1 ? po1 : po2);
  qo2 = (iobj == 1 ? po2 : po1);
  kpar = (iobj == 1 ? 0 : po1->iobj);

  *jstat = 0;


  if (qo1->iobj == SISLCURVE)
    {
      /* Subdivide the curve at the found subdivision parameter value.  */

      /* printf("Subdivide curve. Parameter value = %10.6f \n",spar[0]); */

      s1231 (qo1->c1, spar[0], &(wob[0]->c1), &(wob[1]->c1), &kstat);
      if (kstat < 0)
	goto error;


      if (wob[0]->edg[1] == SISL_NULL)
	{
	  if ((qso = wob[0]->edg[1] = newObject (SISLPOINT)) == SISL_NULL)
	    goto err101;

	  /* Pick out end point from a curve. */

	  s1438 (wob[0]->c1, 1, &(qso->p1), &spar[0], &kstat);
	  if (kstat < 0)
	    goto error;
	}
      else
	qso = wob[0]->edg[1];

      if (po1->iobj + po2->iobj > SISLCURVE)
	{
	  /***** Treating edges on sub problems. *****/

	  /* We first have to transform intersection points to the the new
	     intersection format qintdat. */

	  /* UJK, newi */
	  sh6idget (po1, po2, kpar, spar[0], *pintdat, &qintdat, aepsge, &kstat);


	  /* Making new edge object to sub problems. */

	  if ((iobj == 1 ? qso : po1)->iobj == SISLPOINT)
	    uedge[0] = SISL_NULL;
	  else if ((uedge[0] = newEdge (vedge[0]->iedge - (iobj == 1 ? 2 : 0))) == SISL_NULL)
	    goto err101;
	  if ((iobj == 2 ? qso : po2)->iobj == SISLPOINT)
	    uedge[1] = SISL_NULL;
	  else if ((uedge[1] = newEdge (vedge[1]->iedge - (iobj == 2 ? 2 : 0))) == SISL_NULL)
	    goto err101;

	  /* Update edge intersection on sub problems. */

	  sh6idalledg ((iobj == 1 ? qso : po1), (iobj == 1 ? po2 : qso), qintdat, uedge, &kstat);
	  if (kstat < 0)
	    goto error;

	  /* Examine if the subdividing point intersect the second object. */

	  qso->o1 = qso;

	  sh1762 ((iobj == 1 ? qso : po1), (iobj == 1 ? po2 : qso), aepsge,
		  &qintdat, uedge, &kstat);
	  if (kstat < 0)
	    goto error;

	  if (uedge[0] != SISL_NULL)
	    freeEdge (uedge[0]);
	  if (uedge[1] != SISL_NULL)
	    freeEdge (uedge[1]);
	}
      else
	{
	  sh1761 ((iobj == 1 ? qso : po1), (iobj == 1 ? po2 : qso), aepsge,
		  &qintdat, &kstat);
	  if (kstat < 0)
	    goto error;
	}

      if (kstat)
	{
	  /* Total number of points. */

	  knum = (*pintdat == SISL_NULL ? 0 : (*pintdat)->ipoint) + qintdat->ipoint;

	  *jstat = 1;		/* Mark intersection found. */

	  /* Intersection found and we have to register the intersection
	     points. */

	  /* UJK newi */
	  sh1782 (po1, po2, aepsge, qintdat, kpar, spar[0], pintdat, &idummy, &kstat);
	  if (kstat < 0)
	    goto error;

	  /* UJK newi divide curve */
	  /* UPDATE: ? what about help points from s1782 knum?? */
	  if (qpt != SISL_NULL && (*pintdat)->ipoint == knum)
	    {
	      /* Find the closest poin to qpt. */

	      s6idcpt (*pintdat, qpt, &qp);

	      /* UJK newi, unite the points : */
	      sh6idnewunite (po1, po2, pintdat, &qpt, &qp, (double) 0.5,
			     aepsge, &kstat);
	      if (kstat < 0)
		goto error;
	    }
	}

      if (qintdat != SISL_NULL)
	{
	  freeIntdat (qintdat);
	  qintdat = SISL_NULL;
	}
    }
  else if (qo1->iobj == SISLSURFACE)
    {

      /* Subdivide surface and treat subdivision curves.  */

      for (ki = 0; ki < (idiv < 3 ? 1 : 3); ki++)
	{
	  if (idiv == 1)
	    {
	      /* printf("Subdivide surface. 1. par dir. par = %10.6f \n",spar[0]); */

	      s1711 (qo1->s1, 1, spar[0], &(wob[0]->s1), &(wob[1]->s1), &kstat);
	      if (kstat < 0)
		goto error;

	      if (wob[0]->edg[1] == SISL_NULL)
		{
		  if ((qso = wob[0]->edg[1] = newObject (SISLCURVE)) == SISL_NULL)
		    goto err101;

		  /* Pick out edge curve from a surface. */

		  s1435 (wob[0]->s1, 1, &(qso->c1), spar, &kstat);
		  if (kstat < 0)
		    goto error;
		}
	      else
		qso = wob[0]->edg[1];

	      /* Pick curve from mother object of surface, and make
		 motherobject of curve.                             */

	      s1437(qo1->o1->s1,spar[0],&qcrv,&kstat);
	      if (kstat < 0) goto error;

	      if ((qmotherobj = newObject(SISLCURVE)) == SISL_NULL) goto err101;
	      qmotherobj->c1 = qcrv;
	      qso->o1 = qmotherobj;
	    }
	  else if (idiv == 2)
	    {
	      /* printf("Subdivide surface. 2. par dir. par = %10.6f \n",spar[1]); */

	      s1711 (qo1->s1, 2, spar[1], &(wob[0]->s1), &(wob[1]->s1), &kstat);
	      if (kstat < 0)
		goto error;

	      if (wob[0]->edg[2] == SISL_NULL)
		{
		  if ((qso = wob[0]->edg[2] = newObject (SISLCURVE)) == SISL_NULL)
		    goto err101;

		  /* Pick out edge curve from a surface. */

		  s1435 (wob[0]->s1, 2, &(qso->c1), spar + 1, &kstat);
		  if (kstat < 0)
		    goto error;
		}
	      else
		qso = wob[0]->edg[2];

	      /* Pick curve from mother object of surface, and make
		 motherobject of curve.                             */

	      s1436(qo1->o1->s1,spar[1],&qcrv,&kstat);
	      if (kstat < 0) goto error;

	      if ((qmotherobj = newObject(SISLCURVE)) == SISL_NULL) goto err101;
	      qmotherobj->c1 = qcrv;
	      qso->o1 = qmotherobj;
	    }
	  else if (ki == 0)
	    {
	      if ((qs1 = newObject (SISLSURFACE)) == SISL_NULL)
		goto err101;
	      if ((qs2 = newObject (SISLSURFACE)) == SISL_NULL)
		goto err101;

	      /* printf("Subdivide surface. 1. par dir. par = %10.6f \n",spar[0]); */

	      s1711 (qo1->s1, 1, spar[0], &(qs1->s1), &(qs2->s1), &kstat);
	      if (kstat < 0)
		goto error;

	      if (qs1->edg[1] == SISL_NULL)
		{
		  if ((qso = qs1->edg[1] = newObject (SISLCURVE)) == SISL_NULL)
		    goto err101;

		  /* Pick out edge curve from a surface. */

		  s1435 (qs1->s1, 1, &(qso->c1), spar, &kstat);
		  if (kstat < 0)
		    goto error;
		}
	      else
		qso = qs1->edg[1];

	      /* Pick curve from mother object of surface, and make
		 motherobject of curve.                             */

	      s1437(qo1->o1->s1,spar[0],&qcrv,&kstat);
	      if (kstat < 0) goto error;

	      if ((qmotherobj = newObject(SISLCURVE)) == SISL_NULL) goto err101;
	      qmotherobj->c1 = qcrv;
	      qso->o1 = qmotherobj;
	    }
	  else if (ki == 1)
	    {
	      /* printf("Subdivide surface. 2. par dir. par = %10.6f \n",spar[1]); */

	      s1711 (qs1->s1, 2, spar[1], &(wob[0]->s1), &(wob[1]->s1), &kstat);
	      if (kstat < 0)
		goto error;

	      if (wob[0]->edg[2] == SISL_NULL)
		{
		  if ((qso = wob[0]->edg[2] = newObject (SISLCURVE)) == SISL_NULL)
		    goto err101;

		  /* Pick out edge curve from a surface. */

		  s1435 (wob[0]->s1, 2, &(qso->c1), spar + 1, &kstat);
		  if (kstat < 0)
		    goto error;
		}
	      else
		qso = wob[0]->edg[2];

	      /* Pick curve from mother object of surface, and make
		 motherobject of curve.                             */

	      s1436(qo1->o1->s1,spar[1],&qcrv,&kstat);
	      if (kstat < 0) goto error;

	      if ((qmotherobj = newObject(SISLCURVE)) == SISL_NULL) goto err101;
	      qmotherobj->c1 = qcrv;
	      qso->o1 = qmotherobj;
	    }
	  else
	    /* if (ki == 2) */
	    {
	      /* printf("Subdivide surface. 2. par dir. par = %10.6f \n",spar[1]); */

	      s1711 (qs2->s1, 2, spar[1], &(wob[2]->s1), &(wob[3]->s1), &kstat);
	      if (kstat < 0)
		goto error;

	      if (wob[2]->edg[2] == SISL_NULL)
		{
		  if ((qso = wob[2]->edg[2] = newObject (SISLCURVE)) == SISL_NULL)
		    goto err101;

		  /* Pick out edge curve from a surface. */

		  s1435 (wob[2]->s1, 2, &(qso->c1), spar + 1, &kstat);
		  if (kstat < 0)
		    goto error;
		}
	      else
		qso = wob[2]->edg[2];

	      /* Pick curve from mother object of surface, and make
		 motherobject of curve.                             */

	      s1436(qo1->o1->s1,spar[1],&qcrv,&kstat);
	      if (kstat < 0) goto error;

	      if ((qmotherobj = newObject(SISLCURVE)) == SISL_NULL) goto err101;
	      qmotherobj->c1 = qcrv;
	      qso->o1 = qmotherobj;
	    }

	  /***** Treating edges on sub problems. *****/

	  /* We first have to transform intersection points to the the new
	     intersection format qintdat. */

	  /* UJK, newi */
	  sh6idget ((iobj == 1 ? (ki == 0 ? po1 : (ki == 1 ? qs1 : qs2)) : po1),
		(iobj == 2 ? (ki == 0 ? po2 : (ki == 1 ? qs1 : qs2)) : po2),
		(ki == 0   ? (idiv == 2 ? 1 : 0) : 1) + kpar,
		(ki == 0   ? (idiv == 2 ? spar[1] : spar[0]) : spar[1]),
		    *pintdat, &qintdat, aepsge, &kstat);

	  /* Making new edge object to sub problems. */

	  if ((iobj == 1 ? qso : po1)->iobj == SISLPOINT)
	    uedge[0] = SISL_NULL;
	  else if ((uedge[0] = newEdge (vedge[0]->iedge - (iobj == 1 ? 2 : 0))) == SISL_NULL)
	    goto err101;
	  if ((iobj == 2 ? qso : po2)->iobj == SISLPOINT)
	    uedge[1] = SISL_NULL;
	  else if ((uedge[1] = newEdge (vedge[1]->iedge - (iobj == 2 ? 2 : 0))) == SISL_NULL)
	    goto err101;

	  /* Update edge intersection on sub problems. */

	  sh6idalledg ((iobj == 1 ? qso : po1), (iobj == 1 ? po2 : qso), qintdat, uedge, &kstat);
	  if (kstat < 0)
	    goto error;

	  /* START of update, UJK,jan.93.__________________________ */
	  /* UJK, jan 1993, 1D: test if end pt of curve is intersection pt
	     and not registred on edge.
	     This will very seldom ocurr, boarder line case. */
	  if (qso->c1->idim == 1)
	  {
	     int changes = FALSE;
	     int loop;
	     double endpar;
	     double qt_par[2];
	     SISLPoint  *end_point=SISL_NULL;
	     SISLObject *pt_obj=SISL_NULL;
	     SISLCurve  *pcrv=SISL_NULL;
	     int knum;
	     int ind_missing, ind_kept;
	     SISLIntpt *qt  = SISL_NULL;
	     SISLIntpt *pcl = SISL_NULL;
	     SISLIntpt **up = SISL_NULL;	 /* Array of poiners to intersection point. */

	     /* Get edge points to SUB-problem. */
	     sh6edgpoint (uedge, &up, &knum, &kstat);
	     if (kstat < 0)
		goto error;

	     /* Set up case navigators. */
	     pcrv   = qso->c1;
	     pt_obj = (iobj == 1 ? po2 :po1);
	     ind_missing = (ki == 0 ? (idiv == 2 ? 1 : 0) : 1);
	     ind_kept    = (ki == 0 ? (idiv == 2 ? 0 : 1) : 0);

	     if (knum < 2)
		for (loop = 0; loop < 2; loop++)
		{

		   /* Pick out end point from a curve. */
		   s1438 (pcrv, loop, &end_point, &endpar, &kstat);
		   if (kstat < 0)
		      goto error;
		   if (fabs(end_point->ecoef[0] - pt_obj->p1->ecoef[0]) < aepsge &&
		       (knum == 0 || DNEQUAL(up[0]->epar[0], endpar)))
		   {
		      /* Making intersection point. */
		      double *nullp = SISL_NULL;

		      changes = TRUE;
		      qt_par[ind_kept]    = endpar;
		      qt_par[ind_missing] = spar[ind_missing];
		      qt = hp_newIntpt (2, qt_par, DZERO, SI_ORD,
					SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
					0, 0, nullp, nullp);

		      if (qt == SISL_NULL)
			 goto err101;

		      sh6tohelp (qt,&kstat);
		      if (kstat < 0) goto error;

		      /* get closest point pcl to qt in pintdat. */
		      sh6idnpt(pintdat,&qt,TRUE,&kstat);
		      if (kstat < 0) goto error;
		      kpos=1;
		      if (kstat) goto errinconsis;

		      kpos=2;
		      s6idcpt(*pintdat,qt,&pcl);
		      if (!pcl) goto errinconsis;

		      if (DEQUAL(pcl->epar[ind_kept],qt_par[ind_kept]) &&
			  fabs(pcl->epar[ind_missing] - qt_par[ind_missing])
			  < 0.000001)
		      {
			 qt->epar[ind_missing] = pcl->epar[ind_missing];
			 pcl->epar[ind_missing] = qt_par[ind_missing];
			 sh6tomain (pcl,&kstat);
			 if (kstat < 0) goto error;
			 sh6idcon (pintdat,&qt,&pcl,&kstat);
			 if (kstat < 0) goto error;
		      }
		   }
		   if (end_point)
		      freePoint(end_point);
		   end_point = SISL_NULL;

		}

	     if (changes)
	     {
		/* Clean up and regenerate uedge and qintdat. */
		if (uedge[0] != SISL_NULL)
		   freeEdge (uedge[0]);
		if (uedge[1] != SISL_NULL)
		   freeEdge (uedge[1]);

		if (qintdat != SISL_NULL)
		{
		   freeIntdat (qintdat);
		   qintdat = SISL_NULL;
		}

		sh6idget ((iobj == 1 ? (ki == 0 ? po1 : (ki == 1 ? qs1 : qs2)) : po1),
			  (iobj == 2 ? (ki == 0 ? po2 : (ki == 1 ? qs1 : qs2)) : po2),
			  (ki == 0   ? (idiv == 2 ? 1 : 0) : 1) + kpar,
			  (ki == 0   ? (idiv == 2 ? spar[1] : spar[0]) : spar[1]),
			  *pintdat, &qintdat, aepsge, &kstat);

		/* Making new edge object to sub problems. */

		if ((iobj == 1 ? qso : po1)->iobj == SISLPOINT)
		   uedge[0] = SISL_NULL;
		else if ((uedge[0] = newEdge (vedge[0]->iedge - (iobj == 1 ? 2 : 0))) == SISL_NULL)
		   goto err101;
		if ((iobj == 2 ? qso : po2)->iobj == SISLPOINT)
		   uedge[1] = SISL_NULL;
		else if ((uedge[1] = newEdge (vedge[1]->iedge - (iobj == 2 ? 2 : 0))) == SISL_NULL)
		   goto err101;

		/* Update edge intersection on sub problems. */

		sh6idalledg ((iobj == 1 ? qso : po1), (iobj == 1 ? po2 : qso), qintdat, uedge, &kstat);
		if (kstat < 0)
		   goto error;


	     }
	     if (up) freearray(up);
	  }
	  /* END of update, UJK,jan.93.__________________________ */


	  /* Examine if the subdividing curve intersect the second object. */

	  sh1762 ((iobj == 1 ? qso : po1), (iobj == 1 ? po2 : qso), aepsge,
		  &qintdat, uedge, &kstat);
	  if (kstat < 0)
	    goto error;

	  if (uedge[0] != SISL_NULL)
	    freeEdge (uedge[0]);
	  if (uedge[1] != SISL_NULL)
	    freeEdge (uedge[1]);


	  /* Free mother object of the subdividing curve.  */

	  if (qmotherobj != SISL_NULL) freeObject(qmotherobj);
	  qmotherobj = SISL_NULL;
	  qcrv = SISL_NULL;
	  qso->o1 = qso;

	  /* Examine if there is an intersection point close to the
	     subdivision point. If there is, correct the subdivision point. */

	  /* ALA and UJK 31.10.90 don't change divide point
	     when fixflag is set */
	  if ((fixflag == 0) && idiv == 3 && ki == 0 && qintdat != SISL_NULL)
	    /*  if (idiv == 3 && ki == 0 && qintdat != SISL_NULL) */
	    {
	      tdel = (qso->c1->et[qso->c1->in] -
		      qso->c1->et[qso->c1->ik - 1]) * (double) 0.1;

	      for (kn = 0; kn < qintdat->ipoint; kn++)
		/* UJK, aug.92 Do NOT subdiv in a help point */
		if (sh6ismain(qintdat->vpoint[kn]))
		  if ((fabs (qintdat->vpoint[kn]->epar[kpar]
			     - spar[1]) < fabs (tdel)) &&
		      DNEQUAL (qintdat->vpoint[kn]->epar[kpar],
			       qso->c1->et[qso->c1->in]) &&
		      DNEQUAL (qintdat->vpoint[kn]->epar[kpar],
			       qso->c1->et[qso->c1->ik - 1]))
		    spar[1] = qintdat->vpoint[kn]->epar[kpar];
	    }

	  /*ujk, ala 921218, dont't split very close to a
	     new intersection point */
	  else if ((fixflag) && idiv == 3 && ki == 0 && qintdat != SISL_NULL)
	  {
	     tdel = (qso->c1->et[qso->c1->in] -
		     qso->c1->et[qso->c1->ik - 1]) * (double) 0.000001;

	     for (kn = 0; kn < qintdat->ipoint; kn++)
		if (DNEQUAL(qintdat->vpoint[kn]->epar[kpar],spar[1])
		    && (fabs (qintdat->vpoint[kn]->epar[kpar]
			      - spar[1]) < fabs (tdel))) break;

	     if (kn <  qintdat->ipoint)
		/* Using midpoint */
	     {
		spar[1] = s1792 (qo1->s1->et2, qo1->s1->ik2, qo1->s1->in2);
		fixflag = 0;
	     }
	  }

	  /* TDO,UJK 02.08.89 */
	  /* A possible better strategy for
	     subdividing is to divide in the
	     middle value of the first and last
	     intersection on the subdividing curve.*/
/*
	  if (idiv == 3 && ki == 0 && qintdat != SISL_NULL)
	    {
	      int kn;
	      double tdel,tdiv,tmin,tmax;

	      if(qintdat->ipoint > 0)
	      {
		  tmin = qintdat->vpoint[0]->epar[kpar];
		  tmax = qintdat->vpoint[0]->epar[kpar];

		  for (kn=1;kn<qintdat->ipoint;kn++)
		    {
		      if (tmin < qintdat->vpoint[kn]->epar[kpar])
			tmin = qintdat->vpoint[kn]->epar[kpar];

		      if (tmax > qintdat->vpoint[kn]->epar[kpar])
			tmax = qintdat->vpoint[kn]->epar[kpar];

		    }

		  tdiv = (tmax + tmin)/(double)2.0;

		  if (DNEQUAL(tdiv,qso->c1->et[qso->c1->in]) &&
		      DNEQUAL(tdiv,qso->c1->et[qso->c1->ik-1]))
		    spar[1]=tdiv;

		}
	    }
*/


	  if (kstat)
	    {
	      /* Total number of points. */

	      knum = (*pintdat == SISL_NULL ? 0 : (*pintdat)->ipoint) +
		qintdat->ipoint;

	      *jstat = 1;	/* Mark intersection found. */

	      /* Intersection found and we have to register the intersection
		 points. */

	      /* UJK newi */
	      sh1782 (po1, po2, aepsge, qintdat,
		      (ki == 0 ? (idiv == 2 ? 1 : 0) : 1) + kpar,
		      (ki == 0 ? (idiv == 2 ? spar[1] : spar[0]) : spar[1]),
		      pintdat, &idummy, &kstat);
	      if (kstat < 0)
		goto error;

	      /* UJK newi divide surface */
	      /* UPDATE: ? what about help points from s1782 knum?? */
	      if (qpt != SISL_NULL && (*pintdat)->ipoint == knum)
		{
		  /* Find the closest poin to qpt. */

		  s6idcpt (*pintdat, qpt, &qp);

		  /* UJK newi, unite the points : */
		  sh6idnewunite (po1, po2, pintdat, &qpt, &qp, (double) 0.5,
				 aepsge, &kstat);
		  if (kstat < 0)
		    goto error;
		}

	    }

	  if (qintdat != SISL_NULL)
	    {
	      freeIntdat (qintdat);
	      qintdat = SISL_NULL;
	    }
	}
      if (qs1 != SISL_NULL)
	freeObject (qs1);
      if (qs2 != SISL_NULL)
	freeObject (qs2);
    }
  else
    goto err121;


  goto out;

/* Error. Inconsistency.  */

errinconsis:*jstat = -231;
  s6err ("sh1762_s9div", *jstat, kpos);
  goto out;

/* Error. Kind of object does not exist.  */

err121:*jstat = -121;
  s6err ("sh1762_s9div", *jstat, kpos);
  goto out;

/* Error in space allocation.  */

err101:*jstat = -101;
  s6err ("sh1762_s9div", *jstat, kpos);
  goto out;

/* Error in lower level routine.  */

error:*jstat = kstat;
  s6err ("sh1762_s9div", *jstat, kpos);
  goto out;

out:;
}

#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9update (SISLObject * po1, SISLObject * po2, double aepsge,
		 SISLIntdat ** pintdat, SISLEdge ** vedge[], int *jstat)
#else
static void
sh1762_s9update (po1, po2, aepsge, pintdat, vedge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **pintdat;
     SISLEdge **vedge[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Uppdating intersection data if possible.
*
*
*
* INPUT      : po1      - First object in intersection.
*              po2      - Second object in intersection.
*              aepsge    - Geometry resolution.
*
*
* OUTPUT     : pintdat  - intersection data.
*              vedge[2] - intersection on edges.
*              jstat    - status messages
*                          = 2     : no intersection found
*                          = 1     : Intersection found.
*                          = 0     : not a simple case.
*                          < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*
*********************************************************************
*/
{
  /* UJK newi */
  int ki, no_new;
  SISLIntpt *qt = SISL_NULL;

  int kpos = 0;
  int kstat = 0;
  int kdim;
  SISLObject *qo;
  SISLIntpt **up = SISL_NULL;

  /* Test input.  */

  kdim = (po1->iobj == SISLPOINT ? po1->p1->idim :
	  (po1->iobj == SISLCURVE ? po1->c1->idim : po1->s1->idim));

  if (kdim != (po2->iobj == SISLPOINT ? po2->p1->idim :
	       (po2->iobj == SISLCURVE ? po2->c1->idim : po2->s1->idim)))
    goto err106;

  /* Initiate to no intersection. */

  *jstat = 2;

  if (po1->iobj == SISLPOINT || po2->iobj == SISLPOINT)
    {
      int kturn = 0;
      int knum = 0;

      if (po1->iobj != SISLPOINT)
	{
	  qo = po1;
	  po1 = po2;
	  po2 = qo;
	  kturn = 1;
	}

      knum = (*vedge)[1 - kturn]->ipoint;


      /* UPDATE ALA 010993. Start */
      if (knum == 0 && (*pintdat) != SISL_NULL)
      {
	for (ki = 0; ki < (*pintdat)->ipoint; ki++)
	  if (po2->iobj == SISLCURVE)
	  {
	    if ((*pintdat)->vpoint[ki]->epar[0] > po2->c1->et[po2->c1->ik-1] &&
		(*pintdat)->vpoint[ki]->epar[0] < po2->c1->et[po2->c1->in])
	    {
	      knum = 1;
	    }
	  }
	  else
	  {
	    if ((*pintdat)->vpoint[ki]->epar[0] > po2->s1->et1[po2->s1->ik1-1] &&
		(*pintdat)->vpoint[ki]->epar[0] < po2->s1->et1[po2->s1->in1]   &&
		(*pintdat)->vpoint[ki]->epar[1] > po2->s1->et2[po2->s1->ik2-1] &&
		(*pintdat)->vpoint[ki]->epar[1] < po2->s1->et2[po2->s1->in2])
	    {
	      knum = 1;
	    }
	  }
      }
      /* UPDATE ALA 010993.  End */


      if (knum > 1)
	{
	  /* sh1762_s9edgpoint ((*vedge), &up, &knum, &kstat); */
	  sh6edgpoint ((*vedge), &up, &knum, &kstat);
	  if (kstat < 0)
	    goto error;
	}


      /* We have more than one intersection point on the edges.
	 If the dimension is one and the second object is a point
	 we just connect the point else we kill these points and
	 try to find a new intersection point. */

      if (knum > 1)
	{
	  if (po2->iobj == SISLSURFACE && kdim == 1)
	    {
	      int ksimple;
	      if (po2->o1 == po2)
		ksimple = 0;
	      else
		ksimple = 1;

	      /* UPDATE: UJK, new parameter turn ?? */
	      sh1762_s9edgpscon ((*vedge)[1 - kturn], po1->p1->ecoef[0],
				 po2->s1, ksimple, *pintdat, aepsge, &kstat);
	      if (kstat < 0)
		goto error;
	      else if (kstat)
		*jstat = 0;	/* Not a simple case. */
	    }

	    else  if (po2->iobj == SISLSURFACE && kdim == 2 && knum == 2)
	       {
		  /* 2D point surf, connect */
		  sh6idcon (pintdat, up, up + 1, &kstat);
		  if (kstat < 0)
		     goto error;
	        }

	    else
	    {
	      /* UJK newi */
	      for (ki = 1; ki < knum; ki++)
		{
		   sh6idnewunite (po1, po2, pintdat, &up[0], &up[ki],
				(double) 0.5, aepsge, &kstat);
		  if (kstat < 0)
		    goto error;

		}
	      qt = up[0];

	      ki = (*vedge)[1 - kturn]->iedge;
	      freeEdge ((*vedge)[1 - kturn]);
	      if (((*vedge)[1 - kturn] = newEdge (ki)) == SISL_NULL)
		goto err101;
	      knum = 0;
	    }
	}

      if (knum == 0)
	{
	  double spar[2];

	  if (po2->iobj == SISLCURVE)
	    {
	      double tstart, tend;

	      tstart = po2->c1->et[po2->c1->ik - 1];
	      tend = po2->c1->et[po2->c1->in];
	      spar[0] = (tstart + tend) * (double) 0.5;


	      s1771 (po1->p1, po2->o1->c1, aepsge,
		     tstart, tend, spar[0], spar, &kstat);
	      if (kstat < 0)
		goto error;

	      if (kstat == 1)
		/*Intersection point found. Control edges. */
		if (DEQUAL (spar[0], tstart) || DEQUAL (spar[0], tend))
		  kstat = 0;
	    }
	  else if (po2->iobj == SISLSURFACE)
	    {
	      double sstart[2], send[2];

	      sstart[0] = po2->s1->et1[po2->s1->ik1 - 1];
	      sstart[1] = po2->s1->et2[po2->s1->ik2 - 1];

	      send[0] = po2->s1->et1[po2->s1->in1];
	      send[1] = po2->s1->et2[po2->s1->in2];

	      spar[0] = (sstart[0] + send[0]) * (double) 0.5;
	      spar[1] = (sstart[1] + send[1]) * (double) 0.5;

	      s1773 (po1->p1, po2->o1->s1, aepsge, sstart, send, spar, spar, &kstat);
	      if (kstat < 0)
		goto error;

	      if (kstat == 1)
		/*Intersection point found. Control edges. */
		if (DEQUAL (spar[0], sstart[0]) || DEQUAL (spar[0], send[0])
		|| DEQUAL (spar[1], sstart[1]) || DEQUAL (spar[1], send[1]))
		  kstat = 0;
	    }



	  /* UJK, October 91, 2D crv and surf's may be degenerate,
	     continue when iteration fails */
	  if (po1->p1->idim == 2)
	  {
	     if ((po2->iobj == SISLSURFACE && kstat == 9)||
		 (po2->iobj == SISLCURVE   && kstat != 1))
	     {
		*jstat = 0;
		goto out;
	     }
	  }


	  /* TESTING UJK !!!!!!!!!!!!!!!!!!!!! */
	    /* UJK, August 92, 1D crvs may be "degenerate",
	       continue when iteration fails */
	    if (kstat != 1 && po1->p1->idim == 1)
	    {
	       *jstat = 0;
	       goto out;
	    }


	    if (kstat == 1)	/* Intersection point found. */
	    {
	      *jstat = 1;	/* Mark intersection found.  */

	      /* UJK newi */
	      if (qt)
		{
		  /* We have an instance of a point, use it */
		  for (ki = 0; ki < qt->ipar; ki++)
		    qt->epar[ki] = spar[ki];
		}
	      else
		{
		  /* Making intersection point. */
		  double *nullp = SISL_NULL;
		  qt = hp_newIntpt (po2->iobj, spar, DZERO, SI_ORD,
				    SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
				    0, 0, nullp, nullp);

		  if (qt == SISL_NULL)
		    goto err101;

		  /* Uppdating pintdat. */
		  sh6idnpt (pintdat, &qt, 1, &kstat);
		  if (kstat < 0)
		    goto error;

		  /* Set pretopology */
		  if (po2->iobj == SISLCURVE)
		    {
		      /* Case point, curve */
		      if (po1->p1->idim == 1)
			{
			  /* 1D point curve treated,
			     2D, 3D is set to undef */

			  sh1781 ((kturn ? po2 : po1),
				  (kturn ? po1 : po2),
				  aepsge, pintdat, qt, &no_new, &kstat);
			  if (kstat < 0)
			    goto error;
			}
		      else
			{
			  /* UPDATE (ujk) 1D touch well defined ? */
			  /* Case point, surface */
			  if (po1->p1->idim == 2)
			    {
			      /* 2D point surface is treated,
			         1D, 3D is set to undef */
			      sh1786 ((kturn ? po2 : po1),
				      (kturn ? po1 : po2),
				      aepsge, pintdat, qt, &no_new, &kstat);

			      if (kstat < 0)
				goto error;
			    }
			}

		    }
		}
	    }
	}
    }
  else if (po1->iobj == SISLCURVE || po2->iobj == SISLCURVE)
    {
      int kturn1 = 1, kturn2 = 0;

      if (po1->iobj != SISLCURVE)
	{
	  qo = po1;
	  po1 = po2;
	  po2 = qo;
	  kturn1 = 0;
	  kturn2 = 2;
	}

      if ((*vedge)[0]->ipoint + (*vedge)[1]->ipoint > 1)
	{
	  int knum;

	  /* sh1762_s9edgpoint ((*vedge), &up, &knum, &kstat); */
	  sh6edgpoint ((*vedge), &up, &knum, &kstat);
	  if (kstat < 0)
	    goto error;

	  if (knum > 1)
	    {
	      int ki;

	      /* We have more than one intersection point on the edges.
	         We therefor kill these points and
	         try to find a new intersection point. */

	      /* UJK newi CONNECT */
	      for (ki = 1; ki < knum; ki++)
		{
		   /* sh6idnewunite (po1, po2, pintdat, &up[0], &up[ki],
		      (double) 0.5, aepsge, &kstat); */
		   sh6idcon(pintdat, &up[0], &up[ki], &kstat);
		  if (kstat < 0)
		    goto error;

		}

	      *jstat = 1;
	      goto out;

	      /* qt = up[0];

		 ki = (*vedge)[0]->iedge;
		 freeEdge ((*vedge)[0]);
		 if (((*vedge)[0] = newEdge (ki)) == SISL_NULL)
		 goto err101;
		 ki = (*vedge)[1]->iedge;
		 freeEdge ((*vedge)[1]);
		 if (((*vedge)[1] = newEdge (ki)) == SISL_NULL)
		 goto err101;
		 knum = 0; */
	    }
	}

      if ((*vedge)[0]->ipoint + (*vedge)[1]->ipoint == 0)
	{
	  double spar[3];

          /* UPDATE ALA 010993. Start */
          if ((*pintdat) != SISL_NULL)
          {
	  for (ki = 0; ki < (*pintdat)->ipoint; ki++)
	     if (po2->iobj == SISLCURVE)
	     {
	       if ((*pintdat)->vpoint[ki]->epar[0] > po1->c1->et[po1->c1->ik-1] &&
		   (*pintdat)->vpoint[ki]->epar[0] < po1->c1->et[po1->c1->in] &&
	           (*pintdat)->vpoint[ki]->epar[1] > po2->c1->et[po2->c1->ik-1] &&
		   (*pintdat)->vpoint[ki]->epar[1] < po2->c1->et[po2->c1->in])
	         goto out;
	     }
	     else
	     {
	       if ((*pintdat)->vpoint[ki]->epar[kturn2] > po1->c1->et[po1->c1->ik-1] &&
		   (*pintdat)->vpoint[ki]->epar[kturn2] < po1->c1->et[po1->c1->in] &&
	           (*pintdat)->vpoint[ki]->epar[kturn1] > po2->s1->et1[po2->s1->ik1-1] &&
		   (*pintdat)->vpoint[ki]->epar[kturn1] < po2->s1->et1[po2->s1->in1]   &&
		   (*pintdat)->vpoint[ki]->epar[kturn1+1] > po2->s1->et2[po2->s1->ik2-1] &&
		   (*pintdat)->vpoint[ki]->epar[kturn1+1] < po2->s1->et2[po2->s1->in2])
	         goto out;
	     }
         }
         /* UPDATE ALA 010993.  End */

	 if (po2->iobj == SISLCURVE)
	    {
	      double tstart1, tend1;
	      double tstart2, tend2;

	      tstart1 = po1->c1->et[po1->c1->ik - 1];
	      tend1 = po1->c1->et[po1->c1->in];
	      spar[0] = (tstart1 + tend1) * (double) 0.5;

	      tstart2 = po2->c1->et[po2->c1->ik - 1];
	      tend2 = po2->c1->et[po2->c1->in];
	      spar[1] = (tstart2 + tend2) * (double) 0.5;


	      s1770 (po1->o1->c1, po2->o1->c1, aepsge, tstart1,
		     tstart2, tend1, tend2, spar[0], spar[1], spar, spar + 1, &kstat);
	      if (kstat < 0)
		goto error;

	      if (kstat == 2)
		{
		  /* Search for a better start point for the
		     iteration. */
		  sh6cvvert(po1->c1, po2->c1, spar, spar+1);

		  /* Iterate. */
		  kstat = 0;
		  s1770 (po1->o1->c1, po2->o1->c1, aepsge, tstart1,
			 tstart2, tend1, tend2, spar[0], spar[1],  
			 spar, spar + 1, &kstat);
		  if (kstat < 0)
		    { kpos=__LINE__; goto error; }
		}

	      if (kstat == 1)
		/*Intersection point found. Control edges. */
		if (DEQUAL (spar[0], tstart1) || DEQUAL (spar[0], tend1)
		    || DEQUAL (spar[1], tstart2) || DEQUAL (spar[1], tend2))
		  kstat = 0;
	    }
	  else if (po2->iobj == SISLSURFACE)
	    {
	      double tstart, tend;
	      double sstart[2], send[2];

	      tstart = po1->c1->et[po1->c1->ik - 1];
	      tend = po1->c1->et[po1->c1->in];
	      spar[kturn2] = (tstart + tend) * (double) 0.5;


	      sstart[0] = po2->s1->et1[po2->s1->ik1 - 1];
	      sstart[1] = po2->s1->et2[po2->s1->ik2 - 1];

	      send[0] = po2->s1->et1[po2->s1->in1];
	      send[1] = po2->s1->et2[po2->s1->in2];

	      spar[kturn1] = (sstart[0] + send[0]) * (double) 0.5;
	      spar[kturn1 + 1] = (sstart[1] + send[1]) * (double) 0.5;

	      kstat = 0;
	      s1772 (po1->o1->c1, po2->o1->s1, aepsge, tstart, sstart, tend, send,
		     spar[kturn2], &spar[kturn1],
		     &spar[kturn2], &spar[kturn1], &kstat);
	      if (kstat < 0)
		goto error;

		if (kstat == 3)
		{
		   /* FLAT */
		   *jstat = 0;
		   goto out;
		}

	      /* UJIK, Retry, with better startpoint */
	      if (kstat == 2)
	      {
		/* No intersection point is found. Try again with a new
		   start point to the iteration.  */

		sh6closevert(po1->c1,po2->s1,&spar[kturn2],&spar[kturn1]);
		kstat = 0;
		s1772 (po1->o1->c1, po2->o1->s1, aepsge, tstart, sstart, tend, send,
		       spar[kturn2], &spar[kturn1],
		       &spar[kturn2], &spar[kturn1], &kstat);
		if (kstat < 0)
		  goto error;
	      }

	      if (kstat == 1)
		/*Intersection point found. Control edges. */
		if (DEQUAL (spar[kturn2], tstart) ||
		    DEQUAL (spar[kturn2], tend) ||
		    DEQUAL (spar[kturn1], sstart[0]) ||
		    DEQUAL (spar[kturn1], send[0]) ||
		    DEQUAL (spar[kturn1 + 1], sstart[1]) ||
		    DEQUAL (spar[kturn1 + 1], send[1]))
		  kstat = 0;
	    }

	  if (kstat == 1)	/* Intersection point found. */
	    {

	      *jstat = 1;	/* Mark intersection found.  */

	      /* UJK newi */
	      if (qt)
		{
		  /* We have an instance of a point, use it */
		  for (ki = 0; ki < qt->ipar; ki++)
		    qt->epar[ki] = spar[ki];
		}
	      else
		{
		  /* Making intersection point. */
		  double *nullp = SISL_NULL;
		  qt = hp_newIntpt (po1->iobj + po2->iobj, spar, DZERO, SI_ORD,
				    SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
				    0, 0, nullp, nullp);

		  if (qt == SISL_NULL)
		    goto err101;

		  /* Uppdating pintdat. */
		  sh6idnpt (pintdat, &qt, 1, &kstat);
		  if (kstat < 0)
		    goto error;

		  /* Set pretopology */
		  if (po2->iobj == SISLCURVE)
		    {
		      /* Case curve, curve */
		      if (po1->c1->idim == 2)
			{
			  /* Only 2D is treated */
			  sh1780 (po1, po2,
				  aepsge, pintdat, qt, &no_new, &kstat);
			  if (kstat < 0)
			    goto error;
			}
		    }
		  else if (po2->iobj == SISLSURFACE)
		    {
		      /* Case curve, surface */
		      if (po1->c1->idim == 3)
			{
			  /* Only 3D  */
			  sh1779 ((kturn1 ? po1 : po2),
				  (kturn1 ? po2 : po1),
				  aepsge, pintdat, qt, &no_new, &kstat);
			  if (kstat < 0)
			    goto error;
			}
		    }

		}
	    }
	}
    }
  else if (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE)
    {
      if ((*vedge)[0]->ipoint + (*vedge)[1]->ipoint > 1)
	{
	  /* We have more than one intersection point on the edges,
             we therefor connect these points to each other. */
	  int ksimple;

	  if (po1->psimple == po2)
	    ksimple = 1;
	  else
	    ksimple = 0;

	  sh1762_s9edgsscon ((*vedge), po1->s1, po2->s1, *pintdat, ksimple,
			     aepsge, &kstat);
	  if (kstat < 0)
	    goto error;
	  else if (kstat)
	    *jstat = 0;		/* Not a simple case. */
	}
    }
  else
    goto err121;

  goto out;

/* Error. Kind of object does not exist.  */

err121:*jstat = -121;
  s6err ("sh1762_s9update", *jstat, kpos);
  goto out;

/* Error in input. Conflicting dimensions.  */

err106:*jstat = -106;
  s6err ("s1770", *jstat, kpos);
  goto out;

/* Error in space allocation.  */

err101:*jstat = -101;
  s6err ("sh1762_s9update", *jstat, kpos);
  goto out;

/* Error in lower level routine.  */

error:*jstat = kstat;
  s6err ("sh1762_s9update", *jstat, kpos);
  goto out;

out:if (up != SISL_NULL)
    freearray (up);
}

#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9con (SISLObject * po1, SISLObject * po2, double aepsge,
	      SISLIntdat ** pintdat, SISLEdge * vedge[], int *jstat)
#else
static void
sh1762_s9con (po1, po2, aepsge, pintdat, vedge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **pintdat;
     SISLEdge *vedge[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if one part of an object coincide in one part
*              of an other object. Both object must intersect in at least two
*              ends/edges. This rutine also use a rotated SISLbox tests
*              to test if intersection is possible.
*
*
*
* INPUT      : po1      - The first SISLObject to check.
*              po2      - The second SISLObject to check.
*              aepsge   - Geometrical resolution.
*              vedge[]  - SISLEdge intersection.
*              *jstat    - Flag
*                          = 202 : Complicated point-surface intersection
*                                  in 3D. Perform extra interception test.
*
*
* INPUT/OUTPUT:pintdat  - Intersection data.
*
*
* OUTPUT     :
*              jstat    - status messages
*                           = 2     : No intersection..
*                           = 1     : Coinside found.
*                           = 0     : no coinside.
*                           < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* CALLS      : sh1762_s9intercept - Perform improved box tests to intercept
*                            recursion.
*              sh1762_s9coincide  - Test coincidence between crv/crv and surf/crv objects.
*              sh1762_s9toucharea - Test coincidence between surf/surf objects.
*              sh6edgpoint        - Fetch intersection points on edges.
*
*
* WRITTEN BY : Arne Laksaa,  SI, 89-05.
* REWISED BY : Vibeke Skytt, SI, 91-01.
*              UJK,          SI, 91-06
* REVISED BY : Michael Floater SI, 91-09
		- check for existence of touching areas.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int ki,kj;			/* Counter.                                */
  int knum = 0;			/* Number of intersection points on edges. */
  SISLIntpt **up = SISL_NULL;	/* Intersection points on edges.           */
  SISLdir *qd1, *qd2;		/* Direction cones of objects.             */
  SISLIntpt *qpt;               /* Evt 3. intersection point.              */
  int knpar=po1->iobj+po2->iobj; /* Number of parameter directions.        */
  int kcrv1;                    /* Indicates if 1. object is a curve.      */
  int kcrv2;                    /* Indicates if 2. object is a curve.      */
  int pretop[2][4];
  SISLObject *qobj;
  int ind1, ind2, perm[2], obj, ipar;
  int klist1, klist2;
  int linear = FALSE;
  int kpt,kpt2;                 /* Number of elements in int. list.        */
  int kstat2 = 0;               /* Remember status from s9toucharea.       */
  double mintang1;
  double mintang2;
  double tboxsize1;
  double tboxsize2;
  int kxintercept = (*jstat == 202);  /* Extra interception       */

  /*int loopcount;*/		/* Count up num intpts in a list. */
  int one_edge = 0;             /* Indicates if all intersection points
				   lies on one edge in each surface.     */
  SISLPtedge *qpt1, *qpt2;      /* Pointers used to traverse edge intersections. */

  /* Set kcrv parameters.  */

  kcrv1 = (po1->iobj == SISLCURVE) ? 1 : 0;
  kcrv2 = (po2->iobj == SISLCURVE) ? 1 : 0;

  if ((po1->iobj == SISLPOINT && po1->p1->idim == 1) ||
      (po2->iobj == SISLPOINT && po2->p1->idim == 1))
    *jstat = 0;
  else
    {

       if (po1->iobj == SISLPOINT) qd1 = SISL_NULL;
       else
	  qd1 = (po1->iobj == SISLCURVE ? po1->c1->pdir : po1->s1->pdir);

       if (po2->iobj == SISLPOINT) qd2 = SISL_NULL;
       else
	  qd2 = (po2->iobj == SISLCURVE ? po2->c1->pdir : po2->s1->pdir);

       knum = 0;
       if (vedge[0] != SISL_NULL) knum += vedge[0]->ipoint;
       if (vedge[1] != SISL_NULL) knum += vedge[1]->ipoint;

      if (knum > 0)
	{
	  /* Organize intersection points on an array. */

	  /* sh1762_s9edgpoint (vedge, &up, &knum, &kstat); */
	  sh6edgpoint (vedge, &up, &knum, &kstat);
	  if (kstat < 0)
	    goto error;
	}

      /* We test coincide by linearity. If the two object is liniar
         and have end/edge intersection we just connect these
         intersection points, else we have no internal intersections. */

      if (po1->iobj == SISLCURVE)
      {
	 tboxsize1 = po1->c1->pbox->e2max[2][0] - po1->c1->pbox->e2min[2][0];
	 if (po1->c1->idim > 1)
	    tboxsize1 = MAX(tboxsize1,
			   po1->c1->pbox->e2max[2][1] - po1->c1->pbox->e2min[2][1]);
	 if (po1->c1->idim > 2)
	    tboxsize1 = MAX(tboxsize1,
			   po1->c1->pbox->e2max[2][2] - po1->c1->pbox->e2min[2][2]);
	 mintang1 = aepsge/((double)2*tboxsize1);
      }
      else  if (po1->iobj == SISLSURFACE)
	 mintang1 = ANGULAR_TOLERANCE/(double)10;

      if (po2->iobj == SISLCURVE)
      {
	 tboxsize2 = po2->c1->pbox->e2max[2][0] - po2->c1->pbox->e2min[2][0];
	 if (po2->c1->idim > 1)
	    tboxsize2 = MAX(tboxsize2,
			   po2->c1->pbox->e2max[2][1] - po2->c1->pbox->e2min[2][1]);
	 if (po2->c1->idim > 2)
	    tboxsize2 = MAX(tboxsize2,
			   po2->c1->pbox->e2max[2][2] - po2->c1->pbox->e2min[2][2]);
	 mintang2 = aepsge/((double)2*tboxsize2);
      }
      else if (po2->iobj == SISLSURFACE)
	 mintang2 = ANGULAR_TOLERANCE/(double)10;

      /* if (qd1->igtpi || qd2->igtpi || qd1->aang > ANGULAR_TOLERANCE ||
	  qd2->aang > ANGULAR_TOLERANCE) */
      if (qd1 == SISL_NULL || qd2 == SISL_NULL)
	 *jstat = 0;
      else if (qd1->igtpi || qd2->igtpi || qd1->aang > mintang1 ||
	  qd2->aang > mintang2)
	*jstat = 0;
      else if (knum == 2)
	/* Newi (ujk) When linear and 2 points, we know how to set
           the pretopology for curves, this is done a bit further down */
	linear = TRUE;
      else if (po1->iobj + po2->iobj < 2*SISLSURFACE)
	{
	  if (knum > 1)
	    {
	      /* We have more than one intersection point on the edges.
                 We therefore connect these points. */
	      /* UPDATE (ujk) don't like this connection */

	      for (ki = 0; ki < knum; ki++)
		sh6tomain (up[ki], &kstat);

	      for (ki = 1; ki < knum; ki++)
		{
		  sh6idcon (pintdat, &up[ki - 1], &up[ki], &kstat);
		  if (kstat < 0)
		    goto error;
		}
	      *jstat = 1;
	    }
	  else
	    *jstat = 2;

	  goto out;		/* Test performed.  */
	}

      if (knum >= 2 &&
	  po1->iobj == SISLSURFACE &&
	  po2->iobj == SISLSURFACE)
      {
	 /* VSK. Change test on possibility of coincidence.
	    More than two intersection points on the edges.
	    Check if there is
	    coincidence between the (surface) objects.
	    Fetch all closed loops. Then call s9toucharea to
	    see if the surfaces coincide everywhere inside the loop. */

	 for (kstat2=0, kpt=0; kpt<knum; kpt+=kpt2)
	 {
	    sh6floop(up+kpt,knum-kpt,&kpt2,&kstat);

	    if (kstat == 1)
	    {
	       sh1762_s9toucharea (po1, po2, aepsge, kpt2, up+kpt, &kstat);
	       /*fprintf (stdout, "\n s9_toucharea, kstat=%d", kstat); */
	       if (kstat < 0)
		  goto error;
	       kstat2 = MAX(kstat2,kstat);
	    }
	    else if (kpt == 0 && kpt2 == knum)
	    {
	       /* Only one open edge curve. Test if the entire
		  curve lies on one edge in each surface.  */

	       for (one_edge=1, ki=1; ki<knum; ki++)
	       {
		  sh6comedg(po1, po2, up[ki-1], up[ki], &kstat);
		  if (kstat < 0) goto error;

		  if (kstat != 3) one_edge = 0;  /* Not a common edge. */
	       }
	    }

	 }

	 *jstat = kstat2;

	 if (kstat2 == 1)
	 {
	    /* Do something with the pertopology. */
	    /* fprintf (stdout, "\n Coincidence, kstat=%d", kstat); */
	 }
      }

      if (knum < 2 || (po1->iobj == SISLSURFACE &&
	  po2->iobj == SISLSURFACE && one_edge))
	{
	  /* Number of intersection points on the edges is less than
             two. Try to intercept further subdivision by performing
             improved box tests.  */

	  kstat = (kxintercept) ? 202 : 0;
	  sh1762_s9intercept (po1, po2, aepsge, knum, up, &kstat);
	  if (kstat < 0)
	    goto error;

	  *jstat = kstat;
	}
      else if (knum == 2 && !(po1->iobj == SISLSURFACE &&
			      po2->iobj == SISLSURFACE))
	{
	  /* Two intersection points on the edges. Check if there is
             coincidence between the objects.  */

	  if (linear)
	    kstat = 1;
	  else
	    {
	      sh1762_s9coincide (po1, po2, aepsge, knum, up, &kstat);
	      if (kstat < 0)
		goto error;
	    }

	  *jstat = kstat;

	  if (kstat == 1)
	    {
	      int kstat1 = 0;

	      for (ki = 0; ki < knum; ki++)
		sh6tomain (up[ki], &kstat);

	      sh6idcon (pintdat, &up[0], &up[1], &kstat);
	      if (kstat < 0)
		goto error;
	      /* Newi (ujk) */
	      /*	      for (ind1 = 0; ind1 < 2; ind1++)
		for (ind2 = 0; ind2 < 4; ind2++)
		pretop[ind1][ind2] = SI_UNDEF; */

	      /* Fetch existing pretopology. */
	      sh6gettop (up[0], -1, &pretop[0][0], &pretop[0][1],
			 &pretop[0][2], &pretop[0][3], &kstat1);

	      sh6gettop (up[1],  -1, &pretop[1][0], &pretop[1][1],
			 &pretop[1][2], &pretop[1][3], &kstat1);

	      for (qobj = po1, obj = 0, ipar = 0; obj < 2;
	       qobj = po2, obj++, ipar = ((po1->iobj == SISLCURVE) ? 1 : 2))
		/* Pretopology for curves */
		if (qobj->iobj == SISLCURVE)
		  {
		    if (up[0]->epar[ipar] < up[1]->epar[ipar])
		      {
			perm[0] = 0;
			perm[1] = 1;
			ind1 = 0;
			ind2 = 1;
		      }
		    else
		      {
			perm[0] = 1;
			perm[1] = 0;
			ind1 = 1;
			ind2 = 0;
		      }

		    /* Left point on curve */
		    pretop[ind1][1 + 2 * obj] = SI_ON;
		    /* Point at edge */
		    if (pretop[ind1][2 * obj] != SI_IN &&
			pretop[ind1][2 * obj] != SI_OUT &&
			DEQUAL (up[perm[0]]->epar[ipar],
				qobj->c1->et[qobj->c1->ik - 1]))
		      {
			/* Point at edge */
			pretop[ind1][2 * obj] = SI_AT;
		      }

		    /* Right point of curve */
		    pretop[ind2][2 * obj] = SI_ON;
		    if (pretop[ind2][1 + 2 * obj] != SI_IN &&
			pretop[ind2][1 + 2 * obj] != SI_OUT &&
			DEQUAL (up[perm[1]]->epar[ipar],
				qobj->c1->et[qobj->c1->in]))
		      {
			/* Point at edge */
			pretop[ind2][1 + 2 * obj] = SI_AT;
		      }

		    /*    / Left point on curve /
		    if (DEQUAL (up[perm[0]]->epar[ipar],
				qobj->c1->et[qobj->c1->ik - 1]))
		      {
			* Point at edge *
			pretop[ind1][2 * obj] = SI_AT;
			pretop[ind1][1 + 2 * obj] = SI_ON;
		      }
		    else
		      {
			pretop[ind1][1 + 2 * obj] = SI_ON;
		      }

		    * Right point of curve *
		    if (DEQUAL (up[perm[1]]->epar[ipar],
				qobj->c1->et[qobj->c1->in]))
		      {
			* Point at edge *
			pretop[ind2][2 * obj] = SI_ON;
			pretop[ind2][1 + 2 * obj] = SI_AT;
		      }
		    else
		      {
			pretop[ind2][2 * obj] = SI_ON;
			} */

		  }
	      sh6getlist (up[0], up[1], &klist1, &klist2, &kstat);
	      if (kstat != 0)
		{
		  kstat = -1;
		  goto error;
		}

	      sh6settop (up[0], -1, pretop[0][0], pretop[0][1],
			 pretop[0][2], pretop[0][3], &kstat);
	      if (kstat < 0)
		goto error;

	      sh6settop (up[1], -1, pretop[1][0], pretop[1][1],
			 pretop[1][2], pretop[1][3], &kstat);
	      if (kstat < 0)
		goto error;


	      if (knpar<4 && (*pintdat)->ipoint == 3)
		 {
		    /* There is 3 intersection points. Test if the 3. point
		       lies between the endpoints of the coincidence curve.
		       First fetch 3. point.         */

		    for (kj=0; kj<3; kj++)
		      {
			 qpt = (*pintdat)->vpoint[kj];
			 if (qpt!=up[0] && qpt!=up[1]) break;
		      }

		    /* Check if the point lies inside the current
		       intersection area.  */

		    sh6isinside(po1,po2,qpt,&kstat);
		    if (kstat < 0) goto error;

		    if (kstat == 1)
		    {
		       /* Check parameter value of evt curves. */

		       if ((kcrv1 &&
			   (up[0]->epar[0] < qpt->epar[0] &&
			    qpt->epar[0] < up[1]->epar[0])) ||
			   (up[1]->epar[0] < qpt->epar[0] &&
			    qpt->epar[0] < up[0]->epar[0])) kcrv1 = -1;

		       if ((kcrv2 &&
			   (up[0]->epar[po1->iobj] < qpt->epar[po1->iobj] &&
			    qpt->epar[po1->iobj] < up[1]->epar[po1->iobj])) ||
			   (up[1]->epar[po1->iobj] < qpt->epar[po1->iobj] &&
			     qpt->epar[po1->iobj] < up[0]->epar[po1->iobj]))
			  kcrv2 = -1;

		       if (kcrv1 < 1 && kcrv2 < 1)
		       {

			  /* The point lies inside the coincidence curve.
			     Place it between the endpoints of the curve. */

			  sh6tomain(qpt,&kstat);
			  sh6insertpt(up[0],up[1],qpt,&kstat);
			  if (kstat < 0) goto error;
			}
		     }
		 }

	    }
	}
      else if (knum > 2 && !(po1->iobj == SISLSURFACE &&
			     po2->iobj == SISLSURFACE))
      {
	 /* There is more than two edge intersection and it is no
	    surface - surface intersection. Check if the edge
	    intersections are already connected. */

	 sh6floop(up, knum, &kpt2, &kstat);
	 if (kpt2 == knum)
	    /* All edge intersections lie in one loop.  */

	    *jstat = 1;
	 else
	 {
	    /* Check if (one of) the curve(s) lies entirely in an
	       intersection curve found at an edge.  */

	    if (po1->iobj == SISLCURVE)
	    {
	       for (qpt1=vedge[0]->prpt[0]; qpt1!=SISL_NULL; qpt1=qpt1->pnext)
	       {
		  for (qpt2=vedge[0]->prpt[1]; qpt2!=SISL_NULL; qpt2=qpt2->pnext)
		  {
		     /* UJK, aug 93, oo-loop in sh6isconn, BEOrd20786. */
		     int is_conn,kcount;
		     is_conn = sh6isconnect(SISL_NULL, qpt1->ppt, qpt2->ppt);
		     for (kcount = 0;kcount<(*pintdat)->ipoint;kcount++)
			(*pintdat)->vpoint[kcount]->marker = 0;

		     if (is_conn) break;
		  }
		  if (qpt2 != SISL_NULL) break;
	       }

	       if (qpt1 != SISL_NULL && qpt2 != SISL_NULL) *jstat = 1;
	       else *jstat = 0;
	    }
	    if (*jstat != 1 && po2->iobj == SISLCURVE)
	    {
	       for (qpt1=vedge[1]->prpt[0]; qpt1!=SISL_NULL; qpt1=qpt1->pnext)
	       {
		  for (qpt2=vedge[1]->prpt[1]; qpt2!=SISL_NULL; qpt2=qpt2->pnext)
		  {
		     /* UJK, aug 93, oo-loop in sh6isconn, BEOrd20786. */
		     int is_conn,kcount;
		     is_conn = sh6isconnect(SISL_NULL, qpt1->ppt, qpt2->ppt);
		     for (kcount = 0;kcount<(*pintdat)->ipoint;kcount++)
			(*pintdat)->vpoint[kcount]->marker = 0;

		     if (is_conn) break;
		  }

		  if (qpt2 != SISL_NULL) break;
	       }

	       if (qpt1 != SISL_NULL && qpt2 != SISL_NULL) *jstat = 1;
	       else *jstat = 0;
	    }
	 }
      }
    }

  goto out;

  /* Error in subroutines      */

error:*jstat = kstat;
  s6err ("sh1762_s9con", *jstat, 0);
  goto out;

out:
  if (up != SISL_NULL)
    freearray (up);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9intercept (SISLObject * po1, SISLObject * po2, double aepsge,
		    int inmbpt, SISLIntpt * vintpt[], int *jstat)
#else
static void
sh1762_s9intercept (po1, po2, aepsge, inmbpt, vintpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     int inmbpt;
     SISLIntpt *vintpt[];
     int *jstat;
#endif
 /*
 *********************************************************************
 *
 *********************************************************************
 *
 * PURPOSE    : If there is less than two intersections between two
 *              objects, use rotated box tests to check if more
 *              intersections are possible.
 *
 *
 *
 * INPUT      : po1      - The first object to check.
 *              po2      - The second object to check.
 *              aepsge   - Geometrical resolution.
 *              inmbpt   - Number of intersections found on the edges.
 *              vintpt   - The intersections at the edges.
 *                         Dimension of pointer array is inmbpt.
 *              *jstat    - Flag
*                          = 202 : Complicated point-surface intersection
*                                  in 3D. Perform extra interception test.
*
 *
 *
 * OUTPUT     :
 *              jstat    - status messages
 *                           = 2     : No intersection..
 *                           = 1     : Coinside found.
 *                           = 0     : no coinside.
 *                           < 0     : error
 *
 *
 * METHOD     :
 *
 *
 * REFERENCES :
 *
 * CALLS      : s1221  - Evaluation of B-spline curve.
 *              sh1830  - Improved box test between curve and surface.
 *              sh1834  - Rotated box test.
 *              sh1839  - Improved box test between surface and object.
 *
 * WRITTEN BY : Arne Laksaa, SI, 89-05.
 * REWRITTEN BY : Vibeke Skytt, SI, 91-01.
 *
 *********************************************************************
 */
{
  int kstat = 0;		/* Status variable.               */
  int kdim;			/* Dimension of geometry space.   */
  int kleft = 0;		/* Parameter to curve evaluation. */
  int kleft2 = 0;               /* Parameter to evaluator.        */
  int incr, ind;		/* indexes and loop control       */
  int ratflag = 0;              /* Indicates if rational object.  */
  int kxintercept = (*jstat == 202);  /* Extra interception       */
  double tepsge;                /* Local tolerance in 1D box test. */
  double testpar[2];		/* Par val when treating help p.  */
  double trad;                  /* Radius of geometry object.     */
  double spar[2];               /* Parameter pair of surface.     */
  double scentre[3];            /* Centre of sphere of cylinder.  */
  double sder1[9];		/* Value and derivative of object.  */
  double sder2[9];		/* Pointer to value of second object.*/
  double snorm1[3];             /* Normal to first surface.       */
  double snorm2[3];             /* Normal to second surface.      */
  double splitgeom[16];         /* Matrix description of a sphere
				   or cylinder.                   */
  SISLSurf *qs1=SISL_NULL;           /* B-spline surface put into sphere
				   or cylinder equation.          */
  SISLSurf *qs2=SISL_NULL;           /* B-spline surface put into sphere
				   or cylinder equation.          */
  SISLCurve *qc=SISL_NULL;           /* B-spline curve put into sphere
				   equation.                      */
  SISLCurve *qc2=SISL_NULL;           /* B-spline curve put into sphere
				   equation.                      */
  SISLPoint *pp1=SISL_NULL;
  SISLObject *qobjs;		/* Pointer to surface object.     */
  SISLObject *qobjc;		/* Pointer to curve object.       */

  /*   long time_before;
  long time_used = 0;  */

  /* Test number of found intersection points.  */

  /* VSK, 01/93. if (inmbpt > 1 || inmbpt < 0)
    goto err128; */

  *jstat = 0;

  if (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE)
    {
       kdim = po1->s1->idim;

      /*      rotate_nmb++;
      time_before = clock(); */

       if (inmbpt == 0)
       {
	  /* No intersections at the edges.  */

	  if (xc % 2 == 0)
	  {
	     sh1839 (po1, po2, aepsge, &kstat);
	     if (kstat < 0)
		goto error;
	  }
	  else
	  {
	     sh1839 (po2, po1, aepsge, &kstat);
	     if (kstat < 0)
		goto error;
	  }
	  /*   time_used = clock() - time_before; */

	  if (kstat == 1)
	  {
	     if (xc % 2 == 0)
	     {
		sh6findsplit(po1->s1, po2->s1, aepsge, &kstat);
		if (kstat < 0) goto error;
	     }
	     else
	     {
		sh6findsplit(po2->s1, po1->s1, aepsge, &kstat);
		if (kstat < 0) goto error;
	     }
	  }

	  if (kstat == 0 || kstat == 2)
	  {
	     *jstat = 2;
	     goto out;
	  }

       }
       else
       {
	  /* Evaluate the surfaces in the first intersection point, and
	     use the partial derivatives in this point as rotation axises. */

	  s1421(po1->s1,1,vintpt[0]->epar,&kleft,&kleft2,sder1,snorm1,&kstat);
	  if (kstat < 0) goto error;

	  s1421(po2->s1,1,vintpt[0]->epar+2,&kleft,&kleft2,sder2,snorm2,&kstat);
	  if (kstat < 0) goto error;

	  sh1834(po1,po2,aepsge,kdim,sder1+kdim,sder1+2*kdim,&kstat);
	  if (kstat < 0) goto error;

	  if (kstat == 1 &&
	      fabs(s6ang(sder1+kdim,sder1+2*kdim,kdim) - PIHALF) > ANGULAR_TOLERANCE)
	  {
	     sh1834(po1,po2,aepsge,kdim,sder1+2*kdim,sder1+kdim,&kstat);
	     if (kstat < 0) goto error;
	  }

	  if (kstat == 1 &&
	      s6ang(sder1+kdim,sder2+kdim,kdim) > ANGULAR_TOLERANCE &&
	      s6ang(sder1+2*kdim,sder2+kdim,kdim) > ANGULAR_TOLERANCE)
	  {
	     sh1834(po1,po2,aepsge,kdim,sder2+kdim,sder2+2*kdim,&kstat);
	     if (kstat < 0) goto error;
	  }

	  if (kstat == 1 &&
	      fabs(s6ang(sder2+kdim,sder2+2*kdim,kdim) - PIHALF) > ANGULAR_TOLERANCE &&
	      s6ang(sder1+kdim,sder2+2*kdim,kdim) > ANGULAR_TOLERANCE &&
	      s6ang(sder1+2*kdim,sder2+2*kdim,kdim) > ANGULAR_TOLERANCE)

	  {
	     sh1834(po1,po2,aepsge,kdim,sder2+2*kdim,sder2+kdim,&kstat);
	     if (kstat < 0) goto error;
	  }

	  if (kstat == 0 || kstat == 2)
	  {
	     *jstat = 2;
	     goto out;
	  }

       }
    }
  else if ((po1->iobj == SISLSURFACE && po2->iobj == SISLCURVE) ||
	   (po1->iobj == SISLCURVE && po2->iobj == SISLSURFACE))
    {
      /*We test if intersection is possible using
	 rotated box tests. */

      if (po1->iobj == SISLSURFACE)
	{
	  qobjs = po1;
	  qobjc = po2;
	}
      else
	{
	  qobjs = po2;
	  qobjc = po1;
	}

      /* Perform improved box-test.  */

      /*   rotate_nmb++;
      time_before = clock(); */

      /* Improved box-test based on main tangent of curve and
	 main normal of surface.     */

      sh1830 (qobjs, qobjc, aepsge, &kstat);
      if (kstat < 0)
	goto error;

      if (kstat == 1)
	{
	  sh1839 (qobjs, qobjc, aepsge, &kstat);
	  if (kstat < 0)
	    goto error;
	}
      /*	time_used = clock() - time_before; */
      /*
      if (kstat == 1)
      {
	 Try to separate the objects by a sphere.

	 sh6sepgeom(qobjs->s1, qobjc->c1, aepsge, scentre, &trad, &kstat);
	 if (kstat < 0) goto error;

	  If kstat = 0 is returned, no splitting geometry is found,
	    and no further interception is to be tried.

	 if (kstat > 0)
	 {
	    The splitting geometry object is a sphere.
	       Make a matrix of dimension (idim+1)x(idim+1) describing a hyper
	       sphere as an implicit function.

	    s1321(scentre,trad,qobjc->c1->idim,1,splitgeom,&kstat);
	    if (kstat < 0) goto error;



	    * Put the description of the surface and the curve into the
	    * implicit equation for the sphere.
	    * ----------------------------------------------------------

	    ratflag = (qobjs->s1->ikind == 2 || qobjs->s1->ikind == 4) ? 1 : 0;
	    s1320(qobjs->s1,splitgeom,1,ratflag,&qs1,&kstat);
	    if (kstat < 0) goto error;

	    ratflag = (qobjc->c1->ikind == 2 || qobjc->c1->ikind == 4) ? 1 : 0;
	    s1370(qobjc->c1,splitgeom,qobjc->c1->idim,1,ratflag,&qc,&kstat);
	    if (kstat < 0) goto error;

	     Set up local tolerance.

	    tepsge = (double)2.0*trad*aepsge;

	     Make box of 1D surface.

	    sh1992su(qs1,2,tepsge,&kstat);
	    if (kstat < 0) goto error;

	     Make box of 1D curve.

	    sh1992cu(qc,2,tepsge,&kstat);
	    if (kstat < 0) goto error;

	     Check if the boxes overlap.

	    if (qs1->pbox->e2min[2][0] > qc->pbox->e2max[2][0] ||
		qs1->pbox->e2max[2][0] < qc->pbox->e2min[2][0])
	    {

	       No intersection is possible.

	       *jstat = 2;
	       goto out;
	    }
	    else kstat = 1;   Mark possibility of intersection.
	 }
	 else kstat = 1;   Mark possibility of intersection.
      } */

      if (kstat == 0 || kstat == 2)
	{
	  *jstat = 2;
	  goto out;
	}
    }
  else if (po1->iobj == SISLCURVE && po2->iobj == SISLCURVE &&
	   2*po1->c1->ik >= po1->c1->in && 2*po2->c1->ik >= po2->c1->in)
  {
     double spoint[3];  /* Point in splitting plane. */
     double snorm[3];   /* Normal to splitting plane. */
     double sn1[3], sn2[3];
     int ki;
     double t1, t2;
     int ksign;

     /* Find dimension of geometry space. */

     kdim = po1->c1->idim;
     if (kdim != po2->c1->idim)
	goto err106;

     if (inmbpt == 1)
     {
	/* One intersection point between two curves found. Find splitting
	   plane. */
	/* First allocate space for local arrays.  */
	/* NEWI, (ujk), Lets try to find a help point */

	incr = 0;
	if (DEQUAL (vintpt[0]->epar[0], po1->c1->et[po1->c1->ik - 1]))
	{
	   incr++;
	   testpar[0] = po1->c1->et[po1->c1->in];
	}
	else if (DEQUAL (vintpt[0]->epar[0], po1->c1->et[po1->c1->in]))
	{
	   incr++;
	   testpar[0] = po1->c1->et[po1->c1->ik - 1];
	}

	if (DEQUAL (vintpt[0]->epar[1], po2->c1->et[po2->c1->ik - 1]))
	{
	   incr++;
	   testpar[1] = po2->c1->et[po2->c1->in];
	}
	else if (DEQUAL (vintpt[0]->epar[1], po2->c1->et[po2->c1->in]))
	{
	   incr++;
	   testpar[1] = po2->c1->et[po2->c1->ik - 1];
	}

	if (incr == 2)
	   for (ind = 0; ind < vintpt[0]->no_of_curves; ind++)
	      if (sh6ishelp (vintpt[0]->pnext[ind]) &&
		  DEQUAL (vintpt[0]->pnext[ind]->epar[0], testpar[0]) &&
		  DEQUAL (vintpt[0]->pnext[ind]->epar[1], testpar[1]))
	      {
		 *jstat = 2;
		 goto out;
	      }

	/* Evaluate the curves in the intersection point.  */

	s1221 (po1->c1, 1, vintpt[0]->epar[0], &kleft, sder1, &kstat);
	if (kstat < 0)
	   goto error;

	s1221 (po2->c1, 1, vintpt[0]->epar[1], &kleft, sder2, &kstat);
	if (kstat < 0)
	   goto error;

	/* Normalize derivatives. */

	t1 = s6norm(sder1+kdim, kdim, sder1+kdim, &kstat);
	t2 = s6norm(sder2+kdim, kdim, sder2+kdim, &kstat);
	ksign = (s6scpr(sder1+kdim, sder2+kdim, kdim) >
		 DZERO) ? 1 : -1;
	for (ki=0; ki<kdim; ki++)
	{
	   /* sder1[kdim+ki] *= t2;
	   sder2[kdim+ki] *= t1; */
	   spoint[ki] = (double)0.5*(sder1[ki] + sder2[ki]);
	   sn1[ki] = (double)0.5*(sder1[kdim+ki]+sder2[kdim+ki]);
	}
	if (kdim == 2)
	{
	   snorm[0] = sn1[1];   /* KYS 5/7-94: normal corrected */
	   snorm[1] = -sn1[0];
	}
	else if (kdim == 3)
	{
	   s6crss(sder1+kdim, sder2+kdim, sn2);
	   s6crss(sn1, sn2, snorm);
	}
	(void)s6norm(snorm, kdim, snorm, &kstat);
	if (!kstat) kstat = 1;
	else kstat = 0;
     }
     else if (inmbpt == 0 && po1->c1->pdir->aang < ANGULAR_TOLERANCE &&
	      po2->c1->pdir->aang < ANGULAR_TOLERANCE  &&
	      s6ang(po1->c1->pdir->ecoef,po1->c1->pdir->ecoef,kdim) <
	      (double)10*ANGULAR_TOLERANCE)
     {
	double tpar2;
	SISLPoint *pt = SISL_NULL;
	double *s1, *s2, *s3, *s4;

	s1 = po1->c1->ecoef;
	s2 = po1->c1->ecoef+kdim*(po1->c1->in-1);
	s3 = po2->c1->ecoef;
	s4 = po2->c1->ecoef+kdim*(po2->c1->in-1);

	/* Evaluate midpoint of first curve. */

	/* tpar1 = (double)0.5*(po1->c1->et[po1->c1->ik-1] +
	   po1->c1->et[po1->c1->in]);
	   s1221 (po1->c1, 0, tpar1, &kleft, sder1, &kstat);
	   if (kstat < 0)
	   goto error; */
	if (MIN(s6dist(s1,s3,kdim),s6dist(s1,s4,kdim)) <
	    MIN(s6dist(s2,s3,kdim),s6dist(s2,s4,kdim)))
	   memcopy(sder1,s1,kdim,DOUBLE);
	else
	   memcopy(sder1,s2,kdim,DOUBLE);

	/* Find closest point on the other curve. */

	if ((pt = newPoint(sder1, kdim, 0)) == SISL_NULL) goto err101;

	/* tpar2 = (double)0.5*(po2->c1->et[po2->c1->ik-1] +
	   po2->c1->et[po2->c1->in]); */
	if (s6dist(s3,sder1,kdim) < s6dist(s4,sder1,kdim))
	   tpar2 = po2->c1->et[po2->c1->ik-1];
	else
	   tpar2 = po2->c1->et[po2->c1->in];
	s1771(pt, po2->c1, aepsge, po2->c1->et[po2->c1->ik-1],
	      po2->c1->et[po2->c1->in], tpar2, &tpar2, &kstat);

	if (pt) freePoint(pt);
	if (kstat < 0)
	   goto error;

	s1221 (po1->c1, 1, tpar2, &kleft, sder1, &kstat);
	if (kstat < 0)
	   goto error;
	s1221 (po2->c1, 1, tpar2, &kleft, sder2, &kstat);
	if (kstat < 0)
	   goto error;

	/* Let the splitting plane pass through the midpoint of the
	   points on the two curves and let the medium of the
	   axises of the direction cones of the curves lie in the
	   plane. */

	/* Normalize the tangents. */

	t1 = s6norm(sder1+kdim, kdim, sder1+kdim, &kstat);
	t2 = s6norm(sder2+kdim, kdim, sder2+kdim, &kstat);
	ksign = (s6scpr(sder1+kdim, sder2+kdim, kdim) >
		 DZERO) ? 1 : -1;
	for (ki=0; ki<kdim; ki++)
	{
	   /* sder1[kdim+ki] *= t2;
	   sder2[kdim+ki] *= t1; */
	   spoint[ki] = (double)0.5*(sder1[ki] + sder2[ki]);
	   sn1[ki] = (double)0.5*(sder1[kdim+ki] +
				  (double)ksign*sder2[kdim+ki]);
	}

	if (kdim == 3)
	{
	   s6crss(sder1+kdim, sder2+kdim, sn2);
	   s6crss(sn1, sn2, snorm);
	}
	else
	{
	   snorm[0] = sn1[1]; /* KYS 5/7-94: normal corrected */
	   snorm[1] = -sn1[0];
	}

	(void)s6norm(snorm, kdim, snorm, &kstat);
	if (!kstat) kstat = 1;
	else kstat = 0;
     }
     else kstat = 1;


     /* Try to intercept with the found plane. */

     if (kstat == 0)
     {
	/* nmb_rotated++; */
	sh1831(po1->c1, po2->c1, ksign, spoint, snorm, aepsge, &kstat);
	if (kstat < 0) goto error;
     }

     if (kstat == 0)
     {
	/* nmb_succ_rotated++; */
	*jstat = 2;
	goto out;
     }


     if (kstat == 1 && inmbpt == 0 && po1->c1->idim > 2)
     {
	/* kstat = 1; */			/* Make sure to subdivide further if there
						   is two curves and no intersection point. */
	/* Try to separate the objects by a sphere. */

	   if (xc % 2 == 0)
	   {
	      /* nmb_sep++; */
	      sh6sepcrv(po1->c1, po2->c1, aepsge, scentre, &trad, &kstat);
	      if (kstat < 0) goto error;
	   }
	   else
	   {
	      /* nmb_sep++; */
	      sh6sepcrv(po2->c1, po1->c1, aepsge, scentre, &trad, &kstat);
	      if (kstat < 0) goto error;
	   }

	/* If kstat = 0 is returned, no splitting geometry is found,
	   and no further interception is to be tried.  */

	if (kstat)
	{
	   /* The splitting geometry object is a sphere.
	      Make a matrix of dimension (idim+1)x(idim+1) describing a hyper
	      sphere as an implicit function.      	      */

	   /* nmb_try_sep++; */
	   s1321(scentre,trad,po1->c1->idim,1,splitgeom,&kstat);
	   if (kstat < 0) goto error;


	   /*
	   * Put the description of the surface and the curve into the
	   * implicit equation for the sphere.
	   * ----------------------------------------------------------
	   */

	   ratflag = (po1->c1->ikind == 2 || po1->c1->ikind == 4) ? 1 : 0;
	   s1370(po1->c1,splitgeom,po1->c1->idim,1,ratflag,&qc,&kstat);
	   if (kstat < 0) goto error;

	   ratflag = (po2->c1->ikind == 2 || po2->c1->ikind == 4) ? 1 : 0;
	   s1370(po2->c1,splitgeom,po2->c1->idim,1,ratflag,&qc2,&kstat);
	   if (kstat < 0) goto error;

	   /* Set up local tolerance. */

	   tepsge = (double)2.0*trad*aepsge;

	   /* Make box of 1D surface. */

	   sh1992cu(qc,2,tepsge,&kstat);
	   if (kstat < 0) goto error;

	   /* Make box of 1D curve. */

	   sh1992cu(qc2,2,tepsge,&kstat);
	   if (kstat < 0) goto error;

	   /* Check if the boxes overlap.  */

	   if (qc2->pbox->e2min[2][0] > qc->pbox->e2max[2][0] ||
	       qc2->pbox->e2max[2][0] < qc->pbox->e2min[2][0])
	   {

	      /* No intersection is possible.  */

	      /* numb_succ_sep++; */
	      *jstat = 2;
	      goto out;
	   }
	   else kstat = 1;  /* Mark possibility of intersection.  */
	}
	else kstat = 1;
     }
     else kstat = 1;
  }
  else if ((po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT &&
	   po2->p1->idim == 2) ||
	   (po2->iobj == SISLSURFACE && po1->iobj == SISLPOINT &&
	   po1->p1->idim == 2))
  {
     /* Compute the mid-parameter value of the surface. First set
	pointer to the surface.  */

     if (po1->iobj == SISLSURFACE) qs1 = po1->s1;
     else qs1 = po2->s1;

     spar[0] = (double)0.5*(qs1->et1[qs1->ik1-1] + qs1->et1[qs1->in1]);
     spar[1] = (double)0.5*(qs1->et2[qs1->ik2-1] + qs1->et2[qs1->in2]);

     /* Evaluate the surface in the midpoint. */

     s1421(qs1, 1, spar, &kleft, &kleft2, sder1, snorm1, &kstat);
     if (kstat < 0) goto error;

     if (s6ang(sder1+2, sder1+4, 2) < ANGULAR_TOLERANCE)
     {
	spar[0] = (double)0.5*(sder1[2]+sder1[4]);
	spar[1] = (double)0.5*(sder1[3]+sder1[5]);
	sh1834(po1, po2, aepsge, 2, spar, sder1+4, &kstat);
	if (kstat < 0) goto error;
	if (kstat == 5) kstat = 0;   /* No 45 degree testing for rotated
					box test meens no danger of
					intersection point near corner that
					is not caught by the box test. */
     }
     else kstat = 1;

     qs1 = SISL_NULL;     /* Make sure that the input surface is not freed. */
  }
  else if (((po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT &&
	   po2->p1->idim == 3) ||
	   (po2->iobj == SISLSURFACE && po1->iobj == SISLPOINT &&
	   po1->p1->idim == 3)) && kxintercept && xc > 7 && xc % 2 == 0)
  {
    if (po1->iobj == SISLSURFACE) 
      {
	qs1 = po1->s1;
	pp1 = po2->p1;
      }
    else 
      {
	qs1 = po2->s1;
	pp1 = po1->p1;
      }
    kdim = qs1->idim;

    if (qs1->in1 > qs1->ik1 || qs1->in2 > qs1->ik2)
      kstat = 1;
    else
      {
	int ind1, ind2, ind3;
	int kpt = 0;
	int kcrv = 0;
	double *spar = SISL_NULL;
	SISLIntcurve **ucurve = SISL_NULL;
	double eps = 0.001*aepsge;

	/* Find the closest points between the surface and the point */
	s1954(qs1, pp1->ecoef, qs1->idim, 0.0, eps, &kpt, &spar,
	      &kcrv, &ucurve, &kstat);
	if (kstat < 0)
	  goto error;


	/* Test distance between the closest points on the surface and
	   the point  */
	for (ind1=0; ind1<kpt; ind1++)
	  {
	    s1421(qs1, 0, spar+2*ind1, &kleft, &kleft2, sder1, 
		  snorm1, &kstat);
	    if (s6dist(pp1->ecoef, sder1, kdim) <= aepsge)
	      break;
	  }

	for (ind2=0; ind2<kcrv; ind2++)
	  {
	    for (ind3=0; ind3<ucurve[ind2]->ipoint; ind3++)
	      {
		s1421(qs1, 0, ucurve[ind2]->epar1+2*ind3, &kleft, &kleft2, 
		      sder1, snorm1, &kstat);
		if (s6dist(pp1->ecoef, sder1, kdim) <= aepsge)
		  break;
	      }
	    if (ind3 < ucurve[ind2]->ipoint)
	      break;
	  }

	if (ind1 < kpt || ind2 < kcrv)
	  kstat = 1;
	else kstat = 0;

	/* fprintf(stdout,"%7.13f %7.13f %7.13f %7.13f \n",qs1->et1[0],
		qs1->et1[qs1->in1],qs1->et2[0],qs1->et2[qs1->in2]);
	fprintf(stdout,"Point-srf : kstat = %d\n",kstat); */

	if (spar)
	  freearray(spar);
	if (ucurve)
	  freeIntcrvlist(ucurve, kcrv);
      }
    qs1 = SISL_NULL;
  }  
  else kstat = 1;


  *jstat = (kstat == 0 || kstat == 2) ? 2 : 0;
  goto out;

  /* Error in scratch allocation.  */

  err101: *jstat = -101;
  goto out;

  /* Error in input. Confliciting dimensions.  */

err106:*jstat = -106;
  goto out;

  /* Wrong number of intersection points on edge.  */

  /* err128:*jstat = -128;
  goto out; */

  /* Error in lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
   /* Free scratch used by 1D surfaces. */

   if (qs1 != SISL_NULL) freeSurf(qs1);
   if (qs2 != SISL_NULL) freeSurf(qs2);
   if (qc != SISL_NULL) freeCurve(qc);
   if (qc2 != SISL_NULL) freeCurve(qc2);

  /*	rotate_box_time += time_used;	 */
  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9coincide (SISLObject * po1, SISLObject * po2, double aepsge,
		   int inmbpt, SISLIntpt * vintpt[], int *jstat)
#else
static void
sh1762_s9coincide (po1, po2, aepsge, inmbpt, vintpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     int inmbpt;
     SISLIntpt *vintpt[];
     int *jstat;
#endif
 /*
 *********************************************************************
 *
 *********************************************************************
 *
 * PURPOSE    : Check if one part of an object coincide in one part
 *              of an other object. Both object must intersect in two
 *              ends/edges.
 *
 *
 * INPUT      : po1      - The first SISLObject to check.
 *              po2      - The second SISLObject to check.
 *              aepsge   - Geometrical resolution.
 *              inmbpt   - Number of intersection points on edges (=2).
 *              vintpt   - Intersection points found on edges.
 *                         Dimension of pointer array is inmbpt.
 *              vedge[]  - SISLEdge intersection.
 *
 *
 * INPUT/OUTPUT:pintdat  - Intersection data.
 *
 *
 * OUTPUT     :
 *              jstat    - status messages
 *                           = 2     : No intersection..
 *                           = 1     : Coinside found.
 *                           = 0     : no coinside.
 *                           < 0     : error
 *
 *
 * METHOD     :
 *
 *
 * REFERENCES :
 *
 * CALLS      : s1221  - Evaluation of curve.
 *              s1421  - Evaluation of surface.
 *              s1785  - Test coincidence between curve and surface.
 *              s1786  - Test coincidence between two curves.
 *              sh6getlist - Check if two intersection points are already
 *                           connected.
 *
 *
 * WRITTEN BY : Arne Laksaa, SI, 89-05.
 * REVISED BY : Vibeke Skytt, SI, 91-01.
 * REVISED BY : Michael Floater, SI, 91-08.
 *                   Removed angle test.
 * REVISED BY : Vibeke Skytt, SINTEF Oslo, 94-11. Coincidence marching
 *                                                for 2D point - surface.
 *
 *********************************************************************
 */
{
  int kstat = 0;		/* Status variable.                           */
  int kdim;			/* Dimension of geometry space.               */
  int kcur;			/* Indicates the curve in curve-surface
				   intersection.                              */
  int kn;			/* Counter.                                   */
  int kleft1 = 0, kleft2 = 0;	/* Parameters used in evaluation.           */
  int kind1,kind2;              /* Dummy parameters to sh6getlist.            */
  double tang;			/* Angle between vectors.                     */
  double *snorm;		/* Pointer to surface normal.                 */
  double *sder1 = SISL_NULL;		/* Array containing position etc. of objects. */
  double *sder2;		/* Pointer to position of second object.      */
  SISLSurf *qs;			/* Pointer to surface.                        */
  SISLCurve *qc;		/* Pointer to curve.                          */
  SISLPoint *qp;

  if (inmbpt != 2)
    goto err128;

  if ((po1->iobj == SISLSURFACE && po2->iobj == SISLCURVE) ||
      (po1->iobj == SISLCURVE && po2->iobj == SISLSURFACE))
    {
      /* We test coincidence for curve-surface. */

       /* VSK, 10.92. First check if the points are already connected.     */

       sh6getlist(vintpt[0],vintpt[1],&kind1,&kind2,&kstat);
       if (kstat < 0) goto error;

       if (kstat == 0)
       {
	  /* The points are already connected.  */

	  *jstat = 1;
	  goto out;
       }

      if (po1->iobj == SISLSURFACE)
	{
	  qs = po1->s1;
	  qc = po2->c1;
	  kcur = 0;
	}
      else
	{
	  qs = po2->s1;
	  qc = po1->c1;
	  kcur = 1;
	}

      /* Allocate space for local arrays.  */

      if ((sder1 = newarray (6 * qc->idim, double)) == SISL_NULL)
	goto err101;
      sder2 = sder1 + 2 * qc->idim;
      snorm = sder2 + 3 * qc->idim;

      for (kn = 0; kn < 2; kn++)
	{
	  /* We have to test if the curve and the surface
	     have coinciding derivatives in intersection ponts. */

	  s1221 (qc, 1, vintpt[kn]->epar[(kcur ? 0 : 2)], &kleft1, sder1, &kstat);
	  if (kstat < 0)
	    goto error;

	  s1421 (qs, 1, vintpt[kn]->epar + kcur, &kleft1, &kleft2, sder2, snorm, &kstat);
	  if (kstat < 0)
	    goto error;
	  else if (kstat > 0)
	    {
	      /* Singular point.  */

	      *jstat = 0;
	      goto out;
	    }
/*
	  tang = s6ang (sder1 + qc->idim, snorm, qc->idim);

	  if (PIHALF - tang > ANGULAR_TOLERANCE)
	    {
	      *jstat = 0;
	      goto out;
	    }
*/
	}
      /* Removed the angle test. M.F. 30/8/91.  */
      /* If the first derivatives are equal we call a routine
	 to test further for coincidence.  */

      s1785 (qc, qs, aepsge, vintpt[0]->epar, vintpt[1]->epar, kcur, &kstat);
      if (kstat < 0)
	goto error;

    }
  else if (po1->iobj == SISLCURVE && po2->iobj == SISLCURVE)
    {
      kdim = po1->c1->idim;
      if (kdim != po2->c1->idim)
	goto err106;

      /* Test coincidence between two curves. First allocate
	 space for local arrays.  */

      if ((sder1 = newarray (8 * kdim, double)) == SISL_NULL)
	goto err101;
      sder2 = sder1 + 4 * kdim;

      /* Evaluate the curves in the first intersection point.  */

      s1221 (po1->c1, 1, vintpt[0]->epar[0], &kleft1, sder1, &kstat);
      if (kstat < 0)
	goto error;

      s1221 (po2->c1, 1, vintpt[0]->epar[1], &kleft1, sder2, &kstat);
      if (kstat < 0)
	goto error;

      /* Evaluate the curves in the second intersection point.  */

      s1221 (po1->c1, 1, vintpt[1]->epar[0], &kleft1, sder1 + (2 * kdim), &kstat);
      if (kstat < 0)
	goto error;

      s1221 (po2->c1, 1, vintpt[1]->epar[1], &kleft2, sder2 + (2 * kdim), &kstat);
      if (kstat < 0)
	goto error;

      /* Test if the curves are parallel in the endpoints. */

      tang = s6ang (sder1 + kdim, sder2 + kdim, kdim);

      if (tang > ANGULAR_TOLERANCE)
	{
	  *jstat = 0;
	  goto out;
	}

      tang = s6ang (sder1 + (3 * kdim), sder2 + (3 * kdim), kdim);

      if (tang > ANGULAR_TOLERANCE)
	{
	  *jstat = 0;
	  goto out;
	}

      s1786 (po1->c1, po2->c1, aepsge, vintpt[0]->epar, vintpt[1]->epar, &kstat);
      if (kstat < 0)
	goto error;

    }
  else if ((po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT &&
	    po2->p1->idim >= 2) ||
	   (po2->iobj == SISLSURFACE && po1->iobj == SISLPOINT &&
	    po1->p1->idim >= 2))
  {
     if (po1->iobj == SISLSURFACE)
     {
	qs = po1->s1;
	qp = po2->p1;
     }
     else
     {
	qs = po2->s1;
	qp = po1->p1;
     }

     /* Allocate space for local arrays.  */

     if ((sder1 = newarray (7 * qs->idim, double)) == SISL_NULL)
	goto err101;
     sder2 = sder1 + 3 * qs->idim;
     snorm = sder2 + 3 * qs->idim;

     /* Evaluate the surface in the intersection points at the edges. */

     s1421 (qs, 1, vintpt[0]->epar, &kleft1, &kleft2, sder1, snorm, &kstat);
     if (kstat < 0)
	goto error;

     s1421 (qs, 1, vintpt[1]->epar, &kleft1, &kleft2, sder2, snorm, &kstat);
     if (kstat < 0)
	goto error;

     /* Test if this is a singular situation. */

     if (s6ang(sder1+qs->idim, sder1+2*qs->idim, qs->idim) <= 
	 ANGULAR_TOLERANCE &&
	 s6ang(sder2+qs->idim, sder2+2*qs->idim, qs->idim) <= 
	 ANGULAR_TOLERANCE)
     {
	/* Perform marching to check if there is coincidence between
	   the intersection points. */

	 /* fprintf(stdout,"Try coincidence marching \n"); 
	 fprintf(stdout,"%7.13f %7.13f %7.13f %7.13f \n",qs->et1[0],
		 qs->et1[qs->in1],qs->et2[0],qs->et2[qs->in2]); */

	s1789(qp, qs, aepsge, vintpt[0]->epar, vintpt[1]->epar, &kstat);
	if (kstat < 0) goto error;
	 /* fprintf(stdout,"kstat = %d \n",kstat); */
     }
     else
	kstat = 0;   /* No coincidence. */
  }

  *jstat = kstat;
  goto out;

  /* Error in scratc allocation.  */

err101:*jstat = -101;
  goto out;

  /* Error in input. Conflicting dimensions.  */

err106:*jstat = -106;
  goto out;

  /* Wrong number of edge intersections found.  */

err128:*jstat = -128;
  goto out;

  /* Error in lower level routine.  */

error:*jstat = kstat;
  goto out;

out:

  /* Free scratch occupied by local array.  */

  if (sder1 != SISL_NULL)
    freearray (sder1);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9toucharea (SISLObject * po1, SISLObject * po2, double aepsge,
		    int inmbpt, SISLIntpt * vintpt[], int *jstat)
#else
static void
sh1762_s9toucharea (po1, po2, aepsge, inmbpt, vintpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     int inmbpt;
     SISLIntpt *vintpt[];
     int *jstat;
#endif
 /*
 *********************************************************************
 *
 *********************************************************************
 *
 * PURPOSE    : Check if the intersection of two surfaces
 *              is a surface (a "touch area").
 *              Both objects must intersect in at least
 *              three ends or edges and all the intersection
 *              points must be linked in one closed list.
 *
 *
 * INPUT      : po1      - The first SISLObject to check.
 *              po2      - The second SISLObject to check.
 *              aepsge   - Geometrical resolution.
 *              inmbpt   - Number of intersection points on edges (>=3).
 *              vintpt   - Intersection points found on edges.
 *                         Dimension of pointer array is inmbpt.
 *                         VSK. The array is sorted.
 *              vedge[]  - SISLEdge intersection.
 *
 *
 * INPUT/OUTPUT:pintdat  - Intersection data.
 *
 *
 * OUTPUT     :
 *              jstat    - status messages
 *                           = 2     : No intersection..
 *                           = 1     : Coinside found.
 *                           = 0     : no coinside.
 *                           < 0     : error
 *
 *
 * METHOD     :
 *
 *
 * REFERENCES :
 *
 * CALLS      : s1421  - Evaluation of surface.
 *              s1773  - Closest point to a surface.
 *              s1436  - Pick a const. v curve from surface.
 *              s1437  - Pick a const. u curve from surface.
 *              sh1784 - Test coincidence between curve and surface.
 *
 *
 * WRITTEN BY : Vibeke Skytt, 9403.
 *              Note that this verions of coincidence testing between
 *              surfaces are a prelimenary one, not a final version.
 * Changed by : Per OEyvind Hvidsten, SINTEF, 11-94
 *              Added a goto out; at end (thus skipping err101).
 *
 *********************************************************************
 */
{
   int kstat = 0;         /* Local status variable. */
   int kntest1, kntest2;  /* Number of locations to test coincidence in
			     both parameter directions.                 */
   double tint1, tint2;   /* Parameter interval between testing spots.  */
   int ki,kj;             /* Counters.                                  */
   int kdim = po1->s1->idim; /* Dimension of geometry space.            */
   double spar[2];        /* Parameter of testing spot.                 */
   double sder1[3];       /* Position of first surface.                 */
   double sder2[3];       /* Position of second surface.                */
   double snorm1[3], snorm2[3];  /* Dummy normals of surface.           */
   int kleft11 = 0, kleft12 = 0; /* Pointers into knot arrays of surface. */
   int kleft21 = 0, kleft22 = 0; /* Pointers into knot arrays of surface. */
   int kn11 = po1->s1->in1;
   int kn12 = po1->s1->in2;
   int kk11 = po1->s1->ik1;
   int kk12 = po1->s1->ik2;
   double *st11 = po1->s1->et1;
   double *st12 = po1->s1->et2;
   int kn21 = po2->s1->in1;
   int kn22 = po2->s1->in2;
   int kk21 = po2->s1->ik1;
   int kk22 = po2->s1->ik2;
   double *st21 = po2->s1->et1;
   double *st22 = po2->s1->et2;
   SISLPoint *pt = SISL_NULL;       /* Point in point surface iteration.       */
   double sstart[2], send[2];  /* Parameter boundaries of second surface. */
   double spar2[2];            /* Parameter value of second surface.      */

   /* Set number of locations to test coincidence. */

   kntest1 = 30*(kn11 - kk11 + 1);
   kntest2 = 30*(kn12 - kk12 + 1);
   tint1 = (st11[kn11] - st11[kk11-1])/(double)(kntest1+1);
   tint2 = (st12[kn12] - st12[kk12-1])/(double)(kntest2+1);

   /* Set parameter boundaries and midpoint of second surface. */

   sstart[0] = st21[kk21-1];
   sstart[1] = st22[kk22-1];
   send[0] = st21[kn21];
   send[1] = st22[kn22];
   spar2[0] = (double)0.5*(sstart[0] + send[0]);
   spar2[1] = (double)0.5*(sstart[1] + send[1]);

   for (spar[0]=st11[kk11-1]+tint1, ki=0; ki<kntest1; ki++, spar[0]+=tint1)
   {
      for (spar[1]=st12[kk12-1]+tint2, kj=0; kj<kntest2; kj++, spar[1]+=tint2)
      {
	 /* Evaluate first surface. */

	 s1421(po1->s1, 0, spar, &kleft11, &kleft12, sder1, snorm1, &kstat);
	 if (kstat < 0) goto error;

	 /* Find closest point on the other surface. */

	 if ((pt =  newPoint(sder1, kdim, 0)) == SISL_NULL) goto err101;

	 s1773(pt, po2->s1, aepsge, sstart, send, spar2, spar2, &kstat);
	 if (kstat < 0) goto error;

	 /* Evalutate second surface. */

	 s1421(po2->s1, 0, spar2, &kleft21, &kleft22, sder2, snorm2, &kstat);
	 if (kstat < 0) goto error;

	 if (pt != SISL_NULL) freePoint(pt);
	 pt = SISL_NULL;

	 /* Check distance between the closest points. */

	 if (s6dist(sder1, sder2, kdim) > aepsge) break;  /* Not a coincidence.*/
      }
      if (kj < kntest2) break;  /* Not a coincidence. */
   }

   *jstat = (ki==kntest1 && kj==kntest2) ? 1 : 0;
   goto out;

   err101 : *jstat = -101;    /* Error in scratch allocation. */
   goto out;

   error : *jstat = kstat;    /* Error in lower level function. */
   goto out;

   out:
      if (pt != SISL_NULL) freePoint(pt);

      return;
}


#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9edgsscon (SISLEdge * vedge[], SISLSurf * ps1, SISLSurf * ps2,
		   SISLIntdat * rintdat, int isimple, double aepsge,
		   int *jstat)
#else
static void
sh1762_s9edgsscon (vedge, ps1, ps2, rintdat, isimple, aepsge, jstat)
     SISLEdge *vedge[];
     SISLSurf *ps1;
     SISLSurf *ps2;
     SISLIntdat *rintdat;
     int isimple;
     double aepsge;
     int *jstat;
#endif
 /*
 *********************************************************************
 *
 *********************************************************************
 *
 * PURPOSE    : To connect intersection points on edges on two surfaces
 *              into curves.
 *
 *
 *
 * INPUT      : vedge[]  - SISLEdge intersection.
 *              ps1      - First surface in intersection.
 *              ps2      - Second surface in intersection.
 *              isimple  - If value 1, we use a special connecting
 *                         strategy for 4 points.
 *              aepsge   - Geometry tolerance.
 *
 * OUTPUT     : rintdat  - Intersection dates to be updated.
 *              jstat    - status messages
 *                           = 1     : Not a simple case.
 *                           = 0     : OK.
 *                           < 0     : error
 *
 *
 * METHOD     :
 *
 *
 * REFERENCES :
 *
 *
 * WRITTEN BY : Arne Laksaa, SI, 89-06.
 *
 *********************************************************************
 */
{
  int kstat,kstat1,kstat2;
  int kmarch = 0;               /* Indicates if marching is to be done. */
  int *ldir = SISL_NULL;		/* Local array containing one of the statusvalues for
				 * each point:
				 *  0 - The intersect.curve is parallel to one
				 *      parameter direction.
				 *  1 - The intersect.curve has direction into the
				 *      domain.
				 * -1 - The intersect.curve has direction out of the
				 *      domain.
				 *  2 - The point is singulear.
				 * 10 - The intersect.curve touch one corner of the
				 *      domain.
				 * --------------------------------------------------
				 */

  int lant[2];
  unsigned char *edg = SISL_NULL;
  double *sval1 = SISL_NULL;
  SISLIntpt **uipt = SISL_NULL;
  SISLIntpt **uinewpt = SISL_NULL;

  double *spar = SISL_NULL;		/* Local array with parameter values used as
				   input to s9surmarch. */
  int *lperm = SISL_NULL;		/* Local permutation array after sorting
				   input points to s9surmarch.*/
  int *lpermdir = SISL_NULL;		/* Local array with status values used as
				   input to s9surmarch. */
  int lstatus[4];		/* Local array containing the possible status
				   constants -1,1,0,2. */
  int lnumb[4];			/* Local array containing the number of points
				   with status lstatus. */

  double *sparout = SISL_NULL;	/* Local array with parameter values used as
				   output from s9surmarch.*/
  int *lpar = SISL_NULL;		/* Local array containing the connection information
				   from s9surmarch.*/
  int kpoints;			/* Local integer containing number of points returned
				   from s9surmarch.*/
  double *nullp = SISL_NULL;
  int klist1, klist2;		/* List index in iintpt.   */

  double tdist;                 /* Distance between surfaces in point.   */
  double tref;                  /* Reference value.                      */
  double spos[4];               /* Parameter value of singular point.    */
  double start[4];              /* Start parameter to iteration.         */
  double slimit[8];             /* Limits to the parameter areas.        */
  SISLIntpt *qsing=SISL_NULL;        /* Singular intersection point.          */

  /* Experiment UJK, sept 92 (BEOrd12754) */
    double tolpar = (double) 0.001;

    /* double tolpar = (double) 0.01; */

  if (ps1->idim != 3 || ps2->idim != 3)
    goto err200;

  *jstat = 0;

  if (vedge[0] == SISL_NULL)
    lant[0] = 0;
  else
    lant[0] = vedge[0]->ipoint;

  if (vedge[1] == SISL_NULL)
    lant[1] = 0;
  else
    lant[1] = vedge[1]->ipoint;


  if (lant[0] + lant[1] > 1)
    {
      int kn1, kn, ki, kj, kv, kant, klfs, klft, kdir, kpar;
      double *sval2, *snorm1, *snorm2, *stang, *sdec1, *sdec2;
      SISLPtedge *qpt;
      double tang;

      kant = lant[0] + lant[1];

      /* Allocate array of pointers to the points. */

      if ((uipt = newarray (kant, SISLIntpt *)) == SISL_NULL)
	goto err101;
      if ((edg = new0array (kant, unsigned char)) == SISL_NULL)
	goto err101;
      if ((ldir = new0array (kant, int)) == SISL_NULL)
	goto err101;
      if ((sval1 = newarray (33 * kant, double)) == SISL_NULL)
	goto err101;
      sval2 = sval1 + 9 * kant;
      snorm1 = sval2 + 9 * kant;
      snorm2 = snorm1 + 3 * kant;
      stang = snorm2 + 3 * kant;
      sdec1 = stang + 3 * kant;
      sdec2 = sdec1 + 3 * kant;


      /* UPDATE (ujk) Main point vs helppoints covered ? */
      /* Update the arrays. */

      for (kn1 = 0, kn = 0; kn < 2; kn++)
	if (lant[kn] > 0)
	  for (kj = 0; kj < vedge[kn]->iedge; kj++)
	    for (qpt = vedge[kn]->prpt[kj]; qpt != SISL_NULL; qpt = qpt->pnext)
	      {
		for (ki = 0; ki < kn1; ki++)
		  {
		    if (qpt->ppt == uipt[ki])
		      break;
		  }
		if (ki == kn1)
		  uipt[kn1++] = qpt->ppt;

		edg[ki] |= 1 << (vedge[kn]->iedge * kn + kj);
	      }

      if (kn1 > 1)
	for (ki = 0; ki < kn1; ki++)
	  {
	    kv = 3 * ki;
	    kn = 9 * ki;
	    /* UPDATE (ujk) : getgeom !!  */

	    klfs = klft = 0;
	    s1421 (ps1, 1, uipt[ki]->epar, &klfs, &klft, sval1 + kn, snorm1 + kv,
		   &kstat);
	    if (kstat < 0)
	      goto error;
	    else if (kstat > 0)
	      {
		ldir[ki] = 2;
		continue;
	      }

	    klfs = klft = 0;
	    s1421 (ps2, 1, uipt[ki]->epar + 2, &klfs, &klft, sval2 + kn, snorm2 + kv,
		   &kstat);
	    if (kstat < 0)
	      goto error;
	    else if (kstat > 0)
	      {
		ldir[ki] = 2;
		continue;
	      }

	    tang = s6ang (snorm1 + kv, snorm2 + kv, 3);
	    if (tang < REL_PAR_RES)
	      /*	    if (tang < ANGULAR_TOLERANCE) */
	      {
		ldir[ki] = 2;
		continue;
	      }

	    s6crss (snorm1 + kv, snorm2 + kv, stang + kv);

	    s6decomp (stang + kv, sdec1 + kv, sval1 + kn + 3, sval1 + kn + 6,
		      snorm1 + kv, &kstat);
	    if (kstat < 0)
	      goto error;
	    else if (kstat > 0)
	      {
		ldir[ki] = 2;
		continue;
	      }


	    s6decomp (stang + kv, sdec2 + kv, sval2 + kn + 3, sval2 + kn + 6,
		      snorm2 + kv, &kstat);
	    if (kstat < 0)
	      goto error;
	    else if (kstat > 0)
	      {
		ldir[ki] = 2;
		continue;
	      }

	    for (kpar = 1, kj = 0; kj < 8; kj++)
	      if ((edg[ki] & 1 << kj) == 1 << kj)
		{
		  switch (kj)
		    {
		    case 0:
		      tang = s6ang (stang + kv, sval1 + kn + 3, 3);
		      kdir = (sdec1[kv + 1] > DZERO ? 1 : -1);
		      break;
		    case 4:
		      tang = s6ang (stang + kv, sval2 + kn + 3, 3);
		      kdir = (sdec2[kv + 1] > DZERO ? 1 : -1);
		      break;
		    case 1:
		      tang = s6ang (stang + kv, sval1 + kn + 6, 3);
		      kdir = (sdec1[kv] > DZERO ? -1 : 1);
		      break;
		    case 5:
		      tang = s6ang (stang + kv, sval2 + kn + 6, 3);
		      kdir = (sdec2[kv] > DZERO ? -1 : 1);
		      break;
		    case 2:
		      tang = s6ang (stang + kv, sval1 + kn + 3, 3);
		      kdir = (sdec1[kv + 1] > DZERO ? -1 : 1);
		      break;
		    case 6:
		      tang = s6ang (stang + kv, sval2 + kn + 3, 3);
		      kdir = (sdec2[kv + 1] > DZERO ? -1 : 1);
		      break;
		    case 3:
		      tang = s6ang (stang + kv, sval1 + kn + 6, 3);
		      kdir = (sdec1[kv] > DZERO ? 1 : -1);
		      break;
		    case 7:
		      tang = s6ang (stang + kv, sval2 + kn + 6, 3);
		      kdir = (sdec2[kv] > DZERO ? 1 : -1);
		    }

		  if (tang < tolpar)
		    kdir = 0;
		  /* if (tang < REL_PAR_RES)
		    kdir = 0; */
		  /*  if (tang < ANGULAR_TOLERANCE) kdir = 0;*/
		  if (kdir == 0)
		    kpar = 0;
		  else if (ldir[ki] != kdir)
		    {
		      if (ldir[ki] == 0)
			ldir[ki] = kdir;
		      else
			{
			  ldir[ki] = 10;
			  break;
			}
		    }
		}
	    if (kpar == 0 && ldir[ki] != 10)
	      ldir[ki] = 0;
	  }


      /* When only two points, check if they are connected. */
      if (kn1 == 2)
	{
	  sh6getlist (uipt[0], uipt[1], &klist1, &klist2, &kstat);
	  if (kstat == 0)
	    kn1 = 0;
	}


      for (kv = ki = 0; ki < kn1; ki++)
	if (ldir[ki] < 10)
	  kv++;

      if (kv == 1 && kn1 > 1 )
      {
	 int i,j;

	 if (kn1 > 2)
	 {
	    *jstat = 0;
	    for (i=0;i<kn1;i++)
	       if (ldir[i] == 1 || ldir[i] == -1)
		  *jstat = 1;
	    goto out;
	 }


	 if (ldir[0] == 1 || ldir[0] == -1)
	 {
	    i =  0;
	    j =  1;
	 }
	 else if (ldir[1] == 1 || ldir[1] == -1)
	 {
	    i =  1;
	    j =  0;
	 }
	 else						i = -1;

	 if (i >=0)
	 {
	    /* We have one point, and the direction is in/out.
	       The other point is touching a corner it is probably
	       an error and the point must be close to tangential.*/
	    if (ldir[i] == -1)
	    {
	       int k=j;
	       j  =  i;
	       i  =  k;
	    }

	    sh6tomain (uipt[i], &kstat);
	    sh6tomain (uipt[j], &kstat);

	    sh6idcon (&rintdat, &uipt[i], &uipt[j], &kstat);
	    if (kstat < 0)
	       goto error;
	    sh6setdir (uipt[i], uipt[j], &kstat);
	    if (kstat < 0)
	       goto error;
	    *jstat = 0;
	 }
      }
      else if (kv < 2)
	/* Less than two points, it's a simple case. */
	*jstat = 0;

      else if (kv > 8)
	/* More than eight points, it's not a simple case. */
	*jstat = 1;

      else if (kv == 4 && isimple > 0)
	{
	  /* Connect point in to closest point out. */
	  int kinn = 0, kout = 0;
	  int ki, ki1, ki2, ko1, ko2;

	  for (ki = 0; ki < kn1; ki++)
	    if (ldir[ki] == 1)
	      {
		if (kinn == 0)
		  ki1 = ki;
		else
		  ki2 = ki;
		kinn++;
	      }
	    else if (ldir[ki] == -1)
	      {
		if (kout == 0)
		  ko1 = ki;
		else
		  ko2 = ki;
		kout++;
	      }

	  if (kinn == 2 && kout == 2)
	    {
	      double tdist1, tdist2, tdir1, tdir2;

	      tdist1 = s6dist (sval1 + (ki1 * 9), sval1 + (ko1 * 9), 3);
	      tdist2 = s6dist (sval1 + (ki1 * 9), sval1 + (ko2 * 9), 3);

	      tdir1 = s6scpr (stang + (ki1 * 3), stang + (ko1 * 3), 3);
	      tdir2 = s6scpr (stang + (ki1 * 3), stang + (ko2 * 3), 3);

	      if (tdir2 < DZERO || (tdir1 >= DZERO && tdist1 <= tdist2))
		{
		  /* We can connect the points. */
		  /* UPDATE (ujk) : do we need to make conistency
	             between the points ? or can it be delayed until
	             s6idlis ?? (this comments is valid for all
	             connects in sscon and pscon)*/

		  sh6tomain (uipt[ki1], &kstat);
		  sh6tomain (uipt[ko1], &kstat);
		  sh6tomain (uipt[ki2], &kstat);
		  sh6tomain (uipt[ko2], &kstat);

		  sh6idcon (&rintdat, &uipt[ki1], &uipt[ko1], &kstat);
		  if (kstat < 0)
		    goto error;
		  /* Newi (ujk) Set in/out direction */
		  sh6setdir (uipt[ki1], uipt[ko1], &kstat);
		  if (kstat < 0)
		    goto error;

		  sh6idcon (&rintdat, &uipt[ki2], &uipt[ko2], &kstat);
		  if (kstat < 0)
		    goto error;
		  /* Newi (ujk) Set in/out direction */
		  sh6setdir (uipt[ki2], uipt[ko2], &kstat);
		  if (kstat < 0)
		    goto error;

		}
	      else
		{
		  /* We can connect the points. */

		  sh6tomain (uipt[ki1], &kstat);
		  sh6tomain (uipt[ko1], &kstat);
		  sh6tomain (uipt[ki2], &kstat);
		  sh6tomain (uipt[ko2], &kstat);

		  sh6idcon (&rintdat, &uipt[ki1], &uipt[ko2], &kstat);
		  if (kstat < 0)
		    goto error;
		  /* Newi (ujk) Set in/out direction */
		  sh6setdir (uipt[ki1], uipt[ko2], &kstat);
		  if (kstat < 0)
		    goto error;

		  sh6idcon (&rintdat, &uipt[ki2], &uipt[ko1], &kstat);
		  if (kstat < 0)
		    goto error;
		  /* Newi (ujk) Set in/out direction */
		  sh6setdir (uipt[ki2], uipt[ko1], &kstat);
		  if (kstat < 0)
		    goto error;

		}
	      *jstat = 0;
	    }
	  else
	    *jstat = 1;
	}

      else
	/* Prepare a marching strategy for connection.*/
	{
	  lstatus[0] = 1;
	  lstatus[1] = -1;
	  lstatus[2] = 0;
	  lstatus[3] = 2;

	  lnumb[0] = 0;
	  lnumb[1] = 0;
	  lnumb[2] = 0;
	  lnumb[3] = 0;

	  if ((lperm = new0array (kv, int)) == SISL_NULL)
	    goto err101;
	  if ((lpermdir = new0array (kv, int)) == SISL_NULL)
	    goto err101;
	  if ((spar = newarray (4 * kv, double)) == SISL_NULL)
	    goto err101;

	  for (kn = 0, kj = 0; kn < 4 && kj < kv; kn++)
	    {
	      /* We pass through the points four times, one for each
                 possible statusvalue.*/
	      /* In this way they will be sorted : 111..-1-1-1..000..222.. */

	      for (ki = 0; ki < kn1 && kj < kv; ki++)
		{

		  if (ldir[ki] == lstatus[kn])
		    {
		      lnumb[kn] += 1;
		      spar[4 * kj] = uipt[ki]->epar[0];
		      spar[4 * kj + 1] = uipt[ki]->epar[1];
		      spar[4 * kj + 2] = uipt[ki]->epar[2];
		      spar[4 * kj + 3] = uipt[ki]->epar[3];
		      lperm[kj] = ki;
		      lpermdir[kj] = ldir[ki];
		      kj++;
		    }
		}
	    }

	  /* If the first point is a parallel point, change status. */
	  if (lpermdir[0] == 0)
	    lpermdir[0] = 11;


	  /* UJK, aug. 92 */
	  if (kv == 3)
	    {
	       /* When three points, check if they are connected. */

	       sh6getlist (uipt[lperm[0]], uipt[lperm[1]],
			   &klist1, &klist2, &kstat);
	       if (kstat < 0) goto error;
	       sh6getlist (uipt[lperm[0]], uipt[lperm[2]],
			   &klist1, &klist2, &kstat1);
	       if (kstat1 < 0) goto error;
	       sh6getlist (uipt[lperm[1]], uipt[lperm[2]],
			   &klist1, &klist2, &kstat2);
	       if (kstat2 < 0) goto error;

	       if (kstat + kstat1 + kstat2 <= 1)
		 {
		    *jstat=0;
		    goto out;
		 }
	    }
	  if (kv == 3 &&  lpermdir[2] == 2)
	     {
		/* When three points, one singular, check if they are connected. */

		sh6getlist (uipt[lperm[0]], uipt[lperm[2]],
			    &klist1, &klist2, &kstat);
		sh6getlist (uipt[lperm[1]], uipt[lperm[2]],
			    &klist1, &klist2, &kstat1);

		if (kstat == 0 && kstat1 == 0)
		{
		   *jstat=0;
		   goto out;
		}
		else if (kstat1 == 0 &&
			 abs(lpermdir[0]) == 1 &&
			     lpermdir[1] == 0  &&
			     lpermdir[2] == 2)
		   {
		      /* Sing point is connected to parallell,
			 connect singular point to in(out) point. */
		      sh6tomain (uipt[lperm[0]], &kstat);
		      sh6tomain (uipt[lperm[2]], &kstat);

		      sh6idcon (&rintdat, &uipt[lperm[0]],
				&uipt[lperm[2]], &kstat);
		      if (kstat < 0)
			 goto error;

			 /* Set in/out direction */
			 if (lpermdir[0] == 1)
			    sh6setdir (uipt[lperm[0]],
				       uipt[lperm[2]], &kstat);
			 else
			 sh6setdir (uipt[lperm[2]],
				    uipt[lperm[0]], &kstat);

			 if (kstat < 0)
			    goto error;

			    /* Set status JUNCTION point. */
			    uipt[lperm[2]]->iinter = SI_SING;


		      *jstat=0;
		      goto out;
		   }
	     }


	  /* UJK 18.09.90 Connecting whenever we have two points and
             at least one of them is moving in or out. */
	  if (kv == 2 && (abs (lpermdir[0]) == 1))
/*	  if (kv==2 &&
	      ((lpermdir[0]*lpermdir[1] == -1) ||
	       (abs(lpermdir[0])==1 && lpermdir[1] == 0))) */
	    {
	      /* We have only two points and there has to be a curve,
                 connect. */
	      sh6tomain (uipt[lperm[0]], &kstat);
	      sh6tomain (uipt[lperm[1]], &kstat);

	      sh6idcon (&rintdat, &uipt[lperm[0]], &uipt[lperm[1]], &kstat);
	      if (kstat < 0)
		goto error;

	      /* Newi (ujk) Set in/out direction */
	      if (lpermdir[0] == 1)
		sh6setdir (uipt[lperm[0]],
			   uipt[lperm[1]], &kstat);
	      else
		sh6setdir (uipt[lperm[1]],
			   uipt[lperm[0]], &kstat);

	      if (kstat < 0)
		goto error;

	      /* ALA && UJK 19.09.90 Set status JUNCTION point. */
	      if (lpermdir[1] == 2)
		uipt[lperm[1]]->iinter = SI_SING;


	      *jstat = 0;
	      goto out;
	    }

	  else if (kv > 4 && lpermdir[kv - 1] != -1)
	     {
		/* More than four points, some with status /= +-1;
		   , it's not a simple case. */
		*jstat = 1;
		goto out;
	     }


	  else if (kv == 3 &&
		   lpermdir[0] == 1 &&
		   lpermdir[1] == -1 &&
		   lpermdir[2] == 0)
	    /* Three points, one in, one out, one parallel, connect in-out,
               delete parallel.. */
	    {

	      sh6tomain (uipt[lperm[0]], &kstat);
	      sh6tomain (uipt[lperm[1]], &kstat);
	      sh6idcon (&rintdat, &uipt[lperm[0]], &uipt[lperm[1]], &kstat);
	      if (kstat < 0)
		goto error;

	      /* Newi (ujk) Set in/out direction */
	      sh6setdir (uipt[lperm[0]], uipt[lperm[1]], &kstat);

	      /* UPDATE: (ujk) Avoid killing on edge !!! */
	      /* UJK, Aug.92, can't delete the parallel point
		 if it's in a corner. */

	      /* Corner first surface ? */
	      for (kn1 = kj = 0; kj < 4; kj++)
		if (edg[lperm[2]] & (1<<kj)) kn1++;

	      /* Corner second surface ? */
	      if (kn1 < 2)
		for (kn1 = 0, kj = 4; kj < 8; kj++)
		  if (edg[lperm[2]] & (1<<kj)) kn1++;

	      if (kn1 < 2)
		{
		   sh6idkpt (&rintdat, &uipt[lperm[2]], 1, &kstat);
		   if (kstat < 0)
		     goto error;
		}

	      *jstat = 0;
	      goto out;
	    }


	    else if (kv == 3 && (lpermdir[0] == 1 || lpermdir[0] == -1) &&
		     lpermdir[1] == 0 && lpermdir[2] == 0)
	       {
		  /* Set up parameters for search for singularity.  */

		  slimit[0] = *(ps1->et1 + ps1->ik1 - 1);
		  slimit[1] = *(ps1->et1 + ps1->in1);
		  slimit[2] = *(ps1->et2 + ps1->ik2 - 1);
		  slimit[3] = *(ps1->et2 + ps1->in2);
		  slimit[4] = *(ps2->et1 + ps2->ik1 - 1);
		  slimit[5] = *(ps2->et1 + ps2->in1);
		  slimit[6] = *(ps2->et2 + ps2->ik2 - 1);
		  slimit[7] = *(ps2->et2 + ps2->in2);
		  tref = MAX(MAX(slimit[1]-slimit[0],slimit[3]-slimit[2]),
			     MAX(slimit[5]-slimit[4],slimit[7]-slimit[6]));

		  for (kj=0; kj<4; kj++)
		  {
		     start[kj] = DZERO;
		     for (ki=1; ki<3; ki++) start[kj] += spar[ki*4+kj];
		     start[kj] /= (double)2;
		  }

		  /* Search for singularity.  */

		  shsing(ps1, ps2, slimit, start, spos, &kstat);
		  if (kstat < 0) goto error;

		  if (kstat == 1)
		  {
		     /* A singularity is found. Check if it is an
			intersection point. */

		     s1421(ps1, 0, spos, &klfs, &klft, sval1, snorm1, &kstat);
		     if (kstat < 0) goto error;

		     s1421(ps2, 0, spos+2, &klfs, &klft, sval2, snorm2, &kstat);
		     if (kstat < 0) goto error;

		     tdist = s6dist(sval1, sval2, 3);
		     if (tdist < aepsge)
		     {
			/* A singular intersection point is found. Check if
			   it is identical to any of the existing. */

			for (ki=0; ki<3; ki++)
			{
			   for (kj=0; kj<4; kj++)
			      if (DNEQUAL(spos[kj]+tref,spar[4*ki+kj]+tref)) break;
			   if (kj == 4) break;
			}

			if (ki > 0 && ki < 3)
			{
			   /* Connect the in/out point to the singular,
			      parallel point. */

			   sh6tomain(uipt[lperm[0]],&kstat);
			   sh6tomain(uipt[ki],&kstat);
			   sh6idcon(&rintdat,&uipt[lperm[0]],&uipt[ki],&kstat);
			   if (kstat < 0) goto error;

			   /* Newi Set in/out direction */
			   if (lpermdir[0] == 1)
			      sh6setdir (uipt[lperm[0]], uipt[ki], &kstat);
			   else
			      sh6setdir (uipt[ki], uipt[lperm[0]], &kstat);

			   goto out;
			}

			if (ki == 3)
			{
			   /* The singular intersection point do not exist already.
			      Create intersection point.  */

			   qsing = hp_newIntpt (4, spos, DZERO,
						SI_ORD, SI_UNDEF, SI_UNDEF,
						SI_UNDEF, SI_UNDEF,
						0, 0, nullp, nullp);
			   if (qsing == SISL_NULL) goto err101;

			  /*  Check if it lies on an edge.       */

			   for (kj=0; kj<4; kj++)
			      if (DEQUAL(spos[kj]+tref,slimit[2*kj]+tref) ||
				  DEQUAL(spos[kj]+tref,slimit[2*kj+1]+tref))
				 break;

			   if (kj < 4)
			   {
			      /* Singular intersection point at an edge.
				 Insert between the two parallel points.  */

			      sh6insert(&rintdat,uipt[lperm[1]],uipt[lperm[2]],
					&qsing, &kstat);
			      if (kstat < 0) goto error;

			      /* Connect to in/out point.  */

			      sh6tomain(uipt[lperm[0]],&kstat);
			      sh6tomain(qsing,&kstat);
			      sh6idcon(&rintdat,&uipt[lperm[0]],&qsing,&kstat);
			      if (kstat < 0) goto error;

			      /* Mark point as singular.  */

			      qsing -> iinter = SI_SING;

			      /* Newi Set in/out direction */
			      if (lpermdir[0] == 1)
				 sh6setdir (uipt[lperm[0]], qsing, &kstat);
			      else
				 sh6setdir (qsing, uipt[lperm[0]], &kstat);

			      goto out;
			   }
			   else
			   {
			      /* Singular point in the inner. Connect to
				 all edge points. First put into data structure. */

			      sh6idnpt (&rintdat, &qsing, 1, &kstat);
			      if (kstat < 0)
				 goto error;

			      sh6tomain(uipt[lperm[0]],&kstat);
			      sh6tomain(uipt[lperm[1]],&kstat);
			      sh6tomain(uipt[lperm[2]],&kstat);
			      sh6tomain(qsing,&kstat);
			      sh6idcon(&rintdat,&uipt[lperm[0]],&qsing,&kstat);
			      if (kstat < 0) goto error;
			      sh6idcon(&rintdat,&uipt[lperm[1]],&qsing,&kstat);
			      if (kstat < 0) goto error;
			      sh6idcon(&rintdat,&uipt[lperm[1]],&qsing,&kstat);
			      if (kstat < 0) goto error;

			      /* Mark point as singular.  */

			      qsing -> iinter = SI_SING;

			      /* Newi Set in/out direction */
			      if (lpermdir[0] == 1)
				 sh6setdir (uipt[lperm[0]], qsing, &kstat);
			      else
				 sh6setdir (qsing, uipt[lperm[0]], &kstat);

			      goto out;
			   }
			}
		     }
		  }
		  /* else
		  {
		      Check if only one of the parallel point is connected
			to only one point and the other is connected to two,
			connect to the first parallel point.

		     if ((sh6nmbmain(uipt[lperm[1]],&kstat) +
			  sh6nmbhelp(uipt[lperm[1]],&kstat) == 1) &&
			 (sh6nmbmain(uipt[lperm[2]],&kstat) +
			  sh6nmbhelp(uipt[lperm[2]],&kstat) == 2))
		     {
			sh6tomain(uipt[lperm[0]],&kstat);
			sh6tomain(uipt[lperm[1]],&kstat);
			sh6idcon(&rintdat,&uipt[lperm[0]],&uipt[lperm[1]],&kstat);
			if (kstat < 0) goto error;

			 Newi Set in/out direction
			if (lpermdir[0] == 1)
			   sh6setdir (uipt[lperm[0]], uipt[lperm[1]], &kstat);
			else
			   sh6setdir (uipt[lperm[1]], uipt[lperm[0]], &kstat);

			goto out;
	             }
		     else if ((sh6nmbmain(uipt[lperm[2]],&kstat) +
			       sh6nmbhelp(uipt[lperm[2]],&kstat) == 1) &&
			      (sh6nmbmain(uipt[lperm[1]],&kstat) +
			       sh6nmbhelp(uipt[lperm[1]],&kstat) == 2))
		     {
			sh6tomain(uipt[lperm[0]],&kstat);
			sh6tomain(uipt[lperm[2]],&kstat);
			sh6idcon(&rintdat,&uipt[lperm[0]],&uipt[lperm[2]],&kstat);
			if (kstat < 0) goto error;

			 Newi Set in/out direction
			if (lpermdir[0] == 1)
			   sh6setdir (uipt[lperm[0]], uipt[lperm[2]], &kstat);
			else
			   sh6setdir (uipt[lperm[2]], uipt[lperm[0]], &kstat);

			 Newi Set in/out direction
			if (lpermdir[0] == 1)
			   sh6setdir (uipt[lperm[0]], uipt[lperm[2]], &kstat);
			else
			   sh6setdir (uipt[lperm[2]], uipt[lperm[0]], &kstat);

			goto out;
	             }
		     else
		     {
			 Try to march from the non-parallel point.
			   Set parameter to enter marching.

			kmarch = 1;
		     }
		  }   */
	       }

	  if (kv == 3 &&  abs(lpermdir[2]) == 1)
	  {
	     /* There is three points, all of which pointing in or out.
		Prepare for a search for a singularity.  */

	     slimit[0] = *(ps1->et1 + ps1->ik1 - 1);
	     slimit[1] = *(ps1->et1 + ps1->in1);
	     slimit[2] = *(ps1->et2 + ps1->ik2 - 1);
	     slimit[3] = *(ps1->et2 + ps1->in2);
	     slimit[4] = *(ps2->et1 + ps2->ik1 - 1);
	     slimit[5] = *(ps2->et1 + ps2->in1);
	     slimit[6] = *(ps2->et2 + ps2->ik2 - 1);
	     slimit[7] = *(ps2->et2 + ps2->in2);
	     tref = MAX(MAX(slimit[1]-slimit[0],slimit[3]-slimit[2]),
			MAX(slimit[5]-slimit[4],slimit[7]-slimit[6]));

	     for (kj=0; kj<4; kj++)
	     {
		start[kj] = DZERO;
		for (ki=0; ki<3; ki++) start[kj] += spar[ki*4+kj];
		start[kj] /= (double)3;
	     }

	     /* Search for singularity.  */

	     shsing(ps1, ps2, slimit, start, spos, &kstat);
	     if (kstat < 0) goto error;

	     if (kstat == 1)
	     {
		/* A singularity is found. Check if it is an
		   intersection point. */

		s1421(ps1, 0, spos, &klfs, &klft, sval1, snorm1, &kstat);
		if (kstat < 0) goto error;

		s1421(ps2, 0, spos+2, &klfs, &klft, sval2, snorm2, &kstat);
		if (kstat < 0) goto error;

		tdist = s6dist(sval1, sval2, 3);
		if (tdist < aepsge)
		{
		   /* A singular intersection point is found. Check if
		      it is identical to any of the existing. */

		   for (ki=0; ki<3; ki++)
		   {
		      for (kj=0; kj<4; kj++)
			 if (DNEQUAL(spos[kj]+tref,spar[4*ki+kj]+tref)) break;
		      if (kj == 4) break;
		   }

		   if (ki < 3)
		   {
		      /* Something is wrong. Non of the existing
			 points are classified as singular. Continue
			 subdividing.   */
		   }
		   else
		   {
		      /* The singular intersection point do not exist already.
			 Create intersection point.  */

		      qsing = hp_newIntpt (4, spos, DZERO,
					   SI_ORD, SI_UNDEF, SI_UNDEF,
					   SI_UNDEF, SI_UNDEF,
					   0, 0, nullp, nullp);
		      if (qsing == SISL_NULL) goto err101;

		      /* Connect to all edge points. First put into
			 data structure. */

		      sh6idnpt (&rintdat, &qsing, 1, &kstat);
		      if (kstat < 0)
			 goto error;

		      sh6tomain(uipt[lperm[0]],&kstat);
		      sh6tomain(uipt[lperm[1]],&kstat);
		      sh6tomain(uipt[lperm[2]],&kstat);
		      sh6tomain(qsing,&kstat);
		      sh6idcon(&rintdat,&uipt[lperm[0]],&qsing,&kstat);
		      if (kstat < 0) goto error;
		      sh6idcon(&rintdat,&uipt[lperm[1]],&qsing,&kstat);
		      if (kstat < 0) goto error;
		      sh6idcon(&rintdat,&uipt[lperm[1]],&qsing,&kstat);
		      if (kstat < 0) goto error;

		      /* Mark point as singular.  */

		      qsing -> iinter = SI_SING;

		      /* Newi Set in/out direction */
		      if (lpermdir[0] == 1)
			 sh6setdir (uipt[lperm[0]], qsing, &kstat);
		      else
			 sh6setdir (qsing, uipt[lperm[0]], &kstat);

		      goto out;
		   }
		}
	     }
	  }

	  if (kv == 2 || lnumb[2] < 2 || kmarch == 1)
	  {
	      /* If we got two points, we always try marching as the second
                 strategy. */
	      /* If we got more than two points, we march only when
                 no more than one parallel point is given.*/

	      s9surmarch (ps1, ps2, spar, lpermdir, kv, &sparout,
			  &lpar, &kpoints, &kstat);
	      if (kstat < 0)
		goto error;

	      /* Branch on the given result from s9surmarch. */
	      if (kstat == 0)
		/* The points are not connected, continue subdividing. */
		*jstat = 1;

	      else if (kstat == 2)
		/* Only singulear points, set status ok. */
		*jstat = 0;

	      else
		{
		  if (kpoints > kv)
		    {
		      if ((uinewpt = newarray (kpoints - kv, SISLIntpt *))
			  == SISL_NULL)
			goto err101;
		      for (kj = kv, ki = 0; kj < kpoints; kj++, ki++)
			{
			  /* For each new point returned from s9surmarch,
                             we create an intersection point. */
			  SISLIntpt *qt;
			  double sintpar[4];

			  sintpar[0] = sparout[4 * kj];
			  sintpar[1] = sparout[4 * kj + 1];
			  sintpar[2] = sparout[4 * kj + 2];
			  sintpar[3] = sparout[4 * kj + 3];

			  uinewpt[ki] = qt = hp_newIntpt (4, sintpar, DZERO,
						 SI_ORD, SI_UNDEF, SI_UNDEF,
							  SI_UNDEF, SI_UNDEF,
							0, 0, nullp, nullp);


			  if (qt == SISL_NULL)
			    goto err101;

			  sh6idnpt (&rintdat, &qt, 1, &kstat);
			  if (kstat < 0)
			    goto error;

			}
		    }

		  /* Now we connect the points. */
		  *jstat = 0;

		  for (kj = 0; kj < kv; kj++)
		    {
		      ki = lpar[kj] - 1;

		      if (lpar[kj] > 0)
			{
			  if (lpar[kj] > kv)
			    {
			      /* SISLPoint is to be connected to a new point */

			      sh6tomain (uipt[lperm[kj]], &kstat);
			      sh6tomain (uinewpt[ki - kv], &kstat);

			      sh6idcon (&rintdat, &uipt[lperm[kj]],
					&uinewpt[ki - kv], &kstat);
			      if (kstat < 0)
				goto error;
			      /* Newi (ujk) Set in/out direction */
			      if (lpermdir[kj] == 1)
				sh6setdir (uipt[lperm[kj]],
					   uinewpt[ki - kv], &kstat);
			      else if (lpermdir[kj] == -1)
				sh6setdir (uinewpt[ki - kv],
					   uipt[lperm[kj]], &kstat);

			      if (kstat < 0)
				goto error;
			    }
			  else
			    {
			      /* SISLPoint is to be connected to an old point */
			      sh6tomain (uipt[lperm[kj]], &kstat);
			      sh6tomain (uipt[lperm[ki]], &kstat);
			      sh6idcon (&rintdat, &uipt[lperm[kj]],
					&uipt[lperm[ki]], &kstat);

			      /* Newi (ujk) Set in/out direction */
			      if (lpermdir[kj] == 1 || lpermdir[ki] == -1)
				sh6setdir (uipt[lperm[kj]],
					   uipt[lperm[ki]], &kstat);
			      else if (lpermdir[kj] == -1 ||
				       lpermdir[ki] == 1)
				sh6setdir (uipt[lperm[ki]],
					   uipt[lperm[kj]], &kstat);
			      if (kstat < 0)
				goto error;

			    }
			}
		    }
		}
	    }
	  else
	    /* No marching is to be done. */
	    *jstat = 1;
	}
    }


  goto out;



/* Error in sub rutines.      */

error:*jstat = kstat;
  s6err ("sh1762_s9edgsscon", *jstat, 0);
  goto out;

/* Error in memory allocation.      */

err101:*jstat = -101;
  s6err ("sh1762_s9edgsscon", *jstat, 0);
  goto out;

/* Error dimension.      */

err200:*jstat = -200;
  s6err ("sh1762_s9edgsscon", *jstat, 0);
  goto out;

out:
  if (uipt != SISL_NULL)
    freearray (uipt);
  if (edg != SISL_NULL)
    freearray (edg);
  if (ldir != SISL_NULL)
    freearray (ldir);
  if (sval1 != SISL_NULL)
    freearray (sval1);


  if (uinewpt != SISL_NULL)
    freearray (uinewpt);
  if (spar != SISL_NULL)
    freearray (spar);
  if (sparout != SISL_NULL)
    freearray (sparout);
  if (lpar != SISL_NULL)
    freearray (lpar);
  if (lperm != SISL_NULL)
    freearray (lperm);
  if (lpermdir != SISL_NULL)
    freearray (lpermdir);
}

#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9edgpscon (SISLEdge * pedge, double alevel, SISLSurf * ps,
		   int isimple, SISLIntdat * rintdat, double aepsge, int *jstat)
#else
static void
sh1762_s9edgpscon (pedge, alevel, ps, isimple, rintdat, aepsge, jstat)
     SISLEdge *pedge;
     double alevel;
     SISLSurf *ps;
     int isimple;
     SISLIntdat *rintdat;
     double aepsge;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To connect intersection points on edges on a surface
*              into curves. The problem must be of type point-surface
*              intersection in one dimension.
*
*
*
* INPUT      : pedge    - SISLEdge intersection.
*              alevel   - The point value in intersection.
*              ps       - Surface in intersection.
*              isimple  - Flag telling if we are to try marching in troublesom cases.
*              aepsge   - Geometry tolerance
*
* OUTPUT     : rintdat  - Intersection dates to be updated.
*              jstat    - status messages
*                           = 1     : Not a simple case.
*                           = 0     : OK.
*                           < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-07.
*
*********************************************************************
*/
{
  int kstat,kstat1,kstat2;
  int *ldir = SISL_NULL;		/* Local array containing one of the statusvalues for
			         * each point:
	                         *  0 - The intersect.curve is parallel to one
			         *      parameter direction.
	                         *  1 - The intersect.curve has direction into the
			         *      domain.
	                         * -1 - The intersect.curve has direction out of the
			         *      domain.
	                         *  2 - The point is singulear.
	                         * 10 - The intersect.curve touch one corner of the
			         *      domain.
			         * ---------------------------------------------------
			         */

  unsigned char *edg = SISL_NULL;
  double *sval = SISL_NULL;
  SISLPtedge *qpt = SISL_NULL;
  SISLIntpt **uipt = SISL_NULL;
  SISLIntpt **uinewpt = SISL_NULL;

  double *spar = SISL_NULL;		/* Local array with parameter values used as
				   input to s9conmarch. */
  int *lperm = SISL_NULL;		/* Local permutation array after sorting
				   input points to s9conmarch.*/
  int *lpermdir = SISL_NULL;		/* Local array with status values used as
				   input to s9conmarch. */
  int lstatus[4];		/* Local array containing the possible status
				   constants -1,1,0,2. */
  int lnumb[4];			/* Local array containing the number of points
				   with status lstatus. */

  double *sparout = SISL_NULL;	/* Local array with parameter values used as
				   output from s9conmarch.*/
  int *lpar = SISL_NULL;		/* Local array containing the connection information
				   from s9conmarch.*/
  int kpoints;			/* Local integer containing number of points returned
				   from s9conmarch.*/
  int klist1, klist2;		/* Indices for Intpoints. */

  if (ps->idim != 1)
    goto err200;

  *jstat = 1;

  if (pedge->ipoint > 1)
    {
      int kn1, kn, ki, kj, kv, kant, klfs, klft, kdir, kpar;
      /* UJK 18.09.90  multiplying; treating near singularities as singularities.*/
      double ttol = (double) 1000000.0 * REL_COMP_RES;
      double tolpar = (double) 0.00001;

      double *snorm;
      double tmax;
      double *nullp = SISL_NULL;
      kant = pedge->ipoint;

      /* Allocate array of pointers to the points. */

      if ((uipt = newarray (kant, SISLIntpt *)) == SISL_NULL)
	goto err101;
      if ((edg = new0array (kant, unsigned char)) == SISL_NULL)
	goto err101;
      if ((ldir = new0array (kant, int)) == SISL_NULL)
	goto err101;
      if ((sval = newarray (4 * kant, double)) == SISL_NULL)
	goto err101;
      snorm = sval + 3 * kant;


      /* Update the arrays. */

      for (kn1 = kj = 0; kj < pedge->iedge; kj++)
	for (qpt = pedge->prpt[kj]; qpt != SISL_NULL; qpt = qpt->pnext)
	  {
	    for (ki = 0; ki < kn1; ki++)
	      {
		if (qpt->ppt == uipt[ki])
		  break;
	      }
	    if (ki == kn1)
	      uipt[kn1++] = qpt->ppt;

	    edg[ki] |= 1 << kj;
	  }

      /* NEWI (ujk) some tuning ?! */
      if (kn1 == 2)
	tolpar = ttol;


      if (kn1 > 1)
	for (ki = 0; ki < kn1; ki++)
	  {
	    kn = 3 * ki;

	    klfs = klft = 0;
	    s1421 (ps, 1, uipt[ki]->epar, &klfs, &klft, sval + kn, snorm + ki, &kstat);
	    if (kstat < 0)
	      goto error;
	    else if (kstat > 0)
	      {
		ldir[ki] = 2;
		continue;
	      }

	    tmax = sqrt (sval[kn + 1] * sval[kn + 1] + sval[kn + 2] * sval[kn + 2]);
	    if (tmax < ttol)
	      {
		ldir[ki] = 2;
		continue;
	      }

	    for (kpar = 1, kj = 0; kj < 4; kj++)
	      if ((edg[ki] & 1 << kj) == 1 << kj)
		{
		  switch (kj)
		    {
		    case 0:
		      if (fabs (sval[kn + 1] / tmax) < tolpar)
			kdir = 0;
		      else
			kdir = (sval[kn + 1] > DZERO ? 1 : -1);
		      break;
		    case 1:
		      if (fabs (sval[kn + 2] / tmax) < tolpar)
			kdir = 0;
		      else
			kdir = (sval[kn + 2] > DZERO ? 1 : -1);
		      break;
		    case 2:
		      if (fabs (sval[kn + 1] / tmax) < tolpar)
			kdir = 0;
		      else
			kdir = (sval[kn + 1] > DZERO ? -1 : 1);
		      break;
		    case 3:
		      if (fabs (sval[kn + 2] / tmax) < tolpar)
			kdir = 0;
		      else
			kdir = (sval[kn + 2] > DZERO ? -1 : 1);
		    }

		  if (kdir == 0)
		    kpar = 0;
		  else if (ldir[ki] != kdir)
		    {
		      if (ldir[ki] == 0)
			ldir[ki] = kdir;
		      else
			{
			  ldir[ki] = 10;
			  break;
			}
		    }
		}
	    if (kpar == 0 && ldir[ki] != 10)
	      ldir[ki] = 0;
	  }
      /* End of for ki=0 ..... */


     /* When only two points, check if they are connected. */
      if (kn1 == 2)
	{
	  sh6getlist (uipt[0], uipt[1], &klist1, &klist2, &kstat);
	  if (kstat == 0)
	    kn1 = 0;
	}

      /* Count all the points that are candidates for connection. */
      for (kv = ki = 0; ki < kn1; ki++)
	if (ldir[ki] < 10)
	  kv++;


      if (kv == 1 && kn1 > 1 )
      {
	 int i,j;

	 if (kn1 > 2)
	 {
	    *jstat = 0;
	    for (i=0;i<kn1;i++)
	       if (ldir[i] == 1 || ldir[i] == -1)
		  *jstat = 1;
	    goto out;
	 }


	 if (ldir[0] == 1 || ldir[0] == -1)
	 {
	    i =  0;
	    j =  1;
	 }
	 else if (ldir[1] == 1 || ldir[1] == -1)
	 {
	    i =  1;
	    j =  0;
	 }
	 else						i = -1;

	 if (i >=0)
	 {
	    /* We have one point, and the direction is in/out.
	       The other point is touching a corner it is probably
	       an error and the point must be close to tangential.*/
	    if (ldir[i] == -1)
	    {
	       int k=j;
	       j  =  i;
	       i  =  k;
	    }

	    sh6tomain (uipt[i], &kstat);
	    sh6tomain (uipt[j], &kstat);

	    sh6idcon (&rintdat, &uipt[i], &uipt[j], &kstat);
	    if (kstat < 0)
	       goto error;
	    sh6setdir (uipt[i], uipt[j], &kstat);
	    if (kstat < 0)
	       goto error;
	    *jstat = 0;
	 }
      }
      else if (kv < 2)
	/* Less than two points, it's a simple case. */
	*jstat = 0;

      else if (kv > 8)
	/* More than eight points, it's not a simple case. */
	*jstat = 1;
      else
	/* We have two, three or four points, prepare a marching strategy
           for connection.*/
	{
	  lstatus[0] = 1;
	  lstatus[1] = -1;
	  lstatus[2] = 0;
	  lstatus[3] = 2;

	  lnumb[0] = 0;
	  lnumb[1] = 0;
	  lnumb[2] = 0;
	  lnumb[3] = 0;

	  if ((lperm = new0array (kv, int)) == SISL_NULL)
	    goto err101;
	  if ((lpermdir = new0array (kv, int)) == SISL_NULL)
	    goto err101;
	  if ((spar = newarray (2 * kv, double)) == SISL_NULL)
	    goto err101;

	  for (kn = 0, kj = 0; kn < 4 && kj < kv; kn++)
	    {
	      /* We pass through the points four times, one for each
                 possible statusvalue.*/
	      /* In this way they will be sorted : 111..-1-1-1..000..222.. */

	      for (ki = 0; ki < kn1 && kj < kv; ki++)
		{

		  if (ldir[ki] == lstatus[kn])
		    {
		      lnumb[kn] += 1;
		      spar[2 * kj] = uipt[ki]->epar[0];
		      spar[2 * kj + 1] = uipt[ki]->epar[1];
		      lperm[kj] = ki;
		      lpermdir[kj] = ldir[ki];
		      kj++;
		    }
		}
	    }

	  /* If the first point is a parallel point, change status. */
	  if (lpermdir[0] == 0)
	    lpermdir[0] = 11;

	  /* UJK, aug. 92 */
	  if (kv == 3)
	    {
	       /* When three points, check if they are connected. */

	       sh6getlist (uipt[lperm[0]], uipt[lperm[1]],
			   &klist1, &klist2, &kstat);
	       if (kstat < 0) goto error;
	       sh6getlist (uipt[lperm[0]], uipt[lperm[2]],
			   &klist1, &klist2, &kstat1);
	       if (kstat1 < 0) goto error;
	       sh6getlist (uipt[lperm[1]], uipt[lperm[2]],
			   &klist1, &klist2, &kstat2);
	       if (kstat2 < 0) goto error;

	       if (kstat + kstat1 + kstat2 <= 1)
		 {
		    *jstat=0;
		    goto out;
		 }
	    }

	  if (kv == 3 &&  lpermdir[2] == 2)
	     {
		/* When three points, one singular, check if they are connected. */

		sh6getlist (uipt[lperm[0]], uipt[lperm[2]],
			    &klist1, &klist2, &kstat);
		sh6getlist (uipt[lperm[1]], uipt[lperm[2]],
			    &klist1, &klist2, &kstat1);

		if (kstat == 0 && kstat1 == 0)
		{
		   *jstat=0;
		   goto out;
		}
		else if (kstat1 == 0 &&
			 abs(lpermdir[0]) == 1 &&
			     lpermdir[1] == 0  &&
			     lpermdir[2] == 2)
		   {
		      /* Sing point is connected to parallell,
			 connect singular point to in(out) point. */
		      sh6tomain (uipt[lperm[0]], &kstat);
		      sh6tomain (uipt[lperm[2]], &kstat);

		      sh6idcon (&rintdat, &uipt[lperm[0]],
				&uipt[lperm[2]], &kstat);
		      if (kstat < 0)
			 goto error;

			 /* Set in/out direction */
			 if (lpermdir[0] == 1)
			    sh6setdir (uipt[lperm[0]],
				       uipt[lperm[2]], &kstat);
			 else
			 sh6setdir (uipt[lperm[2]],
				    uipt[lperm[0]], &kstat);

			 if (kstat < 0)
			    goto error;

			    /* Set status JUNCTION point. */
			    uipt[lperm[2]]->iinter = SI_SING;


		      *jstat=0;
		      goto out;
		   }


	     }


	  /* UJK 18.09.90 Connecting whenever we have two points and
             at least one of them is moving in or out. */
	  if (kv == 2 && (abs (lpermdir[0]) == 1))
/*	  if (kv==2 &&
	      ((lpermdir[0]*lpermdir[1] == -1) ||
	       (abs(lpermdir[0])==1 && lpermdir[1] == 0))) */
	    {
	      /* We have only two points and there has to be a curve,
                 connect. */
	      sh6tomain (uipt[lperm[0]], &kstat);
	      sh6tomain (uipt[lperm[1]], &kstat);

	      sh6idcon (&rintdat, &uipt[lperm[0]], &uipt[lperm[1]], &kstat);
	      if (kstat < 0)
		goto error;

	      /* Newi (ujk) Set in/out direction */
	      if (lpermdir[0] == 1)
		sh6setdir (uipt[lperm[0]],
			   uipt[lperm[1]], &kstat);
	      else
		sh6setdir (uipt[lperm[1]],
			   uipt[lperm[0]], &kstat);

	      if (kstat < 0)
		goto error;

	      /* ALA && UJK 19.09.90 Set status JUNCTION point. */
	      if (lpermdir[1] == 2)
		uipt[lperm[1]]->iinter = SI_SING;

	      /* UJK,aug 93, Problems in sh1d_div, simple case may result
		 in connecting to edge pts that is singular in the original
		 problem, create one extra, internal pt if possible. */
	      if (DEQUAL(uipt[lperm[0]]->epar[0],uipt[lperm[1]]->epar[0]) ||
		  DEQUAL(uipt[lperm[0]]->epar[1],uipt[lperm[1]]->epar[1]))
	      {
		 SISLObject *obj=SISL_NULL;
		 SISLIntpt *pint=SISL_NULL;
		 double start;
		 double *result;
		 double coor[2];
		 if ((obj = newObject(SISLCURVE))== SISL_NULL) goto err101;
		 if (fabs(uipt[lperm[0]]->epar[0]-uipt[lperm[1]]->epar[0]) >
		     fabs(uipt[lperm[0]]->epar[1]-uipt[lperm[1]]->epar[1]))
		 {
		    coor[0] = (uipt[lperm[0]]->epar[0]+uipt[lperm[1]]->epar[0])/(double)2.0;
		    result = coor+1;
		    s1437(ps,coor[0],&(obj->c1),&kstat);
		    if (kstat < 0) goto error;
		 }
		 else
		 {
		    coor[1] = (uipt[lperm[0]]->epar[1]+uipt[lperm[1]]->epar[1])/(double)2.0;
		    result = coor;
		    s1436(ps,coor[1],&(obj->c1),&kstat);
		    if (kstat < 0) goto error;
		 }
		 start = s1792(obj->c1->et,obj->c1->ik,obj->c1->in);

		 sh6ptobj(&alevel, obj, aepsge, &start, result, &kstat);
		 if (kstat < 0) goto error;
		 if (kstat == 1)
		 {
		    pint = hp_newIntpt (2, coor, DZERO,
					SI_ORD, SI_UNDEF, SI_UNDEF,
					SI_UNDEF, SI_UNDEF,
					0, 0, SISL_NULL, SISL_NULL);


		    if (pint == SISL_NULL)
		       goto err101;

		    sh6insert (&rintdat, uipt[lperm[0]], uipt[lperm[1]],
			       &pint,&kstat);
		    if (kstat < 0)
		       goto error;


		 }
		 if (obj) freeObject(obj);

	      }
	      /* ______E N D aug.93________ */

	      *jstat = 0;
	    }

	  else if (kv > 4 && lpermdir[kv - 1] != -1)
	    /* More than four points, some with status /= +-1;
               , it's not a simple case. */
	    *jstat = 1;

	  else if (kv == 3 &&
		   lpermdir[0] == 1 &&
		   lpermdir[1] == -1 &&
		   lpermdir[2] == 0)
	    /* Three points, one in, one out, one parallel, connect in-out,
               delete parallel.. */
	    {

	      sh6tomain (uipt[lperm[0]], &kstat);
	      sh6tomain (uipt[lperm[1]], &kstat);

	      sh6idcon (&rintdat, &uipt[lperm[0]], &uipt[lperm[1]], &kstat);
	      if (kstat < 0)
		goto error;
	      sh6setdir (uipt[lperm[0]],
			 uipt[lperm[1]], &kstat);
	      if (kstat < 0)
		goto error;

	      /* UJK, Aug.92, can't delete the parallel point
		 if it's in the corner. */

	      for (kn1 = kj = 0; kj < pedge->iedge; kj++)
		if (edg[lperm[2]] & (1<<kj)) kn1++;

	       if (kn1 < 2)
		 {
		    sh6idkpt (&rintdat, &uipt[lperm[2]], 1, &kstat);
		    if (kstat < 0)
		      goto error;
		 }

	      *jstat = 0;
	    }
	  /* ALA UJK 31.10.90 Some special analysis when simple case
             and two nondirectional points. */
	  else if (kv == 2 && isimple > 0)
	    {
	      if (edg[lperm[0]] & edg[lperm[1]])
		{
		  /* The two points are on the same edge. */
		  sh6getlist (uipt[lperm[0]], uipt[lperm[1]], &klist1, &klist2, &kstat);
		  if (kstat == 0)
		    /* The points are connected, ok. */
		    *jstat = 0;
		  else
		    /* We must subdivide further. */
		    *jstat = 1;
		}
	      else
		{
		  /* The points are not on the same edge, connect. */

		  sh6tomain (uipt[lperm[0]], &kstat);
		  sh6tomain (uipt[lperm[1]], &kstat);

		  sh6idcon (&rintdat, &uipt[lperm[0]], &uipt[lperm[1]], &kstat);
		  if (kstat < 0)
		    goto error;
		  /* Set status JUNCTION point. */
		  if (lpermdir[0] == 2)
		    uipt[lperm[0]]->iinter = SI_SING;
		  if (lpermdir[1] == 2)
		    uipt[lperm[1]]->iinter = SI_SING;

		  *jstat = 0;
		}
	    }


	  /* end 31.10.90 */

	  /* UJK 18.09.90 Marching singular cases is VERY expensiv and leads
             to nothing.  (Number of points is here more than 2)*/
	  else if (lnumb[3] >= 2)
	    *jstat = 1;
	  /*end 18.09.90 */

	  else if (kv == 2 || (isimple > 0 && (lnumb[2] < 2)))
	    {
	      /* If we got two points, we always try marching as the second
                 strategy. */
	      /* If we got more than two points, we march only when the
                 input parameter isimple is set and the no more than one
                 parallel point is given.*/

	      s9conmarch (ps, alevel, spar, lpermdir, kv, &sparout,
			  &lpar, &kpoints, &kstat);
	      if (kstat < 0)
		goto error;

	      /* Branch on the given result from s9conmarch. */
	      if (kstat == 0)
		/* The points are not connected, continue subdividing. */
		*jstat = 1;

	      else if (kstat == 2)
		/* Only singulear points, set status ok. */
		*jstat = 0;

	      else
		{
		  if (kpoints > kv)
		    {
		      if ((uinewpt = newarray (kpoints - kv, SISLIntpt *))
			  == SISL_NULL)
			goto err101;
		      for (kj = kv, ki = 0; kj < kpoints; kj++, ki++)
			{
			  /* For each new point returned from s9conmarch,
                             we create an intersection point. */
			  SISLIntpt *qt;
			  double sintpar[2];

			  sintpar[0] = sparout[2 * kj];
			  sintpar[1] = sparout[2 * kj + 1];

			  uinewpt[ki] = qt = hp_newIntpt (2, sintpar, DZERO,
						 SI_ORD, SI_UNDEF, SI_UNDEF,
							  SI_UNDEF, SI_UNDEF,
							0, 0, nullp, nullp);


			  if (qt == SISL_NULL)
			    goto err101;

			  sh6idnpt (&rintdat, &qt, 1, &kstat);
			  if (kstat < 0)
			    goto error;

			}
		    }

		  /* Now we connect the points. */
		  *jstat = 0;

		  for (kj = 0; kj < kv; kj++)
		    {
		      ki = lpar[kj] - 1;

		      if (lpar[kj] > 0)
			{
			  if (lpar[kj] > kv)
			    {
			      /* SISLPoint is to be connected to a new point */
			      sh6tomain (uipt[lperm[kj]], &kstat);
			      sh6tomain (uinewpt[ki - kv], &kstat);

			      sh6idcon (&rintdat, &uipt[lperm[kj]],
					&uinewpt[ki - kv], &kstat);
			      if (kstat < 0)
				goto error;

			      /* Newi (ujk) Set in/out direction */
			      if (lpermdir[kj] == 1)
				sh6setdir (uipt[lperm[kj]],
					   uinewpt[ki - kv], &kstat);
			      else if (lpermdir[kj] == -1)
				sh6setdir (uinewpt[ki - kv],
					   uipt[lperm[kj]], &kstat);

			      if (kstat < 0)
				goto error;
			    }
			  else
			    {
			      /* SISLPoint is to be connected to an old point */
			      sh6tomain (uipt[lperm[kj]], &kstat);
			      sh6tomain (uipt[lperm[ki]], &kstat);
			      sh6idcon (&rintdat, &uipt[lperm[kj]],
					&uipt[lperm[ki]], &kstat);
			      if (kstat < 0)
				goto error;


			      /* Newi (ujk) Set in/out direction */
			      if (lpermdir[kj] == 1 || lpermdir[ki] == -1)
				sh6setdir (uipt[lperm[kj]],
					   uipt[lperm[ki]], &kstat);
			      else if (lpermdir[kj] == -1 ||
				       lpermdir[ki] == 1)
				sh6setdir (uipt[lperm[ki]],
					   uipt[lperm[kj]], &kstat);
			      if (kstat < 0)
				goto error;
			    }
			}
		    }
		}
	    }
	  else
	    /* The input parameter isimple is not set, no marching
               is to be done. */
	    *jstat = 1;
	}
    }


  goto out;

/* Error in sub rutines.      */

error:*jstat = kstat;
  s6err ("sh1762_s9edgpscon", *jstat, 0);
  goto out;

/* Error in memory allocation.      */

err101:*jstat = -101;
  s6err ("sh1762_s9edgpscon", *jstat, 0);
  goto out;

/* Error dimension.      */

err200:*jstat = -200;
  s6err ("sh1762_s9edgpscon", *jstat, 0);
  goto out;

out:
  if (uipt != SISL_NULL)
    freearray (uipt);
  if (edg != SISL_NULL)
    freearray (edg);
  if (ldir != SISL_NULL)
    freearray (ldir);
  if (sval != SISL_NULL)
    freearray (sval);

  if (uinewpt != SISL_NULL)
    freearray (uinewpt);
  if (spar != SISL_NULL)
    freearray (spar);
  if (sparout != SISL_NULL)
    freearray (sparout);
  if (lpar != SISL_NULL)
    freearray (lpar);
  if (lperm != SISL_NULL)
    freearray (lperm);
  if (lpermdir != SISL_NULL)
    freearray (lpermdir);
}

#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9simple (SISLObject * po1, SISLObject * po2, SISLEdge * vedge[],
		 int *jstat)
#else
static void
sh1762_s9simple (po1, po2, vedge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLEdge *vedge[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : A simle case test were we have edge intersections.
*              In this case we just test if every B-spline bases
*              from the orginal object is divided.
*
*
* INPUT      : po1      - First  object in intersection.
*              po2      - Second object in intersection.
*
*
* OUTPUT     : jstat    - status messages
*                           = 1     : A simple case.
*                           = 0     : No simel case.
*                           < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-06.
*
*********************************************************************
*/
{
  int kstat;
  int kk1, kk2;
  int kn1, kn2;
  int kdim, knum;
  int klin1 = 0, klin2 = 0;
  double *et1, *et2;
  double tstart1, tstart2, tend1, tend2;
  double tvolorg, tvol, tvolfac;
  SISLIntpt **up = SISL_NULL;

  kdim = po1->s1->idim;
  if (kdim != po2->s1->idim)
    goto err106;


  /* Init to no simpel case. */
  *jstat = 0;

  /* Count number of different edge-intersection points. */
  /* sh1762_s9edgpoint (vedge, &up, &knum, &kstat); */
  sh6edgpoint (vedge, &up, &knum, &kstat);
  if (kstat < 0)
    goto error;

  /* No point on edge, continue subdividing. */
  if (knum == 0)
    goto out;

  /* We test if one of the two objects is liniar. */
  if (po1->s1->pdir != SISL_NULL)
    if (po1->s1->pdir->igtpi == 0 && po1->s1->pdir->aang <= ANGULAR_TOLERANCE)
      klin1 = 1;

  if (po2->s1->pdir != SISL_NULL)
    if (po2->s1->pdir->igtpi == 0 && po2->s1->pdir->aang <= ANGULAR_TOLERANCE)
      klin2 = 1;
  /* UJK, This is a coincidence or simple case test ? */
  /* Both objects are linear, stop subdividing. */
  /*  if (klin1 == 1 && klin2 == 1)
     {
     *jstat = 1;
     goto out;
     } */

  /* Set different percent factor depending on number of edge points. */
  if (knum == 1)
    tvolfac = (double) 0.0001;
  else
    tvolfac = (double) 0.00390625;

  /* Get attributes from surface no 1. */
  tstart1 = po1->s1->et1[po1->s1->ik1 - 1];
  tend1 = po1->s1->et1[po1->s1->in1];
  tstart2 = po1->s1->et2[po1->s1->ik2 - 1];
  tend2 = po1->s1->et2[po1->s1->in2];

  et1 = po1->o1->s1->et1;
  et2 = po1->o1->s1->et2;
  kk1 = po1->o1->s1->ik1;
  kk2 = po1->o1->s1->ik2;
  kn1 = po1->o1->s1->in1;
  kn2 = po1->o1->s1->in2;

  tvolorg = (et1[kn1] - et1[kk1 - 1]) * (et2[kn2] - et2[kk2 - 1]);
  tvol = (tend1 - tstart1) * (tend2 - tstart2);

  /* Get attributes from surface no 1. */
  tstart1 = po2->s1->et1[po2->s1->ik1 - 1];
  tend1 = po2->s1->et1[po2->s1->in1];
  tstart2 = po2->s1->et2[po2->s1->ik2 - 1];
  tend2 = po2->s1->et2[po2->s1->in2];

  et1 = po2->o1->s1->et1;
  et2 = po2->o1->s1->et2;
  kk1 = po2->o1->s1->ik1;
  kk2 = po2->o1->s1->ik2;
  kn1 = po2->o1->s1->in1;
  kn2 = po2->o1->s1->in2;

  tvolorg = tvolorg * (et1[kn1] - et1[kk1 - 1]) * (et2[kn2] - et2[kk2 - 1]);
  tvol = tvol * (tend1 - tstart1) * (tend2 - tstart2);

  if (tvol <= tvolorg * tvolfac)
    *jstat = 1;

  goto out;

  /* Error in input. Dimensions of curves conflicting. */

err106:*jstat = -106;
  s6err ("sh1762_s9simple", *jstat, 0);
  goto out;

  /* Error in sub rutines.      */

error:*jstat = kstat;
  s6err ("sh1762_s9simple", *jstat, 0);
  goto out;

out:if (up != SISL_NULL)
    freearray (up);
}

#if 0
#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9reex (SISLObject * po1, SISLObject * po2, SISLEdge * vedge[],
	       double aepsge, SISLIntdat * pintdat, int *jstat)
#else
static void
sh1762_s9reex (po1, po2, vedge, aepsge, pintdat, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLEdge *vedge[];
     double aepsge;
     SISLIntdat *pintdat;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Aa extra last test whether the rest of the rutine
*              have done theire jobb. If no edge points is
*              connected to an internal or another edge point.
*              we just connecet.
*
*
* INPUT      : po1      - First  object in intersection.
*              po2      - Second object in intersection.
*              vedge[2] - SISLEdge intersections.
*              aepsge   - Geometry tolerance.
*
*
* INPUT/OUTPUT: pintdat - Intersections dates.
*
*
* OUTPUT     : jstat    - status messages
*                           = 0     : OK!.
*                           = 1     : Updating is performed.
*                           < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-07.
* REWISED BY : Vibeke Skytt, SI, 92-10. Set status if updating performed.
*
*********************************************************************
*/
{
  int kstat;
  SISLIntpt **up = SISL_NULL;

  *jstat = 0;

  if (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE)
    {
      int knum = 0;

      if (vedge[0]->ipoint + vedge[1]->ipoint > 1)
	{
	  /* The edges might be corrupt caused by killpoint,
             check existance of edgepoints. */

	  int kn1, kn, ki, kj;
	  SISLPtedge *qpt;

	  kn1 = pintdat->ipoint;

	  /* Loop for both objects*/
	  for (kn = 0; kn < 2; kn++)
	    if (vedge[kn] != SISL_NULL)
	      /* Loop for objects edges*/
	      for (kj = 0; kj < vedge[kn]->iedge; kj++)
		/* Loop for all points on edge*/
		for (qpt = vedge[kn]->prpt[kj]; qpt != SISL_NULL; qpt = qpt->pnext)
		  {
		    /* Loop for all points in intersection data*/
		    for (ki = 0; ki < kn1; ki++)
		      if (qpt->ppt == pintdat->vpoint[ki])
			break;

		    /* SISLPoint not found, quit! */
		    if (ki == kn1)
		      goto out;

		  }

	  /* sh1762_s9edgpoint (vedge, &up, &knum, &kstat); */
	  sh6edgpoint (vedge, &up, &knum, &kstat);
	  if (kstat < 0)
	    goto error;
	}

      if (knum > 1)
	{
	  int ki, kj;
	  int klist1, klist2;
	  double *spar;

	  /* We have to examine if any intersection points
             point to an internal or an other edge point. */

	  for (ki = 0; ki < knum; ki++)
	    if (sh6nmbmain (up[ki], &kstat) > 0)
	      {
		for (kj = 0; kj < knum; kj++)
		  {
		    sh6getlist (up[ki], up[kj], &klist1, &klist2, &kstat);
		    if (kstat == 0)
		      goto out;
		  }

		for (kj = 0; kj < up[ki]->no_of_curves; kj++)
		  {
		    spar = up[ki]->pnext[kj]->epar;

		    if (spar[0] > po1->s1->et1[po1->s1->ik1 - 1] &&
			spar[0] < po1->s1->et1[po1->s1->in1] &&
			spar[1] > po1->s1->et2[po1->s1->ik2 - 1] &&
			spar[1] < po1->s1->et2[po1->s1->in2] &&
			spar[2] > po2->s1->et1[po2->s1->ik1 - 1] &&
			spar[2] < po2->s1->et1[po2->s1->in1] &&
			spar[3] > po2->s1->et2[po2->s1->ik2 - 1] &&
			spar[3] < po2->s1->et2[po2->s1->in2])
		      goto out;
		  }

	      }


	  sh1762_s9edgsscon (vedge, po1->s1, po2->s1, pintdat, 0,
			     aepsge, &kstat);
	  if (kstat < 0)
	    goto error;

	  /* Test if any connections is done in sh1762_s9edgsscon. */

	  if (kstat == 0) *jstat = 1;
	}
    }
  else if ((po1->iobj == SISLPOINT &&
	    po2->iobj == SISLSURFACE &&
	    po1->p1->idim == 1) ||
	   (po1->iobj == SISLSURFACE &&
	    po2->iobj == SISLPOINT &&
	    po2->p1->idim == 1))
    {
      int knum = 0;

      if (vedge[0] != SISL_NULL)
	knum = vedge[0]->ipoint;
      if (vedge[1] != SISL_NULL)
	knum += vedge[1]->ipoint;

      if (knum > 1)
	{
	  /* The edges might be corrupt caused by killpoint,
             check existance of edgepoints. */

	  int kn1, kn, ki, kj;
	  SISLPtedge *qpt;

	  kn1 = pintdat->ipoint;

	  /* Loop for both objects*/
	  for (kn = 0; kn < 2; kn++)
	    if (vedge[kn] != SISL_NULL)
	      /* Loop for objects edges*/
	      for (kj = 0; kj < vedge[kn]->iedge; kj++)
		/* Loop for all points on edge*/
		for (qpt = vedge[kn]->prpt[kj]; qpt != SISL_NULL; qpt = qpt->pnext)
		  {
		    /* Loop for all points in intersection data*/
		    for (ki = 0; ki < kn1; ki++)
		      if (qpt->ppt == pintdat->vpoint[ki])
			break;

		    /* SISLPoint not found, quit! */
		    if (ki == kn1)
		      goto out;

		  }

	  /* sh1762_s9edgpoint (vedge, &up, &knum, &kstat); */
	  sh6edgpoint (vedge, &up, &knum, &kstat);
	  if (kstat < 0)
	    goto error;
	}

      if (knum > 1)
	{
	  int ki, kj;
	  int klist1, klist2;
	  double *spar;
	  SISLSurf *qs1;

	  /* We have to examine if any intersection points
             point to an internal or an other edge point. */

	  if (po1->iobj == SISLSURFACE)
	    qs1 = po1->s1;
	  else
	    qs1 = po2->s1;

	  for (ki = 0; ki < knum; ki++)
	    if (sh6nmbmain (up[ki], &kstat) > 0)
	      {
		for (kj = 0; kj < knum; kj++)
		  {
		    sh6getlist (up[ki], up[kj], &klist1, &klist2, &kstat);
		    if (kstat == 0)
		      goto out;
		  }

		for (kj = 0; kj < up[ki]->no_of_curves; kj++)
		  {
		    spar = up[ki]->pnext[kj]->epar;

		    /* ALA and UJK 19.09.90, To treat the problem of
                       junction points, we have introduced equality
                       in this test. */
		    if (spar[0] >= qs1->et1[qs1->ik1 - 1] &&
			spar[0] <= qs1->et1[qs1->in1] &&
			spar[1] >= qs1->et2[qs1->ik2 - 1] &&
			spar[1] <= qs1->et2[qs1->in2])
		      goto out;
		  }

	      }


	  sh1762_s9edgpscon (vedge[(po1->iobj == SISLSURFACE ? 0 : 1)],
		       (po1->iobj == SISLSURFACE ? po2 : po1)->p1->ecoef[0],
			     qs1, 1, pintdat, aepsge, &kstat);
	  if (kstat < 0)
	    goto error;

	  /* Test if any connections is done in sh1762_s9edgpscon. */

	  if (kstat == 0) *jstat = 1;
	}
    }

  goto out;

/* Error in subroutines.      */

error:*jstat = kstat;
  s6err ("sh1762_s9reex", *jstat, 0);
  goto out;

out:if (up != SISL_NULL)
    freearray (up);
}

#endif /* if 0 */
#if defined(SISLNEEDPROTOTYPES)
static void
sh1762_s9ptiter (SISLObject * po1, SISLObject * po2, double aepsge,
	       SISLIntdat ** pintdat, SISLEdge *vedge[], int *jstat)
#else
static void
sh1762_s9ptiter (po1, po2, aepsge, pintdat, vedge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double     aepsge;
     SISLIntdat **pintdat;
     SISLEdge   *vedge[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : In point-object intersection, there is no overlap, but
*              the objects might intersect anyway. Test if an intersection
*              is possible, and in that case, iterate to find it. The
*              dimension of the geometry space is greater than one.
*
*
* INPUT      : po1      - First  object in intersection.
*              po2      - Second object in intersection.
*              aepsge   - Geometry tolerance.
*              vedge    - Intersections found at the ends of the objects.
*
*
* INPUT/OUTPUT: pintdat - Intersections data.
*
*
* OUTPUT     : jstat    - status messages
*                           = 0     : OK!.
*                           = 1     : Updating is performed.
*                           < 0     : error
*
*
* METHOD     : First try to intercept by checking how close the point is
*              to the endpoints/corner points of the other object, and that
*              it is at the same side of the endpoints as the object. If this
*              is not possible, iteration is performed starting from the
*              closest endpoint. If the iteration produced an intersection
*              this is stored, otherwise no intersection is possible.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Vibeke Skytt, SI, 92-10.
*
*********************************************************************
*/
{
   int kstat = 0;
   int kturn;      /* Indicates if the order of the objects is changed.      */
   int kdim;         /* Dimension of geometry space.                         */
   int kcrv;         /* Index of curve in array of edge intersections.       */
   double tptdist1;  /* Distance between point and first endpoint of curve.  */
   double tptdist2;  /* Distance between point and second endpoint of curve. */
   double tcoefd1,tcoefd2;   /* Distance between vertices of curve.          */
   double tstart,tend; /* Endparameters of curve used in iteration.          */
   double tref;        /* Referance value in equality test.                  */
   double tpar,tres;   /* Start- and endparameter of intersection point.     */
   double *sc1;      /* Pointer to coefficient of curve.                     */
   double *sc2;      /* Pointer to coefficient of curve.                     */
   double sdiff1[3]; /* Vector between point and endpoint of object.         */
   double sdiff2[3]; /* Vector between endpoint of object and closest inner
			vertex.                                              */
   double *nullp = SISL_NULL;
   SISLPoint *qpt;   /* Pointer to the point in the intersection.            */
   SISLObject *qobj2; /* Pointer to the other object in the intersection.    */
   SISLCurve *qcrv;  /* Pointer to the curve in the intersection.            */
   SISLIntpt *qt = SISL_NULL; /* Pointer to intersection point.                   */

   /* Find the order of the objects. */

   if (po1->iobj == SISLPOINT)
   {
      qpt = po1->p1;
      qobj2 = po2;
      kturn = 0;
   }
   else if (po2->iobj == SISLPOINT)
   {
      qpt = po2->p1;
      qobj2 = po1;
      kturn = 1;
   }
   else
      goto err122;

   /* Find dimension of geometry space. */

   kdim = qpt->idim;

   if (qobj2->iobj == SISLCURVE)
   {
      /* Curve object intersection. Compute distances between point and
	 endpoints of curve, and between endpoints of curve and closest
	 vertex.        */

      kcrv = 1 - kturn;
      qcrv = qobj2->c1;
      sc1 = qcrv->ecoef;
      sc2 = qcrv->ecoef + kdim*(qcrv->in - 1);

      if (qcrv->in == 2)
      {
	 /* No intersection is possible.  */

	 *jstat = 0;
	 goto out;
      }

      tptdist1 = s6dist(qpt->ecoef,sc1,kdim);
      tptdist2 = s6dist(qpt->ecoef,sc2,kdim);
      tcoefd1 = s6dist(sc1,sc1+kdim,kdim);
      tcoefd2 = s6dist(sc2,sc2-kdim,kdim);

      if (tptdist1 > (double)1.5*tcoefd1 && tptdist2 > (double)1.5*tcoefd2)
      {
	 /* The point is not close to an endpoint of the curve.
	    No intersection. NB! It may be necessary to make this
	    test less strict when we have gained some experience
	    with this routine. */

	 *jstat = 0;
	 goto out;
      }

      if (tptdist1 < tptdist2 && tptdist1 <= (double)1.5*tcoefd1)
      {
	 s6diff(qpt->ecoef,sc1,kdim,sdiff1);
	 s6diff(sc1+kdim,sc1,kdim,sdiff2);
      }
      else
      {
	 s6diff(qpt->ecoef,sc2,kdim,sdiff1);
	 s6diff(sc2-kdim,sc2,kdim,sdiff2);
      }

      /* Check if the point lies on the same side of the closest endpoint
	 of the curve as the curve itself.                                */

      if (s6scpr(sdiff1,sdiff2,kdim) < DZERO)
      {
	 /* No intersection.  */

	 *jstat = 0;
	 goto out;
      }

      /* Iterate to find the intersection point. */

      tstart = qcrv->et[qcrv->ik - 1];
      tend = qcrv->et[qcrv->in];
      tref = tend - tstart;

      if (tptdist1 < tptdist2 && tptdist1 <= (double)1.5*tcoefd1)
      {
	 tpar = tstart;

	 /* Check if there exists an intersection in the start of the
	    curve.  */

	 if (vedge[kcrv]->prpt[0] != SISL_NULL)
	 {
	    /* No iteration is to be performed. */

	    *jstat = 0;
	    goto out;
	 }
      }
      else
      {
	 tpar = tend;

	 /* Check if there exists an intersection in the end of the
	    curve.  */

	 if (vedge[kcrv]->prpt[1] != SISL_NULL)
	 {
	    /* No iteration is to be performed. */

	    *jstat = 0;
	    goto out;
	 }
      }

      s1771 (qpt, qcrv, aepsge,
	     tstart, tend, tpar, &tres, &kstat);
      if (kstat < 0)
	 goto error;

      if (kstat == 1)
	 /*Intersection point found. Control edges. */
	 if (DEQUAL (tres+tref, tstart+tref) || DEQUAL (tres+tref, tend+tref))
	    kstat = 0;

      if (kstat == 1)	/* Intersection point found. */
      {
	 *jstat = 1;	/* Mark intersection found.  */

	 /* Making intersection point. */
	 qt = hp_newIntpt (SISLCURVE, &tres, DZERO, SI_ORD,
			   SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
			   0, 0, nullp, nullp);

	 if (qt == SISL_NULL)
	    goto err101;

	 /* Uppdating pintdat. */
	 sh6idnpt (pintdat, &qt, 1, &kstat);
	 if (kstat < 0)
	    goto error;
      }
   }
   else
   {
      /* The other object is a surface. For the time being, set no
	 intersection. This part of the routine will be implemented
	 later. */

      *jstat = 0;
      goto out;
   }

   goto out;

   /* Error in space allocation.  */

   err101 : *jstat = -101;
   goto out;

   /* None of the objects is a point. */

   err122 : *jstat = -122;
   goto out;

   /* Error in lower level routine.  */

   error : *jstat = kstat;
   goto out;


   out:
      return;
}
