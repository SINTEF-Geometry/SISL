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
 * $Id: s1162.c,v 1.3 2001-03-19 15:58:41 afr Exp $
 *
 */


#define S1162

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void s1162_s9mic(SISLObject *,SISLObject *,SISLIntdat **,
			SISLEdge *[],int *);
static void s1162_s9num(SISLObject *,int *,int *);
static void s1162_s9edge(SISLObject *[],SISLObject *[],int,int,
			 SISLIntdat *,SISLEdge *[],int *);
static void s1162_s9con(SISLObject *,double *,double,SISLIntdat **,
			SISLEdge *[],int *,int *,int *);
static void s1162_s9update(SISLObject *,double *,double,SISLIntdat **,
			   SISLEdge *[2],int *);
static void s1162_s9div(SISLObject *,double *,double,int,int,int,
			SISLObject *[],SISLIntdat **,SISLEdge *[2],int,int *);
#else
static void s1162_s9mic();
static void s1162_s9num();
static void s1162_s9edge();
static void s1162_s9con();
static void s1162_s9update();
static void s1162_s9div();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
s1162(SISLObject *po1,double *cmax,double aepsge,
	   SISLIntdat **pintdat,SISLEdge *vedge[2],
	   int ilevel,int inum,int *jstat)
#else
void s1162(po1,cmax,aepsge,pintdat,vedge,ilevel,inum,jstat)
     SISLObject *po1;
     double *cmax;
     double aepsge;
     SISLIntdat **pintdat;
     SISLEdge   *vedge[2];
     int    ilevel;
     int    inum;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : SISLObject - object maximum. Treat the inner of the
*              object.
*
*
*
* INPUT      : po1       - Pointer to  object
*              aepsge    - Geometry resolution.
*              vedge     - Pointers to structure of edge-maximums.
*                          vedge[1] must be SISL_NULL
*              ilevel    - Debt in recursion with inumb(>2) max. on the edges(if bezier case).
*              inum      - Number of max. on the edges.
*
* INPUT/OUTPUT : pintdat - Pointer to maximum data.
*                cmax    - Level value.
*
*
* OUTPUT     : jstat     - status messages
*                                   = 2 : maximum found equal to level value
*                                   = 1 : maximum found over to level value
*                                   = 0 : no maximum
*                                   < 0 : error
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6err      - Treating error situation.
*              s1119      - Simple Case test for maximums.
*              s1190      - Box test.
*              s1791      - Test if possible to subdivide
*              s1792      - Computing midpoint of parameter intervall.
*              s1770      - Curve/curve iteration.
*              s1770p     - Point/curve iteration.
*              s1771      - Curve/surface iteration.
*              s1771p     - Point/surface iteration.
*              s1231      - Subdivide curve.
*              s1711      - Subdivide surface.
*              s1161      - Object/object maximum.
*              s1435      - Pick an edge curve from a surface.
*              s1438      - Pick an end point from a curve.
*              s6idnpt    -
*              s6idkpt    -
*              s6idcpt    -
*              s6idcon    -
*              s6idput    -
*              s6idedg    -
*              s6idint    -
*              newPoint   -
*              newObject  -
*              newIntpt   -
*              newEdge    -
*              newIntdat  -
*              freeObject - Free space occupied by a given object
*              freeCurve  -
*              freeEdge   -
*              freeIntdat -
*
* WRITTEN BY : Ulf J. Krystad, SI, 89-05.
*
*********************************************************************
*/
{
  int klevel;             /* Local - Debt in recursion with.    */
  int knumedge;           /* Local - Number of max. on the edges*/
  int kpos  = 0;          /* Position of error.                 */
  int kstat = 0;          /* Local error status.                */
  int ksimple = 0;        /* Local simple case status.          */
  int kdiv  = 0;          /* Parameter direction of subdivsion. */
  int knum;               /* Number of edges in subproblems.    */
  int ki;                 /* Counter.                           */
  int kind1,kind2;        /* Index two knots with multiplicity. */
  SISLObject *uob1[4];        /* Pointers to subdivided object.     */
  SISLObject *qdum;           /* Pointer to dummy object.           */
  SISLEdge **uedge=SISL_NULL;      /* Pointer to array (to be allocated)
				    of edges to use in subproblems.    */
  SISLIntpt *up[2];
  SISLPtedge *qpt0,*qpt1;

  /*Init*/
  knumedge   = inum;
  klevel     = ilevel;

  for (ki=0;ki<4;ki++)  uob1[ki] = SISL_NULL;
  if ((qdum = newObject(SISLPOINT)) == SISL_NULL) goto err101;

  /* Initiate no maximum.*/
  *jstat = 0;

  /* Test if maximum is possible (perform box-test).  */
  s1190(po1,cmax,aepsge,&kstat);
  if (kstat < 0) goto error;

  /* We may have four different values on kstat.
     kstat = 1 : The SISLbox is beyond level value.
     kstat = 2 : The object is of constant value.
     kstat = 3 : The object is beyond one of its corners.
     kstat = 0 : No conclusion.*/

  if (kstat == 1);

  /* No max is possible */

  else if (kstat == 2)
    {
      /* The geometry is of constant value. Since it is not taken by
	 the SISLbox test and the edges already are treated in s1161,
	 we just connect the point on the edges. */


      if (vedge[0] != SISL_NULL && vedge[0]->iedge == 2)
	{
	  /* Only curves has to do connect */

	  qpt0=vedge[0]->prpt[0];
	  qpt1=vedge[0]->prpt[1];
	  if (qpt0 != SISL_NULL && qpt1 != SISL_NULL)
	    {

	      up[0] = qpt0->ppt;
	      up[1] = qpt1->ppt;
	      s6idcon(pintdat,&up[0],&up[1],&kstat);
	      if (kstat<0) goto error;
	    }
	}
    }

  else if (kstat == 3);


  /* Maximum for the object is a corner value, it has been found
     while treating the edges. */


  else
    {
      /* Simple Case test (more than one maximum possible?)  */
      if(po1->iobj ==SISLCURVE)

	s1119(po1->c1->ecoef,po1->c1->et,po1->c1->et,
	      po1->c1->ik,po1->c1->in,
	      1,1,&ksimple,&kind1,&kind2,&kstat);

      else
	s1119(po1->s1->ecoef,po1->s1->et1,po1->s1->et2,
	      po1->s1->ik1,po1->s1->in1,
	      po1->s1->ik2,po1->s1->in2,&ksimple,&kind1,&kind2,&kstat);
      if (kstat < 0) goto error;

      /* We may have three different values on ksimple.
	 ksimple = 0 : Not possible with interior max.
	 ksimple = 1 : Simpel case
	 ksimple = 2 : Not simpel case.*/

      if (ksimple == 0)
	*jstat = 0;

      else if (ksimple == 1)
	{
	  /* Simple Case, uppdate maximum list. */

	  s1162_s9update(po1,cmax,aepsge,pintdat,vedge,&kstat);
	  if (kstat < 0) goto error;
	  *jstat = kstat;

	}
      else
	{
	  /* Check for interval maximum.*/

	  s1162_s9con(po1,cmax,aepsge,pintdat,vedge,&klevel,&knumedge,&kstat);
	  if (kstat < 0) goto error;

	  /* We may have two different values on kstat.
	     kstat = 0 : No intervall maximum.
	     kstat = 1 : More than 2 maximum found on the edges.
	                 (bezier case only).
	     kstat = 2 : Intervall maximum found.
	     kstat = 3 : Simple case  */

	  if (kstat == 3)
	    /* Simple Case, uppdate maximum list. */
	    {

	      s1162_s9update(po1,cmax,aepsge,pintdat,vedge,&kstat);
	      if (kstat < 0) goto error;
	      *jstat = kstat;
	    }

	  else if (kstat == 2)

	    *jstat = kstat;     /*Uppdating maximum found. */

	  else
	    {
	      /* Find number of possible subdivision directions.
		 kdiv may have 4 difference values :
		 kdiv = 0 : Subdivision not possible.
		 kdiv = 1 : Subdivision in first parameter direction.
		 kdiv = 2 : Subdivision in second parameter direction.
		 kdiv = 3 : Subdivision in both parameter directions. */

	      s1162_s9num(po1,&kdiv,&kstat);
	      if (kstat < 0) goto error;


	      if(kdiv == 0)
		{
		  /* Microcase in parameter plane.*/

		  s1162_s9mic(po1,qdum,pintdat,vedge,&kstat);
		  if (kstat < 0) goto error;
		  else *jstat = kstat;
		}
	      else
		{
		  /* We do not have simpel case and it is possible to
		     subdivide. We therfor subdivide and uppdate the
		     edge maximum and then do a recurcive call
		     to treat the sub problems. Curves are subdivided
		     into two, surfaces into four. We can therfor get
		     up to four recursive calls.*/

		  /* Computing total number of subobjects in sub problems. */
		  knum = (kdiv<3 ? 2:4);

		  /***** Treating objects on sub problems. *****/

		  if (kdiv > 0) /* New objects for subdivision of po1. */
		    {
		      for (ki=0;ki<knum;ki++)
			{
			  if ((uob1[ki] = newObject(po1->iobj)) == SISL_NULL)
			    goto err101;

			  /*Initiate o1 pointer to point to top level object.*/

			  uob1[ki]->o1 = po1->o1;
			}

		      /* Subdivide the po1 object. */

		      s1162_s9div(po1,cmax,aepsge,kdiv,kind1,kind2,
			    uob1,pintdat,vedge,klevel,&kstat);
		      if (kstat < 0) goto error;
		      *jstat = max(*jstat,kstat);

		    }

		  /***** Treating edges on sub problems. *****/


		  /* Making array of pointers to edge object
		     to the sub problems. */
		  if ((uedge = new0array(2*knum,SISLEdge *)) == SISL_NULL)
		    goto err101;

		  /* Making new edge object to sub problems. */
		  for (ki=0; ki<2*knum; ki+=2)
		    {

		      if ((uedge[ki]   = newEdge(vedge[0]->iedge)) == SISL_NULL)
			goto err101;
		      /* No edge for the dummy point: */
		      uedge[ki+1] = SISL_NULL;

		    }


		  /***** Recursion. *****/
		  for (ki=0;ki<knum;ki+=1)
		    {

		      /* Uppdate edge maximum on sub problems. */
		      s1162_s9edge(uob1+ki, &qdum, 1, 1, *pintdat,
			     uedge+2*ki, &kstat);
		      if (kstat < 0) goto error;

		      s1162(uob1[ki],cmax,aepsge,pintdat,
			    uedge+2*ki,klevel,knumedge,&kstat);
		      if (kstat < 0) goto error;
		      else *jstat = max(*jstat, kstat);
		    }
		}
	    }
	}
    }

  /* Intersections in the inner found.  */

  goto out;

  /* Error in space allocation.         */
 err101: *jstat = -101;
  s6err("s1162",*jstat,kpos);
  goto out;

  /* Error in lower level routine.      */
  error : *jstat = kstat;
  s6err("s1162",*jstat,kpos);
  goto out;

  /* Free the space that is  allocated. */

 out:
  if (qdum != SISL_NULL) freeObject(qdum);

  for (ki=0;ki<4;ki++)
    if (uob1[ki] != SISL_NULL) freeObject(uob1[ki]);

  if (uedge != SISL_NULL)
    {
       /* 26.10.92 UJK/ BEOrd13969 */
       /* for (ki=0;ki<knum;ki++) */
       for (ki=0;ki<2*knum;ki++)
	  if (uedge[ki] != SISL_NULL) freeEdge(uedge[ki]);

      freearray(uedge);
    }
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1162_s9mic(SISLObject *po1,SISLObject *po2,SISLIntdat **rintdat,SISLEdge *vedge[],int *jstat)
#else
static void s1162_s9mic(po1,po2,rintdat,vedge,jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat **rintdat;
     SISLEdge   *vedge[];
     int    *jstat;
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
* INPUT      : vedge[2] - SISLEdge intersection objects to the two
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
  int kpos = 0;                 /* Position of error.                      */
  int kstat=0;                  /* Local error status.                     */
  int kpoint;                   /* Number of intpt on edges.               */
  double *spar = SISL_NULL;          /* Array to store parameter values.        */
  SISLIntpt **up = SISL_NULL;     /* Array of poiners to intersection point. */


  /* Initiate to now new intersection point. */


  *jstat = 0;


  /* Compute number of intersection points on edges. */

  if (vedge[0] == SISL_NULL )
    kpoint = 0;
  else
    kpoint = vedge[0]->ipoint;

  if (vedge[1] != SISL_NULL )
    kpoint += vedge[1]->ipoint;


  if (kpoint == 0 )
    {
      int kpar = 0;
      SISLIntpt *qt;


      /* There is not any intersection points on the edges.
	 We therfor make one new intersection point with parameter
	 values in senter of each object. */


      /* Number of parameter values of object 1. */

      if (po1->iobj == SISLCURVE) kpar = 1;
      else if (po1->iobj == SISLSURFACE) kpar = 2;


      /* Number of parameter values of object 2. */

      if (po2->iobj == SISLCURVE) kpar++;
      else if (po2->iobj == SISLSURFACE) kpar += 2;


      /* Allocate array to store midpoint parameter values. */

      if ((spar = newarray(kpar,double)) == SISL_NULL)
	goto err101;


      /* Compute midpoint parameter values. */

      if (po1->iobj == SISLCURVE)
	{
	  spar[0] = (po1->c1->et[po1->c1->ik - 1] +
		     po1->c1->et[po1->c1->in])*(double)0.5;
	  kpar = 1;
	}
      else if (po1->iobj == SISLSURFACE)
	{
	  spar[0] = (po1->s1->et1[po1->s1->ik1 - 1] +
		     po1->s1->et1[po1->s1->in1])*(double)0.5;
	  spar[1] = (po1->s1->et2[po1->s1->ik2 - 1] +
		     po1->s1->et2[po1->s1->in2])*(double)0.5;
	  kpar = 2;
	}

      if (po2->iobj == SISLCURVE)
	{
	  spar[kpar] = (po2->c1->et[po2->c1->ik - 1] +
			po2->c1->et[po2->c1->in])*(double)0.5;
	  kpar++;
	}
      else if (po2->iobj == SISLSURFACE)
	{
	  spar[kpar] = (po2->s1->et1[po2->s1->ik1 - 1] +
			po2->s1->et1[po2->s1->in1])*(double)0.5;
	  spar[kpar+1] = (po2->s1->et2[po2->s1->ik2 - 1] +
			  po2->s1->et2[po2->s1->in2])*(double)0.5;
	  kpar += 2;
	}

      *jstat = 1;         /* Mark intersection found. */


      /* Makeing intersection point. */

      qt = newIntpt(kpar,spar,DZERO);
      if (qt == SISL_NULL) goto err101;

      /* Uppdating pintdat. */

      s6idnpt(rintdat,&qt,1,&kstat);
      if (kstat < 0) goto error;
    }
  else if (kpoint > 1)
    {
      int kn,kn1,ki,kj;
      SISLPtedge *qpt;


      /* We have more than one intersection point on the edges,
	 we therfor conect these points to each other. */

      /* Allacate array of pointers to these points. */

      if ((up = newarray(kpoint,SISLIntpt *)) == SISL_NULL) goto err101;


      /* Uppdate the array. */

      for (kn=0,kn1=0; kn<2; kn++)
	if (vedge[kn] != SISL_NULL && vedge[kn]->ipoint > 0)
	  for(kj=0; kj<vedge[kn]->iedge; kj++)
	    for(qpt=vedge[kn]->prpt[kj]; qpt != SISL_NULL; qpt=qpt->pnext,kn1++)
	      up[kn1] = qpt->ppt;


      /* Connect the points to each other. */

      for (ki=1; ki<kpoint; ki++)
	{
	  s6idcon(rintdat,&up[ki-1],&up[ki],&kstat);
	  if (kstat<0) goto error;
	}
    }

  goto out;

  /* Error in space allocation.         */

 err101: *jstat = -101;
  s6err("s1162_s9mic",*jstat,kpos);
  goto out;

  /* Error in lower level routine.      */

  error : *jstat = kstat;
  s6err("s1162_s9mic",*jstat,kpos);
  goto out;

 out: if (spar != SISL_NULL) freearray(spar);
  if (up != SISL_NULL)   freearray(up);
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1162_s9num(SISLObject *po,int *jdiv,int *jstat)
#else
static void s1162_s9num(po,jdiv,jstat)
     SISLObject *po;
     int    *jdiv;
     int    *jstat;
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
*
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
* Revised BY : Christophe Rene Birkeland, SINTEF Oslo, May 93
*              *jstat = 0    initialization
*
*********************************************************************
*/
{
  *jstat = 0;
  if (po->iobj == SISLPOINT)                             *jdiv = 0;
  else if (po->iobj == SISLCURVE)
    {
      if(s1791(po->c1->et,po->c1->ik,po->c1->in))    *jdiv = 1;
      else                                           *jdiv = 0;
    }
  else if (po->iobj == SISLSURFACE)
    {
      if(s1791(po->s1->et1,po->s1->ik1,po->s1->in1)) *jdiv = 1;
      else                                           *jdiv = 0;

      if(s1791(po->s1->et2,po->s1->ik2,po->s1->in2)) *jdiv += 2;
    }
  else
    {

      /* Error. Kind of object does not exist.  */

      *jstat = -121;
      s6err("s1162_s9num",*jstat,0);
    }
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1162_s9edge(SISLObject *vob1[],SISLObject *vob2[],
		   int iobj1,int iobj2,SISLIntdat *pintdat,SISLEdge *wedge[],int *jstat)
#else
static void s1162_s9edge(vob1,vob2,iobj1,iobj2,pintdat,wedge,jstat)
     SISLObject  *vob1[];
     SISLObject  *vob2[];
     int     iobj1;
     int     iobj2;
     SISLIntdat  *pintdat;
     SISLEdge    *wedge[];
     int     *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE     : Uppdate edge structure. It may be up to four object
*               in vob1 and vob2, and wedge may therfor contain
*               seexteen elements to be uppdated.
*
*
*
* INPUT      : vob1[]  - Array of pointers to first objects.
*              vob2[]  - Array of pointers to second objects.
*              iobj1   - Number of elements in vob1.
*              iobj2   - Number of elements in vob2.
*              pintdat - intersection data.
*
*
*
* OUTPUT     : wedge[] - SISLEdge structures to uppdate.
*              jstat   - status messages
*                               = 0     : OK!
*                               < 0     : error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-05.
*
*********************************************************************
*/
{
  int kpos = 0;                 /* Position of error.       */
  int kstat=0;                  /* Local error status.      */
  int ki1,ki2,kj,kn;            /* Counters.                */
  int kedg;                     /* Number of edges.         */
  int kpar;                     /* Parameter number.        */
  double tpar;                  /* Parameter value at edge. */


  for (kn=0,ki1=0; ki1<iobj1; ki1++)
    for (ki2=0; ki2<iobj2; ki2++,kn+=2)
      {
        kedg = (vob1[ki1]->iobj == SISLPOINT ?0:(vob1[ki1]->iobj == SISLCURVE ?2:4));

	for (kj=0; kj<kedg; kj++)
	  {
	    if (vob1[ki1]->iobj == SISLCURVE)
	      {
		tpar = (kj == 0 ? vob1[ki1]->c1->et[vob1[ki1]->c1->ik-1] :
			vob1[ki1]->c1->et[vob1[ki1]->c1->in]);
		kpar = 1;
	      }
	    else if (kj == 0)
	      {
		tpar = vob1[ki1]->s1->et2[vob1[ki1]->s1->ik2-1];
		kpar = 2;
	      }
	    else if (kj == 1)
	      {
		tpar = vob1[ki1]->s1->et1[vob1[ki1]->s1->in1];
		kpar = 1;
	      }
	    else if (kj == 2)
	      {
		tpar = vob1[ki1]->s1->et2[vob1[ki1]->s1->in2];
		kpar = 2;
	      }
	    else
	      {
		tpar = vob1[ki1]->s1->et1[vob1[ki1]->s1->ik1-1];
		kpar = 1;
	      }


	    s6idedg(vob1[ki1],vob2[ki2],1,kpar,tpar,pintdat,
		    &(wedge[kn]->prpt[kj]),&(wedge[kn]->ipoint),&kstat);
	    if (kstat < 0) goto error;
	  }

        kedg = (vob2[ki2]->iobj == SISLPOINT ?0:(vob2[ki2]->iobj == SISLCURVE ?2:4));

	for (kj=0; kj<kedg; kj++)
	  {
	    if (vob2[ki2]->iobj == SISLCURVE)
	      {
		tpar = (kj == 0 ? vob2[ki2]->c1->et[vob2[ki2]->c1->ik-1] :
			vob2[ki2]->c1->et[vob2[ki2]->c1->in]);
		kpar = 1;
	      }
	    else if (kj == 0)
	      {
		tpar = vob2[ki2]->s1->et2[vob2[ki2]->s1->ik2-1];
		kpar = 2;
	      }
	    else if (kj == 1)
	      {
		tpar = vob2[ki2]->s1->et1[vob2[ki2]->s1->in1];
		kpar = 1;
	      }
	    else if (kj == 2)
	      {
		tpar = vob2[ki2]->s1->et2[vob2[ki2]->s1->in2];
		kpar = 2;
	      }
	    else
	      {
		tpar = vob2[ki2]->s1->et1[vob2[ki2]->s1->ik1-1];
		kpar = 1;
	      }


	    s6idedg(vob1[ki1],vob2[ki2],2,kpar,tpar,pintdat,
		    &(wedge[kn+1]->prpt[kj]),&(wedge[kn+1]->ipoint),&kstat);
	    if (kstat < 0) goto error;
	  }
      }

  *jstat = 0;

  goto out;

  /* Error in lower level routine.      */

  error : *jstat = kstat;
  s6err("s1162_s9edge",*jstat,kpos);
  goto out;

 out: ;
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1162_s9con(SISLObject *po1,double *cmax,double aepsge,SISLIntdat **pintdat,
		  SISLEdge *vedge[],int *jlevel,int *jnum,int *jstat)
#else
static void s1162_s9con(po1,cmax,aepsge,pintdat,vedge,jlevel,jnum,jstat)
     SISLObject *po1;
     double *cmax;
     double aepsge;
     SISLIntdat **pintdat;
     SISLEdge   *vedge[];
     int    *jlevel;
     int    *jnum;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if we are able to connect the maxpoints of the object
*              found on the edges.
*
*
*
* INPUT      : po1      - The first SISLObject to check.
*              cmax     - The max value found.
*              aepsge   - Geometrical resolution.
*              vedge[]  - SISLEdge intersection.
*
*
* INPUT/OUTPUT:pintdat  - Intersection data.
*              jlevel    - Debt in recursion with inumb(>2) max. on the edges(if bezier case).
*              jnum      - Number of max. on the edges.
*
*
* OUTPUT     :
*              jstat    - status messages
*                           = 3     : Simple Case.
*                           = 2     : Connection done.
*                           = 1     : Level value set to one or increased by one.
*                           = 0     : Level value set to zero.
*                           < 0     : error
*
*
* METHOD     : When  - object is a surface
*                    - the surf is a bezier patch
*                    - the patch has three (or more) max on the edges
*              If jlevel = 4 and number of max on the edges are equal to jnum,
*                      we connect the points treating one of them as singulear.
*              If jlevel in {1,2,3} and number of max
*              on the edges are equal to jnum,
*              we increase the jlevel by one.
*              If jlevel is 0 we set jlevel to one and
*              jnum to number of max on edge.
*              If jnum not equal to number of max on edge,
*              we set jlevel to one and
*              jnum to number of max on edge..
*
*              If the number of points on the edge is > 10, we do nothing.
* REFERENCES :
*
*
* WRITTEN BY : Ulf J Krystad, SI, 89-05.
*
*********************************************************************
*/
{
  SISLIntpt  *qintpt,*up[10];
  SISLPtedge *qpt;

  int kstat;      /* Local status.                */
/*guen  int kpos;   */   /* Local status counter.        */
/*guen changed into:*/
  int kpos = 0;   /* Local status counter.        */
  int kk1;        /* Local SURFACE attribute.     */
  int kk2;        /* Local SURFACE attribute.     */
  int kn1;        /* Local SURFACE attribute.     */
  int kn2;        /* Local SURFACE attribute.     */
  int ki,kj;      /* Local counter.               */
  int kfound;     /* Local flag in loop.          */
  int knum   = 0; /* Local number of max on edge. */
  int klevel = 0; /* Local level.                 */
  int kleft1 = 0; /* Local input parameter s1421  */
  int kleft2 = 0; /* Local input parameter s1421  */
  int kder   = 1; /* Local input parameter s1421  */

  double spar[2];  /* Parameter value              */
  double spar1[2]; /* Parameter value              */
  double smidle[2];/* middle parameter value       */
  double *sval=SISL_NULL;/*  Values from s1421.          */
  double *snorm=SISL_NULL;/* Values from s1421.         */


  kstat = 0;

  if (po1->iobj == SISLSURFACE)
    {
      if ((po1->s1->in1 == po1->s1->ik1) && (po1->s1->in2 == po1->s1->ik2))
	/* Bezier case for surface */
	{

	  /*-------------------------------------------------------*/
	  /* Count number of max on the edges. */
	  kk1 = po1->s1->ik1;
	  kk2 = po1->s1->ik2;
	  kn1 = po1->s1->in1;
	  kn2 = po1->s1->in2;

	  for (kj=0,knum=0;kj<vedge[0]->iedge;kj++)
	    {
	      qpt = vedge[0]->prpt[kj];
	      while(qpt != SISL_NULL)
		{
		  qintpt = qpt->ppt;
		  for (ki=0,kfound=0;ki<knum && kfound == 0;ki++)
		    if (qintpt == up[ki]) kfound = 1;

		  if (kfound == 0)
		    {
		      if (knum > 9) goto out;
		      up[knum]=qintpt;
		      knum++;
		    }

		  qpt = qpt->pnext;
		}

	    }

	  /*---------------------------------------------------------*/

	  if (knum > 0 )
	    /* Number of max on the edges more than 1. */
	    {
	      klevel = *jlevel;

	      if (klevel == 0 || knum !=*jnum)
		/* No continuation of suspected singulear point,
		   start a new one. */
		{
		  kstat = 1;
		  klevel = 1;
		}
	      else if (klevel < 2)
		/* Continuation of suspected singulear point. */
		{
		  kstat = 1;
		  klevel += 1;
		}
	      else if (knum < 2 )
		/* Simple Case */
		{
		  kstat = 3;
		  klevel += 1;
		}
	      else
		{

		  /*--------------------------------------------------*/
		  /* Connection case. */

		  /* Allocate local used memory */

		  sval = newarray(4,double);
		  if (sval == SISL_NULL) goto err101;
		  snorm = sval + 3;

		  for (kj=0;kj<knum-1;kj++)
		    {
		      spar[0] = up[kj]->epar[0];
		      spar[1] = up[kj]->epar[1];

		      for (ki=kj+1;ki<knum;ki++)
			{
			  /* First we linearize. */
			  spar1[0] = up[ki]->epar[0];
			  spar1[1] = up[ki]->epar[1];
			  smidle[0] = (spar[0] + spar1[0])/(double)2.0;
			  smidle[1] = (spar[1] + spar1[1])/(double)2.0;

			  /* Evaluate 0-1.st derivatives of surface */

			  s1421(po1->s1,kder,smidle,&kleft1,&kleft2,
				sval,snorm,&kstat);
			  if (kstat < 0) goto error;
			  if (fabs(sval[0]-*cmax) < aepsge)
			    {
			      /* Connect. */
			      s6idcon(pintdat,&up[kj],&up[ki],&kstat);
			      if (kstat<0) goto error;
			    }

			}
		    }


		  kstat = 2;
		  /*------------------------------------------*/


		}
	    }
	}
    }

  goto out;

  /* Error in allocation */
 err101: kstat = -101;
  s6err("s1162_s9con",kstat,kpos);
  goto out;

  /* Error in lower level routine.  */
  error : s6err("s1162_s9con",kstat,kpos);
  goto out;

 out:    if (sval != SISL_NULL) freearray(sval);
  *jlevel = klevel;
  *jnum   = knum;
  *jstat  = kstat;

}

#if defined(SISLNEEDPROTOTYPES)
static void
s1162_s9update(SISLObject *po1,double *cmax,double aepsge,
		     SISLIntdat **pintdat,SISLEdge *vedge[2],int *jstat)
#else
static void s1162_s9update(po1,cmax,aepsge,pintdat,vedge,jstat)
     SISLObject *po1;
     double *cmax;
     double aepsge;
     SISLIntdat **pintdat;
     SISLEdge   *vedge[2];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To find an inner maxima in an object when
*              it has no more then one.
*
* INPUT      : po1      - The object.
*              aepsge   - Geometry resolution.
*              vedge    - Pointer to edge maximum.
*
*
* INPUT/OUTPUT : cmax  - The level value.
*
* OUTPUT     : pintdat  - Maximum data.
*              jstat    - Status messages
*                          = 2     : New maximum found inside edges.
*                          = 1     : Maximum found inside object.
*                          = 0     : no maximum found.
*                          < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, 89-05.
*
*********************************************************************
*/
{
  int i, kj, ki;      /* Counters.                          */
  int kpos = 0;       /* Position of error.                 */
  int kstat= 0;       /* Local status                       */
  int kk1, kk2, kn1, kn2; /* Local number of knots and vertices.     */
  int kmax, kind1, kind2; /* Indexes of the maximum vertice.         */
  int kleft = 0;          /* Used in s1221 .                         */
  int kleft2 = 0;         /* Used in s1424 .                         */
  int kconn  = 0;         /* Connection flag.                        */
  int knum  = 0;          /* Number of max on the edge.              */
  int kfound = 0;         /* Flag.                                   */

  double tstart, tend;       /* Start, end values for curve parameter.    */
  double sstart[2], send[2]; /* Start, end values for surface parameter.  */
  double tpar;               /* The parameter vallue for
				subdivision of a curve.    */
  double spar[2];            /* The parameter vallue for subdivision
				of a surface.  */
  double tmax, tmin;         /* Local max and min value for the
				vertices of object. */
  double tval;               /* The value of the geometry at the found point.*/


  SISLObject *qop=SISL_NULL,*qcuo = SISL_NULL;/* Help pointers        */
  SISLIntdat *qintdat=SISL_NULL;         /* Local max data.      */
  SISLIntdat *qintdat1=SISL_NULL;        /* Local for double upgrading. */
  SISLIntpt  *qintpt,*up[3];
  SISLPtedge *qpt;

  /* Init */
  *jstat = 0;
  if (po1 == SISL_NULL || po1->iobj == SISLPOINT) goto out;
  if ((qop = newObject(SISLPOINT)) == SISL_NULL) goto err101;

  if (po1->iobj == SISLCURVE)
    {
      kk1   = po1->c1->ik;
      kn1   = po1->c1->in;
      kmax  = po1->c1->pbox->imax;
      tmax  = po1->c1->pbox->emax[0];
      tmin  = po1->c1->pbox->emin[0];

      tstart = po1->c1->et[kk1-1];
      tend   = po1->c1->et[kn1];

      /* Try to find an inner ekstremal point by iteration. */


      /* First get a good starting point for the iteration. */
      tpar = 0;
      for (i=kmax+1;i<kmax+kk1;i++)
	tpar += po1->c1->et[i];

      tpar /= kk1 - 1;

      s1252(po1->c1,aepsge,tpar,&tpar,&kstat);
      if (kstat < 0) goto error;

      /* Test if the found point is at start or end. */
      if(DEQUAL(tpar,tstart)  || DEQUAL(tpar,tend)) goto out;


      /* Evaluate curve at parameter value. */
      kleft = 0;
      s1221(po1->o1->c1,0,tpar,&kleft,&tval,&kstat);
      if (kstat < 0) goto error;

      /* Here we are ready to examine if we really found a new max point.*/
      if ((qop->p1 = newPoint(&tval,1,1)) == SISL_NULL) goto err101;

      s1161(qop,cmax,aepsge,&qintdat,&kstat);
      if (kstat < 0) goto error;

      if (kstat == 2)
	/* New maximum found, delete old ones */
	if (*pintdat != SISL_NULL)
	  {
	    freeIntdat(*pintdat);
	    *pintdat = SISL_NULL;
	  }

      if ( kstat )
	{
	  /* Maximum found, add it to the list */
	  *jstat = max(*jstat,kstat);         /* Mark maximum found. */

	  /* Set parameter parameter value of curve. */
	  s6idput(pintdat,qintdat,0,tpar,&kstat);
	  if (kstat < 0) goto error;
	}

    }
  else if (po1->iobj == SISLSURFACE)
    {


      kk1   = po1->s1->ik1;
      kn1   = po1->s1->in1;
      kk2   = po1->s1->ik2;
      kn2   = po1->s1->in2;
      kmax = po1->s1->pbox->imax;
      tmax = po1->s1->pbox->emax[0];
      tmin = po1->s1->pbox->emin[0];

      sstart[0] = po1->s1->et1[kk1-1];
      sstart[1] = po1->s1->et2[kk2-1];

      send[0]   = po1->s1->et1[kn1];
      send[1]   = po1->s1->et2[kn2];


      /* Get the two dimensional index of the greatest vertice. */
      kind2 = kmax/kk1;
      kind1 = kmax - kind2*kk1;


      /*-----------------------------------------------------------*/
      /* Count number of max on the edges. */

      for (kj=0,knum=0;kj<vedge[0]->iedge&&knum<3;kj++)
	{
	  qpt = vedge[0]->prpt[kj];
	  while(qpt != SISL_NULL && knum<3)
	    {
	      qintpt = qpt->ppt;
	      for (ki=0,kfound=0;ki<knum && kfound == 0;ki++)
		if (qintpt == up[ki]) kfound = 1;

	      if (kfound == 0)
		{
		  up[knum]=qintpt;
		  knum++;
		}

	      qpt = qpt->pnext;
	    }

	}

      /* Try if connection is possible.*/
      if (knum == 2)
	{
	  /* if on same edge, they are be connected before
	     (when in simple case.)*/
	  if ((DEQUAL(up[0]->epar[0],sstart[0]) &&
	       DEQUAL(up[1]->epar[0],sstart[0]))||
	      (DEQUAL(up[0]->epar[0],send[0])   &&
	       DEQUAL(up[1]->epar[0],send[0]))  ||
	      (DEQUAL(up[0]->epar[1],sstart[1]) &&
	       DEQUAL(up[1]->epar[1],sstart[1]))||
	      (DEQUAL(up[0]->epar[1],send[1])   &&
	       DEQUAL(up[1]->epar[1],send[1])))
	    kconn = 0;

	  else
	    {
	      /* Pick out two curves between the parameter
		 value on the edges. */
	      kconn = 0;
	      ki = 0;
	      if (fabs(up[0]->epar[0]-up[1]->epar[0]) <
		  fabs(up[0]->epar[1]-up[1]->epar[1]))
		ki =1;

	      tpar = (double)0.25*up[0]->epar[ki] +
		     (double)0.75*up[1]->epar[ki];
	      if ((qcuo = newObject(SISLCURVE)) == SISL_NULL) goto err101;
	      if (ki==0)
		s1437(po1->s1,tpar,&(qcuo->c1),&kstat);
	      else
		s1436(po1->s1,tpar,&(qcuo->c1),&kstat);
	      if (kstat < 0) goto error;

	      s1161(qcuo,cmax,aepsge,&qintdat,&kstat);
	      if (kstat < 0) goto error;

	      if (kstat == 1)
		{
		  freeCurve(qcuo->c1);
		  qcuo->c1 = SISL_NULL;

		  tpar = (double)0.75*up[0]->epar[ki] +
		         (double)0.25*up[1]->epar[ki];
		  if (ki==0)
		    s1437(po1->s1,tpar,&(qcuo->c1),&kstat);
		  else
		    s1436(po1->s1,tpar,&(qcuo->c1),&kstat);
		  if (kstat < 0) goto error;

		  s1161(qcuo,cmax,aepsge,&qintdat,&kstat);
		  if (kstat < 0) goto error;

		  if (kstat == 1)
		    {
		      /* Connect. */
		      kconn = 1;
		      s6idcon(pintdat,&up[0],&up[1],&kstat);
		      if (kstat<0) goto error;
		    }
		}
	    }
	}


      if (kconn == 0)
	{
	  /* No connection is done. */

	  /* Try to find an inner ekstremal point by iteration. */

	  /* First get a good starting point for the iteration. */
	  spar[0] = 0;
	  for (i=kind1+1;i<kind1+kk1;i++)
	    spar[0] += po1->s1->et1[i];

	  spar[0] /= kk1 - 1;

	  spar[1] = 0;
	  for (i=kind2+1;i<kind2+kk2;i++)
	    spar[1] += po1->s1->et2[i];

	  spar[1] /= kk2 - 1;


	  /* Create a point greater than the surface */
	  if ((qop->p1 = newPoint(&tmax,1,1)) == SISL_NULL) goto err101;

	  /* Iterate using aepsge=tmax-tmin to ensure covergence. */
	  s1173(qop->p1,po1->o1->s1,aepsge,sstart,send,spar,spar,&kstat);
	  if (kstat < 0) goto error;

	  /* Test if the found point is at start or end. */
	  if(DEQUAL(spar[0],sstart[0])  ||
	     DEQUAL(spar[0],send[0])    ||
	     DEQUAL(spar[1],sstart[1])  ||
	     DEQUAL(spar[1],send[1])) goto out;

	  /* Evaluate surface at parameter value. */
	  kleft  = 0;
	  kleft2 = 0;
	  s1424(po1->o1->s1,0,0,spar,&kleft,&kleft2,&tval,&kstat);
	  if (kstat < 0) goto error;

	  /* Here we are ready to examine if we really found a max point.*/
	  freePoint(qop->p1);
	  qop->p1 = SISL_NULL;
	  if ((qop->p1 = newPoint(&tval,1,1)) == SISL_NULL) goto err101;

	  s1161(qop,cmax,aepsge,&qintdat,&kstat);
	  if (kstat < 0) goto error;

	  if (kstat == 2)
	    /* New maximum found, delete old ones */
	    if (*pintdat != SISL_NULL)
	      {
		freeIntdat(*pintdat);
		*pintdat = SISL_NULL;
	      }

	  if ( kstat )
	    {
	      /* Maximum found, add them to the list */

	      *jstat = max(*jstat,kstat);  /* Mark maximum found. */

	      /* Special treatment for putting two
		 new parameters into pintdat from qintdat. */
	      s6idput(&qintdat1,qintdat,0,spar[0],&kstat);
	      if (kstat < 0) goto error;
	      s6idput(pintdat,qintdat1,1,spar[1],&kstat);
	      if (kstat < 0) goto error;

	    }
	}
    }


  goto out;

  /* -------------------ERROR SECTION----------------------------*/

  /* Error in space allocation.  */
 err101: *jstat = -101;
  s6err("s1162_s9update",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */
  error : *jstat = kstat;

 out:
  if (qcuo != SISL_NULL)
    {
      if (qcuo->c1 != SISL_NULL)
	{
	  freeCurve(qcuo->c1);
	  qcuo->c1 = SISL_NULL;
	}
      freeObject(qcuo);
      qcuo = SISL_NULL;
    }

  if (qop != SISL_NULL)
    {
      if (qop->p1 != SISL_NULL)
	{
	  freePoint(qop->p1);
	  qop->p1 = SISL_NULL;
	}
      freeObject(qop);
      qop = SISL_NULL;
    }
  if (qintdat != SISL_NULL)
    {
      freeIntdat(qintdat);
      qintdat = SISL_NULL;
    }
  if (qintdat1 != SISL_NULL)
    {
      freeIntdat(qintdat1);
      qintdat1 = SISL_NULL;
    }
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1162_s9div(SISLObject *po1,double *cmax,double aepsge,int idiv,int iind1,
		  int iind2,SISLObject *wob[],SISLIntdat **pintdat,SISLEdge *vedge[2],
		  int ilevel,int *jstat)
#else
static void s1162_s9div(po1,cmax,aepsge,idiv,iind1,iind2,
		  wob,pintdat,vedge,ilevel,jstat)
     SISLObject *po1;
     double *cmax;
     double aepsge;
     int    idiv;
     int    iind1;
     int    iind2;
     SISLObject *wob[];
     SISLIntdat **pintdat;
     SISLEdge   *vedge[2];
     int    ilevel;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Subdivide an object possible.
*
*
*
* INPUT      : po1      - The object to subdivide.
*              aepsge   - Geometry resolution.
*              idiv     - Subdivision direction.
*                          = 0     : No subdivision.
*		           = 1     : Subdivision in first parameter
*                                                          direction.
*		           = 2     : Subdivision in second parameter
*                                                          direction.
*		           = 3     : Subdivision in first and second
*                                                parameter direction.
*
*             iind1     - Index to first interior knot multiplicity in
*                         first parameter direction.
*             iind2     - Index to first interior knot multiplicity in
*                         second parameter direction.
*             ilevel    - If > 0  we subdivide in middlepoint.(SISLSurface only)
*              vedge    - Pointer to edge maximum.
*
* INPUT/OUTPUT : cmax  - The level value.
*
* OUTPUT     : pintdat  - Maximum data.
*              wob[]    - Pointers to the subdivided objects.
*              jstat    - Status messages
*                          = 2     : New maximum found on new edges.
*                          = 1     : Maximum found on new edges.
*                          = 0     : no maximum found.
*                          < 0     : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Ulf J. Krystad, SI, 89-05.
* REVISED BY : Paal Fugelli, SINTEF, Oslo, 94-07. Fixed memory leaks.
*
*********************************************************************
*/
{

  SISLPtedge *qpt;

  int ki, kj,i; /* Counters.                                           */
  int kpos = 0; /* Position of error.                                  */
  int kstat= 0; /* Local status                                        */
  int kk1, kk2, kn1, kn2; /* Local number of knots and vertices.       */
  int kmax, kind1, kind2; /* Indexes of the maximum vertice.           */
  double tstart, tend;    /* Start,end and middle values for curve parameter.*/
  double sstart[2], send[2],tmidle;/* Start, end values for surface parameter*/
  double tpar, tparold;  /* The parameter vallue for subdivision curve.*/
  double spar[2],sparold[2]; /* The parameter vallue for subdivision
				of a surface.  */
  double *tmax, *tmin;/* Local max and min value for the vertices of object.*/
  double sdiff[2];    /* The length of parameter intervall for surface.      */
  double smin[2];     /* The lower allowed limit in the prameter intervall
			 for subdividing a surface.      */
  double smax[2];    /* The upper allowed limit in the prameter intervall
			for subdividing a surface.      */
  SISLSurf *qs1=SISL_NULL;  /* Help pointers while subdividing        */
  SISLSurf *qs2=SISL_NULL;  /* Help pointers while subdividing        */
  SISLObject *qop = SISL_NULL;/* Help pointers while subdividing      */
  SISLObject *qoc = SISL_NULL;/* Help pointers while subdividing      */
  SISLIntdat *qintdat=SISL_NULL;/* Local max data for the new edges.  */



  /* Init */
  *jstat = 0;
  if ((qop = newObject(SISLPOINT)) == SISL_NULL) goto err101;

  if (po1 == SISL_NULL || po1->iobj == SISLPOINT)
    /* Nothing to do. */
    ;
  else if (po1->iobj == SISLCURVE)
    {
      kk1   = po1->c1->ik;
      kn1   = po1->c1->in;
      kmax  = po1->c1->pbox->imax;
      tmax  = po1->c1->pbox->emax;
      tmin  = po1->c1->pbox->emin;

      tstart = po1->c1->et[kk1-1];
      tend   = po1->c1->et[kn1];

      /* If we got problems with subdiv in max points, remove as comment: */
      /*   kmax = 0;  */


      /* ------------------Determination of sudiv parameter value-----------*/
      if (iind1 != 0)
	/* We subdivide in an interior knot with multiplicity. */
	tpar = po1->c1->et[iind1];

      else if (kmax == 0 || kmax == kn1-1)
	/*The greatest coeff is the first or last, divide in middlepoint. */
	tpar = s1792(po1->c1->et,kk1, kn1);


      else
	/* Try to find an inner subdivision (ekstremal) point by iteration. */
	{

	  /* First get a good starting point for the iteration. */
	  tpar = 0;
	  for (i=kmax+1;i<kmax+kk1;i++)
	    tpar += po1->c1->et[i];

	  tpar /= kk1 - 1;
	  tparold = tpar;

	  /*Iterate using Newton. */
	  s1252(po1->c1,aepsge,tpar,&tpar,&kstat);
	  if (kstat < 0) goto error;

	  /* Test if the found point is at start or end. */
	  if(DEQUAL(tpar,tstart)  || DEQUAL(tpar,tend))
	    /*Try Schoenbergs approximation to max vertice. */
	    {
	      tpar = tparold;

	      if(DEQUAL(tpar,tstart)  || DEQUAL(tpar,tend))
		/*Divide in middlepoint. */
		tpar = s1792(po1->c1->et,kk1,kn1);
	    }
	}
      /* ------------------Subdivision -------------------------------- */

      /* Subdivide the curve at the given parameter value. */
      s1231(po1->c1,tpar,&(wob[0]->c1),&(wob[1]->c1),&kstat);
      if (kstat < 0) goto error;


      /* Pick out end point from a curve. */
      s1438(wob[0]->c1,1,&(qop->p1),&tpar,&kstat);
      if (kstat < 0) goto error;


      /* Examin if the subdividing point is a max. */
      s1161(qop,cmax,aepsge,&qintdat,&kstat);
      if (kstat < 0) goto error;

      if (kstat == 2)
	/* New maximum found, delete old ones */
	if (*pintdat != SISL_NULL)
	  {

	    freeIntdat(*pintdat);
	    *pintdat = SISL_NULL;
	  }

      if (kstat)
	{
	  /* Maximum found, add them to the list */

	  *jstat = max(*jstat,kstat);         /* Mark maximum found. */

	  /* Put maximum found on edges into pintdat. */

	  /* Set parameter border values of object. */
	  s6idput(pintdat,qintdat,0,tpar,&kstat);
	  if (kstat < 0) goto error;

	  if (qintdat != SISL_NULL)
	    {
	      freeIntdat(qintdat);
	      qintdat = SISL_NULL;
	    }
	}
    }
  else if (po1->iobj == SISLSURFACE)
    {
      kk1   = po1->s1->ik1;
      kn1   = po1->s1->in1;
      kk2   = po1->s1->ik2;
      kn2   = po1->s1->in2;
      kmax = po1->s1->pbox->imax;
      tmax = po1->s1->pbox->emax;
      tmin = po1->s1->pbox->emin;

      sstart[0] = po1->s1->et1[kk1-1];
      sstart[1] = po1->s1->et2[kk2-1];

      send[0]   = po1->s1->et1[kn1];
      send[1]   = po1->s1->et2[kn2];

      sdiff[0] = send[0] - sstart[0];
      sdiff[1] = send[1] - sstart[1];
      smin[0]  = sstart[0] + (double)0.01*sdiff[0];
      smin[1]  = sstart[1] + (double)0.01*sdiff[1];
      smax[0]  = send[0] - (double)0.01*sdiff[0];
      smax[1]  = send[1] - (double)0.01*sdiff[0];

      kind2 = kmax/kn1;
      kind1 = kmax - kind2*kn1;


      /* If we got problems with subdiv in max points, remove as comment: */
      /*  kind1 = 0; */

      /* ------------------Determination of sudiv parameter value-------*/
      if (iind1 != 0 || iind2 != 0 || ilevel > 0)
	{
	  if (ilevel > 0)
	    /* We are forced to subdivide in middlepoint. */
	    {
	      spar[0] = s1792(po1->s1->et1,kk1, kn1);
	      spar[1] = s1792(po1->s1->et2,kk2, kn2);
	    }

	  else
	    /*We have knot multiplicity at least in one parameter direction.
	      Subdivide in interior knot multiplicity. If the other parameter
	      direction is without multiplicities, subdivide in middlepoint.*/
	    {
	      if (iind1 != 0)
		spar[0] = po1->s1->et1[iind1];
	      else
		spar[0] = s1792(po1->s1->et1,kk1, kn1);

	      if (iind2 != 0 )
		spar[1] = po1->s1->et2[iind2];
	      else
		spar[1] = s1792(po1->s1->et2,kk2, kn2);
	    }
	}


      else if (kind1 == 0 || kind1 == kn1-1 || kind2 == 0 || kind2 == kn2-1)
	{

	  /*The greatest coeff is on the edge.
	    Examin the edge for max and divide
	    in these parameter values. If more than one max,
	    use the one closest to the middlepoint*/

	  tmidle = s1792(po1->s1->et1,kk1, kn1);
	  spar[0] = sstart[0];

	  for (kj=0;kj<3;kj+=2)
	    {
	      qpt = vedge[0]->prpt[kj];
	      while (qpt != SISL_NULL)
		{
		  if (fabs(qpt->ppt->epar[0] - tmidle) <
		      fabs(spar[0] - tmidle))
		    spar[0] = qpt->ppt->epar[0];
		  qpt = qpt->pnext;
		}
	    }
	  if (DEQUAL(spar[0],sstart[0])  || DEQUAL(spar[0],send[0]))
	    spar[0] = tmidle;

	  tmidle = s1792(po1->s1->et2,kk2, kn2);
	  spar[1] = sstart[1];

	  for (kj=1;kj<4;kj+=2)
	    {
	      qpt = vedge[0]->prpt[kj];
	      while (qpt != SISL_NULL)
		{
		  if (fabs(qpt->ppt->epar[0] - tmidle) <
		      fabs(spar[1] - tmidle))
		    spar[1] = qpt->ppt->epar[1];
		  qpt = qpt->pnext;
		}
	    }
	  if (DEQUAL(spar[1],sstart[1])  || DEQUAL(spar[1],send[1]))
	    spar[1] = tmidle;
	}


      else
	/* Try to find an inner subdivision (ekstremal) point by iteration. */
	{

	  /* First get a good starting point for the iteration. */
	  spar[0] = 0;
	  for (i=kind1+1;i<kind1+kk1;i++)
	    spar[0] += po1->s1->et1[i];

	  spar[0] /= kk1 - 1;
	  sparold[0] = spar[0];

	  spar[1] = 0;
	  for (i=kind2+1;i<kind2+kk2;i++)
	    spar[1] += po1->s1->et2[i];

	  spar[1] /= kk2 - 1;
	  sparold[1] = spar[1];


	  /* Create a point greater than the surface */
	  if ((qop->p1 = newPoint(tmax,1,1)) == SISL_NULL) goto err101;

	  /* Iterate using Newton. */
	  s1173(qop->p1,po1->o1->s1,aepsge,sstart,send,spar,spar,&kstat);
	  freePoint(qop->p1);
	  qop->p1 = SISL_NULL;
	  if (kstat < 0) goto error;

	  /* Test if the found point is near one edge. */
	  if(spar[0] < smin[0] ||spar[0] > smax[0]
	     || spar[1] < smin[1] ||spar[1] > smax[1])
	    {
	      /*Try Schoenberg. */
	      spar[0] = sparold[0];
	      spar[1] = sparold[1];

	      if(spar[0] < smin[0] ||spar[0] > smax[0]
		 || spar[1] < smin[1] ||spar[1] > smax[1])
		{
		  /*Divide in middlepoint. */
		  spar[0] = s1792(po1->s1->et1,kk1,kn1);
		  spar[1] = s1792(po1->s1->et2,kk2,kn2);
		}
	    }


	}


      /* ------------------Subdivision ------------------------------*/
      /* Now we have found the parameters for subdivision, divide! */

      if ((qoc = newObject(SISLCURVE)) == SISL_NULL)
	goto err101;

      for (ki=0; ki<(idiv<3 ? 1:3); ki++)
	{

	  if (idiv == 1)
	    {
	      s1711(po1->s1,1,spar[0],&(wob[0]->s1),&(wob[1]->s1),&kstat);
	      if (kstat < 0) goto error;

	      /* Pick out edge curve from a surface. */

	      s1435(wob[0]->s1,1,&(qoc->c1),spar,&kstat);
	      if (kstat < 0) goto error;
	    }
	  else if (idiv == 2)
	    {
	      s1711(po1->s1,2,spar[1],&(wob[0]->s1),&(wob[1]->s1),&kstat);
	      if (kstat < 0) goto error;

	      /* Pick out edge curve from a surface. */

	      s1435(wob[0]->s1,2,&(qoc->c1),spar+1,&kstat);
	      if (kstat < 0) goto error;
	    }
	  else if (ki == 0)
	    {
	      s1711(po1->s1,1,spar[0],&qs1,&qs2,&kstat);
	      if (kstat < 0) goto error;

	      /* Pick out edge curve from a surface. */

	      s1435(qs1,1,&(qoc->c1),spar,&kstat);
	      if (kstat < 0) goto error;
	    }
	  else if (ki == 1)
	    {
	      s1711(qs1,2,spar[1],&(wob[0]->s1),&(wob[1]->s1),&kstat);
	      if (kstat < 0) goto error;

	      /* Pick out edge curve from a surface. */

	      s1435(wob[0]->s1,2,&(qoc->c1),spar+1,&kstat);
	      if (kstat < 0) goto error;
	    }
	  else   /* if (ki == 2) */
	    {
	      s1711(qs2,2,spar[1],&(wob[2]->s1),&(wob[3]->s1),&kstat);
	      if (kstat < 0) goto error;

	      /* Pick out edge curve from a surface. */

	      s1435(wob[2]->s1,2,&(qoc->c1),spar+1,&kstat);
	      if (kstat < 0) goto error;
	    }


	  /* Examine the new edge for max. */

	  s1161(qoc, cmax, aepsge, &qintdat, &kstat);
	  if (kstat < 0) goto error;

	  freeCurve(qoc->c1);
	  qoc->c1 = SISL_NULL;


	  if (kstat == 2)
	    /* New maximum found, delete old ones */
	    if (*pintdat != SISL_NULL)
	      {
		freeIntdat(*pintdat);
		*pintdat = SISL_NULL;
	      }

	  if (kstat)
	    {
	      /* Maximum found, add them to the list */

	      *jstat = max(kstat,*jstat);         /* Mark maximum found. */

	      /* Put maximum found on edges into pintdat. */

	      /* Test if we can pick the second subdivision parameter
		 from a max on subdiv curve.*/
	      if(ki==0 && qintdat->vpoint[0]->epar[0] > smin[1]
		 && qintdat->vpoint[0]->epar[0] < smax[1] )
		spar[1]=qintdat->vpoint[0]->epar[0];


	      /* Set parameter border values of object. */
	      s6idput(pintdat,qintdat,(ki==0 ? 0:1),spar[(ki==0 ? 0:1)],&kstat);
	      if (kstat < 0) goto error;

	      if (qintdat != SISL_NULL)
		{
		  freeIntdat(qintdat);
		  qintdat = SISL_NULL;
		}

	    }

	  /* End of for (ki=/..............) */
	}

    }
  goto out;

  /* -------------------ERROR SECTION------------------------------------*/

  /* Error in space allocation.  */
 err101: *jstat = -101;
  s6err("s1162_s9div",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */
  error : *jstat = kstat;
  s6err("s1162_s9div",*jstat,kpos);
  goto out;
  /* -------------------END OF ERROR SECTION----------------------------*/

 out:
  if (qop != SISL_NULL) freeObject(qop);
  if (qoc != SISL_NULL) freeObject(qoc);
  if (qs1 != SISL_NULL) freeSurf(qs1);  /* PFU 15/07-94 */
  if (qs2 != SISL_NULL) freeSurf(qs2);  /* PFU 15/07-94 */
}
