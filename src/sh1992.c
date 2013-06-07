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
 * $Id: sh1992.c,v 1.4 2007-08-06 13:09:13 vsk Exp $
 *
 */


#define SH1992

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void sh1992_s9mbox3(double [],int,int,double,double,double [],double []);
static void sh1992_s9mbox2(double [],int,int,double,double,double [],double []);
static void sh1992_s9mbox(double [],int,int,int,double,double,double[],double [],int *);
#else
static void sh1992_s9mbox3();
static void sh1992_s9mbox2();
static void sh1992_s9mbox();
#endif

#if defined(SISLNEEDPROTOTYPES)
void sh1992(SISLObject *po,int itype,double aepsge,int *jstat)
#else
void sh1992(po,itype,aepsge,jstat)
     SISLObject *po;
     int    itype;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make a SISLbox on the control-polygons given by
*              ecoef to the object. If dimension is 2 then one
*	       SISLbox rotated 45 degree is also made. If dimension
*	       is 3 then one SISLbox roteted 45 degree around each
*	       main axes, in all 4 boxes is made.
*
*
*
* INPUT      : po     - Object to treat.
*              itype  - Kind of box to make.
*                       = 0 : Do not expand box.
*                       = 1 : Make a totally expanded box.
*                       = 2 : Make a box expanded in the inner of the
*                             object, and reduced along the edges/endpoints.
*                       If itype>=10, it is interpreted as itype-10, exept
*                       for that no rotation of the box is performed.
*              aepsge - Geometry resolution.
*
* OUTPUT     : jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                                                     
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :  s6existbox - Check if the wanted box exist already.
*               s6newbox   - Make a box of a given type.
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
* MODIFIED BY : Vibeke Skytt, SI, 91-01. Tolerance dependant boxes.
* MODIFIED BY : Vibeke Skytt, SI, 92-10. Points are connected
*                                        non-expanded boxes. Exept for
*                                        1D, boxes are expanded with
*                                        eps/2 and eps.
*
*********************************************************************
*/                                     
{
   int kstat = 0;                       /* Status variable.        */
   int kdim;                       /* Dimension of geometry space. */
   int ktype = itype % 10;              /* Kind of box.            */
   int knum;                            /* Number of sides of box. */
   int k2;                              /* Other box type.         */
   int kbez = 0;                        /* Indicates if Bezier case. */
   double teps_inner;     /* Tolerance with which to expand in the inner. */
   double teps_edge;      /* Tolerance with which to expand at the edge.  */

   /* Set correct tolerances.  */
   
   teps_inner = (ktype == 0) ? DZERO : (double)0.5*aepsge;
   teps_edge = (ktype == 2) ? -teps_inner : teps_inner;
   
   if (po -> iobj == SISLPOINT)
   {
      if (po->p1->pbox == SISL_NULL)
	 if ((po->p1->pbox = newbox(po->p1->idim)) == SISL_NULL) goto err101;
      
      if (s6existbox(po->p1->pbox,ktype,aepsge) < 1)
      {
     	 kdim = po->p1->idim;
	 if (itype < 10 && kdim == 3) knum = 9;
	 else if (itype < 10 && kdim == 2) knum = 4;
	 else knum = kdim;
	   
	 /* The box do not exist already. For a point we always 
	    use non-expanded boxes.  */

	 /* Create the box.  */
	 
	 s6newbox(po->p1->pbox,knum,ktype,aepsge,&kstat);
	 if (kstat < 0) goto error;
	 
	 teps_inner = teps_edge = DZERO;
	 
	 k2 = (ktype == 0) ? 0 : ((ktype == 1) ? 2 : 1);
	 if (ktype > 0 && s6existbox(po->p1->pbox,k2,aepsge))
	    {
	       memcopy(po->p1->pbox->e2min[ktype],po->p1->pbox->e2min[k2],
		       (1+(kdim!=1))*knum,double);
	       memcopy(po->p1->pbox->e2max[ktype],po->p1->pbox->e2max[k2],
		       (1+(kdim!=1))*knum,double);
	    }
	    else
	    {
	       /* Make the requested box. */
	       
	       if (knum == 9) 
		  sh1992_s9mbox3(po->p1->ecoef,1,1,teps_inner,teps_edge,
			  po->p1->pbox->e2max[ktype],po->p1->pbox->e2min[ktype]);
	       else if (knum == 4)
		  sh1992_s9mbox2(po->p1->ecoef,1,1,teps_inner,teps_edge,
			  po->p1->pbox->e2max[ktype],po->p1->pbox->e2min[ktype]);
	       else
	       {
		  sh1992_s9mbox(po->p1->ecoef,1,1,kdim,teps_inner,teps_edge,
			 po->p1->pbox->e2max[ktype],po->p1->pbox->e2min[ktype],
			 &kstat);
		  if (kstat < 0) goto error;
	       }
	    }
      }
   }
   else if (po -> iobj == SISLCURVE)
   {
      if (po->c1->pbox == SISL_NULL)
	 if ((po->c1->pbox = newbox(po->c1->idim)) == SISL_NULL) goto err101;
      
      if (s6existbox(po->c1->pbox,ktype,aepsge) < 1)
      {
     	 kdim = po->c1->idim;
	 if (itype < 10 && kdim == 3) knum = 9;
	 else if (itype < 10 && kdim == 2) knum = 4;
	 else knum = kdim;
	 
	 /* The box do not exist already. In the Bezier case,
	    it is not necessary to expand in the inner of the curve.  */
	 
	 /* Create the box.  */
	 
	 s6newbox(po->c1->pbox,knum,ktype,aepsge,&kstat);
	 if (kstat < 0) goto error;
	 
	 /*if (po->c1->ik == po->c1->in) 
         {
            teps_inner = DZERO;
            kbez = 1;
	    }*/
	 
	 /* Make the requested box. First allocate scratch for
	    box arrays.  */
	 
	 if (knum == 9) 
	    sh1992_s9mbox3(po->c1->ecoef,po->c1->in,1,teps_inner,teps_edge,
		    po->c1->pbox->e2max[ktype],po->c1->pbox->e2min[ktype]);
	 else if (knum == 4)
	    sh1992_s9mbox2(po->c1->ecoef,po->c1->in,1,teps_inner,teps_edge,
		    po->c1->pbox->e2max[ktype],po->c1->pbox->e2min[ktype]);
	 else
	 {
	    sh1992_s9mbox(po->c1->ecoef,po->c1->in,1,kdim,teps_inner,
		   teps_edge,po->c1->pbox->e2max[ktype],
		   po->c1->pbox->e2min[ktype],&kstat);
	    if (kstat < 0) goto error;
         }
      }
   } 
   else if (po -> iobj == SISLSURFACE)
   {
      if (po->s1->pbox == SISL_NULL)
	 if ((po->s1->pbox = newbox(po->s1->idim)) == SISL_NULL) goto err101;
      
      if (s6existbox(po->s1->pbox,ktype,aepsge) < 1)
      {
     	 kdim = po->s1->idim;
	 if (itype < 10 && kdim == 3) knum = 9;
	 else if (itype < 10 && kdim == 2) knum = 4;
	 else knum = kdim;
	 
	 /* The box do not exist already. In the Bezier case, it
	    is not necessary to expand in the inner of the surface.  */
	 
	 /* Create the box.  */
	 
	 s6newbox(po->s1->pbox,knum,ktype,aepsge,&kstat);
	 if (kstat < 0) goto error;
	 
	 /*if (po->s1->ik1 == po->s1->in1 && po->s1->ik2 == po->s1->in2) 
         {
	    teps_inner = DZERO;
            kbez = 1;
	    }*/
	 
	 /* Make the requested box. First allocate scratch for
	    box arrays.  */
	 
	 if (knum == 9) 
	    sh1992_s9mbox3(po->s1->ecoef,po->s1->in1,po->s1->in2,teps_inner,
		    teps_edge,po->s1->pbox->e2max[ktype],
		    po->s1->pbox->e2min[ktype]);
	 else if (knum == 4)
	    sh1992_s9mbox2(po->s1->ecoef,po->s1->in1,po->s1->in2,teps_inner,
		    teps_edge,po->s1->pbox->e2max[ktype],
		    po->s1->pbox->e2min[ktype]);
	 else
	 {
	    sh1992_s9mbox(po->s1->ecoef,po->s1->in1,po->s1->in2,kdim,
		   teps_inner,teps_edge,po->s1->pbox->e2max[ktype],
		   po->s1->pbox->e2min[ktype],&kstat);
	    if (kstat < 0) goto error;
	 }
      }
   }  
  
  *jstat = kbez;
  goto out;

  /* Error in space allocation.  */
  
  err101 : *jstat = -101;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  goto out;
     
 out:
    return;
}

#if defined(SISLNEEDPROTOTYPES)
void sh1992cu(SISLCurve *pc,int itype,double aepsge,int *jstat)
#else
void sh1992cu(pc,itype,aepsge,jstat)
     SISLCurve *pc;
     int   itype;
     double aepsge;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make a SISLbox on the control-polygons given by
*              ecoef to the object. If dimension is 2 then one
*	       SISLbox rotated 45 degree is also made. If dimension
*	       is 3 then one SISLbox roteted 45 degree around each
*	       main axes, in all 4 boxes is made.
*
*
*
* INPUT      : pc     - SISLObject to treat.
*              itype  - Kind of box to make.
*                       = 0 : Do not expand box.
*                       = 1 : Make a totally expanded box.
*                       = 2 : Make a box expanded in the inner of the
*                             object, and reduced along the edges/endpoints.
*                       If itype>=10, it is interpreted as itype-10, exept
*                       for that no rotation of the box is performed.
*              aepsge - Geometry resolution.
*
*
* OUTPUT     : jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                                                     
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :  s6existbox - Check if the wanted box exist already.
*               s6newbox   - Make a box of a given type.
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
* MODIFIED BY : Vibeke Skytt, SI, 91-01. Tolerance dependant boxes.
*
*********************************************************************
*/                                     
{
   int kstat = 0;                       /* Status variable.        */
   int kdim = pc->idim;                 /* Dimension of geometry space. */
   int ktype = itype % 10;              /* Kind of box.            */
   int knum;                            /* Number of sides of box. */
   int kbez = 0;                        /* Indicates if Bezier case. */
   double teps_inner;     /* Tolerance with which to expand in the inner. */
   double teps_edge;      /* Tolerance with which to expand at the edge.  */

   /* Set number of box sides.  */
   
   if (itype < 10 && kdim == 3) knum = 9;
   else if (itype < 10 && kdim == 2) knum = 4;
   else knum = kdim;
   
   /* Set correct tolerances.  */
   
   teps_inner = (ktype == 0) ? DZERO : (double)0.5*aepsge;
   teps_edge = (ktype == 2) ? -teps_inner : teps_inner;
   
   if (pc->pbox == SISL_NULL)
      if ((pc->pbox = newbox(pc->idim)) == SISL_NULL) goto err101;
   
   if (s6existbox(pc->pbox,ktype,aepsge) < 1)
   {
      /* The box do not exist already. In the Bezier case,
	 it is not necessary to expand in the inner of the curve.  */
      
      /* Create the box.  */
      
      s6newbox(pc->pbox,knum,ktype,aepsge,&kstat);
      if (kstat < 0) goto error;
		     
      if (pc->ik == pc->in) 
      {
          teps_inner = DZERO;
          kbez = 1;
      }
      
      /* Make the requested box. First allocate scratch for
	 box arrays.  */
      
      if (knum == 9) 
	 sh1992_s9mbox3(pc->ecoef,pc->in,1,teps_inner,teps_edge,
		 pc->pbox->e2max[ktype],pc->pbox->e2min[ktype]);
      else if (knum == 4)
	 sh1992_s9mbox2(pc->ecoef,pc->in,1,teps_inner,teps_edge,
		 pc->pbox->e2max[ktype],pc->pbox->e2min[ktype]);
      else
      {
	 sh1992_s9mbox(pc->ecoef,pc->in,1,kdim,teps_inner,teps_edge,
		pc->pbox->e2max[ktype],pc->pbox->e2min[ktype],&kstat);
	 if (kstat < 0) goto error;
       }
   }
  
  *jstat = kbez;
  goto out;
  
  /* Error in space allocation.  */
  
  err101 : *jstat = -101;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  goto out;
  
 out:
    return;
}

#if defined(SISLNEEDPROTOTYPES)
void sh1992su(SISLSurf *ps,int itype,double aepsge,int *jstat)
#else
void sh1992su(ps,itype,aepsge,jstat)
     SISLSurf *ps;
     int  itype;
     double aepsge;
     int  *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make a SISLbox on the control-polygons given by
*              ecoef to the object. If dimension is 2 then one
*	       SISLbox rotated 45 degree is also made. If dimension
*	       is 3 then one SISLbox roteted 45 degree around each
*	       main axes, in all 4 boxes is made.
*
*
*
* INPUT      : ps     - SISLObject to treat.
*              itype  - Kind of box to make.
*                       = 0 : Do not expand box.
*                       = 1 : Make a totally expanded box.
*                       = 2 : Make a box expanded in the inner of the
*                             object, and reduced along the edges/endpoints.
*                       If itype>=10, it is interpreted as itype-10, exept
*                       for that no rotation of the box is performed.
*              aepsge - Geometry resolution.
*
*
* OUTPUT     : jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                                                     
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :  s6existbox - Check if the wanted box exist already.
*               s6newbox   - Make a box of a given type.
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
* MODIFIED BY : Vibeke Skytt, SI, 91-01. Tolerance dependant boxes.
*
*********************************************************************
*/                                     
{
   int kstat = 0;                       /* Status variable.        */
   int kdim = ps->idim;                 /* Dimension of geometry space. */
   int ktype = itype % 10;              /* Kind of box.            */
   int knum;                            /* Number of sides of box. */
   int kbez = 0;                        /* Indicates if Bezier.    */
   double teps_inner;     /* Tolerance with which to expand in the inner. */
   double teps_edge;      /* Tolerance with which to expand at the edge.  */

   /* Set correct tolerances.  */
   
   teps_inner = (ktype == 0) ? DZERO : (double)0.5*aepsge;
   teps_edge = (ktype == 2) ? -teps_inner : teps_inner;
   
   /* Set number of box sides.  */
   
   if (itype < 10 && kdim == 3) knum = 9;
   else if (itype < 10 && kdim == 2) knum = 4;
   else knum = kdim;
   
   if (ps->pbox == SISL_NULL)
      if ((ps->pbox = newbox(ps->idim)) == SISL_NULL) goto err101;
   
   if (s6existbox(ps->pbox,ktype,aepsge) < 1)
   {
      /* The box do not exist already. In the Bezier case, it
	 is not necessary to expand in the inner of the surface.  */
      
      /* Create the box.  */
      
      s6newbox(ps->pbox,knum,ktype,aepsge,&kstat);
      if (kstat < 0) goto error;
      
      if (ps->ik1 == ps->in1 && ps->ik2 == ps->in2) 
      {
	 teps_inner = DZERO;
         kbez = 1;
      }
      
      /* Make the requested box. First allocate scratch for
	 box arrays.  */
      
      if (knum == 9) 
	 sh1992_s9mbox3(ps->ecoef,ps->in1,ps->in2,teps_inner,teps_edge,
		 ps->pbox->e2max[ktype],ps->pbox->e2min[ktype]);
      else if (knum == 4)
	 sh1992_s9mbox2(ps->ecoef,ps->in1,ps->in2,teps_inner,teps_edge,
		 ps->pbox->e2max[ktype],ps->pbox->e2min[ktype]);
      else
      {
	 sh1992_s9mbox(ps->ecoef,ps->in1,ps->in2,kdim,
		teps_inner,teps_edge,ps->pbox->e2max[ktype],
		ps->pbox->e2min[ktype],&kstat);
	 if (kstat < 0) goto error;
      }
   }  
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
  err101 : *jstat = -101;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  goto out;
  
 out:
    return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
  sh1992_s9mbox3(double ecoef[],int icoef1,int icoef2,double aeps1,
		 double aeps2,double e2max[],double e2min[])
#else
static void sh1992_s9mbox3(ecoef,icoef1,icoef2,aeps1,aeps2,e2max,e2min)
     double ecoef[];
     int    icoef1;
     int    icoef2;
     double aeps1;
     double aeps2;
     double e2max[];
     double e2min[];
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make 4 boxes on the control-polygons given by
*              ecoef to the object.
*
*
* INPUT      : ecoef   - Control-polygon.
*              icoef1  - Number of vertices in control-polygon in
*                        first parameter direction.
*              icoef2  - Number of vertices in control-polygon in
*                        second parameter direction.
*              aeps1   - Size of expansion to use in the inner of
*                        the object.
*              aeps2   - Size of expansion to use at the edge of the object.
*
*                                                                     
*
* OUTPUT     : e2max   - Array to contain maximum values of the box.
*	       e2min   - Array to contain minimum values of the box.
*
*
* METHOD     : Make 4 boxes. One ordinary and tree boxes rotated
*	       45 degree around the main axes.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
* MODIFIED BY : Vibeke Skytt, SI, 91-01. 
*
*********************************************************************
*/                                     
{
  int ki,kj;             /* Counters.                                 */
  int kant = 9;          /* Number of box sides.                      */
  int kinset = 0;        /* Indicates if inner box is set.            */
  double teps1 = aeps1+aeps1; /* Double tolerance in the inner.       */
  double teps2;               /* Tolerance at edge.                   */
  double teps3;               /* Double tolerance at edge.            */
  double t1,t2,t3,t4;    /* To store elements of the rotation matrix. */
  double *tmin,*tmax;    /* Pointers used to traverse e2min and e2max.  */
  double sminin[9],smaxin[9];   /* Box boundaries in the inner.       */
  double sminedg[9],smaxedg[9]; /* Box boundaries at the edge.        */
  double *sc1,*sc2,*sc3; /* Pointers to coefficients.                 */
  
  /* Set tolerances at edge. If the tolerance is positive or dimension
     is 1D, the input tolerance is used, otherwise we must make sure 
     that the maximum distance from the total box at the edges to the
     reduced box is aeps2.                                             */
  
  if (aeps2 >= DZERO)
     teps2 = aeps2;
  else
     teps2 = (double)0.2767326953*aeps2; 
  teps3 = teps2 + teps2;

  /* Initiate box boundaries of inner box.  */
  
  for (ki=0; ki<kant; ki++)
  {
     sminin[ki] = HUGE;
     smaxin[ki] = -HUGE;
  }
  
  /* Fetch value of first vertex.  */
  
  sc1 = ecoef; sc2 = sc1+1;  sc3 = sc2 + 1;
  t1= ROTM * sc1[0];
  t2= ROTM * sc2[0];
  t3= ROTM * sc3[0];
  
  tmin = sminedg;
  tmax = smaxedg;
  *tmin = *tmax = *sc1;
  tmin++; tmax++;
  *tmin = *tmax = *sc2;
  tmin++; tmax++;
  *tmin = *tmax = *sc3;
  tmin++; tmax++;
  *tmin = *tmax = t2-t3;
  tmin++; tmax++;
  *tmin = *tmax = t2+t3;
  tmin++; tmax++;
  *tmin = *tmax = t1-t3;
  tmin++; tmax++;
  *tmin = *tmax = t1+t3;
  tmin++; tmax++;
  *tmin = *tmax = t1-t2;
  tmin++; tmax++;
  *tmin = *tmax = t1+t2;
  
  /* For each vertice at the edge check and corrigate the box.  */
  
  for (ki=0,sc1+=3,sc2+=3,sc3+=3; ki<icoef2; ki++)
     for (kj=(ki==0); kj<icoef1; kj++,sc1+=3,sc2+=3,sc3+=3)
     {
	
	/* Set correct pointers.  */ 
	
	if (((ki==0 || ki==icoef2-1) && icoef2>1) ||
		  ((kj==0 || kj==icoef1-1) && icoef1>1))
	   tmin = sminedg, tmax = smaxedg;
	else 
	   kinset = 1,  tmin = sminin,  tmax = smaxin;
	
	t1= ROTM * sc1[0];
	t2= ROTM * sc2[0];
	t3= ROTM * sc3[0];

	if(*sc1 < *tmin) *tmin = *sc1;
	if(*sc1 > *tmax) *tmax = *sc1;
	tmin++; tmax++;
	if(*sc2 < *tmin) *tmin = *sc2;
	if(*sc2 > *tmax) *tmax = *sc2;
	tmin++; tmax++;
	if(*sc3 < *tmin) *tmin = *sc3;
	if(*sc3 > *tmax) *tmax = *sc3;
	tmin++; tmax++;
	t4= t2 - t3;
	if(t4 < *tmin) *tmin = t4;
	if(t4 > *tmax) *tmax = t4;
	tmin++; tmax++;
	t4= t2 + t3;
	if(t4 < *tmin) *tmin = t4;
	if(t4 > *tmax) *tmax = t4;
	tmin++; tmax++;
	t4= t1 - t3;
	if(t4 < *tmin) *tmin = t4;
	if(t4 > *tmax) *tmax = t4;
	tmin++; tmax++;
	t4= t1 + t3;
	if(t4 < *tmin) *tmin = t4;
	if(t4 > *tmax) *tmax = t4;
	tmin++; tmax++;
	t4= t1 - t2;
	if(t4 < *tmin) *tmin = t4;
	if(t4 > *tmax) *tmax = t4;
	tmin++; tmax++;
	t4= t1 + t2;
	if(t4 < *tmin) *tmin = t4;
	if(t4 > *tmax) *tmax = t4;
     }
  
  /* Merge the inner and the outer box, and adjust with the
     tolerance.  */
  
  if (!kinset)
  {
     memcopy(sminin,sminedg,kant,DOUBLE);
     memcopy(smaxin,smaxedg,kant,DOUBLE);
  }
  for (ki=0; ki<kant; ki++)
  {
     e2min[ki] = MIN(sminin[ki]-aeps1,sminedg[ki]-teps2);
     e2max[ki] = MAX(smaxin[ki]+aeps1,smaxedg[ki]+teps2);
     e2min[kant+ki] = MIN(sminin[ki]-teps1,sminedg[ki]-teps3);
     e2max[kant+ki] = MAX(smaxin[ki]+teps1,smaxedg[ki]+teps3);
  }
}
 
#if defined(SISLNEEDPROTOTYPES)
static void
  sh1992_s9mbox2(double ecoef[],int icoef1,int icoef2,double aeps1,
		 double aeps2,double e2max[],double e2min[])
#else
static void sh1992_s9mbox2(ecoef,icoef1,icoef2,aeps1,aeps2,e2max,e2min)
     double ecoef[];
     int    icoef1;
     int    icoef2;
     double aeps1;
     double aeps2;
     double e2max[];
     double e2min[];
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make 2 boxes on the control-polygons given by
*              ecoef to the object.
*
*
* INPUT      : ecoef  - Control-polygon.
*              icoef1  - Number of vertices in control-polygon in
*                        first parameter direction.
*              icoef2  - Number of vertices in control-polygon in
*                        second parameter direction.
*              aeps1   - Size of expansion to use in the inner of
*                        the object.
*              aeps2   - Size of expansion to use at the edge of the object.
*
*                                                                     
*
* OUTPUT     : e2max   - Array to contain maximum values of the box.
*	       e2min   - Array to contain minimum values of the box.
*
*
* METHOD     : Make two boxes. One ordinary and one rotated 45 degree.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
* MODIFIED BY : Vibeke Skytt, SI, 91-01.
*
*********************************************************************
*/                                     
{
  int ki,kj;             /* Counters.                                 */
  int kant = 4;          /* Number of box sides.                      */
  int kinset = 0;        /* Indicates if an inner box is found.       */
  double teps1 = aeps1+aeps1; /* Double tolerance in the inner.       */
  double teps2;               /* Tolerance at edge.                   */
  double teps3;               /* Double tolerance at edge.            */
  double t1,t2,t3;       /* To store elements of the rotation matrix. */
  double *tmin,*tmax;    /* Pointers used to traverse e2min and e2max.  */
  double sminin[4],smaxin[4];   /* Box boundaries in the inner.       */
  double sminedg[4],smaxedg[4]; /* Box boundaries at the edge.        */
  double *sc1,*sc2;      /* Pointers into coefficient array.          */
  
  /* Set tolerances at edge. If the tolerance is positive or dimension
     is 1D, the input tolerance is used, otherwise we must make sure 
     that the maximum distance from the total box at the edges to the
     reduced box is aeps2.                                             */
  
  if (aeps2 >= DZERO)
     teps2 = aeps2;
  else
     teps2 = (double)0.38268343*aeps2;   /* aeps2 * sin(PI/8).   */
  teps3 = teps2 + teps2;
 
  /* Initiate box boundaries of inner box.  */
  
  for (ki=0; ki<kant; ki++)
  {
     sminin[ki] = HUGE;
     smaxin[ki] = -HUGE;
  }
  
  /* Fetch value of first vertex.  */
  
  sc1 = ecoef;  sc2 = sc1 + 1;
  t1= ROTM * sc1[0];
  t2= ROTM * sc2[0];
  
  tmin = sminedg;
  tmax = smaxedg;
  *tmin = *tmax = *sc1;
  tmin++; tmax++;
  *tmin = *tmax = *sc2;
  tmin++; tmax++;
  *tmin = *tmax = t1-t2;
  tmin++; tmax++;
  *tmin = *tmax = t1+t2;
  
  /* For each vertex check and corrigate the box.  */
  
  for (ki=0,sc1+=2,sc2+=2; ki<icoef2; ki++)
     /* UJK, writing error */
     /*for (kj=(ki==1); kj<icoef1; kj++,sc1+=2,sc2+=2) */

     for (kj=(ki==0); kj<icoef1; kj++,sc1+=2,sc2+=2)
     {
	/* Set correct box boundaries.  */
	
	if (((ki==0 || ki==icoef2-1) && icoef2>1) ||
		  ((kj==0 || kj==icoef1-1) && icoef1>1))
	   tmin = sminedg,  tmax = smaxedg;
	else
	   kinset = 1,  tmin = sminin,  tmax = smaxin;
	
	t1= ROTM * sc1[0];
	t2= ROTM * sc2[0];

	if(*sc1 < *tmin) *tmin = *sc1;
	if(*sc1 > *tmax) *tmax = *sc1;
	tmin++; tmax++;
	if(*sc2 < *tmin) *tmin = *sc2;
	if(*sc2 > *tmax) *tmax = *sc2;
	tmin++; tmax++;
	t3= t1 - t2;
	if(t3 < *tmin) *tmin = t3;
	if(t3 > *tmax) *tmax = t3;
	tmin++; tmax++;
	t3= t1 + t2;
	if(t3 < *tmin) *tmin = t3;
	if(t3 > *tmax) *tmax = t3;
     }
  
  /* Merge the inner and the outer box, and adjust with the
     tolerance.  */
  
  if (!kinset)
  {
     memcopy(sminin,sminedg,kant,DOUBLE);
     memcopy(smaxin,smaxedg,kant,DOUBLE);
  }
  for (ki=0; ki<kant; ki++)
  {
     e2min[ki] = MIN(sminin[ki]-aeps1,sminedg[ki]-teps2);
     e2max[ki] = MAX(smaxin[ki]+aeps1,smaxedg[ki]+teps2);
     e2min[kant+ki] = MIN(sminin[ki]-teps1,sminedg[ki]-teps3);
     e2max[kant+ki] = MAX(smaxin[ki]+teps1,smaxedg[ki]+teps3);
  }
}

#if defined(SISLNEEDPROTOTYPES)
static void
  sh1992_s9mbox(double ecoef[],int icoef1,int icoef2,int idim,
		double aeps1,double aeps2,double e2max[],
		double e2min[],int *jstat)
#else   
static void sh1992_s9mbox(ecoef,icoef1,icoef2,idim,aeps1,aeps2,
			  e2max,e2min,jstat)
     double ecoef[];
     int    icoef1;
     int    icoef2;
     int    idim;
     double aeps1;
     double aeps2;
     double e2max[];
     double e2min[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make a SISLbox on the control-polygons given by
*              ecoef to the object.
*
*
* INPUT      : ecoef  - Control-polygon.
*              icoef1  - Number of vertices in control-polygon in
*                        first parameter direction.
*              icoef2  - Number of vertices in control-polygon in
*                       second parameter direction.
*              aeps1   - Size of expansion to use in the inner of
*                        the object.
*              aeps2   - Size of expansion to use at the edge of the object.
*
*                                                                     
*
* OUTPUT     : e2max   - Array to contain maximum values of the box.
*	       e2min   - Array to contain minimum values of the box.
*              jstat  - status messages  
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
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
* MODIFIED BY : Vibeke Skytt, SI, 91-01.
*
*********************************************************************
*/                                     
{
  int ki,ki1,kj;       /* Counters.  */
  int kant = idim;     /* Number of box sides.                        */
  int kinset = 0;      /* Indicates if the inner box is set.          */
  double noice = (double)100.0*REL_COMP_RES;   /* Noice killer.       */
  double teps1 = aeps1+aeps1; /* Double tolerance in the inner.       */
  double teps2;               /* Tolerance at edge.                   */
  double teps3;               /* Double tolerance at edge.            */
  double *tmin,*tmax;  /* Pointers into box boundary arrays.          */
  double *sc;          /* Pointer into coefficient array.             */
  double *sminin=SISL_NULL,*smaxin=SISL_NULL;  /* Box boundaries of the inner.  */
  double *sminedg=SISL_NULL,*smaxedg=SISL_NULL; /* Box boundaries of the edge.  */
  
  /* Set tolerances at edge. If the tolerance is positive or dimension
     is 1D, the input tolerance is used, otherwise we must make sure 
     that the maximum distance from the total box at the edges to the
     reduced box is aeps2.                                             */
  
  if (idim == 1 || aeps2 >= DZERO)
     teps2 = aeps2;
  else
     teps2 = aeps2/sqrt((double)idim);
  teps3 = teps2 + teps2;
  
  /* Allocate scratch for intermediate box arrays.  */
  
  if ((sminin = newarray(kant,double)) == SISL_NULL) goto err101;
  if ((smaxin = newarray(kant,double)) == SISL_NULL) goto err101;
  if ((sminedg = newarray(kant,double)) == SISL_NULL) goto err101;
  if ((smaxedg = newarray(kant,double)) == SISL_NULL) goto err101;
  
  /* Initiate box boundaries of inner box.  */
  
  for (ki=0; ki<kant; ki++)
  {
     sminin[ki] = HUGE;
     smaxin[ki] = -HUGE;
  }
  
  /* Fetch value of first vertex.  */
  
  for (ki = 0; ki < idim; ki++) 
     sminedg[ki] = smaxedg[ki] = ecoef[ki];
  
  /* For each vertice check and corrigate the box.  */
  
  for (sc=ecoef+idim, ki=0; ki<icoef2; ki++)
     for (kj=(ki==0); kj<icoef1; kj++)
     {
	/* Set correct box.  */
	
	if (((ki==0 || ki==icoef2-1) && icoef2>1) ||
		  ((kj==0 || kj==icoef1-1) && icoef1>1))
	   tmin = sminedg,  tmax = smaxedg;
	else 
	    kinset = 1,  tmin = sminin,  tmax = smaxin;
	
	for (ki1=0; ki1<idim; ki1++,sc++,tmin++,tmax++)
	{
	   if(*sc < *tmin) *tmin = *sc;
	   if(*sc > *tmax) *tmax = *sc;
	}
     }

  /* Merge the inner and the outer box, and adjust with the
     tolerance.  */
  
  if (!kinset)
  {
     memcopy(sminin,sminedg,kant,DOUBLE);
     memcopy(smaxin,smaxedg,kant,DOUBLE);
  }
  for (ki=0; ki<kant; ki++)
  {
     e2min[ki] = MIN(sminin[ki]-aeps1,sminedg[ki]-teps2);
     e2max[ki] = MAX(smaxin[ki]+aeps1,smaxedg[ki]+teps2);
     if (idim > 1)
     {
	e2min[kant+ki] = MIN(sminin[ki]-teps1,sminedg[ki]-teps3);
	e2max[kant+ki] = MAX(smaxin[ki]+teps1,smaxedg[ki]+teps3);
     }	
  }
  
  /* ALA and UJK 30.10.90, remove noice near by zero.  */
  
  if (idim == 1)
  {
     if (fabs(e2max[0]) < noice) e2max[0] = DZERO;
     if (fabs(e2min[0]) < noice) e2min[0] = DZERO;
  }
  
  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation. */
  
  err101 : *jstat = -101;
  goto out;
  
  out :
  if (sminin != SISL_NULL) freearray(sminin);
  if (smaxin != SISL_NULL) freearray(smaxin);
  if (sminedg != SISL_NULL) freearray(sminedg);
  if (smaxedg != SISL_NULL) freearray(smaxedg);		       
}
 

