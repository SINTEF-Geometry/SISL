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
 * $Id: s1992.c,v 1.2 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S1992

#include "sislP.h"                                                 

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void
s1992_s9mbox3(double [],int,double [],double []);
static void
s1992_s9mbox2(double [],int,double [],double []);
static void
s1992_s9mbox(double [],int,int idim,double[],double []);
#else
static void s1992_s9mbox3();
static void s1992_s9mbox2();
static void s1992_s9mbox();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
s1992(SISLObject *po,int *jstat)
#else
void s1992(po,jstat)
     SISLObject *po;
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
* INPUT      : po     - SISLObject to treat.
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
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
*
*********************************************************************
*/                                     
{
  int kpos = 0;                        /* Position of error.   */

  if (po -> iobj == SISLPOINT)
    {
      if (po->p1->pbox == SISL_NULL)
	{
	  if ((po->p1->pbox = newbox(po->p1->idim))==SISL_NULL)
	    goto err101;
	  
	  if (po->p1->idim == 3) 
	    s1992_s9mbox3(po->p1->ecoef,1,po->p1->pbox->emax,
		    po->p1->pbox->emin);
	  else 
	    if (po->p1->idim == 2)
	      s1992_s9mbox2(po->p1->ecoef,1,po->p1->pbox->emax,
		      po->p1->pbox->emin);
	    else
	      s1992_s9mbox(po->p1->ecoef,1,po->p1->idim,
		     po->p1->pbox->emax,po->p1->pbox->emin);
	}
    }
  else
    if (po -> iobj == SISLCURVE)
      {
	if (po->c1->pbox == SISL_NULL)
	  {
	    if ((po->c1->pbox = newbox(po->c1->idim))==SISL_NULL)
	      goto err101;
	    
	    if (po->c1->idim == 3) 
	      s1992_s9mbox3(po->c1->ecoef,po->c1->in,po->c1->pbox->emax,
		      po->c1->pbox->emin);
	    else 
	      if (po->c1->idim == 2)
		s1992_s9mbox2(po->c1->ecoef,po->c1->in,po->c1->pbox->emax,
			po->c1->pbox->emin);
	      else
		s1992_s9mbox(po->c1->ecoef,po->c1->in,po->c1->idim,
		       po->c1->pbox->emax,po->c1->pbox->emin);
	  }
      }
    else
      if (po -> iobj == SISLSURFACE)
	{
	  if (po->s1->pbox == SISL_NULL)
	    {
	      if ((po->s1->pbox = newbox(po->s1->idim))==SISL_NULL)
		goto err101;
	      
	      if (po->s1->idim == 3) 
		s1992_s9mbox3(po->s1->ecoef,po->s1->in1 * po->s1->in2,
			po->s1->pbox->emax,po->s1->pbox->emin);
	      else 
		if (po->s1->idim == 2)
		  s1992_s9mbox2(po->s1->ecoef,po->s1->in1 * po->s1->in2,
			  po->s1->pbox->emax,po->s1->pbox->emin);
		else
		  s1992_s9mbox(po->s1->ecoef,po->s1->in1 * po->s1->in2,
			 po->s1->idim,po->s1->pbox->emax,po->s1->pbox->emin);
	    }
	}
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1992",*jstat,kpos);
  goto out;
  
 out: ;
}

#if defined(SISLNEEDPROTOTYPES)
void 
s1992cu(SISLCurve *pc,int *jstat)
#else
void s1992cu(pc,jstat)
     SISLCurve *pc;
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
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
*
*********************************************************************
*/                                     
{
  int kpos = 0;                        /* Position of error.   */
  
  if (pc->pbox == SISL_NULL)
    {
      if ((pc->pbox = newbox(pc->idim))==SISL_NULL)
	goto err101;
      
      if (pc->idim == 3) 
	s1992_s9mbox3(pc->ecoef,pc->in,pc->pbox->emax,pc->pbox->emin);
      else if (pc->idim == 2)
	s1992_s9mbox2(pc->ecoef,pc->in,pc->pbox->emax,pc->pbox->emin);
      else
	s1992_s9mbox(pc->ecoef,pc->in,pc->idim,
	       pc->pbox->emax,pc->pbox->emin);
    }
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1992cu",*jstat,kpos);
  goto out;
  
 out: ;
}

#if defined(SISLNEEDPROTOTYPES)
void 
s1992su(SISLSurf *ps,int *jstat)
#else
void s1992su(ps,jstat)
     SISLSurf *ps;
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
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-02.
*
*********************************************************************
*/                                     
{
  int kpos = 0;                        /* Position of error.   */
  
  if (ps->pbox == SISL_NULL)
    {
      if ((ps->pbox = newbox(ps->idim))==SISL_NULL) goto err101;
      
      if (ps->idim == 3) s1992_s9mbox3(ps->ecoef,ps->in1 * ps->in2,
				 ps->pbox->emax,ps->pbox->emin);
      
      else if (ps->idim == 2) s1992_s9mbox2(ps->ecoef,ps->in1 * ps->in2,
				      ps->pbox->emax,ps->pbox->emin);
      else
	s1992_s9mbox(ps->ecoef,ps->in1 * ps->in2,
	       ps->idim,ps->pbox->emax,ps->pbox->emin);
    }
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1992su",*jstat,kpos);
  goto out;
  
 out: ;
}

#if defined(SISLNEEDPROTOTYPES)
static void
  s1992_s9mbox3(double ecoef[],int icoef,double gmax[],double gmin[])
#else
static void s1992_s9mbox3(ecoef,icoef,gmax,gmin)
     double ecoef[];
     int    icoef;
     double gmax[];
     double gmin[];
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
* INPUT      : ecoef  - Control-polygon.
*              icoef - Number of vertices in control-polygon.
*
*                                                                     
*
* OUTPUT     : gmax   - Array to contain maximum values of the box.
*	       gmin   - Array to contain minimum values of the box.
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
*
*********************************************************************
*/                                     
{
  int ki,ki1,ki2,ki3;    /* Counters.                                 */
  double t1,t2,t3,t4;    /* To store elements of the rotation matrix. */
  double *tmin,*tmax;    /* Pointers used to traverse gmin and gmax.  */
  
  /* Fetch value of first vertex.  */
  
  t1= ROTM * ecoef[0];
  t2= ROTM * ecoef[1];
  t3= ROTM * ecoef[2];
  
  tmin = gmin;
  tmax = gmax;
  *tmin = *tmax = ecoef[0];
  tmin++; tmax++;
  *tmin = *tmax = ecoef[1];
  tmin++; tmax++;
  *tmin = *tmax = ecoef[2];
  tmin++; tmax++;
  *tmin = *tmax = ecoef[0];
  tmin++; tmax++;
  *tmin = *tmax = t2-t3;
  tmin++; tmax++;
  *tmin = *tmax = t2+t3;
  tmin++; tmax++;
  *tmin = *tmax = t1-t3;
  tmin++; tmax++;
  *tmin = *tmax = ecoef[1];
  tmin++; tmax++;
  *tmin = *tmax = t1+t3;
  tmin++; tmax++;
  *tmin = *tmax = t1-t2;
  tmin++; tmax++;
  *tmin = *tmax = t1+t2;
  tmin++; tmax++;
  *tmin = *tmax = ecoef[2];
  
  /* For each vertice check and corrigate the box.  */
  
  for (ki=1,ki1=3,ki2=4,ki3=5; ki<icoef; ki++,ki1+=3,ki2+=3,ki3+=3)
    {
      
      
      t1= ROTM * ecoef[ki1];
      t2= ROTM * ecoef[ki2];
      t3= ROTM * ecoef[ki3];
      tmin = gmin;
      tmax = gmax;
      if(ecoef[ki1] < *tmin) *tmin = ecoef[ki1];
      if(ecoef[ki1] > *tmax) *tmax = ecoef[ki1];
      tmin++; tmax++;
      if(ecoef[ki2] < *tmin) *tmin = ecoef[ki2];
      if(ecoef[ki2] > *tmax) *tmax = ecoef[ki2];
      tmin++; tmax++;
      if(ecoef[ki3] < *tmin) *tmin = ecoef[ki3];
      if(ecoef[ki3] > *tmax) *tmax = ecoef[ki3];
      tmin++; tmax++;
      if(ecoef[ki1] < *tmin) *tmin = ecoef[ki1];
      if(ecoef[ki1] > *tmax) *tmax = ecoef[ki1];
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
      if(ecoef[ki2] < *tmin) *tmin = ecoef[ki2];
      if(ecoef[ki2] > *tmax) *tmax = ecoef[ki2];
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
      tmin++; tmax++;
      if(ecoef[ki3] < *tmin) *tmin = ecoef[ki3];
      if(ecoef[ki3] > *tmax) *tmax = ecoef[ki3];
    }
}
 
#if defined(SISLNEEDPROTOTYPES)
static void
  s1992_s9mbox2(double ecoef[],int icoef,double gmax[],double gmin[])
#else
static void s1992_s9mbox2(ecoef,icoef,gmax,gmin)
     double ecoef[];
     int    icoef;
     double gmax[];
     double gmin[];
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
*              icoef - Number of vertices in control-polygon.
*
*                                                                     
*
* OUTPUT     : gmax   - Array to contain maximum values of the box.
*	       gmin   - Array to contain minimum values of the box.
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
*
*********************************************************************
*/                                     
{
  int ki,ki1,ki2;        /* Counters.                                 */
  double t1,t2,t3;       /* To store elements of the rotation matrix. */
  double *tmin,*tmax;    /* Pointers used to traverse gmin and gmax.  */
  
  
  
  /* Fetch value of first vertex.  */
  
  t1= ROTM * ecoef[0];
  t2= ROTM * ecoef[1];
  
  tmin = gmin;
  tmax = gmax;
  *tmin = *tmax = ecoef[0];
  tmin++; tmax++;
  *tmin = *tmax = ecoef[1];
  tmin++; tmax++;
  *tmin = *tmax = t1-t2;
  tmin++; tmax++;
  *tmin = *tmax = t1+t2;
  
  /* For each vertice check and corrigate the box.  */
  
  for (ki=1,ki1=2,ki2=3; ki<icoef; ki++,ki1+=2,ki2+=2)
    {
      
      
      t1= ROTM * ecoef[ki1];
      t2= ROTM * ecoef[ki2];
      tmin = gmin;
      tmax = gmax;
      if(ecoef[ki1] < *tmin) *tmin = ecoef[ki1];
      if(ecoef[ki1] > *tmax) *tmax = ecoef[ki1];
      tmin++; tmax++;
      if(ecoef[ki2] < *tmin) *tmin = ecoef[ki2];
      if(ecoef[ki2] > *tmax) *tmax = ecoef[ki2];
      tmin++; tmax++;
      t3= t1 - t2;
      if(t3 < *tmin) *tmin = t3;
      if(t3 > *tmax) *tmax = t3;
      tmin++; tmax++;
      t3= t1 + t2;
      if(t3 < *tmin) *tmin = t3;
      if(t3 > *tmax) *tmax = t3;
    }
}

#if defined(SISLNEEDPROTOTYPES) 
static void
  s1992_s9mbox(double ecoef[],int icoef,int idim,double gmax[],double gmin[])
#else
static void s1992_s9mbox(ecoef,icoef,idim,gmax,gmin)
     double ecoef[];
     int    icoef;
     int    idim;
     double gmax[];
     double gmin[];
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
*              icoef - Number of vertices in control-polygon.
*
*                                                                     
*
* OUTPUT     : gmax   - Array to contain maximum values of the box.
*	       gmin   - Array to contain minimum values of the box.
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
*
*********************************************************************
*/                                     
{
  int ki,ki1,kj;         /* Counters.  */
  double noice = (double)100.0 * REL_COMP_RES;   /* Noice killer */ 
  
  
  
  /* Fetch value of first vertex.  */
  
  for (ki = 0; ki < idim; ki++) gmin[ki] = gmax[ki] = ecoef[ki];
  
  /* For each vertice check and corrigate the box.  */
  
  for (kj=1; kj<icoef; kj++)
    for (ki1 = 0; ki1 < idim; ki1++,ki++) 
      {
	if(ecoef[ki] < gmin[ki1]) gmin[ki1] = ecoef[ki];
	if(ecoef[ki] > gmax[ki1]) gmax[ki1] = ecoef[ki];
      }

  /* ALA and UJK 30.10.90, remove noice near by zero */
  if (idim == 1)
    {
      if (fabs(gmax[0]) < noice) gmax[0] = DZERO; 
      if (fabs(gmin[0]) < noice) gmin[0] = DZERO; 
    }
  
}
 



