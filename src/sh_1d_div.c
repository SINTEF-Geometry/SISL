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
 * $Id: sh_1d_div.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH_1D_DIV

#include "sislP.h"


/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void sh_1d_div_sh9idnpt(SISLSurf*, SISLPoint*,SISLIntdat **,
			       SISLIntpt **,int, double,int *);
#else
static void sh_1d_div_sh9idnpt();
#endif


#if defined(SISLNEEDPROTOTYPES)
void
  sh_1d_div (SISLObject *po1, SISLObject *po2, double aepsge,
	     SISLIntdat **pintdat,  SISLEdge * vedge[], int *jstat)
#else
void 
   sh_1d_div (po1, po2, aepsge, pintdat, vedge, jstat)
   SISLObject *po1; 
   SISLObject *po2; 
   double aepsge;
   SISLIntdat **pintdat;
   SISLEdge * vedge[];
   int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE     :To check if we can subdivide a 1D surface problem
*              by a factorization method.
*
*
*
*
* INPUT      : po1          - Pointer to surface object.
*              po2          - Pointer to point object.
*              aepsge       - Geometry tolerance.
*              vedge[2]  - Pointers to structure of edge-intersections.
*
*
* INPUT/OUTPUT : pintdat - Pointer to intersection data.
*
* OUTPUT     : jstat     - status messages
*                                           1      : Factorization done
*                                                    together with intersection.
*                                         = 0      : No factorization.
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
*
* WRITTEN BY : Ulf J. Krystad, SI, 92-12.
* MODIFIED BY :
*
**********************************************************************/
{

   int kant;                    /* Number of parameter directions          */
   int cv_dir_1, cv_dir_2;      /* Locals for curve_dir                    */
   int ind1, ind2;              /* Locals indexes                          */
   int kpos = 0;		/* Position of error.                      */
   int kstat = 0;		/* Local error status.                     */
   int knum;                    /* Number of intersection pts on edge      */              
   int knum2;                   /* Number of intersection pts in corners   */              
   int which_end_1=-1;          /* Branch paraeter fro zero edge.          */
   int which_end_2=-1;          /* Branch paraeter fro zero edge.          */
   int kn, kj, ki;              /* Loop control                            */
   int edge_1=0, edge_2=0;      /* No of pts in pt_arr_1[2]                */
   int alloc_1=0, alloc_2=0;    /* Size of pt_arr_1[2]                     */
   SISLIntpt *pcurr = SISL_NULL;	/* Array of poiners to int points.         */
   SISLIntpt **uintpt = SISL_NULL;	/* Array of poiners to int points.         */
   SISLIntpt **up = SISL_NULL;	/* Array of poiners to edge int points.    */
   SISLIntpt **pt_arr_1 = SISL_NULL;	/* Array of poiners to ZERO edge.          */
   SISLIntpt **pt_arr_2 = SISL_NULL;	/* Array of poiners to ZERO edge.          */
   SISLIntpt **up2 = SISL_NULL;	/* Array of poiners to corner int points.  */
   SISLIntdat *qintdat = SISL_NULL;	/* Data structure of sub inters problem    */
   SISLObject *qo1 = SISL_NULL;      /* Pointer to surface in
				   object/point intersection. */
   double *nullp = SISL_NULL;
   /* ____________________________________________________________________ */

   *jstat = 0;
   
   /* Check input */
   if (po1->iobj != SISLSURFACE) goto err150;
   if (po1->s1->idim != 1) goto err150;
   
   /* Bezier case ? */
   if (po1->s1->ik1 != po1->s1->in1 ||
       po1->s1->ik2 != po1->s1->in2) goto out;

   if (po1->s1->ik1 < 3 ||
       po1->s1->ik2 < 3 ) goto out;
   
   sh6edgpoint (vedge, &up, &knum, &kstat);
   if (kstat < 0)
      goto error;
   if (knum < 2) goto out;
   
   /* Find corner points */
   /* Allocate an array for intersection points. */
   if ((up2= newarray (knum, SISLIntpt *)) == SISL_NULL)
      goto err101;
   
   for (knum2=ki=0;ki<knum;ki++)
   {
      sh6isinside (po1, po2, up[ki], &kstat);
      if (kstat < 0 ) goto error;
      
      if (kstat == 3)
      {
	 up2[knum2] = up[ki];
	 knum2++;
      }
   }
   
   if (knum2 < 2) goto out;
   
   /* Find connections */
   for (ki=0;ki<knum2-1;ki++)
      for (kj=1;kj<knum2;kj++)
      {
	 sh6comedg (po1, po2, up2[ki], up2[kj], &kstat);
	 if (kstat < 0) goto error;
	 
	 if (kstat == 1)
	 {
	    /* One edge is zero, find which */
	    if (DEQUAL(up2[ki]->epar[0], up2[kj]->epar[0]))
	    {
	       
	       /* Store the two corner points (sorted). */
	       if (alloc_1 == 0)
	       {
		  alloc_1 = 10;
		  if((pt_arr_1 = newarray(alloc_1,SISLIntpt *))
		     == SISL_NULL) goto err101;
	       }
	       edge_1 = 2;
	       if (up2[ki]->epar[1] < up2[kj]->epar[1])
	       {		   
		  pt_arr_1[0] = up2[ki];
		  pt_arr_1[1] = up2[kj];
	       }
	       else
	       {		   
		  pt_arr_1[1] = up2[ki];
		  pt_arr_1[0] = up2[kj];
	       }
	       
	       if (DEQUAL(up2[ki]->epar[0],po1->s1->et1[0]))
		  which_end_1 = 0;
	       else
		  which_end_1 = 1;
	    }
	    else
	    {
	       
	       /* Store the two corner points (sorted). */
	       if (alloc_2 == 0)
	       {
		  alloc_2 = 10;
		  if((pt_arr_2 = newarray(alloc_2,SISLIntpt *))
		     == SISL_NULL) goto err101;
	       }
	       edge_2 = 2;
	       if (up2[ki]->epar[0] < up2[kj]->epar[0])
	       {		   
		  pt_arr_2[0] = up2[ki];
		  pt_arr_2[1] = up2[kj];
	       }
	       else
	       {		   
		  pt_arr_2[1] = up2[ki];
		  pt_arr_2[0] = up2[kj];
	       }
	       
	       if (DEQUAL(up2[ki]->epar[1],po1->s1->et2[0]))
		  which_end_2 = 0;
	       else
		  which_end_2 = 1;
	    }
	    
	 }
      }
   
   if (which_end_1 >=0 || which_end_2 >=0)
   {
      
      /*
      * Create new object and create surface to object.
      * ------------------------------------------------
      */
      
      if (!(qo1 = newObject (SISLSURFACE)))
	 goto err101;
      qo1->s1 = SISL_NULL;
      qo1->o1 = qo1;
      
      /* Filter coefficients less than aepsge. */
      for (ki=0; ki< po1->s1->in1*po1->s1->in2;ki++)
	 if ( fabs(po1->s1->ecoef[ki]-po2->p1->ecoef[0]) < aepsge)
	    po1->s1->ecoef[ki] = po2->p1->ecoef[0];
      
      sh_div_surf(po1->s1,which_end_1, which_end_2, aepsge, &qo1->s1, &kstat);
      if (kstat < 0) goto error;
      
      sh1761 (qo1, po2, aepsge, &qintdat, &kstat);
      if (kstat < 0)
	 goto error;
      
      /* UJK, JUNE 93: start____________ */
      if (qintdat)
      {
	 
	 /* Kill all help.pts. */
	 for (ki = 0; ki < qintdat->ipoint; ki++)
	 {
	    pcurr = qintdat->vpoint[ki];
	    if(sh6ishelp(pcurr))
	    {
	       sh6idkpt (&qintdat, &pcurr, 0, &kstat);
	       if (kstat < 0) goto error;
	       ki--;
	    }
	 }
      }
      /* UJK, JUNE 93: end____________ */
      
      if (qintdat)
      {
	 /* Intersection found, transfere it to pintdat. */
	 
	 /* Number of parameter direction. */
	 kant = qintdat->vpoint[0]->ipar;
	 
	 /* Allocate an array for intersection points. */
	 if ((uintpt = newarray (qintdat->ipoint, SISLIntpt *)) == SISL_NULL)
	    goto err101;
	 
	 /* Copy all intersection points. */
	 for (ki = 0; ki < qintdat->ipoint; ki++)
	 {
	    
	    uintpt[ki] = hp_newIntpt (kant,  
				      qintdat->vpoint[ki]->epar, 
				      qintdat->vpoint[ki]->adist,
				      qintdat->vpoint[ki]->iinter,
				      qintdat->vpoint[ki]->left_obj_1[0],
				      qintdat->vpoint[ki]->right_obj_1[0],
				      qintdat->vpoint[ki]->left_obj_2[0],
				      qintdat->vpoint[ki]->right_obj_2[0],
				      0, 0,
				      nullp, nullp);
	    
	    if (uintpt[ki] == SISL_NULL)
	       goto err101;
	 }
	 
	 /* Insert all new intersection points in rintdat. */

	 for (ki = 0; ki < qintdat->ipoint; ki++)
	 {
	    sh_1d_div_sh9idnpt (po1->o1->s1,po2->p1,pintdat, &uintpt[ki], 1,aepsge,
				&kstat);
	    if (kstat < 0)
	       goto error;
	 }

	 
	 /* Insert points on edges divided out/(splitting strategy. */
	 for (ki = 0; ki < qintdat->ipoint; ki++)
	 {
	    if (which_end_1 >=0 &&
		DEQUAL(uintpt[ki]->epar[0], pt_arr_1[0]->epar[0]) &&
		DEQUAL(uintpt[ki]->epar[0], pt_arr_1[1]->epar[0]))
	    {
	       for (kj=0; kj < edge_1 - 1; kj++)
		  if (uintpt[ki]->epar[1] > pt_arr_1[kj]->epar[1] &&
		      uintpt[ki]->epar[1] < pt_arr_1[kj+1]->epar[1])
		  {
		     sh6insertpt(pt_arr_1[kj],pt_arr_1[kj+1],uintpt[ki],&kstat);
	             if (kstat < 0) goto error;
		     if (edge_1 >= alloc_1)
		     {
			alloc_1 += 10;
			if ((pt_arr_1 = 
			     increasearray(pt_arr_1,alloc_1,SISLIntpt *))
			    == SISL_NULL) goto err101;
		     }
		     
		     for (kn = edge_1; kn > kj+1; kn--)
			pt_arr_1[kn] = pt_arr_1[kn-1];
		     pt_arr_1[kj+1] = uintpt[ki];
		     edge_1++;
		     
		     break;
		  }
	    }
	    else if (which_end_2 >=0 &&
		DEQUAL(uintpt[ki]->epar[1], pt_arr_2[0]->epar[1]) &&
		DEQUAL(uintpt[ki]->epar[1], pt_arr_2[1]->epar[1]))
	    {
	       for (kj=0; kj < edge_2 - 1; kj++)
		  if (uintpt[ki]->epar[0] > pt_arr_2[kj]->epar[0] &&
		      uintpt[ki]->epar[0] < pt_arr_2[kj+1]->epar[0])
		  {
		     sh6insertpt(pt_arr_2[kj],pt_arr_2[kj+1],uintpt[ki],&kstat);
	             if (kstat < 0) goto error;
		     if (edge_2 >= alloc_2)
		     {
			alloc_2 += 10;
			if ((pt_arr_2 = 
			     increasearray(pt_arr_2,alloc_2,SISLIntpt *))
			    == SISL_NULL) goto err101;
		     }
		     
		     for (kn = edge_2; kn > kj+1; kn--)
			pt_arr_2[kn] = pt_arr_2[kn-1];
		     pt_arr_2[kj+1] = uintpt[ki];
		     edge_2++;
		     
		     break;
		  }
	    }
	    
	 }	      
	    
	    
	 /* Transform the connections. */
	 for (ki = 0; ki < qintdat->ipoint; ki++)
	 {
	    for (kj = ki + 1; kj < qintdat->ipoint; kj++)
	    {
	       sh6getlist (qintdat->vpoint[ki], qintdat->vpoint[kj],
			   &ind1, &ind2, &kstat);
	       if (kstat < 0)
		  goto error;
	       if (kstat == 0 && uintpt[kj] != uintpt[ki])
	       {
                  cv_dir_1 = qintdat->vpoint[ki]->curve_dir[ind1];
                  cv_dir_2 = qintdat->vpoint[kj]->curve_dir[ind2];

		  sh6idcon (pintdat, &uintpt[ki], &uintpt[kj], &kstat);
		  if (kstat < 0)
		    goto error;
	       
		  sh6getlist (uintpt[ki], uintpt[kj],
			   &ind1, &ind2, &kstat);
		  if (kstat != 0) goto error;
		  uintpt[ki]->curve_dir[ind1] |= cv_dir_1;
		  uintpt[kj]->curve_dir[ind2] |= cv_dir_2;
		  
	       }
	    }
	    
	    if (sh6ismain (qintdat->vpoint[ki]) &&
		sh6nmbmain (qintdat->vpoint[ki], &kstat))
	    {
	       sh6tomain (uintpt[ki], &kstat);
	       if (kstat < 0)
		  goto error;
	    }
	 }
	 
	 /* Remains splitting of zero curves inside edges ! */	 
	 
      }

      *jstat = 1;
    }
   
   
   
   goto out;
   /* ______________ ERROR EXITS ______________________________ */
/* Lower level problem. */
error:
   *jstat = kstat;
   s6err("sh_1d_div",*jstat,kpos);
   goto out;
   
/* Space problem. */
err101:
   *jstat = -101;
   s6err("sh_1d_div",*jstat,kpos);
   goto out;
   
/* Input wrong. */
err150:
   *jstat = -150;
   s6err("sh_1d_div",*jstat,kpos);
   goto out;


out:
   if (uintpt) freearray(uintpt);
   if (up) freearray(up);
   if (up2) freearray(up2);
   if (pt_arr_1) freearray(pt_arr_1);
   if (pt_arr_2) freearray(pt_arr_2);
   if (qo1)
      freeObject (qo1);
   if (qintdat)
      freeIntdat (qintdat);
 }


#if defined(SISLNEEDPROTOTYPES)
static void 
  sh_1d_div_sh9idnpt(SISLSurf* surf, SISLPoint* point, SISLIntdat **pintdat,
		     SISLIntpt **pintpt, int itest, double aepsge, int *jstat)
#else
static void
  sh_1d_div_sh9idnpt(surf,point,pintdat,pintpt,itest, aepsge, jstat)
     SISLSurf* surf;
     SISLPoint* point;
     SISLIntdat **pintdat;
     SISLIntpt  **pintpt;
     int    itest;
     double aepsge;
     int    *jstat;
#endif   


/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To insert a new intersection point into pintdat.
*              If pintdat is SISL_NULL a new pintdat is also made.
*              If pintpt is close to an other intersection point
*              the object pintpt is pointing to is freed, and
*              pintpt is set to point to the already inserted point.
*
*
*
* INPUT      : pintpt   - Pointer to a pointer to new intersection point.
*              pintdat  - Pointer to a pointer to intersection data.
*              itest    - Indikate testing equalety.
*                               = 1      : Testing.
*                               = 0      : No testing.
*              surf	- pointer to the surface object in the intersection
*              point    - pointer to the point object in the intersection
*              aepsge  -  Resolution in space.
*
*
* OUTPUT     : jstat  - status messages  
*                               = 2      : Already existing.
*                               = 1      : Already inserted.
*                               = 0      : Intersection point inserted.
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
*              newIntdat  - Create new intdat structure.
*              freeIntpt  - free instant of intpt structure.
*
* WRITTEN BY : Michael Floater, June 91.
* REVISED BY : Kyrre Stroem, Jan-93.
*
*********************************************************************
*/                                     
{
  register int ki;              /* Counters.    */
  double eps_ball_par;
  int kstat;
  
  /* We have to be sure that we have an intdat structure. */
  
  if ((*pintdat) == SISL_NULL)
    {
      if (((*pintdat) = newIntdat()) == SISL_NULL) goto err101;
    }
  
  
  /* Then we have to be sure that we do not have the intersection point
     before or an equal point. */
  
  for (ki=0; ki<(*pintdat)->ipoint; ki++)
    if ((*pintdat)->vpoint[ki] == (*pintpt))
      {
	*jstat = 1;
	goto out;
      }
    else if (itest)
      {
	eps_ball_par = surf->et1[surf->in1]- surf->et1[surf->ik1];
	eps_ball_par = max(eps_ball_par,
			   surf->et2[surf->in2]- surf->et2[surf->ik2])+1;
	eps_ball_par *= 1e-6;

	s6identify(surf,(*pintpt)->epar,
		   (*pintdat)->vpoint[ki]->epar,
		   point->ecoef[0],eps_ball_par,aepsge,&kstat);
	if (kstat < 0 ) goto error;
	if (kstat == 1 ) 
	  {
	    freeIntpt(*pintpt);
	    (*pintpt) = (*pintdat)->vpoint[ki];
	    *jstat = 2;
	    goto out;
	  }
      }
  
  
  /* Then we have to be sure that the array vpoint is great enough. */
  
  if (ki == (*pintdat)->ipmax)
    {
      (*pintdat)->ipmax += 20;
      
      if (((*pintdat)->vpoint = increasearray((*pintdat)->vpoint,
					      (*pintdat)->ipmax,SISLIntpt *)) == SISL_NULL) 
	goto err101;
    }
  
  
  /* Now we can insert the new point. */
  
  (*pintdat)->vpoint[ki] = (*pintpt);
  (*pintdat)->ipoint++;
  *jstat = 0;
  goto out;
  

/* Error in space allocation.  */

err101: *jstat = -101;
        s6err("sh_1d_div_sh9idnpt",*jstat,0);
        goto out;
error: *jstat = kstat;
        s6err("sh_1d_div_sh9idnpt",*jstat,0);
        goto out;

 out: ;
}














