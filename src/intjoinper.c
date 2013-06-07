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
 * $Id: intjoinper.c,v 1.4 2001-03-19 16:13:07 afr Exp $
 *
 */


#define INT_JOIN_PER

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    int_join_per (SISLIntdat ** pintdat,
	    SISLObject * po1,
	    SISLObject * po2,
	    double eimpli[],
	    int ideg,
	    double aepsge,
	    int *jstat)

#else
void
   int_join_per (pintdat,
	    po1,
	    po2,
	    eimpli,
	    ideg,
	    aepsge,
	    jstat)

     SISLIntdat **pintdat;
     SISLObject *po1;
     SISLObject *po2;
     double eimpli[];
     int ideg;
     double aepsge;
     int *jstat;

#endif
/*
*********************************************************************
*
* PURPOSE    : To join intersection curves based on object space data only.
*              Curves are joined if they are G1 object space.
*              Single point with same object space posisition are
*              made to one point.
*
*
* INPUT      : pintdat     - Pointer to pointer to the SISLIntdat data.
*              po1         - Pointer surface object.
*              po2         - Pointer surface object.
*              eimpli      - Array containing descr. of implicit surf
*	       ideg        - Type of impl surf.
              ang_tol     - Angle control tolerance ie ??
*              aepsge      - Absolute tolerance
*
*
* OUTPUT     :  jstat  - status messages
*                       = ?      : ?
*                       = 0      : ok
*                       < 0      : error
*
*
* METHOD     :
*
* USE        : The function is only working i 3-D
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway, July-1990
*
*********************************************************************
*/
{
  int kstat;			/* Local status variable.          */
  int kpos = 0;			/* Position of error.              */
  int ki=0, kj=0;		        /* Counter                         */
  SISLIntpt *pcurr;		/* to traverse list of points.     */
  SISLIntpt *pnext;	        /* to traverse list of points.     */
  SISLIntpt *pother;	/* to traverse list of points.     */
  int logtest=0;                /* used for constant crv test      */

  SISLIntpt *pfirst_1;		/* First point in a list, periodicity     */
  SISLIntpt *pfirst_2;		/* First point in a list, periodicity     */
  SISLIntpt *plast_1;		/* Last point in a list, periodicity      */
  SISLIntpt *plast_2;		/* Last point in a list, periodicity      */
  double *curve_val_3d  = SISL_NULL; /* Pos, tang and curvature ,3d    */
  double *curve_val_2d_1= SISL_NULL; /* Pos, tang and curvature ,2d    */
  double *curve_val_2d_2= SISL_NULL; /* Pos, tang and curvature ,2d    */
  double delta_par[4];          /* Delta vector in par space, periodicity */
  double *delta_1=delta_par;    /* Delta vector in par space, periodicity */
  double *delta_2=delta_par+2;  /* Delta vector in par space, periodicity */
  double dot;                   /* Scalar product, periodicity            */
  double dist;                  /* Distance in 3D, periodicity            */
  double dist_a;                /* Distance in 3D, periodicity            */
  double dist_b;                /* Distance in 3D, periodicity            */
  double dist_c;                /* Distance in 3D, periodicity            */
  double dist_d;                /* Distance in 3D, periodicity            */
  double ang;                   /* Angel between tangents, periodicity    */
  int dimobj;                   /* Dimension 3, periodicity               */
  //int dim2=2;                   /* Dimension 2, periodicity               */
  int join=TRUE;                /* Flag to kill int point, periodicity    */
  int const_crv_1;              /* Constant parameter direction           */
  int const_crv_2;              /* Constant parameter direction           */
  int same_curve;               /* Flag                                   */
  int no_of_main;               /* No of mainpts connected to curr pt     */
  int treat_2d;                 /* Flag, shifting of parameter values is 
				   only done in surf/surf cases.          */
  int cas;                      /* Flag, surf surf, surf analytic, other. */ 
  double epar[2];
  double *nullp = SISL_NULL;
  SISLIntpt *pturn=SISL_NULL;	/* Last point in a list, periodicity      */
  int log_1, log_2;             /* To test on an edge curve lies along
				   the same parameter direction.          */
  int kp,no_par,index_1;
  double min_par[4],max_par[4],legal_min[4],legal_max[4];
  SISLObject *qo=SISL_NULL;
  /* -------------------------------------------------------------------- */ 
       
  /* If we do not have any intersection data we just return. */     
  
  if ((*pintdat) == SISL_NULL)
     goto out;
  if ((*pintdat)->ipoint <=1)
     goto out;
  
  
  if      (po1->iobj == SISLSURFACE) dimobj = po1->s1->idim;
  else if (po1->iobj == SISLCURVE)   dimobj = po1->c1->idim;
  else goto errinp;
  
  /* UJK, 13.08.93 : Check periodicity flag.*/
  no_par = (*pintdat)->vpoint[0]->ipar;
  if (no_par > 4) goto errinp;
  
  for (kp=0;kp<4;kp++)
  {
     legal_min[kp] = -HUGE;
     legal_max[kp] =  HUGE;
  }
  
  for (ki=0,qo=po1,index_1=0;ki<2;ki++,qo=po2)
     if (qo)
     {
	if (qo->iobj == SISLSURFACE) 
	{
	   if (qo->s1->cuopen_1 != SISL_SURF_PERIODIC)
	   {
	      legal_min[index_1] = qo->s1->et1[qo->s1->ik1-1];
	      legal_max[index_1] = qo->s1->et1[qo->s1->in1];
	   }
	   
	   index_1++;
	   if (qo->s1->cuopen_2 != SISL_SURF_PERIODIC)
	   {
	      legal_min[index_1] = qo->s1->et2[qo->s1->ik2-1];
	      legal_max[index_1] = qo->s1->et2[qo->s1->in2];
	   }
	   index_1++;
	}
	else if (qo->iobj == SISLCURVE) 
	{
	   if (qo->c1->cuopen != SISL_CRV_PERIODIC)
	   {
	      legal_min[index_1] = qo->c1->et[qo->c1->ik-1];
	      legal_max[index_1] = qo->c1->et[qo->c1->in];
	   }
	   index_1++;
	}
     }
  
  
  /*________________________________________ */
  /* UJK, TESTING !!!!!!!!!!!!!!!!!!! For the moment: */
  if      (dimobj == 3 && po2->iobj == SISLPOINT) goto out;
  /*________________________________________ */
  
  if (po1->iobj == SISLSURFACE && 
      po2->iobj == SISLSURFACE &&  ideg == 0)
     /* Two B-spline surf's in 3D. */
     cas = 1;
  else if (po1->iobj == SISLSURFACE && ideg != 2000 && ideg !=0) 
     /* B-spline surf vs analytic. */
     cas = 2;
  else 
     cas = 0;
  
  
  
  if (ideg != 2000 && po1->iobj == SISLSURFACE) treat_2d = TRUE;
  else treat_2d = FALSE;
  
  for (ki=1;ki<5;ki++) logtest |= 1<< ki;
  
  /* Ensure that object space information is in place. */
  for (kj = 0; kj < (*pintdat)->ipoint; kj++)
  {
     pcurr = (*pintdat)->vpoint[kj];
     sh6evalint (po1, po2, eimpli, ideg, pcurr, aepsge,
		 &curve_val_3d, &curve_val_2d_1,
		 &curve_val_2d_2, &kstat);
     if (kstat < 0)
	goto error;
  }
  
  /* Check direction for const_crvs. */
  if (cas)
     for (ki = 0; ki < (*pintdat)->ilist; ki++)
     {
	pfirst_1 = (*pintdat)->vlist[ki]->pfirst;
	plast_1  = (*pintdat)->vlist[ki]->plast;
	
	if(((*pintdat)->vlist[ki]->inumb == 2) &&
	   (pfirst_1->curve_dir[(*pintdat)->vlist[ki]->ind_first] & logtest))
	{
	   if (cas == 1)
	   {
	   }
	   else /*if (cas == 2)*/
	   {
	      /* Select midpoint. */
	      epar[0] = (pfirst_1->epar[0] + plast_1->epar[0])/2.0;
	      epar[1] = (pfirst_1->epar[1] + plast_1->epar[1])/2.0;
	      pturn = hp_newIntpt (2, epar, DZERO, SI_ORD,
				   SI_UNDEF, SI_UNDEF, SI_UNDEF, SI_UNDEF,
				   0, 0, nullp, nullp);
	      if (pturn == SISL_NULL)
		 goto err101;
	      
	      sh6evalint (po1, po2, eimpli, ideg, pturn, aepsge,
			  &curve_val_3d, &curve_val_2d_1,
			  &curve_val_2d_2, &kstat);
	      if (kstat < 0)
		 goto error;
	      
	      if (pturn->iinter == SI_ORD)
	      {
		 double dot;
		 dot = (plast_1->epar[0]-pfirst_1->epar[0])*curve_val_2d_1[2]+
		    (plast_1->epar[1]-pfirst_1->epar[1])*curve_val_2d_1[3];
		 if (dot < 0) 
		 {
		    /* Turn direction. */
		    int ind_1, ind_2, dir_1, dir_2;
		    
		    ind_1    = (*pintdat)->vlist[ki]->ind_first;
		    ind_2    = (*pintdat)->vlist[ki]->ind_last;
		    dir_1    = pfirst_1->curve_dir[ind_1];
		    dir_2    = plast_1->curve_dir[ind_2];
		    
		    (*pintdat)->vlist[ki]->pfirst = plast_1;
		    (*pintdat)->vlist[ki]->plast  = pfirst_1;
		    (*pintdat)->vlist[ki]->ind_first = ind_2;
		    (*pintdat)->vlist[ki]->ind_last  = ind_1;
		    pfirst_1->curve_dir[ind_1] = dir_2;
		    plast_1->curve_dir[ind_2]  = dir_1;
		    
		 }
	      }		 
	      if (pturn) freeIntpt(pturn);
	      pturn = SISL_NULL;
	      
	   }
	}
     }
  
  /* Traverse the lists to remove doubly represented curves along
     periodic edges.     */
  
  for (ki = 0; ki < (*pintdat)->ilist; ki++)
  {
     pfirst_1 = (*pintdat)->vlist[ki]->pfirst;
     plast_1  = (*pintdat)->vlist[ki]->plast;
     
     if(((*pintdat)->vlist[ki]->inumb == 2) &&
	(pfirst_1->curve_dir[(*pintdat)->vlist[ki]->ind_first] & logtest))
	const_crv_1 = TRUE;
     else 
	const_crv_1 = FALSE;
     
     
     if (pfirst_1 == plast_1) continue;
     
     for (kj = 0; kj < (*pintdat)->ilist; kj++)
     {
	
	same_curve = (ki == kj); 
	pfirst_2 = (*pintdat)->vlist[kj]->pfirst;
	plast_2 = (*pintdat)->vlist[kj]->plast;
	
	if(((*pintdat)->vlist[kj]->inumb == 2) &&
	   (pfirst_2->curve_dir[(*pintdat)->vlist[kj]->ind_first] & logtest))
	   const_crv_2 = TRUE;
	else 
	   const_crv_2 = FALSE;
	
	/* To treat the case when two curves are on an periodic edge */
	if (const_crv_1 && const_crv_2 && !same_curve)
	{
	   log_1 = pfirst_1->curve_dir[(*pintdat)->vlist[ki]->ind_first];
	   log_1 = log_1>>1;
	   log_1 &= 15;
	   log_2 = pfirst_2->curve_dir[(*pintdat)->vlist[kj]->ind_first];
	   log_2 = log_2>>1;
	   log_2 &= 15;
	   
	   if (log_1 & log_2)
	   {
	      dist_a = s6dist(plast_1->geo_track_3d,
			      pfirst_2->geo_track_3d,
			      dimobj);
	      dist_b = s6dist(plast_1->geo_track_3d,
			      plast_2->geo_track_3d,
			      dimobj);
	      dist_c = s6dist(pfirst_1->geo_track_3d,
			      pfirst_2->geo_track_3d,
			      dimobj);
	      dist_d = s6dist(pfirst_1->geo_track_3d,
			      plast_2->geo_track_3d,
			      dimobj);
	      if ((dist_a <aepsge && dist_d < aepsge) ||	
		  (dist_b <aepsge && dist_c < aepsge)	)
	      {
		 /* Kill the two points */
		 sh6idkpt(pintdat, &pfirst_2, join=FALSE, &kstat);
		 if (kstat < 0) goto error;
		 
		 sh6idkpt(pintdat, &plast_2, join=FALSE, &kstat);
		 if (kstat < 0) goto error;
		 
		 /* Remove the curve list kj, and pack vlist array */
		 freeIntlist ((*pintdat)->vlist[kj]);
		 (*pintdat)->ilist--;
		 (*pintdat)->vlist[kj] =
		    (*pintdat)->vlist[(*pintdat)->ilist];
		 
		 /* Reset to start */
		 ki = -1;
		 break; /* kj loop */
	      }
	      
	   }
	}
     }
  }
  
  for (ki = 0; ki < (*pintdat)->ilist; ki++)
  {
     pfirst_1 = (*pintdat)->vlist[ki]->pfirst;
     plast_1  = (*pintdat)->vlist[ki]->plast;
     
     if(((*pintdat)->vlist[ki]->inumb == 2) &&
	(pfirst_1->curve_dir[(*pintdat)->vlist[ki]->ind_first] & logtest))
	const_crv_1 = TRUE;
     else 
	const_crv_1 = FALSE;
     
     
     if (pfirst_1 == plast_1) continue;
     
     for (kj = 0; kj < (*pintdat)->ilist; kj++)
     {
	
	same_curve = (ki == kj); 
	pfirst_2 = (*pintdat)->vlist[kj]->pfirst;
	plast_2 = (*pintdat)->vlist[kj]->plast;
	
	if(((*pintdat)->vlist[kj]->inumb == 2) &&
	   (pfirst_2->curve_dir[(*pintdat)->vlist[kj]->ind_first] & logtest))
	   const_crv_2 = TRUE;
	else 
	   const_crv_2 = FALSE;
	
	/* Joining of curves */
	if (plast_1 == pfirst_2 ||
	    plast_1->iinter == SI_TRIM ||
	    plast_1->iinter == SI_SING ||
	    plast_1->iinter == SI_TOUCH ||
	    pfirst_2->iinter == SI_TRIM ||
	    pfirst_2->iinter == SI_SING ||
	    pfirst_2->iinter == SI_TOUCH)
	{
	   /* Test on edge ? */
	   continue;
	}
	else
	{
	   dist = s6dist(plast_1->geo_track_3d,
			 pfirst_2->geo_track_3d,
			 dimobj);
	   ang  = s6ang(plast_1->geo_track_3d + dimobj,
			pfirst_2->geo_track_3d + dimobj,
			dimobj);
	   dot  = s6scpr(plast_1->geo_track_3d + dimobj,
			 pfirst_2->geo_track_3d + dimobj,
			 dimobj);
	   
	   if (dist < aepsge &&
	       dot  > 0      &&
	       ang  < ANGULAR_TOLERANCE)
	   { 
	      /* HIT, Join the two lists of curves */
	      if (same_curve)
	      {
		 /*sh6connect(pfirst_1, plast_1, &kstat);
		    if (kstat != 0) goto error;
		    sh6idkpt (pintdat, &plast_1, join=TRUE, &kstat);
		    if (kstat < 0) goto error;
		    
		    (*pintdat)->vlist[ki]->inumb -= 1; 
		    
		    periodic set to true ! 
		    
		    (*pintdat)->vlist[ki]->plast = 
		    pfirst_1; */
	      }
	      else
	      {
		 /* Move 2D values of second part */
		 s6diff(plast_1->epar,
			pfirst_2->epar,
			no_par,
			delta_par);
		 
		 /* Last check to see if we move outside legal area */
		 for(kp=0;kp<no_par;kp++)
		 {
		    min_par[kp] = HUGE;
		    max_par[kp] = -HUGE;
		 }
		 
		 pcurr = pfirst_2;
		 pnext = pfirst_2->pnext[(*pintdat)->vlist[kj]->ind_first];
		 while (pcurr != plast_2 && pcurr && pnext)
		 {
		    for(kp=0;kp<no_par;kp++)
		    {
		       min_par[kp] = min(pnext->epar[kp]+delta_par[kp],min_par[kp]);
		       max_par[kp] = max(pnext->epar[kp]+delta_par[kp],max_par[kp]);
		    }
		    
		    sh6getother (pnext, pcurr, &pother, &kstat);
		    if (kstat && pnext != plast_2) goto errinconsist;
		    pcurr = pnext;
		    pnext = pother;
		 }		     
		 
		 for(kp=0;kp<no_par;kp++)
		 {
		    if (min_par[kp] < legal_min[kp] && 
			DNEQUAL(min_par[kp],legal_min[kp])) break;
		    if (max_par[kp] > legal_max[kp] && 
			DNEQUAL(max_par[kp],legal_max[kp])) break;
		    
		 }
		 
		 
		 /* _______________________ */
		 if (kp == no_par)
		 {
		    pcurr = pfirst_2;
		    pnext = pfirst_2->pnext[(*pintdat)->vlist[kj]->ind_first];
		    while (pcurr != plast_2 && pcurr && pnext)
		    {
		       for(kp=0;kp<no_par;kp++)
			  pnext->epar[kp] += delta_par[kp];
		       
		       if (treat_2d)
		       {
			  pnext->geo_track_2d_1[0] += delta_1[0];
			  pnext->geo_track_2d_1[1] += delta_1[1];
			  if (cas==1)
			  {
			     pnext->geo_track_2d_2[0] += delta_2[0];
			     pnext->geo_track_2d_2[1] += delta_2[1];
			  }
		       }
		       sh6getother (pnext, pcurr, &pother, &kstat);
		       if (kstat && pnext != plast_2) goto errinconsist;
		       pcurr = pnext;
		       pnext = pother;
		    }		     
		 
		 sh6connect(plast_1, pfirst_2, &kstat);
		 if (kstat != 0) goto error;
		 sh6idkpt (pintdat, &pfirst_2, join=TRUE, &kstat);
		 if (kstat < 0) goto error;
		 
		 if (const_crv_1 && const_crv_2)
		 {
		    sh6idkpt (pintdat, &plast_1, join=TRUE, &kstat);
		    if (kstat < 0) goto error;
		    
		    sh6getlist(pfirst_1,plast_2,
			       &((*pintdat)->vlist[ki]->ind_first),
			       &((*pintdat)->vlist[ki]->ind_last),
			       &kstat);
		    kpos = 111;
		    if (kstat != 0) goto errinconsist;
		    (*pintdat)->vlist[ki]->inumb = 2; 
		    
		    (*pintdat)->vlist[ki]->plast = 
		       plast_2;
		    
		 }
		 else 
		 {
		    
		    (*pintdat)->vlist[ki]->inumb += 
		       (*pintdat)->vlist[kj]->inumb -1;
		    
		    (*pintdat)->vlist[ki]->plast = 
		       plast_2;
		    (*pintdat)->vlist[ki]->ind_last = 
		       (*pintdat)->vlist[kj]->ind_last;
		 }
		 
		 /* Remove the curve list kj, and pack vlist array */
		 freeIntlist ((*pintdat)->vlist[kj]);
		 (*pintdat)->ilist--;
		 (*pintdat)->vlist[kj] =
		    (*pintdat)->vlist[(*pintdat)->ilist];
		 
		 
		 /* Reset to start */
		 ki = -1;
		 break; /* kj loop */
		 }
	      } /* End of else not same curve */
	      
	   } /* if hit */
	   
	} /* else */
     } /* kj */
  } /* ki */
  
  
  /* Treat single points */
  for (ki = 0; ki < (*pintdat)->ipoint; ki++)
  {
     kstat = 0;
     pcurr = (*pintdat)->vpoint[ki];
     if (sh6ismain(pcurr)) no_of_main = sh6nmbmain(pcurr, &kstat);
     else no_of_main = -1;
     if (kstat < 0) goto error;
     
     if (no_of_main == 0)		  
	for (kj = 0; kj < (*pintdat)->ipoint; kj++)
	{
	   pother = (*pintdat)->vpoint[kj];
	   if (sh6ismain(pother) && pother != pcurr)
	   {
	      dist = s6dist(pcurr->geo_track_3d,
			    pother->geo_track_3d,
			    dimobj);
	      if (dist < aepsge)
	      {
		 sh6idkpt(pintdat, &pcurr, join = FALSE, &kstat);
		 if (kstat < 0) goto error;
		 ki--; /* New point in array place no ki, redo the job */
		 break; /* kj loop */
	      }
	   }
	}
  }
  
  
  
  *jstat = 0;
  goto out;

  /* _________________________ EXIT ____________________________ */
  /* Error in alloc */
err101:
  *jstat = -101;
  s6err ("int_join_per", *jstat, kpos);
  goto out;

  /* Error inconsistency */
errinconsist:
  *jstat = -500;
  s6err ("int_join_per", *jstat, kpos);
  goto out;

  /* Error in input */
errinp:
     *jstat = -200;
  s6err ("int_join_per", *jstat, kpos);
  goto out;

  /* Error in lower level function */
error:
     *jstat = kstat;
  s6err ("int_join_per", *jstat, kpos);
  goto out;

out:
  ;
}
