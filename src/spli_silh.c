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
 * $Id: spli_silh.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SPLI_SILH

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    spli_silh (SISLIntdat ** pintdat,
	    SISLObject * po1,
	    int *jstat)

#else
void
   spli_silh (pintdat,
	    po1,
	    jstat)

     SISLIntdat **pintdat;
     SISLObject *po1;
     int *jstat;

#endif
/*
*********************************************************************
*
* PURPOSE    : To split intersection curves at knot values where
*              they are only C0, only silhuett case.
*              (easy to change to ordinary intersections when
*              only C0 geometry is allowed.)
*
*
* INPUT      : pintdat     - Pointer to pointer to the SISLIntdat data.
*              po1         - Pointer surface object.
*
*
* OUTPUT     :  jstat  - status messages
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
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway, January-1993
*
*********************************************************************
*/
{
  int kstat;			/* Local status variable.          	  */
  int kpos = 0;			/* Position of error.              	  */
  int ktype = 0;		/* Type of intersection curve.    	  */
  int klist = 0;		/* No of lists in pintdat.         	  */
  int indstart = 0;		/* Index.                      	          */
  int indlast = 0;		/* Index.                      	          */
  int index1 = 0;		/* Index.                      	          */
  int indum = 0;		/* Not used index.            	          */
  int ki=0, kp=0;		/* Loop controls.                         */
  int mu_1=0, mu_2=0;		/* Knot multiplicity.                     */
  int left1=0,left2=0;		/* Knot navigators.                       */
  int ac_no=0;		        /* No of points in first half of split crv*/
  int knum=0;		        /* No of points in crv to split.          */
  int test=0;		        /* Flag input to s6idnpt.                 */
  SISLIntpt *pfirst=SISL_NULL;	/* First point in a list.                 */
  SISLIntpt *plast=SISL_NULL;	/* Last point in a list.                  */
  SISLIntpt *plast_prev=SISL_NULL;	/* Next to last point in a list.          */
  SISLIntpt *pshadow=SISL_NULL;	/* Copy of pnext.                         */
  SISLIntpt *pnext=SISL_NULL;	/* Navigator in  a list.             	  */
  SISLIntpt *prev=SISL_NULL;	 	/* Navigator in  a list.             	  */
  SISLIntpt *pother=SISL_NULL;	/* Navigator in  a list.             	  */
  SISLSurf  *psurf=SISL_NULL;	/* The surface in question.            	  */
  double arr[2];            	/* Delta vector in par space              */
  /* -------------------------------------------------------------------- */ 
       
  *jstat = 0;
  
  /* If we do not have any intersection data we just return. */     
  if ((*pintdat) == SISL_NULL)
    goto out;
  
  if (po1->iobj != SISLSURFACE) goto errinp;
  psurf = po1->s1;
   

  /* First pass, examine for split pts. */
  for (ki = 0; ki < (*pintdat)->ilist; ki++)
  {
     knum     = (*pintdat)->vlist[ki]->inumb;
     pfirst   = (*pintdat)->vlist[ki]->pfirst;
     plast    = (*pintdat)->vlist[ki]->plast;
     indstart = (*pintdat)->vlist[ki]->ind_first;
     indlast  = (*pintdat)->vlist[ki]->ind_last;
     ktype    = (*pintdat)->vlist[ki]->itype;
     prev     = pfirst;
     pnext    = pfirst->pnext[indstart];
     ac_no    = 1;

     /* VSK, 23/2-93. First check if the list describes an area of
	coincidence.  */
     
     if (pfirst->iinter == SI_TRIM) continue;
     
     if (pfirst == plast &&
	 knum > 3)
     {
	mu_1 = s6knotmult(psurf->et1,psurf->ik1,
			   psurf->in1,&left1,pfirst->epar[0], &kstat);
	mu_2 = s6knotmult(psurf->et2,psurf->ik2,
			  psurf->in2,&left2,pfirst->epar[1], &kstat);
	
	/* Always seperate */
	/* if(mu_1 > psurf->ik1 - 3 ||
	   mu_2 > psurf->ik2 - 3) 
	   { */
	
	/* Don't allow closed ring. */
	kpos = 10;
	plast_prev = plast->pnext[indlast];
	pshadow    = hp_copyIntpt(plast);
	if (!pshadow) goto err101;
	sh6idnpt(pintdat, &pshadow, test=FALSE, &kstat);
	if (kstat < 0) goto error;
	
	
	sh6insertpt(plast_prev, plast, pshadow, &kstat);
	if (kstat < 0) goto error;
	
	sh6disconnect(plast, pshadow, &kstat);
	if (kstat < 0) goto error;
	sh6getlist(pshadow, plast_prev, &index1, &indum, &kstat);
	kpos = 11;
	if (kstat < 0 || index1 < 0) goto errinconsis;
	
	plast    = (*pintdat)->vlist[ki]->plast = pshadow;
	indlast  = (*pintdat)->vlist[ki]->ind_last = index1;

	/*UJK, July 93 */
	sh6getlist(pfirst, pnext, &index1, &indum, &kstat);
	kpos = 12;
	if (kstat < 0 || index1 < 0) goto errinconsis;
	indstart = (*pintdat)->vlist[ki]->ind_first = index1;
	
	/*	} */
     }
     
     if (knum > 2)
     {
	while (pnext != plast)
	{
	   ac_no++;
	   kpos = 1;
	   if (ac_no >= knum) goto errinconsis;
	   sh6getother(pnext,prev,&pother,&kstat);
	   if (kstat < 0) goto error;
	   
	   mu_1 = s6knotmult(psurf->et1,psurf->ik1,
			      psurf->in1,&left1,pnext->epar[0], &kstat);
	   mu_2 = s6knotmult(psurf->et2,psurf->ik2,
			      psurf->in2,&left2,pnext->epar[1], &kstat);
	   
	   if(mu_1 > psurf->ik1 - 3 ||
	      mu_2 > psurf->ik2 - 3)
	   {

	      /* Split list */
	      pshadow = hp_copyIntpt(pnext);
	      if (!pshadow) goto err101;
	      sh6idnpt(pintdat, &pshadow, test=FALSE, &kstat);
	      if (kstat < 0) goto error;
	      
	      sh6insertpt(prev, pnext, pshadow, &kstat);
	      if (kstat < 0) goto error;
	      
	      sh6disconnect(pnext, pshadow, &kstat);
	      if (kstat < 0) goto error;
	      
	      /* Update first part of list. */
	      sh6getlist(pshadow, prev, &index1, &indum, &kstat);
	      kpos = 3;
	      if (kstat < 0 || index1 < 0) goto errinconsis;
	      
	      (*pintdat)->vlist[ki]->inumb    = ac_no;
	      (*pintdat)->vlist[ki]->ind_last = index1;
	      (*pintdat)->vlist[ki]->plast    = pshadow;
	      
	      /* Make and update second part of list. */
	      sh6getlist(pnext, pother, &index1, &indum, &kstat);
	      kpos = 4;
	      if (kstat < 0 || index1 < 0) goto errinconsis;
	      
	      klist = (*pintdat)->ilist;
	      kpos = 5;
	      if (klist > (*pintdat)->ilmax) goto errinconsis;
	      if (klist == (*pintdat)->ilmax)
	      {
		 /* Increase size of vlist-array */
		 (*pintdat)->ilmax += 20;
		 
		 if (((*pintdat)->vlist =
		      increasearray ((*pintdat)->vlist, (*pintdat)->ilmax,
				     SISLIntlist *)) == SISL_NULL)
		    goto err101;
	      }
	      
	      
	      if(((*pintdat)->vlist[klist] = newIntlist (pnext, plast, ktype))
		 == SISL_NULL) goto err101;
	      (*pintdat)->vlist[klist]->inumb = knum + 1 - ac_no;
	      (*pintdat)->vlist[klist]->ind_first = index1;
	      (*pintdat)->vlist[klist]->ind_last  = indlast;
	      for (kp = 0; kp < 4; kp++)
		 (*pintdat)->vlist[klist]->pretop[kp] =
		    (*pintdat)->vlist[ki]->pretop[kp];
	      
	      (*pintdat)->ilist++;
	      break;
	   }
	   else
	   {
	      prev = pnext;
	      pnext = pother;
	   }
	   
	} 
	
     }
  }
  
  /* Second pass, select left/right evaluator at start and stop pts. */
  for (ki = 0; ki < (*pintdat)->ilist; ki++)
  {
     pfirst   = (*pintdat)->vlist[ki]->pfirst;
     plast    = (*pintdat)->vlist[ki]->plast;
     indstart = (*pintdat)->vlist[ki]->ind_first;
     indlast  = (*pintdat)->vlist[ki]->ind_last;

     /* Startpoint. */
     pnext    = pfirst->pnext[indstart];
     s6diff(pnext->epar,pfirst->epar,2,arr);
     pfirst->iside_1 = (arr[0]>=0 ? 1 : -1);
     pfirst->iside_2 = (arr[1]>=0 ? 1 : -1);
     
     /* Endpoint. */
     pnext    = plast->pnext[indlast];
     s6diff(pnext->epar,plast->epar,2,arr);
     plast->iside_1 = (arr[0]>=0 ? 1 : -1);
     plast->iside_2 = (arr[1]>=0 ? 1 : -1);
     
     
     /* UJK, July 93: Clean up computed results */
     if (pfirst->geo_data_1)
	freearray (pfirst->geo_data_1);
     if (pfirst->geo_data_2)
	freearray (pfirst->geo_data_2);
     pfirst->geo_data_1 = SISL_NULL;
     pfirst->size_1 = 0;
     pfirst->geo_data_2 = SISL_NULL;
     pfirst->size_2 = 0;
     pfirst->evaluated = FALSE;
     
     if (plast->geo_data_1)
	freearray (plast->geo_data_1);
     if (plast->geo_data_2)
	freearray (plast->geo_data_2);
     plast->geo_data_1 = SISL_NULL;
     plast->size_1 = 0;
     plast->geo_data_2 = SISL_NULL;
     plast->size_2 = 0;
     plast->evaluated = FALSE;
  }  
    
     
  
  *jstat = 0;
  goto out;
  /* ERROR EXITS. ____________________________________ */
  
  /* Error in alloc */
  err101:
  *jstat = -101;
  s6err ("spli_silh", *jstat, kpos);
  goto out;

  /* Error in input */
  errinp:
  *jstat = -200;
  s6err ("spli_silh", *jstat, kpos);
  goto out;

  /* Error data structure */
  errinconsis:
  *jstat = -200;
  s6err ("spli_silh", *jstat, kpos);
  goto out;

  error:
  *jstat = kstat;
  s6err ("spli_silh", *jstat, kpos);
  goto out;

out:
  ;
}
