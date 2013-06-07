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
 * $Id: sh6idlis.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH6IDLIS

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
sh6idlis (SISLObject * po1, SISLObject * po2, SISLIntdat ** pintdat,
	  double aepsge, int *jstat)
#else
void
sh6idlis (po1, po2, pintdat, aepsge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat **pintdat;
     double aepsge;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To update vlist in pintdat.
*
*
* INPUT      : pintdat  - Pointer to a pointer to intersection data.
*
*
* OUTPUT     : jstat  - status messages
*                               = 0      : OK !
*                               < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6err       - Gives error message.
*              newIntlist  - Create new intlist structure.
*              freeIntlist - Free an intlist structure.
*
* WRITTEN BY :Ulf J. Krystad, 06.90.
*
*********************************************************************
*/
{
  int kstat;			/* Local status variable.          */
  int kpos = 0;			/* Position of error.              */
  int list_index = 0;		/* Counter                         */
  int knum = 0;			/* Counter                         */
  int indstart;			/* Indexes used in lists           */
  int indlast;			/* Indexes used in lists           */
  int ind1, ind2;
  int inddum;			/* Indexes used in lists           */
  int no_main = 0;		/* Counter                         */
  int ki1, ki2, ki, kj;		/* Counters                        */
  int r1, r2, l1, l2;		/* Pretopology info.		   */
  int ktype = 0;		/* To indicate type of list.       */
  int direction;		/* Direction of curve              */
  double *geom, *norm1, *norm2;	/* help pointers.		   */
  SISLIntpt *prev, *pcurr;	/* to traverse list of points.     */
  SISLIntpt *pnext, *pstart;	/* to traverse list of points.     */
  SISLIntpt *plast, *pother;	/* to traverse list of points.     */
  int pretop[4];
  int case_2d = 0;		/* Case flag, 2d Sf vs Pnt.        */
  int const_dir;		/* Reduction of internal points
				   along a constant parameter.     */
  int log_1, log_2;
  /* ------------------------------------------------------------- */

  /* If we do not have any intersection data we just return. */

  if ((*pintdat) == SISL_NULL)
    goto out;
  if ((po1->iobj == SISLSURFACE && po1->s1->idim == 2) ||
      (po2->iobj == SISLSURFACE && po2->s1->idim == 2))
    case_2d = TRUE;
  else
    case_2d = FALSE;

  /* We first destroy existing intersection lists. */

  for (kj = 0; kj < (*pintdat)->ilist; kj++)
    freeIntlist ((*pintdat)->vlist[kj]);

  /* Set SI_AT info in topology part */
  sh_set_at (po1, po2, *pintdat, &kstat);
  if (kstat < 0)
    goto error;

  
  
  /* Traverse all intersection points to get pretopology from help
     points */
  
  for (ki1 = 0; *pintdat && ki1 < (*pintdat)->ipoint; ki1++)
  {
     pcurr = (*pintdat)->vpoint[ki1];
     if (sh6ismain (pcurr))
     {
	sh6gettop (pcurr, 0, pretop, pretop + 1, pretop + 2, pretop + 3, &kstat);
	if (kstat < 0)
	   goto error;
	
	for (ki2 = 0; ki2 < pcurr->no_of_curves; ki2++)
	{
	   sh6gettophlp (pcurr->pnext[ki2], pretop, case_2d, &kstat);
	   if (kstat < 0)
	      goto error;
	}
	
	sh6settop (pcurr, 0, pretop[0], pretop[1], pretop[2], pretop[3], &kstat);
	if (kstat < 0)
	   goto error;
	
     }
  }
  
  
  /* Remove all internal points in a list when along a
     constant parameter direction */
  /* Remove all singularpoints that has exactly two 
     singular neighbours. */
  
  if ((po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT
       && po1->s1->idim == 1) ||
      (po2->iobj == SISLSURFACE && po1->iobj == SISLPOINT
       && po2->s1->idim == 1) ||
      (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE
       && po1->s1->idim == 3))
     const_dir = 2;
  else if ((po1->iobj == SISLSURFACE && po2->iobj == SISLCURVE) ||
	   (po1->iobj == SISLCURVE && po2->iobj == SISLSURFACE))
     const_dir = 1;
  else
     const_dir = 0;
  
  for (kj = 0; kj < (*pintdat)->ipoint; kj++)
  {
     
     pcurr = (*pintdat)->vpoint[kj];
     
     sh6getnhbrs (pcurr, &pstart, &plast, &kstat);
     if (kstat < 0)
	goto error;
     
     if (kstat == 0)
     {
	/* Two neighbours, check */
	sh6getlist (pcurr, pstart, &indstart, &inddum, &kstat);
	if (kstat < 0)
	   goto error;		/* Error. */
	if (kstat == 1)
	   goto errinconsist;	/* pcurr and pstart are not linked. */
	
	sh6getlist (pcurr, plast, &indlast, &inddum, &kstat);
	if (kstat < 0)
	   goto error;		/* Error. */
	if (kstat == 1)
	   goto errinconsist;	/* pcurr and plast are not linked. */
	
	log_1 = pcurr->curve_dir[indstart];
	log_1 = log_1>>1;
	log_1 &= 15;
	log_2 = pcurr->curve_dir[indlast];
	log_2 = log_2>>1;
	log_2 &= 15;
	
	
	if (const_dir == 0 || 
	    (log_1 & log_2 ) ||
	    (pcurr->iinter == SI_SING &&
	     pstart->iinter == SI_SING && plast->iinter == SI_SING))
	{
	   sh6idkpt (pintdat, &pcurr, 1, &kstat);
	   if (kstat < 0)
	      goto error;
	   /* Recursive nature : */
	   kj = -1;
	}
	
	
     }
  }
  
  /* -------------------------------------------- */
  if (const_dir > 1)
  {
     /* Split curves at points when change in curve_dir */
     for (kj = 0; kj < (*pintdat)->ipoint; kj++)
     {
	
	pcurr = (*pintdat)->vpoint[kj];
	
	sh6getnhbrs (pcurr, &pstart, &plast, &kstat);
	if (kstat < 0)
	   goto error;
	
	if (pcurr->iinter == SI_ORD && kstat == 0)
	{
	   /* Two neighbours, check */
	   sh6getlist (pcurr, pstart, &indstart, &inddum, &kstat);
	   if (kstat < 0)
	      goto error;		/* Error. */
	   if (kstat == 1)
	      goto errinconsist;	/* pcurr and pstart are not linked. */
	   
	   sh6getlist (pcurr, plast, &indlast, &inddum, &kstat);
	   if (kstat < 0)
	      goto error;		/* Error. */
	   if (kstat == 1)
	      goto errinconsist;	/* pcurr and plast are not linked. */
	   
	   log_1 = pcurr->curve_dir[indstart];
	   log_1 = log_1>>1;
	   log_1 &= 15;
	   log_2 = pcurr->curve_dir[indlast];
	   log_2 = log_2>>1;
	   log_2 &= 15;
	   
	   /* If both curve_dir is set as constant, this must be a singular
	      point, (remember internal points on same edge has been
	      removed! )*/
	   if (log_1 && log_2 ) pcurr->iinter = SI_SING;
	}
	
	
     }
  }
  
  /* -------------------------------------------- */
  
  if (const_dir > 1)
  {
     
     /* Split curves at corner points. This is put in to avoid
	problems for the marching.   */
     
     for (kj = 0; kj < (*pintdat)->ipoint; kj++)
     {
	
	pcurr = (*pintdat)->vpoint[kj];
	
	pcurr->marker = FALSE;
	if (pcurr->iinter == SI_ORD && sh6nmbmain(pcurr,&kstat) == 2)
	   
	{
	   sh6isinside(po1,po2,pcurr,&kstat);
	   
	   /* UJK, February 1993, Sometimes an intersection point on an edge is not
	      identified as singular. Therefore we split the curve at edges, 
	      if this is no natural ending (ie parallel pnt), 
	      int_join_per will join them. */
	   /* if (kstat == 3 || kstat == 4) */
	   
	   if (kstat < 0) goto error;
	   
	   if (kstat > 1) 
	   {
	      /* The point lies on a boarder. Mark it
		 achieve a split.  */
	      
	      /* UJK, aug 93, always mark singular in corners */
	      if (kstat == 3 || kstat == 4) pcurr->iinter = SI_SING;
	      pcurr->marker = TRUE;
	   }
	}
     }
  }
  
  
  if (const_dir > 1)
  {
     /* All previous trim points in a corner, must split it's neighbour if
	this lies on a corner */
     for (kj = 0; kj < (*pintdat)->ipoint; kj++)
     {
	
	pcurr = (*pintdat)->vpoint[kj];
	
	if (pcurr->iinter == SI_TRIM &&
	    sh6nmbmain (pcurr, &kstat) == 1)
	{
	   sh6isinside (po1, po2, pcurr, &kstat);
	   if (kstat == 3 || kstat == 4)
	   {
	      sh6getnhbrs (pcurr, &pstart, &plast, &kstat);
	      if (kstat < 0)
		 goto error;
	      if (pstart->iinter == SI_TRIM)
	      {
		 sh6isinside (po1, po2, pstart, &kstat);
		 if (kstat == 3 || kstat == 4)
		 {
		    
		    pstart->iinter = SI_SING;
		    
		    /* Recursive nature : */
		    kj = -1;
		 }
	      }
	   }
	}
     }
  }
  
  /* April 92, we need one instanse of each end point. This
     because of the storing of geometric data in the points
     that is in some cases contex dependent (mirroring
     in singular situation or translation of parameter space 
     values in periodicity treatment).
     We identify junction points and make a copy
     for each branch. */
  for (ki1 = 0; *pintdat && ki1 < (*pintdat)->ipoint; ki1++)
  {
     pcurr = (*pintdat)->vpoint[ki1];
     if (sh6ismain (pcurr))
     {
	/* Get number of neighbours */
	no_main = sh6nmbmain (pcurr, &kstat);
	if (kstat < 0)
	   goto error;
	
	if (pcurr->marker ||
	    (no_main == 2 && pcurr->iinter == SI_SING) ||
	    no_main > 2)
	{
	   if (pcurr->iinter == SI_ORD && no_main > 2) 
	      pcurr->iinter = SI_SING;
	   sh6idsplit(pintdat, pcurr, &kstat);
	   if (kstat < 0) goto error;
	}
     }
  }
  /* End of split */
  
  /* Traverse all intersection points, mark all main points */
  for (ki1 = 0; *pintdat && ki1 < (*pintdat)->ipoint; ki1++)
     if (sh6ismain ((*pintdat)->vpoint[ki1]))
     {
	/* Get number of neighbours */
	(*pintdat)->vpoint[ki1]->marker =
	   sh6nmbmain ((*pintdat)->vpoint[ki1], &kstat);
	if (kstat < 0)
	   goto error;
     }
     else
	(*pintdat)->vpoint[ki1]->marker = 0;
     
     
     /* Traverse all intersection points to look for
	start points to lists. If a point has only one neighbour,
	or is SI_SING or has more than two neighbours,
	it is a start or end point. */
     
     for (ki1 = 0; *pintdat && ki1 < (*pintdat)->ipoint; ki1++)
	if ((*pintdat)->vpoint[ki1]->marker > 0)
	{
	   /* Get number of neighbours */
	   no_main = sh6nmbmain ((*pintdat)->vpoint[ki1], &kstat);
	   if (kstat < 0)
	      goto error;
	   
	   if (no_main == 1 ||
	       (no_main == 2 && (*pintdat)->vpoint[ki1]->iinter == SI_SING) ||
	       no_main > 2)
	   {
	      pstart = (*pintdat)->vpoint[ki1];
	      
	      for (ki2 = 0; ki2 < pstart->no_of_curves; ki2++)
		 if (sh6ismain (pstart->pnext[ki2]) &&
		     pstart->pnext[ki2]->marker)
		 {
		    pcurr = pstart->pnext[ki2];
		    prev = pstart;
		    knum = 1;
		    
		    /* Get first index */
		    sh6getlist (pstart, pcurr, &indstart, &inddum, &kstat);
		    
		    while (pcurr)
		    {
		       /* Remember index */
		       sh6getlist (prev, pcurr, &inddum, &indlast, &kstat);
		       
		       prev->marker--;
		       pcurr->marker--;
		       sh6getother (pcurr, prev, &pnext, &kstat);
		       if (kstat < 0)
			  goto error;
		       
		       prev = pcurr;
		       pcurr = pnext;
		       knum++;
		    }
		    
		    /* Create list */
		    /* To be sure that list array is big enough. */
		    
		    if (list_index == (*pintdat)->ilmax)
		    {
		       (*pintdat)->ilmax += 20;
		       
		       if (((*pintdat)->vlist =
			    increasearray ((*pintdat)->vlist, (*pintdat)->ilmax,
					   SISLIntlist *)) == SISL_NULL)
			  goto err101;
		    }
		    
		    /* Type setting may be done in s1880? */
		    ktype = 0;
		    
		    /* Making a new list structure. */
		    if (((*pintdat)->vlist[list_index] =
			 newIntlist (pstart, prev, ktype)) == SISL_NULL)
		       goto err101;
		    (*pintdat)->vlist[list_index]->inumb = knum;
		    (*pintdat)->vlist[list_index]->ind_first = indstart;
		    (*pintdat)->vlist[list_index]->ind_last = indlast;
		    list_index++;
		    
		    
		    
		 }
	   }
	}
     
     /* Only closed list left */
     for (ki1 = 0; *pintdat && ki1 < (*pintdat)->ipoint; ki1++)
	if ((*pintdat)->vpoint[ki1]->marker > 0)
	{
	   /* Get number of neighbours */
	   no_main = sh6nmbmain ((*pintdat)->vpoint[ki1], &kstat);
	   if (kstat < 0)
	      goto error;
	   
	   
	   if (no_main == 2)
	   {
	      pstart = prev = (*pintdat)->vpoint[ki1];
	      
	      sh6getnhbrs (prev, &pcurr, &pnext, &kstat);
	      if (kstat < 0 || pcurr == SISL_NULL)
		 goto error;
	      knum = 1;
	      
	      /* Get first index */
	      sh6getlist (prev, pcurr, &indstart, &inddum, &kstat);
	      
	      /* Get last index */
	      sh6getlist (prev, pnext, &indlast, &inddum, &kstat);
	      
	      
	      
	      while (pcurr && pcurr != pstart)
	      {
		 prev->marker = 0;
		 
		 sh6getother (pcurr, prev, &pnext, &kstat);
		 if (kstat < 0)
		    goto error;
		 
		 prev = pcurr;
		 pcurr = pnext;
		 knum++;
	      }
	      prev->marker = 0;
	      if (pcurr == pstart)
	      {
		 /* It really is a closed curve */
		 knum++;
		 prev = pstart;
	      }
	      else
		 goto errinconsist;
	      
	      /* Create list */
	      /* To be sure that list array is big enough. */
	      
	      if (list_index == (*pintdat)->ilmax)
	      {
		 (*pintdat)->ilmax += 20;
		 
		 if (((*pintdat)->vlist =
		      increasearray ((*pintdat)->vlist, (*pintdat)->ilmax,
				     SISLIntlist *)) == SISL_NULL)
		    goto err101;
	      }
	      
	      /* Type setting may be done in s1880? */
	      ktype = 0;
	      
	      /* Making a new list structure. */
	      if (((*pintdat)->vlist[list_index] =
		   newIntlist (pstart, prev, ktype)) == SISL_NULL)
		 goto err101;
	      (*pintdat)->vlist[list_index]->inumb = knum;
	      (*pintdat)->vlist[list_index]->ind_first = indstart;
	      (*pintdat)->vlist[list_index]->ind_last = indlast;
	      list_index++;
	      
	      
	   }
	}
     
     
     (*pintdat)->ilist = list_index;
     
     
     /* If direction of a list is wrong, turn it. */
     
     for (kj = 0; kj < (*pintdat)->ilist; kj++)
     {
	knum = (*pintdat)->vlist[kj]->inumb;
	indstart = (*pintdat)->vlist[kj]->ind_first;
	
	pcurr = (*pintdat)->vlist[kj]->pfirst;
	plast = (*pintdat)->vlist[kj]->plast;
	pnext = (*pintdat)->vlist[kj]->pfirst->pnext[indstart];
	direction = (*pintdat)->vlist[kj]->pfirst->curve_dir[indstart];
	
	(*pintdat)->vlist[kj]->pretop[0] = SI_UNDEF;
	(*pintdat)->vlist[kj]->pretop[1] = SI_UNDEF;
	(*pintdat)->vlist[kj]->pretop[2] = SI_UNDEF;
	(*pintdat)->vlist[kj]->pretop[3] = SI_UNDEF;
	
	
	while (pnext != plast)
	{
	   if (direction)
	      break;
	   
	   sh6getother (pnext, pcurr, &pother, &kstat);
	   if (kstat < 0)
	      goto error;
	   sh6getlist (pnext, pother, &ind1, &ind2, &kstat);
	   if (kstat < 0)
	      goto error;
	   direction = pnext->curve_dir[ind1];
	   pcurr = pnext;
	   pnext = pother;
	}
	
	if (direction < 0)
	{
	   pcurr = (*pintdat)->vlist[kj]->pfirst;
	   (*pintdat)->vlist[kj]->pfirst = (*pintdat)->vlist[kj]->plast;
	   (*pintdat)->vlist[kj]->plast = pcurr;
	   
	   inddum = (*pintdat)->vlist[kj]->ind_first;
	   (*pintdat)->vlist[kj]->ind_first = (*pintdat)->vlist[kj]->ind_last;
	   (*pintdat)->vlist[kj]->ind_last = inddum;
	   indstart = (*pintdat)->vlist[kj]->ind_first;
	}
	
	/* Set pretopology information. */
	/* Traverse the list */
	
	if ((po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT
	     && po1->s1->idim == 1) ||
	    (po2->iobj == SISLSURFACE && po1->iobj == SISLPOINT
	     && po2->s1->idim == 1) ||
	    (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE
	     && po1->s1->idim == 3))
	{
	   pstart = prev = (*pintdat)->vlist[kj]->pfirst;
	   plast = (*pintdat)->vlist[kj]->plast;
	   
	   if ((pstart->edge_1 && plast->edge_1) ||
	       (pstart->edge_2 && plast->edge_2))
	   {
	      l1 = SI_IN;
	      r1 = SI_OUT;
	      l2 = SI_OUT;
	      r2 = SI_IN;
	      
	      if (plast->edge_1 == SI_RIGHT)
		 r1 = SI_AT;
	      else if (plast->edge_1 == SI_LEFT)
		 l1 = SI_AT;
	      
	      if (plast->edge_2 == SI_RIGHT)
		 r2 = SI_AT;
	      else if (plast->edge_2 == SI_LEFT)
		 l2 = SI_AT;
	      
	      (*pintdat)->vlist[kj]->pretop[0] = l1;
	      (*pintdat)->vlist[kj]->pretop[1] = r1;
	      (*pintdat)->vlist[kj]->pretop[2] = l2;
	      (*pintdat)->vlist[kj]->pretop[3] = r2;
	      
	   }
	   
	   else
	   {
	      
	      r1 = r2 = l1 = l2 = 0;
	      pcurr = sh6getnext (prev, indstart);
	      
	      /* UJK,Does not work
		 while (pcurr != plast) */
	      while (0)
	      {
		 if (sh6nmbhelp (pcurr, &kstat))
		 {
		    if (pcurr->left_obj_1[0] == pcurr->right_obj_1[0] &&
			pcurr->left_obj_1[0] != 0 && r1 == 0)
		    {
		       l1 = r1 = pcurr->left_obj_1[0];
		       if (po1->iobj == 2 && po2->iobj == 2 && po1->s1->idim == 3)
		       {
			  sh6getgeom (po1, 1, pcurr, &geom, &norm1, aepsge, &kstat);
			  if (kstat < 0)
			     goto error;
			  sh6getgeom (po2, 2, pcurr, &geom, &norm2, aepsge, &kstat);
			  if (kstat < 0)
			     goto error;
			  
			  if (s6scpr (norm1, norm2, 3) < 0.0)
			     l2 = r2 = l1;
			  else
			     l2 = r2 = (l1 == SI_IN ? SI_OUT : SI_IN);
		       }
		       else
			  l2 = r2 = (l1 == SI_IN ? SI_OUT : SI_IN);
		       
		       break;
		    }
		    
		    if (pcurr->left_obj_2[0] == pcurr->right_obj_2[0] &&
			pcurr->left_obj_2[0] != 0 && r2 == 0)
		    {
		       l2 = r2 = pcurr->left_obj_2[0];
		       if (po1->iobj == 2 && po2->iobj == 2 && po1->s1->idim == 3)
		       {
			  sh6getgeom (po1, 1, pcurr, &geom, &norm1, aepsge, &kstat);
			  if (kstat < 0)
			     goto error;
			  sh6getgeom (po2, 2, pcurr, &geom, &norm2, aepsge, &kstat);
			  if (kstat < 0)
			     goto error;
			  
			  if (s6scpr (norm1, norm2, 3) < 0.0)
			     l1 = r1 = l2;
			  else
			     l1 = r1 = (l2 == SI_IN ? SI_OUT : SI_IN);
		       }
		       else
			  l1 = r1 = (l2 == SI_IN ? SI_OUT : SI_IN);
		       
		       break;
		    }
		 }
		 else if (kstat < 0)
		    goto error;
		 
		 
		 sh6getother (pcurr, prev, &pnext, &kstat);
		 if (kstat < 0)
		    goto error;
		 
		 prev = pcurr;
		 pcurr = pnext;
	      }
	      if (r1 == 0)
	      {
		 l1 = SI_IN;
		 r1 = SI_OUT;
		 l2 = SI_OUT;
		 r2 = SI_IN;
	      }
	      
	      
	      (*pintdat)->vlist[kj]->pretop[0] = l1;
	      (*pintdat)->vlist[kj]->pretop[1] = r1;
	      (*pintdat)->vlist[kj]->pretop[2] = l2;
	      (*pintdat)->vlist[kj]->pretop[3] = r2;
	      
	      prev = pstart;
	      pcurr = sh6getnext (prev, indstart);
	      sh6getlist (prev, pcurr, &ind1, &ind2, &kstat);
	      if (kstat < 0)
		 goto error;
	      
	      for (ki = 0; ki < knum; ki++)
	      {
		 sh6settop (prev, ind1, l1, r1, l2, r2, &kstat);
		 if (!pcurr)
		    break;
		 sh6settop (pcurr, ind2, l1, r1, l2, r2, &kstat);
		 sh6getother (pcurr, prev, &pnext, &kstat);
		 if (kstat < 0)
		    goto error;
		 prev = pcurr;
		 pcurr = pnext;
		 sh6getlist (prev, pcurr, &ind1, &ind2, &kstat);
		 if (kstat < 0)
		    goto error;
	      }
	   }
	}
	
	else if ((po1->iobj == SISLCURVE && po2->iobj == SISLPOINT) ||
		 (po2->iobj == SISLCURVE && po1->iobj == SISLPOINT))
	   
	{
	   /* Curve point cases */
	   
	   (*pintdat)->vlist[kj]->pretop[0] =
	      (*pintdat)->vlist[kj]->pfirst->left_obj_1[0];
	   (*pintdat)->vlist[kj]->pretop[1] =
	      (*pintdat)->vlist[kj]->plast->right_obj_1[0];
	   (*pintdat)->vlist[kj]->pretop[2] =
	      (*pintdat)->vlist[kj]->pfirst->left_obj_2[0];
	   (*pintdat)->vlist[kj]->pretop[3] =
	      (*pintdat)->vlist[kj]->plast->right_obj_2[0];
	}
	
	else if (po1->iobj == SISLCURVE && po2->iobj == SISLCURVE)
	{
	   
	   /* Curve curve cases */
	   if ((*pintdat)->vlist[kj]->pfirst->epar[0] <=
	       (*pintdat)->vlist[kj]->plast->epar[0])
	   {
	      (*pintdat)->vlist[kj]->pretop[0] =
		 (*pintdat)->vlist[kj]->pfirst->left_obj_1[0];
	      (*pintdat)->vlist[kj]->pretop[1] =
		 (*pintdat)->vlist[kj]->plast->right_obj_1[0];
	   }
	   else
	   {
	      (*pintdat)->vlist[kj]->pretop[0] =
		 (*pintdat)->vlist[kj]->plast->left_obj_1[0];
	      (*pintdat)->vlist[kj]->pretop[1] =
		 (*pintdat)->vlist[kj]->pfirst->right_obj_1[0];
	   }
	   
	   if ((*pintdat)->vlist[kj]->pfirst->epar[1] <=
	       (*pintdat)->vlist[kj]->plast->epar[1])
	   {
	      (*pintdat)->vlist[kj]->pretop[2] =
		 (*pintdat)->vlist[kj]->pfirst->left_obj_2[0];
	      (*pintdat)->vlist[kj]->pretop[3] =
		 (*pintdat)->vlist[kj]->plast->right_obj_2[0];
	   }
	   else
	   {
	      (*pintdat)->vlist[kj]->pretop[2] =
		 (*pintdat)->vlist[kj]->plast->left_obj_2[0];
	      (*pintdat)->vlist[kj]->pretop[3] =
		 (*pintdat)->vlist[kj]->pfirst->right_obj_2[0];
	   }
	}
	
	else if ((po1->iobj == SISLSURFACE && po2->iobj == SISLCURVE
		  && po1->s1->idim == 3))
	{
	   
	   /* Suface curve case */
	   
	   (*pintdat)->vlist[kj]->pretop[0] =
	      (*pintdat)->vlist[kj]->pfirst->left_obj_1[0];
	   (*pintdat)->vlist[kj]->pretop[1] =
	      (*pintdat)->vlist[kj]->plast->right_obj_1[0];
	   
	   if ((*pintdat)->vlist[kj]->pfirst->epar[2] <=
	       (*pintdat)->vlist[kj]->plast->epar[2])
	   {
	      (*pintdat)->vlist[kj]->pretop[2] =
		 (*pintdat)->vlist[kj]->pfirst->left_obj_2[0];
	      (*pintdat)->vlist[kj]->pretop[3] =
		 (*pintdat)->vlist[kj]->plast->right_obj_2[0];
	   }
	   else
	   {
	      (*pintdat)->vlist[kj]->pretop[2] =
		 (*pintdat)->vlist[kj]->plast->left_obj_2[0];
	      (*pintdat)->vlist[kj]->pretop[3] =
		 (*pintdat)->vlist[kj]->pfirst->right_obj_2[0];
	   }
	   
	}
	
	else if ((po1->iobj == SISLCURVE && po2->iobj == SISLSURFACE
		  && po1->c1->idim == 3))
	{
	   
	   /* Curve suface case */
	   
	   if ((*pintdat)->vlist[kj]->pfirst->epar[0] <=
	       (*pintdat)->vlist[kj]->plast->epar[0])
	   {
	      (*pintdat)->vlist[kj]->pretop[0] =
		 (*pintdat)->vlist[kj]->pfirst->left_obj_1[0];
	      (*pintdat)->vlist[kj]->pretop[1] =
		 (*pintdat)->vlist[kj]->plast->right_obj_1[0];
	   }
	   else
	   {
	      (*pintdat)->vlist[kj]->pretop[0] =
		 (*pintdat)->vlist[kj]->plast->left_obj_1[0];
	      (*pintdat)->vlist[kj]->pretop[1] =
		 (*pintdat)->vlist[kj]->pfirst->right_obj_1[0];
	   }
	   
	   (*pintdat)->vlist[kj]->pretop[2] =
	      (*pintdat)->vlist[kj]->pfirst->left_obj_2[0];
	   (*pintdat)->vlist[kj]->pretop[3] =
	      (*pintdat)->vlist[kj]->plast->right_obj_2[0];
	   
	}
     }
     
     
     
     
     *jstat = 0;
     goto out;
     
     /* ------------------------------------------------------ */
     errinconsist:
	*jstat = -500;
     s6err ("sh6idlis", *jstat, kpos);
     goto out;
     
     err101:
	*jstat = -101;
     s6err ("sh6idlis", *jstat, kpos);
     goto out;
     
     error:
	*jstat = kstat;
     s6err ("sh6idlis", *jstat, kpos);
     goto out;
     
     out:
	;
}
   
