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
 * $Id: s6idlis.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */
#define S6IDLIS

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void s6idlis_s9ssexamin(SISLSurf *,SISLSurf *,SISLIntdat **,int *);
static void s6idlis_s9psexamin(SISLSurf *,double,SISLIntdat **,int *);
#else
static void s6idlis_s9ssexamin();
static void s6idlis_s9psexamin();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
     s6idlis(SISLObject *po1,SISLObject *po2,SISLIntdat **pintdat,int *jstat)
#else
void s6idlis(po1,po2,pintdat,jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat **pintdat;
     int    *jstat;
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
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  int kstat;                /* Local status variable.          */
  int kpos=0;               /* Position of error.              */
  int kj,ki1,ki2;           /* Counters                        */
  int ktype;                /* To indicate type of list.       */
  SISLIntpt   *pt;       /* to traverse list of points.     */
  
  *jstat = 0;
  
  /* If we do not have any intersection data we just return. */
  
  if ((*pintdat) == SISL_NULL) goto out;
  
  /* We first destroy existing intersection lists. */
  
  for (kj=0; kj<(*pintdat)->ilist; kj++) freeIntlist((*pintdat)->vlist[kj]);
  
  
  /* Then we split lists with internal junction points. We have to
     be sure that all junction points are end points in the lists. */
  
  for (kj=0; kj<(*pintdat)->ipoint; kj++)
    if ((*pintdat)->vpoint[kj]->iinter == 2)
      {
	if ((*pintdat)->vpoint[kj]->pcurve != SISL_NULL)
	  {
	    for (ki1=0; ki1<(*pintdat)->ipoint; ki1++)
	      if ((*pintdat)->vpoint[ki1]->pcurve == (*pintdat)->vpoint[kj])
		break;
	    
	    if (ki1<(*pintdat)->ipoint)
	      {
		pt = copyIntpt((*pintdat)->vpoint[kj]);
		
		s6idnpt(pintdat,&pt,0,&kstat);
		if (kstat < 0) goto error;
		
		pt->pcurve = (*pintdat)->vpoint[kj]->pcurve;
		
		(*pintdat)->vpoint[kj]->pcurve = SISL_NULL;
	      }
	  }
      }
  
  
  /* At least we can traverse all intersection points to look for
     start points to lists. If a point have a next point
     and no other point pointing on itself. It is a start point. */
  
  for (ki1=0,ki2=0; ki1 < (*pintdat)->ipoint; ki1++)
    if ((*pintdat)->vpoint[ki1]->pcurve != SISL_NULL)
      {
	for (kj=0; kj<(*pintdat)->ipoint; kj++)
	  if ((*pintdat)->vpoint[kj]->pcurve == (*pintdat)->vpoint[ki1])
	    break;
	
	if (kj == (*pintdat)->ipoint)
	  {
	    /* To be sure that list array is big enough. */
	    
	    if (ki2 == (*pintdat)->ilmax)
	      {
		(*pintdat)->ilmax += 20;
		
		if (((*pintdat)->vlist = increasearray((*pintdat)->vlist,
						       (*pintdat)->ilmax,SISLIntlist *)) == SISL_NULL)
		  
		  goto err101;
	      }
	    
	    
	    /* Finding the last point in the list, and number of points. */
	    
	    kj = 0;
	    for (pt=(*pintdat)->vpoint[ki1];pt->pcurve!=SISL_NULL;
		 pt=pt->pcurve,kj++);
	    
	    
	    /* Computing type of point, junctions in the end points. */
	    
	    ktype = 0;
	    
	    if ((*pintdat)->vpoint[ki1]->iinter == 2)
	      ktype = 2;
	    
	    if (pt->iinter == 2)
	      ktype = (ktype == 2 ? 4 : 3);
	    
	    
	    /* Making a new list structure. */
	    
	    if (((*pintdat)->vlist[ki2] = newIntlist((*pintdat)->vpoint[ki1],
						     pt,ktype)) == SISL_NULL) goto err101;
	    
	    (*pintdat)->vlist[ki2]->inumb = kj + 1;
	    ki2++;
	    
	  }
      }
  
  /*------------------------------------------------------------------*/
  
  /* We also have to find closed lists.    */
  
  /* Mark found list elements: */
  for (ki1=0; ki1 < ki2; ki1++)
    for (pt=(*pintdat)->vlist[ki1]->pfirst;pt!=SISL_NULL;pt=pt->pcurve)
      pt->iinter += 10;
  
  /* Now travers the point array untill we find an unmarked point.
     This point has to be a single (unconnected) one or a member 
     of a closed connection. Mark points in the closed connection and 
     establish a new list. */
  
  for (ki1=0; ki1 < (*pintdat)->ipoint; ki1++)
    {
      if ((*pintdat)->vpoint[ki1]->iinter>=10)
	/* Unmark point. */
	(*pintdat)->vpoint[ki1]->iinter -= 10;
      
      else  if ((*pintdat)->vpoint[ki1]->pcurve != SISL_NULL)
	{
	  /* It has to be a closed connection, travers all elements. */
	  kj = 1;
	  for (pt=(*pintdat)->vpoint[ki1]->pcurve;pt!=(*pintdat)->vpoint[ki1];
	       pt=pt->pcurve)
	    {	
	      if (pt == SISL_NULL) goto err105;
	      /* Mark found list elements: */
	      pt->iinter += 10;
	      kj++;
	    }
	  
	  /*Create new list element. */
	  
	  /* To be sure that list array is big enough. */
	  if (ki2 == (*pintdat)->ilmax)
	    {
	      (*pintdat)->ilmax += 20;
	      
	      if (((*pintdat)->vlist = increasearray((*pintdat)->vlist,
						     (*pintdat)->ilmax,SISLIntlist *)) == SISL_NULL) 
		goto err101;
	    }
	  
	  /* Closed curves will have no singularities: */
	  ktype = 1;
	  
	  /* Making a new list structure. */
	  if (((*pintdat)->vlist[ki2] = 
	       newIntlist((*pintdat)->vpoint[ki1]->pcurve,
			  (*pintdat)->vpoint[ki1],ktype)) == SISL_NULL) 
	    goto err101;
	  (*pintdat)->vlist[ki2]->inumb = kj;
	  ki2++;
	  
	}
    }
  
  (*pintdat)->ilist = ki2;
  
  /*------------------------------------------------------------------*/
  
  /* A final check if the geometry found is ok. */
  
  if (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE && po1->s1->idim == 3)
    {
       s6idlis_s9ssexamin(po1->s1,po2->s1,pintdat,&kstat);
      if (kstat < 0) goto error;
    }
  else if (po1->iobj == SISLPOINT && po2->iobj == SISLSURFACE && po1->p1->idim == 1)
    {
       s6idlis_s9psexamin(po2->s1,po1->p1->ecoef[0],pintdat,&kstat);
      if (kstat < 0) goto error;
    }
  else if (po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT && po2->p1->idim == 1)
    {
       s6idlis_s9psexamin(po1->s1,po2->p1->ecoef[0],pintdat,&kstat);
      if (kstat < 0) goto error;
    }
  
  goto out;
  
  /* Error in space allocation.  */
  
  err101: *jstat = -101;
  s6err("s6idlis",*jstat,kpos);
  goto out;
  
  /* Error in vpoint array.  */
  
  err105: *jstat = -105;
  s6err("s6idlis",*jstat,kpos);
  goto out;
  
  /* Error in sub function.  */
  
  error:  *jstat = kstat;
  s6err("s6idlis",*jstat,kpos);
  goto out;
  
  out: ;
  
}

#if defined(SISLNEEDPROTOTYPES)
   static
      void
            s6idlis_s9ssexamin(SISLSurf *ps1,SISLSurf *ps2,
			       SISLIntdat **rintdat,int *jstat)
#else
static void s6idlis_s9ssexamin(ps1,ps2,rintdat,jstat)
     SISLSurf   *ps1;
     SISLSurf   *ps2;
     SISLIntdat **rintdat;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To examin if the intersection beetween two surfaces
*              is correct. If something wrong is found the error is
*              is removed.
*
*
* INPUT      : ps1      - First surface in intersection.
*              ps2      - Second surface in intersection.
*
*
* OUTPUT     : rintdat  - Intersection dates to be examin.
*              jstat    - status messages  
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
*              Ulf J. Krystad
*********************************************************************
*/                                     
{
  int kstat;
  int kdirstat;
  
  unsigned char edg=0;
  SISLIntpt **uipt=SISL_NULL;
  SISLIntlist **uilst=SISL_NULL;
  
  int ki,kj,kv,klnr,kpnr,klfs,klft,kdir,kpar;
  SISLIntpt  *qpt1,*qpt2;
  SISLIntpt  *qipt, *qp;
  
  double tepsge, tang,sedg[8];
  double sval1[9],sval2[9],snorm1[3],snorm2[3];
  double stang[3],sdec1[3],sdec2[3];
  double epar1[4],epar2[4];
  
  tepsge = (double)0.001;
  
  sedg[0] = ps1->et1[ps1->ik1-1];
  sedg[1] = ps1->et1[ps1->in1];
  sedg[2] = ps1->et2[ps1->ik2-1];
  sedg[3] = ps1->et2[ps1->in2];
  sedg[4] = ps2->et1[ps2->ik1-1];
  sedg[5] = ps2->et1[ps2->in1];
  sedg[6] = ps2->et2[ps2->ik2-1];
  sedg[7] = ps2->et2[ps2->in2];
  
  /* Init */  
  if (!(*rintdat)) goto out;
  if (ps1->idim != 3 || ps2->idim != 3) goto err200;
  *jstat = 0;
  
  
  /* SISLCurve analyse section --------------------------------------------------*/
  if ((*rintdat)->ilist != 0)
    {
      /* Allocate array of pointers to the lists. */
      klnr = (*rintdat)->ilist;
      if ((uilst = newarray(klnr,SISLIntlist *)) == SISL_NULL) goto err101;
      
      /* Update the list array. */
      
      /* Get all open curves into array. */
      for (kv=ki=0; ki<klnr; ki++)
        if ((*rintdat)->vlist[ki]->itype != 1)
	  uilst[kv++] = (*rintdat)->vlist[ki];
      
      /* Correct number of curves.*/
      klnr = kv;
      
      /* Remove all open curves with endpoints on edges or singular 
	 endpoints from array.(they are ok.) */
      for (ki=0; ki<klnr; ki++)
	{
	  for (kj=0; kj<8; kj++)
	    if (DEQUAL(uilst[ki]->pfirst->epar[(kj/2)],sedg[kj]))
	      break;
	  
	  
	  if (kj == 8)
	    {
	      /* Start point is NOT on edge, test if the point 
		 is a singular point. */
	      
	      klfs=klft=0;
	      s1421(ps1,1,uilst[ki]->pfirst->epar,&klfs,&klft,sval1,
		    snorm1,&kstat);
	      if (kstat < 0) goto error;
	      else if (kstat > 0 ) kj--;
	      
	      klfs=klft=0;
	      s1421(ps2,1,uilst[ki]->pfirst->epar+2,&klfs,&klft,sval2,
		    snorm2,&kstat);
	      if (kstat < 0) goto error;
	      else if (kstat > 0 ) kj--;
	      
	      tang = s6ang(snorm1,snorm2,3);
	      if (tang < ANGULAR_TOLERANCE) kj--;
	    }	    
	  
	  
	  if (kj<8)
	    {
	      /* Start point is ok, test end point*/
	      for (kj=0; kj<8; kj++)
		if (DEQUAL(uilst[ki]->plast->epar[(kj/2)],sedg[kj]))
		  break;
	      
	      
	      if (kj == 8)
		{
		  /* End point is NOT on edge, test if the point 
		     is a singular point. */
		  klfs=klft=0;
		  s1421(ps1,1,uilst[ki]->plast->epar,&klfs,&klft,sval1,
			snorm1,&kstat);
		  if (kstat < 0) goto error;
		  else if (kstat > 0 ) kj--;
		  
		  klfs=klft=0;
		  s1421(ps2,1,uilst[ki]->plast->epar+2,&klfs,&klft,sval2,
			snorm2,&kstat);
		  if (kstat < 0) goto error;
		  else if (kstat > 0 ) kj--;
		  
		  tang = s6ang(snorm1,snorm2,3);
		  if (tang < ANGULAR_TOLERANCE) kj--;
		}	    
	      
	      
	      if (kj<8)
		{
		  /* Start point and end point is ok, remove it from the array*/
		  klnr--;
		  if (ki<klnr)
		    {
		      uilst[ki] = uilst[klnr];
		      ki--;
		    }
		  
		}
	    }
	}
      
      /* Now we only have curves with bad endpoints in the array. */
      
      for (ki=0; ki< klnr; ki++)
	{
	  /* Now we kill all the points in the list except the
	     end point that is an internal point . */
	  
	  for (kj=0; kj<8; kj++)
	    if (DEQUAL(uilst[ki]->pfirst->epar[(kj/2)],sedg[kj]))
	      break;
	  
	  if (kj<8)
	    {
	      
	      /* The first point is on the edge, keep the last. */
	      qipt = uilst[ki]->pfirst;
	      for (qp=qipt->pcurve; qipt != uilst[ki]->plast;qipt=qp,qp=qp->pcurve)
		{
		  s6idkpt(rintdat,&qipt,&qpt1,&qpt2,&kstat);
		  if (kstat < 0) goto error;
		}
	      
	      qipt = uilst[ki]->plast;
	      uilst[ki]->pfirst = uilst[ki]->plast =  qipt;
	      uilst[ki]->inumb = 1;
	    }
	  else
	    {
	      /* The first point is not on the edge, keep it. */
	      
	      for (qipt = uilst[ki]->pfirst->pcurve; qipt != SISL_NULL;qipt=qp)	      
		{
		  s6idkpt(rintdat,&qipt,&qpt1,&qp,&kstat);
		  if (kstat < 0) goto error;
		}
	      
	      qipt = uilst[ki]->pfirst;
	      uilst[ki]->pfirst = uilst[ki]->plast =  qipt;
	      uilst[ki]->inumb = 1;
	      
	    }
	  
	  /* March from point qipt */
	  s1788(ps1,ps2,tepsge,qipt->epar,epar1,epar2,&kstat);
          if (kstat<0) goto error;
	  
	  if (kstat == 0)
            {
	      /* No succes. */
	      /* Kill point and the list */
	      s6idklist(rintdat,uilst[ki],&kstat);
	      if (kstat<0) goto error;
	      klnr--;
	      if (ki < klnr)
		{
		  uilst[ki] = uilst[klnr];
		  ki--;
		}
	      
	    }
	  else if (kstat == 11 || kstat == 12 || kstat == 13 ||
	           kstat == 14 || kstat == 21 || kstat == 22 || kstat == 24 )
	    {
	      /* Making a new open curve with endpoint in epar1 and epar2.*/
	      
	      uilst[ki]->pfirst = newIntpt(4,epar1,DZERO);
	      if (uilst[ki]->pfirst == SISL_NULL) goto err101;
	      
	      s6idnpt(rintdat,&uilst[ki]->pfirst,0,&kstat);
	      if (kstat < 0)goto error;
	      
	      uilst[ki]->plast = newIntpt(4,epar2,DZERO);
	      if (uilst[ki]->plast == SISL_NULL) goto err101;
	      
	      s6idnpt(rintdat,&uilst[ki]->plast,0,&kstat);
	      if (kstat < 0)goto error;
	      
	      uilst[ki]->pfirst->pcurve = qipt;
	      qipt->pcurve = uilst[ki]->plast;
	      uilst[ki]->inumb = 3;
              uilst[ki]->itype = 4;
	    }
	  
	  else if (kstat == 16 || kstat == 17 || kstat == 26 || kstat == 27)
	    {
	      /* Making a new closed curve with pfirst and plast
                 pointing on qipt.*/
	      
              uilst[ki]->pfirst = uilst[ki]->plast = qipt;
	      uilst[ki]->inumb = 1;
              uilst[ki]->itype = 1;
	      qipt->pcurve = qipt;
	      
	    }
	}
      /* Now we check equality between the remaining curves in the array. */
      for (ki=0;ki<klnr-1;ki++)
	for (kj=ki+1;kj<klnr;kj++)
	  {
	    if ((s6dist(uilst[ki]->pfirst->epar,uilst[kj]->pfirst->epar,4)
		 < REL_COMP_RES &&
		 s6dist(uilst[ki]->plast->epar,uilst[kj]->plast->epar,4)
		 < REL_COMP_RES)  ||
		(s6dist(uilst[ki]->pfirst->epar,uilst[kj]->plast->epar,4)
		 < REL_COMP_RES &&
		 s6dist(uilst[ki]->plast->epar,uilst[kj]->pfirst->epar,4)
		 < REL_COMP_RES))
	      /* The two curves has common start+end, remove the last one of them. */
	      {
		
		s6idklist(rintdat,uilst[kj],&kstat);
		if (kstat<0) goto error;
		
		klnr--;
		if (kj < klnr)
		  {
		    uilst[kj] = uilst[klnr];
		    kj--;
		  }
		
	      }
	  }
      
    }
  
  
  /* End of curve analyse section -------------------------------------------*/
  
  /* SISLPoint analyse section --------------------------------------------------*/
  
  if ((*rintdat) && (*rintdat)->ipoint != 0)
    /* Update the point array. */
    {     
      kpnr = (*rintdat)->ipoint;
      if ((uipt = newarray(kpnr,SISLIntpt *)) == SISL_NULL) goto err101;
      
      
      for (kv=ki=0; ki<kpnr; ki++)
        if ((*rintdat)->vpoint[ki]->pcurve == SISL_NULL)
	  uipt[kv++] = (*rintdat)->vpoint[ki];
      
      for (ki=0; ki<kpnr; ki++)
        for (kj=0; kj<kv; kj++)
          if ((*rintdat)->vpoint[ki]->pcurve == uipt[kj])
            {
	      kv--;
      	      uipt[kj] = uipt[kv];
              break;
            }
      
      /* All single points found. */
      kpnr = kv;
      
      
      /* Sorting out and killing all points but single touch points. */
      
      for (ki=0; ki<kpnr; ki++)
	{
	  klfs=klft=0;
	  s1421(ps1,1,uipt[ki]->epar,&klfs,&klft,sval1,snorm1,&kstat);
	  if (kstat < 0) goto error;
	  else if (kstat > 0 ) continue;
	  
	  if (s6length(snorm1,3,&kstat) <= REL_COMP_RES) continue;
	  
	  klfs=klft=0;
	  s1421(ps2,1,uipt[ki]->epar+2,&klfs,&klft,sval2,snorm2,&kstat);
	  if (kstat < 0) goto error;
	  else if (kstat > 0 ) continue;
	  
	  if (s6length(snorm2,3,&kstat) <= REL_COMP_RES) continue;
	  
	  tang = s6ang(snorm1,snorm2,3);
	  if (tang < ANGULAR_TOLERANCE) continue;	    
	  
	  
	  /* All singular points or degenerated points is ok. We
	     then remove all other internal points. */
	  
	  for (kj=0,edg=0; kj<8; kj++)
	    if (DEQUAL(uipt[ki]->epar[(kj/2)],sedg[kj]))
	      edg |= 1<<kj;
	  
	  if (edg == 0)
	    {
	      /* The point is removed. */
	      
	      s6idkpt(rintdat,&uipt[ki],&qpt1,&qpt2,&kstat);
	      if (kstat < 0) goto error;
	      
	      kpnr--;
	      
	      if (ki < kpnr)
		{
		  uipt[ki] = uipt[kpnr];
		  ki--;
		}
	      
	      continue;
	    }
	  
	  
	  s6crss(snorm1,snorm2,stang);
	  
	  s6decomp(stang,sdec1,sval1+3,sval1+6,snorm1,&kstat);
	  if (kstat < 0) goto error;
	  else if (kstat > 0 ) continue;
	  
	  s6decomp(stang,sdec2,sval2+3,sval2+6,snorm2,&kstat);
	  if (kstat < 0) goto error;
	  else if (kstat > 0 ) continue;
	  
	  for (kpar=1,kdir=kdirstat=0,kj=0; kj<8; kj++)
	    if ((edg & 1<<kj) == 1<<kj)
	      {
		switch (kj) 
		  {
		  case 0: tang = s6ang(stang,sval1+3,3);
		    kdir = (sdec1[1] > DZERO ?  1 : -1);
		    break;
		  case 4: tang = s6ang(stang,sval2+3,3);
		    kdir = (sdec2[1] > DZERO ?  1 : -1);
		    break;
		  case 1: tang = s6ang(stang,sval1+6,3);
		    kdir = (sdec1[0] > DZERO ?  -1 : 1);
		    break;
		  case 5: tang = s6ang(stang,sval2+6,3);
		    kdir = (sdec2[0] > DZERO ?  -1 : 1);
		    break;
		  case 2: tang = s6ang(stang,sval1+3,3);
		    kdir = (sdec1[1] > DZERO ?  -1 : 1);
		    break;
		  case 6: tang = s6ang(stang,sval2+3,3);
		    kdir = (sdec2[1] > DZERO ?  -1 : 1);
		    break;
		  case 3: tang = s6ang(stang,sval1+6,3);
		    kdir = (sdec1[0] > DZERO ?  1 : -1);
		    break;
		  case 7: tang = s6ang(stang,sval2+6,3);
		    kdir = (sdec2[0] > DZERO ?  1 : -1);
		  }
		
		if (tang < ANGULAR_TOLERANCE) kdir = 0;
		
		if (kdir == 0)
		  kpar = 0;
		else if (kdirstat != kdir)
		  {
		    if (kdirstat == 0)
		      kdirstat = kdir;
		    else
		      {
			kdirstat = 10;
			break;
		      }
		  }
	      }
	  
	  if (kpar == 0 && kdirstat != 10) kdirstat = 0;
	  
	  /* Test if the point is to be removed.*/
	  if (kdirstat == 1 || kdirstat == -1)
	    {
	      /* The point is removed. */
	      
	      s6idkpt(rintdat,&uipt[ki],&qpt1,&qpt2,&kstat);
	      if (kstat < 0) goto error;
	      
	      kpnr--;
	      
	      if (ki < kpnr)
		{
		  uipt[ki] = uipt[kpnr];
		  ki--;
		}
	      
	      continue;
	    }
	  
	} /* End of for ki= .... */
      
    }
  /* End of point analyse section ---------------------------------------------*/ 
  
  
  goto out;
  
  /* Error in sub rutines.      */
  
  error : *jstat = kstat;
  s6err("s6idlis_s9ssexamin",*jstat,0);
  goto out;                       
  
  /* Error in memory allocation.      */
  
  err101 : *jstat = -101;
  s6err("s6idlis_s9ssexamin",*jstat,0);
  goto out;                       
  
  /* Error dimention.      */
  
  err200 : *jstat = -200;
  s6err("s6idlis_s9ssexamin",*jstat,0);
  goto out;                       
  
  out:
  if (uipt != SISL_NULL)  freearray(uipt);
  if (uilst != SISL_NULL)   freearray(uilst);
}

#if defined(SISLNEEDPROTOTYPES)
   static
      void
            s6idlis_s9psexamin(SISLSurf *ps1,double alevel,
			       SISLIntdat **rintdat,int *jstat)
#else
static void s6idlis_s9psexamin(ps1,alevel,rintdat,jstat)
     SISLSurf   *ps1;
     double alevel;
     SISLIntdat **rintdat;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To examin if the intersection beetween a surface
*              and a point in dim 1 is correct. If something wrong 
*              is found the error is
*              is removed.
*
*
* INPUT      : ps1      - First surface in intersection.
*              alevel   - The point value.
*
*
*INPUT/OUTPUT: rintdat  - Intersection dates to be examin.
*
* OUTPUT     : jstat    - status messages  
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
* WRITTEN BY : Ulf J. Krystad, SI, 10.08.89
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              jstat set to 0 in start.
*
*********************************************************************
*/                                     
{
  int kstat;
  int kdirstat;
  
  unsigned char edg=0;
  SISLIntpt **uipt=SISL_NULL;
  SISLIntlist **uilst=SISL_NULL;
  
  int ki,kj,kv,klnr,kpnr,klfs,klft,kdir,kpar;
  SISLIntpt *qpt1,*qpt2;
  SISLIntpt *qipt, *qp;
  
  double tepsge, tmax,sedg[4];
  double sval1[9],snorm1[3];
  double epar1[2],epar2[2];
/* ALA && UJK 19.09.90 */
  double ttol= 10000.0 * REL_COMP_RES;
  
  tepsge = (double)0.001;
  
  sedg[0] = ps1->et1[ps1->ik1-1];
  sedg[1] = ps1->et1[ps1->in1];
  sedg[2] = ps1->et2[ps1->ik2-1];
  sedg[3] = ps1->et2[ps1->in2];
  
  
  /* Init */ 
  if (!(*rintdat)) goto out;
  if (ps1->idim != 1) goto err200;
  *jstat = 0;
  
  
  /* SISLCurve analyse section --------------------------------------------------*/
  
  if ((*rintdat)->ilist != 0)
    {
      /* Allocate array of pointers to the lists. */
      klnr = (*rintdat)->ilist;
      if ((uilst = newarray(klnr,SISLIntlist *)) == SISL_NULL) goto err101;
      
      /* Update the list array. */
      
      /* Get all open curves into array. */
      for (kv=ki=0; ki<klnr; ki++)
        if ((*rintdat)->vlist[ki]->itype != 1)
	  uilst[kv++] = (*rintdat)->vlist[ki];
      
      /* Correct number of curves.*/
      klnr = kv;
      
      /* Remove all open curves with endpoints on edges or singular 
	 endpoints from array.(they are ok.) */
      for (ki=0; ki<klnr; ki++)
	{
	  for (kj=0; kj<4; kj++)
	    if (DEQUAL(uilst[ki]->pfirst->epar[(kj/2)],sedg[kj]))
	      break;
	  
	  
	  if (kj == 4)
	    {
	      /* Start point is NOT on edge, test if the point 
		 is a singular point. */
	      
	      klfs=klft=0;
	      s1421(ps1,1,uilst[ki]->pfirst->epar,&klfs,&klft,sval1,
		    snorm1,&kstat);
	      if (kstat < 0) goto error;
	      else if (kstat > 0 ) kj--;
	      
	      tmax = sqrt(sval1[1]*sval1[1] + sval1[2]*sval1[2]);
	      if ( tmax < ttol ) kj--;
	      
	    }	    
	  
	  
	  if (kj<4)
	    {
	      /* Start point is ok, test end point*/
	      for (kj=0; kj<4; kj++)
		if (DEQUAL(uilst[ki]->plast->epar[(kj/2)],sedg[kj]))
		  break;
	      
	      
	      if (kj == 4)
		{
		  /* End point is NOT on edge, test if the point 
		     is a singular point. */
		  klfs=klft=0;
		  s1421(ps1,1,uilst[ki]->plast->epar,&klfs,&klft,sval1,
			snorm1,&kstat);
		  if (kstat < 0) goto error;
		  else if (kstat > 0 ) kj--;
		  
		  tmax = sqrt(sval1[1]*sval1[1] + sval1[2]*sval1[2]);
		  if ( tmax < ttol ) kj--;
		  
		}	    
	      
	      
	      if (kj<4)
		{
		  /* Start point and end point is ok, remove it from the array*/
		  klnr--;
		  if (ki<klnr)
		    {
		      uilst[ki] = uilst[klnr];
		      ki--;
		    }
		  
		}
	    }
	}
      
      /* Now we only have curves with bad endpoints in the array. */
      
      for (ki=0; ki< klnr; ki++)
	{
	  /* Now we kill all the points in the list except the
	     end point that is an internal point . */
	  
	  for (kj=0; kj<4; kj++)
	    if (DEQUAL(uilst[ki]->pfirst->epar[(kj/2)],sedg[kj]))
	      break;
	  
	  if (kj<4)
	    {
	      
	      /* The first point is on the edge, keep the last. */
	      qipt = uilst[ki]->pfirst;
	      for (qp=qipt->pcurve; qipt != uilst[ki]->plast;qipt=qp,qp=qp->pcurve)
		{
		  s6idkpt(rintdat,&qipt,&qpt1,&qpt2,&kstat);
		  if (kstat < 0) goto error;
		}
	      
	      qipt = uilst[ki]->plast;
	      uilst[ki]->pfirst = uilst[ki]->plast =  qipt;
	      uilst[ki]->inumb = 1;
	    }
	  else
	    {
	      /* The first point is not on the edge, keep it. */
	      
	      for (qipt = uilst[ki]->pfirst->pcurve; qipt != SISL_NULL;qipt=qp)	      
		{
		  s6idkpt(rintdat,&qipt,&qpt1,&qp,&kstat);
		  if (kstat < 0) goto error;
		}
	      qipt = uilst[ki]->pfirst;
	      uilst[ki]->pfirst = uilst[ki]->plast =  qipt;
	      uilst[ki]->inumb = 1;
	    }
	  
	  /* March from point qipt */
	  s1787(ps1,alevel,tepsge,qipt->epar,epar1,epar2,&kstat);
          if (kstat<0) goto error;
	  
	  if (kstat == 0)
            {
	      /* No succes. */
	      /* Kill point and the list */
	      s6idklist(rintdat,uilst[ki],&kstat);
	      if (kstat<0) goto error;
	      
	      klnr--;
	      if (ki < klnr)
		{
		  uilst[ki] = uilst[klnr];
		  ki--;
		}
	      
	    }
	  else if (kstat == 11 || kstat == 12 || kstat == 13 ||
	           kstat == 14 || kstat == 21 || kstat == 22 || kstat == 24 )
	    {
	      /* Making a new open curve with endpoint in epar1 and epar2.*/
	      
	      uilst[ki]->pfirst = newIntpt(2,epar1,DZERO);
	      if (uilst[ki]->pfirst == SISL_NULL) goto err101;
	      
	      s6idnpt(rintdat,&uilst[ki]->pfirst,0,&kstat);
	      if (kstat < 0)goto error;
	      
	      uilst[ki]->plast = newIntpt(2,epar2,DZERO);
	      if (uilst[ki]->plast == SISL_NULL) goto err101;
	      
	      s6idnpt(rintdat,&uilst[ki]->plast,0,&kstat);
	      if (kstat < 0)goto error;
	      
	      uilst[ki]->pfirst->pcurve = qipt;
	      qipt->pcurve = uilst[ki]->plast;
	      uilst[ki]->inumb = 3;
              uilst[ki]->itype = 4;
	    }
	  
	  else if (kstat == 16 || kstat == 17 || kstat == 26 || kstat == 27)
	    {
	      /* Making a new closed curve with pfirst and plast
                 pointing on qipt.*/
	      
              uilst[ki]->pfirst = uilst[ki]->plast = qipt;
	      uilst[ki]->inumb = 1;
              uilst[ki]->itype = 1;
	      qipt->pcurve = qipt;
	    }
	}
      /* Now we check equality between the remaining curves in the array. */
      for (ki=0;ki<klnr-1;ki++)
	for (kj=ki+1;kj<klnr;kj++)
	  {
	    if ((s6dist(uilst[ki]->pfirst->epar,uilst[kj]->pfirst->epar,2)
		 < REL_COMP_RES &&
		 s6dist(uilst[ki]->plast->epar,uilst[kj]->plast->epar,2)
		 < REL_COMP_RES)  ||
		(s6dist(uilst[ki]->pfirst->epar,uilst[kj]->plast->epar,2)
		 < REL_COMP_RES &&
		 s6dist(uilst[ki]->plast->epar,uilst[kj]->pfirst->epar,2)
		 < REL_COMP_RES))
	      /* The two curves has common start+end, remove the last one of them. */
	      {
		
		s6idklist(rintdat,uilst[kj],&kstat);
		if (kstat<0) goto error;
		
		klnr--;
		if (kj < klnr)
		  {
		    uilst[kj] = uilst[klnr];
		    kj--;
		  }
		
	      }
	  }
      
    }
  
  
  /* End of curve analyse section -------------------------------------------*/
  
  /* SISLPoint analyse section --------------------------------------------------*/
  
  if ((*rintdat) && (*rintdat)->ipoint != 0)
    
    /* Update the point array. */
    {     
      kpnr = (*rintdat)->ipoint;
      if ((uipt = newarray(kpnr,SISLIntpt *)) == SISL_NULL) goto err101;
      
      
      for (kv=ki=0; ki<kpnr; ki++)
        if ((*rintdat)->vpoint[ki]->pcurve == SISL_NULL)
	  uipt[kv++] = (*rintdat)->vpoint[ki];
      
      for (ki=0; ki<kpnr; ki++)
        for (kj=0; kj<kv; kj++)
          if ((*rintdat)->vpoint[ki]->pcurve  == uipt[kj])
            {
	      kv--;
      	      uipt[kj] = uipt[kv];
              break;
            }
      
      /* All single points found. */
      kpnr = kv;
      
      
      /* Sorting out and killing all points but single touch points. */
      
      for (ki=0; ki<kpnr; ki++)
	{
	  klfs=klft=0;
	  s1421(ps1,1,uipt[ki]->epar,&klfs,&klft,sval1,snorm1,&kstat);
	  if (kstat < 0) goto error;
	  else if (kstat > 0 ) continue;
	  
	  tmax = sqrt(sval1[1]*sval1[1] + sval1[2]*sval1[2]);
	  if ( tmax < ttol ) continue;
	  
	  
	  /* All singular points or degenerated points is ok. We
	     then remove all other internal points. */
	  
	  for (kj=0,edg=0; kj<4; kj++)
	    if (DEQUAL(uipt[ki]->epar[(kj/2)],sedg[kj]))
	      edg |= 1<<kj;
	  
	  if (edg == 0)
	    {
	      /* The point is removed. */
	      
	      s6idkpt(rintdat,&uipt[ki],&qpt1,&qpt2,&kstat);
	      if (kstat < 0) goto error;
	      
	      kpnr--;
	      
	      if (ki < kpnr)
		{
		  uipt[ki] = uipt[kpnr];
		  ki--;
		}
	      
	      continue;
	    }
	  
	  /* Now we remove all edge points with in/out component. */
	  for (kpar=1,kj=0,kdir=kdirstat=0; kj<4; kj++)
	    if ((edg & 1<<kj) == 1<<kj)
	      {
		switch (kj) 
		  {
		  case 0:
		    if (fabs(sval1[1]/tmax) < ttol)
		      kdir = 0;
		    else
		      kdir = (sval1[1] > DZERO ?  1 : -1);
		    break;
		  case 1:
		    if (fabs(sval1[2]/tmax) < ttol)
		      kdir = 0;
		    else
		      kdir = (sval1[2] > DZERO ?  1 : -1);
		    break;
		  case 2:
		    if (fabs(sval1[1]/tmax) < ttol)
		      kdir = 0;
		    else
		      kdir = (sval1[1] > DZERO ?  -1 : 1);
		    break;
		  case 3:
		    if (fabs(sval1[2]/tmax) < ttol)
		      kdir = 0;
		    else
		      kdir = (sval1[2] > DZERO ?  -1 : 1);
		  }
		
		if (kdir == 0)
		  kpar = 0;
		else if (kdirstat != kdir)
		  {
		    if (kdirstat == 0)
		      kdirstat = kdir;
		    else
		      {
			kdirstat = 10;
			break;
		      }
		  }
	      }
	  
	  if (kpar == 0 && kdirstat != 10) kdirstat = 0;
	  
	  /* Test if the point is to be removed.*/
	  if (kdirstat == 1 || kdirstat == -1)
	    {
	      /* The point is removed. */
	      
	      s6idkpt(rintdat,&uipt[ki],&qpt1,&qpt2,&kstat);
	      if (kstat < 0) goto error;
	      
	      kpnr--;
	      
	      if (ki < kpnr)
		{
		  uipt[ki] = uipt[kpnr];
		  ki--;
		}
	      
	      continue;
	    }
	  
	} /* End of for ki= .... */
      
    }
  /* End of point analyse section ---------------------------------------------*/ 
    
  goto out;
  
  /* Error in sub rutines.      */
  
  error : 
    *jstat = kstat;
    s6err("s6idlis_s9psexamin",*jstat,0);
    goto out;                       
  
  /* Error in memory allocation.      */
  
  err101 : 
    *jstat = -101;
    s6err("s6idlis_s9psexamin",*jstat,0);
    goto out;                       
  
  /* Error dimention.      */
  
  err200 : 
    *jstat = -200;
    s6err("s6idlis_s9psexamin",*jstat,0);
    goto out;                       
  
  out:
    if (uipt != SISL_NULL)  freearray(uipt);
    if (uilst != SISL_NULL) freearray(uilst);
}
