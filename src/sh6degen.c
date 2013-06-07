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
 * $Id: sh6degen.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6DEGEN

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static void sh6degen_geom(SISLObject *,SISLObject *,double [],double [],
			  int *);
#else
static void sh6degen_geom();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  sh6degen(SISLObject * po1, SISLObject * po2, SISLIntdat ** pintdat,
	   double aepsge, int *jstat)
#else
void sh6degen(po1, po2, pintdat, aepsge, jstat)
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
* PURPOSE    : To replace any curves which are degenerate (to a
*              point) by points and to remove any subsequent
*              duplicate points.
*
*
* INPUT      : pintdat  - Pointer to a pointer to intersection data.
*
*
* OUTPUT     : jstat  - status messages
*                               = 0      : OK
*                               < 0      : error
*
*
* METHOD     : Run through the Intlists. If an Intlist is degenerate,
*              then it is a point. If this same point occurs elsewhere
*              in pintdat (as an Intpt or an end point of another
*              Intlist) then delete it. Otherwise add it to the
*              array of Intpoints.
*              At the end, the Intpoint array may have increased
*              and the Intlist array may have decresed.
*
*
* REFERENCES :
*
*-
* CALLS      : sh6getnext       - Get next point in Intlist.
*              sh6getother      - Get next point in Intlist.
*              sh6idkpt         - Remove intersection point.
*              s6length         - Lenght of vector.
*              s6dist           - Distance between points.
*              sh6degen_geom    - Fetch geometry.
*              newIntpt    - Create new Intpt structure.
*              freeIntlist - Free an intlist structure.
*
* WRITTEN BY :Mike Floater, 20/2/92.
* REWISED BY: Vibeke Skytt, SI, 06.92.
*
*********************************************************************
*/
{
  int kstat;		/* Local status variable.          */
  int kpos = 0;		/* Position of error.              */
  int ki,kj, kl;        /* Loop variables.        */
  int kpoint;           /* Number of points in list.       */
  int kstoremark;       /* Save marker of point.   */
  int knumbpt;          /* Number of points in intersection list. */
  int knpar;            /* Number of parameter directions.        */
  int kexact;           /* Indicates which parameter direction is exact. */
  SISLIntpt** vpoint;   /* Array of Int points.    */
  int ipoint;           /* Number of Intpoints.            */
  int ipmax;            /* Size of vpoint array.           */
  SISLIntlist** vlist;  /* Array of Intlists.   */
  int ilist;            /* Number of Intlists.             */
  int ilmax;            /* Size of vlist array.           */
  int isdegen;          /* Flag. Is curve degenerate? */
  double spar[4];       /* Parameter value on intersection curve of
			   the  objects.                         */
  SISLIntpt *pprev, *pcurr, *pnext;    /*    Pointers to Intpts. */    
  SISLIntpt *pfirst, *plast;           /*    Pointers to Intpts. */    
  SISLIntpt *qpt;  
  SISLObject *qo2=SISL_NULL; /* Either po2 or dummy point. */
  int newilist;         /* New number of Intlists. */
  int idim;              /* Dimension of space. */
  double geomfirst[3]; /* geometric info --- position etc. */
  double geomcurr[3];  /* geometric info --- position etc. */
  double geommid[3];   /* geometric info --- position etc. */
  double len;    /* Distance between two vectors. */
  
  knpar = po1->iobj + po2->iobj;
  *jstat = 0;

  if(*pintdat == SISL_NULL) goto out;

  /* Set up local variables. */

  vpoint = (*pintdat) -> vpoint;
  ipoint = (*pintdat) -> ipoint;
  ipmax = (*pintdat) -> ipmax;
  vlist = (*pintdat) -> vlist;
  ilist = (*pintdat) -> ilist;
  ilmax = (*pintdat) -> ilmax;

  if(ipoint == 0 && ilist == 0) goto out;

  /* Make dimension of geometry space. */

  if (po1->iobj == SISLPOINT) idim = po1->p1->idim;
  else if (po1->iobj == SISLCURVE) idim = po1->c1->idim;
  else idim = po1->s1->idim;

  /* UJK, feb 93, often po1 and po2 is the same geometry ie po2
     is a dummy, so let create a dummy point in those cases. */
  if ((*pintdat) -> vpoint[0]->ipar == po1->iobj)
  {
     if ((qo2 = newObject(SISLPOINT)) == SISL_NULL) goto err101;
     knpar = po1->iobj;
  }
  else qo2 = po2;
  
     
     for(ki=0; ki<ilist; ki++)
  {
      /* Find out whether the curve vlist[ki] is degenerate.
      We assume that if all the guide points in vlist[ki]
      are equal then it is degenerate. This should be
      correct if SISL has done its job properly before
      reaching this stage of the intersection algorithm.
      Otherwise it clearly is not degenerate. */


      isdegen = TRUE;

      pfirst = vlist[ki]->pfirst;
      plast  = vlist[ki]->plast;
      knumbpt = vlist[ki]->inumb;

      /* VSK, 02-93. Check if the list describes a trimming area.
	 In that case no reduction is to be performed.            */
      
      if (pfirst->iinter == SI_TRIM) continue;
      
      kstat = 0;
      sh6degen_geom(po1,qo2,pfirst->epar,geomfirst,&kstat);
      if(kstat < 0) goto error;

      pprev = pfirst;

      pcurr = sh6getnext(pprev,vlist[ki]->ind_first);

      /* Loop through the guide points and compare with pfirst. */

      while(pcurr != plast)
      {
	 kstat = 0;
	  sh6degen_geom(po1,qo2,pcurr->epar,geomcurr,&kstat);
	  if(kstat < 0) goto error;

	  len = s6dist(geomcurr,geomfirst,idim);
	  if(len > aepsge)
	  {
	      isdegen = FALSE;
	      break;
	  }

	  sh6getother(pcurr,pprev,&pnext,&kstat);
	  if(kstat < 0) goto error;
  
	  pprev = pcurr;
	  pcurr = pnext;
      }
      
      /* Compare plast with pfirst. */

      if(isdegen == TRUE)
      {
	 kstat = 0;
	  sh6degen_geom(po1,qo2,plast->epar,geomcurr,&kstat);
	  if(kstat < 0) goto error;

	  len = s6dist(geomcurr,geomfirst,idim);
	  if(len > aepsge)
	  {
	      isdegen = FALSE;
	  }
      }
      
      /* Test if it is only 2 points in the intersection list. In that
	 case it may be a cyclic constant parameter curve, and it is
	 necessary to check in an internal point. */
      
      if (isdegen && knumbpt == 2)
      {
	 for (kexact=0, kj=0; kj<knpar; kj++) 
	 {
	    if (pfirst->curve_dir[vlist[ki]->ind_first] &
		(1 << (kj + 1))) 
	       kexact = kj + 1;
	    spar[kj] = (double)0.5*(pfirst->epar[kj]+plast->epar[kj]);
	 
	 }
	 
	 /* UJK, 930811:
	    kstat = (kexact == 0 || kexact > po1->iobj) ? 1 : 0; */
	 kstat = (kexact > po1->iobj) ? 1 : 0;
	 sh6degen_geom(po1,qo2,spar,geommid,&kstat);
	 
	 len = s6dist(geommid,geomfirst,idim);
	 if(len > aepsge)
	 {
	    isdegen = FALSE;
	 }
	 else
	 {
	    len = s6dist(geommid,geomcurr,idim);
	    if(len > aepsge)
	    {
	       isdegen = FALSE;
	    }
	 }
	 
      }

      if(isdegen == TRUE)
      {
	  /* Curve vlist[ki] is degenerate, i.e. it's just a point.  
	     Free the list in pintdat, but make sure that there is one
	     point left. If there is a point that is an endpoint of an
	     other list, and that are equal to this degenerate list,
	     this endpoint should represent these point.
	     
	     First pass through the list and mark all points as removed.  */

	 pcurr = pfirst;
	 kstoremark = pfirst->marker;
	 pnext = sh6getnext(pcurr,vlist[ki]->ind_first);
	 kpoint = vlist[ki]->inumb;
	 
	 for (kl=0; kl<kpoint; kl++)
	   {
	      /* Mark the point as removed.  */
	      
	      pcurr->marker = -99;
		 
	       pprev = pcurr;
	       pcurr = pnext;
	       
	       if (kl < kpoint-1)
	       {
		  sh6getother(pcurr,pprev,&pnext,&kstat);
		  if (kstat < 0) goto error;
	       }
	    }
	 
	 /* Free the list.  */
	 
	 freeIntlist(vlist[ki]);
	 vlist[ki] = SISL_NULL;
	 
	  /* Check the intersection points. If no point represent the
	     degenerated curve,restore the first point of the curve. */

          for(kj=0; kj<ipoint; kj++)
	  {
	     qpt = vpoint[kj];
	     if (qpt->marker == -99) continue;
	     kstat = 0;
	    sh6degen_geom(po1,qo2,qpt->epar,geomcurr,&kstat);
	    if(kstat < 0) goto error;
    
	    len = s6dist(geomfirst,geomcurr,idim);
	    if(len <= aepsge) break;
	  }
	  if (kj >= ipoint) pfirst->marker = kstoremark;

      }
  }

  /*  Tidy up the vlist array afterwards since some of the
      lists have been freed. */

  newilist = 0;

  for(ki=0; ki<ilist; ki++)
  {
      if(vlist[ki] != SISL_NULL)
      {
	  newilist++;
	  if(newilist < ki+1)
	  {
	      vlist[newilist-1] = vlist[ki];
	      vlist[ki] = SISL_NULL;
	  }
      }
  }

  ilist = newilist;

  /* Pass through the points and remove all marked points. */
  
  for (kj=0; kj<(*pintdat) -> ipoint; )
  {
     qpt = (*pintdat)->vpoint[kj];
     if (qpt->marker == -99)
     {
	/* Remove point.  */
	
	sh6idkpt (pintdat, &qpt, 0, &kstat);
	if (kstat < 0) goto error;
     }
     else kj++;
  }
  
  /* Reset the integer pintdat variables with the local
  ones in case they have changed. */

  (*pintdat)->ipmax =  ipmax;
  (*pintdat)->ilist =  ilist;
  (*pintdat)->ilmax =  ilmax;


  /* Degeneracies removed. Finish. */

  *jstat = 0;
  goto out;

  /* Error in alloc. */
err101:
  *jstat = -101;
  s6err ("sh6degen", *jstat, kpos);
  goto out;

  /* Error in lower level routine. */
error:
  *jstat = kstat;
  s6err ("sh6degen", *jstat, kpos);
  goto out;

out:
   if (qo2 && qo2 != po2) freeObject(qo2);

}

   
#if defined(SISLNEEDPROTOTYPES)
static void
  sh6degen_geom(SISLObject *po1,SISLObject *po2,double epar[],
		double geom[],int *jstat)
#else
static void sh6degen_geom(po1,po2,epar,geom,jstat)
      SISLObject *po1;
      SISLObject *po2;
      double epar[];
      double geom[];
      int *jstat;
#endif      
/*
*********************************************************************
*
* PURPOSE    : Evaluate object in an intersection point between the
*              two objects, po1 and po2. 
*
*
* INPUT      : po1  - First object in intersection.
*              po2  - Second object, in the case of a surface/curve - analytic
*                     intersection, po2 = po1 = object pointing at
*                     the surface or curve.
*              epar - Parameter values of intersection point inherited
*                     from the actual intersection point.
*              jstat - = 0 : Evaluate first object.
*                      = 1 : Evaluate second object.
*              
*
* OUTPUT     : geom - Position of intersection point in geometry space.
*              jstat   - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
*********************************************************************
*/
{
   int kstat = 0;             /* Local status variable.             */
   int kdim;                  /* Dimension.                         */
   int kder = 0;              /* Evaluate position.                 */
   int kleft1 = 0;            /* Parameter to the evaluator.        */
   int kleft2 = 0;            /* Parameter to the evaluator.        */
   int keval = *jstat;        /* Indicates which object to evaluate. */
   double sder[3];            /* Position and derivatives.          */
   double snorm[3];           /* Dummy normal in surf. evalutation. */

   *jstat = 0;
   if (keval == 0)
   {
      if (po1->iobj == SISLSURFACE)
      {
	 /* Evaluate the surface.  */
	 
	 kdim = po1->s1->idim;
	 s1421(po1->s1,kder,epar,&kleft1,&kleft2,sder,snorm,&kstat);
	 if (kstat < 0) goto error;
      }
      
      else if (po1->iobj == SISLCURVE)
      {
	 
	 /* Evaluate the curve.  */
	 
	 kdim = po1->c1->idim;
	 s1221(po1->c1,kder,epar[0],&kleft1,sder,&kstat);
	 if (kstat < 0) goto error;
      }
      
      else if (po1->iobj == SISLPOINT)
      {
	 
	 /* Copy the value of the point. */
	 
	 kdim = po1->p1->idim;
	 memcopy(sder,po1->p1->ecoef,kdim,DOUBLE);
      }
   }
   else if (keval == 1)
   {
      if (po2->iobj == SISLSURFACE)
      {
	 /* Evaluate the surface.  */
	 
	 kdim = po2->s1->idim;
	 s1421(po2->s1,kder,epar+po1->iobj,&kleft1,&kleft2,sder,snorm,&kstat);
	 if (kstat < 0) goto error;
      }
      
      else if (po2->iobj == SISLCURVE)
      {
	 
	 /* Evaluate the curve.  */
	 
	 kdim = po2->c1->idim;
	 s1221(po2->c1,kder,epar[po1->iobj],&kleft1,sder,&kstat);
	 if (kstat < 0) goto error;
      }
      
      else if (po2->iobj == SISLPOINT)
      {
	 
	 /* Copy the value of the point. */
	 
	 kdim = po2->p1->idim;
	 memcopy(sder,po1->p1->ecoef,kdim,DOUBLE);
      }
   }
      

   /* Copy result to output.  */
   
   memcopy(geom,sder,kdim,DOUBLE);
   
   *jstat = 0;
   goto out;
   
   /* Error in lower level function. */
   
   error :
      *jstat = kstat;
   goto out;
   
   out :
      return;
}
