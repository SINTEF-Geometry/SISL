/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl-copyright.h"

/*
 *
 * $Id: sh6red.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH6RED

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static int sh6red_help(SISLObject *po1, SISLObject *po2,
		       SISLIntpt *pt, SISLIntpt *pmain, int *jstat);
static void sh6red2help (SISLObject * po1, SISLObject * po2, 
			 SISLIntdat * pintdat, int *jstat);
#else
static int sh6red_help;
static void sh6red2help;
#endif

#if defined(SISLNEEDPROTOTYPES)
void
sh6red (SISLObject * po1, SISLObject * po2,
	SISLIntdat * pintdat, int *jstat)
#else
void
sh6red (po1, po2, pintdat, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat *pintdat;
     int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Reduse main int point to help point if they are not legal.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pintdat  - Intersection data structure.
*
*
* OUTPUT     : jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      :
*
* REFERENCES :
*
* WRITTEN BY : Arne Laksaa, SI, Oslo, sep. -91.
*********************************************************************
*/
{
  int kstat, i, j;
  double tepsge = (double)10000.0*REL_COMP_RES;
  double weight = (double) 0.5;
  int changed, doreduce;
  SISLIntpt *pcurr,*pstart,*plast;	/* to traverse list of points.     */
  int indstart,indlast,inddum;		/* Indexes used in lists           */
  int log_1, log_2;
  int dim;
  SISLIntpt *qmain;
  if (po1->iobj == SISLSURFACE)
    dim = po1->s1->idim;
  else if (po1->iobj == SISLCURVE)
    dim = po1->c1->idim;
  else if (po1->iobj == SISLPOINT)
    dim = po1->p1->idim;

  /* Remove all internal points in a list when along a
     constant parameter direction */
  
  if (((po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT
        && dim == 1) ||
       (po2->iobj == SISLSURFACE && po1->iobj == SISLPOINT
        && dim == 1) ||
       (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE
        && dim == 3) ||
      (po1->iobj + po2->iobj == SISLSURFACE + SISLCURVE
        && dim == 3) ||
      (po2->iobj == SISLSURFACE && po1->iobj == SISLCURVE
       && dim == 3)) &&
       pintdat != SISL_NULL)
     for (j = 0; j < pintdat->ipoint; j++)
     {
	
	pcurr = pintdat->vpoint[j];
	sh6isinside (po1, po2, pcurr, &kstat);
	if (kstat < 0)
	   goto error;
	
	/* VSK && ALA. 01.93. Do not remove points at corners. */
	if (kstat != 1 && kstat != 2) continue;
	
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
	   	   
	   if (log_1 & log_2 )
	   {
	      sh6idkpt (&pintdat, &pcurr, 1, &kstat);
	      if (kstat < 0)
		 goto error;
	      /* Recursive nature : */
	      j = -1;
	   }
	   
	   
	}
     }
   
  
  if (pintdat != SISL_NULL)
    {
      /* Weight value in 3D sf vs sf case is one */
      if (pintdat->vpoint[0]->ipar == 4)
	weight = (double) 1.0;

      /* Reduce an illegal trim_curve to one point. */

      for (i = 0; i < pintdat->ipoint; i++)
	{
	  if (pintdat->vpoint[i]->iinter == SI_TRIM)
	    {
	      SISLIntpt **trim = SISL_NULL;
	      int no_trim = 0;
	      int no_alloc = 0;
	      sh6trimlist (pintdat->vpoint[i], &trim, &no_trim, &no_alloc);
	      for (j = 0; j < no_trim; j++)
		{
		  sh6isinside (po1, po2, trim[j], &kstat);
		  if (kstat < 0)
		    goto error;
		  if (kstat != 1)
		    break;
		}
	      if (j == no_trim)
		{
		  /* Internal trim area. */
		  for (j = 1; j < no_trim; j++)
		    {
		       /* sh6idunite (&pintdat, &trim[0], &trim[j], weight, &kstat);
			  */
		       /* VSK. 01.93. */
		       sh6idnewunite(po1, po2, &pintdat, &trim[0], &trim[j], 
				     weight, tepsge, &kstat);
		      if (kstat < 0)
			goto error;

		      /* We now need to correct the intpoint. */
		    }
		  trim[0]->iinter = SI_SING;
		}
	      if (trim)
		freearray (trim);
	    }
	}

      /* Reduse ilegal main points to help points. */
      do
	{
	  changed = 0;
	  //for (i = 0; i < pintdat->ipoint; i++)
	  for (i = pintdat->ipoint-1; i >= 0; i--)
	    {
	      sh6isinside (po1, po2, pintdat->vpoint[i], &kstat);
	      if (kstat < 0)
		goto error;
	      if (kstat == 1 || (kstat == 2 && po1->iobj+po2->iobj == 2*SISLSURFACE))
		{
		  if (sh6ismain (pintdat->vpoint[i]) &&
		      sh6nmbmain (pintdat->vpoint[i], &kstat) == 1)
		    {
		      if (pintdat->vpoint[i]->fromhelp > 0)
			doreduce = 1;
		      else
			{
			  /* Check if a reduction is recommended */
			  qmain = SISL_NULL;
			  for (j=0; j<pintdat->vpoint[i]->no_of_curves; ++j)
			    if (pintdat->vpoint[i]->pnext[j]->iinter > 0)
			      {
				qmain = pintdat->vpoint[i]->pnext[j];
				break;
			      }

			  if (qmain && po1->iobj+po2->iobj == 2*SISLSURFACE)
			    {
			      /* Check if the points are connected along an
				 edge */
			      sh6comedg(po1, po2, pintdat->vpoint[i], qmain, &kstat);
			      if (kstat < 0)
				goto error;
			      if (kstat == 0)
				qmain = SISL_NULL;   /* Not a common edge */
			    }
			  doreduce = (qmain == SISL_NULL) ? 0 :
			    sh6red_help(po1, po2, pintdat->vpoint[i], qmain, &kstat);
			  if (kstat < 0)
			    goto error;
			}
		      if (doreduce)
			{
			  sh6tohelp (pintdat->vpoint[i], &kstat);
			  if (kstat < 0)
			    goto error;
			  changed = 1;
			}
		    }
		}
	    }
      } while (changed);
 
      /*UJK, 12.08.93 */
      /* Disconnect trim pts with 3 neighbours */
      do
	{
	   int ind_1,ind_2;
	   SISLIntpt *p_neighb[3];
	   int log_check[3];
	   changed = 0;
	   for (i = 0; i < pintdat->ipoint; i++)
	   {
	      pcurr = pintdat->vpoint[i];
	      sh6isinside (po1, po2, pcurr, &kstat);
	      if (kstat < 0)
		 goto error;
	      if (kstat &&
		  pcurr->iinter == SI_TRIM &&
		  sh6nmbmain (pcurr, &kstat) == 3)
	      {
		 for (ind_1=ind_2=0;ind_1<pcurr->no_of_curves;ind_1++)
		    if (pcurr->pnext[ind_1]->iinter == SI_TRIM)
		    {
		       sh6isinside (po1, po2, pcurr->pnext[ind_1], &kstat);
		       if (kstat < 0)
			  goto error;
		       if (kstat)
		       {
			  p_neighb[ind_2]  = pcurr->pnext[ind_1];
			  log_check[ind_2] = pcurr->curve_dir[ind_1];
			  log_check[ind_2] = log_check[ind_2]>>1;
			  log_check[ind_2] &= 15;
			  ind_2++;
		       }
		    }
		 
		 if (ind_2 == 3)
		 {
		    if (log_check[0] & log_check[1])
		       ind_2 = 2;
		    else if (log_check[0] & log_check[2])
		       ind_2 = 1;
		    else if (log_check[1] & log_check[2])
		       ind_2 = 0;
		    
		    if (ind_2 < 3)
		    {
		       changed = TRUE;
		       sh6disconnect(pcurr,p_neighb[ind_2],&kstat);
		       if (kstat < 0) goto error;
		       /* afr: Changed line below from an empty if-statement. */
		       sh6nmbmain (p_neighb[ind_2], &kstat);
		       if (kstat < 0) goto error;
		       sh6idkpt (&pintdat, &p_neighb[ind_2], 0, &kstat);
		       if (kstat < 0) goto error;
		    }
		 }
	      }
	   }
	} while (changed);
    }

  /* VSK, 03-2018. Take another iteration to remove internal points in a
     list after the help point reduction is performed.  */
  if (pintdat != SISL_NULL)
    {
     for (j = 0; j < pintdat->ipoint; j++)
     {
	
	pcurr = pintdat->vpoint[j];
	sh6isinside (po1, po2, pcurr, &kstat);
	if (kstat < 0)
	   goto error;
	
	/* Do not remove points at corners. */
	if (kstat != 1 && kstat != 2) continue;
	
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
	   	   
	   if ((log_1 & log_2) || 
	       (dim == 3 && po1->iobj+po2->iobj < 4))
	   {
	      sh6idkpt (&pintdat, &pcurr, 1, &kstat);
	      if (kstat < 0)
		 goto error;
	      /* Recursive nature : */
	      j = -1;
	   }
	}
     }
    }

  /* Reduce illegal main points to help points. */
  /* sh6red2help(po1, po2, pintdat, &kstat); */
  /* if (kstat < 0) */
  /*   goto error; */

   /* Reduction done. */

  (*jstat) = 0;
  goto out;

errinconsist:
  *jstat = -500;
  s6err ("sh6red", *jstat, 0);
  goto out;
  
error:(*jstat) = kstat;
  s6err ("sh6red", *jstat, 0);
  goto out;

out:
  return;
}

#if defined(SISLNEEDPROTOTYPES)
static int sh6red_help(SISLObject *po1, SISLObject *po2,
		       SISLIntpt *pt, SISLIntpt *pmain, int *jstat)
#else
  static int sh6red_help(po1, po2, pt, pmain, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt *pt;
     SISLIntpt *pmain;
     int *jstat;
#endif
{
  int kstat = 0;
  int ki, k1;
  int kleft1 = 0, kleft2 = 0, kleft3;
  double tpar1, tpar2;
  double *sst[4];   /* Knot vectors in the problem, maximum 4. */
  int lkn[4];       /* Number of coefficients in each parameter direction. */
  int lkk[4];       /* Order in each parameter direction. */
  int  doreduce = 1;

  *jstat = 0;
  if (po1->iobj == SISLSURFACE)
    {
      sst[0] = po1->s1->et1;
      sst[1] = po1->s1->et2;
      lkn[0] = po1->s1->in1;
      lkn[1] = po1->s1->in2;
      lkk[0] = po1->s1->ik1;
      lkk[1] = po1->s1->ik2;
      k1 = 2;
    }
  else if (po1->iobj == SISLCURVE)
    {
      sst[0] = po1->c1->et;
      lkn[0] = po1->c1->in;
      lkk[0] = po1->c1->ik;
      k1 = 1;
    }
  else
    k1 = 0;

  if (po2->iobj == SISLSURFACE)
    {
      sst[k1] = po2->s1->et1;
      sst[k1+1] = po2->s1->et2;
      lkn[k1] = po2->s1->in1;
      lkn[k1+1] = po2->s1->in2;
      lkk[k1] = po2->s1->ik1;
      lkk[k1+1] = po2->s1->ik2;
    }
  else if (po2->iobj == SISLCURVE)
    {
      sst[k1] = po2->c1->et;
      lkn[k1] = po2->c1->in;
      lkk[k1] = po2->c1->ik;
    }
      
  for (ki=0; ki<pt->ipar; ki++)
    {
      tpar1 = pt->epar[ki];
      s1219(sst[ki], lkk[ki], lkn[ki], &kleft1, tpar1, &kstat);
      if (kstat < 0)
	goto error;
      if (tpar1 >= sst[ki][lkn[ki]])
	kleft1 = lkn[ki];

      tpar2 = pmain->epar[ki];
      s1219(sst[ki], lkk[ki], lkn[ki], &kleft2, tpar2, &kstat);
      if (kstat < 0)
	goto error;
      if (tpar2 >= sst[ki][lkn[ki]])
	kleft2 = lkn[ki];

      kleft3 = (kleft2 > kleft1) ? kleft1 + 1 :
	((kleft1 == kleft2) ? kleft1 : kleft1 - 1);
      if (DEQUAL(sst[ki][kleft1],tpar1) ||
	  (kleft1 != kleft2 && DNEQUAL(sst[ki][kleft3],sst[ki][kleft2])))
	{
	  doreduce = 0;
	  break;
	}
    }
  goto out;

 error:
  *jstat = kstat;
  goto out;

 out:
  return doreduce;
}
	

#if defined(SISLNEEDPROTOTYPES)
static void sh6red2help (SISLObject * po1, SISLObject * po2, 
			 SISLIntdat * pintdat, int *jstat)
#else
static void sh6red2help (po1, po2, pintdat, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat *pintdat;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Reduse main int point to help point if they are not legal.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pintdat  - Intersection data structure.
*
*
* OUTPUT     : jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*

*********************************************************************
*/
{
  int change = 0;
  int kstat = 0, gstat = 0;
  int i;
  SISLIntpt *pt1 = SISL_NULL;
  SISLIntpt *pt2 = SISL_NULL;
  SISLIntpt *pcurr = SISL_NULL;

  if (po1->iobj + po2->iobj == 2*SISLSURFACE && pintdat != SISL_NULL)
    {
  do
    {
      change = 0;
      //for (i = 0; i < pintdat->ipoint; i++)
      for (i = pintdat->ipoint-1; i >= 0; i--)
	{
	  if (sh6ismain (pintdat->vpoint[i]))
	    {
	      sh6isinside (po1, po2, pintdat->vpoint[i], &kstat);
	      if (kstat < 0)
		goto error;
	      
	      /* ALA and VSK. Test if the point lies on edge in 
		 one or two objects.         */
	      if (kstat == 2 || kstat == 5)
		{
		  sh6getnhbrs (pintdat->vpoint[i], &pt1, &pt2, &gstat);
		  if (gstat == 1)
		    {
		      sh6comedg (po1, po2, pintdat->vpoint[i], pt1, &gstat);
		      
		      /* ALA and VSK. Test if the points lie on the same
			 edge in both objects if it lies on an edge in
			 both objects.                  */
		      if ((kstat == 2 && gstat > 0) ||
			  (kstat == 5 && gstat == 3))
			{
			  sh6tohelp (pintdat->vpoint[i], &kstat);
			  if (kstat < 0)
			    goto error;
			  change = 1;
			}
		    }
		}
	    }
	}
    } while (change);
    }
  
  /* Reduction done. */

  (*jstat) = 0;
  goto out;

error:(*jstat) = kstat;
  s6err ("sh6red2help", *jstat, 0);
  goto out;

out:
  return;
}
