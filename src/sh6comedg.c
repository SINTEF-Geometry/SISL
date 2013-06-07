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
 * $Id: sh6comedg.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6COMEDG

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
sh6comedg (SISLObject * po1, SISLObject * po2, SISLIntpt *pt1, SISLIntpt *pt2,
	   						int *jstat)
#else
void
sh6comedg (po1, po2, pt1, pt2, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt *pt1;
     SISLIntpt *pt2;
     int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Test if two intersection points is on a common edge and
*		connected along this edge on one
*		of the objects.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pt1      - Intersection point.
*              pt2      - Intersection point.
*
*
* OUTPUT     : jstat    - status messages
*                        < 0  : Error.
*			   0  : Not on a common edge.
*			   1  : On a common edge in first object.
*			   2  : On a common edge in second object.
*			   3  : On a common edge in bouth objects.
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
 int kstat=0;
 double minpar[4];
 double maxpar[4];
 int nrpar, np1, np2;
 int common_edg, i, j;
 int is_inside = 1;
 int on_edge1 = 0;
 int on_edge2 = 0;
 int ind1, ind2;
 /* ---------------------------------------------------------------- */
 
 *jstat = 0;
   
 if (pt1 != SISL_NULL && pt2 != SISL_NULL)
 {
    /* Making the parametric boarders */
 
    if (po1->iobj == SISLSURFACE)
    {
       nrpar = 2;
       np1 = 4;

       minpar[0] = po1->s1->et1[po1->s1->ik1-1];
       minpar[1] = po1->s1->et2[po1->s1->ik2-1];
       maxpar[0] = po1->s1->et1[po1->s1->in1];
       maxpar[1] = po1->s1->et2[po1->s1->in2];
    }
    else if (po1->iobj == SISLCURVE)
    {
       nrpar = 1;
       np1 = 2;

       minpar[0] = po1->c1->et[po1->c1->ik-1];
       maxpar[0] = po1->c1->et[po1->c1->in];
    }
    else    /* SISLPOINT */
       np1 = nrpar = 0;

    if (po2->iobj == SISLSURFACE)
    {
       minpar[nrpar]   = po2->s1->et1[po2->s1->ik1-1];
       minpar[nrpar+1] = po2->s1->et2[po2->s1->ik2-1];
       maxpar[nrpar]   = po2->s1->et1[po2->s1->in1];
       maxpar[nrpar+1] = po2->s1->et2[po2->s1->in2];
       nrpar += 2;
       np2 = 4;
    }
    else if (po2->iobj == SISLCURVE)
    {
       minpar[nrpar] = po2->c1->et[po2->c1->ik-1];
       maxpar[nrpar] = po2->c1->et[po2->c1->in];
       nrpar++;
       np2 = 2;
    }
    else   np2 = 0;     /* SISLPOINT */

    /* Testing. */
    
    /* UJK, aug.92 */
    /* for (i = 0; i < nrpar || !is_inside; i++) */
    for (i = 0; i < nrpar && is_inside; i++)
      {
	 if (pt1->epar[i] <= maxpar[i] + REL_PAR_RES &&
	     pt1->epar[i] >= minpar[i] - REL_PAR_RES)
	   {
	      /* pt1 is inside. */
	      
	      if (pt1->epar[i] >= maxpar[i] - REL_PAR_RES)
		on_edge1 +=  (1 << (2*i));	/* On edge/end */
	      if (pt1->epar[i] <= minpar[i] + REL_PAR_RES)
		on_edge1 +=  (1 << (2*i+1));	/* On edge/end */
	      
	   }
	 else  is_inside = 0;
	 
	 if (pt2->epar[i] <= maxpar[i] + REL_PAR_RES &&
	     pt2->epar[i] >= minpar[i] - REL_PAR_RES)
	   {
	      /* pt2 is inside. */
	      
	      if (pt2->epar[i] >= maxpar[i] - REL_PAR_RES)
		on_edge2 +=  (1 << (2*i));	/* On edge/end */
	      if (pt2->epar[i] <= minpar[i] + REL_PAR_RES)
		on_edge2 +=  (1 << (2*i+1));	/* On edge/end */
	      
	   }
	 else  is_inside = 0;
      }
    
    common_edg = on_edge1 & on_edge2;
    (*jstat) = 0;
    
    if(is_inside && common_edg)
    {
       if (np1 > 0)
       {
	  j = (15>>(4-np1));
	  if (common_edg & j)
	  {
	     sh6getlist(pt1,pt2,&ind1,&ind2,&kstat);
             if (kstat < 0) goto err106;
             if (kstat == 0)
	     {
		if (common_edg & 3) i = 2;
		else i = 0;
		if (common_edg & (3<<2)) i+= 4;
		if (pt1->curve_dir[ind1] & i)  (*jstat) = 1;
	     }
	  }
       }
       if (np2 > 0)
       {
	  j = (15>>(4-np2));
	  j <<= np1;
	  if (common_edg & j)
	  {
	     sh6getlist(pt1,pt2,&ind1,&ind2,&kstat);
             if (kstat < 0) goto err106;
             if (kstat == 0)
	     {
		if (common_edg & (3<<np1)) i = 8;
		else i = 0;
		if (common_edg & (3<<(np1+2))) i+= 16;
		if (pt1->curve_dir[ind1] & i)  (*jstat) += 2;
	     }
	  }
       }
    }
    else
       (*jstat) = 0;
 }
 else goto err108;

 
  /* Done. */

  goto out;

  /* Error in input. Conflicting dimensions.  */

err106:*jstat = -106;
  s6err("sh6comedg",*jstat,0);
  goto out;

  /* Error in input. No points */

err108:*jstat = -108;
  s6err("sh6comedg",*jstat,0);
  goto out;

out:
  return;
}

