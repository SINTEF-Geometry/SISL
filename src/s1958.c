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
 * $Id: s1958.c,v 1.2 2001-03-19 15:58:58 afr Exp $
 *
 */


#define S1958

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1958(SISLSurf *psurf,double epoint[],int idim,double aepsco,double aepsge,
           double gpar[],double *dist,int *jstat)
#else
void s1958(psurf,epoint,idim,aepsco,aepsge,gpar,dist,jstat)
     SISLSurf    *psurf;
     double   epoint[];
     int      idim;
     double   aepsco;
     double   aepsge;
     double   gpar[];
     double   *dist;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find the closest point between the surface psurf
*              and the point epoint.
*              The method is fast and should work well
*              in clear cut cases but does not guarantee finding
*              the right solution. As long as it doesn't fail,
*              it will find exactly one point.
*
*
*
* INPUT      : psurf  - Pointer to the surface in the closest point problem.
*              epoint - The point in the closest point problem.
*              idim   - Dimension of the space in which epoint lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : gpar   - Array(2) containing the parameter values of the
*                       closest point in the parameter space
*                       of the surface.
*              dist   - The closest distance between point and the surface.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : Solution in interior 
*                                         = 1      : Solution at an edge 
*                                         = 2      : Solution at a corner 
*                                         < 0      : error
*
*
* METHOD     : Find an initial guess solution by finding, essentially,
*              the closest control point to the given point
*              and estimating the corresponding parameter pair, s1959.
*              This pair is then the starting point for a Newton
*              iteration in s1773. The distance of this solution
*              is then compared with the distance of the edge curves
*              (from the given point) found by s1957
*              and the minimum is returned.
*
*
* REFERENCES :
*
*- 
* CALLS      : s1960,s1773,newPoint,freePoint,
*              s1957,freeCurve,s1436,s1437.
*
* WRITTEN BY : Michael Floater, SI, 91-10.
*
*********************************************************************
*/                                                               
{                                                                     
  double clspt[3];          /* Coeffs. of potential closest point.       */
  double newdist,cldist;    /* Distances of edges from epoint.           */
  double enext[2];          /* Initial guess for iteration               */
  double estart[2],eend[2]; /* Parameter area for Newton iteration.      */
  double *et1=SISL_NULL;         /* Knot vector in u direc.                   */
  double *et2=SISL_NULL;         /* Knot vector in v direc.                   */
  int ik1;                  /* Order of curve in u direction.            */
  int ik2;                  /* Order of curve in v direction.            */
  int in1;                  /* No. control points of curve in u direec.  */
  int in2;                  /* No. control points of curve in v direec.  */
  int kleft1=0,kleft2=0;    /* Dummies used in evualation.               */
  int kstat = 0;            /* Local status variable.                    */
  int clkstat = 0;          /* Local status variable.                    */
  int kpos = 0;             /* Error position.                           */
  double gpos[2];           /* Parameters of closest point on surface.   */
  double clgpos[2];         /* Current parameters of cl. pt. on surface. */
  double crvpar;            /* Parameter of closest point on an edge.    */
  SISLPoint *ppoint=SISL_NULL;   /* epoint in SISLPoint form.                 */
  SISLCurve *pcurve=SISL_NULL;   /* An edge of the surface.                   */
  
  /* Test input.  */
  
  if (idim != 3) goto err105;
  if (psurf->idim != idim) goto err106;


  /* Set up local variables. */

  ik1 = psurf->ik1;
  ik2 = psurf->ik2;
  in1 = psurf->in1;
  in2 = psurf->in2;
  et1 = psurf->et1;
  et2 = psurf->et2;

  /* First of all, find the closest point on the boundary.
     This is done by calling s1957 on each of the four edges. */

  /* Find closest point of edge 0. */
  
  s1437(psurf,et1[ik1-1],&pcurve,&kstat);
  if (kstat < 0) goto error;

  s1957(pcurve,epoint,idim,aepsco,aepsge,&crvpar,&cldist,&clkstat);
  if (kstat < 0) goto error;
  if(pcurve != SISL_NULL)
  {
    freeCurve(pcurve);
    pcurve = SISL_NULL;  
  }
  clgpos[0] = et1[ik1-1];
  clgpos[1] = crvpar;
  clkstat++;

  /* Find closest point of edge 1. */
  
  s1437(psurf,et1[in1],&pcurve,&kstat);
  if (kstat < 0) goto error;

  s1957(pcurve,epoint,idim,aepsco,aepsge,&crvpar,&newdist,&kstat);
  if (kstat < 0) goto error;
  if(pcurve != SISL_NULL)
  {
    freeCurve(pcurve);
    pcurve = SISL_NULL;
  }
  if(newdist < cldist)
  {
      cldist = newdist;
      clgpos[0] = et1[in1];
      clgpos[1] = crvpar;
      clkstat = kstat+1;
  }

  /* Find closest point of edge 2. */
  
  s1436(psurf,et2[ik2-1],&pcurve,&kstat);
  if (kstat < 0) goto error;

  s1957(pcurve,epoint,idim,aepsco,aepsge,&crvpar,&newdist,&kstat);
  if (kstat < 0) goto error;
  if(pcurve != SISL_NULL)
  {
    freeCurve(pcurve);
    pcurve = SISL_NULL;
  }

  if(newdist < cldist)
  {
      cldist = newdist;
      clgpos[0] = crvpar;
      clgpos[1] = et2[ik2-1];
      clkstat = kstat+1;
  }

  /* Find closest point of edge 3. */
  
  s1436(psurf,et2[in2],&pcurve,&kstat);
  if (kstat < 0) goto error;

  s1957(pcurve,epoint,idim,aepsco,aepsge,&crvpar,&newdist,&kstat);
  if (kstat < 0) goto error;
  if(pcurve != SISL_NULL) freeCurve(pcurve);

  if(newdist < cldist)
  {
      cldist = newdist;
      clgpos[0] = crvpar;
      clgpos[1] = et2[in2];
      clkstat = kstat+1;
  }



  /* Next, try the interior of the surface by Newton iteration. */

  ppoint = newPoint(epoint,idim,1);
  if(ppoint == SISL_NULL) goto err101;

  /* Find a good guess point based on finding the closest control
     point and its corresponding parameter values. */

  s1960(ppoint,psurf,enext,&kstat);
  if(kstat < 0) goto error;

  /* Do the Newton iteration. */
    
  estart[0] = et1[ik1-1];
  estart[1] = et2[ik2-1];
  eend[0] = et1[in1];
  eend[1] = et2[in2];

  s1773(ppoint,psurf,aepsge,estart,eend,enext,gpos,&kstat);
  if(kstat >= 0)
  {
    
      /* Closest point found. */
      /* Find distance and compare with edges. */
    
      s1424(psurf,0,0,gpos,&kleft1,&kleft2,clspt,&kstat);
      if (kstat < 0) goto error;
      
      newdist = s6dist(epoint,clspt,idim);

      if(newdist < cldist)
      {
          cldist = newdist;
          clgpos[0] = gpos[0];
          clgpos[1] = gpos[1];
          clkstat = 0;
      }

    
  }

  
  /* Closest point found. Return point is in interior, on an edge,
     or at a corner. */

  *dist = cldist;
  gpar[0] = clgpos[0];
  gpar[1] = clgpos[1];

  *jstat = clkstat;

  goto out;
 


  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1958",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 3.  */
  
 err105: *jstat = -105;
  s6err("s1958",*jstat,kpos);
  goto out;
  
  /* Dimensions conflicting.  */
  
 err106: *jstat = -106;
  s6err("s1958",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1958",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free allocated space.  */
  
  if (ppoint != SISL_NULL) freePoint(ppoint);
  
  return;
}                                               
                                           
                       

