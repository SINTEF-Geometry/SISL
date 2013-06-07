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
 * $Id: s1302.c,v 1.2 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1302

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1302(SISLCurve *pc,double aepsge,double angle,double ep[],double eaxis[],
	   SISLSurf **rs,int *jstat)
#else
void s1302(pc,aepsge,angle,ep,eaxis,rs,jstat)
     SISLCurve  *pc;
     double aepsge;
     double angle;
     double ep[];
     double eaxis[];
     SISLSurf   **rs;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To create a B-spline rotational surface by rotating
*              the curve *pc around the axisdefined by e}[] and eaxis[]
*              the given angle. The maximal deviation between the true
*              rotational surface and the generated surface allowed 
*              is controlled by aepsge.
*             
*
* INPUT      : pc     - Pointer to curve to be rotated
*              aepsge - Maximal deviation allowed between true rotational
*                       surface and generated surface.
*              angle  - The rotational angle. Counter clockwise around axis.
*                       If the absolute value of the angle is greater than
*                       2 PI then a rotational surface closed in the
*                       rotation direction is made.
*              ep     - SISLPoint on rotational axis
*              eaxis  - Direction of rotational axis
*
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              rs     - Pointer to the surface produced
*
* METHOD     : First the maximal distance between the curve and the
*              rotational axis is determined. Then by comparing this
*              with aepsge the allowed relative error is found. This
*              relative error and the rotational angle is used for
*              generating a normalized circle segment spanning the
*              actual angle. This circle is then translated to generate
*              the actual rows of control vertices of the surface
*
* REFERENCES :
*
*-                                      
* CALLS      : s6norm, s6scpr, s1301, s6rotax, s6mvec,
*              newSurf, test_cyclic_knots, s6err
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 24. May 1988
* REVISED BY : Christophe Birkeland, SI, Oslo, Norway. 13. August 1992
*              (Generates Nurbs-surface if aepsge = 0)
* REVISED BY : Christophe Rene Birkeland, SINTEF, Oslo, May 1993.
*              (jstat = kstat after call to s1520)
*********************************************************************
*/
{
  double *st1;            /* Pointer to knot vector of circle segment   */
  double *scoef1;         /* Pointer to vertices of circle segment      */
  int    kn1;             /* Number of vertice of circle segment        */
  int    kk1;             /* Order of circle segment                    */
  double *st2;            /* Pointer to knot vector of input curve      */
  double *scoef2;         /* Pointer to vertices of input curve         */
  int    kn2;             /* Number of vertice of input curve           */
  int    kk2;             /* Order of input curve                       */
  int    kdim;            /* Dimension of space in which curve lies     */
  double sdiff[3];        /* Array for storing differences              */
  double saxis[3];        /* Array for storing nomalize eaxis           */
  int    kl;              /* Pointer into array                         */
  int    kj;              /* Control variable in loop                   */
  int    ki;              /* Control variable in loop                   */
  double tlength;         /* Variable used for length calculation       */
  double tmaxl;           /* Variable used for accumulating max lengths */
  double treler;          /* Variable used for relative error           */
  double *sucof=SISL_NULL;     /* Pointer to vertex array for surface        */
  SISLCurve *pnorm=SISL_NULL;  /* Pointer to normalized circle segment       */
  double smat[16];        /* Transformation matrix                      */
  int    kstat;           /* Status variable                            */
  double tfak;            /* Value of cross product                     */
  double *srow;           /* Pointer to row of vertices in surface      */
  
  int    kpos = 1;       /* Position of error                     */
  

  /* If aepsge = 0.0, a nurbs surface is generated  */
  
  if (aepsge < REL_COMP_RES)
  {
     s1520(pc,angle,ep,eaxis,rs,&kstat);
     if (kstat < 0)
	goto error;
     *jstat = kstat;
     goto out;
  }
	     
     
  /* Make local pointers to description of curve */
  
  st2    = pc -> et;
  scoef2 = pc -> ecoef;
  kn2    = pc -> in;
  kk2    = pc -> ik;
  kdim  = pc -> idim;
  
  /* The routine is only working for dimension=3 */
  
  if (kdim != 3) goto err104;
  
  
  /* Normalize axis direction */                
  (void)s6norm(eaxis,kdim,saxis,&kstat);
  if (kstat<0) goto error;
  
  
  /* Find maximal distance between axis and vertices of curve */
  
  tmaxl = (double)0.0;
  
  for (ki=0;ki<kn2;ki++)
    {
      
      /*  Make difference between vertex and point on axis */
      kl = ki*kdim;
      for (kj=0;kj<3;kj++)
        {
	  sdiff[kj] = scoef2[kl] - ep[kj];
	  kl++;
        }
      tfak = s6scpr(sdiff,saxis,kdim);
      
      /*  Find vector normal to axis going to vertex by subtracing the
       *   component of the difference vector along the axis                */
      for (kj=0;kj<3;kj++)
        {
	  sdiff[kj] = sdiff[kj] - tfak*saxis[kj];
        }
      
      /*  Find length of this vector */
      tlength = s6norm(sdiff,kdim,sdiff,&kstat);
      if (kstat<0) goto error;
      tmaxl = MAX(tmaxl,tlength);
    }
  
  /* Calculate relative error, if this is <=0 then the whole curve
   *  lies on the axis. */
  
  if (tmaxl <= (double)0.0) goto err127;
  treler = aepsge/tmaxl;
  
  
  /* Calculate normalized circle */
  
  s1301(treler,angle,kdim,&pnorm,&kstat);
  if (kstat<0) goto error;
  
  /* Make local variables for curve description */
  
  st1    = pnorm -> et;
  scoef1 = pnorm -> ecoef;
  kn1    = pnorm -> in;
  kk1    = pnorm -> ik;
  
  /* Allocate vertex array for surface */
  sucof = newarray(kn1*kn2*kdim,DOUBLE);
  if (sucof == SISL_NULL) goto err101;
  
  /* Make the surface vertices circle segment by circle segment */
  
  for (ki=0;ki<kn2;ki++)
    {
      
      /*  Make transformation matrix for first vertex on curve to be rotated */
      
      s6rotax(ep,eaxis,&scoef2[ki*kdim],smat,&kstat);
      if (kstat<0) goto error;
      
      /*  Transform the vertices of this row into right position */
      
      srow = sucof + ki*3*kn1;
      s6mvec(smat,scoef1,kn1,srow);
    }    
  
  /* Create the surface */
  
  *rs =  newSurf(kn1,kn2,kk1,kk2,st1,st2,sucof,1,kdim,1);
  
  /* Check if the surface is cyclic in the first parameter direction. */
  
  test_cyclic_knots(st1,kn1,kk1,&kstat);
  if (kstat < 0) goto error;
  if (kstat == 2) (*rs)->cuopen_1 = SISL_SURF_PERIODIC;
  
  /* Copy periodicity flag from curve in second parameter direction. */
  
  (*rs)->cuopen_2 = pc->cuopen;
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1302",*jstat,kpos);
  goto out;
  
  /* Error in input, dimension not equal to 3 */
  
 err104: *jstat = -104;
  s6err("s1302",*jstat,kpos);
  goto out;
  
  /* Error in input, whole curve lies on axis */
  
 err127: *jstat = -127;
  s6err("s1302",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;     
  s6err("s1302",*jstat,kpos);
  goto out;
  
  
 out:
  
  /* Free allocated arrays */
  
  if (sucof != SISL_NULL) freearray(sucof);
  if (pnorm != SISL_NULL) freeCurve(pnorm);
  
  return;
}
