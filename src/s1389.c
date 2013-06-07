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
 * $Id: s1389.c,v 1.2 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1389

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1389(SISLCurve *pc1,double *gcubic[],int *jnumb,int *jdim,int *jstat)
#else
void s1389(pc1,gcubic,jnumb,jdim,jstat)
     SISLCurve  *pc1;
     double *gcubic[];
     int    *jnumb;
     int    *jdim;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Convert a B-spline curve of order up to 4 to a sequence
*              of cubic segment with uniform parametrization
*
*                                            
* INPUT      : pc1    - Pointer to the curve to be converted
*
*
* OUTPUT     : gcubic - Array containing the sequence of cubic segments.
*                       Each segment is represented by the start point,
*                       followed by the start tangent, end point and end
*                       tangent. Each segment takes 4*jdim doubles.
*                       This array is allocated inside the function and must
*                       be released by the calling function.            
*
*              jstat  - status messages  
*                                         = 1      : Order too high, curve
*                                                    interpolated
*                                         = 0      : ok
*                                         < 0      : error
*                      
*
* METHOD     : For each polynomial segment end, the position and first
*              derivative is calculated.
*
* USE:         int knumb,kdim,kstat;
*              double *scubic;
*              SISLSurf *qc1;
*               .
*               .
*              s1389(qc1,&scubic,&knumb,&kdim,&kstat);
*               .
*              If one of the order of the curve is greater than four
*              (*jstat==1) then degree reduction (order reduction) should
*              be applied before using this routine to get a satisfactory
*              representation of the curve by Coons patches.
*              The degree reduction routine is s1343.
*
*
*
* REFERENCES :
*
*-
* CALLS      : s1221, s1227, s6err
*
* WRITTEN BY : Tor Dokken, SI, Norway, 1988-11
*
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* Position of error.                              */
  int kdim;           /* Dimension of the space in which the surface lies. */
  int kder=1;         /* Calculate all first derivatives                 */
  int kleft=0;        /* Pointer into knot vector                        */
  int kdumlft;        /* Temporary pointer into knot vector              */
  int ksize;          /* Number of doubles to store a cubic segment      */
  
  int ki;             /* Control variables in for loop                   */
  int kn;             /* Number of vertices                              */
  int kk;             /* Polynomial order                                */
  double *st;         /* Knots                                           */
  double tpar;        /* Current parameter value                         */
  double tparx;       /* Temporary parameter value                       */
  double tdiff1;      /* Length of parameter interval                    */
  double *scorn1;     /* Pointer to first end of segment                 */
  double *scorn2;     /* Pointer to second end of segment                */
  
  kn  = pc1->in;
  kk  = pc1->ik;
  kdim = pc1 -> idim;
  st  = pc1 -> et;
  
  
  /* Calculate number of doubles to store a cubic segment*/
  
  ksize = kdim*4;
  
  
  /* Allocate array for storage of the coefficients */
  
  *gcubic = newarray((kn*4*kdim),DOUBLE);
  if (*gcubic == SISL_NULL) goto err101;
  
  kleft = kk - 1;
  
  *jnumb = 0;
  
  scorn1 = *gcubic;
  
  while (kleft < kn)
    {
      
      /* Set pointers to the end of the segment */
      
      scorn2 = scorn1 + 2*kdim;
      
      tpar = st[kleft];
      
      /* The parameter describes the left corner of the segment. By 
	 evaluating at tpar we get the kleft to point to the parameter
	 interval. The other end is at st[kleft+1].
	 */
      
      /* Calulate start of segement */
      
      s1221(pc1,kder,tpar,&kleft,scorn1,&kstat);
      if (kstat<0) goto error;
      
      /* Find length of aprameter intervals */
      
      tdiff1 = st[kleft+1] - st[kleft];
      
      /* Calculate end of segment, us left derivative */
      
      tparx = st[kleft+1];
      kdumlft = kleft;
      
      s1227(pc1,kder,tparx,&kdumlft,scorn2,&kstat);
      if (kstat<0) goto error;
      
      /* Scale derivatives to match uniform parametrization */
      
      for (ki=kdim;ki<2*kdim;ki++)
        {
	  scorn1[ki] *= tdiff1;
	  scorn2[ki] *= tdiff1;
        }
      
      kleft += 1;
      *jnumb +=1;
      scorn1 += kdim*4;
    }
  
  /* The array is probably too big for the Coons patches, decrease the
     array size */
  
  
  /* Allocate array for storage of the coefficients */
  
  *gcubic = increasearray(*gcubic,((*jnumb)*4*kdim),DOUBLE);
  if (*gcubic == SISL_NULL) goto err101;
  
  
  /* Test if order to high */
  
  *jdim = kdim;
  
  if (kk>4) goto war01;
  
  *jstat = 0;
  
  goto out;
  
  /* Orders too high */
  
  war01:*jstat=1;
    goto out;
  
  /* Error in scratch allocation */
  
  err101: 
    *jstat = -101;
    s6err("s1389",*jstat,kpos);
    goto freeout;
  
  /* Error in lower level function.  */
  
  error : *jstat = kstat;
    s6err("s1389",*jstat,kpos); 
    goto freeout;
  
  /* Some error has occured free allocated space */
  
  freeout:
    if (*gcubic != SISL_NULL) freearray(*gcubic);
    goto out;
  
  out:
    return;
}

