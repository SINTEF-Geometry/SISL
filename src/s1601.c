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
 * $Id: s1601.c,v 1.2 2001-03-19 15:58:51 afr Exp $
 *
 */


#define S1601

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1601(SISLSurf *psurf,double epoint[],double enorm[],int idim,SISLSurf **rsurf,int *jstat)
#else
void s1601(psurf,epoint,enorm,idim,rsurf,jstat)
     SISLSurf   *psurf;
     double epoint[];
     double enorm[];
     int    idim;
     SISLSurf   **rsurf;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To mirror a B-spline surface about a plane
*             
*
* INPUT      : psurf  - The input B-spline surface   
*              epoint - a point in the plane
*              enorm  - normal vector to the plane  
*              idim   - The dimension of the space
*
* OUTPUT     : rsurf  - Pointer to the mirrored surface
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : The vertices are mirrored. All other surface data
*              are copied into the new surface.  
*
* EXAMPLE OF USE:
*
* REFERENCES :
*
*-                                                 
* CALLS      : s6norm,s6diff,s6scpr,newSurf,s6err
*
*
* WRITTEN BY : Qyvind Hjelle, SI, Oslo, Norway. 10. Nov 1988
* REVISED BY : Johannes Kaasa, SI, Oslo, Norway April 1992 (Introduced NURBS)
*
*********************************************************************
*/
{
  int kstat=0;        /* Status variable                                 */
  int kiv;            /* Counter for vertices                            */
  int kid;            /* Counter for dimension                           */
  
  int kn1;            /* Number of vertices in first direction of psurf  */
  int kn2;            /* Number of vertices in second direction of psurf */
  int kk1;            /* Order in first direction of input psurf         */
  int kk2;            /* Order in second direction of input psurf        */
  int kind;           /* Type of surface, 2 and 4 are rational surfaces. */
  int kdim;           /* The dimension of the space in which the curve
			 lies. Equivalently, the number of components
			 of each B-spline coefficient.                   */
  int kvert;          /* Counter for position in mirrored vertex array   */
  int kpos=0;         /* Position of error                               */
  
  double *snorm=SISL_NULL; /* Pointer to normalized normal                    */
  double *svecpv=SISL_NULL;/* Pointer to vector from pointin plane to vertex  */
  double *smirr=SISL_NULL; /* Pointer to mirrored vertex array                */
  double *st1;        /* First knot vector is psurf                       */
  double *st2;        /* Second knot vector is psurf                      */
  double *svert;      /* Pointer to a vertex                             */
  double *scoef;      /* Pointer to the first element of the curve's
			 B-spline coefficients. This is assumed to be an
			 array with [kn*kdim] elements stored in the
			 following order:
			 First the kdim components of the first B-spline
			 coefficient, then the kdim components of the
			 second B-spline coefficient and so on.          */
  
  double tdist;       /* Distance */
  
  /* Make local pointers */
  
  kn1   = psurf -> in1;
  kn2   = psurf -> in2;
  kk1   = psurf -> ik1;
  kk2   = psurf -> ik2;
  kdim  = psurf -> idim;
  st1   = psurf -> et1;
  st2   = psurf -> et2;
  kind  = psurf -> ikind;
  if (kind == 2 || kind == 4)
    scoef = psurf -> rcoef;
  else
    scoef = psurf -> ecoef;
  
  /* Check if conflicting dimensions */  
  
  if (kdim != idim) goto err106;
  if (kind == 2 || kind == 4)
    kdim++;
  
  /* Allocate space for normalized normalvector , vector from point in
     plane to vertex and new vertex array */
  
  snorm = newarray(idim,DOUBLE);
  if (snorm == SISL_NULL) goto err101;
  
  svecpv = newarray(idim,DOUBLE);
  if (svecpv == SISL_NULL) goto err101;
  
  smirr = newarray(kdim*kn1*kn2,DOUBLE);
  if (svecpv == SISL_NULL) goto err101;
  
  
  /* Normalize normal vector */
  
  (void)s6norm(enorm,idim,snorm,&kstat);     
  if (kstat < 0) goto error;
  
  /* Do for all vertices */
  
  kvert = 0;
  for (kiv=0; kiv<kn1*kn2; kiv++) 
    {
      
      /* Find vector from point in plane to vertex */
      
      svert = scoef + kiv*kdim;
      s6diff(svert,epoint,idim,svecpv);
      
      /* Find distance between the vertex and the plane */
      
      tdist =  s6scpr(svecpv,snorm,idim);
      tdist *=(double)2.0;
      
      /* Find the mirrored vertex */
      
      for (kid=0; kid<idim; kid++,kvert++)
	{
	  smirr[kvert] = scoef[kvert] - tdist*snorm[kid];
	}
      if (kind == 2 || kind == 4)
	{
	   smirr[kvert] = scoef[kvert];
	   kvert++;
	}
    }
  
  
  /* Make the mirrored surface. */
  
  *rsurf = SISL_NULL;              
  *rsurf = newSurf(kn1,kn2,kk1,kk2,st1,st2,smirr,psurf->ikind,idim,1);
  if (*rsurf == SISL_NULL) goto err101;                
  
  /* Set periodicity flag. */
  
  (*rsurf)->cuopen_1 = psurf->cuopen_1;
  (*rsurf)->cuopen_2 = psurf->cuopen_2;
  
  *jstat = 0;
  goto out;
  
  /* Error in memory allocation */
  
 err101: 
  *jstat = -101;
  s6err("s1601",*jstat,kpos);
  goto out;
  
  /* Error in input, conflicting dimensions */
  
 err106: 
  *jstat = -106;
  s6err("s1601",*jstat,kpos);
  goto out;
  
  /* Error in lower level function */  
  
 error:  
  *jstat = kstat;
  s6err("s1601",*jstat,kpos); 
  goto out;
  
 out:
  if (snorm  != SISL_NULL) freearray(snorm);
  if (svecpv != SISL_NULL) freearray(svecpv);
  if (smirr  != SISL_NULL) freearray(smirr);
  return;
}    
                    
