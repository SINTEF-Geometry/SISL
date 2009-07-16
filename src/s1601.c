/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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
                    
