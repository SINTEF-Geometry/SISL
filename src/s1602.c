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
 * $Id: s1602.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1602

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1602(double estapt[],double endpt[],int ik,int idim,double astpar,
	   double *cendpar,SISLCurve **rc,int *jstat)
#else
void s1602(estapt,endpt,ik,idim,astpar,cendpar,rc,jstat)
     double estapt[];
     double endpt[];
     int    ik;
     int    idim;
     double astpar;
     double *cendpar;
     SISLCurve  **rc;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Convert a straight line to a B-spline described curve.
*             
*
* INPUT      : estapt - start point of the straight line 
*              endpt  - end point of the straight line
*              ik     - the order of the B-spline curve to be found
*              idim   - The dimension of the space
*              astpar - start value of parameterization of the curve
*             
* OUTPUT
*            : cendpar - parameter used at the end of the curve
*              rc      - Pointer to the found curve
*              jstat  - status messages
  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : First the knots are found with ik knots in the start
*              and ik knots in the end of the curve. Then a modified
*              Marsdens identity is used to find ik vertices equally 
*              spaced on the straight line.
*
* EXAMPLE OF USE:
*
* REFERENCES :
*
*-                                                 
* CALLS      : s6dist,newCurve,s6err
*
*
* WRITTEN BY : Qyvind Hjelle, SI, Oslo, Norway. 10. Nov 1988
*
*********************************************************************
*/
{
  int kit;            /* Loop control                                    */
  int kit2;           /* Loop contero                                    */
  int kvert;          /* Counter for position in vertex array            */
  int kpos=0;         /* Position of error                               */
  
  double *st=SISL_NULL;    /* Pointer to the first element of the knot vector
			 of the curve.                                   */
  double *scoef=SISL_NULL; /* Pointer to the first element of the curve's
			 B-spline coefficients.                          */
  double tdist;       /* Distance                                        */
  double tdel;        /* Delta x, y , ....                               */
  
  /* Check input */          
  
  if (idim <  1) goto err102;
  if (ik   <  2) goto err109;
  
  /* Find distance between start nd end point */
  tdist = s6dist(estapt,endpt,idim);
  
  
  /* Make knots. First allocate space */
  
  st = newarray(ik*2,DOUBLE);
  if (st == SISL_NULL) goto err101;
  
  for (kit=0; kit<ik; kit++) 
    {
      st[kit]    = astpar;
      st[kit+ik] = astpar + tdist;
    }
  
  /* calculate the vertices. First allocate space */ 
  
  /* First allocate space for vertices */ 
  
  scoef = newarray(ik*idim,DOUBLE);
  if (scoef == SISL_NULL) goto err101;
  
  /* Find first and last vertex. */ 
  
  kvert = (ik-1) * idim;
  for (kit=0; kit<idim; kit++,kvert++) 
    {
      scoef[kit]   = estapt[kit];
      scoef[kvert] = endpt[kit];
    }
  
  /* Find other vertices */ 
  
  for (kit=0; kit<idim; kit++)
    {   
      tdel = (endpt[kit] - estapt[kit])/(ik - 1);
      for (kit2=2; kit2<ik; kit2++)
	scoef[(kit2-1)*idim + kit] = scoef[(kit2-2)*idim + kit] + tdel; 
    }
  
  /* Make the curve */
  
  *rc = SISL_NULL;              
  *rc = newCurve(ik,ik,st,scoef,1,idim,1);
  if (*rc == SISL_NULL) goto err101;                
  
  *cendpar = st[ik];
  *jstat = 0;
  goto out;
  
  /* Error in memory allocation */
  
 err101: 
  *jstat = -101;
  s6err("s1602",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension less than 1 */
  
 err102: 
  *jstat = -102;
  s6err("s1602",*jstat,kpos);
  goto out;                          
  /* Error in input. Order less than 2 */
  
 err109: 
  *jstat = -109;
  s6err("s1602",*jstat,kpos);
  goto out;                          
    
 out:
  if (st     != SISL_NULL) freearray(st);
  if (scoef  != SISL_NULL) freearray(scoef);
  return;
}          
                    
