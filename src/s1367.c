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
 * $Id: s1367.c,v 1.2 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1366

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1367(SISLSurf *ps,double aoffset,double aepsge,int idim,double epar[],
	   int ider,int *ilfs,int *ilft,double eder[],int *jstat)
#else
void s1367(ps,aoffset,aepsge,idim,epar,ider,ilfs,ilft,eder,jstat)
     SISLSurf   *ps;
     double aoffset;
     double aepsge;
     int    idim;
     double epar[];
     int    ider;
     int    *ilfs;
     int    *ilft;
     double eder[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Evaluate the surface at the parameter pair value epar.
*
*
*
* INPUT      : ps     - The input B-spline surface.
*              aoffset- The offset distance.
*                       If idim=2 a positive signe on this value put the
*                       offset on the side of the positive normal vector,
*                       and a negative sign puts the offset on the sign
*                       of the negative normal vector.
*                       If idim=3 the offset is determined by the cross
*                       product of the tangent vector and the anorm vector.
*                       The offset distance is multiplied by this vector.
*              aepsge - Maximal deviation allowed between true offset surface
*                       and the approximated offset surface.
*              idim   - The dimension of the space (2 or 3).
*              epar   - Parameter pair value at which to compute position
*                       and derivatives.
*              ider   - The number of derivatives to compute.
*                       < 0: Error
*                       = 0: Compute position
*                       = 1: Compute position and first derivative
*                       etc.
* INPUT/OUTPUT:
*              ilfs   - Pointer to the interval in the knot vector in 1.
*                       parameter direction.
*              ilft   - Pointer to the interval in the knot vector in 2.
*                       parameter direction.
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         = 3      : Surface is degenerate
*                                                    at the point, normal
*                                                    found by s6degnorm
*                                         = 2      : Surface is degenerate
*                                                    at the point, normal
*                                                    has zero length
*                                         = 1      : Surface is close to
*                                                    degenerate at the point
*                                                    Angle between tangents,
*                                                    less than angular tolerance
*                                         = 0      : ok
*                                         < 0      : error
*              eder   - Array of dimension [(ider+1)*(ider+1)*idim]
*
* METHOD     : Points, derivatives, cross derivatives and normals are calcu-
*              lated at the selected point on the surface. The selected point
*              corresponds to the input parameter-pair (epar).
*              The offset values are calculated in the following way:
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*                    u          u                           u
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*                    v          v                           v
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*                    uv         uv                          uv
*
*
* EXAMPLE OF USE:
*
*
* REFERENCES :
*
*
*
* CALLS      : s1424     - Evaluate the surface at the parameter pair value
*                          (epar).
*              s6length, s6crss, s6scpr, s6err
*
* WRITTEN BY : Per Evensen,  SI, 89-4.
* REVISED BY : Michael Floater, SI, 3/2/92.
*              If an edge of the surface is degenerate, a normal is
*              still returned via s6degnorm.
*              Derivatives are set to null.
*
*********************************************************************
*/
{
  int kstat = 0;     /* Local status variable.                           */
  int kpos = 0;      /* Position of error.                               */
  int kdim;          /* Dimension of the space in which the surface lies.*/
  int ki;            /* Loop controller.                                 */
  
  int kder = ider+1; /* Number of derivatives to calculate.              */
  int knmb;          /* Number of array-elements.                        */
  double snorm[3];   /* Pointer to array containing the normal.          */
  double utang[3];   /* Pointer to array containing the u tangent. */
  double vtang[3];   /* Pointer to array containing the v tangent. */
  double sdum[27];   /* Array used instead of allocation.                */
  double *sder=SISL_NULL; /* Pointer to array containing the derivatives.     */
  
  if (idim !=3) goto err105;
  
  /* Only allocate space if sdum is too small. */
  
  kdim = ps->idim;
  knmb = kdim*kder*kder;
  if (knmb>27)
    sder = newarray(knmb,DOUBLE);
  else
    sder = &sdum[0];
 

  /* Set the answer to zero in case nothing is found,
     e.g. if ps is degenerate at the given point. */

  for (ki=0;ki<knmb;ki++) eder[ki] = DZERO;


  if (DEQUAL(aoffset,DZERO))
    {
      /* Evaluate the surface (ps) at the parameter value spar.  */
      
      s1424(ps,kder,kder,epar,ilfs,ilft,sder,&kstat);
      if (kstat<0) goto error;
      
      /* Copy relevant information to eder. */
      
      for (ki=0;ki<idim;ki++)
        {
	  eder[ki] = sder[ki];
	  eder[ki+kdim] = sder[ki+kdim];
	  eder[ki+2*kdim] = sder[ki+3*kdim];
	  eder[ki+3*kdim] = sder[ki+4*kdim];
        }
    }
  else
    {
      
      double tlen1,tlen2,tnorm,tang=(double)0.0;
      
      /* Evaluate the surface (ps) at the parameter value spar.  */
      
      s1424(ps,kder,kder,epar,ilfs,ilft,sder,&kstat);
      if (kstat<0) goto error;
      
      /*
       *   Let P(s,t) be a point on the B-spline surface, and N(s,t) the
       *   normal in P.
       *
       *   The normal is found by the expression: N = P  x P .
       *                                               s    t
       *
       *   The length of the normal N is found by:
       *
       *           1/2
       *   n = (N*N)   .
       */
      
      /* Calculate surface normal, make cross products of tangents. */
      
      s6crss(sder+kdim,sder+3*kdim,snorm);
      
      /*  Make length of tangents and normal */
      
      tlen1 = s6length(sder+kdim,kdim,&kstat);
      tlen2 = s6length(sder+3*kdim,kdim,&kstat);
      tnorm = s6length(snorm,kdim,&kstat);
      
      /*  Calculate angle between tangents */
      
      if (tlen1 != DZERO && tlen2 != DZERO && tnorm != DZERO)
        tang = tnorm/(tlen1*tlen2);

      
      if (tang == (double)0.0) *jstat = 2;
      else if (tang <= ANGULAR_TOLERANCE) *jstat = 1;   
      else *jstat = 0;
      

      /* I just chose a reasonable tolerance here.
	 It should be relative to the size of the surface.
	 10e-8 ought to work OK if the size of ps is about
	 1 and we're
	 working in double precision. */

      if (tnorm < 0.00000001)
      {
	  /* Find a normal if possible and then leave. */

	  s6degnorm(ps,2,epar,sder,utang,vtang,snorm,&kstat);
	  if(kstat < 0) goto error;
	  if(kstat == 0)
	  {
              tnorm = s6length(snorm,kdim,&kstat);
	      *jstat = 3;
              for (ki=0;ki<idim;ki++)
		       eder[ki]=sder[ki]+aoffset*snorm[ki]/tnorm;
	  }
	  else
	  {
              for (ki=0;ki<idim;ki++) eder[ki]=sder[ki];
	  }
	  goto out;
      }

      /*
       *   The offset point, O(s,t,o), where o is the offset distance, is
       *   found by the expression:
       *
       *   O = P + o*(N/n)    
       */
      
      /* Calculate position of offset point. */
      
      for (ki=0;ki<idim;ki++) eder[ki]=sder[ki]+aoffset*snorm[ki]/tnorm;
      
      /* Calculate derivatives only if we are at a degenerate point. */

      if (ider > 0)
	{
	  
	  double *spnt;      /* Pointer to the derivatives.                */
	  double tder1,tder2;/* Lebgth of tangents.                        */
	  double snoru1[3];  /* The first derivative of the surface normal
				in first parameter direction, 1. part.     */
	  double snoru2[3];  /* The first derivative of the surface normal
				in first parameter direction, 2. part.     */
	  double snoru[3];   /* The first derivative of the surface normal
				in first parameter direction.              */
	  double snorv1[3];  /* The first derivative of the surface normal
				in second parameter direction, 1. part.    */
	  double snorv2[3];  /* The first derivative of the surface normal
				in second parameter direction, 2. part.    */
	  double snorv[3];   /* The first derivative of the surface normal
				in second parameter direction.             */
	  double snoruv1[3]; /* The cross derivative of the surface normal,
				1. part.                                   */
	  double snoruv2[3]; /* The cross derivative of the surface normal,
				2. part.                                   */
	  double snoruv3[3]; /* The cross derivative of the surface normal,
				3. part.                                   */
	  double snoruv[3];  /* The cross derivative of the surface normal.*/
	  double tnun;       /* The scalar product; snorm snoru.           */
	  double tnvn;       /* The scalar product; snorm snorv.           */
	  double tunvn;      /* The scalar product; snoru snorv.           */
	  double tnuvn;      /* The scalar product; snorm snoruv.          */
	  double snorus[3];  /* The first derivative of the surface normal
				divided by the length of the normal
				in first parameter direction.              */
	  double snorvs[3];  /* The first derivative of the surface normal
				divided by the length of the normal
				in second parameter direction.             */
	  double snoruvs[3]; /* The cross derivative of the surface normal
				divided by the length of the normal.       */
	  
	  /* Calculate length of tangents. */
	  
	  tder1 = s6length(sder+kdim,kdim,&kstat);
	  tder2 = s6length(sder+3*kdim,kdim,&kstat);
	  
	  /* The tangent length might be very different from 1. Scale it and
	     higher order derivatives to give tangent length one. */
	  
	  for (ki=0;ki<idim;ki++)
            {
	      spnt = sder+kdim+ki;
	      *spnt /= tder1;
	      
	      spnt = sder+2*kdim+ki;
	      *spnt /= (tder1+tder1);
	      
	      spnt = sder+3*kdim+ki;
	      *spnt /= tder2;
	      
	      spnt = sder+4*kdim+ki;
	      *spnt /= (tder1+tder2);
	      
	      spnt = sder+5*kdim+ki;
	      *spnt /= (tder1+tder1+tder2);
	      
	      spnt = sder+6*kdim+ki;
	      *spnt /= (tder2+tder2);
	      
	      spnt = sder+7*kdim+ki;
	      *spnt /= (tder1+tder2+tder2);
	      
	      spnt = sder+8*kdim+ki;
	      *spnt /= (tder1+tder1+tder2+tder2);
            }
	  
	  /*
	   *   The first derivative of the surface normal in first parameter direction
	   *   is calculated by the expression:
	   *
	   *
	   *   N  = P  x P   + P x P
	   *    u    uu   v     u   uv
	   */
	  
	  /* Calculate first derivative of the surface normal in first
	     parameter direction. */
	  
	  s6crss(sder+2*kdim,sder+3*kdim,snoru1);
	  s6crss(sder+kdim,sder+4*kdim,snoru2);
	  
	  for (ki=0;ki<idim;ki++) snoru[ki]=snoru1[ki]+snoru2[ki];
	  
	  /*
	   *   The first derivative of the surface normal in second parameter direction
	   *   is calculated by the expression:
	   *
	   *
	   *   N  = P  x P   + P x P
	   *    v    uv   v     u   vv
	   */
	  
	  /* Calculate first derivative of the surface normal in second
	     parameter direction. */
	  
	  s6crss(sder+4*kdim,sder+3*kdim,snorv1);
	  s6crss(sder+kdim,sder+6*kdim,snorv2);
	  
	  for (ki=0;ki<idim;ki++) snorv[ki]=snorv1[ki]+snorv2[ki];
	  
	  /*
	   *   The cross derivative of the surface normal is calculated by the 
	   *   expression:
	   *
	   *
	   *   N   = P   x P   + 2*(P  x P  ) + P  x P
	   *    uv    uuv   v        uu   vv     u    uvv
	   */
	  
	  /* Calculate cross derivative of the surface normal. */
	  
	  s6crss(sder+5*kdim,sder+3*kdim,snoruv1);
	  s6crss(sder+2*kdim,sder+6*kdim,snoruv2);
	  s6crss(sder+kdim,sder+7*kdim,snoruv3);
	  
	  for (ki=0;ki<idim;ki++) 
	    snoruv[ki]=snoruv1[ki]+(double)2.0*snoruv2[ki]+snoruv3[ki];
	  
	  
	  /*
	   *   The first derivative of the surface normal divided by the length of
	   *   the normal in first parameter direction is calculated by the 
	   *   expression:
	   *
	   *                             3
	   *   (N/n)  = N /n - (N * N )/n  N
	   *        u    u           u
	   */
	  
	  /* Calculate first derivative of the surface normal divided
	     by the length of the normal in first parameter direction. */
	  
	  tnun = s6scpr(snorm,snoru,idim);
	  
	  for (ki=0;ki<idim;ki++) 
	    snorus[ki]=
	      snoru[ki]/tnorm-(tnun/(tnorm*tnorm*tnorm))*snorm[ki];
	  
	  /*
	   *   The first derivative of the surface normal divided by the length of
	   *   the normal in second parameter direction is calculated by the 
	   *   expression:
	   *
	   *                             3
	   *   (N/n)  = N /n - (N * N )/n  N
	   *        v    v           v
	   */
	  
	  /* Calculate first derivative of the surface normal divided
	     by the length of the normal in second parameter direction. */
	  
	  tnvn = s6scpr(snorm,snorv,idim);
	  
	  for (ki=0;ki<idim;ki++) 
	    snorvs[ki]=
	      snorv[ki]/tnorm-(tnvn/(tnorm*tnorm*tnorm))*snorm[ki];
	  
	  /*
	   *   The cross derivative of the surface normal divided by the length of
	   *   the normal is calculated by expression:
	   *
	   *                               3                3
	   *   (N/n)   = N  /n - (N * N )/n  N  - (N * N )/n  N  -  ==>
	   *        uv    uv           v      u         u      v
	   *
	   *                        3                3                        5
	   *             (N  * N )/n  N - (N * N  )/n  N + 3(N * N )(N * N )/n  N
	   *               u    v               uv                u       v
	   */
	  
	  /* Calculate cross derivative of the surface normal divided
	     by the length of the normal. */
	  
	  tunvn = s6scpr(snoru,snorv,idim);
	  tnuvn = s6scpr(snorm,snoruv,idim);
	  
	  for (ki=0;ki<idim;ki++) 
	    snoruvs[ki]=
	      snoruv[ki]/tnorm-(tnvn/(tnorm*tnorm*tnorm))*snoru[ki]-
		(tnun/(tnorm*tnorm*tnorm))*snorv[ki]-
		  (tunvn/(tnorm*tnorm*tnorm))*snorm[ki]-
		    (tnuvn/(tnorm*tnorm*tnorm))*snorm[ki]+
		      (((double)3.0*tnun*tnvn)/(tnorm*tnorm*tnorm*tnorm*tnorm))*
			snorm[ki];
	  
	  /*
	   *   The first derivative in first parameter direction offset point, 
	   *   where o is the offset distance, is found by the expression:
	   *
	   *           
	   *   O  = P  + o*(N/n)    
	   *    u    u          u
	   */
	  
	  /* Calculate position of first derivative in first parameter direction
	     offset point. */
	  
	  for (ki=0;ki<idim;ki++) 
            eder[ki+kdim]=sder[ki+kdim]+aoffset*snorus[ki];
	  
	  /*
	   *   The first derivative in second parameter direction offset point, 
	   *   where o is the offset distance, is found by the expression:
	   *
	   *
	   *   O  = P  + o*(N/n)    
	   *    v    v          v
	   */
	  
	  /* Calculate position of first derivative in second parameter direction
	     offset point. */
	  
	  for (ki=0;ki<idim;ki++) 
            eder[ki+2*kdim]=sder[ki+3*kdim]+aoffset*snorvs[ki];
	  
	  /*
	   *   The cross derivative offset point, 
	   *   where o is the offset distance, is found by the expression:
	   *
	   *
	   *   O   = P   + o*(N/n)    
	   *    uv    uv          uv
	   */
	  
	  /* Calculate position of cross derivative offset point. */
	  
	  for (ki=0;ki<idim;ki++) 
            eder[ki+3*kdim]=sder[ki+4*kdim]+aoffset*snoruvs[ki];
	  
	  
	  
	  /* Scale the derivatives to match the original parametrization. */
	  
	  for (ki=0;ki<idim;ki++)
            {
	      spnt = eder+kdim+ki;
	      *spnt *= tder1;
	      
	      spnt = eder+2*kdim+ki;
	      *spnt *= tder2;
	      
	      spnt = eder+3*kdim+ki;
	      *spnt *= (tder1+tder2);
            }
	  
	}
    }
  
  /* SISLPoint and derivatives calculated, with a possible offset distance. */
  
  *jstat = 0;
  goto out;
  
  /* Error in input, dimension not equal to 3 */
  
  err105: *jstat = -105;
  s6err("s1367",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1367",*jstat,kpos);
  goto out;
  
  out: return;
}
