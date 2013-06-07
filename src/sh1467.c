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
 * $Id: sh1467.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH1467

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static void sh1467_s9fac1(double [],double [],double [],int,int,int,double *,
			  double *,double *);
static void sh1467_s9fac2(double [],double [],double [],int,int,int,double *,
			  double *,double *,double *,double *,double *);
#else
  static void  sh1467_s9fac1(); 
  static void  sh1467_s9fac2(); 
#endif  

#if defined(SISLNEEDPROTOTYPES)
void
  sh1467(SISLCurve *ecurve[],double etwist[],int ider,double ebar[],
	 double eval[],int *jstat)
#else	   
void sh1467(ecurve,etwist,ider,ebar,eval,jstat)
     double etwist[],ebar[],eval[];
     SISLCurve *ecurve[];
     int ider,*jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given the barycentric coordinates of a point in a 5-sided
*              vertex region, evaluate the value of the ideal blending 
*              surface of the vertex region in this point.
*
*
*
* INPUT      : ecurve - Position and cross-tangent curves around the vertex
*                       region. For each edge of the region position and cross-
*                       tangent curves are given. The curves follow each other
*                       around the region and are oriented counter-clock-wise.
*                       The dimension of the array is 10.
*              etwist - Twist-vectors of the corners of the vertex region. The
*                       first element of the array is the twist in the corner
*                       before the first edge, etc. The dimension of the array
*                       is 5*kdim.
*              ider   - Number of derivatives to compute. Directions of 
*                       differentiation is that of the two first barycentric
*                       coordinates. 0 <= ider <= 2.
*              ebar   - Generalized barycentric coordinates of the point to be 
*                       evaluated. The dimension of the array is 5.
*                       
*
* OUTPUT     : eval   - Value and derivatives of ideal blending surface in the 
*                       parameter directions of the two first coordinates in the 
*                       given point. Dimension of the array is 3*(1+..+(ider+1)).
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : A functional description of the ideal surface is given as
*              a blend between five surfaces, each of which fulfill the
*              continuity requirements over two edges.
*
* REFERENCES : Gregory and Charrot : A pentagonal surface patch for
*                                    computer aided geometric design
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1221 - Evaluate curve at a given parameter value.
*              sh1467_s9fac1() - Compute value and derivative of expression 
*			         in the barycentric coordinates.  
*              sh1467_s9fac2() - Compute value, first and second derivative.
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
  int kstat=0;         /* Status variable.                                */
  int ki,kj,kk,kh;     /* Counters.                                       */
  int kder;            /* Number of derivatives to evaluate.              */
  int kleft = 0;       /* Local parameter used in s1221.                  */
  int kdim = 3;        /* Dimension of geometry.                          */
  int kwarn = 0;       /* Indicates if a warning is to be sendt.          */
  int knmb;            /* Number of doubles pr derivative.                */
  int kl = 0;          /* Number of derivatives in the output array.      */
  double tdum;         /* Product of three barycentric coordinates.       */
  double tdenom;       /* Denominator of weight function.                 */
  double salpha[5];    /* Weigth of actual blending surface.              */
  double sp[15];       /* Value of function interpolating two sides.      */
  double sstart[5];    /* Start-parameters of edge-curves.                */
  double send[5];      /* End-parameters of edge-curves.                  */
  double sint[5];      /* Parameter intervals of edge-curves.             */
  double sx1[5];       /* First local coordinates of the blending patches.  */
  double sx2[5];       /* Second local coordinates of the blending patches. */
  double spar[5];      /* Parameter values of edges.                      */
  double spos[45];     /* Position of edge-curve.                         */
  double sder[45];     /* Value of oss tangent.                           */
  double scorn[15];    /* Position of edge-curves in actual corner of 
			  vertex region.                                  */
  double scornder1[15]; /* Tangent in corner along next edge.             */
  double scornder2[15]; /* Tangent in corner along previous edge.         */

  /* Test input.  */

  if (ider > 2) kwarn = 1;
  
  /* Initialise.  */

  kder = ider;
  knmb = kdim*(ider+1);
  for (ki=0; ki<ider; ki++) kl += ider + 1;
  
  /* Initiate output array to zero.  */

  for (kh=0; kh<kl*kdim; kh++) eval[kh] = (double)0.0;
        
  for (ki=0; ki<5; ki++)
    {

      /* Get endpoints of parameter intervals of edge curves.  */

      sstart[ki] = *(ecurve[2*ki] -> et + ecurve[2*ki] -> ik - 1);
      send[ki] = *(ecurve[2*ki] -> et + ecurve[2*ki] -> in);
      sint[ki] = send[ki] - sstart[ki];

      /* Compute local coordinates corresponding to one blending patch. */
      
      kj = (ki + 2) % 5;
      kk = (ki > 0) ? ki-1 : 4;
      sx1[ki] = ebar[kj]/(ebar[kk] + ebar[kj]);
      
      kj = (ki > 1) ? ki-2 : 3+ki;
      kk = (ki + 1) % 5;
      sx2[ki] = ebar[kj]/(ebar[kk] + ebar[kj]);
      
      /* Compute parameter value corresponding to the current edge. */

      spar[kj] = ((double)1.0 - sx1[ki])*sstart[ki] + sx1[ki]*send[ki];

      /* Evaluate position and cross-tangent curves at points on
	 the edges needed when evaluating surface.        */
 
      /* Evaluate position and cross-tangent curves at parameter value.  */

      s1221(ecurve[2*ki],kder,spar[kj],&kleft,spos+knmb*kj,&kstat);
      if (kstat < 0) goto error;
      
      s1221(ecurve[2*ki+1],kder,spar[kj],&kleft,sder+knmb*kj,&kstat);
      if (kstat < 0) goto error;
 
      /* Evaluate position and both cross-tangents at the corner nr ki.  */

      kk = (ki > 0) ? ki-1 : 4;
      s1221(ecurve[2*ki],0,sstart[ki],&kleft,scorn+ki*kdim,&kstat);
      if (kstat < 0) goto error;

      s1221(ecurve[2*kk+1],0,send[kk],&kleft,scornder1+ki*kdim,&kstat);
      if (kstat < 0) goto error;

      s1221(ecurve[2*ki+1],0,sstart[ki],&kleft,scornder2+ki*kdim,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* Compute position of the Gregory Charrot patch at the given parameter
     value.  First compute factors in the expression of the weight function. */

  tdenom = (double)0.0;
  for (kh=0; kh<5; kh++)
    {
      kj = (kh + 1) % 5;
      kk = (kh > 0) ? kh-1 : 4;	  

      tdum = ebar[kk]*ebar[kh]*ebar[kj];
      tdenom += tdum*tdum;
    }

  for (ki=0; ki<5; ki++)
    {

      kj = (ki + 1) % 5;
      kk = (ki > 0) ? ki-1 : 4;

      /* Compute the weigth of the actual blending surface.  */

      tdum = ebar[kk]*ebar[ki]*ebar[kj];
      salpha[ki] = tdum*tdum/tdenom;

      /* Compute the blending surface.  */

      kj = (ki + 2) % 5;
      kk = (ki > 1) ? ki-2 : 3+ki;

      for (kh=0; kh<kdim; kh++)
	{
	  sp[ki*kdim+kh] = spos[knmb*kj+kh] + sx1[ki]*sder[knmb*kj+kh] +
	    spos[knmb*kk+kh] + sx2[ki]*sder[knmb*kk+kh] -
	      scorn[ki*kdim+kh] - sx1[ki]*scornder1[ki*kdim+kh] 
		- sx2[ki]*scornder2[ki*kdim+kh] -
		  sx1[ki]*sx2[ki]*etwist[ki*kdim+kh];
	  
	  /* Add the contribution of the value of this blending surface to the
	     value of the ideal surface.   */

	  eval[kh] += salpha[ki]*sp[ki*kdim+kh];
	}
    }
  
  if (ider >= 1)
    {
      /* Compute first derivatives of the Gregory Charrot function. */

      double tconst = (sqrt((double)5.0) - (double)1.0)/(double)2.0; /* Constant. */
      double tdum1,tdum2; /* Factor in expression for derivative of weight function.*/
      double tnom1,tnom2; /* Factor in expression for derivative of weight function.*/
      double sd1alpha[5],sd2alpha[5]; /* 1. derivative of weight functions. */
      double sd1p[15],sd2p[15];       /* 1. derivative of blending functions. */
      double sd1bar[5],sd2bar[5];     /* 1. derivative of barycentric coordinates.*/
      double sd1x1[5],sd2x1[5];       /* 1. derivative of local coordinates.  */      
      double sd1x2[5],sd2x2[5];       /* 1. derivative of local coordinates.  */      
      double sd1par[5],sd2par[5];     /* 1. derivative of parameter value at edge. */

      /* Differentiate barycentric coordinates.  */

      sd1bar[0] = (double)1.0;
      sd1bar[1] = (double)0.0;
      sd1bar[2] = -(double)1.0;
      sd1bar[3] = -tconst;
      sd1bar[4] = tconst;

      sd2bar[0] = (double)0.0;
      sd2bar[1] = (double)1.0;
      sd2bar[2] = tconst;
      sd2bar[3] = -tconst;
      sd2bar[4] = -(double)1.0;

      for (ki=0; ki<5; ki++)
	{
	  /* Differentiate local coordinates.  */      

	  kj = (ki + 2) % 5;
	  kk = (ki > 0) ? ki-1 : 4;
	  tdum = ebar[kk] + ebar[kj];
	  sd1x1[ki] = sd1bar[kj]/tdum
	    - ebar[kj]*(sd1bar[kk] + sd1bar[kj])/(tdum*tdum);
	  sd2x1[ki] = sd2bar[kj]/tdum
	    - ebar[kj]*(sd2bar[kk] + sd2bar[kj])/(tdum*tdum);

	  kj = (ki > 1) ? ki-2 : 3+ki;
	  kk = (ki + 1) % 5;
	  tdum = ebar[kk] + ebar[kj];
	  sd1x2[ki] = sd1bar[kj]/tdum
	    - ebar[kj]*(sd1bar[kk] + sd1bar[kj])/(tdum*tdum);
	  sd2x2[ki] = sd2bar[kj]/tdum
	    - ebar[kj]*(sd2bar[kk] + sd2bar[kj])/(tdum*tdum);

	  /* Compute derivative of parameter value corresponding to
	     the current edge.  */

	  sd1par[kj] = sint[ki]*sd1x1[ki];
	  sd2par[kj] = sint[ki]*sd2x1[ki];
	}
      
      /* Compute factors used in the expression for the 1. derivatives of
	 the weight function.  */

      tnom1 = tnom2 = (double)0.0;
      for (kh=0; kh<5; kh++)
	{
	  kj = (kh + 1) % 5;
	  kk = (kh > 0) ? kh-1 : 4;	  
	      
	  sh1467_s9fac1(ebar,sd1bar,sd2bar,kk,kh,kj,&tdum,&tdum1,&tdum2);

	  tnom1 += (double)2.0*tdum*tdum1;
	  tnom2 += (double)2.0*tdum*tdum2;
	}
      

      for (ki=0; ki<5; ki++)
	{
	  /* Compute the 1. derivatives of the weight functions.  */

	  kj = (ki + 1) % 5;
	  kk = (ki > 0) ? ki-1 : 4;	  
	  
	  sh1467_s9fac1(ebar,sd1bar,sd2bar,kk,ki,kj,&tdum,&tdum1,&tdum2);    

	  sd1alpha[ki] = tdum*((double)2.0*tdum1 - tdum*tnom1/tdenom)/tdenom;
	  sd2alpha[ki] = tdum*((double)2.0*tdum2 - tdum*tnom2/tdenom)/tdenom;
      
	  /* Compute 1. derivatives of the functions which blends two sides
	     of the region.  */

	  kj = (ki + 2) % 5;
	  kk = (ki > 1) ? ki-2 : 3+ki;	  

	  for (kh=0; kh<kdim; kh++)
	    {
	      sd1p[ki*kdim+kh] = spos[kj*knmb+kdim+kh]*sd1par[kj]
		+ sd1x1[ki]*sder[kj*knmb+kh] 
		  + sx1[ki]*sder[kj*knmb+kdim+kh]*sd1par[kj]
		    + spos[kk*knmb+kdim+kh]*sd1par[kk]
		      + sd1x2[ki]*sder[kk*knmb+kh]
			+ sx2[ki]*sder[kk*knmb+kdim+kh]*sd1par[kk]
			  - sd1x1[ki]*scornder1[ki*kdim+kh]
			    - sd1x2[ki]*scornder2[ki*kdim+kh]
			      - (sd1x1[ki]*sx2[ki] 
				 + sx1[ki]*sd1x2[ki])*etwist[ki*kdim+kh];

	      sd2p[ki*kdim+kh] = spos[kj*knmb+kdim+kh]*sd2par[kj]
		+ sd2x1[ki]*sder[kj*knmb+kh] 
		  + sx1[ki]*sder[kj*knmb+kdim+kh]*sd2par[kj]
		    + spos[kk*knmb+kdim+kh]*sd2par[kk]
		      + sd2x2[ki]*sder[kk*knmb+kh]
			+ sx2[ki]*sder[kk*knmb+kdim+kh]*sd2par[kk]
			  - sd2x1[ki]*scornder1[ki*kdim+kh]
			    - sd2x2[ki]*scornder2[ki*kdim+kh]
			      - (sd2x1[ki]*sx2[ki] 
				 + sx1[ki]*sd2x2[ki])*etwist[ki*kdim+kh];

	      /* Add the contribution from this blending surface to the first 
		 derivative of the Gregory Charrot function. */

	      eval[kdim+kh] += sd1alpha[ki]*sp[ki*kdim+kh]
		+ salpha[ki]*sd1p[ki*kdim+kh]; 

	      eval[2*kdim+kh] += sd2alpha[ki]*sp[ki*kdim+kh]
		+ salpha[ki]*sd2p[ki*kdim+kh]; 
	    }
	}
      
      if (ider >= 2)
	{
	  double sd11alpha[5],sd12alpha[5],sd22alpha[5];  /* 2. derivatives of
							     weight function.  */
	  double sd11p[15],sd12p[15],sd22p[15];  /* 2. derivatives of blending
						    surfaces.                  */
	  double sd11x1[5],sd12x1[5],sd22x1[5];  /* 2. derivatives of local
						    coordinates.               */
	  double sd11x2[5],sd12x2[5],sd22x2[5];  /* 2. derivatives of local
						    coordinates.               */
	  double sd11par[5],sd12par[5],sd22par[5];  /* 2. derivatives of
						       parameter value at edge. */
	  double tdum11,tdum12,tdum22;  /* Factor in expression for 2.
					   derivative of weight function.*/
	  double tnom11,tnom12,tnom22;  /* Factor in expression for 2.
					   derivative of weight function.*/
	  
	  /* Compute second derivatives of the Gregory Charrot function. */

	  for (ki=0; ki<5; ki++)
	    {
	      /* Differentiate local coordinates.  */      

	      kj = (ki + 2) % 5;
	      kk = (ki > 0) ? ki-1 : 4;
	      tdum = ebar[kk] + ebar[kj];
	      tdum1 = sd1bar[kk] + sd1bar[kj];
	      tdum2 = sd2bar[kk] + sd2bar[kj];
	  
	      sd11x1[ki] = (double)2.0*(ebar[kj]*tdum1*tdum1/tdum
					- sd1bar[kj]*tdum1)/(tdum*tdum);
	      sd12x1[ki] = ((double)2.0*ebar[kj]*tdum1*tdum2/tdum
			    - sd1bar[kj]*tdum2 - sd2bar[kj]*tdum1)/(tdum*tdum);
	      sd22x1[ki] = (double)2.0*(ebar[kj]*tdum2*tdum2/tdum
					- sd2bar[kj]*tdum2)/(tdum*tdum);

	      kj = (ki > 1) ? ki-2 : 3+ki;
	      kk = (ki + 1) % 5;
	      tdum = ebar[kk] + ebar[kj];
	      tdum1 = sd1bar[kk] + sd1bar[kj];
	      tdum2 = sd2bar[kk] + sd2bar[kj];
	  
	      sd11x2[ki] = (double)2.0*(ebar[kj]*tdum1*tdum1/tdum
					- sd1bar[kj]*tdum1)/(tdum*tdum);
	      sd12x2[ki] = ((double)2.0*ebar[kj]*tdum1*tdum2/tdum
			    - sd1bar[kj]*tdum2 - sd2bar[kj]*tdum1)/(tdum*tdum);
	      sd22x2[ki] = (double)2.0*(ebar[kj]*tdum2*tdum2/tdum
					- sd2bar[kj]*tdum2)/(tdum*tdum);

	      /* Compute derivative of parameter value corresponding to
		 the current edge.  */

	      sd11par[kj] = sint[ki]*sd11x1[ki];
	      sd12par[kj] = sint[ki]*sd12x1[ki];
	      sd22par[kj] = sint[ki]*sd22x1[ki];
	    }

      /* Compute factors used in the expression for the 2. derivatives of
	 the weight function.  */

	  tnom11 = tnom12 = tnom22 = (double)0.0;
	  for (kh=0; kh<5; kh++)
	    {
	      kj = (kh + 1) % 5;
	      kk = (kh > 0) ? kh-1 : 4;	  
	      
              sh1467_s9fac2(ebar,sd1bar,sd2bar,kk,kh,kj,&tdum,&tdum1,&tdum2,
		     &tdum11,&tdum12,&tdum22);
	      
	      tnom11 += (double)2.0*(tdum1*tdum1 + tdum*tdum11);
	      tnom12 += (double)2.0*(tdum1*tdum2 + tdum*tdum12);
	      tnom22 += (double)2.0*(tdum2*tdum2 + tdum*tdum22);

	    }

	  for (ki=0; ki<5; ki++)
	    {
	      
	      /* Compute the 2. derivatives of the weight functions.  */

	      kj = (ki + 1) % 5;
	      kk = (ki > 0) ? ki-1 : 4;	  
	      
	      sh1467_s9fac2(ebar,sd1bar,sd2bar,kk,ki,kj,&tdum,&tdum1,&tdum2,
		     &tdum11,&tdum12,&tdum22);
	      
	      sd11alpha[ki] = (double)2.0*(tdum1*tdum1 + tdum*tdum11)/tdenom
		- (double)4.0*tdum*tdum1*tnom1/(tdenom*tdenom)
		  - tdum*tdum*tnom11/(tdenom*tdenom) 
                     + (double)2.0*tdum*tdum*tnom1*tnom1/(tdenom*tdenom*tdenom);

	      sd12alpha[ki] = (double)2.0*(tdum1*tdum2 + tdum*tdum12)/tdenom
		- (double)2.0*tdum*(tdum1*tnom2 + tdum2*tnom1)/(tdenom*tdenom)
		  - tdum*tdum*tnom12/(tdenom*tdenom)
		    + (double)2.0*tdum*tdum*tnom1*tnom2/(tdenom*tdenom*tdenom);

	      sd22alpha[ki] = (double)2.0*(tdum2*tdum2 + tdum*tdum22)/tdenom
		- (double)4.0*tdum*tdum2*tnom2/(tdenom*tdenom)
		  - tdum*tdum*tnom22/(tdenom*tdenom)
		    + (double)2.0*tdum*tdum*tnom2*tnom2/(tdenom*tdenom*tdenom);
	      
	      /* Compute 2. derivatives of the functions which blends two sides
		 of the region.  */

	      kj = (ki + 2) % 5;
	      kk = (ki > 1) ? ki - 2 : 3 + ki;
	      
	      for (kh=0; kh<kdim; kh++)
		{
		  sd11p[ki*kdim+kh] = spos[kj*knmb+2*kdim+kh]*sd1par[kj]*sd1par[kj]
		    + spos[kj*knmb+kdim+kh]*sd11par[kj]
		      + sd11x1[ki]*sder[kj*knmb+kh]
			+ (double)2.0*sd1x1[ki]*sder[kj*knmb+kdim+kh]*sd1par[kj]
			  + sx1[ki]*sder[kj*knmb+2*kdim+kh]*sd1par[kj]*sd1par[kj]
			    + sx1[ki]*sder[kj*knmb+kdim+kh]*sd11par[kj]
			      + (double)2.0*sd1x2[ki]*sder[kk*knmb+kdim+kh]*sd1par[kk]
				+ sx2[ki]*sder[kk*knmb+2*kdim+kh]*sd1par[kk]*sd1par[kk]
				  + spos[kk*knmb+2*kdim+kh]*sd1par[kk]*sd1par[kk]
				    + spos[kk*knmb+kdim+kh]*sd11par[kk]
				      + sd11x2[ki]*sder[kk*knmb+kh]
					+ sx2[ki]*sder[kk*knmb+kdim+kh]*sd11par[kk]
					  - sd11x1[ki]*scornder1[ki*kdim+kh]
					    - sd11x2[ki]*scornder2[ki*kdim+kh]
					      - (sd11x1[ki]*sx2[ki] 
						 + (double)2.0*sd1x1[ki]*sd1x2[ki]
						 + sx1[ki]*sd11x2[ki])
						*etwist[ki*kdim+kh];
	      
		  sd12p[ki*kdim+kh] = spos[kj*knmb+2*kdim+kh]*sd1par[kj]*sd2par[kj]
		    + spos[kj*knmb+kdim+kh]*sd12par[kj]
		      + sd12x1[ki]*sder[kj*knmb+kh]
			+ sd1x1[ki]*sder[kj*knmb+kdim+kh]*sd2par[kj]
			  + sd2x1[ki]*sder[kj*knmb+kdim+kh]*sd1par[kj]
			    + sx1[ki]*sder[kj*knmb+2*kdim+kh]*sd1par[kj]*sd2par[kj]
			      + sx1[ki]*sder[kj*knmb+kdim+kh]*sd12par[kj]
				+ sx2[ki]*sder[kk*knmb+2*kdim+kh]*sd1par[kk]*sd2par[kk]
				  + spos[kk*knmb+2*kdim+kh]*sd1par[kk]*sd2par[kk]
				    + spos[kk*knmb+kdim+kh]*sd12par[kk]
				      + sd12x2[ki]*sder[kk*knmb+kh]
					+ sd1x2[ki]*sder[kk*knmb+kdim+kh]*sd2par[kk]
					  + sd2x2[ki]*sder[kk*knmb+kdim+kh]*sd1par[kk]
					    + sx2[ki]*sder[kk*knmb+kdim+kh]*sd12par[kk]
					      - sd12x1[ki]*scornder1[ki*kdim+kh]
						- sd12x2[ki]*scornder2[ki*kdim+kh]
						  - (sd12x1[ki]*sx2[ki] 
						     + sd1x1[ki]*sd2x2[ki]
						     + sd2x1[ki]*sd1x2[ki]
						     + sx1[ki]*sd12x2[ki])
						    *etwist[ki*kdim+kh];

		  sd22p[ki*kdim+kh] = spos[kj*knmb+2*kdim+kh]*sd2par[kj]*sd2par[kj]
		    + spos[kj*knmb+kdim+kh]*sd22par[kj]
		      + sd22x1[ki]*sder[kj*knmb+kh]
			+ (double)2.0*sd2x1[ki]*sder[kj*knmb+kdim+kh]*sd2par[kj]
			  + sx1[ki]*sder[kj*knmb+2*kdim+kh]*sd2par[kj]*sd2par[kj]
			    + sx1[ki]*sder[kj*knmb+kdim+kh]*sd22par[kj]
			      + (double)2.0*sd2x2[ki]*sder[kk*knmb+kdim+kh]*sd2par[kk]
				+ sx2[ki]*sder[kk*knmb+2*kdim+kh]*sd2par[kk]*sd2par[kk]
				  + spos[kk*knmb+2*kdim+kh]*sd2par[kk]*sd2par[kk]
				    + spos[kk*knmb+kdim+kh]*sd22par[kk]
				      + sd22x2[ki]*sder[kk*knmb+kh]
					+ sx2[ki]*sder[kk*knmb+kdim+kh]*sd22par[kk]
					  - sd22x1[ki]*scornder1[ki*kdim+kh]
					    - sd22x2[ki]*scornder2[ki*kdim+kh]
					      - (sd22x1[ki]*sx2[ki] 
						 + (double)2.0*sd2x1[ki]*sd2x2[ki]
						 + sx1[ki]*sd22x2[ki])
						*etwist[ki*kdim+kh];
	      

		  /* Add the current contribution to the 2. derivative of the Gregory 
		     Charrot function. */

		  eval[3*kdim+kh] += sd11alpha[ki]*sp[ki*kdim+kh]
		    + (double)2.0*sd1alpha[ki]*sd1p[ki*kdim+kh]
		      + salpha[ki]*sd11p[ki*kdim+kh];
		  
		  eval[4*kdim+kh] += sd12alpha[ki]*sp[ki*kdim+kh]
		    + sd1alpha[ki]*sd2p[ki*kdim+kh] 
		      + sd2alpha[ki]*sd1p[ki*kdim+kh]
			+ salpha[ki]*sd12p[ki*kdim+kh];
		  
		  eval[5*kdim+kh] += sd22alpha[ki]*sp[ki*kdim+kh]
		    + (double)2.0*sd2alpha[ki]*sd2p[ki*kdim+kh]
		      + salpha[ki]*sd22p[ki*kdim+kh];
		}
	    }
	}
    }
  
  /* Ideal surface evaluated.  */

  *jstat = kwarn;
  goto out;


  /* Error in lower level function.  */

  error :
    *jstat = kstat;
  goto out;

  out :
    return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
  sh1467_s9fac1(double ebar[],double ed1bar[],double ed2bar[],
		int i1,int i2,int i3,double *cfac,double *cfac1,
		double *cfac2)
#else	       
static void sh1467_s9fac1(ebar,ed1bar,ed2bar,i1,i2,i3,cfac,cfac1,cfac2)
     int i1,i2,i3;
     double ebar[],ed1bar[],ed2bar[],*cfac,*cfac1,*cfac2;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute value and 1. derivative of product of barycentric
*              coordinates.
*
*
*
* INPUT      : ebar   - Array containing value of barycentric coordinates.
*              ed1bar - Array containing derivative in first parameter 
*                       direction of barycentric coordinates.
*              ed2bar - Array containing derivative in second parameter 
*                       direction of barycentric coordinates.
*              i1     - Index of first barycentric coordinate in product.
*              i2     - Index of second barycentric coordinate in product.
*              i3     - Index of third barycentric coordinate in product.
*                       
*
* OUTPUT     : cfac   - Product of barycentric coordinates.
*              cfac1  - Derivative in first parameter direction of product.
*              cfac2  - Derivative in second parameter direction of product.
*
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
  *cfac = ebar[i1]*ebar[i2]*ebar[i3];
  *cfac1 = ed1bar[i1]*ebar[i2]*ebar[i3] + ebar[i1]*ed1bar[i2]*ebar[i3]
    + ebar[i1]*ebar[i2]*ed1bar[i3];
  *cfac2 = ed2bar[i1]*ebar[i2]*ebar[i3] + ebar[i1]*ed2bar[i2]*ebar[i3]
    + ebar[i1]*ebar[i2]*ed2bar[i3];

  return;
}

  
#if defined(SISLNEEDPROTOTYPES)
static void
  sh1467_s9fac2(double ebar[],double ed1bar[],double ed2bar[],
		int i1,int i2,int i3,double *cfac,double *cfac1,
		double *cfac2,double *cfac11,double *cfac12,
		double *cfac22)
#else	       
static void sh1467_s9fac2(ebar,ed1bar,ed2bar,i1,i2,i3,cfac,cfac1,cfac2,
			  cfac11,cfac12,cfac22)
     int i1,i2,i3;
     double ebar[],ed1bar[],ed2bar[],*cfac,*cfac1,*cfac2,*cfac11;
     double *cfac12,*cfac22;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute value and 1. and 2. derivatives of product of 
*              barycentric coordinates.
*
*
*
* INPUT      : ebar   - Array containing value of barycentric coordinates.
*              ed1bar - Array containing derivative in first parameter 
*                       direction of barycentric coordinates.
*              ed2bar - Array containing derivative in second parameter 
*                       direction of barycentric coordinates.
*              i1     - Index of first barycentric coordinate in product.
*              i2     - Index of second barycentric coordinate in product.
*              i3     - Index of third barycentric coordinate in product.
*                       
*
* OUTPUT     : cfac   - Product of barycentric coordinates.
*              cfac1  - Derivative in first parameter direction of product.
*              cfac2  - Derivative in second parameter direction of product.
*              cfac11 - 2. derivative in first parameter direction of product.
*              cfac11 - Mixed derivative of product.
*              cfac22 - 2. derivative in second parameter direction of product.
*
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
  
  *cfac = ebar[i1]*ebar[i2]*ebar[i3];
  *cfac1 = ed1bar[i1]*ebar[i2]*ebar[i3] + ebar[i1]*ed1bar[i2]*ebar[i3]
    + ebar[i1]*ebar[i2]*ed1bar[i3];
  *cfac2 = ed2bar[i1]*ebar[i2]*ebar[i3] + ebar[i1]*ed2bar[i2]*ebar[i3]
    + ebar[i1]*ebar[i2]*ed2bar[i3];
  *cfac11 = ed1bar[i1]*(ed1bar[i2]*ebar[i3] + ebar[i2]*ed1bar[i3])
    + ed1bar[i2]*(ed1bar[i1]*ebar[i3] + ebar[i1]*ed1bar[i3])
      + ed1bar[i3]*(ed1bar[i1]*ebar[i2] + ebar[i1]*ed1bar[i2]);
  *cfac12 = ed1bar[i1]*(ed2bar[i2]*ebar[i3] + ebar[i2]*ed2bar[i3])
    + ed1bar[i2]*(ed2bar[i1]*ebar[i3] + ebar[i1]*ed2bar[i3])
      + ed1bar[i3]*(ed2bar[i1]*ebar[i2] + ebar[i1]*ed2bar[i2]);
  *cfac22 = ed2bar[i1]*(ed2bar[i2]*ebar[i3] + ebar[i2]*ed2bar[i3])
    + ed2bar[i2]*(ed2bar[i1]*ebar[i3] + ebar[i1]*ed2bar[i3])
      + ed2bar[i3]*(ed2bar[i1]*ebar[i2] + ebar[i1]*ed2bar[i2]);
}



