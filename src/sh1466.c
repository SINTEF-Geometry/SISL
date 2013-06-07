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
 * $Id: sh1466.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH1466

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      sh1466(SISLCurve *ecurve[],double etwist[],int ider,double ebar[],
	     double eval[],int *jstat)
#else	 
void sh1466(ecurve,etwist,ider,ebar,eval,jstat)
     double etwist[],ebar[],eval[];
     SISLCurve *ecurve[];
     int ider,*jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given the barycentric coordinates of a point in a 3-sided
*              vertex region, evaluate the value of the ideal blending 
*              surface of the vertex region in this point.
*
*
*
* INPUT      : ecurve - Position and cross-tangent curves around the vertex
*                       region. For each edge of the region position and cross-
*                       tangent curves are given. The curves follow each other
*                       around the region and are oriented counter-clock-wise.
*                       The dimension of the array is 6.
*              etwist - Twist-vectors of the corners of the vertex region. The
*                       first element of the array is the twist in the corner
*                       before the first edge, etc. The dimension of the array
*                       is 3*kdim.
*              ider   - Number of derivatives to compute. Directions of 
*                       differentiation is that of the two first barycentric
*                       coordinates. 0 <= ider <= 2.
*              ebar   - Barycentric coordinates of the point to be evaluated.
*                       The dimension of the array is 3.
*                       
*
* OUTPUT     : eval   - Value and derivatives of ideal blending surface in the 
*                       given point. Dimension of the array is 3*(1+..+(ider+1)).
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : A functional description of the ideal surface is given as
*              a blend between three surfaces, each of which fulfill the
*              continuity requirements over two edges.
*
* REFERENCES : Gregory and Charrot : A C1 Triangular Interpolation Patch for
*                                    Computer Aided Geometric Design
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1221 - Evaluate curve at a given parameter value.
*
* WRITTEN BY : Vibeke Skytt, SI, Dec 89.
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
  double tpar1;        /* Parameter value of edge between actual corner 
			  and next corner.                                */
  double tpar2;        /* Parameter value of edge between actual corner
			  and previous corner.                            */
  double tlambi;       /* First barycentric coordinate.                   */
  double tlambj;       /* Second barycentric coordinate.                  */
  double tlambk;       /* Third barycentric coordinate.                   */
  double salpha[3];    /* Weight of actual blending surface.              */
  double sp[9];        /* Value of function interpolating two sides.      */
  double sstart[3];    /* Start-parameters of edge-curves.                */
  double send[3];      /* End-parameters of edge-curves.                  */
  double sint[3];      /* Parameter intervals of edge-curves.             */
  double spos1[27];    /* Position of edge-curve in tpar1.                */
  double sder1[27];    /* Cross tangent in tpar1.                         */
  double spos2[27];    /* Position of edge-curve in tpar2.                */
  double sder2[27];    /* Cross tangent in tpar2.                         */
  double scorn[9];     /* Position of edge-curves in actual corner of 
			  vertex region.                                  */
  double scornder1[9]; /* Tangent in corner along next edge.              */
  double scornder2[9]; /* Tangent in corner along previous edge.          */

  /* Test input.  */

  if (ider > 2) kwarn = 1;
  
  /* Initialise.  */

  kder = ider;
  knmb = kdim*(ider+1);
  for (ki=0; ki<ider; ki++) kl += ider + 1;
  
  /* Initiate output array to zero.  */

  for (kh=0; kh<kl*kdim; kh++) eval[kh] = (double)0.0;
        
  /* Get endpoints of parameter intervals of edge curves.  */

  for (ki=0; ki<3; ki++)
    {
      sstart[ki] = *(ecurve[2*ki] -> et + ecurve[2*ki] -> ik - 1);
      send[ki] = *(ecurve[2*ki] -> et + ecurve[2*ki] -> in);
      sint[ki] = send[ki] - sstart[ki];
    }

  /* Evaluate position and cross-tangent curves at points on
     the edges needed when evaluating surface.        */
  
  for (ki=0; ki<3; ki++)
    {
      kj = (ki+1) % 3;
      kk = (ki+2) % 3;

      /* Copy barycentric coordinates to local variables.  */

      tlambi = ebar[ki];      
      tlambj = ebar[kj];
      tlambk = ebar[kk];
      
      /* Find parameter values of points on edges used to evaluate
	 actual blending surface.  */

      tpar1 = ((double)1.0 - tlambj)*sstart[ki] + tlambj*send[ki];
      tpar2 = tlambk*sstart[kk] + ((double)1.0 - tlambk)*send[kk];
 
      /* Evaluate position and cross-tangent curves at first
	 found parameter value.  */

      s1221(ecurve[2*ki],kder,tpar1,&kleft,spos1+knmb*ki,&kstat);
      if (kstat < 0) goto error;
      
      s1221(ecurve[2*ki+1],kder,tpar1,&kleft,sder1+knmb*ki,&kstat);
      if (kstat < 0) goto error;

      /* Evaluate position and cross-tangent curves at second
	 found parameter value.  */
      
      s1221(ecurve[2*kk],kder,tpar2,&kleft,spos2+knmb*ki,&kstat);
      if (kstat < 0) goto error;
      
      s1221(ecurve[2*kk+1],kder,tpar2,&kleft,sder2+knmb*ki,&kstat);
      if (kstat < 0) goto error;

      /* Evaluate position and both cross-tangents at the corner nr ki.  */

      s1221(ecurve[2*ki],0,sstart[ki],&kleft,scorn+ki*kdim,&kstat);
      if (kstat < 0) goto error;

      s1221(ecurve[2*kk+1],0,send[kk],&kleft,scornder1+ki*kdim,&kstat);
      if (kstat < 0) goto error;

      s1221(ecurve[2*ki+1],0,sstart[ki],&kleft,scornder2+ki*kdim,&kstat);
      if (kstat < 0) goto error;

      /* Compute the weigth of the actual blending surface.  */

      salpha[ki] = tlambi*tlambi*((double)3.0 - (double)2.0*tlambi +
				  (double)6.0*tlambj*tlambk);
      
      /* Add the contribution of the value of this blending surface to the
	 value of the ideal surface.   */

      for (kh=0; kh<kdim; kh++)
	{
	  sp[ki*kdim+kh] = spos1[knmb*ki+kh] + tlambk*sder1[knmb*ki+kh] +
	    spos2[knmb*ki+kh] + tlambj*sder2[knmb*ki+kh] -
	      scorn[ki*kdim+kh] - tlambj*scornder1[ki*kdim+kh] 
		- tlambk*scornder2[ki*kdim+kh] -
		  tlambj*tlambk*etwist[ki*kdim+kh];
	  
	  eval[kh] += salpha[ki]*sp[ki*kdim+kh];
	}
    }
  
  if (ider >= 1)
    {
      /* Compute first derivatives of the Gregory Charrot function. */

      double tl1,tl2,tl3;            /* Barycentric coordinates.  */
      double sd1alpha[3],sd2alpha[3];  /* 1. derivative of weight functions. */
      double sd1p[9],sd2p[9];        /* 1. derivative of blending functions. */
      
      /* Copy barycentric coordinates to local variables.  */

      tl1 = ebar[0];
      tl2 = ebar[1];
      tl3 = ebar[2];
      
      /* Compute the 1. derivatives of the weight functions.  */

      sd1alpha[0] = (double)6.0*tl1*((double)1.0 - tl1 - tl1*tl2 + 
				  (double)2.0*tl2*tl3);
      sd2alpha[0] = (double)6.0*tl1*tl1*(tl3 - tl2);

      sd1alpha[1] = (double)6.0*tl2*tl2*(tl3 - tl1);
      sd2alpha[1] = (double)6.0*tl2*((double)1.0 - tl2 - tl1*tl2 + 
				  (double)2.0*tl1*tl3);

      sd1alpha[2] = (double)6.0*tl3*(-(double)1.0 + tl3 + tl2*tl3 - 
				  (double)2.0*tl1*tl2);
      sd2alpha[2] = (double)6.0*tl3*(-(double)1.0 + tl3 + tl1*tl3 - 
				  (double)2.0*tl1*tl2);
      
      /* Compute 1. derivatives of the functions which blends two sides
	 of the region.  */

      for (kh=0; kh<kdim; kh++)
	{
	  sd1p[kh] = (spos2[kdim+kh] + tl2*sder2[kdim+kh])*sint[2]
	    - sder1[kh] + scornder2[kh] + tl2*etwist[kh];
	  sd2p[kh] = (spos2[kdim+kh] + tl2*sder2[kdim+kh])*sint[2] 
	    + (spos1[kdim+kh] + tl3*sder1[kdim+kh])*sint[0]
	      + sder2[kh] - sder1[kh] - scornder1[kh]
		+ scornder2[kh] + (tl2 - tl3)*etwist[kh];

	  sd1p[kdim+kh] = -sder2[knmb+kh] + sder1[knmb+kh]
	    - (spos2[knmb+kdim+kh] + tl3*sder2[knmb+kdim+kh])*sint[0] 
	      - (spos1[knmb+kdim+kh] + tl1*sder1[knmb+kdim+kh])*sint[1] 
		+ scornder1[kdim+kh] - scornder2[kdim+kh]
		  + (tl1 - tl3)*etwist[kdim+kh];
	  sd2p[kdim+kh] = - sder2[knmb+kh]
	    - (spos1[knmb+kdim+kh] + tl1*sder1[knmb+kdim+kh])*sint[1]
	      + scornder1[kdim+kh] + tl1*etwist[kdim+kh];

	  sd1p[2*kdim+kh] = sder2[2*knmb+kh]
	    + (spos1[2*knmb+kdim+kh] + tl2*sder1[2*knmb+kdim+kh])*sint[2]
	      - scornder1[2*kdim+kh] - tl2*etwist[2*kdim+kh];
	  sd2p[2*kdim+kh] = sder1[2*knmb+kh]
	     - (spos2[2*knmb+kdim+kh] + tl1*sder2[2*knmb+kdim+kh])*sint[1]
	      - scornder2[2*kdim+kh] - tl1*etwist[2*kdim+kh];
	  
	  /* Compute the first derivative of the Gregory Charrot function. */

	  for (ki=0; ki<3; ki++)
	    {
	      eval[kdim+kh] += sd1alpha[ki]*sp[ki*kdim+kh]
		+ salpha[ki]*sd1p[ki*kdim+kh];

	      eval[2*kdim+kh] += sd2alpha[ki]*sp[ki*kdim+kh]
		+ salpha[ki]*sd2p[ki*kdim+kh];
	    }
	}
      
      if (ider >= 2)
	{
	  double sd11alpha[3],sd12alpha[3],sd22alpha[3];  /* 2. derivatives of
                                                             blending patches. */
	  double sd11p[9],sd12p[9],sd22p[9];  /* 2. derivatives of weight function. */
	  
	  /* Compute second derivatives of the Gregory Charrot function. */

	  /* Compute the 2. derivatives of the weight functions.  */

	  sd11alpha[0] = (double)6.0 - (double)12.0*tl1 
	    - (double)24.0*tl1*tl2 + (double)12.0*tl2*tl3;
	  sd12alpha[0] = tl1*((double)12.0 - (double)18.0*tl1
			      - (double)24.0*tl2);
	  sd22alpha[0] = -(double)12.0*tl1*tl1;
	  
	  sd11alpha[1] = -(double)12.0*tl2*tl2; 
	  sd12alpha[1] = tl2*((double)12.0 - (double)18.0*tl2
			      - (double)24.0*tl1);
	  sd22alpha[1] = (double)6.0 - (double)12.0*tl2 
	    - (double)24.0*tl1*tl2 + (double)12.0*tl1*tl3;

	  sd11alpha[2] = (double)6.0 - (double)12.0*tl3 
	    - (double)24.0*tl2*tl3 + (double)12.0*tl1*tl2;
	  sd12alpha[2] = (double)6.0 + (double)6.0*tl3*tl3
	    + (double)12.0*(-tl3 + tl1*tl2 - tl1*tl3 - tl2*tl3);
	  sd22alpha[2] = (double)6.0 - (double)12.0*tl3 
	    - (double)24.0*tl1*tl3 + (double)12.0*tl1*tl2;

	  /* Compute 2. derivatives of the functions which blends two sides
	     of the region.  */

	  for (kh=0; kh<kdim; kh++)
	    {
	      sd11p[kh] = (spos2[2*kdim+kh]+ tl2*sder2[2*kdim+kh])*sint[2]*sint[2];
	      sd12p[kh] = ((spos2[2*kdim+kh] + tl2*sder2[2*kdim+kh])*sint[2]
			   + sder2[kdim+kh]*sint[2] - sder1[kdim+kh])*sint[0] 
			      + etwist[kh];
	      sd22p[kh] = ((spos2[2*kdim+kh] + tl2*sder2[2*kdim+kh])*sint[2]
			   + (double)2.0*sder2[kdim+kh])*sint[2]
			     + ((spos1[2*kdim+kh] + tl3*sder1[2*kdim+kh])*sint[0]
				- (double)2.0*sder1[kdim+kh])*sint[0]
				  + (double)2.0*etwist[kh];
	      
	      sd11p[kdim+kh] = ((spos2[knmb+2*kdim+kh] 
				 + tl3*sder2[knmb+2*kdim+kh])*sint[0]
				+ (double)2.0*sder2[knmb+kdim+kh])*sint[0]
				  + ((spos1[knmb+2*kdim+kh] 
				      + tl1*sder1[knmb+2*kdim+kh])*sint[1]
				     - (double)2.0*sder1[knmb+kdim+kh])*sint[1]
				       + (double)2.0*etwist[kdim+kh];
	      sd12p[kdim+kh] = ((spos1[knmb+2*kdim+kh] 
				 + tl1*sder1[knmb+2*kdim+kh])*sint[1]
				- sder1[knmb+kdim+kh])*sint[1]
				  + sder2[knmb+kdim+kh]*sint[0]
				  + etwist[kdim+kh];
	      sd22p[kdim+kh] = (spos1[knmb+2*kdim+kh]
				+ tl1*sder1[knmb+2*kdim+kh])*sint[1]*sint[1];
	      
	      sd11p[2*kdim+kh] = (spos1[2*knmb+2*kdim+kh]
				  + tl2*sder1[2*knmb+2*kdim+kh])*sint[2]*sint[2];
	      sd12p[2*kdim+kh] = -sder2[2*knmb+kdim+kh]*sint[1] 
		+ sder1[2*knmb+kdim+kh]*sint[2] - etwist[2*kdim+kh];
	      sd22p[2*kdim+kh] = (spos2[2*knmb+2*kdim+kh]
				+ tl1*sder2[2*knmb+2*kdim+kh])*sint[1]*sint[1];

	      /* Compute the 2. derivative of the Gregory Charrot function. */

	      for (ki=0; ki<3; ki++)
		{
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


			      
	     


