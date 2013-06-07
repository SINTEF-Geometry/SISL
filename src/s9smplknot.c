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
 * $Id: s9smplknot.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S9SIMPLE_KNOT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s9simple_knot(SISLSurf* surf, int idiv, double epar[], 
		   int *fixflag, int *jstat)
#else
void s9simple_knot(surf, idiv, epar, fixflag, jstat)
     SISLSurf* surf;
     int idiv; 
     double epar[];
     int *fixflag;
     int *jstat;
#endif
/*
***************************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To check if the spline surface surf has more than one inner 
*              breakpoint in direction idiv (could be both directions).
*
*
*
* INPUT      : surf      - Pointer to the spline surface
*              idiv      - direction to check, possible values:
*			                1 first parameter direction
*			                2 second parameter direction
*			                3 both parameter directions
*
*
*
* OUTPUT     : jstat     - status messages  
*                                         = 1      : at most one inner 
*						     breakpoint found in 
*						     direction(s) idiv
*                                         = 0      : more than one inner 
*						     breakpoint found in 
*						     direction idiv, (in one
*						     of them if idiv=3)
*                                         < 0      : error
*              epar[i]   - parameter value of simple inner knot in diretion i 
*			   (i=0 first direc, i=1 second direc.),
*			   midpoint if no inner simple knot found in direc i,
*                          not set if more than one inner breakpt in direc i or
*			   direc i not indicated by idiv.
*			   
*              fixflag   - indicates simple inner breakpt found in 
*			   any of the direcs given by idiv:
*			   1 simple knot found in first direc,
*			   2 simple knot found in second direc,
*			   3 simple knot found in both direc.
*                          0 no or more than one simple breakpt found in all
*                            the directions indicated by idiv.
*
* ASSUMPTIONS: k-tuple knots at the ends of the knot vector
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6err      - Treating error situation.
*              s6knotmult - Counting multiplicity of a knot
*
* WRITTEN BY : Kyrre Strom, SI, 93-01.
*
*
****************************************************************************
*/
{
  int k1,k2,kstat,mult;

  k1 = k2 = *fixflag = 0;

  if ( idiv < 1 || idiv > 3 ) goto err202;
  if (idiv == 1 || idiv == 3) /* Check in first parameter direction */
    {
      if ( surf->in1 == surf->ik1 )
	{
	  epar[0] = (surf->et1[0] + surf->et1[surf->in1+surf->ik1-1])/2.0;
	  k1 = 1;
	}
      else 
	{
	  int left = surf->ik1;
	  mult = s6knotmult(surf->et1,surf->ik1,surf->in1, &left,
			    surf->et1[surf->ik1],&kstat);
	  if (kstat < 0 ) goto error;
	  if ( surf->ik1+mult == surf->in1 )
	    {
	      epar[0] = surf->et1[surf->ik1];
	      k1 = 1;
	      *fixflag += 1;
	    }
	}
    }

  if (idiv == 2 || idiv == 3)
    {
      if ( surf->in2 == surf->ik2 )
	{
	  epar[1] = (surf->et2[0] + surf->et2[surf->in2+surf->ik2-1])/2.0;
	  k1 += 2;
	}
      else 
	{
	  int left = surf->ik2;
	  mult = s6knotmult(surf->et2,surf->ik2,surf->in2, &left,
			    surf->et2[surf->ik2],&kstat);
	  if (kstat < 0 ) goto error;
	  if ( surf->ik2+mult == surf->in2 )
	    {
	      epar[1] = surf->et2[surf->ik2];
	      k1 += 2;
	      *fixflag += 2 ;
	    }
	}
    }
  
  *jstat = ((idiv == k1 && (*fixflag)) ? 1 : 0);
  goto out;

 error : *jstat = kstat;
         s6err("s9simple_knot",*jstat,0);
  	 goto out;

 err202 : *jstat = -202;
         s6err("s9simple_knot",*jstat,0);

 out:  return;         
 }
