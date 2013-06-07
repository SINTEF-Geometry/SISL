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
 * $Id: s1301.c,v 1.2 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1301

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1301(double areler,double angle,int idim,SISLCurve **pc,int *jstat)
#else
void s1301(areler,angle,idim,pc,jstat)
     double areler;
     double angle;
     int    idim;
     SISLCurve  **pc;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To create a B-spline curve approximating a circular
*              arc defined by radius=1 and the input angle. The
*              arc will start in the point (1,0) in the xy-plane
*              and the rotational direction is positive. If the opening
*              angle is less than TWOPI then a k-regular basis is made,
*              if the opening angle describes a full circle then a cyclic
*              basis with double knots are made.
*                
*
*             
*
* INPUT      : areler - The relative error the circular arc is to be
*                       produced within.
*              angle  - The angle defining how much of the circle is
*                       to be made (in radians)
*              idim   - The dimension of the curve to be produced
*
*
*
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              pc     - Pointer to the curve produced
*
* METHOD     : First the number of polynomial segment to be produced
*              are found by using an error estimat formula
*              The tangent lengths of the interpolation points are
*              determined. Finally the actual knot vector and and
*              vertices are made.
*
*
* REFERENCES :
*
*-                                      
* CALLS      : pow, cos, sin. sqrt, fabs, newCurve
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, May 1988
* Reviced by : Tor Dokken, SI, Oslo, Norway, April 1992 Introduction of cyclic basis.
*
*********************************************************************
*/
{
  int    kant;           /* Number of segments                     */
  int    kl;             /* Pointer into array                     */
  int    kj;             /* Pointer into array                     */
  int    ki;             /* Control variable in loop               */
  int    kn;             /* Number of vertice                      */
  int    kk = 4;         /* Polynomial order of curve              */
  int    kpos = 1;       /* Position of error                      */
  int    kstart,kstop;   /* Loop control variable                  */
  double tangle;         /* Absolute value of angle                */
  double talfa;          /* Variable used for different angles     */
  double tcos,tsin;      /*Varable used for cos and sin of talfa   */
  double *st = SISL_NULL;     /* Pointer to knot vecttor                */
  double *scoef = SISL_NULL;  /* Pointer to vertex array                */
  double ta,tb,tc;       /* Temporary variables                    */
  double tl;             /* Length of tangents                     */
  double tconst;         /* Internal constant                      */
  
  /* Check that the relative error wanted is positive */
  if (areler<=(double)0.0) goto err120;
  
  /* Check that the dimesnion given is greater or equal to two */
  if (idim<2) goto err103; 
  
  /* Make first the curve as if the angle is only positive */
  
  tangle = fabs(angle);
  
  /*  If the absolute value of the rotation angle is greater than 2*PI
   *  then make the angle to 2*PI and update angle to have absolute
   *  value 2*PI 
   */
  
  
  if (tangle >= TWOPI)
    {
      tangle = TWOPI;
    }
  
  /*  Estimat the opening angle of the segments based on the error 
   *   formula. */                                                  
  ta = (double)1.0/(double)6.0;
  talfa = PI*pow(areler,ta)/((double)0.4879);
  
  /*   Find number of segments actually to be produced */
  
  kant = (int) ((tangle/talfa)+1);
  
  /*   Calculate actual angle to be used
   *   ---------------------------------
   */
  
  talfa = tangle/kant;
  tcos = cos(talfa);
  tsin = sin(talfa);
  
  /*  Calculate length of tangents   
   *  tconst = (3-2sqrt(2))**1/3 + (3+2sqrt(2))**1/3 - 0.5 
   */
  
  tconst = (double)1.85530139760811990992528773586425;
  ta     = (double)0.6*tconst - (double)0.9*tcos;
  tb     = ((double)0.4*tconst+(double)1.8)*tsin;
  tc     = ((double)0.4*tconst+(double)1.0)*tcos -
    (double)0.4*tconst - (double)1.0;
  tl     = (-tb+sqrt(tb*tb-(double)4.0*ta*tc))/((double)2.0*ta);
  
  /*   Allocate space for vertices and knot vector; */
  
  kn = 2*kant + 2;
  scoef = newarray(kn*idim,DOUBLE);
  st    = newarray(kn+kk,DOUBLE);                       
  if (scoef == SISL_NULL || st == SISL_NULL) goto err101;
  
  
  /*   Calculate vertices and make knots. */
  
  kl = 0;
  kj = 2;
  
  /*   Initate the whole vertex array to zero. */
  
  kstop = idim*kn;
  for (ki=0;ki<kstop;ki++)
    {
      scoef[ki] = (double)0.0;
    }
  
  /* If the opening angle is less than TWOPI the two first and two lastrcle is
   *                       to be made (in radians)
   vertices and knots are to be made in a special way to give a k-regular basis */
  
  if (tangle >= TWOPI)
    {
      st[0] = (double)-1;  
      st[1] = (double)-1;
      kstart = 0;
      kstop = kant + 1;
      
      /* Next knot is to be made at index 2 */
      
      kj=2;
      
    }
  
  else
    {
      
      
      /* Make the two first vertices, and four first knots */
      
      tcos = (double)1.0;
      tsin = (double)0.0;
      scoef[0] = (double)1.0;
      scoef[1] = (double)0.0;
      scoef[idim] = (double)1.0;
      scoef[idim+1] = tl;
      st[0] = (double)0.0;  
      st[1] = (double)0.0;
      st[2] = (double)0.0;
      st[3] = (double)0.0;
      kstart = 1;
      kstop = kant;
      
      /* Next knot is to be made at index 2 */
      
      kj = 4;
      
    }  
  /* Make internal vertices, remember that only double knots exist */
  
  for (ki=kstart;ki<kstop;ki++)
    {
      kl = 2*ki*idim;
      talfa = ki*tangle/kant;
      
      /* To make a stable calculation of the sin and cos values, we use
       * arguments only in the interval [0,PI] 
       */
      
      
      tcos = cos(talfa);
      tsin = sin(talfa);
      
      /*
	if (talfa<=PIHALF)
        {
        tcos = cos(talfa);
        tsin = sin(talfa);
        }
	else if ( talfa<=PI)
        {
        talfa = PI - talfa;
        tcos = -cos(talfa);
        tsin = sin(talfa);
        }
	else if (talfa<=THREEPIHALF)
        {
        talfa = talfa - PI;
        tcos = -cos(talfa);
        tsin = -sin(talfa);
        }
	else
        {
        talfa = TWOPI - talfa;
        tcos = cos(talfa);
        tsin = -sin(talfa);
        }                      
	*/
      scoef[kl]        = tcos + tl*tsin;
      scoef[kl+1]      = tsin - tl*tcos;
      scoef[kl+idim]   = tcos - tl*tsin;
      scoef[kl+idim+1] = tsin + tl*tcos;
      st[kj++] = ki;
      st[kj++] = ki;
    }                     
  
  if (tangle >= TWOPI)
    {
      st[kn+2] = kant+1;
      st[kn+3] = kant+1;
      
      
      /*  Make sure that the two first and to last vertices correspond */
      
      ki = (kn-2)*idim;
      scoef[ki] = scoef[0];
      scoef[ki+1] = scoef[1];
      scoef[ki+2] = scoef[2];
      scoef[ki+3] = scoef[3];
      
      
    }
  else
    {
      
      /*   Make two last vertices and four last knots  */
      
      tcos = cos(tangle);
      tsin = sin(tangle);
      kl = 2*kant*idim;;
      scoef[kl]        = tcos + tl*tsin;
      scoef[kl+1]      = tsin - tl*tcos;
      scoef[kl+idim]   = tcos;
      scoef[kl+idim+1] = tsin;
      st[kn]   = kant;
      st[kn+1] = kant;
      st[kn+2] = kant;
      st[kn+3] = kant;
    }
  
  
  /* If negative rotation angle make opposite sign of y-variable */
  
  if (angle<(double)0.0)
    {
      kl = 1;
      for (ki=0;ki<kn;ki++)
	{
	  scoef[kl] = -scoef[kl];
	  kl+=idim;
	}
    }
  
  /* Make curve, copy input arrays */
  
  *pc = newCurve(kn,kk,st,scoef,1,idim,1);
  if (*pc == SISL_NULL) goto err101;

  if (tangle >= TWOPI)
    (*pc)->cuopen = SISL_CRV_PERIODIC;
  
  /* Normal circle segment made */
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1301",*jstat,kpos);
  goto out;
  
  /* Error in input, negative relative tolerance given */
  
 err120: *jstat = -120;
  s6err("s1301",*jstat,kpos);
  goto out;
  
  /* Error in input, dimension less than 2 given */
  
 err103: *jstat = -103;
  s6err("s1301",*jstat,kpos);
  goto out;
  
  /* Free allocated arrays */
 out:
  
  if (st != SISL_NULL)    freearray(st);
  if (scoef != SISL_NULL) freearray(scoef);
  
  return;
  
}
    
