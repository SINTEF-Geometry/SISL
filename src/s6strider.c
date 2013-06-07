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
 * $Id: s6strider.c,v 1.3 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6STRIDER

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6strider(double eder[],int idim,int ider,double gder[],int *jstat)
#else
void s6strider(eder,idim,ider,gder,jstat)
     double eder[];
     int    idim;
     int    ider;
     double gder[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the value and ider*ider derivatives of 
*              a rational B-spline surface.
*
* INPUT      : eder    - Double array of dimenson [(ider+1)*(ider+2)*(idim+1)/2]
*                        containing the position and the derivative vectors
*                        of the homogeneous surface at the point with parameter value
*                        (epar[0],epar[1]).
*                        (idim+1 is the number of components of each B-spline
*                        coefficient, i.e. the dimension of the homogemous
*                        space in which the surface lies.)
*                        These vectors are stored in the following order:
*                        First the idim+1 components of the position vector,
*                        then the idim+1 components of the D(1,0) vector,
*                        then the idim+1 components of the D(0,1) vector,
*                        then the idim+1 components of the D(2,0) vector,
*                        followed by D(1,1), D(0,2)
*                        and so on up to the idim+1 components of the D(0,ider).
*              idim    - The dimension of the non homogenous space
*              ider    - The number of input derivatives with respect
*                        to both parameter directions.
*
*
* OUTPUT     : jstat   - Status message
*                                        >0      : Warning
*                                        =0      : ok
*                                        <0      : Error
*              gder    - Double array of dimension [(ider+1)*(ider+2)*idim/2]
*                        containing the position and the derivative vectors
*                        of the surface at the point with parameter value
*                        (epar[0],epar[1]).
*                        (idim is the number of components of each B-spline
*                        coefficient, i.e. the dimension of the Euclidean
*                        space in which the surface lies.)
*                        These vectors are stored in the following order:
*                        First the idim components of the position vector,
*                        then the idim components of the D(1,0) vector,
*                        then the idim components of the D(0,1) vector,
*                        then the idim components of the D(2,0) vector,
*                        followed by D(1,1), D(0,2)
*                        and so on up to the idim components of the D(0,ider).
*
*
* METHOD     :  The surface P(u,v) can be written as the quotient
*               P(u,v) = T(u,v) / w(u,v) where T and w are ordinary splines.
*               The dimensions of T and w are idim and 1
*               respectively. The array eder contains position
*               and derivatives of the idim+1 dimensional surface
*               (T(u,v),w(u,v)).
*
*               Now, since wP = T, we find, by the Leibnitz formula,
*
*      k   l
*                  k!       l!     (k-i,l-j) (i,j)         (k,l)
*     sum sum   -------- -------- w         P         =   T       .
*               i!(k-i)! j!(l-j)!
*     i=0 j=0
*
*               Therefore
*               
*
*              --            k   l                                     --
*      (k,l)   |   (k,l)                k!       l!     (k-i,l-j) (i,j) |    
*     P      = |  T      -  sum sum  -------- -------- w         P      | / w .
*              |                     i!(k-i)! j!(l-j)!                  |
*              --           i=0 j=0                                    --
*                               i+j<k+l
*
*               This formula is applied recursively to evaluate P's derivatives.
*
*                                                          MF.
*
*
*
* CALLS      :
*
* WRITTEN BY : Michael Floater, SI, 3.9.92.
*                Essentially the same as s6sratder
*                except that we work with triangular matrices
*                ((0,0), (1,0), (0,1), (2,0), (1,1), ...)
*                instead of rectangular ones.
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              Error message corrected
*
*********************************************************************
*/
{
  int kpos=0;          /* Position of error.                     */
  double w0;           /* The denominator.                       */
  int ki;              /* Count through dimensions.              */
  int idu;             /* Count through derivatives in u.        */
  int idv;             /* Count through derivatives in v.        */
  int *binom=SISL_NULL;     /* Array for binomial coefficients.       */
  int *binomu=SISL_NULL;    /* Pointer to binomial coefficients in u. */
  int *binomv=SISL_NULL;    /* Pointer to binomial coefficients in v. */
  double *sum1=SISL_NULL;   /* Leibnitz expansion in u                */
  double *sum2=SISL_NULL;   /* Leibnitz expansion in u and v.         */
  double sumdum1[4];   /* Fixed space for sum1.                  */
  double sumdum2[4];   /* Fixed space for sum2.                  */
  int idimp1;          /* idim + 1.                              */
  int iw;              /* Pointer to a weight.                   */
  int igder;           /* Pointer to already calculated derivs.  */
  int i,iu,iv,j,k;     /* Counters.                              */
  int iderp1;          /* ider + 1.                              */
  int igrow;           /* (ider+1) * idim.                       */  
  int iwrow;           /* (ider+1) * idimp1.                     */  
  int iutemp,ivtemp;   /* Used to find next weight in the sum.   */
  int tot,temp1;       /* Temporary variables.                   */
  int bidum[10];       /* Array for storing binomial coeffs.     */
  double temp;         /* Temporary multiple.                    */
  
  if (ider<0) goto err178;
  if (idim<1) goto err102;
  
  *jstat = 0;

  /* Find denominator. */ 
  
  w0 = eder[idim];
  if (DEQUAL(w0,DZERO)) w0 = (double)1.0;

  /* If we're only asked for position, we'll do it
     now and exit for the sake of speed. */

  if(ider == 0)
  {
    for(ki=0; ki<idim; ki++)
      gder[ki] = eder[ki] / w0;

    goto out;
  }

  /* Set up some constants. */

  idimp1  = idim + 1;
  iderp1 = ider + 1;
  igrow   = iderp1 * idim;
  iwrow   = igrow + iderp1;  /* = iderp1 * idimp1 */

  /* Set up  binomial coefficients.
     Use new array only when ider > 3. */

  if (ider > 3)
  { 
    binom = newarray((iderp1*(iderp1+1)) >> 1, INT);
    if(binom == SISL_NULL) goto err179;
  }
  else
  { 
    binom = bidum;
  }

  for(j=0,k=0; j<=ider; j++,k+=j)
  {
      /* Calculate the new row of binomial coefficients. */
  
      binom[k] = 1;
  
      for(i=k+1; i<k+j; i++)
      {
          binom[i] = binom[i-j-1] + binom[i-j];
      }

      binom[k+j] = 1;
  }

  /* Set up space for sum1 and sum2 if necessary.
     Use new arrays only when idim > 4. */

  if (idim > 4)
  { 
    sum1 = newarray(idim, DOUBLE);
    if(sum1 == SISL_NULL) goto err179;
    sum2 = newarray(idim, DOUBLE);
    if(sum2 == SISL_NULL) goto err179;
  }
  else
  { 
    sum1=sumdum1;
    sum2=sumdum2;
  }
  
  /* Loop through derivatives in u and v. */

  for(idv=0,binomv=binom; idv<=ider; idv++,binomv+=idv)
  {
    for(idu=0,binomu=binom; idu<=ider-idv; idu++,binomu+=idu)
    {
      if(idu == 0 && idv == 0)
      {
          /* Position is a special case. */
    
          for(ki=0; ki<idim; ki++)
            gder[ki] = eder[ki] / w0;
      }
      else
      {
          /* Calculate indices in eder and gder. */
    
          tot = idu + idv;
          temp1 = ((tot * (tot+1)) >> 1) + idv;
    
          j = temp1 * idim;
          k = j + temp1;

          /* Calculating each coefficient of the (idu,idv)'th
	     derivative of the rational surface (in gder).
        
  	     This requires calculating the Liebnitz sum from
  	     the subarray of gder (0,..,idu, 0,...,idv) and
             the subarray of eder (0,..,idu, 0,...,idv). */

          /* Calculate the Leibnitz sum. */

          for(ki=0; ki<idim; ki++)
            sum2[ki] = (double)0.0;        

          for(iv=0; iv<=idv; iv++)
          {
            for(ki=0; ki<idim; ki++)
               sum1[ki] = (double)0.0;	               
            ivtemp = idv-iv;

            for(iu=0; iu<=idu; iu++)
            {
                tot = iu + iv;
                temp1 = ((tot * (tot+1)) >> 1) + iv;

	        igder = temp1 * idim;
                iutemp = idu-iu;

                tot = iutemp + ivtemp;
	        temp1 = ((tot * (tot+1)) >> 1) + ivtemp;

                iw   = temp1 * idimp1 + idim;

     	      /* Add the next Leibnitz term unless we
       		 have reached the last one (the unknown). */
  
                if(iu<idu || iv<idv)
  	        {
  	       	  /* If iu=0 or iu=idu, the u binomial
  	       	     coefficient is 1 so don't multiply. */
  
  	            if(iu>0 && iu<idu)
  	       	    {
  		      temp = (double)binomu[iu] * eder[iw];
                      for(ki=0; ki<idim; ki++,igder++)
  	       	         sum1[ki] += temp * gder[igder];
  		     }
  		     else
                       for(ki=0; ki<idim; ki++,igder++)
  		         sum1[ki] += eder[iw] * gder[igder];
                }
            }
  
  	    /* If iv=0 or iv=idv, the v binomial
  	       coefficient is 1 so don't multiply. */
  
  	    if(iv>0 && iv<idv)
              for(ki=0; ki<idim; ki++)
  	          sum2[ki] += (double)binomv[iv] * sum1[ki];
  	    else
              for(ki=0; ki<idim; ki++)
		 sum2[ki] += sum1[ki];		    
          }
          for(ki=0; ki<idim; ki++,j++,k++)
            gder[j] = (eder[k] - sum2[ki]) / w0;
      }
    }
  }  

  /* Free arrays. */

  if (ider > 3 && binom != SISL_NULL)
     freearray(binom);
  
  if (idim > 4)
  { 
     if(sum1 != SISL_NULL) freearray(sum1);
     if(sum2 != SISL_NULL) freearray(sum2);
  }

  /* Done. */

  goto out;

  /* idim less than 1. */

  err102: 
    *jstat = -102;
    s6err("s6strider",*jstat,kpos);
    goto out;

  /* Derivative negative */

  err178: 
    *jstat = -178;
    s6err("s6strider",*jstat,kpos);
    goto out;

  /* Not enough memory */

  err179: 
    *jstat = -179;
    s6err("s6strider",*jstat,kpos);
    goto out;

  out:
    return;
}

