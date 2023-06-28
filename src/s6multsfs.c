/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl-copyright.h"


#define S6MULTSFS

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
s6multsfs(double *c1, int order11,int order12,
		double* c2, int order21,int order22,
                double *Pascal,
		double *newsurf, int *order_newsurf1,int *order_newsurf2)
#else
   void s6multsfs(c1, order11, order12, c2, order21, order22, Pascal,
		  newsurf, order_newsurf1, order_newsurf2)
      double *c1;
      int order11;
      int order12;
      double *c2;
      int order21;
      int order22;
      double *Pascal;
      double *newsurf;
      int *order_newsurf1;
      int *order_newsurf2;
#endif      
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Multiply two Bezier surfaces of dimension 1.
*
*
*
* INPUT      : c1      - Coefficients of first Bezier surface.
*              order11 - Order of first surface in first parameter dir.
*              order12 - Order of first surface in second parameter dir.
*              c2      - Coefficients of second Bezier surface
*              order21 - Order of second surface in first parameter dir.
*              order22 - Order of second surface in first parameter dir.
*              Pascal  - Factors used in surface multiplication.
*
*
*
* OUTPUT     : newsurf - Coefficients of product surface.
*              order_newsurf1 - Number of coefficients of product surface
*                               in first parameter direction.
*              order_newsurf2 - Number of coefficients of product surface
*                               in second parameter direction.
*
*
* REFERENCES :
*
* NOTE       : Factors necessary for the multiplication must be
*              precomputed. If the multiplication routine is to be
*              calles severatl times, it speeds up the computations
*              if this factors are called only once.
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SINTEF, 1993. 
* COMMENTED BY : Vibeke Skytt, SINTEF, 06.94.
*
*********************************************************************
*/                                     
{
  int p1,p2,r01,r02,r1,r2,h,h1;
  int kgrad11=order11-1;
  int kgrad12=order12-1;
  int kgrad21=order21-1;
  int kgrad22=order22-1;
  int kgrad1 = kgrad11 + kgrad21;
  int kgrad2 = kgrad12 + kgrad22;
  int kstop2=order12+order22-1;
  int kstop1=order11+order21-1;
  double *psl_kgrad11=Pascal+kgrad11*(kgrad11+1)/2;
  double *psl_kgrad12=Pascal+kgrad12*(kgrad12+1)/2;
  double *psl_kgrad21=Pascal+kgrad21*(kgrad21+1)/2;
  double *psl_kgrad22=Pascal+kgrad22*(kgrad22+1)/2;
  double *psl_kgrad1 =Pascal+kgrad1 *(kgrad1 +1)/2;
  double *psl_kgrad2 =Pascal+kgrad2 *(kgrad2 +1)/2;
  double *qsc1,*qsc2;
  double tsum, sumi;
  double *temp;
  double tdiv, t2;
  int kstop3, kstop4;
  double scratch[190];
  double *binom = NULL;

  if (kstop1*order11 < 190)
      binom = scratch;
  else
  {
      binom = newarray(kstop1*order11, double);
  }
      

      
  for (p1=0, h=0; p1<kstop1; p1++)
  {
      kstop3 = MIN(p1,kgrad11);
      for (r1=MAX(0,p1-kgrad21); r1<=kstop3; r1++)
	  binom[h++] = psl_kgrad11[r1]*psl_kgrad21[p1-r1];
  }
      	  
  

  for(p2=0,temp=newsurf;p2<kstop2;p2++)
  {
      kstop4 = MIN(p2,kgrad12);
      r02 = MAX(0,p2-kgrad22);
      for(p1=0, h=h1=0; p1<kstop1; p1++, temp++)
      {
	  tdiv = psl_kgrad1[p1]*psl_kgrad2[p2];
	  kstop3 = MIN(p1,kgrad11);
	  r01 = MAX(0,p1-kgrad21);
	  for(r2=r02, tsum=(double)0; r2<=kstop4; r2++)
	  {
	      t2 = psl_kgrad12[r2]*psl_kgrad22[p2-r2];
	      for(r1=r01, h=h1, qsc1=c1+r2*order11, qsc2=c2+(p2-r2)*order21, sumi=0.0;
		  r1<=kstop3; r1++, h++)
		  sumi += binom[h]*qsc1[r1]*qsc2[p1-r1];
	      tsum += t2*sumi;
	  }
	  h1 += (kstop3+1-r01);
	  tsum /= tdiv;
	 
	  *temp = tsum;
      }
  }

  *order_newsurf1=kstop1;
  *order_newsurf2=kstop2;
 
  if (binom && binom != scratch)
      freearray(binom);
}

