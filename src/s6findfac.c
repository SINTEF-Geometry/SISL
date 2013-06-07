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
 * $Id: s6findfac.c,v 1.4 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S6FINDFAC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6findfac(double evecu[],double evecv[],double evecw[],double etang[],
               int idim,int isign,double *coef1,double *coef2,double *coef3,int *jstat)

#else
void s6findfac(evecu,evecv,evecw,etang,idim,isign,coef1,coef2,coef3,jstat)
     double evecu[];
     double evecv[];
     double evecw[];
     double etang[];
     int    idim;
     int    isign;
     double *coef1;
     double *coef2;
     double *coef3;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given four vectors, evecu, evecv, evecw and etang, find 
*              the factors, coef1, coef2 and coef3, such that the vector
*              coef1*evecu + coef2*evecv + coef3*evecw = isign*etang.
*              
*
* INPUT      : evecu      - First vector.
*              evecv      - Second vector.
*              evecw      - Third vector.
*              etang      - Vector to approximate.
*              idim       - Dimension of geometry space.
*              isign      - Sign with wich etang is to be multiplied.
*
*
* OUTPUT     : coef1      - First factor.
*              coef2      - Second factor.
*              coef3      - Third factor.
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Minimize the square of the expression 
*                   dist(coef1*evecu+coef2*evecv,isign*etang)
*              over coef1 and coef2.
*              The expression is differentiated and set equal to
*              zero. Then this equation system of 2 equations
*              with two unknowns is solved. 
*              If the three vectors evecu, evecv and evecw span
*              3D, |coef3| = |d|/|evecw| where 
*              d = isign*etang - coef1*evecu - coef2*evecv
*
*********************************************************************
*/
{

  int kstat = 0;           /* Status variable.                    */
  int ki;                  /* Counter.                            */
  double tdotuu;           /* Scalar product of evecu and evecu.  */
  double tdotuv;           /* Scalar product of evecu and evecv.  */
  double tdotutang;        /* Scalar product of evecu and etang.  */
  double tdotvv;           /* Scalar product of evecv and evecv.  */
  double tdotvtang;        /* Scalar product of evecv and etang.  */
  double tdiv;             /* Determinant of equation system.     */
  double sdum[3];          /* Help vector.     */

  *jstat = 0;
  
  /* Test input.  */

  /* if (idim != 3) goto err104; */
  
  /* Set output to zero. */

  *coef1 = (double)0.0;
  *coef2 = (double)0.0;
  
  /* Compute coefficients of equation system.  */

  tdotuu = s6scpr(evecu,evecu,idim);
  tdotuv = s6scpr(evecu,evecv,idim);
  tdotutang = (double)isign*s6scpr(evecu,etang,idim);
  tdotvv = s6scpr(evecv,evecv,idim);
  tdotvtang = (double)isign*s6scpr(evecv,etang,idim);

  tdiv = tdotuv*tdotuv - tdotuu*tdotvv;
  if (DEQUAL(tdiv,DZERO))
    {
      if (DEQUAL(tdotuu,DZERO) && DEQUAL(tdotvv,DZERO));
      else if (DEQUAL(tdotuu,DZERO))
	  *coef2 = s6length(etang,idim,&kstat)/sqrt(tdotvv);
      else
	*coef1 = s6length(etang,idim,&kstat)/sqrt(tdotuu);
      goto out;
    }
  
  /* Compute the first two output factors.  */

  *coef1 = (tdotvtang*tdotuv - tdotutang*tdotvv)/tdiv;
  *coef2 = (tdotutang*tdotuv - tdotvtang*tdotuu)/tdiv;

  /* Find third output factor.  */

  for (ki=0; ki<idim; ki++) 
    sdum[ki] = (double)isign*etang[ki] - *coef1*evecu[ki] - *coef2*evecv[ki];
  *coef3 = s6length(sdum,idim,&kstat)/s6length(evecw,idim,&kstat);
  
  if (s6scpr(sdum,evecw,idim) < DZERO) (*coef3) *= -(double)1.0;

  goto out;


  out :
    return;
}
