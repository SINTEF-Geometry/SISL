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
 * $Id: s6crvcheck.c,v 1.3 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6CRVCHECK

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6crvcheck(SISLCurve *pc,int *jstat)
#else
void s6crvcheck(pc,jstat)
     SISLCurve *pc;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To check a curve descripiton and remove unneccessary
*              knots and vertices. Such that no continuous curve will
*              have knots with more than the order minus one in 
*              multiplicity.
*
*
*
* INPUT/OUTPUT:pc     - The curve identifcation
*
* OUTPUT     : kstat  - Status variable
*                        < 0 - Error
*                        = 0 - SISLCurve object not changed
*                        = 1 - SISLCurve object changed
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6dist  - Distance between two points.
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway,  August 1989
* REWISED BY : Vibeke Skytt, SI, August 1990
* REVISED BY : Johannes Kaasa, SI, April 1992 (Intoduced NURBS)
* REVISED BY : Christophe Rene Birkeland, SINTEF, May 1993 (Status variable)
*
*********************************************************************
*/                                                               
{
  int kstat = 0;              /* Status variable.                 */
  int ki,kj;                  /* Counter.                         */
  int kdim;                   /* Dimension of space               */
  int rdim;                   /* Rational dimension.              */
  int kn;                     /* Number of knots                  */
  int kk;                     /* Number of vertices               */
  int kmark;                  /* Indicates if k-tupple knots      */
  int knnew;                  /* New number of vertices           */
  int kind;                   /* Type of curve, 2 and 4 rational. */
  double *snt=SISL_NULL;           /* Compressed knot vector           */
  double *sncoef=SISL_NULL;        /* Compressed vertex vector         */
  double *srcoef=SISL_NULL;        /* Compressed vertex vector         */
  double *st;                 /* Knots                            */
  double *scoef;              /* Vertices                         */
  double *rcoef;              /* Rational vertices.               */
  
  *jstat = 0;

  if (pc == SISL_NULL) goto out;
  
  kk    = pc -> ik;
  kn    = pc -> in;
  kdim  = pc -> idim;
  rdim  = kdim + 1;
  kind  = pc -> ikind;
  st    = pc -> et;
  scoef = pc -> ecoef;
  rcoef = pc -> rcoef;
  
  /* Run through all knots to detect if st[ki]=st[ki+kk-1] e.g. that we
     have at least kk-tupple internal knots */
  
  kmark = 0;
  for (ki=1 ; ki < kn-1 ; ki++)
    if (st[ki] == st[ki+kk-1] && 
	DEQUAL(s6dist(scoef+(ki-1)*kdim,scoef+ki*kdim,kdim),DZERO))
      {
        kmark = 1;
        break;
      }
  
  if (kmark == 0) goto out;
  
  /* We have at least kk-tupple knots, remove not necessary knots and vertices */
  
  if((snt = newarray(kn+kk,DOUBLE)) == SISL_NULL) goto err101;  
  if((sncoef = newarray(kn*kdim,DOUBLE)) == SISL_NULL) goto err101;

  if (kind == 2 || kind == 4)
    {
      srcoef = newarray(kn*rdim,DOUBLE);
      if (srcoef == SISL_NULL) goto err101;
      for (ki=0,kj=0 ; ki < kn ; ki ++)
        if (ki == 0 || ki == kn-1 || st[ki] < st[ki+kk-1] || 
	  DNEQUAL(s6dist(rcoef+(ki-1)*rdim,rcoef+ki*rdim,rdim),DZERO))
          {
            snt[kj] = st[ki];
            memcopy(sncoef+kdim*kj,scoef+kdim*ki,kdim,DOUBLE);
            memcopy(srcoef+rdim*kj,rcoef+rdim*ki,rdim,DOUBLE);
            kj++;
          }
    }
  else
    {
      for (ki=0,kj=0 ; ki < kn ; ki ++)
        if (ki == 0 || ki == kn-1 || st[ki] < st[ki+kk-1] || 
	  DNEQUAL(s6dist(scoef+(ki-1)*kdim,scoef+ki*kdim,kdim),DZERO))
          {
            snt[kj] = st[ki];
            memcopy(sncoef+kdim*kj,scoef+kdim*ki,kdim,DOUBLE);
            kj++;
          }
    }
  
  for (ki=kn ; ki<kn+kk ; ki++,kj++)
    snt[kj] = st[ki];
  
  knnew = kj - kk;
  
  /* An additional end knot might have been left */
  
  if (snt[knnew-1] == snt[knnew+kk-1]) knnew--;
  
  /* Put compressed description back to curve object */      
  
  if (pc->icopy > 0)
    {
      pc -> in = knnew;
      memcopy(pc->et,snt,knnew+kk,DOUBLE);
      memcopy(pc->ecoef,sncoef,knnew*kdim,DOUBLE);
      if (kind == 2 || kind == 4)
        memcopy(pc->rcoef,srcoef,knnew*rdim,DOUBLE);
      kstat = 1;
    }
  
  /* Task done. */
  
  *jstat = kstat;
  goto out;
  
  /* Error in space allocation. */
  
  err101: 
    *jstat = -101;
    goto out;
  
  out:
    if (snt != SISL_NULL) freearray(snt);
    if (sncoef != SISL_NULL) freearray(sncoef);
}
