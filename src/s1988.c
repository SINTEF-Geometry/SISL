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
 * $Id: s1988.c,v 1.2 2001-03-19 15:58:58 afr Exp $
 *
 */


#define S1988

#include "sislP.h"                                                 


#if defined(SISLNEEDPROTOTYPES)
void s1988(SISLCurve *pc,double **emax,double **emin,int *jstat)
#else
void s1988(pc,emax,emin,jstat)
     SISLCurve *pc;
     double **emax;
     double **emin;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find the bounding box of the SISLCurve. NB. The geometric
*              bounding box is returned also in the rational case, that
*              is the box in homogenous coordinates is NOT computed.
*
*
* INPUT      : pc        - SISLCurve to treat.
*
* OUTPUT     : emin      - Array of dimension idim containing
*                          the minimum values of the bounding box,
*                          i.e. down-left corner of the box.
*              emax      - Array of dimension idim containing
*                          the maximum values of the bounding box,
*                          i.e. top-right corner of the box.
*              jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                                                     
*
*
* METHOD     :
*
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF Oslo, July 1993.
*
*********************************************************************
*/                                     
{
  int i,j;                          /* Loop control variables    */
  int kpos = 0;                     /* Position of error.        */
  int bsdim;
  int len;
  int in = pc->in;
  double *coeff;
  double *minim=SISL_NULL;
  double *maxim=SISL_NULL;

  /* initialize variables */

  bsdim = pc->idim;
  coeff = pc->ecoef;
  len = bsdim;

  minim = newarray(bsdim, DOUBLE);
  maxim = newarray(bsdim, DOUBLE);
  if(minim == SISL_NULL || maxim == SISL_NULL) goto err101;

  for(j=0; j<bsdim; j++)
    {
      minim[j] = coeff[j];
      maxim[j] = coeff[j];
    }
  for(i=1, len=bsdim; i<in; i++, len+=bsdim)
    for(j=0; j<bsdim; j++)
      {
	minim[j] = MIN(minim[j], coeff[len+j]);
	maxim[j] = MAX(maxim[j], coeff[len+j]);
      }
  *emin = minim;
  *emax = maxim;      

  /* Success ! */
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
  err101: 
    *jstat = -101;
    s6err("s1988",*jstat,kpos);
    goto out;
  
  out: 
    return;
}
