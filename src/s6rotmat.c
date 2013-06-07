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
 * $Id: s6rotmat.c,v 1.2 2005-02-28 09:04:50 afr Exp $
 *
 */


#define S6ROTMAT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6rotmat (double eorigo[], double exaxis[], double enorm[],
	  double ematrix[], int *jstat)
#else
void 
s6rotmat (eorigo, exaxis, enorm, ematrix, jstat)
     double eorigo[];
     double exaxis[];
     double enorm[];
     double ematrix[];
     int *jstat;
#endif
/*
*************************************************************************
*
* Purpose: To compute the transformation matrix (dim = 4X4), given
*	   points (of dim = 3), in the local coordinate system.
*
* Input:
*        Eorigo  - Origo of the local coordinate system.
*        Exaxis  - Point on (local) positive x-axis.
*        Enorm   - Normal to the xy-plane. (z-axis).
*
* Output:
*        Ematrix - Transformation metrix to the local coordinate system.
*        Jstat   - Output status.
*
* Calls: s6scpr, s6crss, s6err.
*
* Written by: A.M. Ytrehus, SI Oslo Nov.91.
* After FORTRAN (P6INN).
*****************************************************************
*/
{
  double syaxis[3], sxaxis[3];
  double teps = (double) 0.000001;
  double tdum, tdum1, tdum2, tdum3;
  int kp;
  int idim = 3;

  int kpos = 0;

  *jstat = 0;


  /* First column will be the relative X-axis, second column
     the relative Y-axis, third column the relative Z-axis, and
     fourth column the local origo. The fourth line will be 0,0,0,1. */


  /* Find the local X-axis. */

  for (kp = 0; kp < idim; kp++)
    sxaxis[kp] = exaxis[kp] - eorigo[kp];


  /* Calculate the local Y-axis = (Z) X (X). */

  s6crss (enorm, sxaxis, syaxis);


  /* Normalise the three axis-vectors, and put into the matrix. */

  tdum = s6scpr (sxaxis, sxaxis, idim);
  tdum1 = sqrt (tdum);
  if (tdum1 < teps)
    goto err166;

  tdum = s6scpr (syaxis, syaxis, idim);
  tdum2 = sqrt (tdum);
  if (tdum2 < teps)
    goto err166;

  tdum = s6scpr (enorm, enorm, idim);
  tdum3 = sqrt (tdum);
  if (tdum3 < teps)
    goto err166;


  for (kp = 0; kp < idim; kp++)
    {
      ematrix[4 * kp] = sxaxis[kp] / tdum1;
      ematrix[4 * kp + 1] = syaxis[kp] / tdum2;
      ematrix[4 * kp + 2] = enorm[kp] / tdum1;
      ematrix[4 * kp + 3] = eorigo[kp];
      ematrix[12 + kp] = (double) 0.0;
    }

  ematrix[15] = (double) 1.0;

  goto out;

  /* Impossible to create matrix. */

err166:
  *jstat = -166;
  s6err ("s6rotmat", *jstat, kpos);
  goto out;

out:
  return;
}
