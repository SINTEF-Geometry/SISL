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
 * $Id: s6inv4.c,v 1.2 2005-02-28 09:04:50 afr Exp $
 *
 */



#define S6INV4

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s6inv4 (double em[], double einv[], int *jstat)
#else
void
s6inv4 (em, einv, jstat)
     double em[];
     double einv[];
     int *jstat;
#endif
/*
*************************************************************************
*
* Purpose: To invert a 4*4-matrix.
*
* Input:
*        Em      - Array (length 16) containing the matrix to be inverted.
*
* Output:
*        Einv    - The inverted matrix.
*
* Calls: s6err.
*
* Written by: A.M. Ytrehus, SI Oslo Feb.92.
*
*****************************************************************
*/
{
  int ki;
  double det;

  *jstat = 0;


  /* Calculate the determinant of the matrix em. */

  det = em[0] * em[5] * (em[10] * em[15] - em[14] * em[11])
    - em[0] * em[6] * (em[9] * em[15] - em[13] * em[11])
    + em[0] * em[7] * (em[9] * em[14] - em[13] * em[10])
    - em[1] * em[4] * (em[10] * em[15] - em[14] * em[11])
    + em[1] * em[6] * (em[8] * em[15] - em[12] * em[11])
    - em[1] * em[7] * (em[8] * em[14] - em[12] * em[10])
    + em[2] * em[4] * (em[9] * em[15] - em[13] * em[11])
    - em[2] * em[5] * (em[8] * em[15] - em[12] * em[11])
    + em[2] * em[7] * (em[8] * em[13] - em[12] * em[9])
    - em[3] * em[4] * (em[9] * em[14] - em[13] * em[10])
    + em[3] * em[5] * (em[8] * em[14] - em[12] * em[10])
    - em[3] * em[6] * (em[8] * em[13] - em[12] * em[9]);


  /* Calculate the inverse matrix. */

  einv[0] = em[5] * (em[10] * em[15] - em[14] * em[11])
    - em[6] * (em[9] * em[15] - em[13] * em[11])
    + em[7] * (em[9] * em[14] - em[13] * em[10]);

  einv[4] = -em[4] * (em[10] * em[15] - em[14] * em[11])
    + em[6] * (em[8] * em[15] - em[12] * em[11])
    - em[7] * (em[8] * em[14] - em[12] * em[10]);

  einv[8] = em[4] * (em[9] * em[15] - em[13] * em[11])
    - em[5] * (em[8] * em[15] - em[12] * em[11])
    + em[7] * (em[8] * em[13] - em[12] * em[9]);

  einv[12] = -em[4] * (em[9] * em[14] - em[13] * em[10])
    + em[5] * (em[8] * em[14] - em[12] * em[10])
    - em[6] * (em[8] * em[13] - em[12] * em[9]);


  einv[1] = -em[1] * (em[10] * em[15] - em[14] * em[11])
    + em[2] * (em[9] * em[15] - em[13] * em[11])
    - em[3] * (em[9] * em[14] - em[13] * em[10]);

  einv[5] = em[0] * (em[10] * em[15] - em[14] * em[11])
    - em[2] * (em[8] * em[15] - em[12] * em[11])
    + em[3] * (em[8] * em[14] - em[12] * em[10]);

  einv[9] = -em[0] * (em[9] * em[15] - em[13] * em[11])
    + em[1] * (em[8] * em[15] - em[12] * em[11])
    - em[3] * (em[8] * em[13] - em[12] * em[9]);

  einv[13] = em[0] * (em[9] * em[14] - em[13] * em[10])
    - em[1] * (em[8] * em[14] - em[12] * em[10])
    + em[2] * (em[8] * em[13] - em[12] * em[9]);


  einv[2] = em[1] * (em[6] * em[15] - em[14] * em[7])
    - em[2] * (em[5] * em[15] - em[13] * em[7])
    + em[3] * (em[5] * em[14] - em[13] * em[6]);

  einv[6] = -em[0] * (em[6] * em[15] - em[14] * em[7])
    + em[2] * (em[4] * em[15] - em[12] * em[7])
    - em[3] * (em[4] * em[14] - em[12] * em[6]);

  einv[10] = em[0] * (em[5] * em[15] - em[13] * em[7])
    - em[1] * (em[4] * em[15] - em[12] * em[7])
    + em[3] * (em[4] * em[13] - em[12] * em[5]);

  einv[14] = -em[0] * (em[5] * em[14] - em[13] * em[6])
    + em[1] * (em[4] * em[14] - em[12] * em[6])
    - em[2] * (em[4] * em[13] - em[12] * em[5]);


  einv[3] = -em[1] * (em[6] * em[11] - em[10] * em[7])
    + em[2] * (em[5] * em[11] - em[9] * em[7])
    - em[3] * (em[5] * em[10] - em[9] * em[6]);

  einv[7] = em[0] * (em[6] * em[11] - em[10] * em[7])
    - em[2] * (em[4] * em[11] - em[8] * em[7])
    + em[3] * (em[4] * em[10] - em[8] * em[6]);

  einv[11] = -em[0] * (em[5] * em[11] - em[9] * em[7])
    + em[1] * (em[4] * em[11] - em[8] * em[7])
    - em[3] * (em[4] * em[9] - em[8] * em[5]);

  einv[15] = em[0] * (em[5] * em[10] - em[9] * em[6])
    - em[1] * (em[4] * em[10] - em[8] * em[6])
    + em[2] * (em[4] * em[9] - em[8] * em[5]);


  for (ki = 0; ki < 16; ++ki)
    einv[ki] /= det;

  return;
}
