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
 * $Id: s6drawseq.c,v 1.3 1994-12-19 16:58:19 pfu Exp $
 *
 */


#define S6DRAWSEQ

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
extern void s6move(DOUBLE[]);
extern void s6line(DOUBLE[]);
#else
extern void s6move();
extern void s6line();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  s6drawseq(double epoint[],int ipoint)
#else
void s6drawseq(epoint,ipoint)
     double epoint[];
     int ipoint;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Draw a broken line as a sequence of straight lines
*              described by the array epoint.
*
*
*
* INPUT      : epoint - Array describing the corners of the broken line
*                       to draw. Each corner is given as a 3-dimensional
*                       point. The corners is placed continuous in the
*                       array.
*              ipoint - Number of points in epoint.
*
*
* OUTPUT     :
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* REMARK     : This routine is machine-dependant and the internal of
*              this routine has to be reprogrammed if the functions
*              s6move() and s6line() are not defined for your graphics
*              sub-system.
*
*
*-
* CALLS      : s6move - Place pen at given position (empty dummy routine).
*              s6line - Draw a line from the current position to the given one
*                       (empty dummy routine).
*
* WRITTEN BY :
*
*********************************************************************
*/
{
  int ki;          /* Counter.                                    */
  double *spoint;  /* Pointer to corner point in the broken line. */

  /* Position pen at start of the broken line.  */

  s6move(epoint);

  /* Draw sequence of line-segments.  */

  for (ki=1,spoint=epoint+3; ki<ipoint; ki++,spoint+=3)
    s6line(spoint);

  return;
}
