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

#define SH6SETSEG

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void sh6setseg(SISLSurf *ps, int idir, double *segmentation, int nseg,
	       int type, int *jstat)
#else
  void sh6setseg(ps, idir, segmentation, nseg, type, jstat)
     SISLSurf *ps;
     int idir;
     double *segmentation;
     int nseg;
     int type;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : 
*
*
*
* INPUT      : ps     - Surface to aquire segmentation information
*              idir   - Parameter direction corresponding to segmentation
*              segmentation - Segmentation parameter values
*              nseg   - Number of segmentation values
*              type   - Segmentation type
*                                                                     
*
* OUTPUT     : jstat  - status messages  
*                            = 1      : Tangential curve found
*			     = 0      : The curve is not tangential or
*			                do not end in corners
*                            < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 2018-02
*
*********************************************************************
*/                                     
{
  int kstat = 0;
  int ki;
  SISLSegmentation *pseg = NULL;    /* Current segmentation direction */
  int *seg_type = NULL;

  seg_type = newarray(nseg, INT);
  if (seg_type == NULL)
    goto error;
  
  for (ki=0; ki<nseg; ++ki)
    seg_type[ki] = type;
  pseg = (idir == 0) ? ps->seg1 : ps->seg2;
  if (idir == 0 && ps->seg1 == NULL)
    {
      /* Create segmentation array */
      ps->seg1 = newSegmentation(segmentation, seg_type, nseg);
    }
  else if (idir == 0)
    {
      /* Merge arrays */
    }
  if (idir == 1 && ps->seg2 == NULL)
    {
      /* Create segmentation array */
      ps->seg2 = newSegmentation(segmentation, seg_type, nseg);
    }
  else /* if (idir == 1) */
    {
      /* Merge arrays */
    }

  *jstat = 0;
  goto out;
  
 error:
  *jstat = kstat;
  goto out;

 out:
  if (seg_type != NULL) freearray(seg_type);

  return;
}
