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
 * $Id: s1439.c,v 1.2 1994-12-05 15:46:49 pfu Exp $
 *
 */


#define S1439

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1439(SISLSurf *ps1,double apar,int idirec,SISLCurve **rcurve,int *jstat)
#else
void s1439(ps1,apar,idirec,rcurve,jstat)
     SISLSurf   *ps1;
     double apar;
     int idirec;
     SISLCurve  **rcurve;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Pick a curve along a constant parameter line in a NURBS
*              surface.
*              The constant parameter value used is apar and is in the
*              idirec parameter direction.
*              This routine replaces s1436() and s1437().
*
*
*
* INPUT      : ps1    - Surface.
*              apar   - Parameter value to use when picking out constant
*                       parameter curve.
*              idirec - Parameter direction in which to pick (must be 1 or 2)
*
*
* OUTPUT     : rcurve - Constant parameter curve.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      : s1436, s1437 - These two routines do the job, which
*                             one is called depends on what parameter
*                             direction to pick from.
*
* WRITTEN BY : Christophe Rene Birkelan, SINTEF Oslo, July 1993.
*
*********************************************************************
*/
{
  int kpos = 0;      /* Position of error.                            */

  if(idirec == 1)
    {
      s1437(ps1, apar, rcurve, jstat);
      if(*jstat < 0) goto error;
    }
  else if(idirec == 2)
    {
      s1436(ps1, apar, rcurve, jstat);
      if(*jstat < 0) goto error;
    }
  else
    goto err115;

  /* Success !  Curve picked */

  goto out;


  /* Error in input parameter idirec.  */

  err115:
    *jstat = -115;
    s6err("s1439",*jstat,kpos);
    goto out;

  /* Error in lower level routine.  */

  error:
    s6err("s1439",*jstat,kpos);
    goto out;

  out:
    return;
}
