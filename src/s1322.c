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
 * $Id: s1322.c,v 1.3 2001-03-19 15:58:44 afr Exp $
 *
 */


#define S1322

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1322(double epoint[],double edirec[],double aradiu,int idim,
	   int inumb,double carray[],int *jstat)
#else
void s1322(epoint,edirec,aradiu,idim,inumb,carray,jstat)
     double epoint[];
     double edirec[];
     double aradiu;
     int    idim;
     int    inumb;
     double carray[];
     int    *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To make a matrix of dimension (idim+1)x(idim+1)
*              describing a cylinder as an implicit function.
*
*
* INPUT      : epoint - SISLPoint on the cylinder axis
*              edirec - Direction of cylinder axis
*              aradiu - Radius of hyper sphere
*              idim   - The dimension of the space the cylinder lies
*              inumb  - The number of copies that are to be made of the
*                       matrix.
*
*
*
* OUTPUT     : carray - The description of the cylinder. Outside
*                       this function the space for this array must be
*                       allocated. The need is (idim+1)*(idim+1)*inumb
*                       dimension 4x4 (xinarr)
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*     The matrix is described in the following way: (X0,Y0,Z0) point
*     on cylinder axis, (WX,WY,WZ) direction of cylinder axis and
*     R radius of cylinder:
*
*          I-                                 -I
*          I   1-WX*WX  -WX*WY   -WX*WZ    A   I
*          I                                   I
*          I   -WX*WY   1-WY*WY  -WY*WZ    B   I
*          I                                   I
*          I   -WX*WZ   -WY*WZ   1-WZ*WZ   *   I
*          I                                   I
*          I      A        B        C      D   I
*          I-                                 -I
*
*     where
*
*         A = X0*(WX*WX-1)+WX*(Y0*WY+Z0*WZ)
*         B = Y0*(WY*WY-1)+WY*(Z0*WZ+X0*WX)
*         C = Z0*(WZ*WZ-1)+WZ*(X0*WX+Y0*WY)
*         D = X0*X0+Y0*Y0+Z0*Z0-X0*X0*WX*WX-Y0*Y0*WY*WY-Z0*Z0*WZ*WZ
*                -2*X0*Y0*WX*WY-2*Y0*Z0*WY*WZ-2*Z0*X0*WZ*WX-R*R
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 28-June-1988
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Changed loop
*              building diagonal elements to avoid over-running 'sdirec'.
*
*********************************************************************
*/
{
  int kdimp1;         /* Dimension of matrix kdimp1 = idim + 1         */
  int kdimp2;         /* idim + 2                                      */
  int kstop;          /* Stop condition for for loop                   */
  int ki,kj,kl;       /* Running variables in loop                     */
  int kpos=0;         /* Position of error                             */
  double twx,twy,twz; /* Local version of normalized direction vector  */
  double tx0,ty0,tz0; /* Local version of point on axis                */
  double temp;        /* Temporary storage variable                    */
  double tsum;        /* Varaible used for summation                   */
  double sdirec[3];   /* Normalized direction vector                   */


  /* Test i legal input */
  if (inumb <1 ) inumb = 1;
  if (idim != 3 ) goto err104;

  kdimp1 = idim + 1;
  kdimp2 = idim + 2;
  kstop  = kdimp1*kdimp1;

  for (ki=0;ki<kstop;ki++)
    {
      carray[ki] = DZERO;
    }

  /* Normalize direction vector */

  tsum = DZERO;

  for (ki=0;ki<idim;ki++)
    {
      temp = edirec[ki];
      tsum += (temp*temp);
    }

  tsum = sqrt(tsum);
  if (DEQUAL(tsum,DZERO)) goto err173;

  for (ki=0;ki<idim;ki++)
    {
      sdirec[ki] = edirec[ki]/tsum;
    }

  /* Make diagonal elements */

  for (ki=0,kl=0 ; ki<kstop-1 ; kl++,ki+=kdimp2)   /* (PFU 14/11-1994) */
    {
      temp = sdirec[kl];
      carray[ki] = (double)1.0 - temp*temp;
    }
  carray[ki] = (double) 1.0;  /* (PFU 14/11-1994) */

  /* Make element 1,...,idim of last row and 1,...,idim of last column */

  tsum = DZERO;
  twx = sdirec[0];
  twy = sdirec[1];
  twz = sdirec[2];
  tx0 = epoint[0];
  ty0 = epoint[1];
  tz0 = epoint[2];

  /* Make element (1,4) and (4,1) */

  temp = tx0*(twx*twx-(double)1.0) + twx*(ty0*twy+tz0*twz);

  carray[3]  = temp;
  carray[12] = temp;

  /* Make element (2,4) and (4,2) */

  temp = ty0*(twy*twy-(double)1.0) + twy*(tz0*twz+tx0*twx);

  carray[7]  = temp;
  carray[13] = temp;

  /* Make element (3,4) amd (4,3) */

  temp = tz0*(twz*twz-(double)1.0) + twz*(tx0*twx+ty0*twy);

  carray[11] = temp;
  carray[14] = temp;

  /* Make element (4,4) */

  temp = tx0*tx0*((double)1.0-twx*twx) + ty0*ty0*((double)1.0-twy*twy)
         + tz0*tz0*((double)1.0-twz*twz) - (double)2.0*tx0*ty0*twx*twy
         - (double)2.0*ty0*tz0*twy*twz - (double)2.0*tz0*tx0*twz*twx
	 - aradiu*aradiu;

  carray[15] = temp;

  /* Make element (1,2) and (2,1) */

  temp = -twx*twy;
  carray[1] = temp;
  carray[4] = temp;

  /* Make element (1,3) and (3,1) */

  temp = -twx*twz;
  carray[2] = temp;
  carray[8] = temp;

  /* Make element (2,3) and (3,2) */
  temp = -twy*twz;
  carray[6] = temp;
  carray[9] = temp;

  /* Make extra copies of cylinder */

  kj = kstop;
  for (ki=1;ki<inumb;ki++)
    {
      for (kl=0;kl<kstop;kl++,kj++)
        {
	  carray[kj] = carray[kl];
        }
    }

  *jstat = 0;
  goto out;

  /* Dimension less than 1 */
 err104: *jstat = -104;
  s6err("s1322",*jstat,kpos);
  goto out;

  /* Direction vector of length 0 */
 err173: *jstat = -173;
  s6err("s1322",*jstat,kpos);
  goto out;
 out:
  return;
}
