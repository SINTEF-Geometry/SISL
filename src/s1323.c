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
 * $Id: s1323.c,v 1.2 2001-03-19 15:58:44 afr Exp $
 *
 */


#define S1323

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1323(double etop[],double eaxis[],double econe[],int idim,
	   int inumb,double carray[],int *jstat)
#else
void s1323(etop,eaxis,econe,idim,inumb,carray,jstat)
     double etop[];
     double eaxis[];
     double econe[];
     int    idim;
     int    inumb;
     double carray[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make a matrix of dimension (idim+1)x(idim+1)
*              describing a cone as an implicit function.
*
*
* INPUT      : etop   - The top point of the cone 
*              edirec - Direction of cylinder axis
*              econe  - A point on the cone surface different from the
*                       top point
*              idim   - The dimension of the space the cylinder lies
*              inumb  - The number of copies that are to be made of the
*                       matrix.
*
*
*
* OUTPUT     : carray - The description of the cone. Outside 
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
*     If the top point of the cone is denoted (X0,Y0,Z0), the direction
*     vector of the cone axis is denoted (WX,WY,WZ) and COS(T) is
*     the cosine of the opining angle of the cone, the matrix describing
*     the cone is:
*
*     I-                                                          -I
*     I         WX*WX           -WX*WY            -WX*WZ           I
*     I   1 - ----------      -----------       -----------    A   I
*     I       (COS(T)**2)     (COS(T)**2)       (COS(T)**2)        I
*     I                                                            I
*     I       -WX*WY               WY*WY         -WY*WZ            I
*     I     -----------     1 - -----------     -----------    B   I
*     I     (COS(T)**2)         (COS(T)**2)     (COS(T)**2)        I
*     I                                                            I
*     I       -WX*WZ            -WY*WZ              WZ*WZ          I
*     I     -----------       -----------    1 - -----------   C   I
*     I     (COS(T)**2)       (COS(T)**2)        (COS(T)**2)       I
*     I                                                            I
*     I          A                 B                 *         D   I
*     I-                                                          -I
*
*     WHERE
*         A = (X0*WX*WX+WX*(Y0*WY+Z0*WZ))/(COS(T)**2)-X0
*         B = (Y0*WY*WY+WY*(Z0*WZ+X0*WX))/(COS(T)**2)-Y0
*         C = (Z0*WZ*WZ+WZ*(X0*WX+Y0*WY))/(COS(T)**2)-Z0
*         D = X0*X0+Y0*Y0+Z0*Z0-(X0*X0*WX*WY+Y0*Y0*WY*WY+Z0*Z0*WZ*WZ
*             +2*X0*Y0*WX*WY+2*Y0*Z0*WY*WZ+2*Z0*X0*WZ*WX)/(COS(T)**2)
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
*                                    
*********************************************************************
*/
{
  int kdimp1;         /* Dimension of matrix kdimp1 = idim + 1         */
  int kdimp2;         /* idim + 2                                      */
  int kstop;          /* Stop condition for for loop                   */
  int ki,kj,kl;       /* Running variables in loop                     */
  int kpos=0;         /* Position of error                             */
  int kstat;          /* Local status variable                         */
  double twx,twy,twz; /* Local version of normalized direction vector  */
  double tx0,ty0,tz0; /* Local version of point on axis                */
  double temp;        /* Temporary storage variable                    */
  double tcost2;      /* The square of the cosine of the opening angle */
  double sdirec[3];   /* Normalized direction of cone axis             */
  double sdcone[3];   /* Normalized vector from top point to cone point*/
  
  
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
  
  /* Normalize direction vector of axis */
  
  s6diff(etop,eaxis,idim,sdirec);
  (void)s6norm(sdirec,idim,sdirec,&kstat);
  
  /* Normalize vector from top point to point on conic surface */        
  
  s6diff(etop,econe,idim,sdcone);
  (void)s6norm(sdcone,idim,sdcone,&kstat);
  
  /* Make cosinus of angle between the two normalized vectors */
  
  temp = s6scpr(sdirec,sdcone,idim);
  tcost2 = temp*temp;
  
  /* Test if cone degenerate */
  if (DEQUAL(tcost2,DZERO)) goto err174;
  
  /* Make diagonal elements */
  
  for (ki=0,kl=0 ; ki<kstop ; kl++,ki+=kdimp2)
    {
      temp = sdirec[kl];
      carray[ki] = (double)1.0 - temp*temp/tcost2;
    }                                                                          
  
  twx = sdirec[0];
  twy = sdirec[1];
  twz = sdirec[2];
  tx0 = etop[0];
  ty0 = etop[1];
  tz0 = etop[2];
  
  
  /* Make element (1,4) and (4,1) */
  
  temp = (tx0*twx*twx + twx*(ty0*twy+tz0*twz))/tcost2 - tx0;
  
  carray[3]  = temp;
  carray[12] = temp;
  
  
  /* Make element (2,4) and (4,2) */
  
  temp = (ty0*twy*twy + twy*(tz0*twz+tx0*twx))/tcost2 - ty0;
  
  carray[7]  = temp;
  carray[13] = temp;
  
  
  /* Make element (3,4) and (4,3) */
  
  temp = (tz0*twz*twz + twz*(tx0*twx+ty0*twy))/tcost2 - tz0;
  
  carray[11] = temp;
  carray[14] = temp;
  
  /* Make element (4,4) */
  
  temp = tx0*tx0 + ty0*ty0 + tz0*tz0
         - ( tx0*tx0*twx*twx + ty0*ty0*twy*twy + tz0*tz0*twz*twz
         + (double)2.0*tx0*ty0*twx*twy + (double)2.0*ty0*tz0*twy*twz
         + (double)2.0*tz0*tx0*twz*twx )/tcost2;
  carray[15] = temp;
  
  
  /* Make element (1,2) and (2,1) */
  
  temp = -twx*twy/tcost2;
  carray[1] = temp;
  carray[4] = temp;
  
  
  /* Make element (1,3) and (3,1) */
  
  temp = -twx*twz/tcost2;
  carray[2] = temp;
  carray[8] = temp;
  
  
  /* Make element (2,3) and (3,2) */
  
  temp = -twy*twz/tcost2;
  carray[6] = temp;
  carray[9]= temp;
  
  
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
  s6err("s1323",*jstat,kpos);
  goto out;
  
  /* Degenerate cond */
 err174: *jstat = -174;
  s6err("s1323",*jstat,kpos);
  goto out;
 out:
  return;
}
