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
 * $Id: s1500.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1500

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1500(double base[],double norm[],double axisA[],double alpha,
         double ratio,int idim,int inumb,double carray[],int *jstat)
#else
void s1500(base,norm,axisA,alpha,ratio,idim,inumb,carray,jstat)
     double base[];
     double norm[];
     double axisA[];
     double alpha;
     double ratio;
     int    idim;
     int    inumb;
     double carray[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make a matrix of dimension (idim+1)x(idim+1)
*              describing an elliptic cone as an implicit function.
*
*
* INPUT      : base   - The base point of the cone
*                                (the centre of the base ellipse) 
*              norm  - Direction of cone axis
*              axisA  - One of the axes of the ellipse
*              alpha  - The opening angle of axisA
*              ratio  - The ratio of axisA to axisB 
*              idim   - The dimension of the space the cylinder lies in
*              inumb  - The number of copies that are to be made of the
*                       matrix.
*
*
*
* OUTPUT     : carray - The description of the elliptic cone. Outside 
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
*     If the base point of the cone is denoted p=(p1,p2,p3),
*     the axis vector of the cone is denoted n=(n1,n2,n3),
*     the direction vector of one ellipse axis is denoted A=(a1,a2,a3),
*     the second B = (n X A) / |n x A|  is B=(b1,b2,b3),
*     alpha is the opening angle of the cone at the axis A,
*     then the matrix describing
*     the cone is A:
*
*
*     WHERE
*        A11 = cos2*a1*a1+cos2*r*b1*b1-sin2*n1*n1
*  A21 = A12 = cos2*a1*a2+cos2*r*b1*b2-sin2*n1*n2
*  A31 = A13 = cos2*a1*a3+cos2*r*b1*b3-sin2*n1*n3
*  A41 = A14 = -cos2*(A.p)*a1-cos2*r*(B.p)*b1+a*cos*sin*n1+sin2*(n.p)*n1
*        A22 = cos2*a2*a2+cos2*r*b2*b2-sin2*n2*n2
*  A32 = A23 = cos2*a2*a3+cos2*r*b2*b3-sin2*n2*n3
*  A42 = A24 = -cos2*(A.p)*a2-cos2*r*(B.p)*b2+a*cos*sin*n2+sin2*(n.p)*n2
*        A33 = cos2*a3*a3+cos2*r*b3*b3-sin2*n3*n3
*  A43 = A34 = -cos2*(A.p)*a3-cos2*r*(B.p)*b3+a*cos*sin*n3+sin2*(n.p)*n3
*        A44 = cos2*(A.p)**2+cos2*r*(B.p)**2-(cos*a+sin*(n.p))**2
*
*     aND
*
*        cos=cos(alpha), sin=sin(alpha), cos2=cos*cos, sin2=sin*sin
*        a=|A|, r=(|A|/|B|)**2,
*        A.p=a1*p1+a2*p2+a3*p3 (dot product), similarly for B.p, n.p,
*   and  B=n*A/|nxA| (cross product).
*
*        
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Mike Floater, SI, Oslo, Norway, 16-October-1988
*                                    
*********************************************************************
*/
{
  int kdimp1;         /* Dimension of matrix kdimp1 = idim + 1         */
  int kstop;          /* Stop condition for for loop                   */
  int ki,kj,kl;       /* Running variables in loop                     */
  int kpos=0;         /* Position of error                             */
  int kstat;          /* Local status variable                         */
  double a1,a2,a3;    /* Coordinates of A vector                       */
  double b1,b2,b3;    /* Coordinates of B vector                       */
  double n1,n2,n3;    /* Coordinates of ndirec vector                  */
  double temp;        /* Temporary storage variable                    */
  double cosT;        /* cos(alpha)                                    */
  double cosT2;       /* The square of cos(alpha)                      */
  double sinT;        /* sin(alpha)                                    */
  double sinT2;       /* The square of sin(alpha)                      */
  double cosTsinT;    /* cos(alpha)*sin(alpha)                         */
  double ndirec[3];   /* Normalized direction of cone axis             */
  double A[3];        /* Normalized elliptic axisA                     */
  double B[3];        /* Normalized elliptic axisB                     */
  double asize;       /* Length of axisA                               */
  double r;           /* ratio*ratio                                   */
  double Adotp;       /* scalar product of A and base                  */
  double Bdotp;       /* scalar product of B and base                  */
  double ndotp;       /* scalar product of n and base                  */
  
  
  /* Test i legal input */
  if (inumb <1 ) inumb = 1;
  if (idim != 3 ) goto err104;
  
  kdimp1 = idim + 1;
  kstop  = kdimp1*kdimp1;
  
  for (ki=0;ki<kstop;ki++)
    {
      carray[ki] = DZERO;
    }
  
  /* Normalise direction vector of cone axis */
  
  (void)s6norm(norm,idim,ndirec,&kstat);
  
  /* Test if norm degenerate */
  if(kstat == 0) goto err174;

  
  /* Normalise elliptic axisA and find its length */        
  
  asize=s6norm(axisA,idim,A,&kstat);
  
  /* Test if norm degenerate */
  if(kstat == 0) goto err174;

  /* Calculate the other elliptic axis from ndirec and A       */        
  /* Note that norm and axisA are assumed to be orthogonal     */        
  
  (void)s6crss(ndirec,A,B);
  
  /* Calculate cosine and sine of angle alpha */
  
  cosT=cos(alpha);
  sinT=sin(alpha);
  cosT2=cosT*cosT;
  sinT2=sinT*sinT;
  cosTsinT=cosT*sinT;
  
  /* Find the square of the ratio of lengths of axisA and axisB  */

  r=ratio*ratio;

  
  /* set up coordinates of the normalised axes  */

  a1 = A[0];
  a2 = A[1];
  a3 = A[2];
  b1 = B[0];
  b2 = B[1];
  b3 = B[2];
  n1 = ndirec[0];
  n2 = ndirec[1];
  n3 = ndirec[2];
  
  /* set up dot (scalar) products of axes and base point */

  Adotp = s6scpr(A,base,idim);
  Bdotp = s6scpr(B,base,idim);
  ndotp = s6scpr(ndirec,base,idim);

  
  /* Make element (1,1)  */
  
  carray[0] = cosT2*(a1*a1+r*b1*b1)-sinT2*n1*n1;
  
  /* Make element (1,1)  */
  
  carray[5] = cosT2*(a2*a2+r*b2*b2)-sinT2*n2*n2;
  
  /* Make element (1,1)  */
  
  carray[10] = cosT2*(a3*a3+r*b3*b3)-sinT2*n3*n3;
  
  /* Make element (1,2) and (2,1) */
  
  temp = cosT2*(a1*a2+r*b1*b2)-sinT2*n1*n2;
  carray[1] = temp;
  carray[4] = temp;
  
  
  /* Make element (1,3) and (3,1) */
  
  temp = cosT2*(a1*a3+r*b1*b3)-sinT2*n1*n3;
  carray[2] = temp;
  carray[8] = temp;
  
  
  /* Make element (2,3) and (3,2) */
  
  temp = cosT2*(a2*a3+r*b2*b3)-sinT2*n2*n3;
  carray[6] = temp;
  carray[9] = temp;
  
  
  
  /* Make element (1,4) and (4,1) */
  
  temp = -cosT2*(Adotp*a1+r*Bdotp*b1)+(asize*cosTsinT+sinT2*ndotp)*n1;
  carray[3]  = temp;
  carray[12] = temp;
  
  
  /* Make element (2,4) and (4,2) */
  
  temp = -cosT2*(Adotp*a2+r*Bdotp*b2)+(asize*cosTsinT+sinT2*ndotp)*n2;
  carray[7]  = temp;
  carray[13] = temp;
  
  
  /* Make element (3,4) and (4,3) */
  
  temp = -cosT2*(Adotp*a3+r*Bdotp*b3)+(asize*cosTsinT+sinT2*ndotp)*n3;
  carray[11] = temp;
  carray[14] = temp;
  
  /* Make element (4,4) */
  
  temp = cosT*asize+sinT*ndotp;
  carray[15] = cosT2*(Adotp*Adotp+r*Bdotp*Bdotp)-temp*temp;
  
  
  /* Make extra copies of cone */
  
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
  s6err("s1500",*jstat,kpos);
  goto out;
  
  /* norm axis or elliptic axis is zero */
 err174: *jstat = -174;
  s6err("s1500",*jstat,kpos);
  goto out;
  
 out:
  return;
}
