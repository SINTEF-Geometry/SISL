/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s1321.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1321

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1321(double ecentr[],double aradiu,int idim,int inumb,
	   double carray[],int *jstat)
#else
void s1321(ecentr,aradiu,idim,inumb,carray,jstat)
     double ecentr[];
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
*              describing a hyper sphere as an implicit function.
*
*
* INPUT      : ecentr - Center of the hyper sphere
*              aradiu - Radius of hyper sphere
*              idim   - The dimension of the space the hyper sphere lies
*              inumb  - The number of copies that are to be made of the
*                       matrix.
*
*
*
* OUTPUT     : carray - The description of the super sphere. Outside
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
*     The matrix is described in the following way (x,y,z)-center and
*     r radius:
*
*            I-                                  -I
*            I   1    0    0        -x            I      
*            I   0    1    0        -y            I
*            I   0    0    1        -z            I
*            I  -x   -y   -z   x*x+y*y+z*z-r*r    I
*            I-                                  -I
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
  double temp;        /* Temporary storage variable                    */
  double tsum;        /* Varaible used for summation                   */
  
  
  
  /* Test i legal input */
  if (inumb <1 ) inumb = 1;
  if (idim < 1 ) goto err102;
  
  kdimp1 = idim + 1;
  kdimp2 = idim + 2;
  kstop  = kdimp1*kdimp1;
  
  for (ki=0;ki<kstop;ki++)
    {
      carray[ki] = (double)0.0;
    }
  
  /* Make diagonal elements */
  
  for (ki=0;ki<kstop;ki+=kdimp2)
    {
      carray[ki] = (double)1.0;
    }
  
  /* Make element 1,...,idim of last column and element 1,...,idim of last
   *  row */
  
  tsum = (double)0.0;
  for (kl=0,ki=idim,kj=idim*kdimp1;kl<idim;kl++,kj++,ki+=kdimp1)
    {
      temp = -ecentr[kl];
      carray[ki] = temp;
      carray[kj] = temp;
      tsum +=(temp*temp);                                                
    }
  
  /* Make lower right corner element */
  
  carray[kstop-1] = tsum - aradiu*aradiu;
  
  /* Make extra copies of hyper sphere */
  
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
 err102: *jstat = -102;
  s6err("s1321",*jstat,kpos);
  goto out;
 out:
  return;
}
