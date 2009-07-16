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
 * $Id: s6rotax.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6ROTAX

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6rotax(double ep[],double eaxis[],double expnt[],double emat[],int *jstat)
#else
void s6rotax(ep,eaxis,expnt,emat,jstat)
     double ep[];
     double eaxis[];
     double expnt[];
     double emat[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make a matrix that transforms the z-axis ont to
*              the axis specified by ep and eaxis, and transforms
*              the point (1,0,0) onto the point expnt. The matrix is
*              prepared for post multiplication of vectors.
*
* INPUT      : ep      - SISLPoint on axis
*              eaxis   - Direction of axis 
*              expnt   - The point (1,0,0) is trnasformed on to.
*
* OUTPUT     : jstat   - Status message
*                                        >0      : Warning
*                                        =0      : ok
*                                        <0      : Error
*              emat    - The 4x4 matrix produced.
*
*
* METHOD     : 
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-23
*
*********************************************************************
*/
{
  double sdiff[3];            /* Vector for storage of  expnt-ep        */
  double saxis[3];            /* Storage of normalized eaxis            */
  double tradius;             /* Distance from expnt to axis            */
  double sxdir[3];            /* Transformation direction of x-axis     */
  double sydir[3];            /* Transformation direction of y-axis     */
  double strans[3];           /* Translation vector for the origin      */
  double zfak;                /* Length of projection                   */
  int    kdim=3;              /* The dimension of the space we work in  */
  int    kstat;               /* Local status varaible                  */
  int    ki;                  /* Variable used in for loop              */
  
  
  /* Normalize direction of axis */
  
  (void)s6norm(eaxis,kdim,saxis,&kstat);
  
  
  /*  Make difference of expnt  and ep  */
  
  for (ki=0;ki<3;ki++)
    sdiff[ki] = expnt[ki] - ep[ki];
  
  /* Find projection of expnt-ep onto saxis */
  
  zfak = s6scpr(sdiff,saxis,kdim);
  
  
  /* Make transformation of the vector (1,0,0) */
  
  for (ki=0;ki<3;ki++)
    sxdir[ki] = sdiff[ki] - zfak*saxis[ki];

  /* Normalize sxdir */
  
  tradius = s6norm(sxdir,kdim,sxdir,&kstat);
  
  /* Make the vector (0,1,0) is to be projected onto */
  
  s6crss(saxis,sxdir,sydir);
  (void)s6norm(sydir,kdim,sydir,&kstat);
  
  /* Make translation of origo */
  
  for (ki=0;ki<3;ki++)
    strans[ki] = ep[ki] + zfak*saxis[ki];
  
  /* Build transformation matrix for post multiplication of vectors */
  
  emat[0]  = tradius*sxdir[0];
  emat[1]  = tradius*sxdir[1];
  emat[2]  = tradius*sxdir[2];
  emat[3]  = (double)0.0;
  emat[4]  = tradius*sydir[0];
  emat[5]  = tradius*sydir[1];
  emat[6]  = tradius*sydir[2];
  emat[7]  = (double)0.0;
  emat[8]  = tradius*saxis[0];
  emat[9]  = tradius*saxis[1];
  emat[10] = tradius*saxis[2];
  emat[11] = (double)0.0;
  emat[12] = strans[0];
  emat[13] = strans[1];
  emat[14] = strans[2];
  emat[15] = (double)1.0;
  
  /* Matrix made */
  
  *jstat = 0;
  goto out;
  
 out: 
  return;
}
