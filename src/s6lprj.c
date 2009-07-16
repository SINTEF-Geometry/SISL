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
 * $Id: s6lprj.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */
#define S6LPRJ

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
s6lprj(double e1[],double e2[],int idim)
#else
double s6lprj(e1,e2,idim)
     double e1[];
     double e2[];
     int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the length of the projection of a vector on 
*              an arbitrary axis.
*
* INPUT      : e1      - The vector to be projected
*              e2      - The arbitrary axis vector
*              dim     - Number of dimensions in the space the vectors lie
*
* OUTPUT     : s6lprj  - The length of the projection vector
*
* METHOD     : The length of the projection vector is calculated as:
*                       __     __
*                       e1 dot e2
*              ||e3|| = ---------
*                       __     __  1/2
*                      (e2 dot e2)
*-
* CALLS      : s6scpr, s6length
*
* WRITTEN BY : Per Evensen, SI, Oslo, Norway. 1991-aug-16
*                                  
*********************************************************************
*/
{
  int kstat;
  double scpr1,scpr2,lproj; 

  scpr1 = s6scpr(e1,e2,idim);
  scpr2 = s6length(e2,idim,&kstat);
  if (kstat == 0) scpr2=0.000001;
  
  lproj = scpr1/scpr2;
  return(lproj);
}

