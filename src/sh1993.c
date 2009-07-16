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
 * $Id: sh1993.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */

#define SH1993

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void sh1993(SISLCurve *c1,double aepsge,int *jstat)
#else
void sh1993(c1,aepsge,jstat)
     SISLCurve *c1;
     double aepsge;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a point-curve intersection in one dimention
*              is a simple case,
*              i.e. the intersection will result in one single point.
*
*
*
* INPUT      : c1    - SISLCurve in the intersection problem.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jstat  - status messages  
*                                         = 1      : simpel case.
*                                         = 0      : not simpel case.
*                                         < 0      : error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-06.
* Revised by : ALA and UJK 01.11.90, totale rewritten, same strategy
*              as in s1994 (surface case).
* REWISED BY : Vibeke Skytt, SI, 91-02.
*
*********************************************************************
*/
{
  register int ki,kj;

  int kk,kn;
  int kbez;
  double tmax;
  double tmin;
  double tdiff;
  double *scoef=SISL_NULL;
  /* ----------------------------------------------------------- */
  
  /* Init to  simple case. */
  *jstat = 1;
  
  tmax = - HUGE;
  tmin =   HUGE;
  
  /* Get curve attributes. */
  kk  = c1->ik;
  kn  = c1->in;
  kbez = (kk == kn);
  
  /* Run through vertices to find
     intervall of first derivative. */
  
  for (tdiff=DZERO, ki=1, scoef=c1->ecoef; ki<kn; ki+=kj, scoef+=kj)
  {
     for (kj=1; ki+kj<=kn; kj++)
     {
	if (tdiff*(*(scoef+kj) - *(scoef+kj-1)) < DZERO)
	   {
	      scoef += (kj-1);
	      ki += (kj-1);
	      kj = 1;
	   }
	   tdiff = *(scoef + kj) - *scoef;
	   if (fabs(tdiff) >= aepsge) break;
     }
     if (ki+kj > kn) break;
     
     tmin = min(tmin,tdiff);
     tmax = max(tmax,tdiff);
  }
  
  
  /* Simple case when no genuin zero's of first derivative. */
  if (kbez && (tmin*tmax >=DZERO)) 
    *jstat = 1;
  else if (tmin*tmax > DZERO) 
    *jstat = 1;
  else if (tmin == tmax)
    *jstat = 1;
  else
    /* Not a simple case. */
    *jstat = 0;

}



