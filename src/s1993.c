/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: s1993.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1993

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1993(SISLCurve *c1,int *jstat)
#else
void s1993(c1,jstat)
     SISLCurve *c1;
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
*********************************************************************
*/
{
  register int ki;

  int kk,kn;
  int kbez;
  double tmax;
  double tmin;
  double tdiff;
  double *scoef=NULL;
  double noice = (double)100.0 * REL_COMP_RES;   /* Noice killer */ 
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
  
  for (ki = 1,scoef= c1->ecoef;ki < kn; ki++,scoef++ )
      {
	tdiff = *(scoef + 1) - *scoef;
	tmin = min(tmin,tdiff);
	tmax = max(tmax,tdiff);
      }
  
  if (fabs(tmin) < noice) tmin = DNULL; 
  if (fabs(tmax) < noice) tmax = DNULL; 
  
  
  /* Simple case when no genuin zero's of first derivative. */
  if (kbez && (tmin*tmax >=DNULL)) 
    *jstat = 1;
  else if (tmin*tmax > DNULL) 
    *jstat = 1;
  else if (tmin == tmax)
    *jstat = 1;
  else
    /* Not a simple case. */
    *jstat = 0;

}



