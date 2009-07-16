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
 * $Id: sh1994.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH1994

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void sh1994(SISLSurf *s1,double aepsge,int *jstat)
#else
void sh1994(s1,aepsge,jstat)
     SISLSurf *s1;
     double aepsge;
     int  *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a point-surface intersection in one dimention
*              is a simple case,
*              i.e. the intersection will result in one single point.
*
*
*
* INPUT      : s1     - Surface in the intersection problem.
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
* WRITTEN BY : TDO SI, 89-08.
* REWISED BY : Vibeke Skytt, SI, 91-02.
*              UJK, SI,91-10 Bug in first direction when a row
*                            ends with 2 or more equal numbers.
*                            Also added a test to include
*                            tmax < tmin (both HUGE)
*********************************************************************
*/
{
  register int ki,kj,kh;
  int kk1, kk2, kn1, kn2;
  int kbez;
  
  double tmaxt, tmaxs;
  double tmint, tmins;
  double tdiff;
  double *scoef=SISL_NULL;
  
  /* Init to  simple case. */
  *jstat = 1;
  
  tmaxt = tmaxs = - HUGE;
  tmint = tmins =   HUGE;
  
  /* Get surface attributes. */
  kk1  = s1->ik1;
  kk2  = s1->ik2;
  kn1  = s1->in1;
  kn2  = s1->in2;
  kbez = (kk1 == kn1) && (kk2 == kn2); 
  
  
  /* If the surface is linear in some direction it is simpel case. */
  if ((kk1 == 2 && kn1 == 2) || (kk2 == 2 && kn2 == 2)) goto out;
  
  
  /* Run through vertices in first parameter direction to find
     intervall of first derivative. */
  
  /* UJK, 91-10 */
  /* for (kj=0, scoef=s1->ecoef; kj<kn2; kj++,scoef++) */
  for (kj=0, scoef=s1->ecoef; kj<kn2; kj++,scoef=s1->ecoef+kn1*kj)
     for (tdiff=DZERO, ki=1; ki<kn1; ki+=kh, scoef+=kh)
     {
	for (kh=1; ki+kh<=kn1; kh++)
	{
	   if (tdiff*(*(scoef+kh) - *(scoef+kh-1)) < DZERO)
	      {
		 scoef += (kh-1);
		 ki += (kh-1);
		 kh = 1;
	      }
	      tdiff = *(scoef + kh) - *scoef;
	      if (fabs(tdiff) >= aepsge) break;
	}
	if (ki+kh > kn1) break;
	
	tmint = min(tmint,tdiff);
	tmaxt = max(tmaxt,tdiff);
     }
  
  /* Run through vertices in second parameter direction to find
     intervall of first derivative. */
  
  for (ki=0; ki<kn1; ki++)
     for (tdiff=DZERO, kj=1, scoef=s1->ecoef+ki; kj<kn2; kj+=kh, scoef+=kh*kn1)
     {
	for (kh=1; kj+kh<=kn2; kh++)
	{
	   if (tdiff*(*(scoef+kh*kn1) - *(scoef+(kh-1)*kn1)) < DZERO)
	      {
		 scoef += (kh-1)*kn1;
		 kj += (kh-1);
		 kh = 1;
	      }
	      tdiff = *(scoef + kh*kn1) - *scoef;
	      if (fabs(tdiff) >= aepsge) break;
	}
	if (kj+kh > kn2) break;
	
	tmins = min(tmins,tdiff);
	tmaxs = max(tmaxs,tdiff);
     }

  /* UJK, 91-10, maybe parameters not set */
  if (tmint > tmaxt || tmins > tmaxs)
  {
     *jstat = 1;
     goto out;
  }
  
  /* The first derivatives decide directions of possible intersection curves. */
  if (kbez && (tmint*tmaxt >=DZERO || tmins*tmaxs >=DZERO))
    *jstat = 1;
  else if (tmint*tmaxt > DZERO || tmins*tmaxs > DZERO) 
    *jstat = 1;
  else if (tmint == tmaxt  || tmins == tmaxs) 
    *jstat = 1;
  else
    /* Not a simple case. */
    *jstat = 0;
  
  goto out;
 out: ;
}

