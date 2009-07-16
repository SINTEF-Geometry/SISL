/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"


#define S1951

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1951(double etau[], double ecoef[], int in, int ik, int idim, 
	 int ilend, int irend, int incont, double efac[])
#else
void s1951(etau, ecoef, in, ik, idim, ilend, irend, incont, efac)
   double etau[];
   double ecoef[];
   int    in;
   int    ik;
   int    idim;
   int    ilend;
   int    irend;
   int    incont;
   double efac[];
#endif     
/*
*********************************************************************
* 
* PURPOSE    : Multiply the coefficients by dtau(-1/2) and express 
*              the incont last coefficients as a weighted sum
*              of the incont first coeffecients. The weights are given
*              in efac.
* 
* 
* INPUT      : ecoef  - Coefficients of spline curve.
*              in     - Number of coefficients.p in
*              idim   - Dimension of geometry space.
*              incont - Number of continuity conditions, i.e. number of
*                       coefficients at the end to be expressed by
*                       coefficients at the start.
*              efac   - Factors.
*              
*
* 
* OUTPUT     : ecoef  - Coefficients of spline curve.
*             
* 
* METHOD     : 
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      :   
*
* WRITTEN BY : Vibeke Skytt,  SINTEF Oslo, 01.95.
*
*********************************************************************
*/
{
   int ki, kj, kr;   /* Counters.  */
   int kstop;
   double tw;

   /* Multiply the part of ec pointed to by kstart and kstop by the
      corresponding parts of the square matrix dtau(-1/2).   */
   
   for (kstop=in-MAX(incont,irend), ki=ilend; ki<kstop; ki++)
     {
	tw = sqrt((double)ik/(etau[ki+ik] - etau[ki]));
	for (kj=0; kj<idim; kj++)
	  ecoef[ki*idim+kj] *= tw;
     }  
   
   /* Express the incont last coefficients by the incont first ones given
      the factors stored in efac. See s1947. */
   
   for (ki=0; ki<incont; ki++)
   {
      for (kr=0; kr<idim; kr++)
      {
	 ecoef[(in-ki-1)*idim+kr] = DZERO;
	 for (kj=0; kj<=ki; kj++)
	    ecoef[(in-ki-1)*idim+kr] += ecoef[kj*idim+kr]*efac[ki*incont+kj];
      }
   }
   
}
   
