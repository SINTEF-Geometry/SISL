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
 * $Id: sh1927.c,v 1.3 2001-03-19 15:59:07 afr Exp $
 *
 */

#define SH1927

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      sh1927(double etau[],int ik,int in,int idim,SISLCurve *pcurve,
	     int ilend,int irend,double ec[],int *jstat)
#else
void sh1927(etau,ik,in,idim,pcurve,ilend,irend,ec,jstat)
   double etau[];
   int ik;
   int in;
   int idim;
   SISLCurve *pcurve;
   int ilend;
   int irend;
   double ec[];
   int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : To compute the ilend first and irend last coefficients
*              of a spline on the knot vector etau with 0 - ilend-1
*              derivatives at etau[ik-1] equal to the corresponding derivates
*              of the spline given by et and ed at the point et[ik-1], and
*              similarily for the 0 - irend-1 derivative at etau[in] and
*              et[im] respectively. et denotes the knotvector of pcurve and im
*              is the number of knots.
* 
* INPUT      : etau   - Real array of length (in+ik) containing the  
*                       knot vector of the approximating space.
*	       ik     - The order of the spline space.
*              in     - The dimension of the spline space corresponding
*                       to etau.
*              idim   - The dimension of the geometry space.
*              pcurve - The B-spline to be approximated.
*              ilend  - Parameter giving the number of derivatives to be
*                       kept fixed at etau[ik-1], ie.e at the left end of the
*                       curve. This means that the 0 - ilend-1 derivatives at
*                       etau[ik-1], of the spline given by etau and ec will
*                       be the same as the corresponding derivatives at 
*                       et[ik-1], of the spline given by et and ed. 
*                       NB! If ilend is negative it will be set to zero, if
*                       it is greater than ik, it will be set equal to ik.
*              irend  - Parameter giving the number of derivatives to be
*                       kept fixed at etau[in], ie.e at the right end of the
*                       curve. This means that the 0 - irend-1 derivatives at
*                       etau[in], of the spline given by etau and ec will
*                       be the same as the corresponding derivatives at 
*                       et[im], of the spline given by et and ed. 
*                       NB! If irend is negative it will be set to zero, if
*                       it is greater than ik, it will be set equal to ik.
*
* 
* OUTPUT     : ec     - Real array of lengh (in*idim) containing B-spline
*                       coefficients of the new spline with knot vector etau
*                       with derivatives similar to those of the spline given
*                       by et and ed at the beginning and end of the curve.
*                       Only the first ilend and the last irend sets of 
*                       coefficients (each set consists of idim coefficients)
*                       receive a value in this routine.
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
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
* WRITTEN BY : Vibeke Skytt, SI, 05.92, on the basis of a routine
*              written by Tom Lyche and Knut Moerken, 12.85.
*
*********************************************************************
*/
{ 
   int kstat;
   int kn,km;
   int ki,kj,kr,kp;
   int krh,kjh,khc,kjhc;
   int kstart,kstop;
   int kleft = ik - 1;
   double tm;
   double th;
   double *sder=SISL_NULL;
   double *sbder=SISL_NULL;
   double *ssum = SISL_NULL;
   double *shc = SISL_NULL;

   /* Test input.  */

   if (pcurve->ik != ik) goto err106;
   if (pcurve->idim != idim) goto err103;
  
   /* Allocate scratch for local arrays.  */
   
   if ((sder = newarray(idim*MAX(ilend,irend),DOUBLE)) == SISL_NULL)
      goto err101;
   if ((sbder=newarray(ik*ik,DOUBLE)) == SISL_NULL) goto err101;
   if ((ssum = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;
   if ((shc = newarray(MAX(1,irend)*idim,DOUBLE)) == SISL_NULL) goto err101;
   
   /* Check the input and adjust if necessary.  */
   
   if (ilend < 0) ilend = 0;
   if (ilend > ik) ilend = ik;
   if (irend < 0) irend = 0;
   if (irend > ik) irend = ik;
   
   /* Set output array to zero.  */
   
   memzero(ec,in*idim,DOUBLE);
   
   /* Treat left end if any coefficients are kept fixed.  */
   
   if (ilend > 0)
   {
      /* Compute the necessary derivatives. */
      
      s1221(pcurve,ilend-1,pcurve->et[ik-1],&kleft,sder,&kstat);
      if (kstat < 0) goto error;
      
      s1220(etau,ik,in,&kleft,etau[ik-1],ilend-1,sbder,&kstat);
      if (kstat < 0) goto error;

      /* Compute the multiplicity of etau[ik-1] (from the left).  */
		     
      for (km=1; ; km++)
	if (km == ik || etau[ik-km-1] < etau[ik-km]) break;
      
      /* The first diagonal block is known to be of dimension kn=ik-km. */
      
      kn = ik - km;
	
      /* Perform LU-factorization of the transpose of this first block
	 which is a dense, square matrix.   */
      
      if (km < ik)
      {
	 for (kr=0; kr<kn-1; kr++)
	   for (ki=kr+1; ki<kn; ki++)
	     {
		tm = sbder[kr*ilend+ki]/sbder[kr*ilend+kr];
		for (kj=kr+1; kj<kn; kj++)
		  sbder[kj*ilend+ki] -= tm*sbder[kj*ilend+kr];
		for (kp=0; kp<idim; kp++)
		  sder[ki*idim+kp] -= tm*sder[kr*idim+kp];
	     }
	 
	 /* Find the value of ec[0],...,ec[kn-1] by back substitution.  */
	 
	 tm = (double)1.0/sbder[(kn-1)*ilend+kn-1];
	 for (kp=0; kp<idim; kp++)
	   ec[(kn-1)*idim+kp] = sder[(kn-1)*idim+kp]*tm;
	 for (kr=kn-2; kr>=0; kr--)
	   {
	      memzero(ssum,idim,DOUBLE);
	      for (kj=kr+1; kj<kn; kj++)
		{
		   tm = sbder[kj*ilend+kr];
		   for (kp=0; kp<idim; kp++)
		     ssum[kp] += ec[kj*idim+kp]*tm;
		}
	      
	      tm = (double)1.0/sbder[kr*ilend+kr];
	      for (kp=0; kp<idim; kp++)
		ec[kr*idim+kp] = (sder[kr*idim+kp] - ssum[kp])*tm;
	   }
      }
      
      /* Determine ec[kn],...ec[ilend-1]. This is simple since the
	 second block is upper triangular.  */
      
      for (kr=kn; kr<ilend; kr++)
	{
	   memzero(ssum,idim,DOUBLE);
	   for (kj=0; kj<kr; kj++)
	     {
		tm = sbder[kj*ilend+kr];
		for (kp=0; kp<idim; kp++)
		  ssum[kp] += ec[kj*idim+kp]*tm;
	     }
	   
	   tm = (double)1.0/sbder[kr*ilend+kr];
	   for (kp=0; kp<idim; kp++)
	     ec[kr*idim+kp] = (sder[kr*idim+kp] - ssum[kp])*tm;
	}
      
      /* Finished with the left end.  */
   }
      
   /* Treat right end if any coefficients are kept fixed. */
   
   if (irend > 0)
   {
      /* Compute the derivatives of the original spline at et[im],
	 and the derivatives of the last ik B-splines associated with
	 etau at the point etau[in].  */
      
      kleft = pcurve->in;
      s1221(pcurve,irend-1,pcurve->et[pcurve->in],&kleft,sder,&kstat);
      if (kstat < 0) goto error;
		     
      kleft = in;		     
      s1220(etau,ik,in,&kleft,etau[in],irend-1,sbder,&kstat);
      if (kstat < 0) goto error;
				    
      /* Compute the multiplicity of etau[ik-1] (from the right).  */
				    
      for (km=1; ; km++)
	if (km == ik || etau[in+km-1] < etau[in+km]) break;
      
      /* The coefficient matrix is now block upper triangular (really
	 lower triangular but we are working with the transpose of the 
	 matrix) along the bi-diagonal, with a square block of dimension 
	 (ik-im)*(ik-im) in the upper right corner, and a upper triangular
	 matrix of dimensions kn*km in the lower left corner (in case 
	 irend=ik).
	 The solution of this system is stored in an extra array shc
	 so that the average of the two solutions (the one at the beginning
	 and the one at the end) can be easily calculated where the 
	 overlap.       */

      khc = idim*(irend-1);

      /* if pcurve->in = ik, the whole system is triangular.  */
      
      if (km < ik)
      {
	 /* Perform LU-factorizing of the square matrix in the upper right
	    corner and perform the corresponding forward substitution
	    on sder.    */
	 
	 for (krh=0, kstop=ik-km, kr=km+1; kr<ik-1; kr++)
	   for (krh++, ki=krh+1; ki<kstop; ki++)
	     {
		tm = sbder[kr*irend+ki]/sbder[kr*irend+krh];
		for (kj=kr+1; kj<ik; kj++)
		  sbder[kj*irend+ki] -= tm*sbder[kj*irend+krh];
		for (kp=0; kp<idim; kp++)
		  sder[ki*idim+kp] -= tm*sder[krh*idim+kp];
	     }
	 
	 /* Compute the last ik-km unknowns by back substitution.  */
	 
	 krh = ik - km;
	 tm = (double)1.0/sbder[ik*irend+krh];
	 for (kp=0; kp<idim; khc++, kp++)
	   shc[khc] = sder[krh*idim+kp]*tm;
	 for (kstop=in-ik+km+1, kr=in-1; kr>=kstop; kr--)
	   {
	      khc -= 2*idim;
	      krh--;
	      memzero(ssum,idim,DOUBLE);
	      for (kjhc=khc+idim, kjh=km+krh, kj=kr+1; kj<kn; kj++)
		{
		   tm = sbder[kjh*irend+krh];
		   for (kp=0; kp<idim; kjhc++, kp++)
		     ssum[kp] += shc[kjhc]*tm;
		}
	      
	      tm = (double)1.0/sbder[(krh+km)*irend+krh];
	      for (kp=0; kp<idim; khc++, kp++)
		shc[khc] = (sder[krh*idim+kp] - ssum[kp])*tm;
	   }
	 khc -= 2*idim;
      }


   /* Compute the remaining unknowns by straight forward elimination. */
   
   for (krh=ik-km, kstart=in-ik+km-1, kstop=in-irend, kr=kstart;
    kr>=kstop; krh++, kr--)
     {
	memzero(ssum,idim,DOUBLE);
	for (kjhc=khc+idim, kjh=ik-krh, kj=kr+1; kj<in; kjh++, kj++)
	  {
	     tm = sbder[kjh*irend+krh];
	     for (kp=0; kp<idim; kjhc++, kp++)
	       ssum[kp] += shc[kjhc]*tm;
	  }
	
	tm = (double)1.0/sbder[(ik-krh-1)*irend+krh];
	for (kp=0; kp<idim; khc++, kp++)
	  shc[khc] = (sder[krh*idim+kp] - ssum[kp])*tm;
	khc -= 2*idim;
     }
   
   /* Put the last irend coefficients into ec, but set ec equal to the
      average of the two solutions where they overlap.  */
   
   for (khc=0, ki=in-irend; ki<in; ki++)
     for (kp=0; kp<idim; khc++, kp++)
       {
	  if (ki < ilend) th = ec[ki*idim+kp];
	  else th = shc[khc];
	  ec[ki*idim+kp] = (th + shc[khc])/(double)2.0;
       }
   }
   
   *jstat = 0;
   goto out;
   
   /* Error in space allocation.  */
   
   err101: *jstat = -101;
   goto out;

   /* Conflicting dimensions. */
   
   err103: *jstat = -103;
   goto out;
   
   /* Conflicting orders of curves. */
   
   err106: *jstat = -106;
   goto out;
   
   /* Error in lower level routine.  */
   
   error: *jstat = kstat;
   goto out;
   
   out:
      /* Free scratch used for local array.  */
      
      if (sder != SISL_NULL) freearray(sder);
      if (sbder != SISL_NULL) freearray(sbder);
      if (ssum != SISL_NULL) freearray(ssum);
      if (shc != SISL_NULL) freearray(shc);

      return;
}
   
