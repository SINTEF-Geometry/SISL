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
 * $Id: s1384.c,v 1.3 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1384

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1384(SISLCurve *pcurve,SISLSurf *psurf,int idim,int iside,double ax,
	   int *ileftc,int *ilefts1,int *ilefts2,double eder[],
	   double ederu[],double ederv[],double edern[],int *jstat)
#else
void s1384(pcurve,psurf,idim,iside,ax,ileftc,ilefts1,ilefts2,
           eder,ederu,ederv,edern,jstat)
     SISLCurve  *pcurve;
     SISLSurf   *psurf;
     int    idim;
     int    iside;
     double ax;
     int    *ileftc;
     int    *ilefts1;
     int    *ilefts2;
     double eder[];
     double ederu[];
     double ederv[];
     double edern[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the value and first and second derivatives of
*              a curve traced out by a curve in the parameter plane
*              at the point with parameter value ax. In addition to 
*              calculate posision and first derivative of dp(u(s),v(s))/du,
*              ddp(u(s),v(s))/duds, dp(u(s),v(s))/dv, ddp(u(s),v(s))/dvds
*              and similar values for the normal vector to the trace curve
*
*
*
*
* INPUT      : pcurve - Pointer to the curve in the parameter plane
*              psurf  - The surface from which the curve is traced out
*              idim   - Dimension of the space the curve lies in 2 or 3
*              iside  - Calculate derivative from right or left
*                        ileft ==  -1 from left
*                        ileft!=  -1 from right
*              ax     - The parameter value at which to compute
*                       position and derivatives.
*
*                
*
* INPUT/OUTPUT : ileftc- Pointer to the interval in the knot vector
*                        where ax is located. If et is the knot vector,
*                        the relation
*                          
*                          et[ileft] <= ax < et[ileft+1]
* 
*                        should hold. (If ax == et[in] then ileft should
*                        be in-1. Here in is the number of B-spline
*                        coefficients.)
*                        If ileft does not have the right value upon
*                        entry to the routine, its value will be changed
*                        to the value satisfying the above condition.
*              ilefts1- Similar pointer to first knot vector of surface 
*              ilefts2- Similar pointer to second knot vector of surface 
*
*
* OUTPUT     : eder   - Double array of dimension [3*idim]
*                       containing the position and derivative vectors.
*                       (idim is the number of components of each B-spline
*                       coefficient, i.e. the dimension of the Euclidean
*                       space in which the curve lies.)
*                       These vectors are stored in the following order:
*                       First the idim components of the position vector,
*                       then the idim components of the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
*                       (The C declaration of eder as a two dimensional array
*                       would therefore be eder[kder+1,idim].)
*
*              ederu  - Double array of dimension [2*idim] containing
*                        dp(u(s),v(s))/du and ddp(u(s),v(s))/duds.
*              ederv  - Double array of dimension [2*idim] containing
*                       dp(u(s),v(s))/dv, ddp(u(s),v(s))/dvds
*              edern  - Double array of dimension [2*idim] conatining
*                       the value and derivative of the normal to the
*                       tangent curve. NOT IMPLEMENTED
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : First the point and relevant derivatives are calculated on
*              the point in the parameter plane, then relevant derivatives
*              are calculated in the surface and the information is combined
*              to produce the specified derivatives.
*
*
*
* REFERENCES :
*
*-
* CALLS      : s1221, s6length, s6crss, s6scpr, s6err
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, October 1988
* REVISED BY : Christophe rene Birkeland, SINTEF Oslo, May 1993
*              Testing for SISL_NULL (array sders and sderc)
*
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* The position of the error.                      */
  int kder=3;         /* Derivatives necessary to calculate              */
  int ki;             /* Loop variable                                   */
  int knum1,knum2;    /* Numbers of points and derivatives               */
  int kdimc;          /* Dimension of curve                              */
  int kdims;          /* Dimension of surface                            */
  double sdumc[8];    /* Values from curve calculation                   */
  double sdums[30];   /* Values from surface calcualtion                 */
  double snorm[3];    /* Normal to tangent vector                        */
  double *sders=SISL_NULL; /* Values and derivatives on surface               */
  double *sderc=SISL_NULL; /* Values and derivatives on curve                 */
  double tduds,tdvds; /* Derivatives of curve                            */
  double tddudss,tddvdss; /* Derivatives of curve                        */
  double *sdpdu,*sdpdv;/* Pointers to derivatives of surface             */
  double *sddpduu;    /* Pointers to derivatives of surface              */
  double *sddpduv;    /* Pointers to derivatives of surface              */
  double *sddpdvv;    /* Pointers to derivatives of surface              */
  double *ps;         /* Pointer in loop                                 */
  double tlsn;        /* Length of sn                                    */
  
  kdimc = pcurve -> idim;
  kdims = psurf  -> idim;
  
  if (kdimc !=2 && kdims != idim) goto err105;
  
  knum1 = kdimc*(kder+1);
  knum2 = kdims*(kder+1)*(kder+2)/2;
  
  if (knum1>8)
  {
    if((sderc = newarray(knum1,DOUBLE)) == SISL_NULL) goto err101;
  }
  else  
    sderc = sdumc;
  
  if (knum2>30)
  {
    if((sders = newarray(knum2,DOUBLE)) == SISL_NULL) goto err101;
  }
  else   
    sders = sdums;
  
  /* Calculate point values and derivatives on curve */ 
  
  
  if (iside != -1)
    {
      /*  Calculated derivatives from the right */
      s1221(pcurve,kder,ax,ileftc,sderc,&kstat);
      if (kstat<0) goto error;
    }
  else
    {
      /*  Calculated derivatives from the left */
      s1227(pcurve,kder,ax,ileftc,sderc,&kstat);
      if (kstat<0) goto error;
    }
  
  /* Calculate position and derivatives of the surface */ 
  
  s1421(psurf,kder,sderc,ilefts1,ilefts2,sders,snorm,&kstat);            
  if (kstat<0) goto error;
  
  /* Make some pointers that facilitate the futher calculations */
  
  tduds  = sderc[2];
  tdvds  = sderc[3];
  tddudss = sderc[4];
  tddvdss = sderc[5];
  sdpdu   = sders + idim;
  sdpdv   = sdpdu + idim;
  sddpduu = sdpdv + idim;
  sddpduv = sddpduu + idim;
  sddpdvv = sddpduv + idim;
  
  /* FIRST MAKE CALCULATION OF POSITION, DERIVATIVE AND SECOND DERIVATIVE
     OF TRACE CURVE */
  
  /* Calculate position of point */
  
  for (ki=0;ki<idim;ki++) eder[ki] = sders[ki];
  
  /* Calculate first derivatives */
  
  ps = eder + idim;
  for (ki=0;ki<idim;ki++)
    *(ps++) = sdpdu[ki]*tduds + sdpdv[ki]*tdvds;
  
  /*  Calculate second derivatives of */
  
  for (ki=0;ki<idim;ki++)
    *(ps++) = sddpduu[ki]*tduds*tduds + sdpdu[ki]*tddudss +
      sddpdvv[ki]*tdvds*tdvds + sdpdv[ki]*tddvdss +
	(double)2.0*sddpduv[ki]*tduds*tdvds;
  
  /* MAKE CALCULATION OF VALUE (POSITION) AND DERIVATIVES OF dp/du WITH
     RESPECT TO THE TRACE CURVE */
  
  ps = ederu;
  for (ki=0;ki<idim;ki++)
    *(ps++) = sdpdu[ki];
  
  for (ki=0;ki<idim;ki++)
    *(ps++) = sddpduu[ki]*tduds + sddpduv[ki]*tdvds;
  
  
  /* MAKE CALCULATION OF VALUE (POSITION) AND DERIVATIVES OF dp/dv WITH
     RESPECT TO THE TRACE CURVE */
  
  ps = ederv;
  for (ki=0;ki<idim;ki++)
    *(ps++) = sdpdv[ki];
  
  for (ki=0;ki<idim;ki++)
    *(ps++) = sddpduv[ki]*tduds + sddpdvv[ki]*tdvds;
  
  
  /* MAKE POSITION AND DERIVATIVE OF THE NORMAL VECTOR:
     
     NOT IMPLEMENTED..................
     d                 d
     ( -- p(u(s),v(s)) x -- p(u(s),v(s)) )
     dp(u(s),v(s))     du                dv
     ------------- x ----------------------------------
     ds              length of numerator
     
     */
  
  /* First make the normal vector as cross product of the u and v derivative,
     the length of the normal vector and the two components of the derivative
     with respect to s.
     */
  
  s6crss(ederu,ederv,snorm);
  
  tlsn = s6length(snorm,3,&kstat);
  
  if (kstat<0) goto error;
  
  if (DEQUAL(tlsn,DZERO)) tlsn = (double)1.0;
  
  s6crss(eder+idim,snorm,edern);
  
  ps = edern;
  
  for (ki=0;ki<idim;ki++,ps++)
    *(ps) /= tlsn;
  
  /*
    
    To calculate the derivative of the normal we maks:
    
    d                 d
    ( -- p(u(s),v(s)) x -- p(u(s),v(s)) )
    d    dp(u(s),v(s))   du                dv
    --  (------------- x ---------------------------------- ) =
    ds       ds              length of numerator
    
    
    2             
    d p(u(s),v(s))   sn      dp(u(s),v(s))    d            
    -------------- x ---   + ------------- x  --sn/lsn -
    dsds         lsn         ds           ds
    
    
    d p(u(s),v(s))   sn       d
    -------------- x ---- (sn --sn)
    ds              3     ds
    lsn
    
    
    d                 d 
    sn   =   -- p(u(s),v(s)) x -- p(u(s),v(s))  
    du                dv
    
    
    lsn  =   length of(sn)
    
    
    2
    d                   d
    sn1  =   --  p(u(s),v(s))) x -- p(u(s),v(s)) )
    duds                dv
    
    2                   2
    d                   d
    sn2  =   --  p(u(s),v(s))) x -- p(u(s),v(s)) )
    du                  dvds
    
    
    
    d   
    --sn =   sn1 + sn2
    ds
    
    
    Make the components of the derivative 
    
    s6crss(ederu+idim,ederv,sn1);
    s6crss(ederu,ederv+idim,sn2);
    
    for (ki=0;ki<3;ki++)
    sdnds[ki] = sn1[ki] + sn2[ki];
    
    
    Make the cross products of the three components 
    
    s6crss(eder+6,sn,sdum1);
    s6crss(eder+3,sdnds,sdum2);
    s6crss(eder+3,sn,sdum3);
    
    Make the necessary factor for the last component 
    
    tfak = s6scpr(sn,sdnds,3);
    
    Make the first derivative of the normal vector
    
    ps= edern+3;
    
    tlsn3 = tlsn*tlsn*tlsn;
    
    for (ki=0;ki<idim;ki++)
    *(ps++) = (sdum1[ki] + sdum2[ki])/tlsn + sdum3[ki]*tfak/tlsn3;
    
    
    END OF NOT IMPLEMENTED */
  
  /* Successful computations.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in allocations */

  err101: *jstat = -101;
    s6err("s1384",*jstat,kpos);
    goto out;
  
  /* idim not 2 0r 3 */

  err105: *jstat = -105;
    s6err("s1384",*jstat,kpos);
    goto out;
  
  /* Error in lower level routine.  */
  
  error:  *jstat = kstat;
    s6err("s1384",*jstat,kpos);
    goto out;
  
  out:
  
    if (knum1>8)
      { 
        if (sderc != SISL_NULL) freearray(sderc);
      }
    if (knum2>30) 
      {
        if (sders != SISL_NULL) freearray(sders);
      }
    return;
}
