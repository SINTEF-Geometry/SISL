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
 * $Id: s1609.c,v 1.2 2001-03-19 15:58:51 afr Exp $
 *
 */


#define S1609

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1609(SISLCurve *pc1,SISLCurve *pc2,double aepsge,
	   double eps1[],double epf[],double eps2[],double aradius,
	   double enorm[],int itype,int idim,int ik,SISLCurve **rc,
	   double *ct11,double *ctf1,double *ct21,double *ctf2,
	   int *jstat)
#else
void s1609(pc1,pc2,aepsge,eps1,epf,eps2,aradius,enorm,itype,idim,ik,rc,
           ct11,ctf1,ct21,ctf2,jstat)
     SISLCurve  *pc1;
     SISLCurve  *pc2;
     double aepsge;
     double eps1[];
     double epf[];
     double eps2[];
     double aradius;
     double enorm[];
     int    itype;
     int    idim;
     int    ik;
     SISLCurve  **rc;
     double *ct11;
     double *ctf1;
     double *ct21;
     double *ctf2;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate a fillet curve given radius between two curves.
*
* INPUT      : pc1    - The first input curve.   
*              pc2    - The second input curve.   
*              aepsge - Geometry resolution.         
*              eps1   - SISLPoint telling that the fillet should be put on 
*                       the side of curve 1 where eps1 is situated.
*              epf    - SISLPoint indicating where the fillet curve should go. 
*                       eps1 together with epf indicates the direction of 
*                       the start tangent of the curve, while epf together
*                       with eps2 indicates the direction of the end tangent
*                       of the curve. If more than one position of the fillet 
*                       curve is possible, the closest curve to epf is chosen.
*              eps2   - SISLPoint telling that the fillet should be put on the 
*                       side of curve 2 where eps2 is situated.
*              aradius - The radius to be used on the fillet if a circular
*                       fillet possible, else a conic or a quadratic 
*                       polynomial curve is used, approximating the circular
*                       fillet.
*              enorm  - Direction normal to the plane the fillet curve 
*                       should lie close to. (only used in 3-d fillet 
*                       calculations).
*              itype  - Indicator of type of fillet.
*                     = 1  - Circle, interpolating tangent on first curve,
*                            not on curve 2.
*                     = 2  - Conic if possible
*                     else - Polynomial segment
*              idim   - Dimension of space.  
*              ik     - Order of fillet curve.
*
* OUTPUT     : rc     - Fillet curve produced
*              ct11   - Parameter value of the end of curve 1 not affected by 
*                       the fillet.
*              ctf1   - Parameter value of the point on curve 1 where the 
*                       fillet starts.
*              ct21   - Parameter value of the end of curve 2 not affected by 
*                       the fillet.
*              ctf2   - Parameter value of the point on curve 2 where the 
*                       fillet ends.
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                        
* METHOD     : 
*
* USE        : 
*             
* REFERENCES :
*                   
*-                                                 
* CALLS      : 
*
*
* WRITTEN BY : Qyvind Hjelle, SI, Oslo, Norway. 28. Nov 1988
* Reviced by : Tor Dokken, SI, Oslo, Norway, August 1989
*              
*********************************************************************
*/
{
  SISLIntcurve **qic1=SISL_NULL;  /* SISLObject containing intervals if any     */
  SISLIntcurve **qic2=SISL_NULL;  /* SISLObject containing intervals if any     */
  SISLIntcurve **qic3=SISL_NULL;  /* SISLObject containing intervals if any     */
  SISLCurve    *qc1=SISL_NULL;    /* Offset of first curve                  */
  SISLCurve    *qc2=SISL_NULL;    /* Offset of second curve                 */
  
  int kstat;          /* Status variable                                  */
  int kpos=0;         /* Position of error                                */
  int kleft=0;        /* Pointer in knot vector                           */
  int kn;             /* The number of B-splines vertices                 */
  int kk;             /* The polynomial order of the curve.               */
  double *st;         /* Pointer to the first element of the knot vector
			 of the curve. The knot vector has [kn+kk]
			 elements.                                        */
  
  int kcrv1,kcrv2;
  int kcrv3    ;      /* Number of intervals                              */
  int kpt1,kpt2,kpt3; /* Number of points                                 */
  int ki;             /* Loop variable                                    */
  double tpar1,tpar2; /* Parameter values                                 */
  double spnt1[6];    /* SISLPoint and derivate                               */
  double spnt2[6];    /* SISLPoint and derivate                               */
  double *spar1=SISL_NULL; /* Pointer to parameter values                      */
  double *spar2=SISL_NULL; /* Pointer to parameter values                      */
  double *spar3=SISL_NULL; /* Pointer to parameter values                      */
  double *spar4=SISL_NULL; /* Pointer to parameter values                      */
  double snorm[3];    /* Normal vector                                    */
  double sdiff[3];    /* Difference vector                                */
  double tprod;       /* Scalar product                                   */
  double tdir1,tdir2; /* Direction of curves                              */
  double trad1,trad2; /* Offsets with correct direction                   */
  double tmin;        /* Minimum distance                                 */
  double tdist;       /* Distance                                         */
  
  /* Check dimensions */
  
  if (idim != 2 && idim != 3) goto err105;
  if (pc1->idim != pc2->idim) goto err106;
  
  /* Check if curves are  correct */
  
  s1707(pc1,&kstat);
  if (kstat < 0) goto error;
  
  s1707(pc2,&kstat);
  if (kstat < 0) goto error;
  
  /* Calculate closest point to eps1 */
  
  s1953(pc1,eps1,idim,REL_COMP_RES,aepsge,&kpt1,&spar1,&kcrv1,&qic1,&kstat);
  if (kstat < 0) goto error;
  
  /* Remember closest point */
  
  if (kpt1  > 0)
    tpar1 = spar1[0];
  else if (kcrv1>0)
    {
      SISLIntcurve *q1= *qic1;
      if (q1->ipar1 ==1)
        tpar1 = q1 -> epar1[0];
      else if (q1->ipar2 ==1)
        tpar1 = q1 -> epar2[0];
      else
        goto errxxx;
    }
  
  /* Calculate point and tangent on curve at tpar1 */
  
  s1221(pc1,1,tpar1,&kleft,spnt1,&kstat);
  if (kstat < 0) goto error;
  
  
  /* Calculate closest point to eps2 */
  
  s1953(pc2,eps2,idim,REL_COMP_RES,aepsge,&kpt2,&spar2,&kcrv2,&qic2,&kstat);
  if (kstat < 0) goto error;
  
  
  /* Remember closest point */
  
  if (kpt2  > 0)
    tpar2 = spar2[0];
  else if (kcrv1>0)
    {
      SISLIntcurve *q1= *qic2;
      if (q1->ipar1 ==1)
        tpar2 = q1 -> epar1[0];
      else if (q1->ipar2 ==1)
        tpar2 = q1 -> epar2[0];
      else
        goto errxxx;
    }
  
  /* Calculate point and tangent on curve at tpar1 */
  
  s1221(pc2,1,tpar2,&kleft,spnt2,&kstat);
  if (kstat<0) goto error;
  
  
  /* Determine if offset to be used on curve 1 is positive or negative */
  
  if (idim==2)
    {
      snorm[0] = -spnt1[3];
      snorm[1] =  spnt1[2];
      s6diff(eps1,spnt1,idim,sdiff);
      tprod = s6scpr(snorm,sdiff,idim);
      s6diff(epf,spnt1,idim,sdiff);
      tdir1 = s6scpr(sdiff,spnt1+idim,idim);
    }
  else
    {
      s6diff(epf,spnt1,idim,sdiff);
      tdir1 = s6scpr(spnt1+idim,sdiff,idim);
      s6crss(sdiff,spnt1+idim,snorm);
      tprod = s6scpr(snorm,enorm,idim);
    }
  
  /* Set direction of offset of curve 1 */
  
  if (tprod >= DZERO)
    trad1 = aradius;
  else
    trad1 = -aradius;
  
  /* Determine if offset to be used on curve 1 is positive or negative */
  
  if (idim==2)
    {
      snorm[0] = -spnt2[3];
      snorm[1] =  spnt2[2];
      s6diff(eps2,spnt2,idim,sdiff);
      tprod = s6scpr(snorm,sdiff,idim);
      s6diff(epf,spnt2,idim,sdiff);
      tdir2 = s6scpr(sdiff,spnt2+idim,idim);
    }
  else
    {
      s6diff(epf,spnt2,idim,sdiff);
      tdir2 = s6scpr(spnt2+idim,sdiff,idim);
      s6crss(sdiff,spnt2+idim,snorm);
      tprod = s6scpr(snorm,enorm,idim);
    }
  
  
  /* Determine if offset to be used on curve 1 is positive or negative */
  
  if (tprod >= DZERO)
    trad2 = aradius;
  else
    trad2 = -aradius;
  
  /* Make curve offset from curve 1 */
  
  s1360(pc1,trad1,aepsge,enorm,(double)0.0,idim,&qc1,&kstat);
  if (kstat<0) goto error;
  
  
  /* Make curve offset from curve 2 */
  
  s1360(pc2,trad2,aepsge,enorm,(double)0.0,idim,&qc2,&kstat);
  if (kstat<0) goto error;
  
  /* Find closest point between the two curves, 
     if intersection does not succeed
     use closest point between curves */
  
  s1857(qc1,qc2,REL_COMP_RES,aepsge,&kpt3,&spar3,&spar4,&kcrv3,&qic3,&kstat);
  if (kstat<0) goto error;
  
  if (kpt3 == 0 && kcrv3 == 0)
    {
      s1955(qc1,qc2,REL_COMP_RES,aepsge,
	    &kpt3,&spar3,&spar4,&kcrv3,&qic3,&kstat);
      if (kstat<0) goto error;
    }
  
  /* Find point closest to ept */
  
  tmin=HUGE;
  
  if (kpt3>0)
    {
      /*  Find intersection point closest to epf */
      
      for (ki=0;ki<kpt3;ki++)
        {
	  s1221(pc1,0,spar3[ki],&kleft,spnt1,&kstat);
	  if (kstat<0) goto error;
	  
	  tdist = s6dist(spnt1,epf,idim);
	  if (tdist<=tmin)
            {
	      *ctf1 = spar3[ki];
	      *ctf2 = spar4[ki];
	      tmin = tdist;
            }
        }
    }
  else if (kcrv3>0)
    {
      SISLIntcurve *q1= *qic3;
      for (ki=0; ki < q1->ipoint ; ki++)
        {
	  s1221(pc1,0,q1->epar1[ki],&kleft,spnt1,&kstat);
	  if (kstat<0) goto error;
	  tdist = s6dist(spnt1,epf,idim);
	  if (tdist<=tmin)
            {
	      *ctf1 = q1->epar1[ki];
	      *ctf2 = q1->epar2[ki];
	      tmin = tdist;
            }
        }
    }
  
  if (tdist == HUGE) goto errxxx;
  
  /* Set parameter values indicating which part of curves remains after
     the fillet */
  
  st = pc1->et;
  kk = pc1->ik;
  kn = pc1->in;
  
  if (tdir1>=0)
    *ct11 = st[kk-1];
  else
    *ct11 = st[kn];
  
  /* Set parameter values indicating which part of curve two remains after
     the fillet */
  
  st = pc2->et;
  kk = pc2->ik;
  kn = pc2->in;
  
  if (tdir2>=0)
    *ct21 = st[kk-1];
  else
    *ct21 = st[kn];
  
  s1607(pc1,pc2,aepsge,*ct11,*ctf1,*ct21,*ctf2,itype,idim,ik,rc,&kstat);
  if (kstat<0) goto error;
  
  *jstat = 0;
  
  goto out;
  
  
  /* Error in input, conflicting dimensions */
  
 err106: 
  *jstat = -106;
  s6err("s1609",*jstat,kpos);
  goto out;
  
  /* Dimension nmot equal to 2 or 3 */
  
 err105: 
  *jstat = -105;
  s6err("s1609",*jstat,kpos);
  goto out;

  
  /* No fillet produced */
  
 errxxx: 
  *jstat = -1;
  goto out;
  
  /* Error in lower level function */  
  
 error:  
  *jstat = kstat;
  s6err("s1609",*jstat,kpos); 
  goto out;
  
 out:
  if (qc1   != SISL_NULL) freeCurve(qc1);
  if (qc2   != SISL_NULL) freeCurve(qc2);
  if (qic1  != SISL_NULL) freeIntcrvlist(qic1,kcrv1);
  if (qic2  != SISL_NULL) freeIntcrvlist(qic2,kcrv2);
  if (qic3  != SISL_NULL) freeIntcrvlist(qic3,kcrv3);
  if (spar1 != SISL_NULL) freearray(spar1);
  if (spar2 != SISL_NULL) freearray(spar2);
  if (spar3 != SISL_NULL) freearray(spar3);
  if (spar4 != SISL_NULL) freearray(spar4);
  
  return;
}       
