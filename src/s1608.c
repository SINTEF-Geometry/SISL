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
 * $Id: s1608.c,v 1.2 2001-03-19 15:58:51 afr Exp $
 *
 */


#define S1608

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1608(SISLCurve *pc1,SISLCurve *pc2,double aepsge,
	   double ep11[],double epf1[],double ep21[],double epf2[],
	   int itype,int idim,int ik,SISLCurve **rc,
	   double *ct11,double *ctf1,double *ct21,double *ctf2, int *jstat)
#else
void s1608(pc1,pc2,aepsge,ep11,epf1,ep21,epf2,itype,idim,ik,rc,
           ct11,ctf1,ct21,ctf2,jstat) 
	   SISLCurve  *pc1;
	   SISLCurve  *pc2;
	   double aepsge;
	   double ep11[];
	   double epf1[];
	   double ep21[];
	   double epf2[];
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
* PURPOSE    : To calculate a fillet curve between two curves.
*              Points indicate between which points on the input curve 
*              the fillet is to be produced.
*
* INPUT      : pc1    - The first input curve.   
*              pc2    - The second input curve.   
*              aepsge - Geometry resolution.         
*              ep11   - SISLPoint close to curve 1 telling that the part of the 
*                       curve lying on this side of epf1 is not to be 
*                       replaced by the fillet.
*              epf1   - SISLPoint close to curve 1, indicating where the fillet is
*                       to start. The tangent at the start of the fillet will
*                       have the same orientation as the curve from ep11 
*                       to epf1.
*              ep21   - SISLPoint close to curve 2 telling that the part of the
*                       curve lying on this side of epf2 is not to be 
*                       replaced by the fillet.
*              epf2   - SISLPoint close to curve two, indicating where the fillet
*                       is to end. the tangent at the end of the fillet will
*                       have the same orientation as the curve from epf2
*                       to ep21.
*              itype  - Indicator of type of fillet.
*                     = 1  - Circle, interpolating tangent on first curve,
*                            not on curve 2.
*                     = 2  - Conic if possible
*                     else - Polynomial segment
*              idim   - Dimension of space.  
*              ik     - Order of fillet curve.
*
* OUTPUT     : rc     - Fillet curve produced
*              ct11   - parameter value of point ep11 on curve 1.
*              ctf1   - parameter value of point epf1 on curve 1.
*              ct21   - parameter value of point ep21 on curve 2.
*              ctf2   - parameter value of point epf1 on curve 2.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                         =-1      : No fillet produced
*
* METHOD     : First the parameter values  at11, atf1, at21, atf2
*              corresponding respectively to ep11, epf1, ep21 and epf2
*              are found. Then S1607 is used for calculating the actual fillet.
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
  SISLIntcurve **qic1=SISL_NULL;  /* SISLObject containing intervals if any      */
  SISLIntcurve **qic2=SISL_NULL;  /* SISLObject containing intervals if any      */
  SISLIntcurve **qic3=SISL_NULL;  /* SISLObject containing intervals if any      */
  SISLIntcurve **qic4=SISL_NULL;  /* SISLObject containing intervals if any      */
  
  int kstat;          /* Status variable                                   */
  int kpos=0;         /* Position of error                                 */
  
  int kcrv1,kcrv2;
  int kcrv3,kcrv4;    /* Number of intervals                               */
  int kpt1,kpt2,kpt3; /* Number of points                                  */
  int kpt4;           /* Loop variable                                     */
  double *spar1=SISL_NULL; /* Pointer to parameter values                       */
  double *spar2=SISL_NULL; /* Pointer to parameter values                       */
  double *spar3=SISL_NULL; /* Pointer to parameter values                       */
  double *spar4=SISL_NULL; /* Pointer to parameter values                       */
  
  /* Check dimensions */
  
  if (idim != 2 && idim != 3) goto err105;
  if (pc1->idim != pc2->idim) goto err106;
  
  /* Check if curves are  correct */
  
  s1707(pc1,&kstat);
  if (kstat < 0) goto error;
  
  s1707(pc2,&kstat);
  if (kstat < 0) goto error;
  
  /* Calculate closest point to ep11 */
  
  s1953(pc1,ep11,idim,REL_COMP_RES,aepsge,&kpt1,&spar1,&kcrv1,&qic1,&kstat);
  if (kstat < 0) goto error;
  
  /* Remember closest point */
  
  if (kpt1  > 0)
    *ct11 = spar1[0];
  else if (kcrv1>0)
    {
      SISLIntcurve *q1= *qic1;
      if (q1->ipar1 ==1)
        *ct11 = q1 -> epar1[0];
      else if (q1->ipar2 ==1)
        *ct11 = q1 -> epar2[0];
      else
        goto errxxx;
    }
  
  /* Calculate closest point to epf1 */
  
  s1953(pc1,epf1,idim,REL_COMP_RES,aepsge,&kpt2,&spar2,&kcrv2,&qic2,&kstat);
  if (kstat < 0) goto error;
  
  /* Remember closest point */
  
  if (kpt2  > 0)
    *ctf1 = spar2[0];
  else if (kcrv2>0)
    {
      SISLIntcurve *q2= *qic2;
      if (q2->ipar1 ==1)
        *ctf1 = q2 -> epar1[0];
      else if (q2->ipar2 ==1)
        *ctf1 = q2 -> epar2[0];
      else
        goto errxxx;
    }
  
  /* Calculate closest point to ep21 */
  
  s1953(pc2,ep21,idim,REL_COMP_RES,aepsge,&kpt3,&spar3,&kcrv3,&qic3,&kstat);
  if (kstat < 0) goto error;
  
  /* Remember closest point */
  
  if (kpt3  > 0)
    *ct21 = spar3[0];
  else if (kcrv3>0)
    {
      SISLIntcurve *q3= *qic3;
      if (q3->ipar1 ==1)
        *ct21 = q3 -> epar1[0];
      else if (q3->ipar2 ==1)
        *ct21 = q3 -> epar2[0];
      else
        goto errxxx;
    }
  
  /* Calculate closest point to epf2 */
  
  s1953(pc2,epf2,idim,REL_COMP_RES,aepsge,&kpt4,&spar4,&kcrv4,&qic4,&kstat);
  if (kstat < 0) goto error;
  
  /* Remember closest point */
  
  if (kpt4 > 0)
    *ctf2 = spar4[0];
  else if (kcrv1>0)
    {
      SISLIntcurve *q4= *qic4;
      if (q4->ipar1 ==1)
        *ctf2 = q4 -> epar1[0];
      else if (q4->ipar2 ==1)
        *ctf2 = q4 -> epar2[0];
      else
        goto errxxx;
    }
  
  s1607(pc1,pc2,aepsge,*ct11,*ctf1,*ct21,*ctf2,itype,idim,ik,rc,&kstat);
  if (kstat<0) goto error;
  
  *jstat = 0;
  
  goto out;
  
  /* Error in memory allocation */
  
  /* Error in input, conflicting dimensions */
  
 err106: *jstat = -106;
  s6err("s1608",*jstat,kpos);
  goto out;
  
  /* Dimension nmot equal to 2 or 3 */
  
 err105: *jstat = -105;
  s6err("s1608",*jstat,kpos);
  goto out;
  
  
  /*      No fillet produced */
  
 errxxx: *jstat = -1;
  goto out;
  
  /* Error in lower level function */  
  
 error:  *jstat = kstat;
  s6err("s1608",*jstat,kpos); 
  goto out;
  
 out:
  if (qic1  != SISL_NULL) freeIntcrvlist(qic1,kcrv1);
  if (qic2  != SISL_NULL) freeIntcrvlist(qic2,kcrv2);
  if (qic3  != SISL_NULL) freeIntcrvlist(qic3,kcrv3);
  if (qic4  != SISL_NULL) freeIntcrvlist(qic4,kcrv4);
  if (spar1 != SISL_NULL) freearray(spar1);
  if (spar2 != SISL_NULL) freearray(spar2);
  if (spar3 != SISL_NULL) freearray(spar3);
  if (spar4 != SISL_NULL) freearray(spar4);
  
  return;
}       
