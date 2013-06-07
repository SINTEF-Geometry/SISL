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
 * $Id: s1450.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1450

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1450(SISLSurf *ps1,double aepsge,int *jclos1,int *jclos2,
	   int *jdgen1,int *jdgen2,int *jdgen3,int *jdgen4,int *jstat)
#else
void s1450(ps1,aepsge,jclos1,jclos2,jdgen1,jdgen2,jdgen3,jdgen4,jstat)
     SISLSurf   *ps1;
     double aepsge;
     int    *jclos1;
     int    *jclos2;
     int    *jdgen1;
     int    *jdgen2;
     int    *jdgen3;
     int    *jdgen4;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To decide if a surface is closed or degenerate along
*              its boundaries.
*
* INPUT      : ps1    - Pointer to the surface to be checked
*              aepsge - Tolerance for testing
*
* OUTPUT     : jclos1 - Closed indicator in first parameter direction
*                        0 - Surface open if first parameter direction
*                        1 - Surface closed if first parameter direction
*            : jclos2 - Closed indicator in second parameter direction
*                        0 - Surface open if second parameter direction
*                        1 - Surface closed if second parameter direction
*              jdgen1 - Degenerate indicator along standard edge 1
*                       (v=et2[ik2-1])
*                        0 - SISLEdge is not degenerate
*                        1 - SISLEdge is degenerate
*              jdgen2 - Degenerate indicator along standard edge 2
*                       (u=et1[in1]
*                        0 - SISLEdge is not degenerate
*                        1 - SISLEdge is degenerate
*              jdgen3 - Degenerate indicator along standard edge 3
*                       (v=et2[in2]
*                        0 - SISLEdge is not degenerate
*                        1 - SISLEdge is degenerate
*              jdgen4 - Degenerate indicator along standard edge 4
*                       (u=et1[ik1-1])
*                        0 - SISLEdge is not degenerate
*                        1 - SISLEdge is degenerate
*
*              jstat  - status messages  
*                                         > 0      : Warning
*                                         = 0      : ok
*                                         < 0      : error
*                      
*
* METHOD     : All standard edges are picked out and then tested
*              to find if the surface is closed in one of the parameter
*              directions or if it degenerate along some edge.
*
* USE:        int kclos1,kclos2,kdgen1,kdgen2,kdgen3,kdgen4,kstat;
*              double tepsge;
*              SISLSurf *qs1;
*               .
*               .
*              s1450(qs1,tepsge,&kclos1,&kclos2,&kdgen1,&kdgen2,&kdgen3,
*                    &kdgen4,&kstat);
*               .
*
*
* REFERENCES :
*
*-
* CALLS      : s1436, s1437, s1451
*
* WRITTEN BY : Tor Dokken, SI, Norway, 1988-11
* REVISED BY : Johannes Kaasa, SI, Sep 1991 (Introduced NURBS)
*
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* Position of error.                              */
  int kdim;           /* Dimension of the space in which the surface lies. */
  int rdim;           /* Dimension of the rational space.                */
  int kn1;            /* Number of vertices in first parameter direction */
  int kn2;            /* Number of vertices in first parameter direction */
  int kk1;            /* Order in first parameter driection              */
  int kk2;            /* Order in first parameter driection              */
  int ki;             /* Control variables in for loop                   */
  double *scoef;      /* B-spline vertices of surface                    */
  double *st1,*st2;   /* Pointers to knot vectors                        */
  double *sv1,*sv2,*sv3,*sv4;      /* Pointers to vertices               */
  SISLCurve *q1=SISL_NULL;    /* Pointer to boundary curve        */
  SISLCurve *q2=SISL_NULL;    /* Pointer to boundary curve        */
  SISLCurve *q3=SISL_NULL;    /* Pointer to boundary curve        */
  SISLCurve *q4=SISL_NULL;    /* Pointer to boundary curve        */
  
  
  if (aepsge < (double)0.0) goto err184;
  
  kn1 = ps1->in1;
  kn2 = ps1->in2;
  kk1 = ps1->ik1;
  kk2 = ps1->ik2;
  kdim = ps1 ->idim;
  scoef = ps1 -> ecoef;
  st1   = ps1 -> et1;
  st2   = ps1 -> et2;
  if (ps1->ikind == 2 || ps1->ikind == 4)
    rdim = kdim + 1;
  else
    rdim = kdim;
  
  /* Initiate output variables */
  
  *jclos1 = 1;
  *jclos2 = 1;
  *jdgen1 = 1;
  *jdgen2 = 1;
  *jdgen3 = 1;
  *jdgen4 = 1;
  
  /* Pick standard edge 1 */
  
  s1436(ps1,st2[kk2-1],&q1,&kstat);
  if (kstat < 0) goto error;
  
  /* Pick standard edge 2 */
  
  s1437(ps1,st1[kn1],&q2,&kstat);
  if (kstat < 0) goto error;
  
  /* Pick standard edge 3 */
  
  s1436(ps1,st2[kn2],&q3,&kstat);
  if (kstat < 0) goto error;
  
  
  /* Pick standard edge 4 */
  
  s1437(ps1,st1[kk1-1],&q4,&kstat);
  if (kstat < 0) goto error;
  
  /* Check if close in first parameter direction i.e. if the
     vertices of standard edge 2 and 4 correspond */
  
  if (rdim == kdim)
    {
      sv2 = q2->ecoef;
      sv4 = q4->ecoef;
    }
  else
    {
      sv2 = q2->rcoef;
      sv4 = q4->rcoef;
    }
  for (ki=0; ki<kn2 ; ki++,sv2+=rdim,sv4+=rdim)
    if(s6dist(sv2,sv4,kdim)>aepsge)
      {
        *jclos1 = 0;
        ki = kn2;
      }
  
  /* Check if close in second parameter direction i.e. if the
     vertices of standard edge 1 and 3 correspond */
  
  if (rdim == kdim)
    {
      sv1 = q1->ecoef;
      sv3 = q3->ecoef;
    }
  else
    {
      sv1 = q1->rcoef;
      sv3 = q3->rcoef;
    }
  for (ki=0; ki<kn1; ki++,sv1+=rdim,sv3+=rdim)
    if(s6dist(sv1,sv3,kdim)>aepsge)
      {
        *jclos2 = 0;
        ki = kn1;
      }
  
  /* Check if standard edge 1 degenerate */
  
  s1451(q1,aepsge,jdgen1,&kstat);
  if (kstat<0) goto error;
  
  /* Check if standard edge 2 degenerate */
  
  s1451(q2,aepsge,jdgen2,&kstat);
  if (kstat<0) goto error;
  
  /* Check if standard edge 3 degenerate */
  
  s1451(q3,aepsge,jdgen3,&kstat);
  if (kstat<0) goto error;
  
  /* Check if standard edge 4 degenerate */
  
  s1451(q4,aepsge,jdgen4,&kstat);
  if (kstat<0) goto error;
  
  
  
  *jstat = 0;
  
  goto out;
  
  /* Negative absolute tolerance.   */
  
 err184: *jstat = -184;
  s6err("s1450",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1450",*jstat,kpos);
  goto out;
  
 out:
    if (q1 != SISL_NULL) freeCurve(q1);
    if (q2 != SISL_NULL) freeCurve(q2);
    if (q3 != SISL_NULL) freeCurve(q3);
    if (q4 != SISL_NULL) freeCurve(q4);
  return;
}

                     
