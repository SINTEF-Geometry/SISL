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
 * $Id: s1304.c,v 1.2 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1304

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1304(double ep[],double eq[],double eparp[],double eparq[],double egeo3d[],
	   double egeop[],double egeoq[],int *jstat)
#else
void s1304(ep,eq,eparp,eparq,egeo3d,egeop,egeoq,jstat)
     double ep[];
     double eq[];
     double eparp[];
     double eparq[];
     double egeo3d[];
     double egeop[];
     double egeoq[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the radius of curvature at the intersection
*              point between two curves. The first and second derivatives
*              of the intersection point in the two surfaces are give as
*              input. All coordinates are assumed to be in 3-D.
*
*
*
* INPUT      : ep      - 0-2 order derivatives of first surface.
*                        The sequence is position, first derivative in first
*                        parameter direction, first derivative in second
*                        parameter direction, (2,0) derivative, (1,1)
*                        derivative, (0,2) derivative and normal. (21 numbers)
*                        Compatible with output of s1421
*              eq      - 0-2 order derivatives of second surface. Same
*                        sequence as ep.
*              eparp   - Parameter pair in first surface of point
*              eparq   - Parameter pair in second surface of point
*
* OUTPUT     : 
*              jstat  - status messages  
*                         = 1      : Curvature radius infinit
*                         = 0      : ok, curvature radius
*                         < 0      : error
*              egeo3d - 3-D geometry description of the intersection. The
*                       array contains: position, unit tangent, curvature
*                       and radius of curvature. (A total of 10 numbers)
*                       A radius of curvature =-1, indicates that the radius
*                       of curvature is infinit.
*              egeop  - Description of the intersection in the parameter plane
*                       of the first surface. The array contains: position,
*                       unit tangent, curvature and radius of curvature.
*                       (A total of 7 numbers)
*              egeoq  - Description of the intersection in the parameter plane
*                       of the second surface. The array contains: position,
*                       unit tangent, curvature and radius of curvature.
*                       (A total of 7 numbers)
*
* METHOD     : First the most lineary independent selection of 3 vectors
*              from the derivative vectors are found. Then equation systems
*              are made to determine the derivatives of the parameter values
*              with respect to the parameter value not present in the
*              selection of the three vectors. Then the double derivatives
*              of the parameter direction are found. This information is used
*              for expressing the 3-D tangent, curvature and the radius
*              of curvature of the intersection point in question.
*              Corresponding values in both parameter planes are also found.
*
*
* REFERENCES : 
*-
* CALLS      : s6norm,s6scpr,sqrt,fabs,s6length,s6crss
*              
*
* WRITTEN BY : Tor Dokken, SI, Oslo , Norway, 30 May 1988
*
*********************************************************************
*/
{
  int kdim = 3;           /* Dimension of 3-D space                      */
  int k2dim = 2;          /* Dimension of the parameter planes           */
  int kstat = 0;          /* Local status variable                       */
  int ki;                 /* Control variable in loop                    */
  double snpu[3];         /* Nomalized version of epu                    */
  double snpv[3];         /* Nomalized version of epv                    */
  double spn[3];          /* Vector snpu x snpv                          */
  double snqs[3];         /* Nomalized version of eqs                    */
  double snqt[3];         /* Nomalized version of eqt                    */
  double sqn[3];          /* Vector snqs x snqt                          */
  double sright[3];       /* Right hand side when finding s"             */
  double sdc[3];          /* Derivative of intersection curve by w       */
  double sddc[3];         /* Second Derivative of intersection curve by w*/
  double *sqs;            /* Pointer to first row of matrix              */
  double *sqt;            /* Pointer to second row of matrix             */
  double *spu;            /* Pointer to third row of matrix              */
  double *spw;            /* Pointer to fourth row of matrix             */
  double *spuu;           /* Pointer to renamed (2,0)-derivative         */
  double *spuw;           /* Pointer to renamed (1,1)-derivative         */
  double *spww;           /* Pointer to renamed (0,2)-derivative         */
  double *sqss;           /* Pointer to renamed (2,0)-derivative         */
  double *sqst;           /* Pointer to renamed (1,1)-derivative         */
  double *sqtt;           /* Pointer to renamed (0,2)-derivative         */
  double tt;              /* Value of det(snpu,snpv,snqs)                */
  double ts;              /* Value of det(snpu,snpv,snqt)                */
  double tv;              /* Value of det(snqs,snqt,snpu)                */
  double tu;              /* Value of det(snqs,snqt,snpv)                */
  double tlpu;            /* Length of epu                               */
  double tlpv;            /* Length of epv                               */
  double tlqs;            /* Length of eqs                               */
  double tlqt;            /* Length of eqt                               */
  double tmax1;           /* Variable used for maximal value             */
  double tmax2;           /* Variable used for maximal value             */
  double tdum;            /* Dummy variable                              */
  double tdom;            /* The denominator in an equation              */
  double tds;             /* ds/dw                                       */
  double tdt;             /* dt/dw                                       */
  double tdu;             /* du/dw                                       */
  double tddu;            /* ddu/dwdw                                    */
  double tdds;            /* dds/dwdw                                    */
  double tddt;            /* ddt/dwdw                                    */
  double twds;            /* ds/dw after renaming variable second time   */
  double twdt;            /* dt/dw after renaming variable second time   */
  double twdds;           /* dds/dwdw after renaming variable second time*/
  double twddt;           /* ddt/dwdw after renaming variable second time*/
  double twdu;            /* du/dw after renaming variable second time   */
  double twdv;            /* dv/dw after renaming variable second time   */
  double twddu;           /* ddu/dwdw after renaming variable second time*/
  double twddv;           /* ddv/dwdw after renaming variable second time*/
  
  
  /* Get input values into output. */
  egeoq[0] = eparq[0];
  egeoq[1] = eparq[1];
  egeop[0] = eparp[0];
  egeop[1] = eparp[1];
  
  for (ki=2;ki<7;ki++)
    {
      egeop[ki] = DZERO;
      egeoq[ki] = DZERO;
    }
  
  /* Make position of intersection */
  
  for (ki=0;ki<3;ki++)
    {
      egeo3d[ki] = (double)0.5 * (ep[ki]+eq[ki]);
    }
  
  for (ki=3;ki<10;ki++) egeo3d[ki] = DZERO;
  
  
  /* Nomalize derivative vectors */
  
  tlpu = s6norm(ep+3,kdim,snpu,&kstat);
  tlpv = s6norm(ep+6,kdim,snpv,&kstat);
  tlqs = s6norm(eq+3,kdim,snqs,&kstat);
  tlqt = s6norm(eq+6,kdim,snqt,&kstat);
  
  /* Make normal vector for both derivative pairs */
  
  s6crss(snpu,snpv,spn);
  s6crss(snqs,snqt,sqn);                                 
  
  /* Make four scalar product to decide which of the 3 vectors snpu, snpv, 
   * snqs, snqt spans the 3-D space best. (Have the biggest determinant)
   * Remember (axb)c = det(a,b,c). The naming convention is that the
   * name of the variable not present on the left hand side is used for
   * the naming of the determinants. The determinants tt, ts, tu and tv tells
   * which direction is most linearly dependent on the other directions  
   */
  
  tt = fabs(s6scpr(spn,snqs,kdim));
  ts = fabs(s6scpr(spn,snqt,kdim));
  tv = fabs(s6scpr(sqn,snpu,kdim));
  tu = fabs(s6scpr(sqn,snpv,kdim));
  
  /* We want to use the parameter direction names s, t and u on the left
   * hand side of the equation system, thus we want to express all derivatives
   * of the curve as functions of a new parameter  value w, which is chosen
   * to be the parameter direction with partial first derivative most
   * lineary dependent of the other parameter directions 
   */
  
  tmax1 = MAX(tt,ts);
  tmax2 = MAX(tv,tu);
  
  if (tmax1 > tmax2)
    {
      
      /*  The s or t variable should not be used on the left hand side of
       *  the equation system                                                
       */
      
      if (ts>tt) 
        {
	  
	  /*  The s variable should not be used on the left hand side of the 
	   *  equation system
	   *  The renaming of variables is as follows s->w, t->u, u->s, v->t 
	   */
	  
	  spu  = eq+6;
	  spw  = eq+3;
	  spuu = eq+15;
	  spuw = eq+12;
	  spww = eq+9;
	  sqs  = ep+3;
	  sqt  = ep+6;
	  sqss = ep+9;
	  sqst = ep+12;
	  sqtt = ep+15;
        }
      else
        {
	  
	  /* The t variable should not be used on the left hand side of the 
	   * equation system
	   *  The renaming of variables is as follows s->u, t->w, u->s, v->t 
	   */
	  
	  spu  = eq+3;
	  spw  = eq+6;
	  spuu = eq+9;
	  spuw = eq+12;
	  spww = eq+15;
	  sqs  = ep+3;
	  sqt  = ep+6;
	  sqss = ep+9;
	  sqst = ep+12;
	  sqtt = ep+15;
        }
    }
  else
    {
      
      /* The u or v variable should not be used on the left hand side of
       * the equation system                                                
       */
      
      if (tu>tv)
	
        {
	  
	  /* The u variable should not be used on the left hand side of the
	   * equation system.
	   * The renaming of variables is as follows s->s, t->t, u->w, v->u
	   */
	  
	  spu  = ep+6;
	  spw  = ep+3;
	  spuu = ep+15;
	  spuw = ep+12;
	  spww = ep+9;
	  sqs  = eq+3;
	  sqt  = eq+6;
	  sqss = eq+9;
	  sqst = eq+12;
	  sqtt = eq+15;
	  
        }
      
      else
	
        {
	  /* The v variable should not be used on the left hand side of the
	   * equation system.
	   * The renaming of variables is as follows s->s, t->t, u->u, v->w 
	   */
	  
	  spu  = ep+3;
	  spw  = ep+6;
	  spuu = ep+9;
	  spuw = ep+12;
	  spww = ep+15;
	  sqs  = eq+3;
	  sqt  = eq+6;
	  sqss = eq+9;
	  sqst = eq+12;
	  sqtt = eq+15;
	  
        }
    }
  
  /* Now we can solve the equation systems for finding
   *  ds/dw, dt/dw and du/dw and afterwards for
   *  dds/(dwdw), ddt/(dwdw) and ddu/(dwdw), using Cramers Rule.
   *
   *  This equation system is derived in the following way:
   *  Our problem is defined as P(u,w) - Q(s,t) = 0. By taking the derivative
   *  of this equation with repsect to w, we get:
   *
   *  dP(u,w) du   dP(u,w)   dQ(s,t) ds   dQ(s,t) dt
   *  ------- -- + ------- - ------- -- - ------- -- = 0
   *  du      dw   dw        ds      dw   dt      dw    
   *
   *  By using a simplified notation this can be written:
   *
   *  P u' + P  - Q s' - Q t' = 0
   *   u      w    s      t
   *
   *  We can thus set up the equation system:
   *
   *              s'
   *  (Q  Q -P ) (t') = P
   *    s t   u   u'     w
   * 
   *
   *
   *  By making one futher derivative we get an equation systen for s",t" and
   *  u".
   *
   *               s"         2                       2                   2
   *  (Q  Q  -P ) (t") = P  u'  + 2P  u' + P   - Q  s'  - 2Q  s't' - Q  t' 
   *    s  t   u   u"     uu        uw      ww    ss        st        tt
   *
   *
   */
  
  /* Prepare normal vectors for determinants */
  
  s6crss(spu,spw,spn);
  s6crss(sqs,sqt,sqn);
  
  tdom =  -s6scpr(sqn,spu,kdim);
  
  if (DEQUAL(tdom,(double)0.0)) goto war101;
  
  /*  Lineary dependent vectors on left hand side if tdom = 0.0 */
  
  tds = -s6scpr(spn,sqt,kdim)/tdom;
  tdt =  s6scpr(spn,sqs,kdim)/tdom;
  tdu =  s6scpr(sqn,spw,kdim)/tdom;
  
  for (ki=0;ki<3;ki++)
    {
      sright[ki] = (spuu[ki]*tdu + (double)2.0*spuw[ki])*tdu + spww[ki]
	           - (sqss[ki]*tds + sqst[ki]*tdt)*tds
	           - (sqtt[ki]*tdt + sqst[ki]*tds)*tdt;
    }
  
  /* Calculate second derivatives of parameter direction with respect to
   * the w-direction 
   */
  
  tddu = s6scpr(sqn,sright,kdim)/tdom;
  
  /* Use sqn for temporary storage of cross products */
  
  s6crss(sright,sqt,sqn);
  tdds = -s6scpr(sqn,spu,kdim)/tdom;
  s6crss(sqs,sright,sqn);
  tddt = -s6scpr(sqn,spu,kdim)/tdom;
  
  /* We will now express the intersection curve locally as a function
   *  of the w-parameter.
   *
   *  c(w) = p(u(w),w)
   *
   *  This gives the derivative
   *                         
   *   
   *  c' = P u' + P                                                   
   *        u      w
   *
   *  And the second derivative
   *
   *            2
   *  c" = P  u'  + 2P  u' + P   + P u"
   *        uu        uw      ww    u
   *
   *  The curvature vector is defined as the derivative of the unit tangent
   *  vector with respect to the arc length a(w):
   *
   *         d         d    c'(w)    dw   d    c'(w)      da
   *  k(a) = -- T(a) = -- ---------- -- = -- ---------- / --
   *         da        dw sqrt(c'c') da   dw sqrt(c'c')   dw
   *
   *
   *         d       c'(w)                c"        c' (c'c'')
   *         -- ----------------- =   ---------- - ------------- 
   *         dw sqrt(c'(w) c'(w))     sqrt(c'c')   sqrt(c'c')**3
   *
   *
   *
   *         da
   *         -- = sqrt(c'c')
   *         dw 
   */
  for (ki=0;ki<3;ki++)
    {
      sdc[ki] = spu[ki]*tdu + spw[ki];                
      sddc[ki] = (spuu[ki]*tdu+(double)2.0*spuw[ki])*tdu +
	spu[ki]*tddu + spww[ki];
    }
  
  /* To simplify futher calculations we want to normalize the tangent vector
   *  and correspondingly divide the second derivative by the tangent length
   */
  
  tlpu = s6norm(sdc,kdim,egeo3d+3,&kstat);
  
  if (DEQUAL(tlpu,(double)0.0)) goto war101;
  
  for (ki=0;ki<3;ki++)
    {
      sddc[ki] = sddc[ki]/tlpu;
    }
  
  /* Make curvature vector */
  
  tdum = s6scpr(sddc,egeo3d+3,kdim);
  for (ki=0;ki<3;ki++)
    {
      egeo3d[ki+6] = (sddc[ki] - egeo3d[ki+3]*tdum)/tlpu;
    }
  
  /* TO CALCULATE UNIT TANGENT CURVATURE AND RADIUS OF CURVATURE OF THE
   * INTERSECTION POINT IN THE PARAMETER PLANES OF THE TWO SURFACES, WE
   * NOW WANT TO CALCULATE THE TRUE VALUES OF ds/dw, dt/dw, dds/dw,
   * ddt/dwdw, du/dw, dv/dw, ddu/dwdw, ddv/dwdw, where w is the parameter
   * direction we have chosen the other directions to be expressed in.
   * THUS UNDO the changing of parameter directions 
   */
  
  if (tmax1 > tmax2)
    {
      
      /* First and second row of the surface were originally interchanged
       * Thus change sequence of these back again 
       */ 
      
      if (ts>tt) 
        {
	  /* We used the following renaming of variables:
	   * s->w, t->u, u->s, v->t, now express the behavior in the parameter
	   * plane with the original ordering 
	   */
	  
	  twds  = (double)1.0;
	  twdt  = tdu;
	  twdds = (double)0.0;
	  twddt = tddu;
	  twdu  = tds;
	  twdv  = tdt;
	  twddu = tdds;
	  twddv = tddt;
	  
        }
      else
        {
	  
	  /* We used the following renaming of variables:
	   * s->u, t->w, u->s, v->t, now express the behavior in the parameter
	   * plane with the original ordering 
	   */
	  
	  /* The renaming of variables is as follows s->v, t->u, u->s, v->t */
	  
	  twds  = tdu;
	  twdt  = (double)1.0;
	  twdds = tddu;
	  twddt = (double)0.0;
	  twdu  = tds;
	  twdv  = tdt;
	  twddu = tdds;
	  twddv = tddt;
	  
        }
    }
  else
    {
      
      /*  Keep the sequence of surfaces */
      
      if (tu>tv)
	
        {
	  
	  /* We used the following renaming of variables:
	   * s->s, t->t, u->w, v->u, now express the behavior in the parameter
	   * plane with the original ordering 
	   */
	  
	  twds  = tds;
	  twdt  = tdt;
	  twdds = tdds;
	  twddt = tddt;
	  twdu  = (double)1.0;
	  twdv  = tdu;
	  twddu = (double)0.0;
	  twddv = tddu;
	  
        }                  
      
      else
	
        {
	  
	  /* We used the following renaming of variables:
	   * s->s, t->t, u->u, v->w, now express the behavior in the parameter
	   * plane with the original ordering */
	  
	  twds  = tds;
	  twdt  = tdt;
	  twdds = tdds;
	  twddt = tddt;
	  twdu  = tdu;
	  twdv  = (double)1.0;
	  twddu = tddu;
	  twddv = (double)0.0;
	  
        }
    }
  
  /* Now the variable twds, twdt, twdu, twdv contains derivatives of the
   * parameter directions with respect to the w-variable. Correspondingly
   * the second derivatives with respect to w are contained in twdds, twddt,
   * twddu and twddv.
   *
   * THE UNIT TANGENT, CURVATURE VECTOR AND RADIUS OF CURVATURE CAN NOW
   * BE CALCULATED IN BOTH PARAMETER PLANES                                
   */
  
  /* Make description of intersection curve in 
   * parameter plane of first patch
   */
  
  
  tdum = sqrt(twdu*twdu + twdv*twdv);
  if (DEQUAL(tdum,(double)0.0))
    {
      egeop[2] = (double)0.0;
      egeop[3] = (double)0.0;
      egeop[4] = (double)0.0;
      egeop[5] = (double)0.0;
      egeop[6] = (double)0.0;
    }
  else
    {
      
      /* Make unit tangent        */
      
      egeop[2] = twdu/tdum;
      egeop[3] = twdv/tdum;
      
      /* Make curvature vector    */
      
      tdom = egeop[2]*twddu + egeop[3]*twddv;
      egeop[4] = (twddu/tdum - egeop[2]*tdom/tdum)/tdum;
      egeop[5] = (twddv/tdum - egeop[3]*tdom/tdum)/tdum;
    }
  
  /* Make radius of curvature in parameter plane 1 */
  
  tdum = s6length(egeop+4,k2dim,&kstat);
  if (DNEQUAL(tdum,(double)0.0))
    {
      egeop[6] = (double)1.0/tdum;
    }
  else
    {
      egeop[6] = (double)-1.0;
    }
  
  /* Make description of intersection curve in parameter 
   * plane of second patch
   */
  tdum = sqrt(twds*twds + twdt*twdt);
  if (DEQUAL(tdum,(double)0.0))
    {
      egeoq[2] = (double)0.0;
      egeoq[3] = (double)0.0;
      egeoq[4] = (double)0.0;
      egeoq[5] = (double)0.0;
      egeoq[6] = (double)0.0;
    }
  else
    {
      
      /* Make unit tangent  */
      
      egeoq[2] = twds/tdum;
      egeoq[3] = twdt/tdum;
      
      /* Make curvature vector */
      
      tdom     = egeoq[2]*twdds + egeoq[3]*twddt;
      egeoq[4] = (twdds/tdum - egeoq[2]*tdom/tdum)/tdum;
      egeoq[5] = (twddt/tdum - egeoq[3]*tdom/tdum)/tdum;
    }
  
  /* Make radius of curvature in parameter plane 2 */
  
  tdum = s6length(egeoq+4,k2dim,&kstat);
  if (DNEQUAL(tdum,(double)0.0))
    {
      egeoq[6] = (double)1.0/tdum;
    }
  else
    {
      egeoq[6] = (double)-1.0;
    }
  
  /* Make 3-D radius of curvature */
  
  tdum = s6length(egeo3d+6,kdim,&kstat);
  
  if (DNEQUAL(tdum,(double)0.0))
    {
      egeo3d[9] = (double)1.0/tdum;
    }
  else
    {
      egeo3d[9] = (double)-1.0;
      goto war101;
    }
  
  
  *jstat = 0;
  goto out;
  
  /* Infinit radius of curvature */
  
 war101: *jstat=1;
  goto out;
  
 out:
  return;
}
                    
