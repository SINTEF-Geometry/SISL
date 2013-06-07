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
 * $Id: s1331.c,v 1.2 2001-03-19 15:58:45 afr Exp $
 *
 */


#define S1331

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1331(double ep[],double eimpli[],int ideg,int ider,
	   double gder[],double gnorm[],int *jstat)
#else
void s1331(ep,eimpli,ideg,ider,gder,gnorm,jstat)
     double ep[];
     double eimpli[];
     int    ideg;
     int    ider;
     double gder[];
     double gnorm[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the position and derivatives up to second
*              order of a point (and corresponding derivatives) put into
*              the equation of an implicit surface.
*
*   
* INPUT      : ep     - 0-2 order derivatives of B-spline surface.
*                       For ideg=1,2 and 1001 the sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2) derivative
*                       and normal. (21 numbers)
*                       For ideg=1003,1004,1005 the second derivatives are followed
*                       by the third derivatives and the normal (33 numbers)
*                       Compatible with output of s1421
*              eimpli - Description of the implicit surface
*              ideg   - The degree of the implicit surface
*                        ideg=1:    Plane              
*                        ideg=2;    Quadric surface
*                        ideg=1001: Torus surface
*                        ideg=1003: Silhouette line parallel projection
*                        ideg=1004: Silhouette line perspective projection
*                        ideg=1005: Silhouette line circular projection
*              ider   - Derivatives to be calculated
*                        ider=-1:   Normal
*                        ider=0:    Value + normal
*                        ider=1:    Value + first derivatives + normal
*                        ider=2:    Value + 1.st + 2.nd + normal
*                       Note if ideg=1003,1004,1005 then ider is set to max(ider,1).
*
*
* OUTPUT     : 
*              jstat  - status messages
*                         = 0      : ok, curvature radius
*                         < 0      : error
*
*              gder   - Specified derivatives in order: value, 
*                       1.st derivatives and 2.nd derivatives
*              gnorm  - Normal of the implicit surface.
*                       When ideg=1003,1004,1005 i.e. for silhouette curves
*                       gnorm = f P - f P  which is the direction of
*                                v u   u v
*                       the silhouette curve along the surface P.
*
* METHOD 
*
* USE:         This function is only working for 3-D input
*
* REFERENCES : 
*-
* CALLS      : s6scpr,fabs,s1307,s6err
*              
*
* WRITTEN BY : Tor Dokken, SI, Oslo , Norway, 10 October 1988
* REVISED BY : Mike Floater, SI, Oslo , Norway, 31 January 1991 
*               Added perspective silhouette (ideg=1004) and
*                     circular silhouette (ideg=1005).
*               Introduced a proper useable definition of gnorm for silhouettes.
*               It's a little different from the usual cases. See s1331
*               in order to see how it is used.
*
*********************************************************************
*/
{
  int kstat = 0;          /* Local status variable                       */
  int ki,kj,kl;           /* Control variables in loop                   */
  int kpos = 0;           /* Position of error                           */
  int ksize;              /* Number of doubles for storage of derivateves
			     and normal vector */
  int ksizem3;            /* ksize - 3                               */
  double *spu;            /* Pointer to dP/ds                            */
  double *spv;            /* Pointer to dP/dt                            */
  double *spuu;           /* Pointer to ddP/(dsds)                       */
  double *spuv;           /* Pointer to ddP/(dsdt)                       */
  double *spvv;           /* Pointer to ddP/(dtdt)                       */
  
  
  /* If ideg=1,2 or 1001 then only derivatives up to second order
     are calculated, then 18 doubles for derivatives and 3 for the
     normal vector are to be used for calculation of points in the
     spline surface. For ideg=1003,1004,1005 we have a silhouette curve and
     derivatives up to the third are to be calculated,
     thus 30 +3 a total of 33 doubles may be calculated.
     Note, gnorm in these cases is replaced by f P - f P
                                                v u   u v
     so we need ideg to be at least 1. */
  
  if (ideg == 1003 || ideg == 1004 || ideg == 1005)
    {
      ksize = 33;
      ider=max(ider,1);
    }
  else
    {
      ksize = 21;
    }

  ksizem3 = ksize -3;
  
  /* The calculation of the derivatives of one parameter direction with
   *  respect to the other is dependent on the degree of the implicit equation
   */
  
  /* Set local pointers */
  
  spu = ep + 3;
  spv = ep + 6;
  spuu = ep + 9;
  spuv = ep + 12;
  spvv = ep + 15;
  
  if (ideg == 1)
    {
      /*  First degree implicit geometry. 
       *
       *   Let A = (eimpli[0],eimpli[1],eimpli[2],eimpli[3]) and
       *       Q= (P(s,t),1), 
       *   then putting  Q into the implicit equation gives
       *
       *    f(u,v) = A Q
       *
       *    df         dQ
       *    --     = A -- = N P , where N is the normal vector of the plane.
       *    du         du      u
       *
       *    df         dQ
       *    --     = A -- = N P 
       *    dv         dv      v
       *
       *     2          2
       *    d f        d Q
       *    --     = A --- = N P 
       *      2          2      uu
       *    du         du     
       *
       *     2          2
       *    d f        d Q
       *    --     = A --- = N P 
       *    dudv       dudv     uv
       *
       *     2          2
       *    d f        d Q
       *    --     = A --- = N P 
       *      2          2      vv
       *    dv         dv     
       *
       */
      if (ider>=0)
        {
	  gder[0] = s6scpr(ep,eimpli,3) + eimpli[3];
	  if (ider>=1)
            {
	      gder[1] = s6scpr(spu,eimpli,3);
	      gder[2] = s6scpr(spv,eimpli,3);
	      if (ider>=2)
                {
		  gder[3] = s6scpr(spuu,eimpli,3); 
		  gder[4] = s6scpr(spuv,eimpli,3);
		  gder[5] = s6scpr(spvv,eimpli,3); 
                }
            }
        }
      memcopy(gnorm,eimpli,3,DOUBLE);
    }
  else if (ideg == 2)
    {   
      double sq[4],sduq[4],sdvq[4];    /*vectors*/
      double tsum1,tsum2,tsum3;        /*Accumulation of matrix product */ 
      
     /*  Second degree implicit geometry.
     *   Denote the 4x4 matrix representing the implicit surface a A.*
     *
     *   We can now calculate the u and v derivatives of the point put into
     *   the left hand side of the implicit surface equation:
     *
     *            T              T
     *      f(u,v) = Q A Q = sq Q 
     *
     *      d        d        T      dQ    T        dQ
     *      - f    = -- (Q A Q ) = 2 -- A Q  = 2 sq --
     *      du       du              du             du
     *
     *      d        d        T      dQ    T        dQ
     *      - f    = -- (Q A Q ) = 2 -- A Q  = 2 sq --
     *      dv       dv              dv             dv
     *
     *       2        2               2                  T        2
     *      d        d        T      d Q    T     dQ   dQ       dQ        dQ
     *      -- f   = -- (Q A Q ) = 2 --- A Q  + 2 -- A -- = 2(sq--- + sduq--)
     *        2        2               2          du   du         2       du
     *      du       du              du                         du
     *                                                          
     *       2        2                 2               T        2
     *      d        d          T      d Q    T    dQ dQ        d Q        dQ
     *      --  f  = --   (Q A Q ) = 2 ----AQ  + 2 --A--  = 2(sq---- + sduq--)
     *      dudv     dudv              dudv        du dv        dudv       dv
     *
     *       2        2               2                  T        2
     *      d        d        T      d Q    T     dQ   dQ       dQ        dQ
     *      -- f   = -- (Q A Q ) = 2 --- A Q  + 2 -- A -- = 2(sq--  + sdvq--)
     *        2        2               2          dv   dv         2       dv
     *      dv       dv              dv                         dv
     *
     *                                       dQ             dQ
     *      First calculate sq = QA, sduq =  --A and sdvq = --A
     *                                      du             dv
     *
     */
      for (ki=0;ki<4;ki++)
        {
	  tsum1 = eimpli[ki+12];
	  tsum2 = DZERO;
	  tsum3 = DZERO;
	  for (kj=0,kl=ki;kj<3;kj++,kl+=4)
            {
	      tsum1 += eimpli[kl]*ep[kj];
	      if (ider>=1)
                {
		  tsum2 += eimpli[kl]*spu[kj];
		  tsum3 += eimpli[kl]*spv[kj];
                }
            }
	  sq[ki]   = tsum1;
	  sduq[ki] = tsum2;
	  sdvq[ki] = tsum3;
        }
      
      /*  Make value and partial derivatives */
      
      if (ider>=0)
        {
	  gder[0] = s6scpr(sq,ep,3) + sq[3];
	  if(ider>=1)
            {
	      gder[1] = (double)2.0*s6scpr(sq,spu,3);
	      gder[2] = (double)2.0*s6scpr(sq,spv,3);
	      if (ider>1)
                {
		  gder[3] = (double)2.0*(s6scpr(sq,spuu,3) + 
					 s6scpr(sduq,spu,3));
		  gder[4] = (double)2.0*(s6scpr(sq,spuv,3) + 
					 s6scpr(sduq,spv,3));
		  gder[5] = (double)2.0*(s6scpr(sq,spvv,3) + 
					 s6scpr(sdvq,spv,3));
                }
            }
        }
      memcopy(gnorm,sq,3,DOUBLE);
    }
  
  else if (ideg == 1001)
    {
      double sy[3],sz[3],szu[3],szv[3],szuu[3],szuv[3],szvv[3];/* Derivatives */
      double tzn,tzun,tzvn,tzuun,tzuvn,tzvvn;            /* Scalar products */
      double tzduz,tzdvz,tzduduz,tzdudvz,tzdvdvz;        /* Scalar products */
      double tduzduz,tdvzdvz,tduzdvz;                    /* Scalar products */
      double tlenz,tlenz3,tlenz5;                        /* Vector lengths  */
      double sp[3],sdup[3],sdvp[3];                      /* temporary vectrs*/ 
      double sdudup[3],sdudvp[3],sdvdvp[3];              /* temporary vectrs*/ 
      double *scentr,*saxis,tbigr,tsmalr;                /* Torus descript  */

      /*  Torus surface. */
      
      scentr = eimpli;
      saxis  = eimpli + 3;
      tbigr  = *(eimpli+6);
      tsmalr = *(eimpli+7);
      
      
/*     Intersection of implicit function of 1 and implicit torus surface
*      and the the surface. Make the following temporary variables.
*
*      y = p - scentr
*      z = y - (y saxis) saxis
*
*      Put the point and accompanying derivatives into the implicit
*      representation of the torus
*
*                                   2    2
*      f(u,v) = (y - R z/sqrt(z z) )  - r
* 
*
*                                                                
*                                    R z          (z z )
*      df            R z                u      R z    u
*      -- = 2(y - ---------)(y  -  --------- + ----------)
*      du         sqrt(z z)   u    sqrt(z z)            3
*                                              sqrt(z z)
*
*         = 2 sp sdup
*
*
*                                 R z            (z z ) 
*      df            R z             v        R z    v
*      -- = 2(y - ---------)(y  - --------- + ----------)
*      dv         sqrt(z z)   v   sqrt(z z)            3
*                                             sqrt(z z)
*
*         = 2 sp sdvp
*
*
*
*                                                                
*       2             R z            (z z )
*      d f               u      R z      u  2
*      --  = 2(y  - --------- + -----------)  + 
*        2      u   sqrt(z z)             3
*      du                        sqrt(z z)
*
*
*                                    R z           z (z z )
*                     R z               uu      2R  u    u
*            2(y - ---------)(y   - --------- + ----------- +
*                  sqrt(z z)   uu   sqrt(z z)             3
*                                                sqrt(z z)
*
*
*                                                       2
*                            R z(z z +zz  )       z(zz )
*                                 u u   uu      R     u
*                            -------------- - 3 --------- )
*                                        3               5
*                               sqrt(z z)       sqrt(z z)
*
*          = 2 sdup sdup + sp sdudup 
*
*                                                                
*       2               R z           (z z )          R z           (z z )
*      d  f                u      R z     u              v      R z     v   
*      ----  = 2(y  - --------- + -----------)(y  - --------- + -----------) +
*                 u   sqrt(z z)             3   v   sqrt(z z)             3
*      dudv                        sqrt(z z)                    sqrt(z z) 
*
*
*                                   R z           z (z z )       z (z z )
*                    R z               uv       R  u    v      R  v    u
*           2(y - ---------)(y   - --------- + ------------ + ------------ +
*                 sqrt(z z)   uv   sqrt(z z)             3              3
*                                               sqrt(z z)      sqrt(z z) 
*
*
*                                                        
*                            R z(z z +zz  )       z(zz )(zz )
*                                 v u   uv      R     u    v
*                            -------------- - 3 -------------)
*                                        3                 5
*                               sqrt(z z)         sqrt(z z)
*
*
*          = 2 sdup sdvp + sp sdudvp 
*
*                                                                
*       2             R z           (z z )
*      d f               v      R z     v   2
*      --  = 2(y  - --------- + -----------)  + 
*        2      v   sqrt(z z)             3
*      dv                        sqrt(z z)
*
*
*                                    R z            z (z z )
*                     R z               vv       2R  v    v
*            2(y - ---------)(y   - --------- + ------------ +
*                  sqrt(z z)   vv   sqrt(z z)             3
*                                                sqrt(z z)
*
*
*                                                       2
*                             R z(z z +zz  )       z(zz )
*                                  v v   vv      R     v
*                             -------------- - 3 --------- )
*                                         3               5
*                                sqrt(z z)       sqrt(z z)
*
*
*          = 2 sdvp sdvp + sp sdvdvp 
*
*
*                          
*      y  = (p - scentr)  = p
*       u               u    u
*
*      y  = (p - scentr)  = p
*       v               v    v
*
*      y  = (p - scentr)  = p
*       uu              uu   uu
*
*      y  = (p - scentr)  = p
*       uv              uv   uv
*
*      y  = (p - scentr)  = p
*       vv              vv   vv
*
*      z  = (y - (y N)N)  = p  - (p N)N          N = saxis (Torus axis)
*       u               u    u     u
*
*      z  = (y - (y N)N)  = p  - (p N)N
*       v               v    v     v
*
*      z  = (y - (y N)N)  = p  - (p  N)N
*       uu              uu   uu    uu
*
*      z  = (y - (y N)N)  = p  - (p  N)N 
*       uv              uv   uv    uv
*
*      z  = (y - (y N)N)  = p  - (p  N)N
*       vv              vv   vv    vv
*
*/
      for (ki=0;ki<3;ki++)
        sy[ki] = ep[ki] - scentr[ki];
      
      tzn   = s6scpr(sy  ,saxis,3);
      tzun  = s6scpr(spu ,saxis,3);
      tzvn  = s6scpr(spv ,saxis,3);
      if (ider>1)
        {
	  tzuun = s6scpr(spuu,saxis,3);
	  tzuvn = s6scpr(spuv,saxis,3);
	  tzvvn = s6scpr(spvv,saxis,3);
        }
      
      /*  Make z and necessary derivatives of z */ 
      for (ki=0;ki<3;ki++)
        {    
	  sz[ki]   = sy[ki]   - tzn*saxis[ki];
	  szu[ki]  = spu[ki]  - tzun*saxis[ki];
	  szv[ki]  = spv[ki]  - tzvn*saxis[ki];
	  if (ider>1)
            {
	      szuu[ki] = spuu[ki] - tzuun*saxis[ki];
	      szuv[ki] = spuv[ki] - tzuvn*saxis[ki];
	      szvv[ki] = spvv[ki] - tzvvn*saxis[ki];
            }
        }
      
      /*  Make a number of necessary scalar products */
      
      tzduz   = s6scpr(sz,szu,3);
      tzdvz   = s6scpr(sz,szv,3);
      tduzduz = s6scpr(szu,szu,3);
      tduzdvz = s6scpr(szu,szv,3);
      tdvzdvz = s6scpr(szv,szv,3);
      if (ider>1)
        {
	  tzduduz = s6scpr(sz,szuu,3);
	  tzdudvz = s6scpr(sz,szuv,3);
	  tzdvdvz = s6scpr(sz,szvv,3);
        }      
      
      /*  Find lengt of sz */
      
      tlenz = s6length(sz,3,&kstat);
      
      if (kstat<0) goto error;
      if (DEQUAL(tlenz,DZERO)) tlenz =(double)1.0;
      tlenz3 = tlenz*tlenz*tlenz;
      tlenz5 = tlenz3*tlenz*tlenz;
      
      /* Make a number of necessary vectors */
      
      for (ki=0;ki<3;ki++)                    
	{
	  sp[ki]   = sy[ki]  - tbigr*sz[ki]/tlenz;
	  if (ider>=1)
	    {
	      sdup[ki] = spu[ki] - tbigr*(szu[ki]/tlenz - sz[ki]*tzduz/tlenz3);
	      sdvp[ki] = spv[ki] - tbigr*(szv[ki]/tlenz - sz[ki]*tzdvz/tlenz3);
	      if (ider>=2)
		{
		  sdudup[ki] = spuu[ki] -
		    tbigr*(szuu[ki]/tlenz -
			   ((double)2.0*szu[ki]*tzduz+
			    sz[ki]*(tduzduz+tzduduz))/tlenz3 +
			   (double)3.0*sz[ki]*tzduz*tzduz/tlenz5);
		  
		  sdudvp[ki] = spuv[ki] -
		    tbigr*(szuv[ki]/tlenz -
			   (szu[ki]*tzdvz+szv[ki]*tzduz+
			    sz[ki]*(tduzdvz+tzdudvz))/tlenz3 +
			   (double)3.0*sz[ki]*tzduz*tzdvz/tlenz5);
		  
		  sdvdvp[ki] = spvv[ki] -
		    tbigr*(szvv[ki]/tlenz -
			   ((double)2.0*szv[ki]*tzdvz+
			    sz[ki]*(tdvzdvz+tzdvdvz))/tlenz3 +
			   (double)3.0*sz[ki]*tzdvz*tzdvz/tlenz5);
		}
	    }
	}
      
      /*  Make the derivatives */
      if (ider>=0)
        {
	  gder[0] = s6scpr(sp,sp,3) - tsmalr*tsmalr;
	  if (ider>=1)
            {
	      gder[1] = (double)2.0*s6scpr(sp,sdup,3);
	      gder[2] = (double)2.0*s6scpr(sp,sdvp,3);
	      if (ider>=2)
                {
		  gder[3] = (double)2.0*(s6scpr(sdup,sdup,3)+
					 s6scpr(sp,sdudup,3)); 
		  gder[4] = (double)2.0*(s6scpr(sdup,sdvp,3)+
					 s6scpr(sp,sdudvp,3));
		  gder[5] = (double)2.0*(s6scpr(sdvp,sdvp,3)+
					 s6scpr(sp,sdvdvp,3));
                }
            }
        }
      memcopy(gnorm,sp,3,DOUBLE);
    }
  
  else if (ideg == 1003)
    
    {
      
     /* Silhouette curve, the first three elements of eimpli describes
	the viewing direction.
	
	The silhouette line is descibed by the implicit equation, when
	Q(u,v) is the description of the surface:
	
	f(u,v) = (Q x Q ) view = 0
	u   v          
	
	f      = (Q  x Q ) view + (Q x Q  ) view
	u         uu   v           u   uv
	
	
	f      = (Q  x Q ) view + (Q x Q  ) view
	v         uv   v           u   vv
	
	
	
	f      = (Q    x Q ) view + 2 (Q  x Q  ) view + (Q  x Q   ) view
	uu        uuu    v             uu   uv           u    uuv
	
	
	f      = (Q    x Q )view + (Q  x Q  )view + (Q  x Q  )view + (Q xQ)view
	uv        uuv    v          uu   vv          uv   uv          u  uvv
	
	
	f      = (Q    x Q ) view + 2 (Q  x Q  ) view + (Q  x Q   ) view
	vv        uvv    v             uv   vv           u    vvv
	
	
	Note that Q  x Q   has zero length.
	uv   uv
      */
      double *spuuu = ep + 18;
      double *spuuv = ep + 21;
      double *spuvv = ep + 24;
      double *spvvv = ep + 27;
      double sdum1[3],sdum2[3],sdum3[3];
      
      if (ider>=0)
        {
	  s6crss(spu,spv,sdum1);
	  gder[0] = s6scpr(sdum1,eimpli,3);
	  if (ider>=1)
            {
	      s6crss(spuu,spv,sdum1);
	      s6crss(spu,spuv,sdum2);
	      gder[1] = s6scpr(sdum1,eimpli,3) + s6scpr(sdum2,eimpli,3);
	      s6crss(spuv,spv,sdum1);
	      s6crss(spu,spvv,sdum2);
	      gder[2] = s6scpr(sdum1,eimpli,3) + s6scpr(sdum2,eimpli,3);
	      if (ider>=2)
                {
		  s6crss(spuuu,spv,  sdum1);
		  s6crss(spuu, spuv, sdum2);
		  s6crss(spu,  spuuv,sdum3);
		  gder[3] = s6scpr(sdum1,eimpli,3) +
		    (double)2.0*s6scpr(sdum2,eimpli,3)
		      + s6scpr(sdum3,eimpli,3);;
		  s6crss(spuuv,spv,  sdum1);
		  s6crss(spuu, spvv, sdum2);
		  s6crss(spu,  spuvv,sdum3);
		  gder[4] = s6scpr(sdum1,eimpli,3) + s6scpr(sdum2,eimpli,3)
		    + s6scpr(sdum3,eimpli,3);;
		  s6crss(spuvv,spv,  sdum1);
		  s6crss(spuv, spvv, sdum2);
		  s6crss(spu,  spvvv,sdum3);
		  gder[5] = s6scpr(sdum1,eimpli,3) +
		    (double)2.0*s6scpr(sdum2,eimpli,3)
		      + s6scpr(sdum3,eimpli,3);;
                }
            }
        }
    }
  
  else if (ideg == 1004)
    
    {
      
     /* Perspective silhouette curve, the first three elements of eimpli
        describe the eye point E.
	
	The silhouette line is descibed by the implicit equation, when
	P(u,v) is the description of the surface and
        N(u,v) = P x P  is the normal to the surface:
                  u   v
	
	f(u,v) = N(u,v) . (P(u,v) - E) = 0
	
        Differentiating N gives:

        N = P  x  P   + P  x  P
         u   uu    v     u     uv

        N = P  x  P   + P  x  P
         v   uv    v     u     vv

        N  = P   x  P   + 2 P  x  P   + P x P
         uu   uuu    v       uu    uv    u   uuv

        N  = P   x  P   +  P  x  P   + P x P
         uv   uuu    v      uu    vv    u   uuv

        N  = P   x  P   + 2 P  x  P   + P x P
         vv   uvv    v       uv    vv    u   vvv



        Differentiating f gives:

	f = N  . (P - E)
         u   u   

	f = N  . (P - E)
         v   v   

	f  = N   . (P - E) + N  . P
         uu   uu              u    u 

	f  = N   . (P - E) + N  . P
         uv   uv              u    v 

	f  = N   . (P - E) + N  . P
         vv   vv              v    v 

	
      */
      double *spuuu = ep + 18;
      double *spuuv = ep + 21;
      double *spuvv = ep + 24;
      double *spvvv = ep + 27;
      double norm[3],normu[3],normv[3],normuu[3],normuv[3],normvv[3];
      double cprod1[3],cprod2[3],cprod3[3];
      double pediff[3];
      

      if (ider>=0)
        {
	  s6diff(ep,eimpli,3,pediff);

	  s6crss(spu,spv,norm);
	  gder[0] = s6scpr(norm,pediff,3);

	  if (ider>=1)
            {
	      s6crss(spuu,spv,cprod1);
	      s6crss(spu,spuv,cprod2);
              for (ki=0;ki<3;ki++)
              {    
                  normu[ki] = cprod1[ki] + cprod2[ki];
              }    
	      s6crss(spuv,spv,cprod1);
	      s6crss(spu,spvv,cprod2);
              for (ki=0;ki<3;ki++)
              {    
                  normv[ki] = cprod1[ki] + cprod2[ki];
              }    

	      gder[1] = s6scpr(normu,pediff,3);
	      gder[2] = s6scpr(normv,pediff,3);

	      if (ider>=2)
                {
		  s6crss(spuuu,spv,cprod1);
		  s6crss(spuu,spuv,cprod2);
		  s6crss(spu,spuuv,cprod3);
		  for (ki=0;ki<3;ki++)
		  {    
		  normuu[ki] = cprod1[ki] + 2.0 * cprod2[ki] + cprod3[ki];
		  }    

           /*    Note that spuv x spuv = 0 in the formula for normuv  */
		  s6crss(spuuv,spv,cprod1);
		  s6crss(spuu,spvv,cprod2);
		  s6crss(spu,spuvv,cprod3);
		  for (ki=0;ki<3;ki++)
		  {    
		  normuv[ki] = cprod1[ki] + cprod2[ki] + cprod3[ki];
		  }    

		  s6crss(spuvv,spv,cprod1);
		  s6crss(spuv,spvv,cprod2);
		  s6crss(spu,spvvv,cprod3);
		  for (ki=0;ki<3;ki++)
		  {    
		  normvv[ki] = cprod1[ki] + 2.0 * cprod2[ki] + cprod3[ki];
		  }    

		  gder[3] = s6scpr(normuu,pediff,3) + s6scpr(normu,spu,3);
		  gder[4] = s6scpr(normuv,pediff,3) + s6scpr(normu,spv,3);
		  gder[5] = s6scpr(normvv,pediff,3) + s6scpr(normv,spv,3);

                }
            }
        }
    }
  
  else if (ideg == 1005)
    
    {
      
     /* Perspective silhouette curve, the first three elements of eimpli
        describe Q, the next three  describe B.
	
	The silhouette line is descibed by the implicit equation, when
	P(u,v) is the description of the surface and
        N(u,v) = P x P  is the normal to the surface:
                  u   v
	
	f(u,v) = N(u,v) x (P(u,v) - Q) . B = 0
	
        Differentiating N gives:

        N = P  x  P   + P  x  P
         u   uu    v     u     uv

        N = P  x  P   + P  x  P
         v   uv    v     u     vv

        N  = P   x  P   + 2 P  x  P   + P x P
         uu   uuu    v       uu    uv    u   uuv

        N  = P   x  P   +  P  x  P   + P x P
         uv   uuu    v      uu    vv    u   uuv

        N  = P   x  P   + 2 P  x  P   + P x P
         vv   uvv    v       uv    vv    u   vvv



        Differentiating f gives:

	f = { N  x (P - Q) + N x P } . B
         u     u                  u

	f = { N  x (P - Q) + N x P } . B
         v     v                  v

	f  = { N   x (P - Q) + 2 N x P  + N x P  } . B
         uu     uu                u   u        uu

	f  = { N   x (P - Q) + N x P  + N x P  + N x P  } . B
         uv     uv              u   v    v   u        uv

	f  = { N   x (P - Q) + 2 N x P  + N x P  } . B
         vv     vv                v   v        vv

	
      */
      double *spuuu = ep + 18;
      double *spuuv = ep + 21;
      double *spuvv = ep + 24;
      double *spvvv = ep + 27;
      double *bvec  = eimpli + 3;
      double norm[3],normu[3],normv[3],normuu[3],normuv[3],normvv[3];
      double cprod1[3],cprod2[3],cprod3[3],cprod4[3];
      double pqdiff[3],sum[3];
      

      if (ider>=0)
        {
	  s6diff(ep,eimpli,3,pqdiff);

	  s6crss(spu,spv,norm);
	  s6crss(norm,pqdiff,cprod1);
	  gder[0] = s6scpr(cprod1,bvec,3);

	  if (ider>=1)
            {
	      s6crss(spuu,spv,cprod1);
	      s6crss(spu,spuv,cprod2);
              for (ki=0;ki<3;ki++)
              {    
                  normu[ki] = cprod1[ki] + cprod2[ki];
              }    
	      s6crss(spuv,spv,cprod1);
	      s6crss(spu,spvv,cprod2);
              for (ki=0;ki<3;ki++)
              {    
                  normv[ki] = cprod1[ki] + cprod2[ki];
              }    

	      s6crss(normu,pqdiff,cprod1);
	      s6crss(norm,spu,cprod2);
              for (ki=0;ki<3;ki++)
              {    
                  sum[ki] = cprod1[ki] + cprod2[ki];
              }    
	      gder[1] = s6scpr(sum,bvec,3);

	      s6crss(normv,pqdiff,cprod1);
	      s6crss(norm,spv,cprod2);
              for (ki=0;ki<3;ki++)
              {    
                  sum[ki] = cprod1[ki] + cprod2[ki];
              }    
	      gder[2] = s6scpr(sum,bvec,3);

	      if (ider>=2)
                {
		  s6crss(spuuu,spv,cprod1);
		  s6crss(spuu,spuv,cprod2);
		  s6crss(spu,spuuv,cprod3);
		  for (ki=0;ki<3;ki++)
		  {    
		  normuu[ki] = cprod1[ki] + 2.0 * cprod2[ki] + cprod3[ki];
		  }    

           /*    Note that spuv x spuv = 0 in the formula for normuv  */
		  s6crss(spuuv,spv,cprod1);
		  s6crss(spuu,spvv,cprod2);
		  s6crss(spu,spuvv,cprod3);
		  for (ki=0;ki<3;ki++)
		  {    
		  normuv[ki] = cprod1[ki] + cprod2[ki] + cprod3[ki];
		  }    

		  s6crss(spuvv,spv,cprod1);
		  s6crss(spuv,spvv,cprod2);
		  s6crss(spu,spvvv,cprod3);
		  for (ki=0;ki<3;ki++)
		  {    
		  normvv[ki] = cprod1[ki] + 2.0 * cprod2[ki] + cprod3[ki];
		  }    

	          s6crss(normuu,pqdiff,cprod1);
	          s6crss(normu,spu,cprod2);
	          s6crss(norm,spuu,cprod3);
                  for (ki=0;ki<3;ki++)
                  {    
                      sum[ki] = cprod1[ki] + 2.0 * cprod2[ki] + cprod3[ki];
                  }    
	          gder[3] = s6scpr(sum,bvec,3);

	          s6crss(normuv,pqdiff,cprod1);
	          s6crss(normu,spv,cprod2);
	          s6crss(normv,spu,cprod3);
	          s6crss(norm,spuv,cprod4);
                  for (ki=0;ki<3;ki++)
                  {    
                      sum[ki] = cprod1[ki] + cprod2[ki] + cprod3[ki] + cprod4[ki];
                  }    
	          gder[4] = s6scpr(sum,bvec,3);

	          s6crss(normvv,pqdiff,cprod1);
	          s6crss(normv,spv,cprod2);
	          s6crss(norm,spvv,cprod3);
                  for (ki=0;ki<3;ki++)
                  {    
                      sum[ki] = cprod1[ki] + 2.0 * cprod2[ki] + cprod3[ki];
                  }    
	          gder[5] = s6scpr(sum,bvec,3);

                }
            }
        }
    }

  if(ideg == 1003 || ideg == 1004 || ideg == 1005)
    {
      /*  The normal vector has no meaning in this case
          so we calculate the direction D(u,v) of the solution curve
          on the surface P(u,v)
          of  f(u,v) = N(u,v) . d = 0  (ideg = 1003)
          or  f(u,v) = N(u,v) . (P(u,v) - E) = 0  (ideg = 1004)
          or  f(u,v) = N(u,v) x (P(u,v) - Q) . B = 0  (ideg = 1005)
          (through the present solution point)
          directly. This direction D is in fact  f P  - f P
                                                  v u    u v
          where P(u,v) is the parameterised
          surface, N is the normal P  x P  , and d is the 
                                    u    v
	  view direction. The direction D is used by s1331.  */
     
      for(ki=0; ki<3; ki++)
        {
          gnorm[ki] = gder[2] * spu[ki] - gder[1] * spv[ki];
        }
      (void)s6norm(gnorm,3,gnorm,&kstat);
    }

  /* Everyting is ok */
  
  *jstat = 0;
  goto out;
  
  
  /* Error in lower level function */
 error:
  *jstat = kstat;
  s6err("s1331",*jstat,kpos);
  goto out;
  
 out:
  return;
}
