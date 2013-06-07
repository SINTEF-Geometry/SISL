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
 * $Id: s1310.c,v 1.8 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1310

#include "sislP.h"
/*
* Forward declarations.
* ---------------------
*/
#if defined(SISLNEEDPROTOTYPES)
static  void s1310_s9constline(SISLSurf *,SISLSurf *,SISLIntcurve *,
                               double,int,int,int *);
#else
static  void s1310_s9constline();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
     s1310(SISLSurf *psurf1,SISLSurf *psurf2,SISLIntcurve *pinter,
	   double aepsge,double amax,int icur,int igraph,int *jstat)
#else
void s1310(psurf1,psurf2,pinter,aepsge,amax,icur,igraph,jstat)
     SISLSurf     *psurf1;
     SISLSurf     *psurf2;
     SISLIntcurve *pinter;
     double   aepsge;
     double   amax;
     int      icur;
     int      igraph;
     int      *jstat;
#endif
/*

*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To march an intersection curve between two B-splines
*              surfaces. The intersection curve is described by guide
*              parameter pairs in an intersection curve object.
*
*
*
* INPUT      : psurf1 - Pointer to first surface
*              psurf2 - Pointer to second surface
*              pinter - Pointer to intersection curve.
*                       The guide parameter pairs refered by the object are
*                       used for guiding the marching.
*              aepsge - Absolute tolerance
*              amax   - Not used.
*              icur   - Indicator telling if a 3-D curve is to be made
*                        0 - Don't make 3-D curve
*                        1 - Make 3-D curve
*                        2 - Make 3-D curve and curves in parameter plane
*              igraph - Indicator telling if the curve is to be outputted
*                       through function calls:
*                        0 - don't output curve through function call
*                        1 - output as straight line segments through s1line
*
*
*
* OUTPUT     : pinter - Pointer to intersection curve. The routine
*                       adds intersection curve and curve in the parameter
*                       planes to the SISLIntcurve object according to the
*                       values of i3Dcur and iplane
*                       If these curves have already been generated in the
*                       topology part of the intersections, they will first
*                       be free'ed.  This makes it possible to generate
*                       curves for both parameter planes if required.
*                       The geometry will have been generated in the case when
*                       the intersection curve represents a constant parameter
*                       line in the parmeter plane of the surface.

*              jstat  - status messages
*                         = 3      : Iteration stopped due to singular
*                                    point or degenerate surface. A part
*                                    of intersection curve may have been
*                                    traced out. If no curve is traced out
*                                    the curve pointers in the Intcurve
*                                    object point to SISL_NULL.
*                         = 0      : ok
*                         < 0      : error
*                         = -185   : No points produced on intersection curve.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 30 june 1988
* Revised by : Tor Dokken, SI, Oslo, Norway, 24-feb-1989
*              Handles degenerate points
* Revised by : Tor Dokken, SI, Oslo, Norway, 03-April-1988
*              Maximal step length calculation, new strategy around
*              singular points, error correction
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Dec. 1994.  Added check for
*              SISL_NULL 'pinter' and to avoid re-generating the geometry when it has
*              already been generated in the topology part (constant curve, type 9).
*              This fixes memory problems.
*
*********************************************************************
*/
{
  int ki,kj,kl;            /* Control variables in for loops            */
  int kcont;               /* Stop condition for loop                   */
  int kk,kn;               /* Dummy variables                           */
  int kstpch;              /* Status of iteration step                  */
  int kpoint;              /* Number of points in guide curve           */
  int kpar1;               /* Number of parameter direction in 1st. obj */
  int kpar2;               /* Number of parameter direction in 2nd. obj */
  int kpar;                /* Indicater tellin if s1359 shall make
			      parametrization or use parametrization
			      in spar                                   */
  int ktype;               /* Type of intersection curve                */
  int klfu=0;              /* Pointers into knot vectors                */
  int klfv=0;              /* Pointers into knot vectors                */
  int klfs=0;              /* Pointers into knot vectors                */
  int klft=0;              /* Pointers into knot vectors                */
  int kder = 2;            /* Calculate up to second derivatives        */
  int kdim = 3;            /* The dimension of the space we work in     */
  int kfirst = 0;          /* Indicator telling if first guide point
			      degenerate */
  int klast = 0;           /* Indicator telling if last guide point
			      degenerate */
  int kpos = 0;            /* Position of error                         */
  int kstat,kstat1;        /* Status variable returned form routine     */
  int kmaxinf=0;           /* Number of entries object that can be stored
			      in s3dinf, sp1inf, sp2inf                 */
  int knbinf=0;            /* Number of entries stored so far on s3dinf,
			      sp1inf and sp2inf                         */
  int kstart;              /* Start point for iteration among guide pnts*/
  int kguide;              /* Current guide point                       */
  int kdir;                /* March direction                           */
  int kgdir;               /* Direction we march guide point vector     */
  int krem,krem1,krem2;    /* REmember status of boundary crossing      */
  int kbound;              /* Whci boundary is crossed                  */
  int koutside_resolution; /* Flag telling if current seg. outside res. */
  int kdiv=0;              /* Flag telling if iteration diverged        */
  double tlnorm=DZERO;     /* Length of normal vector                   */
  double tltan1=DZERO;     /* Length of tangents                        */
  double tltan2=DZERO;     /* Length of tangents                        */
  double tang1,tang2;      /* Angles                                    */
  int knb1=0;              /* Remember number of points after marching
			      in first marching direction               */
  int kgd1=0;              /* Remeber last guide point used in first
			      marching direction                        */
  double *scorpnt=SISL_NULL;    /* Corrected marching points                 */
  double *scorpr1=SISL_NULL;    /* Corrected marching parameter values in ps1*/
  double *scorpr2=SISL_NULL;    /* Corrected marching parameter values in ps2*/
  double smidd[6];         /* Description of midpoint and tangent of
			      current Bezier segment                    */
  double tcurstep;         /* Current step length                       */
  double tdist;            /* Error at middle of current Bezier segement*/
  double tang;             /* Angle error at midpoint Bezier segement   */
  double tnew;             /* Candidate for new step length             */
  double tfak;             /* How much is the step length to be reduced */
  double *start;           /* Pointer to start of current segment       */
  double *st;              /* Pointer to knot vector                    */
  double tstep;            /* Iteration step length                     */
  double tmax;             /* Local maximal step length                 */
  double tstartstp;        /* Start step length                         */
  double trad;             /* Radius of curvature                       */
  double spar1[2];         /* Parameter pair of current point surface 1 */
  double spar2[2];         /* Parameter pair of current point surface 2 */
  double sparmid1[2];      /* Parameter values at middle of Bezier segm */
  double sparmid2[2];      /* Parameter values at middle of Bezier segm */
  double sipar1[2];        /* Parameter pair iteration point surface  1 */
  double sipar2[2];        /* Parameter pair iteration point surface 2  */
  double simiddpnt[10];    /* Middle point and tangent of segment       */
  double simiddpar1[7];    /* Parameter value at middle point of segment*/
  double simiddpar2[7];    /* Parameter value at middle point of segment*/
  double startg[3];        /* Tangent of start point of iteration       */
  double *sgpar1=SISL_NULL;     /* Parameter pairs of guide point in surf 1  */
  double *sgpar2=SISL_NULL;     /* Parameter pairs of guide point in surf 2  */
  double *sgpara=SISL_NULL;     /* Parameter pairs of guide point in surf 1  */
  double *sgparb=SISL_NULL;     /* Parameter pairs of guide point in surf 2  */
  double *sgd1 = SISL_NULL;     /* 0-2 derivative of guide point + normal
			      of first object                           */
  double *sgd2 = SISL_NULL;     /* 0-2 derivative of guide point + normal
			      of second object                          */
  double spnt1[21];        /* Info on current point in first surface    */
  double spnt2[21];        /* Info on current point in second surface   */
  double sipnt1[21];       /* Info on iteration point in first surface  */
  double sipnt2[21];       /* Info on iteration point in second surface */
                           /* For spnt1, spnt2, sipnt1, sipnt2,         */
                           /* the information is stored 3-tuppels       */
                           /* in the following sequence                 */
                           /* Position, (1,0)-der, (0,1)-der,           */
                           /* (2,0)-der, (1,1)-der (0,2)-der and normal */
                           /* This is compatible with output of s1421   */
  double spntend1[21];     /* End values of candidate end point in surf1*/
  double spntend2[21];     /* End values of candidate end point in surf1*/
  double sparend1[2];      /* Parameter value at candidate end point    */
  double sparend2[2];      /* Parameter value at candidate end point    */
  double *snxt1;           /* SISLPoint in psurf1 we have accepted          */
  double *snxt2;           /* SISLPoint in psurf2 we have accepted          */
  double *snxp1;           /* Parameter value belonging to snxt1        */
  double *snxp2;           /* Parameter value belonging to snxt2        */
  double *s3dinf=SISL_NULL;     /* Pointer to array used for storing 3-D position
			      tangent, curvature and radius of curvature found
			      during the marching process */
  double *sp1inf=SISL_NULL;     /* Pointer to array used for storing position
			      tangent, curvature and radius of curvature found
			      in the first parameter plane during the
			      marching process */
  double *sp2inf=SISL_NULL;     /* Pointer to array used for storing position
			      tangent, curvature and radius of curvature found
			      in the first parameter plane during the
			      marching process */
  double *spar=SISL_NULL;       /* Parametrization of points                 */
  double sval1[2];         /* Limits of parameter plane in first SISLdir    */
  double sval2[2];         /* Limits of parameter plane in second SISLdir   */
  double sval3[2];         /* Limits of parameter plane in third SISLdir    */
  double sval4[2];         /* Limits of parameter plane in fourth SISLdir   */
  double tref1,tref2;      /* Reference values for knot vectors         */
  double tref3,tref4;      /* Reference values for knot vectors         */
  double start1[21];       /* Description of start point in psurf1      */
  double start2[21];       /* Description of start point in psurf1      */
  double stpar1[2];        /* Parameter pair belonging to start1        */
  double stpar2[2];        /* Parameter pair belonging to start2        */
  double sdum1[3],sdum2[3];/* Dummy vectors                             */
  double tdum,tdump1,tdump2;/*Dummy variable                            */
  double *sp1=SISL_NULL;        /* Pointer used when moving information      */
  double *sp2=SISL_NULL;        /* Pointer used when moving information      */
  double stdum[10];        /* Dummy array used when moving information  */
  double *stang;           /* Pointer to tangent of current point       */
  double *stangp1;         /* Pointer to tangent of current point in pp1*/
  double *stangp2;         /* Pointer to tangent of current point in pp2*/
  double *spoint;          /* Pointer to current point                  */
  double t1distgd,t2distgd;/* Distances to guide points                 */
  SISLCurve *q3dcur=SISL_NULL;/* Pointer to 3-D curve                     */
  SISLCurve *qp1cur=SISL_NULL;/* Pointer to curve in first parameter plane*/
  SISLCurve *qp2cur=SISL_NULL;/* Pointer to curve in 2.nd  parameter plane*/


  *jstat = 0;

  if ( pinter == SISL_NULL )  goto err150;


  /* Check if the geometry already has been generated in the topology part.
     This will be the case if the geometry is along a constant parameter line.
     Freeing the geometry her makes it possible to generate curves for both
     parameter planes if required (the pointers will be set to SISL_NULL further
     down, i.e. would cause a memory leak if they weren't free'ed here. */

  if (pinter->itype == 9)
  {
    if (pinter->pgeom)  freeCurve(pinter->pgeom);
    if (pinter->ppar1)  freeCurve(pinter->ppar1);
    if (pinter->ppar2)  freeCurve(pinter->ppar2);
  }


  /* Make maximal step length based on box-size of surface */

  sh1992su(psurf1,0,aepsge,&kstat);
  if (kstat < 0) goto error;

  tmax = MAX(psurf1->pbox->e2max[0][0] - psurf1->pbox->e2min[0][0],
	     psurf1->pbox->e2max[0][1] - psurf1->pbox->e2min[0][1]);
  tmax = MAX(tmax,psurf1->pbox->e2max[0][2] - psurf1->pbox->e2min[0][2]);

  sh1992su(psurf2,0,aepsge,&kstat);
  if (kstat < 0) goto error;

  tmax = MAX(tmax,psurf2->pbox->e2max[0][0] - psurf2->pbox->e2min[0][0]);
  tmax = MAX(tmax,psurf2->pbox->e2max[0][1] - psurf2->pbox->e2min[0][1]);
  tmax = MAX(tmax,psurf2->pbox->e2max[0][2] - psurf2->pbox->e2min[0][2]);

  if (amax>DZERO) tmax = MIN(tmax,amax);

  /* Find a none singular start point for the marching process */

  kpoint = pinter->ipoint;
  kpar1  = pinter->ipar1;
  kpar2  = pinter->ipar2;
  sgpara = pinter->epar1;
  sgparb = pinter->epar2;
  ktype  = pinter->itype;


  /* To support closed curve the first guide point must be copied after
     the last guide point */

  if((sgpar1=newarray(2*kpoint+2,DOUBLE)) == SISL_NULL) goto err101;
  if((sgpar2=newarray(2*kpoint+2,DOUBLE)) == SISL_NULL) goto err101;
  memcopy(sgpar1,sgpara,2*kpoint,DOUBLE);
  memcopy(sgpar2,sgparb,2*kpoint,DOUBLE);

  if (ktype ==2 || ktype == 3)
    {
      /*Closed curve copy first guide point to end of string of guide points */
      memcopy(sgpar1+2*kpoint,sgpara,2,DOUBLE);
      memcopy(sgpar2+2*kpoint,sgparb,2,DOUBLE);
      kpoint = kpoint + 1;
    }

  /* Initiate pointers to intersection curve and intersection curve in
     parameter plane */

  pinter -> pgeom = SISL_NULL;
  pinter -> ppar1 = SISL_NULL;
  pinter -> ppar2 = SISL_NULL;

  /* Initiate parameter direction boundaries */

  kk    = psurf1 -> ik1;
  kn    = psurf1 -> in1;
  st    = psurf1 -> et1;
  sval1[0] = st[kk-1];
  sval1[1] = st[kn];
  tref1 = (double)3.0*MAX(fabs(*sval1),fabs(*(sval1+1)));
  kk    = psurf1 -> ik2;
  kn    = psurf1 -> in2;
  st    = psurf1 -> et2;
  sval2[0] = st[kk-1];
  sval2[1] = st[kn];
  tref2 = (double)3.0*MAX(fabs(*sval2),fabs(*(sval2+1)));
  kk    = psurf2 -> ik1;
  kn    = psurf2 -> in1;
  st    = psurf2 -> et1;
  sval3[0] = st[kk-1];
  sval3[1] = st[kn];
  tref3 = (double)3.0*MAX(fabs(*sval3),fabs(*(sval3+1)));
  kk    = psurf2 -> ik2;
  kn    = psurf2 -> in2;
  st    = psurf2 -> et2;
  sval4[0] = st[kk-1];
  sval4[1] = st[kn];
  tref4 = (double)3.0*MAX(fabs(*sval4),fabs(*(sval4+1)));



  /* Test the both objects have 2 parameter directions */

  if (kpar1 != 2 || kpar2 != 2) goto err123;

  /*THE POINTS , TANGENT, CURVATURE AND RADIUS OF CURVATURE FOUND DURING
    THE MARCHING PROCESS SHOULD ALL BE STORED IN ARRAYS. ALLOCATE ONE ARRAY
    FOR 3-D INFORMATION , ONE ARRAY FOR INFORMATION IN FIRST PARAMETER PLANE
    AND ONE ARRAY FOR INFORMATION IN SECOND PARAMETER PLANE. THESE ARRAYS
    ARE GIVEN AN INITIAL CAPACITY OF STORING 100 POINTS WITH OTHER INFORMATION.
    IF THEY ARE TO SHORT THEY WILL BE REALLOCATED AT A LATER STAGE.

    SINCE THE STEPPING WILL GO IN BOTH DIRECTIONS WE WILL HAVE TO TURN THE
    INFORMATION FOUND WHEN MARCHING IN NEGATIVE DIRECTION, SO THAT IT CAN
    BE COMBINED WITH THE INFORMATION FOUND WHEN WE ARE MARCHING IN POSITVE
    DIRECTION.
    */

  kmaxinf = 100;
  s3dinf = newarray(10*kmaxinf,DOUBLE);
  if (s3dinf == SISL_NULL) goto err101;
  sp1inf = newarray(7*kmaxinf,DOUBLE);
  if (sp1inf == SISL_NULL) goto err101;
  sp2inf = newarray(7*kmaxinf,DOUBLE);
  if (sp2inf == SISL_NULL) goto err101;



  /* Evaluate 0-1-2nd. derivative + normal of all guide points in both
     surfaces, first allocate arrays for storing the information */

  sgd1 = newarray(21*kpoint,DOUBLE);
  if (sgd1==SISL_NULL) goto err101;
  sgd2 = newarray(21*kpoint,DOUBLE);
  if (sgd2==SISL_NULL) goto err101;

  kpos = 5;

  /* Initiate kstart to point at no point */

  kstart = 0;

  for (ki=0,kj=0,kl=0 ; ki<kpoint ; ki++,kj+=2,kl+=21)
    {
      s1421(psurf1,kder,&sgpar1[kj],&klfu,&klfv,&sgd1[kl],&sgd1[kl+18],&kstat);
      if (kstat<0) goto error;

      /*  Find length of normal vector and tangent vectors */

      tlnorm = s6length(&sgd1[kl+18],kdim,&kstat);
      tltan1 = s6length(&sgd1[kl+ 3],kdim,&kstat);
      tltan2 = s6length(&sgd1[kl+ 6],kdim,&kstat);

      /* The cross product satisifes the following conditions:
	 length(axb) = length(a) length(b) sin(angle(a,b)).
	 Thus the angle between the two vectors can be found, close to 0
	 sin(a) is a good approximation of a */

      if (tlnorm == DZERO || tltan1 ==DZERO || tltan2 == DZERO)
        tang1 = DZERO;
      else
        tang1 = tlnorm/(tltan1*tltan2);

      s1421(psurf2,kder,&sgpar2[kj],&klfs,&klft,&sgd2[kl],&sgd2[kl+18],&kstat);
      if (kstat<0) goto error;

      /*  Find length of normal vector and tangent vectors */

      tlnorm = s6length(&sgd2[kl+18],kdim,&kstat);
      tltan1 = s6length(&sgd2[kl+ 3],kdim,&kstat);
      tltan2 = s6length(&sgd2[kl+ 6],kdim,&kstat);

      /*  The cross product satisifes the follwing conditions:
	  length(axb) = length(a) length(b) sin(angle(a,b)).
	  Thus the angle between the two vectors can be found, close to 0
	  sin(a) is a good approximation of a */

      if (tlnorm == DZERO || tltan1 ==DZERO || tltan2 == DZERO)
        tang2 = DZERO;
      else
        tang2 = tlnorm/(tltan1*tltan2);


      if (tang1 >= ANGULAR_TOLERANCE && tang2 >= ANGULAR_TOLERANCE)
        {
	  /* Make tangent of intersection curve */

	  s6crss(&sgd1[kl+18],&sgd2[kl+18],sdum1);

	  tlnorm = s6length(sdum1,kdim,&kstat);

	  /* Remember if start, internal or end point */

	  if (tlnorm != DZERO)
	    {
	      if (ki == 0)
		kfirst = 1;
	      else if (ki == kpoint-1)
		klast = kpoint;
	      else
		kstart = ki+1;
	    }
        }
    }


  /* Check if only degenerate points or singularities exist on the
     intersection curve */

  if (kstart == 0)
    {
      /*  No internal nondegenerate point exits, start marching from first
	  or last point if possible */

      if (kfirst != 0 && ktype != 5 && ktype != 7) kstart = kfirst;
      else if (klast != 0 && ktype != 6 &&
	       ktype != 7 && ktype != 3) kstart = klast;
      else if (kfirst != 0) kstart = kfirst;
      else if (klast != 0) kstart = klast;
      else goto interpolate;
    }


  /* To speed up the marching process when many guide points are given,
     remove guide points that are not at the start, end or the start point */

  if (kpoint >2 && (kstart==1 || kstart==kpoint) )
    {
      /*  No internal guide point necessary, copy last point to second point */
      memcopy(sgd1+21,sgd1+21*(kpoint-1),21,DOUBLE);
      memcopy(sgpar1+2,sgpar1+2*(kpoint-1),2,DOUBLE);
      memcopy(sgd2+21,sgd2+21*(kpoint-1),21,DOUBLE);
      memcopy(sgpar2+2,sgpar2+2*(kpoint-1),2,DOUBLE);

      if (kstart ==  kpoint) kstart = 2;
      kpoint = 2;
    }
  else if (kpoint>2)
    {
      /*  Internal guide point exists, copy this to second position and
	  copy end point to third position */

      memcopy(sgd1+21,sgd1+21*(kstart-1),21,DOUBLE);
      memcopy(sgpar1+2,sgpar1+2*(kstart-1),2,DOUBLE);
      memcopy(sgd2+21,sgd2+21*(kstart-1),21,DOUBLE);
      memcopy(sgpar2+2,sgpar2+2*(kstart-1),2,DOUBLE);

      memcopy(sgd1+2*21,sgd1+21*(kpoint-1),21,DOUBLE);
      memcopy(sgpar1+4,sgpar1+2*(kpoint-1),2,DOUBLE);
      memcopy(sgd2+2*21,sgd2+21*(kpoint-1),21,DOUBLE);
      memcopy(sgpar2+4,sgpar2+2*(kpoint-1),2,DOUBLE);

      kpoint = 3;
      kstart = 2;
    }

  /* Remember description of start point in both surfaces,
     copy point indicated by kstart into spnt1,spnt2,spar1,spar2 */

  memcopy(spnt1,sgd1+21*(kstart-1),21,DOUBLE);
  memcopy(spnt2,sgd2+21*(kstart-1),21,DOUBLE);
  memcopy(spar1,sgpar1+2*(kstart-1),2,DOUBLE);
  memcopy(spar2,sgpar2+2*(kstart-1),2,DOUBLE);

  /* Make position, unit tangent, curvature and radius of curvature for
     start point of iteration, store them in the arrays just allocated */

  kpos = 10;
  s1304(spnt1,spnt2,spar1,spar2,s3dinf,sp1inf,sp2inf,&kstat);

  if (kstat<0) goto error;

  /* Remember start tangent */

  memcopy(startg,s3dinf+3,3,DOUBLE);


  /* Iterate intersection point down to the intersection curve */

  tstep = DZERO;
  s9iterate(s3dinf,spnt1,spnt2,spar1,spar2,psurf1,psurf2,tstep,
	    aepsge,sipnt1,sipnt2,sipar1,sipar2,&kstat);
  if (kstat < 0) goto error;

  /* Copy result of iteration into spnt1,spnt2,spar1,spar2 */


  if (kstat==0 &&
      (s6dist(spnt1,sipnt1,3) > aepsge || s6dist(spnt2,sipnt2,3) > aepsge))
    {
      /*  Copy result of iteration of convergence to no singular point */

      memcopy(spnt1,sipnt1,21,DOUBLE);
      memcopy(spnt2,sipnt2,21,DOUBLE);
      memcopy(spar1,sipar1,2,DOUBLE);
      memcopy(spar2,sipar2,2,DOUBLE);
    }

  if (kstat==0)
    {
      memcopy(start1,sipnt1,21,DOUBLE);
      memcopy(start2,sipnt2,21,DOUBLE);
      memcopy(stpar1,sipar1,2,DOUBLE);
      memcopy(stpar2,sipar2,2,DOUBLE);
    }
  else
    {
      memcopy(start1,spnt1,21,DOUBLE);
      memcopy(start2,spnt2,21,DOUBLE);
      memcopy(stpar1,spar1,2,DOUBLE);
      memcopy(stpar2,spar2,2,DOUBLE);
    }

  /* Make position, unit tangent, curvature and radius of curvature for
     start point of iteration, store them in the arrays just allocated */

  kpos = 10;
  s1304(start1,start2,stpar1,stpar2,s3dinf,sp1inf,sp2inf,&kstat);

  if (kstat<0) goto error;

  /* Test if singular point reached */

  if (kstat == 2) goto war03;

  /* Remember that start point is already stored */

  knbinf = 1;

  /* Make step length based on 3-D radius of curvature, tolerances and
     maks step length */

  kpos = 20;
  tstep = s1311(s3dinf[9],aepsge,tmax,&kstat);
  if (kstat<0) goto error;
  tstartstp = tstep;

  /* STEP IN BOTH DIRECTIONS FROM THE FOUND START POINT */

  /* Indicate that direction in guide point array not determined */

  kguide = kstart;
  kgdir = 0;

  for (kdir=1;kdir<3;kdir++)
    {

      if (kdir == 2)
        {
	  /* Remember result of marching in first direction */

	  knb1 = knbinf;
	  kgd1 = kguide;

	  /* If the previous step direction made no points then knbinf==0. To
	     enable the marching we start from the same start point as the
	     previous step direction, thus in this case knbinf should be 1. */

	  knbinf = MAX(1,knbinf);

	  /* We now step in the second step direction. Turn the sequence of
	     the points found as well as change tangent directions */

	  /* First interchange 3-D info */

	  for (sp1=s3dinf,sp2=s3dinf+10*(knbinf-1) ; sp1<sp2 ; sp1+=10,sp2-=10)
            {
	      memcopy(stdum,sp1,  10,DOUBLE);
	      memcopy(sp1  ,sp2,  10,DOUBLE);
	      memcopy(sp2  ,stdum,10,DOUBLE);
            }

	  for (sp1=s3dinf+3;sp1<s3dinf+10*knbinf;sp1+=10)
            {
	      sp1[0] = - sp1[0];
	      sp1[1] = - sp1[1];
	      sp1[2] = - sp1[2];
            }

	  /* Then interchange info in first parameter plane */

	  for (sp1=sp1inf,sp2=sp1inf+7*(knbinf-1) ; sp1<sp2 ; sp1+=7,sp2-=7)
            {
	      memcopy(stdum,sp1  ,7,DOUBLE);
	      memcopy(sp1  ,sp2  ,7,DOUBLE);
	      memcopy(sp2  ,stdum,7,DOUBLE);
            }

	  for (sp1=sp1inf+2;sp1<sp1inf+7*knbinf;sp1+=7)
            {
	      sp1[0] = - sp1[0];
	      sp1[1] = - sp1[1];
	      sp1[2] = - sp1[2];
            }

	  /* Then interchange info in second parameter plane */

	  for (sp1=sp2inf,sp2=sp2inf+7*(knbinf-1) ; sp1<sp2 ; sp1+=7,sp2-=7)
            {
	      memcopy(stdum,sp1  ,7,DOUBLE);
	      memcopy(sp1  ,sp2  ,7,DOUBLE);
	      memcopy(sp2  ,stdum,7,DOUBLE);
            }

	  for (sp1=sp2inf+2;sp1<sp2inf+7*knbinf;sp1+=7)
            {
	      sp1[0] = - sp1[0];
	      sp1[1] = - sp1[1];
	      sp1[2] = - sp1[2];
            }

	  /* Turn direction of remembered start tangent */

	  for (ki=0;ki<3;ki++)
	    startg[ki] = -startg[ki];


	  /* Update spnt1, spnt2, spar1 and spar2 to
	     have the start point values */

	  memcopy(spnt1,start1,21,DOUBLE);
	  memcopy(spnt2,start2,21,DOUBLE);
	  memcopy(spar1,stpar1,2,DOUBLE);
	  memcopy(spar2,stpar2,2,DOUBLE);

	  /* Turn the direction we march the guide point vector,
	     and set current guide point to kstart */

	  kgdir  = -kgdir;
	  kguide = kstart;

	  /* Update step length */

	  tstep = tstartstp;
        }

      kpos = 30;

      /* Step direction ok, perform marching until stop condition reached */

      kcont = 1;

      while (kcont)
        {


	  /* We must make sure that we are not stepping past a guide point.
	   * Thus if we get close to a guide point, make sure that we step
	   * through this. The direction we travers the guide point array
	   * might not have been determined yet. Thus we have to test in
	   * both directions in guide point array.
	   *
	   * Remember how we step in the varaible kstpch:
	   *       kstpch = -1 : Try to step to previous guide point
	   *       kstpch =  0 : Try not to step through guide point
	   *       kstpch =  1 : Try to step to next guide point
	   *       kstpch =  3 : Step to start point and stop marching
	   *       kstpch =  4 : Don't step through guide point, candidate
	   *                     end point of segement found in iteration loop
	   */


	  kstpch = 0;
	  stang = s3dinf + 10*(knbinf-1) + 3;
	  stangp1 = sp1inf + 7*(knbinf-1) + 2;
	  stangp2 = sp2inf + 7*(knbinf-1) + 2;

	  if (kgdir >=0)
            {

	      /* We are stepping in positive direction in guide point vector
	       * calculate distance to next guide point. If the guide point
	       * is lying closer than the step length to the current point
	       * we should step directly to this point provided that the cross
	       * product of the normal vectors at current point and at the
	       * guide point have the same direction, e.g. that their scalar
	       * product is positiv
	       */

	      t1distgd = (double)2.0*tstep;
	      if (kguide < kpoint)
                {
		  /* Decide if we should step through the guide point */

		  kpos = 40;
		  t1distgd = s9adstep(spnt1,spar1,spnt2,spar2,&sgd1[kguide*21],
				      &sgpar1[kguide*2],&sgd2[kguide*21],
				      &sgpar2[kguide*2],stang,
				      stangp1,stangp2,tstep,&kstat);
		  if (kstat<0) goto error;
		  if (kstat == 1)
                    {
		      /* Step through guide point remember this */

		      kstpch = 1;
		      snxt1 = sgd1 + 21*kguide;
		      snxt2 = sgd2 + 21*kguide;
		      snxp1 = sgpar1 + 2*kguide;
		      snxp2 = sgpar2 + 2*kguide;
                      tstep = MIN(tstep,t1distgd);
		    }
                }
            }

	  if (kgdir <=0)
            {

	      /* We are stepping in negative direction in guide point vector
	       * calculate distance to previous guide point. If the guide point
	       * is lying closer than the step length to the current point
	       * we should step directly to this point provided that the cross
	       * product of the normal vectors at current point and at the
	       * guide point have the same direction, e.g. that their scalar
	       * product is positiv
	       */

	      if (1 < kguide)
                {
		  /* Decide if we should step through the guide point */

		  kpos = 50;
		  t2distgd = s9adstep(spnt1,spar1,spnt2,spar2,
				    &sgd1[(kguide-2)*21],&sgpar1[(kguide-2)*2],
				    &sgd2[(kguide-2)*21],&sgpar2[(kguide-2)*2],
				    stang,stangp1,stangp2,tstep,&kstat);

		  if (kstat<0) goto error;
		  if ((kstat == 1 &&kstpch == 0) ||
		      (kstat == 1 && kstpch == 1 && t2distgd < t1distgd))
                    {
		      /* Step through guide point remember this */

		      kstpch = -1;
		      snxt1 = sgd1 + 21*(kguide-2);
		      snxt2 = sgd2 + 21*(kguide-2);
		      snxp1 = sgpar1 + 2*(kguide-2);
		      snxp2 = sgpar2 + 2*(kguide-2);
                      tstep = MIN(tstep,t2distgd);
                    }
                }
            }

	  /* Check if we step through the start point, should only be necessary
	     if at least 3 points found in this marching direction */

	  if ((kdir==1 && knbinf>3) || (kdir==2 && knbinf>knb1+2))
            {
	      kpos = 60;
	      tdum = s9adstep(spnt1,spar1,spnt2,spar2,start1,stpar1,start2,
			      stpar2,stang,stangp1,stangp2,tstep,&kstat);
	      if (kstat<0) goto error;
	      if (kstat == 1)
		{
		  /* Step to start point remember this */

		  kstpch = 3;
		  snxt1 = start1;
		  snxt2 = start2;
		  snxp1 = stpar1;
		  snxp2 = stpar2;
                  tstep = MIN(tstep,tdum);

		}
	    }

	  /* At this stage kstpch=0 if we have not reached a guide point or
	     if we have not reached the start point of the iteration.

	     Now we want to find a Bezier segement that is lying within the
	     geometric tolerance that is approximating the intersection curve.
	     If a guide point is reached (kstpch=-1 or 1), then we have a
	     candidate for the end point of the Bezier segement. If the start
	     point is reached (kstpch=3) then we also have a candidate
	     end point for the segment.

	     The next loop use kstpch to indicate if we have a candidate
	     end point for the segment:

	     kstpch==0  :  No candidate end point exists
	     kstpch!=0  :  Candidate end point exists


	     To indicate if the segement is within the resolution we
	     use koutside_resolution:

	     koutside_resolution==0 : Segment outside resolution
	     koutside_resolution!=0 : Segment inside resolution
	     */

	  koutside_resolution = 0;


	  /* Make sure that there is enough space for one more point */

	  if (knbinf>=kmaxinf)
            {
	      kmaxinf = kmaxinf + 100;
	      s3dinf = increasearray(s3dinf,((3*kdim+1)*kmaxinf),DOUBLE);
	      if (s3dinf==SISL_NULL) goto err101;
	      sp1inf = increasearray(sp1inf,7*kmaxinf,DOUBLE);
	      if (sp1inf==SISL_NULL) goto err101;
	      sp2inf = increasearray(sp2inf,7*kmaxinf,DOUBLE);
	      if (sp2inf==SISL_NULL) goto err101;
            }


	  /*Make description of candidate endpoint if it exists and store it */

	  if (kstpch != 0)
            {
	      s1304(snxt1,snxt2,snxp1,snxp2,s3dinf+10*knbinf,
		    sp1inf+7*knbinf,sp2inf+7*knbinf,&kstat);
	      if (kstat<0) goto error;

	      /* It is allowed to jump on to a singular point
		 Make sure that the tangents of previous and the new point
		 point in the same direction */

	      if (knbinf>0)
		tdum = s6scpr(s3dinf+10*(knbinf-1)+3,
			      s3dinf+10*knbinf+3,kdim);
	      else
		tdum = s6scpr(startg,s3dinf+3,kdim);


	      if (tdum < DZERO)
                {
		  /* Change tangent direction 3-D and in parameter plane */
		  sp1 = s3dinf + 10*knbinf + 3;
		  sp1[0] = -sp1[0];
		  sp1[1] = -sp1[1];
		  sp1[2] = -sp1[2];
		  sp1 = sp1inf + 7*knbinf + 2;
		  sp1[0] = -sp1[0];
		  sp1[1] = -sp1[1];
		  sp1 = sp2inf + 7*knbinf + 2;
		  sp1[0] = -sp1[0];
		  sp1[1] = -sp1[1];
                }

	      /* Copy the candidate point to spntend1, spntend2,sparend1
		 and sparend2 */

	      memcopy(spntend1,snxt1,21,DOUBLE);
	      memcopy(sparend1,snxp1,2,DOUBLE);
	      memcopy(spntend2,snxt2,21,DOUBLE);
	      memcopy(sparend2,snxp2,2,DOUBLE);
            }


	  while (kstpch == 0 || koutside_resolution == 0)
            {
	      if (kstpch!=0)
                {
		  /* Candidate end point exist, iterate to find point close
		     to the midpoint of the Bezier segement */


		  /* Decide if Hermit shape acceptable and find position and
		     tangent at midpoint of segment */

		  start = s3dinf + 10*(knbinf-1);

		  s1361(start,start+10,3,smidd,smidd+3,&kstat);
		  if (kstat<0) goto error;

		  tcurstep = DZERO;
		  spoint = smidd;
                }
	      else
                {

		  /* Iterate to find end point of segment */

		  /* ITERATE by intersecting the two surface and the plane
		     defined by current point (s3dinf), the tangent (s3dinf+3)
		     and the step length */

		  spoint = s3dinf + 10*(knbinf-1);
		  tcurstep = tstep;
                }

	      /* Perform the actual iteration */

	      kpos = 70;
	      s9iterate(spoint,spnt1,spnt2,spar1,spar2,psurf1,psurf2,tcurstep,
			aepsge,sipnt1,sipnt2,sipar1,sipar2,&kstat);
	      if (kstat < 0) goto error;

	      /* Initiate distance between midpoint and iteration point
		 to -1 to enable detection of divergence */

	      tdist = (double)-1.0;

	      /* Check if iteration has converged */

	      if (kstat == 2)
                {
		  /* Iteration has diverged, half step length if possible,
		     find new endpoint of segement. */

		  kstpch = 0;
		  koutside_resolution = 0;
                }
	      else if(kstat == 1 && kstpch != 0)
                {
		  /* The point found is closer to the input point than
		     the relative computer resolution or is a singular point.
		     We stop the marching in this direction here
		     Half step length if possible, find new endpoint of
		     segement. */

		  kstpch = 0;
		  koutside_resolution = 0;
                }
	      else if (kstpch!=0)
                {


		  /* Make description of intersection point */

		  s1304(sipnt1,sipnt2,sipar1,sipar2,simiddpnt,simiddpar1,
			simiddpar2,&kstat);
		  if (kstat<0) goto error;


		  if (kstat != 2)
                    {
		      /* We iterated to find midpoint of segment, test if
			 it is within resolution */

		      tdist = s6dist(simiddpnt,smidd,3);
		      tang  = s6ang(simiddpnt+3,smidd+3,3);
                    }

		  /* If point is singular or not within resolution a new
		     Hermit segment has to be made */

		  if (kstat == 2 || (fabs(tdist) > aepsge ||
				     (fabs(tang) > ANGULAR_TOLERANCE &&
				      tstep      > aepsge)))
                    {
		      kstpch = 0;
		      koutside_resolution = 0;
                    }
		  else
                    {
		     /*Segment within tolerance.
		       Check that the relationship between the two surfaces
		       has not been interchanged, by making the cross product
		       of the normal vectors in current point and the point
		       found by iteration. Then make the scalar product of
		       these vectors. If the scalar product is negative then
		       we have either jumped to another branch or passed a
		       singularity,iterprete this as the iteration has diverged
		       In addition we don't want the direction of the tangents
		       change to much. We set a limit of approximately PI/3
		       Make normal vectors in implicit surface for both points
		       Make also sure that the curve in the parameter plane
		       does not turn more than 90 degrees.
		       by testing on a cosin value of 0.5
		      */

		      s6crss(spnt1+18,spnt2+18,sdum1);
		      (void)s6norm(sdum1,kdim,sdum1,&kstat);
		      if (kstat < 0) goto error;

		      s6crss(spntend1+18,spntend2+18,sdum2);
		      (void)s6norm(sdum2,kdim,sdum2,&kstat);
		      if (kstat < 0) goto error;

		      tdum = s6scpr(sdum1,sdum2,kdim);

                      s6diff(sipar1,spar1,2,sdum1);
                      tdump1 = s6scpr(sdum1,sp1inf+7*(knbinf-1)+2,2);

                      s6diff(sipar2,spar2,2,sdum1);
                      tdump2 = s6scpr(sdum1,sp2inf+7*(knbinf-1)+2,2);

		      if (tdum == DZERO)
                        {
			  double tl1,tl2;

			  /* If one of the tangents have zero length,
			     accept segment */

			  tl1 = s6length(sdum1,kdim,&kstat);
			  tl2 = s6length(sdum2,kdim,&kstat);

			  if (tl1 == DZERO || tl2 == DZERO)
			    koutside_resolution = 1;
			  else
                            {
			      /* Find new end point of segment */

			      koutside_resolution = 0;
			      kstpch = 0;
                            }

                        }
		      else if (tdum <= (double)0.5 || tdump1 <= DZERO
                               || tdump2 <= DZERO)
                        {
			  /*Find new end point of segment */

			  koutside_resolution = 0;
			  kstpch = 0;
                        }
		      else
                        {
			  koutside_resolution = 1;
                        }
                    }
                }
	      else
                {
		  /* We iterated to find end point of segment,
		     update pointer */

		  memcopy(spntend1,sipnt1,21,DOUBLE);
		  memcopy(sparend1,sipar1,2,DOUBLE);
		  memcopy(spntend2,sipnt2,21,DOUBLE);
		  memcopy(sparend2,sipar2,2,DOUBLE);

		  s1304(sipnt1,sipnt2,sipar1,sipar2,s3dinf+10*knbinf,
			sp1inf+7*knbinf,sp2inf+7*knbinf,&kstat);
		  if (kstat<0) goto error;

		  /* Make sure that the tangents of previous and the new point
		     point in the same direction, singular end point allowed' */

		  if (knbinf>0)
		    tdum = s6scpr(s3dinf+10*(knbinf-1)+3,
				  s3dinf+10*knbinf+3,kdim);
		  else
		    tdum = s6scpr(startg,s3dinf+3,kdim);


		  if (tdum < DZERO)
                    {
		      /* Change tangent direction 3-D and in parameter plane */

		      sp1 = s3dinf + 10*knbinf + 3;
		      sp1[0] = -sp1[0];
		      sp1[1] = -sp1[1];
		      sp1[2] = -sp1[2];
		      sp1 = sp1inf + 7*knbinf + 2;
		      sp1[0] = -sp1[0];
		      sp1[1] = -sp1[1];
		      sp1 = sp2inf + 7*knbinf + 2;
		      sp1[0] = -sp1[0];
		      sp1[1] = -sp1[1];
                    }
		  /* Indicate that end point accepted */

		  kstpch = 4;
		  koutside_resolution = 0;
                }
	      /* It the segment is acceptable clip to the boundary */

	      if (kstpch != 0 && koutside_resolution == 1)
		{

		  /* Check if the curve between the start and end point
		     cross the boundary */

		  memcopy(sparmid1,sipar1,2,double);
		  memcopy(sparmid2,sipar2,2,double);

		  s1330(spar1,spar2,sparend1,sparend2,sval1,sval2,sval3,sval4,
			&kbound,sipar1,sipar2,&kstat);
		  if (kstat<0) goto error;


		  /* In case of kstat==4 (we go from the boundary and out)
		     or kstat==0 and the start is within computer resolution
		     from the boundary, make sure that the tangent points out
		     in both parameter planes.
		     If not set status to 1 e.g, we are inside the patch */

		  if(kstat==0 || kstat==4)
		    {
		      /* Set pointer to tangents at start point */
		      ki = 7*(knbinf-1)+2;

		      if(((DEQUAL(spar1[1]+tref2,sval2[0]+tref2) &&
			   sp1inf[ki+1]>DZERO) ||
			   (DEQUAL(spar1[1]+tref2,sval2[1]+tref2) &&
			    sp1inf[ki+1]<DZERO) ||
			   (DEQUAL(spar1[0]+tref1,sval1[0]+tref1) &&
			    sp1inf[ki  ]>DZERO) ||
			   (DEQUAL(spar1[0]+tref1,sval1[1]+tref1) &&
			    sp1inf[ki  ]<DZERO)
			   ) &&
			 ((DEQUAL(spar2[1]+tref4,sval3[0]+tref4) &&
			   sp2inf[ki+1]>DZERO) ||
			  (DEQUAL(spar2[1]+tref4,sval3[1]+tref4) &&
			   sp2inf[ki+1]<DZERO) ||
			  (DEQUAL(spar2[0]+tref3,sval2[0]+tref3) &&
			   sp2inf[ki  ]>DZERO) ||
			  (DEQUAL(spar2[0]+tref3,sval2[1]+tref3) &&
			   sp2inf[ki  ]<DZERO)))
			kstat = 1;
		    }
		  krem1 = kstat;

		  /* Check if the curve between the start and midpoint cross
		     the boundary */

		  s1330(spar1,spar2,sparmid1,sparmid2,sval1,sval2,sval3,sval4,
			&kbound,sipar1,sipar2,&kstat);
		  if (kstat<0) goto error;

		  krem2 = kstat;

		  /* We now have the following cases:
		     kstat == 0 :
		     Line between (spar1,spar2) and (sparend1,sparend2)
		     outside. If this happens when kdir=1, then
		     just forget the start point. If it happens
		     when kdir=2, then we just stop the marching.
		     kstat == 1 : Line between epar1 and epar2 inside.
		     Continue iteration.
		     kstat == 2 : We step out of the patch. Clip to the edge
		     of the patch. Update start point.
		     kstat == 3 : We step into the patch. Clip to the edge
		     of the patch. Update endpoint
		     kstat == 4 : We go from the boundary and out. Try next
		     iteration direction.
		     */

		  if (krem1 == 0 || krem2 == 0)
		    {
		      if (kdir==1) knbinf--;
		      goto nextdir;
		    }
		  else if ((krem1 !=1 || krem2 !=1) &&
			   krem1 != 4 && krem2 != 4)
		    {

		      /* If we clip to the boundary,
			 forget any guide point identified */

		      kstat1 = 0;
		      if (krem2 == 2 || krem2 == 3)
			{
			  s9clipit(spar1,spar2,sparmid1,sparmid2,psurf1,psurf2,
				   sval1,sval2,sval3,sval4,aepsge,
				   sipnt1,sipnt2,sipar1,sipar2,&kstat);
			  if (kstat<0) goto error;
			  if (krem2==3 && kstat==1) kstpch = 4;
			  kstat1 = kstat;
			  krem = krem2;
			}
		      if (kstat1 !=1 && (krem1 == 2 || krem1 == 3))
			{
			  s9clipit(spar1,spar2,sparend1,sparend2,psurf1,psurf2,
				   sval1,sval2,sval3,sval4,aepsge,
				   sipnt1,sipnt2,sipar1,sipar2,&kstat);
			  if (kstat<0) goto error;
			  if (krem1==3 && kstat==1) kstpch = 4;
			  kstat1 = kstat;
			  krem = krem1;
			}

		      if (kstat1 == 1)
			{
			  /*Check that the relationship between the two
			    surfaces has not been interchanged,
			    by making the cross product
			    of the normal vectors in current point and
			    the point found by iteration.
			    Then make the scalar product of these vectors.
			    If the scalar product is negative then
			    we have either jumped to another branch or passed a
			    singularity, iterprete this as the iteration
			    has diverged.
			    In addition we don't want the direction of the
			    tangents change to much. We set a limit of
			    approximately PI/3 by testing on a
			    cosin value of 0.5
			    Make normal vectors in implicit surface for both
			    points Make also sure that the curves in the
			    parameter plane does not turn more than 90 degrees.
			    */

			  s6crss(spnt1+18,spnt2+18,sdum1);
			  (void)s6norm(sdum1,kdim,sdum1,&kstat);
			  if (kstat < 0) goto error;

			  s6crss(sipnt1+18,sipnt2+18,sdum2);
			  (void)s6norm(sdum2,kdim,sdum2,&kstat);
			  if (kstat < 0) goto error;

			  tdum = s6scpr(sdum1,sdum2,kdim);

			  /*Check that sipar1 lies on the same side of spar1 as
			    the tangent at spar1 */

			  s6diff(sipar1,spar1,2,sdum1);
			  tdump1 = s6scpr(sdum1,sp1inf+7*(knbinf-1)+2,2);

			  s6diff(sipar2,spar2,2,sdum1);
			  tdump2 = s6scpr(sdum1,sp2inf+7*(knbinf-1)+2,2);
			}

		      /* An intersection point has only been
			 found when kstat==1 */

		      if ( kstat1==1 && tdump1 >= DZERO &&
			  tdump1 >= DZERO && tdum > (double)0.5)

			{
			  /* If krem=3 we step into the patch,
			     if krem=2 we step
			     out of the patch */

			  if (krem==2 || krem==3)
			    {
			      /* If krem==3 we step into the patch, make new
				 start point of segment */

			      if (krem==3) knbinf--;

			      memcopy(spntend1,sipnt1,21,DOUBLE);
			      memcopy(sparend1,sipar1,2,DOUBLE);
			      memcopy(spntend2,sipnt2,21,DOUBLE);
			      memcopy(sparend2,sipar2,2,DOUBLE);

			      s1304(sipnt1,sipnt2,sipar1,sipar2,
				    s3dinf+10*knbinf,
				    sp1inf+7*knbinf,
				    sp2inf+7*knbinf,&kstat);
			      if (kstat<0) goto error;

			      /* Make sure that the tangents of previous
				 and the new point
				 point in the same direction */

			      if (knbinf>0)
				tdum = s6scpr(s3dinf+10*(knbinf-1)+3,
					      s3dinf+10*knbinf+3,kdim);
			      else
				tdum = s6scpr(startg,s3dinf+3,kdim);


			      if (tdum < DZERO)
				{
				  /* Change tangent direction 3-D and in
				     parameter plane */

				  sp1 = s3dinf + 10*knbinf + 3;
				  sp1[0] = -sp1[0];
				  sp1[1] = -sp1[1];
				  sp1[2] = -sp1[2];
				  sp1 = sp1inf + 7*knbinf + 2;
				  sp1[0] = -sp1[0];
				  sp1[1] = -sp1[1];
				  sp1 = sp2inf + 7*knbinf + 2;
				  sp1[0] = -sp1[0];
				  sp1[1] = -sp1[1];
				}
			      /* If the new end point tangent points out go to
				 next direction */

			      ki = 7*knbinf;
			      if((sp1inf[ki+1] <= sval2[0] &&
				  sp1inf[ki+3] < DZERO) ||
				 (sp1inf[ki+1] >= sval2[1] &&
				  sp1inf[ki+3] > DZERO) ||
				 (sp1inf[ki  ] <= sval1[0] &&
				  sp1inf[ki+2] < DZERO) ||
				 (sp1inf[ki  ] >= sval1[1] &&
				  sp1inf[ki+2] > DZERO) ||
				 (sp2inf[ki+1] <= sval4[0] &&
				  sp2inf[ki+3] < DZERO) ||
				 (sp2inf[ki+1] >= sval4[1] &&
				  sp2inf[ki+3] > DZERO) ||
				 (sp2inf[ki  ] <= sval3[0] &&
				  sp2inf[ki+2] < DZERO) ||
				 (sp2inf[ki  ] >= sval3[1] &&
				  sp2inf[ki+2] > DZERO))
				{
				  knbinf++;
				  goto nextdir;
				}
			      else if (krem == 2 &&
				       ((sp1inf[ki+1] <= sval2[0] &&
					 sp1inf[ki+3] >= DZERO) ||
					(sp1inf[ki+1] >= sval2[1] &&
					 sp1inf[ki+3] <= DZERO) ||
					(sp1inf[ki  ] <= sval1[0] &&
					 sp1inf[ki+2] >= DZERO) ||
					(sp1inf[ki  ] >= sval1[1] &&
					 sp1inf[ki+2] <= DZERO) ||
					(sp2inf[ki+1] <= sval4[0] &&
					 sp2inf[ki+3] >= DZERO) ||
					(sp2inf[ki+1] >= sval4[1] &&
					 sp2inf[ki+3] <= DZERO) ||
					(sp2inf[ki  ] <= sval3[0] &&
					 sp2inf[ki+2] >= DZERO) ||
					(sp2inf[ki  ] >= sval3[1] &&
					 sp2inf[ki+2] <= DZERO)))
				{
				  /* We were marching out of the patch
				     but the tangent
				     is pointing in half step length */
				  kstpch = 0;
				}



			    }
			}
		      else
			{
			  /* Divergence or point on wrong side in the parameter
			     plane or 3-d */
			  kstpch = 0;
			  koutside_resolution = 0;
			}
		    }
		  else if (kstat==4)
		    goto nextdir;
		}

	      /* Update step length if new endpoint is to be found */

	      if (kstpch==0)
                {
		  if (tdist<DZERO)
                    {
		      tnew = tstep/(double)10.0;
                    }
		  else
                    {
		      tfak = MAX(tdist/aepsge,(double)1.0);
		      tfak = (double)2.0*pow(tfak,ONE_FOURTH);
		      tnew = MIN(tstep/(double)2.0,tstep/tfak);
                    }
		  if (DEQUAL(tmax+tnew,tmax+tstep)) goto nextdir;
		  tstep = tnew;
                }
            }

	  /* If kstpch= -1,1,3 or 4 then a point is accepted and
	     snxt1 points to the position and derivatives
	     of the accepted point. */


	  /* Update number of intersection points */

	  knbinf++;

	  /* Copy point and parameter pair descriptions */

	  memcopy(spnt1,spntend1,21,DOUBLE);
	  memcopy(spar1,sparend1,2,DOUBLE);
	  memcopy(spnt2,spntend2,21,DOUBLE);
	  memcopy(spar2,sparend2,2,DOUBLE);

	  /* Update guide point pointers */


	  if (kstpch ==  1)
            {
	      kguide++;
	      kgdir   = 1;

	      /* Test if end of guide point array reached */

	      if (kguide >= kpoint) goto nextdir;

            }
	  if (kstpch == -1)
            {
	      kguide--;
	      kgdir   = -1;

	      /* Test if start of guide point array reached */

	      if (1 >= kguide) goto nextdir;
            }

	  /* Make new radius of curvature */

	  trad = *(s3dinf + 10*knbinf - 1);
	  tstep = s1311(trad,aepsge,tmax,&kstat);
	  if (kstat<0) goto error;

	  /* Test if start point reached, e.g. that the curve is closed */

	  if (kstpch == 3)
            {
	      /* Closed curve found */

	      goto finished;
            }


	  /*  End while loop */
        }

    nextdir:;

      /*  End two step directions */
    }

 finished:

  /* In certain cases too many marched point may be found. These cases are:

     - Open curve and start of marching first guide point
     - Open curve and start of marching last guide point
     - Closed curve and this found in second marching direction

     In these cases some of the found points have to be discarded */

  scorpnt = s3dinf;
  scorpr1 = sp1inf;
  scorpr2 = sp2inf;

  if (kstpch !=3 && kpoint>1)
    {

      /*  Open curve */

      if ( (kstart==1 && kgd1 == kpoint) ||
	  (kstart==kpoint && kgd1==1)      )
        {
	  /* First marching direction traced curve */

	  knbinf = knb1;
        }
      else if ( (kstart==1 && kguide==kpoint) ||
	       (kstart==kpoint && kguide==1)    )
        {
	  /* Second marching direction traced curve */

	  scorpnt = scorpnt + 10*(knb1-1);
	  scorpr1 = scorpr1 +  7*(knb1-1);
	  scorpr2 = scorpr2 +  7*(knb1-1);
	  knbinf  = knbinf - knb1 + 1;
        }
    }
  else if (kpoint>1)
    {
      /*  Closed curve, correct if result of second marching direction */

      if (kdir != 1)
        {
	  /* Second marching direction, disc ard result of first direction */

	  scorpnt = scorpnt + 10*(knb1-1);
	  scorpr1 = scorpr1 +  7*(knb1-1);
	  scorpr2 = scorpr2 +  7*(knb1-1);
	  knbinf  = knbinf - knb1 + 1;
        }
    }

 interpolate:

  if (knbinf>1)
    {
      if (igraph == 1 && knbinf > 1)
	{
	  /* Output curve through s6line and s6move */

	  s6move(scorpnt);
	  for (ki=1,sp1=scorpnt+10;ki<knbinf;ki++,sp1+=10)
	    s6line(sp1);
	}

      /* A curve is traced out only if at least two points are found */

      if (icur > 0 && knbinf > 1)
	{

	  /*  Make 3-D representation of intersection curve */

	  kpar = 0;

	  spar = newarray(knbinf,DOUBLE);
	  if (spar == SISL_NULL) goto err101;
	  s1359(scorpnt,aepsge,kdim,knbinf,kpar,spar,&q3dcur,&kstat);
	  if (kstat < 0) goto error;

	  /*  Set pointer in intcurve object to 3-D curve */

	  pinter -> pgeom = q3dcur;

	  if (icur == 2)
	    {
	      /* Make curves in parameter planes */

	      kdim = 2;
	      kpar = 1;
	      s1359(scorpr1,aepsge,kdim,knbinf,kpar,spar,&qp1cur,&kstat);
	      if (kstat < 0) goto error;


	      s1359(scorpr2,aepsge,kdim,knbinf,kpar,spar,&qp2cur,&kstat);
	      if (kstat < 0) goto error;

	      /* Set pointers in intcurve object to curves in parameter plane*/

	      pinter -> ppar1 = qp1cur;
	      pinter -> ppar2 = qp2cur;
	    }
	}
    }
  else if( pinter->ipoint > 1)
    {
      /* If no points produced on intersection curve */

       s1310_s9constline(psurf1,psurf2,pinter,aepsge,icur,igraph,&kstat);
      if (kstat<0) goto error;
      if (kstat==0) goto err185;
    }
  else
    goto err185;

  if (kdiv == 1) goto war03;
  *jstat = 0;

  goto out;

  /* Iteration can not continue */
 war03:  *jstat = 3;
  goto out;

  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1310",*jstat,kpos);
  goto out;


  /* Error in surface description parameter direction does not exist */
 err123: *jstat = -123;
  s6err("s1310",*jstat,kpos);
  goto out;


/* Error - SISL_NULL pointer was given */
  err150 :
    *jstat = -150;
    s6err("s1310",*jstat,kpos);
    goto out;

  /* Only degenerate or singular guide points */
 err185: *jstat = -185;
  goto out;

  /* Error in lower leve function */
 error:
  *jstat = kstat;
  s6err("s1310",*jstat,kpos);
  goto out;

 out:

  /* Free allocated space */

  if (sgd1   != SISL_NULL) freearray(sgd1);
  if (sgd2   != SISL_NULL) freearray(sgd2);
  if (s3dinf != SISL_NULL) freearray(s3dinf);
  if (sp1inf != SISL_NULL) freearray(sp1inf);
  if (sp2inf != SISL_NULL) freearray(sp2inf);
  if (spar   != SISL_NULL) freearray(spar);
  if (sgpar1 != SISL_NULL) freearray(sgpar1);
  if (sgpar2 != SISL_NULL) freearray(sgpar2);


  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void s1310_s9constline(SISLSurf *ps1,SISLSurf *ps2,SISLIntcurve *pintcr,
			      double aepsge,int icur,int igraph,int *jstat)
#else
static void s1310_s9constline(ps1,ps2,pintcr,aepsge,icur,igraph,jstat)
     SISLSurf     *ps1;
     SISLSurf     *ps2;
     SISLIntcurve *pintcr;
     double   aepsge;
     int      icur;
     int      igraph;
     int      *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To check if the parameter pairs describe an intersection
*              curve that is a constant parameter line in the parameter
*              plane of a surface and to produce the description of
*              the curve according to the specifications.
*
*
* INPUT      : ps1    - Pointer to first surface.
*              ps2    - Pointer to second surface.
*              aepsge - Geometry resolution.
*              icur   - Indicator telling if a 3-D curve is to be made
*                        0 - Don't make 3-D curve
*                        1 - Make 3-D curve
*                        2 - Make 3-D curve and curves in parameter plane
*              igraph - Indicator telling if the curve is to be outputted
*                       through function calls:
*                        0 - don't output curve through function call
*                        1 - output as straight line segments through
*                            s6move and s6line
*
*
*
* INPUT/OUTPUT:pintcr - The intersection curve. When comming as input
*                       only parameter values in the parameter plane
*                       exist. When comming as output the 3-D geometry
*                       and possibly the curve in the parameter plane
*                       of the surface are added.
*
* OUTPUT:      jstat  - status messages
*                         = 1      : Constant parameter line is intersection
*                         = 0      : No intersection along constant parameter
*                                    line.
*                         < 0      : error
*
*
* METHOD     :
* REFERENCES :
*
*
*-
* CALLS      :
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 2. July 1989
*
*********************************************************************
*/
{
  int kguide1,kguide2,kguide3,kguide4; /* Pointers to guide points       */
  /*int kderc=2;         Number of derivatives to be claculated on curve */
  int kders=1;        /* Number of derivatives to be calculated on surface*/
  int ki,kj,kl;            /* Control variables in for loops            */
  int kk,kn,kk1,kn1,kk2,kn2;/* Orders and numbers of knots               */
  int       kk3,kn3,kk4,kn4;/* Orders and numbers of knots               */
  int kpoint;              /* Number of points in guide curve           */
  int kleft = 0;           /* Pointer into knot vector                  */
  int kpar1;               /* Number of parameter direction in 1st. obj */
  int kpar2;               /* Number of parameter direction in 2st. obj */
  int ktype;               /* Type of intersection curve                */
  int kpos=0;              /* Position of error                         */
  int kstat;               /* Status variable returned form routine     */
  int kdir=0;              /* constant parameter line direction         */
  int knbpnt;              /* Number of points on constant parameter line */
  int kleft1=0,kleft2=0;   /* Pointers into knot vectors                */
  int kstop;               /* Stop value in loop                        */
  double *sp=SISL_NULL;         /* Array for storage of points in
			      parameter plane */
  double *sv=SISL_NULL;         /* Array for storage of tangents in
			      parameter plane*/
  double *spar=SISL_NULL;       /* Array for storage of parameter values     */
  double *stp,*stv,*stpar; /* Pointers to sp,sv and spar                */
  double tdistp,tdistc;    /* Distances between points                  */
  double tfak;             /* Scaling factor                            */
  double sstart[4];        /* Lower boundary of parameter intervals     */
  double send[4];          /* Upper bounadry of parameter intervals     */
  double snext[3];         /* Existing iteration point on  surface      */
  double tmax1,tmin1;      /* Minimum and maximum of 1.rst comp of
			      guide points */
  double tmax2,tmin2;      /* Minimum and maximum of 2.nd. comp of
			      guide points */
  double tmax3,tmin3;      /* Minimum and maximum of 3.rd. comp of
			      guide points */
  double tmax4,tmin4;      /* Minimum and maximum of 4.th  comp of
			      guide points */
  double tmax;             /* Maximum 3-D SISLbox side                      */
  double tdist,tang;       /* Distance and angle error                  */
  double *st,*st1,*st2;    /* Pointers to knot vectors                  */
  double     *st3,*st4;    /* Pointers to knot vectors                  */
  double *spoint;          /* Pointer to points on constant parameter line */
  double *sp1;             /* Pointer into array                        */
  double tsize1,tsize2;    /* Length of knot intervals                  */
  double tsize3,tsize4;    /* Length of knot intervals                  */
  double sval1[2];         /* Limits of parameter plane in first SISLdir    */
  double sval2[2];         /* Limits of parameter plane in second SISLdir   */
  double sval3[2];         /* Limits of parameter plane in third SISLdir    */
  double sval4[2];         /* Limits of parameter plane in fourth SISLdir   */
  double sderc[6];         /* SISLPoint and derivative on curve             */
  double sders[9];         /* SISLPoint and derivative on curve             */
  double snorm[3];         /* Normal on implicit surface                */
  double tsumold,tsum,tval;/* Parameter values                    */
  double ta11,ta12,ta22;   /* Coefficients in equation system           */
  double tb1,tb2,tdom;     /* Left side of eq.syst. and determinant     */
  double t1,t2;            /* Derivatives in parameter plane            */
  double *sgpar1=SISL_NULL;     /* Parameter pairs of guide point in surf 1  */
  double *sgpar2=SISL_NULL;     /* Parameter pairs of guide point in surf 2  */
  SISLCurve *qc1=SISL_NULL;  /* Pointer to 3-D curve                     */
  SISLCurve *qc2=SISL_NULL;  /* Pointer to 3-D curve                     */

  SISLCurve *qp1cur=SISL_NULL;/* Pointer to curve in first parameter plane*/
  SISLCurve *qp2cur=SISL_NULL;/* Pointer to curve in second parameter plane*/
  SISLSurf  *qsurf =SISL_NULL;
  SISLPoint *qpoint=SISL_NULL;


  /* Make maximal step length based on box-size of surface */

  sh1992su(ps1,0,aepsge,&kstat);
  if (kstat < 0) goto error;

  tmax = MAX(ps1->pbox->e2max[0][0] - ps1->pbox->e2min[0][0],
	     ps1->pbox->e2max[0][1] - ps1->pbox->e2min[0][1]);
  tmax = MAX(tmax,ps1->pbox->e2max[0][2] - ps1->pbox->e2min[0][2]);

  sh1992su(ps2,0,aepsge,&kstat);
  if (kstat < 0) goto error;

  tmax = MAX(tmax,ps1->pbox->e2max[0][0] - ps1->pbox->e2min[0][0]);
  tmax = MAX(tmax,ps1->pbox->e2max[0][1] - ps1->pbox->e2min[0][1]);
  tmax = MAX(tmax,ps1->pbox->e2max[0][2] - ps1->pbox->e2min[0][2]);

  /* Find a none singular start point for the marching process */

  kpoint = pintcr->ipoint;
  kpar1  = pintcr->ipar1;
  kpar2  = pintcr->ipar2;
  sgpar1 = pintcr->epar1;
  sgpar2 = pintcr->epar2;
  ktype  = pintcr->itype;


  /* Initiate pointers to intersection curve and intersection curve in
     parameter plane */

  pintcr -> pgeom = SISL_NULL;
  pintcr -> ppar1 = SISL_NULL;
  pintcr -> ppar2 = SISL_NULL;



  /* Test that both objects has 2 parameter direction */

  if (kpar1 == 2 && kpar2 == 2)
    {
      /*  Everything is ok */
      ;
    }
  else
    {
      goto err123;
    }


  /* Run through the parameter pairs to decide if a constant parameter line
     is possible */

  tmax1 = tmin1 = sgpar1[0];
  tmax2 = tmin2 = sgpar1[1];
  tmax3 = tmin3 = sgpar2[0];
  tmax4 = tmin4 = sgpar2[1];

  /* Remember which guide point have minimum value in a specific parameter
     direction */

  kguide1 = kguide2 = kguide3 = kguide4 = 0;

  for (ki=1,kj=2,kl=3 ; ki < kpoint ; ki++,kj+=2,kl+=2)
    {
      if (tmin1>sgpar1[kj])
        {
	  tmin1 = sgpar1[kj];
	  kguide1 = ki;
        }
      if (tmin2>sgpar1[kl])
        {
	  tmin2 = sgpar1[kl];
	  kguide2 = ki;
        }
      if (tmin3>sgpar2[kj])
        {
	  tmin3 = sgpar2[kj];
	  kguide3 = ki;
        }
      if (tmin4>sgpar2[kl])
        {
	  tmin4 = sgpar2[kl];
	  kguide4 = ki;
        }
      tmax1 = MAX(tmax1,sgpar1[kj]);
      tmax2 = MAX(tmax2,sgpar1[kl]);
      tmax3 = MAX(tmax3,sgpar2[kj]);
      tmax4 = MAX(tmax4,sgpar2[kl]);
    }

  /* Initiate parameter direction boundaries */
  kk1    = ps1 -> ik1;
  kn1    = ps1 -> in1;
  st1    = ps1 -> et1;
  sval1[0] = st1[kk1-1];
  sval1[1] = st1[kn1];
  kk2    = ps1 -> ik2;
  kn2    = ps1 -> in2;
  st2    = ps1 -> et2;
  sval2[0] = st2[kk2-1];
  sval2[1] = st2[kn2];

  /* Initiate parameter direction boundaries */
  kk3    = ps2 -> ik1;
  kn3    = ps2 -> in1;
  st3    = ps2 -> et1;
  sval3[0] = st3[kk3-1];
  sval3[1] = st3[kn3];
  kk4    = ps2 -> ik2;
  kn4    = ps2 -> in2;
  st4    = ps2 -> et2;
  sval4[0] = st4[kk4-1];
  sval4[1] = st4[kn4];

  tsize1 = st1[kn1] - st1[kk1-1];
  tsize2 = st2[kn2] - st2[kk2-1];
  tsize3 = st3[kn3] - st3[kk3-1];
  tsize4 = st4[kn4] - st4[kk4-1];

  /* Check if constant parameter value within tolerance */

  if (DEQUAL((tmin1+tsize1),(tmax1+tsize1)) )
    {
      /* Intersection possible constant parameter line with first parameter
	 constant value constant.

	 1. Pick out curve from surface
	 2. Pick out relevant part of curve */

      kdir = 1;

      s1437(ps1,((tmin1+tmax1)/(double)2.0),&qc1,&kstat);
      if (kstat < 0) goto error;

      s1712(qc1,tmin2,tmax2,&qc2,&kstat);
      if (kstat < 0) goto error;

      /* Copy start point of iteration in surface */

      memcopy(snext,sgpar2+2*kguide1,2,DOUBLE);
    }
  else if (DEQUAL((tmin2+tsize2),(tmax2+tsize2)) )
    {
      /* Intersection possible constant parameter line with first parameter
	 constant value constant.

	 1. Pick out curve from surface
	 2. Pick out relevant part of curve
	 */

      kdir = 2;

      s1436(ps1,((tmin2+tmax2)/(double)2.0),&qc1,&kstat);
      if (kstat < 0) goto error;

      s1712(qc1,tmin1,tmax1,&qc2,&kstat);
      if (kstat < 0) goto error;

      /* Copy start point of iteration in surface */

      memcopy(snext,sgpar2+2*kguide2,2,DOUBLE);

    }
  else if (DEQUAL((tmin3+tsize3),(tmax3+tsize3)) )
    {
      /* Intersection possible constant parameter line with first parameter
	 constant value constant from second surface.

	 1. Pick out curve from surface
	 2. Pick out relevant part of curve */

      kdir = 3;

      s1437(ps2,((tmin3+tmax3)/(double)2.0),&qc1,&kstat);
      if (kstat < 0) goto error;

      s1712(qc1,tmin4,tmax4,&qc2,&kstat);
      if (kstat < 0) goto error;

      /* Copy start point of iteration in surface */

      memcopy(snext,sgpar1+2*kguide3,2,DOUBLE);

    }
  else if (DEQUAL((tmin4+tsize4),(tmax4+tsize4)) )
    {
      /* Intersection possible constant parameter line with first parameter
	 constant value constant.

	 1. Pick out curve from surface
	 2. Pick out relevant part of curve */

      kdir = 4;

      s1436(ps2,((tmin4+tmax4)/(double)2.0),&qc1,&kstat);
      if (kstat < 0) goto error;

      s1712(qc1,tmin3,tmax3,&qc2,&kstat);
      if (kstat < 0) goto error;

      /* Copy start point of iteration in surface */

      memcopy(snext,sgpar1+2*kguide4,2,DOUBLE);

    }
  else
    goto war00;

  st = qc2 -> et;
  kk = qc2 -> ik;
  kn = qc2 -> in;

  /* Set boundaries of surface to be iterated in */

  if (kdir==3 || kdir ==4)
    {
      sstart[0] = sval1[0];
      sstart[1] = sval2[0];
      send[0]   = sval1[1];
      send[1]   = sval2[1];
      qsurf = ps1;

    }
  else
    {
      sstart[0] = sval3[0];
      sstart[1] = sval4[0];
      send[0]   = sval3[1];
      send[1]   = sval4[1];
      qsurf = ps2;
    }

  /* Allocate array for storage of points, tangents and parameter values of
     curve in the parameter plane of the surface we test */

  if ((sp=newarray(4*kn,DOUBLE)) == SISL_NULL) goto err101;
  if ((sv=newarray(4*kn,DOUBLE)) == SISL_NULL) goto err101;
  if ((spar=newarray(2*kn,DOUBLE)) == SISL_NULL) goto err101;


  /* Run through 2*kn points of the curve and check that they lie in the
     implicit surface by calculating the 2*kn points. */

  tsumold = st[kk-1];

  for (ki=0,stp=sp,stv=sv,stpar=spar ; ki <kn ; ki++)
    {
      if (kk>1)
        {
	  /* Make parameter value to use for calculation of curve point */

	  for (kl=1,kj=ki+1,tsum=DZERO ; kl<kk ; kl++)
            tsum += st[kj++];

	  tsum = tsum/(double)(kk-1);
        }
      else
        tsum = st[ki];

      tval = (tsum+tsumold)/(double)2.0;

      for (kj=0 ; kj<2 ; kj++,stp+=2,stv+=2,stpar++)
        {

	  /* Calculate point on curve */

	  s1221(qc2,1,tval,&kleft,sderc,&kstat);
	  if (kstat < 0) goto error;

	  /* Remember the parameter value */

	  *stpar = tval;

	  /* Find closest point on surface to sderc */

	  qpoint = newPoint(sderc,3,0);
	  if (qpoint==SISL_NULL) goto err101;

	  /* Calculate closest point to surface */

	  s1773(qpoint,qsurf,aepsge,sstart,send,snext,stp,&kstat);
	  if(kstat<0) goto error;

	  freePoint(qpoint);

	  /* Calculate point and derivatives in surface */

	  s1421(qsurf,kders,stp,&kleft1,&kleft2,sders,snorm,&kstat);

	  if (kstat<0) goto error;

	  /* Find tangent of curve in parameter plane */

	  ta11 = s6scpr(sders+3,sders+3,3);
	  ta12 = s6scpr(sders+3,sders+6,3);
	  ta22 = s6scpr(sders+6,sders+6,3);
	  tb1  = s6scpr(sders+3,sderc+3,3);
	  tb2  = s6scpr(sders+6,sderc+3,3);

	  tdom = ta11*ta22 - ta12*ta12;
	  if (tdom != DZERO)
            {
	      t1 = (ta22*tb1-ta12*tb2)/tdom;
	      t2 = (ta11*tb2-ta12*tb1)/tdom;

	      tdom = sqrt(t1*t1+t2*t2);
	      if (tdom != DZERO)
                {
		  t1 /= tdom;
		  t2 /= tdom;
                }
            }
	  else
            t1 = t2 = DZERO;

	  /* Remember the tangent */

	  *stv     = t1;
	  *(stv+1) = t2;


	  /*Both the position of the two points should be within the relative
	    computer resolution for the point to be accepted.
	    Correspondingly the direction of the intersection curve and the
	    constant parameter line should be within the computer
	    resolution to be accepted. */

	  tdist = s6dist(sders,sderc,3);

	  if (DNEQUAL(tdist+tmax,tmax))
            goto war00;

	  /* Distance within tolerance, check that the angle between surface
	     normal and curve tangent is PIHALF, if both these vectors have a
	     nonzero length. */

	  if (s6length(snorm,3,&kstat) != DZERO &&
	      s6length(sderc+3,3,&kstat) != DZERO  )
	    {
	      tang = s6ang(snorm,sderc+3,3);
	      if (DNEQUAL(fabs(tang),PIHALF) ) goto war00;
	    }
	  tval = tsum;

	  /* Remember start point of iteration */

	  memcopy(snext,stp,2,DOUBLE);
        }
      tsumold = tsum;
    }

  /* Intersection curve along constant parameter line, make right actions
     concerning drawing and/or creation of the curve */

  if (igraph == 1)
    {
      /* Draw curve, first break into straight line segments */

      s1605(qc2,aepsge,&spoint,&knbpnt,&kstat);
      if (kstat < 0) goto error;

      if (knbpnt>1)
        {
	  /* Draw curve */

	  s6move(spoint);
	  for (ki=1,sp1=spoint+3 ; ki<knbpnt ; ki++,sp1+=3)
            s6line(sp1);
        }
      freearray(spoint);
    }

  if (icur >= 1)
    {
      /* Set pointer to 3-D curve */

      pintcr -> pgeom = qc2;
      qc2 = SISL_NULL;
    }

  if (icur == 2)
    {
      /* Make curves in parameter planes */

      double svert[4],sknot[4];


      /* Adjust the tangent lengths to match distance between adjacent points,
	 remember that first and second points are equal and that first point
	 is not used futher on */

      tdistp = s6dist(sp+2,sp+4,2);
      *(sv+2) *= tdistp;
      *(sv+3) *= tdistp;

      for (ki=2,stp=sp+4,stv=sv+4,kstop=kn+kn-1 ; ki < kstop ;
	   ki++,stp+=2,stv+=2)
	{
	  tdistc = s6dist(stp,stp+2,2);
	  tfak = (tdistp+tdistc)/(double)2.0;
	  *stv     *= tfak,
	  *(stv+1) *= tfak;
	  tdistp = tdistc;
	}
      *stv     *= tdistp;
      *(stv+1) *= tdistp;


      /* The first parameter pair is doubly represented */

      stp = sp+2;
      stv = sv+2;
      stpar = spar+1;

      if (kdir==1)
        {
	  svert[0] = svert[2] = (tmin1+tmax1)/(double)2.0;
	  svert[1] = tmin2;
	  svert[3] = tmax2;
	  sknot[0] = sknot[1] = tmin2;
	  sknot[2] = sknot[3] = tmax2;
	  qp1cur = newCurve(2,2,sknot,svert,1,2,1);
	  if (qp1cur==SISL_NULL) goto err101;
	  s1379(stp,stv,stpar,2*kn-1,2,&qp2cur,&kstat);
	  if (kstat<0) goto error;
	}
      else if (kdir==2)
        {
	  svert[0] = tmin1;
	  svert[2] = tmax1;
	  svert[1] = svert[3] = (tmin2+tmax2)/(double)2.0;
	  sknot[0] = sknot[1] = tmin1;
	  sknot[2] = sknot[3] = tmax1;
	  qp1cur = newCurve(2,2,sknot,svert,1,2,1);
	  if (qp1cur==SISL_NULL) goto err101;
	  s1379(stp,stv,stpar,2*kn-1,2,&qp2cur,&kstat);
	  if (kstat<0) goto error;
        }
      else if (kdir==3)
        {
	  svert[0] = svert[2] = (tmin3+tmax3)/(double)2.0;
	  svert[1] = tmin4;
	  svert[3] = tmax4;
	  sknot[0] = sknot[1] = tmin4;
	  sknot[2] = sknot[3] = tmax4;
	  qp2cur = newCurve(2,2,sknot,svert,1,2,1);
	  if (qp2cur==SISL_NULL) goto err101;
	  s1379(stp,stv,stpar,2*kn-1,2,&qp1cur,&kstat);
	  if (kstat<0) goto error;
        }
      else
        {
	  svert[0] = tmin3;
	  svert[2] = tmax3;
	  svert[1] = svert[3] = (tmin4+tmax4)/(double)2.0;
	  sknot[0] = sknot[1] = tmin3;
	  sknot[2] = sknot[3] = tmax3;
	  qp2cur = newCurve(2,2,sknot,svert,1,2,1);
	  if (qp2cur==SISL_NULL) goto err101;
	  s1379(stp,stv,stpar,2*kn-1,2,&qp1cur,&kstat);
	  if (kstat<0) goto error;
        }
      pintcr -> ppar1 = qp1cur;
      pintcr -> ppar2 = qp2cur;
    }

  *jstat = 1;
  goto out;

  /* Iteration can not continue */
 war00:  *jstat = 0;
  goto out;


  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1310_s9constline",*jstat,kpos);
  goto out;

  /* Error in surface description parameter direction does not exist */
 err123: *jstat = -123;
  s6err("s1310_s9constline",*jstat,kpos);
  goto out;

  /* Error in lower leve function */
 error:
  *jstat = kstat;
  s6err("s1310_s9constline",*jstat,kpos);
  goto out;

 out:;
  if (qc1 != SISL_NULL) freeCurve(qc1);
  if (qc2 != SISL_NULL) freeCurve(qc2);
  if (sp  != SISL_NULL) freearray(sp);
  if (sv  != SISL_NULL) freearray(sv);
  if (spar  != SISL_NULL) freearray(spar);
}
