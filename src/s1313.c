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
 * $Id: s1313.c,v 1.12 2005-12-09 14:08:45 afr Exp $
 *
 */


#define S1313

#include "sislP.h"

/*
 * Forward declarations.
 * ---------------------
 */
#if defined(SISLNEEDPROTOTYPES)
static void s1313_s9constline(SISLSurf *,double [],int,double,
			      SISLIntcurve *,int,int,int *);
#else
static void s1313_s9constline();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
s1313(SISLSurf *ps1,double eimpli[],int ideg,double aepsco,double aepsge,
      double amax,SISLIntcurve *pintcr,int icur,int igraph,int *jstat)
#else
     void s1313(ps1,eimpli,ideg,aepsco,aepsge,amax,pintcr,icur,igraph,jstat)
     SISLSurf     *ps1;
     double   eimpli[];
     int      ideg;
     double   aepsco;
     double   aepsge;
     double   amax;
     SISLIntcurve *pintcr;
     int      icur;
     int      igraph;
     int      *jstat;
#endif
     /*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To march an intersection curve described by parameter pairs
*              in the intersection curve, a B-spline surface and an
*              implicit surface.
*
*
* INPUT      : ps1    - Pointer to surface.
*              eimpli - Description of the implicit surface
*              ideg   - The degree of the implicit surface
*                        ideg=1: Plane
*                        ideg=2; Quadric surface
*                        ideg=1001: Torus surface
*                        ideg=1003: Silhouette line parallel projection
*                        ideg=1004: Silhouette line perspective projection
*                        ideg=1005: Silhouette line circular projection
*              aepsco - Not used.
*              aepsge - Geometry resolution.
*              amax   - Not used.
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
*                       If the curves has already been generated in the
*                       topology part of the intersections, nothing will
*                       be done (i.e. not required).  This will be the
*                       case when the intersection curve represents a
*                       constant parameter line in the parmeter plane
*                       of the surface.
*
* OUTPUT:      jstat  - status messages
*                         = 3      : Iteration stopped due to singular
*                                    point or degenerate surface. A part
*                                    of intersection curve may have been
*                                    traced out. If no curve is traced out
*                                    the curve pointers in the Intcurve
*                                    object point to SISL_NULL.
*                         = 3      : Marching not succeded
*                         = 0      : ok
*                         < 0      : error
*                         = -185   : No points produced on intersection curve.
*
*
* METHOD     :
* REFERENCES :
*
* HOW TO EXTEND THIS FUNCTION TO TREAT NEW PROBLEMS.
*
* This function is built as a general function for treating the combination
* of an implicit description and a parametric surface, when introducing
* new problems to be solved this structure can be utilized:
*
* 1. Define the implicit degree of the problem. If it is a special problem
*    utilize ideg>1000.
* 2. Determine how many derivatives have to be calculated to support the
*    function. For ideg=1,2,1001 derivatives up two have been calculated.
*    for ideg=1003,1004,1005, derivatives up to and including 3 have been calculated.
*    The number of derivatives influence the local variables kder, ksize,
*    ksizem3 that are found in the functions:
*     s1306,s1309,s1313,s1331,s9clipimp,s9iterimp,s9adsimp,s9boundimp
*
* 3. The function s1309 branch on itype to decide the distance in the
*    current iteration step. This function has to be updated to support
*    the new type of problem. For ideg=1,2 and 1001 this function
*    currently calculates distance. For ideg=1003,1004,1005 it calculates an angle.
*
* 4. The function s1331 that calculates derivatives of the combination of
*    the implicit problem and the surface have to support the new problem
*
* 5. Look also into s9iterimp and s9boundimp to update the convergence
*    branching on problem type.
*
*
*-
* CALLS      :
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 2. July 1988
* Revised by : Tor Dokken, SI, Oslo, Norway, 22. January 1988,
*               Test for degeneracy and singularities included
* Revised by : Tor Dokken, SI, Oslo, Norway, March 1989
*              Automatic generation of maximal step length, improved
*              marching close to singularities
* Revised by : Correction of error testing
* Revised by : Mike Floater, SI, 1991-01
*                   Improved the routine for parallel silhouettes (ideg=1003) and
*                   added perspective and circular silhouettes (ideg=1004,ideg=1005)
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Dec. 1994.  Added check for
*              SISL_NULL 'pintcr' and to avoid re-generating the geometry when it has
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
  int kpar2;               /* Number of parameter direction in 2st. obj */
  int kpar;                /* Indicater tellin if S1359 shall make
			      parametrization or use parametrization
			      in spar                                   */
  int ktype;               /* Type of intersection curve                */
  int klfu=0;              /* Pointers into knot vectors                */
  int klfv=0;              /* Pointers into knot vectors                */
  int kder = 2;            /* Calculate up to second derivatives        */
  int kdim = 3;            /* The dimension of the space we work in     */
  int kfirst = 0;          /* Indicator telling if first guide point degenerate */
  int klast = 0;           /* Indicator telling if last guide point degenerate */
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
  int krem,krem1,krem2;    /* Remember if we step in or out of patch    */
  int kbound;              /* Dummy variiable                           */
  int koutside_resolution; /* Flag telling if current seg. outside res. */
  int ksize;               /* Number of doubles for storage of derivateves
			      and normal vector */
  int ksizem3;             /* ksize - 3                                 */
  int kdiv=0;              /* Flag remembering if iteration diverged    */
  int knb1=0;              /* Remember number of points after marching
			      in first marching direction               */
  int kgd1=0;              /* Remeber last guide point used in first
			      marching direction                        */
  double *scorpnt=SISL_NULL;    /* Corrected marching points                 */
  double *scorpar=SISL_NULL;    /* Corrected marching parameter values       */
  double smidd[6];         /* Description of midpoint and tangent of
			      current Bezier segment                    */
  double tcurstep;         /* Current step length                       */
  double tdist;            /* Error at middle of current Bezier segement*/
  double tang;             /* Angle error at midpoint Bezier segement   */
  double tnew;             /* Candidate for new step length             */
  double tfak;             /* How much is the step length to be reduced */
  double *start;           /* Pointer to start of current segment       */
  double *st;              /* Pointer to knot vector                    */
  double sval1[2];         /* Limits of parameter plane in first SISLdir    */
  double tref1,tref2;      /* Reference values for knot vectors         */
  double sval2[2];         /* Limits of parameter plane in second SISLdir   */
  double tstep;            /* Iteration step length                     */
  double tmax;             /* Local maximal step length                 */
  double tstartstp;        /* Start step length                         */
  double trad;             /* Radius of curvature                       */
  double tval[6];             /* Dummy array in s1331                  */
  double *spar=SISL_NULL;       /* Parametrization of points in Hermit interp*/
  double spar1[2];         /* Parameter pair of current point surface 1 */
  double spar2[2];         /* Parameter pair of boundarypoint surface 1 */
  double siparmid[2];      /* Parameter value at middle of Bezier segment*/
  double sipar1[2];        /* Parameter pair iteration point surface  1 */
  double simiddpnt[10];    /* Middle point and tangent of segment       */
  double simiddpar[7];     /* Parameter value at middle point of segment*/
  double *sgpar1=SISL_NULL;     /* Parameter pairs of guide point in surf 1  */
  double *sgpar2=SISL_NULL;     /* Parameter pairs of guide point in surf 2  */
  double *sgpar=SISL_NULL;      /* guide points used                         */
  double *sgd1 = SISL_NULL;     /* 0-2 derivative of guide point + normal
				   of first object                           */
  double spnt1[33];        /* Info on current point in first surface    */
  double spnt2[33];        /* Info on boundary point in first surface   */
  double sipnt1[33];       /* Info on iteration point in first surface  */
  /* For spnt1, sipnt1, the                    */
  /* information is stored 3-tuppels in the    */
  /* following sequence                        */
  /* Position, (1,0)-der, (0,1)-der,           */
  /* (2,0)-der, (1,1)-der, (0,2)-der and normal*/
  /* This is compatible with output of s1421   */
  double snorm1[3];        /* Normal vector of implicit surface         */
  double snorm2[3];        /* Normal vector of implicit surface         */
  double startg[3];        /* Start tangent of iteration                */

  double *snxt1;           /* SISLPoint in ps1 we have accepted          */
  double *snxp1;           /* Parameter value belonging to snxt1        */
  double *s3dinf=SISL_NULL;     /* Pointer to array used for storing 3-D position
				   tangent, curvature and radius of curvature found
				   during the marching process if possible */
  double *sp1inf=SISL_NULL;     /* Pointer to array used for storing position
				   tangent, curvature and radius of curvature found
				   in the first parameter plane during the
				   marching process */
  double start1[33];       /* Description of start point in ps1      */
  double stpar1[2];        /* Parameter pair belonging to start1        */
  double sdum1[3],sdum2[3];/* Dummy vectors                             */
  double tdum,tdump;       /* Dummy variable                            */
  double *sp1=SISL_NULL;        /* Pointer used when moving information      */
  double *sp2=SISL_NULL;        /* Pointer used when moving information      */
  double stdum[10];        /* Dummy array used when moving information  */
  double *stang;           /* Pointer to tangent of current point       */
  double *sptang;          /* Pointer to tangent in parameter plane     */
  double *spoint;          /* Pointer to current point                  */
  double t1distgd,t2distgd;/* Distances to guide points                 */
  double tlnorm;           /* Length of normal vector                   */
  double tltan1,tltan2;    /* Tangent lengths                           */
  SISLCurve *q3dcur=SISL_NULL;/* Pointer to 3-D curve                     */
  SISLCurve *qp1cur=SISL_NULL;/* Pointer to curve in first parameter plane*/
  double sdiffcur[3];        /* Difference between current and previous point found */
  double sdiffprev[3];        /* Difference between previous point and the one before that */


  *jstat = 0;

  if (pintcr == SISL_NULL)  goto err150;


  /* Check if the geometry already has been generated in the topology part.
     This will be the case if the geometry is along a constant parameter line. */

  if (pintcr->itype == 9)  goto out;



  /* Make maximal step length based on box-size of surface */

  sh1992su(ps1,0,aepsge,&kstat);
  if (kstat < 0) goto error;

  tmax = MAX(ps1->pbox->e2max[0][0] - ps1->pbox->e2min[0][0],
	     ps1->pbox->e2max[0][1] - ps1->pbox->e2min[0][1]);
  tmax = MAX(tmax,ps1->pbox->e2max[0][2] - ps1->pbox->e2min[0][2]);

  if (amax>DZERO) tmax = MIN(tmax,amax);


  /* If ideg=1,2 or 1001 then only derivatives up to second order
     are calculated, then 18 doubles for derivatives and 3 for the
     normal vector are to be used for calculation of points in the
     spline surface. For ideg=1003,1004,1005 we have a silhouette curve and
     derivatives up to the third are to be calculated,
     thus 30 +3 a total of 33 doubles are to be calculated */

  if (ideg==1003 || ideg==1004 || ideg==1005)
    {
      kder = 3;
      ksize = 33;
    }
  else
    {
      ksize = 21;
      kder =2;
    }
  ksizem3 = ksize -3;


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


  /* Initiate parameter direction boundaries */
  kk    = ps1 -> ik1;
  kn    = ps1 -> in1;
  st    = ps1 -> et1;
  sval1[0] = st[kk-1];
  sval1[1] = st[kn];
  tref1 = (double)3.0*MAX(fabs(*sval1),fabs(*(sval1+1)));
  kk    = ps1 -> ik2;
  kn    = ps1 -> in2;
  st    = ps1 -> et2;
  sval2[0] = st[kk-1];
  sval2[1] = st[kn];
  tref2 = (double)3.0*MAX(fabs(*sval2),fabs(*(sval2+1)));

  /* Test that first object has 2 parameter direction and second object 0*/

  if (kpar1 == 2 && kpar2 == 0)
    {
      /*  Everithing is ok */
      ;
    }
  else if (kpar1 == 0 && kpar2 == 2)
    {
      sgpar1 = sgpar2;
    }
  else
    {
      goto err123;
    }

  /* To support closed curve the first guide point must be copied after
     the last guide point */

  if((sgpar=newarray(2*kpoint+2,DOUBLE)) == SISL_NULL) goto err101;
  memcopy(sgpar,sgpar1,2*kpoint,DOUBLE);
  if (ktype ==2 || ktype == 3)
    {
      /* Closed curve copy first guide point to end of string of guide points */
      memcopy(sgpar+2*kpoint,sgpar1,2,DOUBLE);
      kpoint = kpoint + 1;
    }

  /* THE POINTS , TANGENT, CURVATURE AND RADIUS OF CURVATURE FOUND DURING
   *  THE MARCHING PROCESS SHOULD ALL BE STORED IN ARRAYS. ALLOCATE ONE ARRAY
   *  FOR 3-D INFORMATION AND ONE ARRAY FOR INFORMATION IN THE PARAMETER PLANE.
   *  THESE ARRAYS ARE GIVEN AN INITIAL CAPACITY OF STORING 100 POINTS WITH
   *  OTHER INFORMATION.
   *  IF THEY ARE TOO SHORT THEY WILL BE REALLOCATED AT A LATER STAGE.
   *
   *  SINCE THE MARCHING WILL GO IN BOTH DIRECTIONS WE WILL HAVE TO TURN THE
   *  INFORMATION FOUND WHEN MARCHING IN NEGATIVE DIRECTION, SO THAT IT CAN
   *  BE COMBINED WITH THE INFORMATION FOUND WHEN WE ARE MARCHING IN POSITVE
   *  DIRECTION.
   */

  kmaxinf = 100;
  s3dinf = newarray(10*kmaxinf,DOUBLE);
  if (s3dinf == SISL_NULL) goto err101;
  sp1inf = newarray(7*kmaxinf,DOUBLE);
  if (sp1inf == SISL_NULL) goto err101;

  /* Evaluate 0-1-2nd. derivative + normal of all guide points in the surface,
     first allocate arrays for storing the information, check that the points
     have a defined normal, and that the combination of the implicit surface
     and the surface defines a tangent direction in the curve */

  sgd1 = newarray(ksize*kpoint,DOUBLE);
  if (sgd1==SISL_NULL) goto err101;

  kpos = 5;

  /* Initiate kstart to point at no point */

  kstart = 0;

  for (ki=0,kj=0,kl=0 ; ki<kpoint ; ki++,kj+=2,kl+=ksize)
    {
      s1421(ps1,kder,&sgpar[kj],&klfu,&klfv,&sgd1[kl],&sgd1[kl+ksizem3],&kstat);
      if (kstat<0) goto error;



      if (ideg == 1003 || ideg == 1004 || ideg == 1005)
	{
	  /*  Find length of normal vector and tangent vectors */

	  tlnorm = s6length(&sgd1[kl+ksizem3],kdim,&kstat);
	  tltan1 = s6length(&sgd1[kl+kdim],kdim,&kstat);
	  tltan2 = s6length(&sgd1[kl+kdim+kdim],kdim,&kstat);

	  /*  The cross product satisifes the follwing conditions:
	      length(axb) = length(a) length(b) sin(angle(a,b)).
	      Thus the angle between the two vectors can be found, close to 0
	      sin(a) is a good approximation of a
	  */
	  if (tlnorm == (double)0.0 || tltan1 ==(double)0.0 || tltan2 == (double)0.0)
	    tang = (double)0.0;
	  else
	    tang = tlnorm/(tltan1*tltan2);

	  /* Silhouette line calculation no normal can be found in implicit surface
	     accept point as candidate start point if tang greater tha angular
	     resolution */

	  if (tang >= ANGULAR_TOLERANCE) kstart = ki+1;
	}
      else
	{
	  /* Make direction of intersection curve */

	  s1306(&sgd1[kl],&sgpar[kj],eimpli,ideg,s3dinf,sp1inf,&kstat);
	  if (kstat < 0) goto error;


	  /* Remember if start, internal or end point */

	  if (kstat != 2)
	    {
	      if (ki == 0) kfirst = 1;
	      else if (ki == kpoint-1) klast = kpoint;
	      else  kstart = ki+1;
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
      /*  No internal guide point necessary, copy last point
	  to second point */

      memcopy(sgd1+ksize,sgd1+ksize*(kpoint-1),ksize,DOUBLE);
      memcopy(sgpar+2,sgpar+2*(kpoint-1),2,DOUBLE);

      if (kstart ==  kpoint) kstart = 2;
      kpoint = 2;
    }
  else if (kpoint>2)
    {
      /*Internal guide point exists, copy this to second position and
	copy end point to third position */

      memcopy(sgd1+ksize,sgd1+ksize*(kstart-1),ksize,DOUBLE);
      memcopy(sgpar+2,sgpar+2*(kstart-1),2,DOUBLE);

      memcopy(sgd1+2*ksize,sgd1+ksize*(kpoint-1),ksize,DOUBLE);
      memcopy(sgpar+4,sgpar+2*(kpoint-1),2,DOUBLE);

      kpoint = 3;
      kstart = 2;
    }



  /* Remember description of start point in both surfaces,
   *  copy point indicated by kstart into spnt1,spar1 */

  memcopy(spnt1,sgd1+ksize*(kstart-1),ksize,DOUBLE);
  memcopy(spar1,sgpar+2*(kstart-1),2,DOUBLE);

  /* Make position, unit tangent, curvature and radius of curvature for
   *  start point of iteration, store them in the arrays just allocated */

  kpos = 10;
  s1306(spnt1,spar1,eimpli,ideg,s3dinf,sp1inf,&kstat);

  if (kstat<0) goto error;

  /* Test if singular point reached */

  if (kstat == 2) goto war03;

  /* Remember start tangent */
  memcopy(startg,s3dinf+3,3,DOUBLE);


  /* Iterate intersection point down to the intersection curve */

  tstep = DZERO;
  s9iterimp(s3dinf,spnt1,spar1,ps1,eimpli,ideg,tstep,
	    aepsge,sipnt1,sipar1,&kstat);
  if (kstat < 0) goto error;

  if (kstat==0 && s6dist(spnt1,sipnt1,3) > aepsge)
    {

      /* Nonsingular point found and adjustment greater than aepsge,
	 copy result of iteration into spnt1,spar1 */

      memcopy(spnt1,sipnt1,ksize,DOUBLE);
      memcopy(spar1,sipar1,2,DOUBLE);
    }

  /* Remember start point */

  if (kstat==0)
    {
      memcopy(start1,sipnt1,ksize,DOUBLE);
      memcopy(stpar1,sipar1,2,DOUBLE);
    }
  else
    {
      memcopy(start1,spnt1,ksize,DOUBLE);
      memcopy(stpar1,spar1,2,DOUBLE);
    }

  /* Make position, unit tangent, curvature and radius of curvature for
     start point of iteration, store them in the arrays just allocated */

  kpos = 10;
  s1306(start1,stpar1,eimpli,ideg,s3dinf,sp1inf,&kstat);

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

	  /* Then interchange info in parameter plane */

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
/* 	      sp1[2] = - sp1[2]; */
	    }

	  /* Turn direction of remembered start tangent */

	  for (ki=0;ki<3;ki++)
	    startg[ki] = -startg[ki];


	  /* Update spnt1 and  spar1 to have the start point values */

	  memcopy(spnt1,start1,ksize,DOUBLE);
	  memcopy(spar1,stpar1,2,DOUBLE);

	  /* Turn the direction we march the guide point vector, and set current
	     guide point to kstart */

	  kgdir  = -kgdir;
	  kguide = kstart;

	  /*      Update step length */
	  tstep = tstartstp;
	}


      stang = s3dinf + 10*(knbinf-1) + 3;
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
	   *   kstpch = -1 : Try to step to previous guide point
	   *   kstpch =  0 : Try not to step through guide point
	   *   kstpch =  1 : Try to step to next guide point
	   *   kstpch =  3 : Step to start point and stop marching
	   *   kstpch =  4 : Don't step through guide point, candidate
	   *                 end point of segement found in iteration loop
	   */


	  kstpch = 0;
	  stang = s3dinf + 10*(knbinf-1) + 3;
	  sptang = sp1inf + 7*(knbinf-1) + 2;
	  if (kgdir >=0)
	    {

	      /*  We are stepping in positive direction in guide point vector
	       *  calculate distance to next guide point. If the guide point
	       *  is lying closer than the step length to the current point
	       *  we should step directly to this point provided that the cross
	       *  product of the normal vectors at current point and at the
	       *  guide point have the same direction, e.g. that their scalar
	       *  product is positiv                          */

	      t1distgd = (double)2.0*tstep;
	      if (kguide < kpoint)
		{
		  /* Decide if we should step through the guide point */

		  kpos = 40;
		  t1distgd = s9adsimp(spnt1,spar1,eimpli,ideg,
				      &sgd1[kguide*ksize],
				      &sgpar[kguide*2],
				      stang,sptang,tstep,&kstat);
		  if (kstat < 0) goto error;
		  if (kstat == 1)
		    {
		      /* Step through guide point remember this */

		      kstpch = 1;
		      snxt1 = sgd1 + ksize*kguide;
		      snxp1 = sgpar + 2*kguide;
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
	       * product is positiv                          */

	      if (1 < kguide)
		{
		  /* Decide if we should step through the guide point */

		  kpos = 50;
		  t2distgd = s9adsimp(spnt1,spar1,eimpli,ideg,
				      &sgd1[(kguide-2)*ksize],&sgpar[2*kguide-4],
				      stang,sptang,tstep,&kstat);
		  if (kstat<0) goto error;
		  if ((kstat == 1 &&kstpch == 0) ||
		      (kstat == 1 && kstpch == 1 && t2distgd < t1distgd))
		    {
		      /* Step through guide point remember this */

		      kstpch = -1;
		      snxt1 = sgd1 + ksize*(kguide-2);
		      snxp1 = sgpar + 2*(kguide-2);
		      tstep = MIN(tstep,t2distgd);
		    }
		}
	    }

	  /* Check if we step through the start point. Should only be done if at least
	     three points found in this marching direction, a full closed curve will
	     require 6 segments. */
	  if ((kdir==1 && knbinf>3) || (kdir==2 && knbinf>knb1+3))
	    {

	      kpos = 60;
	      tdum = s9adsimp(spnt1,spar1,eimpli,ideg,start1,stpar1,stang,sptang,
			      tstep,&kstat);
	      if (kstat < 0) goto error;
	      if (kstat == 1)
		{
		  /* Step to start point remember this */

		  kstpch = 3;
		  snxt1 = start1;
		  snxp1 = stpar1;
		  tstep = MIN(tstep,tdum);
		}
	    }


	  /* At this stage kstpch=0 if we have not reached a guide point or
	     if we have not reached the start point of the iteration.

	     Now we want to find a Bezier segement that is lying within the
	     geometric tolerance that is approximating the intersection curve.
	     If a guide point is reached (kstpch=-1 or 1), then we have a
	     candidate for the end point of the Bezier segement. If the start
	     point is reached (kstpch=3) then we also have a candidate end point
	     for the segment.

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
	    }

	  /* Make description of candidate endpoint if it exists and store it */

	  if (kstpch != 0)
	    {
	      s1306(snxt1,snxp1,eimpli,ideg,s3dinf+10*knbinf,
		    sp1inf+7*knbinf,&kstat);
	      if (kstat<0) goto error;

	      /* It is allowed to jump on to a singular point */

	      /* Make sure that the tangents of previous and the
		 candidate end point point in the same direction */

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
		}

	      /* Copy the candidate point to spnt2 and spar2 */

	      memcopy(spnt2,snxt1,ksize,DOUBLE);
	      memcopy(spar2,snxp1,2,DOUBLE);
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
		   * defined by current point (s3dinf), the tangent (s3dinf+3)
		   * and the step length */

		  spoint = s3dinf + 10*(knbinf-1);
		  tcurstep = tstep;
		}

	      /*  Perform the actual iteration */

	      kpos = 70;
	      s9iterimp(spoint,spnt1,spar1,ps1,eimpli,ideg,tcurstep,
			aepsge,sipnt1,sipar1,&kstat);
	      if (kstat < 0) goto error;

	      /*             Initiate distance between midpoint and iteration point
			     to -1 to enable detection of divergence */

	      tdist = -1;

	      /*  Check if iteration has converged        */

	      if (kstat == 2)
		{
		  /*  Iteration has diverged, half step length if possible,
		      find new endpoint of segement. */

		  kstpch = 0;
		  koutside_resolution = 0;
		}
	      else if(kstat == 1 && kstpch != 0)
		{
		  /* The point found is closer to the input point than
		     the relative computer resolution or is a singular point.
		     Half step length if possible, find new endpoint of
		     segement. */

		  kstpch = 0;
		  koutside_resolution = 0;
		}
	      else if (kstpch!=0)
		{


		  /*  Make description of intersection point */

		  s1306(sipnt1,sipar1,eimpli,ideg,simiddpnt,simiddpar,&kstat);
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
				     (fabs(tang) > ANGULAR_TOLERANCE && tstep>aepsge)))
		    {
		      kstpch = 0;
		      koutside_resolution = 0;
		    }
		  else
		    {
		      /* Segment within tolerance.
		       * Check that the relationship between the two surfaces
		       * has not been interchanged, by making the cross product
		       * of the normal vectors in current point and the point
		       * found by iteration. Then make the scalar product of
		       * these vectors. If the scalar product is negative then
		       * we have either jumped to another branch or passed a
		       * singularity, iterprete this as the iteration has diverged
		       * In addition we don't want the direction of the tangents
		       * change to much. We set a limit of approximately 41 degrees
		       * by testing on a cosin value of 0.75
		       * Make normal vectors in implicit surface for both points
		       * Make also sure that the curve in the parameter plane
		       * does not turn more than 90 degrees.
		       */

		      s1331(spnt1,eimpli,ideg,-1,tval,snorm1,&kstat);
		      if (kstat < 0) goto error;
		      s1331(spnt2,eimpli,ideg,-1,tval,snorm2,&kstat);
		      if (kstat < 0) goto error;

		      if(ideg == 1003 || ideg == 1004 || ideg == 1005)
			{
			  tdum = s6scpr(snorm1,snorm2,kdim);
			}
		      else
			{
			  s6crss(spnt1+ksizem3,snorm1,sdum1);
			  (void)s6norm(sdum1,kdim,sdum1,&kstat);
			  if (kstat < 0) goto error;
			  s6crss(sipnt1+ksizem3,snorm2,sdum2);
			  (void)s6norm(sdum2,kdim,sdum2,&kstat);
			  if (kstat < 0) goto error;
			  tdum = s6scpr(sdum1,sdum2,kdim);
			}

		      s6diff(spar2,spar1,2,sdum1);
		      tdump = s6scpr(sdum1,sp1inf+7*(knbinf-1)+2,2);
		      /* The test below was added to detect the case when the
			 normal of the parametric surface turns */
		      if(s6scpr(spnt1+ksizem3,sipnt1+ksizem3,kdim) < 0.0)
			{
			  tdum = -tdum;
			}
		      if (tdum == DZERO)
			{
			  double tl1,tl2;

			  /* If one of the tangents have zero length, accept
			     segment */

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
		      else if (tdum <= (double)0.75 || tdump <= (double)0.0)
			{
			  /* Find new end point of segment */
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
		  /* We iterated to find end point of segment, update pointer */

		  memcopy(spnt2,sipnt1,ksize,DOUBLE);
		  memcopy(spar2,sipar1,2,DOUBLE);

		  s1306(sipnt1,sipar1,eimpli,ideg,s3dinf+10*knbinf,
			sp1inf+7*knbinf,&kstat);
		  if (kstat<0) goto error;

		  /* Make sure that the tangents of previous and the new point
		     point in the same direction, singular end point allowed */

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
		    }

		  /* Indicate that end point accepted */

		  kstpch = 4;
		  koutside_resolution = 0;
		}

	      /* If the segment is accepted. check if we cross the
		 boundary of the patch  */

	      if (kstpch !=0 && koutside_resolution == 1)
		{

		  /*               Check if curve between start and endpoint cross the
				   boundary */

		  memcopy(siparmid,sipar1,2,double);

		  s1305(spar1,spar2,sval1,sval2,&kbound,sipar1,&kstat);
		  if (kstat<0) goto error;


		  if(kstat==0 || kstat==4)
		    {
		      /* Set pointer to tangents at start point */
		      ki = 7*(knbinf-1)+2;

		      if( (DEQUAL(spar1[1]+tref2,sval2[0]+tref2) && sp1inf[ki+1]>DZERO) ||
			  (DEQUAL(spar1[1]+tref2,sval2[1]+tref2) && sp1inf[ki+1]<DZERO) ||
			  (DEQUAL(spar1[0]+tref1,sval1[0]+tref1) && sp1inf[ki]>DZERO) ||
			  (DEQUAL(spar1[0]+tref1,sval1[1]+tref1) && sp1inf[ki]<DZERO)

			  )
			kstat = 1;
		    }
		  krem1 = kstat;



		  /*               Check if curve between start and  midd point cross the
				   boundary */

		  s1305(spar1,siparmid,sval1,sval2,&kbound,sipar1,&kstat);
		  if (kstat<0) goto error;
		  krem2 = kstat;

		  /* We now have the following cases:
		     kstat == 0 : Line between epar1 and epar2 outside,
		     If this happens when kdir=1, then
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
		  if (krem1==0 || krem2==0)
		    {
		      if (kdir==1) knbinf--;
		      goto nextdir;
		    }
		  else if ((krem1 !=1 || krem2 !=1) && krem1 != 4 && krem2 !=4)
		    {

		      /* If we clip to the boundary, forget any guide point
			 identified. The actual action is dependent on which
			 of the points spar2 or sipar1 indicates the crossing */

		      kstat1 = 0;
		      if (krem2==2 || krem2==3)
			{
			  s9clipimp(spar1,siparmid,ps1,eimpli,ideg,sval1,sval2,
				    aepsge,sipnt1,sipar1,&kstat);
			  if (kstat<0) goto error;
			  if (krem2==3 && kstat==1) kstpch = 4;
			  kstat1 = kstat;
			  krem = krem2;
			}
		      if (kstat1 != 1 && (krem1 ==2 || krem1==3))
			{
			  s9clipimp(spar1,spar2,ps1,eimpli,ideg,sval1,sval2,
				    aepsge,sipnt1,sipar1,&kstat);
			  if (kstat<0) goto error;
			  if (krem1==3 && kstat==1) kstpch = 4;
			  kstat1 = kstat;
			  krem = krem1;
			}

		      if (kstat1 == 1)
			{
			  /* Check that the relationship between the two surfaces
			   * has not been interchanged, by making the cross product
			   * of the normal vectors in current point and the point
			   * found by iteration. Then make the scalar product of
			   * these vectors. If the scalar product is negative then
			   * we have either jumped to another branch or passed a
			   * singularity, iterprete this as the iteration has diverged
			   * In addition we don't want the direction of the tangents
			   * change to much. We set a limit of approximately 41 degrees
			   * by testing on a cosin value of 0.75
			   * Make normal vectors in implicit surface for both points
			   * Make also sure that the curve in the parameter plane
			   * does not turn more than 90 degrees.
			   */

			  s1331(spnt1,eimpli,ideg,-1,tval,snorm1,&kstat);
			  if (kstat < 0) goto error;
			  s1331(sipnt1,eimpli,ideg,-1,tval,snorm2,&kstat);
			  if (kstat < 0) goto error;

			  if(ideg == 1003 || ideg == 1004 || ideg == 1005)
			    {
			      tdum = s6scpr(snorm1,snorm2,kdim);
			    }
			  else
			    {
			      s6crss(spnt1+ksizem3,snorm1,sdum1);
			      (void)s6norm(sdum1,kdim,sdum1,&kstat);
			      if (kstat < 0) goto error;
			      s6crss(sipnt1+ksizem3,snorm2,sdum2);
			      (void)s6norm(sdum2,kdim,sdum2,&kstat);
			      if (kstat < 0) goto error;

			      tdum = s6scpr(sdum1,sdum2,kdim);
			    }


			  /* Check that sipar1 lies on the same side of spar1 as
			     the tangent at spar1 */

			  s6diff(sipar1,spar1,2,sdum1);
			  tdump = s6scpr(sdum1,sp1inf+7*(knbinf-1)+2,2);
			  /* The test below was added to detect the case when the
			     normal of the parametric surface turns */
			  if(s6scpr(spnt1+ksizem3,sipnt1+ksizem3,kdim) < 0.0)
			    {
			      tdum = -tdum;
			    }
			}

		      /* An intersection point has only been found when kstat==1
		       */
		      if (kstat1==1 && tdump >= (double)0.0 && tdum > (double)0.75)
			{
			  /* If krem=3 we step into the patch, if krem=2 we step
			     out of the patch */

			  if (krem==2 || krem==3)
			    {
			      /* Since we clip, set kstpch=0, no guide point reached */

			      /* If krem==3 we step into the patch, make new
				 start point of segment */

			      if (krem==3) knbinf--;

			      memcopy(spnt2,sipnt1,ksize,DOUBLE);
			      memcopy(spar2,sipar1,2,DOUBLE);

			      s1306(sipnt1,sipar1,eimpli,ideg,s3dinf+10*knbinf,
				    sp1inf+7*knbinf,&kstat);
			      if (kstat<0) goto error;


			      /* Make sure that the tangents of previous and the new point
				 point in the same direction */

			      if (knbinf>0)
				tdum = s6scpr(s3dinf+10*(knbinf-1)+3,
					      s3dinf+10*knbinf+3,kdim);
			      else
				tdum = s6scpr(startg,s3dinf+3,kdim);


			      if (tdum < DZERO)
				{
				  /* Change tangent direction 3-D and in
				     parameter plane
				  */
				  sp1 = s3dinf + 10*knbinf + 3;
				  sp1[0] = -sp1[0];
				  sp1[1] = -sp1[1];
				  sp1[2] = -sp1[2];
				  sp1 = sp1inf + 7*knbinf + 2;
				  sp1[0] = -sp1[0];
				  sp1[1] = -sp1[1];
				}

			      /* If the new end point tangent points out go
				 to next direction */

			      ki = 7*knbinf;
			      if( (sp1inf[ki+1] <= sval2[0] && sp1inf[ki+3] < DZERO) ||
				  (sp1inf[ki+1] >= sval2[1] && sp1inf[ki+3] > DZERO) ||
				  (sp1inf[ki  ] <= sval1[0] && sp1inf[ki+2] < DZERO) ||
				  (sp1inf[ki  ] >= sval1[1] && sp1inf[ki+2] > DZERO)   )
				{
				  knbinf++;
				  goto nextdir;
				}
			      else if(krem == 2 &&
				      ((sp1inf[ki+1] <= sval2[0] && sp1inf[ki+3] >= DZERO) ||
				       (sp1inf[ki+1] >= sval2[1] && sp1inf[ki+3] <= DZERO) ||
				       (sp1inf[ki  ] <= sval1[0] && sp1inf[ki+2] >= DZERO) ||
				       (sp1inf[ki  ] >= sval1[1] && sp1inf[ki+2] <=DZERO)    ))
				{
				  /* We were marching out ou the patch, but the tangent
				     points in, half step length */
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
		  else if (krem1==4 || krem2==4)
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
	   * snxt1 points to the position and derivatives
	   * of the accepted point. */

	  /* If we have accepted a segment pointing in the opposite
	   * direction of the previous segment, something very wrong
	   * has happened, and we go out with an error. */
	  if (knbinf >= 2)
	    {
	      s6diff(s3dinf + 10*knbinf, s3dinf + 10*(knbinf - 1), kdim, sdiffcur);
	      s6diff(s3dinf + 10*(knbinf-1), s3dinf + 10*(knbinf - 2), kdim, sdiffprev);
	      if (s6scpr(sdiffcur, sdiffprev ,kdim) < DZERO)
		{
		  /* We have a problem with degeneracy, quit now. */
		  goto war03;
		}
	      /* printf("%7.13f\n", s6scpr(sdiffcur, sdiffprev ,kdim)); */
	    }

	  /* Update number of intersection points */

	  knbinf++;

	  /* Copy point and parameter pair descriptions */

	  memcopy(spnt1,spnt2,ksize,DOUBLE);
	  memcopy(spar1,spar2,2,DOUBLE);

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
	      /*             Closed curve found */
	      goto finished;
	    }


	  /*         End while loop */
	}
      /* End inside */
      /* INSIDE TEST REMOVED BECAUSE OF CHANGED STRATEGY
	 }
      */
    nextdir:;
      /* End two step directions */
    }

 finished:

  /* In certain cases too many marched point may be found. These cases are:

  - Open curve and start of marching first guide point
  - Open curve and start of marching last guide point
  - Closed curve and this found in second marching direction

  In these cases some of the found points have to be discarded */

  scorpnt = s3dinf;
  scorpar = sp1inf;

  if (kstpch !=3 && kpoint>1)
    {

      /*  Open curve */

      if ( (kstart==1 && kgd1 == kpoint) ||
	   (kstart==kpoint && kgd1==1)      )
	{
	  /*  First marching direction traced curve */

	  knbinf = knb1;
	}
      else if ( (kstart==1 && kguide==kpoint) ||
		(kstart==kpoint && kguide==1)    )
	{
	  /* Second marching direction traced curve */

	  scorpnt = scorpnt + 10*(knb1-1);
	  scorpar = scorpar +  7*(knb1-1);
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
	  scorpar = scorpar +  7*(knb1-1);
	  knbinf  = knbinf - knb1 + 1;
	}
    }

  /* A curve is traced out only if at least two points are found, if less
     points found try to pick out constant parameter line */

 interpolate:

  if (knbinf>1)
    {
      if (igraph == 1 && knbinf > 1)
	{
	  /*  Output curve through s6line and s6move */

	  s6move(scorpnt);
	  for (ki=1,sp1=scorpnt+10;ki<knbinf;ki++,sp1+=10)
	    s6line(sp1);
	}

      if (icur > 0 && knbinf >1)
	{

	  /*  Make 3-D representation of intersection curve */

	  kpar = 0;

	  /*  We allocate space for parametrization array */


	  spar = newarray(knbinf,DOUBLE);
	  if (spar == SISL_NULL) goto err101;

	  s1359(scorpnt,aepsge,kdim,knbinf,kpar,spar,&q3dcur,&kstat);
	  if (kstat < 0) goto error;

	  /*  Set pointer in intcurve object to 3-D curve */
	  pintcr -> pgeom = q3dcur;

	  if (icur == 2)
	    {
	      /* Make curve in parameter plane */

	      kdim = 2;
	      kpar = 1;
	      s1359(scorpar,aepsge,kdim,knbinf,kpar,spar,&qp1cur,&kstat);
	      if (kstat < 0) goto error;


	      /* Set pointersin intcurve object to curves in parameter plane */
	      if (kpar1 ==2)
		{
		  pintcr -> ppar1 = qp1cur;
		}
	      else
		{
		  pintcr -> ppar2 = qp1cur;
		}
	    }
	}
    }
  /*  Dont use s9constline for the silhouette curves -- it won't work! */
  else if(pintcr->ipoint > 1 && ideg < 1003)
    {

      /* If no points produced on intersection curve */

      s1313_s9constline(ps1,eimpli,ideg,aepsge,pintcr,
			icur,igraph,&kstat);
      if (kstat<0) goto error;
      if (kstat==0) goto err185;
    }
  else
    goto err185;

  if (kdiv==1) goto war03;
  *jstat = 0;

  goto out;

  /* Iteration can not continue */
 war03:  *jstat = 3;
  goto out;

  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1313",*jstat,kpos);
  goto out;

  /* Error in surface description parameter direction does not exist */
 err123: *jstat = -123;
  s6err("s1313",*jstat,kpos);
  goto out;


  /* Error - SISL_NULL pointer was given */
  err150 :
    *jstat = -150;
  s6err("s1313",*jstat,kpos);
  goto out;


  /* Only degenerate or singular guide points */
 err185: *jstat = -185;
  goto out;


  /* Error in lower leve function */
 error:
  *jstat = kstat;
  s6err("s1313",*jstat,kpos);
  goto out;

 out:

  /* Free allocated space */

  if (sgpar  != SISL_NULL) freearray(sgpar);
  if (sgd1   != SISL_NULL) freearray(sgd1);
  if (s3dinf != SISL_NULL) freearray(s3dinf);
  if (sp1inf != SISL_NULL) freearray(sp1inf);
  if (spar   != SISL_NULL) freearray(spar);


  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
s1313_s9constline(SISLSurf *ps1,double eimpli[],int ideg,
		  double aepsge,SISLIntcurve *pintcr,int icur,
		  int igraph,int *jstat)
#else
     static void s1313_s9constline(ps1,eimpli,ideg,aepsge,pintcr,
				   icur,igraph,jstat)
     SISLSurf     *ps1;
     double   eimpli[];
     int      ideg;
     double   aepsge;
     SISLIntcurve *pintcr;
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
*              DO NOT USE IT FOR SILHOUETTE CURVES -- IT WON'T WORK!
*
*
* INPUT      : ps1    - Pointer to surface.
*              eimpli - Description of the implicit surface
*              ideg   - The degree of the implicit surface
*                        ideg=1: Plane
*                        ideg=2; Quadric surface
*                        ideg=1001: Torus surface
*                        ideg=1003: Silhouette line parallel projection
*                        ideg=1004: Silhouette line perspective projection
*                        ideg=1005: Silhouette line circular projection
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*              amax   - Maximal allowed step length. Not used.
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
* Revised by : Mike Floater, SI, 1991-01
*                Tried to add perspective and circular silhouettes (ideg=1004,ideg=1005)
*                but more work is needed. After improving s1313, s1309, and s1331
*                for all silhouette
*                types -- ideg=1003,1004,1005 this routine is no longer
*                useable even for ideg=1003. The problem is that there is not
*                enough information in qc2 and its first derivative for the
*                calls to s1331 and s1309. Perhaps the solution is to write a
*                new routine specially for silhouette curve intersections.
*                Until this problem is fixed S1313_S9CONSTLINE MUST NOT BE
*                CALLED WHEN ideg=1003, ideg=1004, ideg=1005 (see s1313).
*
*********************************************************************
*/
{
  int ki,kj,kl;            /* Control variables in for loops            */
  int kk,kn,kk1,kn1,kk2,kn2;/* Orders and numbers of knots               */
  int kpoint;              /* Number of points in guide curve           */
  int kleft = 0;           /* Pointer into knot vector                  */
  int kpar1;               /* Number of parameter direction in 1st. obj */
  int kpar2;               /* Number of parameter direction in 2st. obj */
  int ktype;               /* Type of intersection curve                */
  int kpos=0;              /* Position of error                         */
  int kstat;               /* Status variable returned form routine     */
  int kdir=0;              /* constant parameter line direction         */
  int knbpnt;              /* Number of points on constant parameter line */
  double tmax1,tmin1;      /* Minimum and maximum of first comp of guide points */
  double tmax2,tmin2;      /* Minimum and maximum of first comp of guide points */
  double tmax;             /* Maximum 3-D SISLbox side                      */
  double tdist,tang;       /* Distance and angle error                  */
  double *st,*st1,*st2;    /* Pointers to knot vectors                  */
  double *spoint;          /* Pointer to points on constant parameter line */
  double *sp1;             /* Pointer into array                        */
  double tsize1,tsize2;    /* Length of knot intervals                  */
  double sval1[2];         /* Limits of parameter plane in first SISLdir    */
  double sval2[2];         /* Limits of parameter plane in second SISLdir   */
  double sder[6];          /* SISLPoint and derivative on curve             */
  double sider[3];         /* SISLPoint on implicit surface                 */
  double snorm[3];         /* Normal on implicit surface                */
  double tsumold,tsum,tval;/* Parameter values                    */
  double *sgpar1=SISL_NULL;     /* Parameter pairs of guide point in surf 1  */
  double *sgpar2=SISL_NULL;     /* Parameter pairs of guide point in surf 2  */
  SISLCurve *qc1=SISL_NULL;         /* Pointer to 3-D curve                     */
  SISLCurve *qc2=SISL_NULL;         /* Pointer to 3-D curve                     */

  SISLCurve *qp1cur=SISL_NULL;      /* Pointer to curve in first parameter plane*/



  /* Make maximal step length based on box-size of surface */

  sh1992su(ps1,0,aepsge,&kstat);
  if (kstat < 0) goto error;

  tmax = MAX(ps1->pbox->e2max[0][0] - ps1->pbox->e2min[0][0],
	     ps1->pbox->e2max[0][1] - ps1->pbox->e2min[0][1]);
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


  /* Test that first object has 2 parameter
     direction and second object 0 */

  if (kpar1 == 2 && kpar2 == 0)
    {
      /*  Everithing is ok */
      ;
    }
  else if (kpar1 == 0 && kpar2 == 2)
    {
      sgpar1 = sgpar2;
    }
  else
    {
      goto err123;
    }


  /* Run through the parameter pairs to decide if a constant parameter line
     is possible */

  tmax1 = tmin1 = sgpar1[0];
  tmax2 = tmin2 = sgpar1[1];

  for (ki=1,kj=2,kl=3 ; ki < kpoint ; ki++,kj+=2,kl+=2)
    {
      tmin1 = MIN(tmin1,sgpar1[kj]);
      tmax1 = MAX(tmax1,sgpar1[kj]);
      tmin2 = MIN(tmin2,sgpar1[kl]);
      tmax2 = MAX(tmax2,sgpar1[kl]);
    }

  tsize1 = st1[kn1] - st1[kk1-1];
  tsize2 = st2[kn2] - st2[kk2-1];

  /* Check if constant parameter value within tolerance */

  if (DEQUAL((tmin1+tsize1),(tmax1+tsize1)) )
    {
      /*  Intersection possible constant parameter line with first parameter
	  constant value constant.

	  1. Pick out curve from surface
	  2. Pick out relevant part of curve */

      kdir = 1;

      s1437(ps1,(tmin1+tmax1)/2.0,&qc1,&kstat);
      if (kstat < 0) goto error;

      s1712(qc1,tmin2,tmax2,&qc2,&kstat);
      if (kstat < 0) goto error;
    }
  else if (DEQUAL((tmin2+tsize2),(tmax2+tsize2)) )
    {
      /*  Intersection possible constant parameter line with first parameter
	  constant value constant.

	  1. Pick out curve from surface
	  2. Pick out relevant part of curve */

      kdir = 2;

      s1436(ps1,(tmin2+tmax2)/2.0,&qc1,&kstat);
      if (kstat < 0) goto error;

      s1712(qc1,tmin1,tmax1,&qc2,&kstat);
      if (kstat < 0) goto error;
    }
  else
    goto war00;

  st = qc2 -> et;
  kk = qc2 -> ik;
  kn = qc2 -> in;

  /* Run through 2*kn points of the curve and check that they lie in the
     implicit surface by calculating the 2*kn points. */

  tsumold = st[kk-1];

  for (ki=0 ; ki <kn ; ki++)
    {
      if (kk>1)
	{
	  /* Make parameter value to use for calculation of curve point */

	  for (kl=1,kj=ki+1,tsum=(double)0.0 ; kl<kk ; kl++)
	    tsum += st[kj++];

	  tsum = tsum/(kk-1);
	}
      else
	tsum = st[ki];

      tval = tsum;

      for (kj=0 ; kj<2 ; kj++)
	{
	  /* Calculate point on curve */

	  s1221(qc2,1,tval,&kleft,sder,&kstat);
	  if (kstat < 0) goto error;

	  /* Calculate normal to implicit surface */

	  s1331(sder,eimpli,ideg,-1,sider,snorm,&kstat);
	  if (kstat < 0) goto error;

	  /* Project point onto implicit surface */

	  tdist = fabs(s1309(sder,snorm,eimpli,ideg,&kstat));
	  if (kstat<0) goto error;
	  if (kstat==2) goto war00;

	  /* Both the position of the two points should be within the relative
	     computer resolution for the point to be accepted. Correspondingly the
	     direction of the intersection curve and the constant parameter line
	     should be within the computer resolution to be accepted. */

	  if (DNEQUAL(tdist+tmax,tmax))
	    goto war00;

	  /* Distance within tolerance, check that the angle between surface
	     normal and curve tangent is PIHALF, if both these vectors have a
	     nonzero length. */

	  if (s6length(snorm,3,&kstat) != (double)0.0 &&
	      s6length(sder+3,3,&kstat) != (double)0.0  )
	    {
	      tang = s6ang(snorm,sder+3,3);
	      if (DNEQUAL(fabs(tang),PIHALF) ) goto war00;
	    }
	  tval = (tsumold+tsum)/(double)2.0;
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
      /* Make curve in parameter plane */

      double svert[4],sknot[4];

      if (kdir==1)
	{
	  svert[0] = svert[2] = (tmin1+tmax1)/(double)2.0;
	  svert[1] = tmin2;
	  svert[3] = tmax2;
	  sknot[0] = sknot[1] = tmin2;
	  sknot[2] = sknot[3] = tmax2;
	}
      else
	{
	  svert[0] = tmin1;
	  svert[2] = tmax1;
	  svert[1] = svert[3] = (tmin2+tmax2)/(double)2.0;
	  sknot[0] = sknot[1] = tmin1;
	  sknot[2] = sknot[3] = tmax1;
	}
      qp1cur = newCurve(2,2,sknot,svert,1,2,1);
      if (qp1cur==SISL_NULL) goto err101;
      pintcr -> ppar1 = qp1cur;
    }

  /* war01: */
  *jstat = 1;
  goto out;

  /* Iteration can not continue */
 war00:  *jstat = 0;
  goto out;


  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1313",*jstat,kpos);
  goto out;

  /* Error in surface description parameter direction does not exist */
 err123: *jstat = -123;
  s6err("s1313_s9constline",*jstat,kpos);
  goto out;

  /* Error in lower leve function */
 error:
  *jstat = kstat;
  s6err("s1313_s9constline",*jstat,kpos);
  goto out;

 out:;
  if (qc1 != SISL_NULL) freeCurve(qc1);
  if (qc2 != SISL_NULL) freeCurve(qc2);
}
