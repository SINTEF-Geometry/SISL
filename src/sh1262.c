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

#define SH1262

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static void sh1262_s9blendcoef(double [],double [],double [],int,int,double *,
			       double *,int *);
static void sh1262_s9blendder(double [],double  [],int,int,int,double,double,
			      double,double,double [],double [],double [],
			      double [],int *);
static void sh1262_s9hermit(double [],int,double,double,int,int *);
#else
static void sh1262_s9blendcoef();
static void sh1262_s9blendder();
static void sh1262_s9hermit();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
      sh1262(SISLCurve *vcurve[],int iedge,int inmbx,double ecoef[],int *jstat)
#else	 
void sh1262(vcurve,iedge,inmbx,ecoef,jstat)
     int iedge,inmbx,*jstat;
     double ecoef[];
     SISLCurve *vcurve[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Find coefficients of cubic blending functions used to
*              blend two derivative curves along each edge of a vertex
*              region into a cross derivative curve pr edge.
*
*
*
* INPUT      : qcurve     - Array containing pointers to edge curves 
*                           of the surfaces around the vertex region.
*                           For each edge the position edge curve and
*                           the derivative curves in the two parameter
*                           directions of the surface, is given. The
*                           sequence is the following : Position curve
*                           along first edge, cross derivative curve along
*                           first edge, derivative curve along first
*                           edge, position curve along second edge, etc.
*                           Dimension of the array is 3*iedge.
*              iedge      - Number of edges of the vertex region.
*                           iedge = 3, 4, 5 or 6.
*              inmbx      - Number of coefficients of all the blending
*                           functions.
* 
*
* OUTPUT     : ecoef      - Coefficients of blending functions. There is
*                           two blending functions for each edge and the
*                           functions are represented as one-dimensional 
*                           cubic Bezier curves. Then the number of coefficients 
*                           of each function is 4. First all the blending
*                           functions corresponding to the first derivative
*                           curve are represented counter clockwise around the
*                           region, the all the functions corresponding to 
*                           the second derivative curve are stored. The dimension
*                           of the array is inmbx, i.e. 2*4*iedge.
*                           
*              jstat      - status messages  
*                                         = 2      : Requirements on tangent
*                                                    not satisfied completely.
*                                         = 1      : Requirements on twist
*                                                    not satisfied completely.
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
* CALLS      : s1221          - Evaluator of B-spline curve.   
*              s6crss         - Cross product of two vectors.  
*              s6scpr         -  Scalar product of two vectors. 
*              sh1262_s9blendcoef    - 
*              sh1262_s9blendder     -
*              sh1262_s9hermit       -
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
  int kstat = 0;     /* Status variable.     */
  int kstat1 = 0;    /* Status variable used to check if the requirements
			to input is satisfied.         */
  int ki,kj,kk;      /* Counters.            */
  int ki4,ki12;      /* ki4 = 4*ki, and ki12 = 12*ki.  */
  int kj4,kj12;      /* kj4 = 4*kj, and kj12 = 12*kj.  */
  int kdim = 3;      /* Dimension of geometry space.   */
  int kleft = 0;     /* Parameter used in curve evaluation.  */
  int kder = 1;      /* Number of derivatives of curve to evaluate. */
  int kncurve = 3*iedge;  /* Number of input curves.   */
  int kncond = 4*iedge;   /* Number of coeffecients of blending
                             functions corresponding to derivative
                             curves in one parameter direction.    */
  int ksign;              /* Indicates if the tangent along a
                             position curve and the cross tangent
                             of the adjacent edge in a corner is
                             to have the same sign.                */
  int krang;              /* Rang of equation system used to find derivatives
                             in endpoints  of blending functions. */
  double ang_tol = (double)0.00001;   /* Tolerance used in computing the
					 space spanned by tang. vectors. */
  double tfac;            /* Used in linearity test of blending function. */
  double tref=(double)0.0; /* Refereance value in equality test.  */
  double tpar1,tpar2;     /* Endpoints of parameter interval of curve. */
  double thelp1,thelp2,thelp3,thelp4;   /* Help variables in computation 
					   of rang.    */
  double tnorm1,tnorm2;   /* Lenght of normals.        */
  double *spar1 = SISL_NULL;   /* Startpoints of parameter intervals of 
                             input curves.       */
  double *spar2 = SISL_NULL;   /* Endpoints of parameter intervals of 
                             input curves.       */
  double *sder = SISL_NULL;    /* Result of curve evaluation. The values are
                             stored as follows : Value and first derivative
                             of all curves corresponding to an edge in
                             first endpoint, the same values in the second
                             endpoint. All curves of one edge is treated
                             before the curve of the next edge is treated. */
  double sa[12];          /* Matrix of equation system used to find
                             derivatives of blending functions.     */
  double sb[3];           /* Right side of equation system.         */
  double snorm1[3];       /* Normal of tangent plane at a corner of
                             the vertex region.                     */
  double snorm2[3];       /* Normal of tangent plane at a corner of
                             the vertex region.                     */
  SISLCurve *qc;              /* Pointer to curve.                      */
  
  /* Test input.  */

  for (ki=0; ki<iedge; ki++)
    if (vcurve[ki]->idim != kdim) goto err104;
  
  /* Allocate space for array containing results of curve evaluation.  */
  
  if ((sder = newarray(12*kdim*iedge,DOUBLE)) == SISL_NULL) goto err101;
  if ((spar1 = newarray(kncurve,DOUBLE)) == SISL_NULL) goto err101;
  if ((spar2 = newarray(kncurve,DOUBLE)) == SISL_NULL) goto err101;  
  
  /* Evaluate boundary curves in the  endpoints.  */

  for (ki=0; ki<kncurve; ki++)
    {
      qc=vcurve[ki];
      
      /* Fetch parameter values at endpoints.  */

      spar1[ki] = tpar1 = *(qc->et + qc->ik - 1);
      spar2[ki] = tpar2 = *(qc->et + qc->in);
      
      kj = ki/3;
      kk = ki % 3;
      
      /* Evaluate curve in startpoint.  */

      s1221(qc,kder,tpar1,&kleft,sder+2*(6*kj+kk)*kdim,&kstat);
      if (kstat < 0) goto error;
      
      /* Evaluate curve in endpoint.  */

      s1221(qc,kder,tpar2,&kleft,sder+2*(6*kj+kk+3)*kdim,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* Compute first and last coefficients of blending functions, i.e find
     value of blending functions in the endpoints.  */

  for (ki=0; ki<iedge; ki++)
    {
      /* For each edge, make sure that the endpoints of the cross tangent 
	 curve is equal (possibly expect from a sign) to the tangents of the
         position curve along the adjacent edges.  */

      kj = (ki > 0) ? ki-1 : iedge-1;
      ki4 = 4*ki, ki12 = 12*ki;
      kj4 = 4*kj, kj12 = 12*kj;
      
      /* Treat startpoint of edge.  */

      ksign = -1;
      sh1262_s9blendcoef(sder+(ki12+2)*kdim,sder+(ki12+4)*kdim,sder+(kj12+7)*kdim,
		  kdim,ksign,ecoef+ki4,ecoef+kncond+ki4,&kstat);
      if (kstat < 0) goto error;
      kstat1 = MAX(kstat,kstat1);
      
      /* Treat endpoint of edge.  */

      ksign = 1;
      sh1262_s9blendcoef(sder+(ki12+8)*kdim,sder+(ki12+10)*kdim,
		  sder+(((ki+1)%iedge)*12+1)*kdim,
		  kdim,ksign,ecoef+ki4+3,ecoef+kncond+ki4+3,&kstat);
      if (kstat < 0) goto error;
      kstat1 = MAX(kstat,kstat1);
    }
  
  for (tref=(double)0.0, ki=0; ki<iedge; ki++)
    {
      /* For each edge, set up an equation system to find the derivatives
	 of the blending functions excisting at the corner lying at the 
	 startpoint of the edge (also the blending functions corresponding
	 to the previous edge). The equation system represent the conditions
	 to be satisfied to have a consistent twist vector in this corner. */

      kj = (ki > 0) ? ki-1 : iedge-1;
      ki4 = 4*ki, ki12 = 12*ki;
      kj4 = 4*kj, kj12 = 12*kj;
      
      for (kk=0; kk < kdim; kk++)
	{
          sa[kk*4] = sder[(ki12+2)*kdim+kk];
	  sa[kk*4+1] = sder[(ki12+4)*kdim+kk];
	  sa[kk*4+2] = sder[(kj12+8)*kdim+kk];
	  sa[kk*4+3] = sder[(kj12+10)*kdim+kk];

	  sb[kk] = - ecoef[ki4]*sder[(ki12+3)*kdim+kk]
	    - ecoef[kncond+ki4]*sder[(ki12+5)*kdim+kk]
	    - ecoef[kj4+3]*sder[(kj12+9)*kdim+kk]
		- ecoef[kncond+kj4+3]*sder[(kj12+11)*kdim+kk];
	  
	  tref = MAX(tref,MAX(sa[kk*4],sa[kk*4+1]));
	  tref = MAX(tref,MAX(sa[kk*4+2],sa[kk*4+3]));
	}

      /* Find rang of equation system, i.e. check wether the 4
	 tangent vectors are able to span a plane.   */

      s6crss(sder+(ki12+2)*kdim,sder+(ki12+4)*kdim,snorm1);
      tnorm1 = s6length(snorm1,kdim,&kstat);
      
      s6crss(sder+(kj12+8)*kdim,sder+(kj12+10)*kdim,snorm2);
      tnorm2 = s6length(snorm2,kdim,&kstat);
      
      thelp1 = s6ang(snorm1,sder+(kj12+8)*kdim,kdim);
      thelp2 = s6ang(snorm1,sder+(kj12+10)*kdim,kdim);

      thelp3 = s6ang(snorm2,sder+(ki12+2)*kdim,kdim);
      thelp4 = s6ang(snorm2,sder+(ki12+4)*kdim,kdim);

      if ((DNEQUAL(tnorm1+tref,tref) && 
	   fabs(thelp1-PIHALF) < ang_tol &&
	   fabs(thelp2-PIHALF) < ang_tol) ||
	  (DNEQUAL(tnorm2+tref,tref) &&
	   fabs(thelp3-PIHALF) < ang_tol &&
	   fabs(thelp4-PIHALF) < ang_tol) ||
	  (DEQUAL(tnorm1+tref,tref) &&
	   (fabs(thelp1-PIHALF) < ang_tol ||
 	    fabs(thelp2-PIHALF) < ang_tol)) ||
	  (DEQUAL(tnorm2+tref,tref) &&
	   (fabs(thelp3-PIHALF) < ang_tol ||
 	    fabs(thelp4-PIHALF) < ang_tol)) ||
	  (DEQUAL(tnorm1+tref,tref) && DEQUAL(tnorm2+tref,tref)))
	krang = 2;
      else
	krang = 3;
      
      /* Estimate derivatives in the current endpoints of the blending 
	 functions existing in this corner.  */

	sh1262_s9blendder(sa,sb,4,kdim,krang,spar1[3*ki],spar2[3*ki],spar1[3*kj],
		 spar2[3*kj],ecoef+ki4,ecoef+kncond+ki4,ecoef+kj4,
		 ecoef+kncond+kj4,&kstat);
      if (kstat < 0) goto error;
      kstat1 = MAX(kstat,kstat1);
    }
   
  for (ki=0; ki<iedge; ki++)
    {
      /* Interpolate blending curves.  */

       sh1262_s9hermit(ecoef+4*ki,4,spar1[3*ki],spar2[3*ki],1,&kstat);
      if (kstat < 0) goto error;

      sh1262_s9hermit(ecoef+kncond+4*ki,4,spar1[3*ki],spar2[3*ki],1,&kstat);
      if (kstat < 0) goto error;
      
      /* Test blending curves. If the coeffecients deviate much from
	 those of a linear curve, they are linearized. In addition, the
	 coefficients of the blending curve corresponding to the input
	 cross derivative curve is not allowed to be negative.      */
      
      tfac = MAX(MAX(fabs(ecoef[4*ki+3]),fabs(ecoef[4*ki])),
		  MAX(fabs(ecoef[kncond+4*ki+3]),fabs(ecoef[kncond+4*ki])));
      for (kj=1; kj<3; kj++)
      {
	 if (ecoef[4*ki+kj] < MIN(ecoef[4*ki+3],ecoef[4*ki])-tfac ||
	     ecoef[4*ki+kj] > MAX(ecoef[4*ki+3],ecoef[4*ki])+tfac ||
	     ecoef[4*ki+kj] < DZERO)
	 {
	    ecoef[4*ki+kj] = ((double)(3-kj)*ecoef[4*ki] +
	       (double)kj*ecoef[4*ki+3])/(double)3.0;
	    kstat1 = MAX(kstat1,1);
	 }
	 
	 if (ecoef[kncond+4*ki+kj] < 
	     MIN(ecoef[kncond+4*ki+3],ecoef[kncond+4*ki])-tfac ||
	     ecoef[kncond+4*ki+kj] > 
	     MAX(ecoef[kncond+4*ki+3],ecoef[kncond+4*ki])+tfac)
	 {
	    ecoef[kncond+4*ki+kj] = ((double)(3-kj)*ecoef[kncond+4*ki] +
	       (double)kj*ecoef[kncond+4*ki+3])/(double)3.0;
	    kstat1 = MAX(kstat1,1);
	 }
      }
    }
  
  
  /* Coefficients of blending functions found.  */

  *jstat = kstat1;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  /* Error in input. Dimension not equal to 3.  */

  err104 :
    *jstat = -104;
  goto out;
  
  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  goto out;
  
  out :
    
    /* Free space occupied by local arrays.  */

    if (sder != SISL_NULL) freearray(sder);
  if (spar1 != SISL_NULL) freearray(spar1);
  if (spar2 != SISL_NULL) freearray(spar2);

  return;
}

					
#if defined(SISLNEEDPROTOTYPES)
   static
      void
	    sh1262_s9blendcoef(double evecu[],double evecv[],double etang[],
			       int idim,int isign,double *coef1,
			       double *coef2,int *jstat)
#else
static void sh1262_s9blendcoef(evecu,evecv,etang,idim,isign,coef1,coef2,jstat)
     int idim,isign,*jstat;
     double evecu[],evecv[],etang[],*coef1,*coef2;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given three vectors, evecu, evecv and etang, find the 
*              coefficients, coef1 and coef2, such that the vector
*              coef1*evecu + coef2*evecv, is as close as possible
*              to the vector etang. If the three vectors lie in a
*              plane, coef1*evecu + coef2*evecv = etang.
*              
*
* INPUT      : evecu      - First vector.
*              evecv      - Second vector.
*              etang      - Vector to approximate.
*              idim       - Dimension of geometry space.
*              isign      - Sign with wich etang is to be multiplied.
*
*
* OUTPUT     : coef1      - First factor.
*              coef2      - Second factor.
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Minimize the square of the expression 
*                   dist(coef1*evecu+coef2*evecv,isign*etang)
*              over coef1 and coef2.
*              The expression is differentiated and set equal to
*              zero. Then this equation system of 2 equations
*              with two unknowns is solved.
*
*********************************************************************
*/
{
  
  int kstat = 0;           /* Status variable.                    */
  int ki;                  /* Counter.                            */
  double tdotuu;           /* Scalar product of evecu and evecu.  */
  double tdotuv;           /* Scalar product of evecu and evecv.  */
  double tdotutang;        /* Scalar product of evecu and etang.  */
  double tdotvv;           /* Scalar product of evecv and evecv.  */
  double tdotvtang;        /* Scalar product of evecv and etang.  */
  double tdiv;             /* Determinant of equation system.     */
  
  /* Set output to zero. */

  *coef1 = (double)0.0;
  *coef2 = (double)0.0;
  
  /* Compute coefficients of equation system.  */

  tdotuu = s6scpr(evecu,evecu,idim);
  tdotuv = s6scpr(evecu,evecv,idim);
  tdotutang = (double)isign*s6scpr(evecu,etang,idim);
  tdotvv = s6scpr(evecv,evecv,idim);
  tdotvtang = (double)isign*s6scpr(evecv,etang,idim);

  tdiv = tdotuv*tdotuv - tdotuu*tdotvv;
  if (DEQUAL(tdiv,DZERO))
    {
      if (DEQUAL(tdotuu,DZERO) && DEQUAL(tdotvv,DZERO));
      else if (DEQUAL(tdotuu,DZERO))
	  *coef2 = s6length(etang,idim,&kstat)/sqrt(tdotvv);
      else
	*coef1 = s6length(etang,idim,&kstat)/sqrt(tdotuu);
      goto out;
    }
  
  /* Compute output factors.  */

  *coef1 = (tdotvtang*tdotuv - tdotutang*tdotvv)/tdiv;
  *coef2 = (tdotutang*tdotuv - tdotvtang*tdotuu)/tdiv;

  /* Test result.  */

  for (ki=0; ki<idim; ki++)
    if (fabs((*coef1)*evecu[ki]+(*coef2)*evecv[ki]-isign*etang[ki]) > REL_PAR_RES) 
      break;
  
  if (ki < idim) goto warn1;

  *jstat = 0;
  goto out;
  
  /* Equality not achieved.  */

  warn1 :
    *jstat = 2;
  goto out;
  
  out :
    return;
}


#if defined(SISLNEEDPROTOTYPES)
   static
      void
	    sh1262_s9hermit(double econd[],int icond,double astart,
			    double aend,int idim,int *jstat)
#else
static void sh1262_s9hermit(econd,icond,astart,aend,idim,jstat)
     int icond,idim,*jstat;
     double econd[],astart,aend;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Hermite interpolation of position derivative in the
*              two endpoints represented as a Bezier curve on the interval 
*              [astart,aend]. icond is expected to be 4.
*              
*
*
*
* INPUT      : icond      - Number of interpolation conditions. 
*                           icond = 4.
*              astart     - Start of parameter interval.
*              aend       - End of parameter interval
*              idim       - Dimension of geometry space.
*
*
* INPUT/OUTPUT : econd    - Interpolation conditions as input, Bezier coefficients
*                           as output. The dimension is icond*idim.
*                       
*
* OUTPUT     : jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
*********************************************************************
*/
{
  int ki;    /* Index.  */

  /* Test input. The number of conditions has to be 4.  */

  if (icond != 4) goto err001;
  
  /* Hermit interpolation with Bezier curve of order 4.  */

  for (ki=0; ki<idim; ki++)
    {
      econd[idim+ki] = ONE_THIRD*(aend-astart)*econd[idim+ki] + econd[ki];
      econd[2*idim+ki] = -ONE_THIRD*(aend-astart)*econd[2*idim+ki] + econd[3*idim+ki];
    }
  
  *jstat = 0;
  goto out;
  
  /* Error in input. Number of coefficients not 4.  */

  err001 :
    *jstat = -1;
  goto out;
  
  out :
    return;
}


#if defined(SISLNEEDPROTOTYPES)
   static
      void
	    sh1262_s9blendder(double ea[],double eb[],int ix,int ieq,
			      int irang,double astarti,double aendi,
			      double astartj,double aendj,double ealfai[],
			      double ebetai[],double ealfaj[],
			      double ebetaj[],int *jstat)
#else
static void sh1262_s9blendder(ea,eb,ix,ieq,irang,astarti,aendi,astartj,aendj,ealfai,
		       ebetai,ealfaj,ebetaj,jstat)
     int ix,ieq,irang,*jstat;
     double ea[],eb[],astarti,aendi,astartj,aendj,ealfai[],ebetai[],
       ealfaj[],ebetaj[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Estimate the derivatives in one endpoint of the blending 
*              functions meeting in a corner of a vertex region.
*
*
* INPUT      : ea         - Matrix containg coefficients in equation
*                           system representing conditions on the 
*                           derivatives.
*              eb         - Rigth side of equation system.
*              ix         - Number of coefficients. ix = 4.
*              ieq        - Number of equations. ieq = 3.
*              irang      - Rang of equation system.
*
*
* INPUT/OUTPUT : ealfai   - Position and derivatives or first blending 
*                           function along first adjacent edge.
*                           Estimate derivative of first endpoint.
*                ebetai   - Position and derivatives or second blending 
*                           function along first adjacent edge.
*                           Estimate derivative of first endpoint.
*                ealfaj   - Position and derivatives or first blending 
*                           function along second adjacent edge.
*                           Estimate derivative of second endpoint.
*                ebetaj   - Position and derivatives or second blending 
*                           function along second adjacent edge.
*                           Estimate derivative of second endpoint.
*                
* 
* OUTPUT     : jstat      - status messages  
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
* CALLS      : s6lufacp - LU-factorizing of matrix in equation system. 
*              s6lusolp - Solve equtaion system. Matrix is LU-factorized. 
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
  int kstat = 0;          /* Status variable.  */
  int ki,kj,kk,kh;        /* Counters.         */
  int l[2][2];            /* Pointer to lines of equation system to be used
                             when solving with respect to two unknowns.   */
  int lpiv[3];            /* Pivoting array. Used in LU-factorizing.      */
  int kcount;             /* Number of iterations used to find derivatives. */
  double teps = 0.001;    /* Local tolerance.                             */
  double tdiff;           /* Difference between derivatives of linear blending
                             function and current blending function.      */
  double tdum;            /* Used in testing of result.                   */
  double tdet;            /* Current determinant of equation system.      */
  double sdet[2];         /* Determinant of equation system.              */
  double sx[4];           /* Derivative of linear blending funtion.       */
  double sxx[4];          /* Derivative of current cubic blending funtion. */
  double sb[3];           /* Rigth side of local equation system.         */
  double smat[9];         /* Matrix of local eqution system.              */
  
  /* Test input.  */

  if (ix != 4 || ieq != 3) goto err001;
  
  /* Set up derivatives of linear blend.  */

  sx[0] = (ealfai[3] - ealfai[0])/(aendi - astarti);
  sx[1] = (ebetai[3] - ebetai[0])/(aendi - astarti);
  sx[2] = (ealfaj[3] - ealfaj[0])/(aendj - astartj);
  sx[3] = (ebetaj[3] - ebetaj[0])/(aendj - astartj);  
  
  if (irang == 2)
    {
      /* Rang of equation system equal to 2. */
      /* Find equations to be used, i.e. find the lines in the
	 equation system which gives the most stable 2x2 equation
	 system when finding the 2 last derivatives when the 2
	 first are set.  */

      sdet[0] = sdet[1] = (double)0.0;
      
      for (ki=0; ki<ieq; ki++)
	for (kj=ki+1; kj<ieq; kj++)
	  for (kk=0,kh=0; kh<2; kk+=2,kh++)
	    {
	      tdet = ea[4*ki+kk]*ea[4*kj+kk+1] - ea[4*ki+kk+1]*ea[4*kj+kk];
	      if (fabs(tdet) > fabs(sdet[kh]))
		{
		  sdet[kh] = tdet;
		  l[kh][0] = ki, l[kh][1] = kj;
		}
	    }
      for (kh=0; kh<2; kh++)
	if (DEQUAL(sdet[kh],DZERO)) sdet[kh] = (double)1.0;
      
      kcount = 0;
      
      /* Set the two first derivatives equal to those of a corresponding
	 linear blending function.   */

      ki = 0, kk = 2;
      sxx[ki] = sx[ki];
      sxx[ki+1] = sx[ki+1];
      
      while (1)
	{
	  /* Iterate until the conditions on the derivatives is satisfied
	     and the difference from the linear case is eqal for all 
	     derivatives or the number of iterations is equal to 4.   */

	  /* Find right side of local equation system for finding the 
	     remaining two derivatives.   */

	  for (kj=0; kj<ieq; kj++)
	    sb[kj] = eb[kj] - ea[4*kj+ki]*sxx[ki] - ea[4*kj+ki+1]*sxx[ki+1];
	  
	  /* Compute the remaining two derivatives.  */

	  kj = kk/2;
	  sxx[kk] = (-ea[4*l[kj][0]+kk+1]*sb[l[kj][1]] 
		       + ea[4*l[kj][1]+kk+1]*sb[l[kj][0]])/sdet[kj];
	  sxx[kk+1] = (ea[4*l[kj][0]+kk]*sb[l[kj][1]] 
		       - ea[4*l[kj][1]+kk]*sb[l[kj][0]])/sdet[kj];
	  
	  /* Check if the difference between current derivatives and 
	     derivatives of a linear blend is small.  */

	  tdiff = fabs(sxx[0] - sx[0]);
	  for (kh=1; kh<4; kh++)
	    if (fabs(tdiff-fabs(sxx[kh]-sx[kh])) > teps) break;
	  
	  if (kh == 4 || kcount == 4) break;   /* Stop iteration.  */
	  
	  /* Initiate next iteration step.  */

	  sxx[kk] = (sx[kk] + sxx[kk])/(double)2.0;
	  sxx[kk+1] = (sx[kk+1] + sxx[kk+1])/(double)2.0;
	  
	  ki = (ki + 2) % 4;
	  kk = (kk + 2) % 4;

	  kcount++;
	}
    }
  else if (irang == 3)
    {
      /* Rang of equation system equal to 3. */
      ki = 0;
      kcount = 0;
      
      while (1)
	{
	  /* Iterate until the conditions on the derivatives is satisfied
	     and the difference from the linear case is eqal for all 
	     derivatives or the number of iterations is equal to 4.   */

	  /* Set the current derivative equal to those of a corresponding
	     linear blending function.   */

	  sxx[ki] = sx[ki];
	  
	  /* Set up local equation system.  */

	  for (kj=0; kj<ieq; kj++)
	    {
	      sb[kj] = eb[kj] - ea[kj*4+ki]*sxx[ki];
	      for (kk=0,kh=0; kk<ix; kk++)
		{
		  if (kk == ki) continue;
		  smat[kj*(ix-1)+kh] = ea[kj*ix+kk];
		  kh++;
		}
	    }
	  
	  /* Compute the remaining three derivatives.  */

	  s6lufacp(smat,lpiv,ieq,&kstat);
	  if (kstat < 0) goto error;
	  
	  s6lusolp(smat,sb,lpiv,ieq,&kstat);
	  if (kstat < 0) goto error;
	  
	  for (kk=0,kh=0; kk<ix; kk++)
	    {
	      if (kk == ki) continue;
	      sxx[kk] = sb[kh++];
	    }
	  
	  /* Check if the difference between current derivatives and 
	     derivatives of a linear blend is small.  */

	  tdiff = fabs(sxx[0] - sx[0]);
	  for (kh=1; kh<4; kh++)
	    if (fabs(tdiff-fabs(sxx[kh]-sx[kh])) > teps) break;
	  
	  if (kh == 4 || kcount == 4) break;   /* Stop iteration.  */
	  
	  /* Initiate next iteration step.  */

	  kj = (ki + 1) % 4;
	  sxx[kj] = (sx[kj] + sxx[kj])/(double)2.0;
	  
	  ki = kj;
	  kcount++;
	}
      
    }
  
  /* Copy derivatives into output array.  */ 

  ealfai[1] = sxx[0];
  ebetai[1] = sxx[1];
  ealfaj[2] = sxx[2];
  ebetaj[2] = sxx[3];
  
  /* Test result.  */

  for (ki=0; ki<ieq; ki++)
    {
      tdum = (double)0.0;
      for (kj=0; kj<ix; kj++)
	tdum += ea[ki*ix+kj]*sxx[kj];
      if (fabs(tdum - eb[ki]) > REL_PAR_RES) break;
    }
  if (ki < ieq) goto warn1;
  
  *jstat = 0;
  goto out;
  
  /* Twist requirement not satisfied completely.  */

  warn1 :
    *jstat = 1;
  goto out;
  
  /* Error in input.  */

  err001 :
    *jstat = -1;
  goto out;
  
  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  goto out;
  
  out :
    return;
}
