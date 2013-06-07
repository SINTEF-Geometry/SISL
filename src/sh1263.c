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

#define SH1263

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static void sh1263_s9makec0(SISLCurve **,int,int *);
static void sh1263_s9checkcrtan(SISLCurve *,int *); 
#else
static void sh1263_s9makec0();
static void sh1263_s9checkcrtan(); 
#endif

#if defined(SISLNEEDPROTOTYPES)
void
      sh1263(SISLCurve *vcurve[],int iedge,SISLCurve *vboundc[],int *jstat)
#else
void sh1263(vcurve,iedge,vboundc,jstat)
   SISLCurve *vcurve[];
     int iedge;
   SISLCurve *vboundc[];
     int *jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Prepare interpolation conditions for blending an
*              n-sided vertex region.
*
*
*
* INPUT      : vcurve     - Array containing pointers to edge curves 
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
*                       
*
* OUTPUT     : vboundc    - Boundary conditions suitable for vertex blend.
*                           For each edge a position curve and cross
*                           derivative curve is given. The sequence is
*                           position curve of first edge, cross derivative
*                           curve of first edge, position curve of second
*                           edge, etc. Dimension of the array is 2*iedge.
*              jstat      - status messages  
*                                         = 1      : Requirements on cross
*                                                    tangent curve not satisfied
*                                                    completely.
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : If necessary reparametrize the edge-curve in order to avoid
*              tangent lengths of the position curves which very long.
*              Blend the derivative curves to produce the cross tangent
*              curves in such a way that :
*              1. In a corner of the region, the cross tangent belonging to
*                 one of the adjacent edges is equal to the tangent belonging
*                 to the other adjacent edge (except for a sign).
*              2. The derivatives of the two cross tangent curves meeting in a
*                 corner are equal, i.e. consistent twist.
*              Check that the cross tangent curves will not imply intersecting
*              parameter curves of the blending surfaces.
*
*
* REFERENCES : 
*              
*
* USE        : 3D geometry only.
*
*-
* CALLS      : sh1260 - Reparametrize if tangents in endpoints
*                       of position curve are long.            
*              sh1261 - Compute cross tangent curve given 
*  		        blending functions.                    
*              sh1262 - Compute coefficients of cubic blending
*		        functions.                             
*              freeCurve - Free space occupied by curve object.   
*              newCurve  - Create new curve object.  
*              sh1263_s9checkcrtan() - Routine to be implemented. 
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
  int kstat = 0;        /* Status variable.  */
  int kstat1 = 0;       /* Status on prepare. */
  int ki;               /* Counter.          */
  int kcopy = 1;        /* Indicates that data arrays is to be copied
			   when creating a new curve.                 */
  int knmbx = 8*iedge;  /* Number of coefficients of blending functions. */

  double tonethird = (double)1.0/(double)3.0;  /* Constant used to check
						  if tangent of position
						  curve is too long.     */
  double *scoef = SISL_NULL;     /* Array containing coefficients fo blending
			       functions.                                */
  double (*s1)[4],(*s2)[4]; /* Pointers to coefficients of one blending
			       function.                                 */
  
  SISLCurve* (*qc)[3];          /* Pointer to curves corresponding to one 
			       edge of the vertex region.                */
  SISLCurve **qcurve = SISL_NULL;    /* Array of pointers to reparametrized curves.*/
  SISLCurve *qpt;               /* Pointer to current curve.                 */
  
  /* Allocate scratch for arrays used to store equation system.  */

  if ((scoef = newarray(knmbx,DOUBLE)) == SISL_NULL) goto err101;
  
  /* Copy input curves to local arrays. First allocate scratch for
     pointer array.  */

  if ((qcurve = newarray(3*iedge,SISLCurve*)) == SISL_NULL) goto err101;
  
  for (ki=0; ki<3*iedge; qpt++,ki++)
    {
      qpt=vcurve[ki];

      /* Copy curve.  */

      if ((qcurve[ki] = newCurve(qpt->in,qpt->ik,qpt->et,qpt->ecoef,qpt->ikind,
				 qpt->idim,kcopy)) == SISL_NULL) goto err101;
     } 

  /* Make sure that all curves are represented with k-tupple knots in
     the endpoints.  */
  
  s1349(3*iedge,qcurve,&kstat);
  if (kstat < 0) goto error;
		 
  /* Check if reparametrization is necessary, and in that case perform
     the reparametrization.  */

  for (ki=0, qc=(SISLCurve*(*)[3])qcurve; ki<iedge; ki++,qc++)
    {
       sh1260(tonethird,*qc,3,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* Find the coeffecients of the blending functions.   */

  sh1262(qcurve,iedge,knmbx,scoef,&kstat);
  if (kstat < 0) goto error;

  kstat1 = kstat;
  
  /* kstat may have the following non-negative values :
            0 - The criterions on the blending functions are satisfied.
	    1 - The criterion is not satisfied regarding twist.
	    2 - The tangent vectors on some corners are not consistent. */
  
  /* Move position curves to output array. Compute and check cross 
     tangent curves.   */

  for (ki=0,s1=(double(*)[4])scoef,s2=s1+iedge; ki<iedge; ki++,s1++,s2++)
    {
      /* Set pointer to position curve.  */

      vboundc[2*ki] = qcurve[3*ki];
      
      /* Compute cross tangent curve.  */

      sh1261(qcurve[3*ki+1],qcurve[3*ki+2],*s1,4,*s2,4,vboundc+2*ki+1,&kstat);
      if (kstat < 0) goto error;
      
      /* Check cross tangent curve. Change curve if necessary. */

      sh1263_s9checkcrtan(vboundc[2*ki+1],&kstat);
      if (kstat < 0) goto error;
    }
  
  if (kstat1 == 2)
  {
     /* The tangent conditions in some corner(s) of the vertex region
	are not satisfied. Thus, the blending will not be C0 in some
	areas. Make sure that the blend is C0, sacrifising the G1-
	condition.      */
     
     sh1263_s9makec0(vboundc,iedge,&kstat);
     if (kstat < 0) goto error;
  }

  /* Boundary conditions for vertex blend prepared.  */

  *jstat = kstat1;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  goto out;
  
  out :
    
    /* Free space occupied by local arrays and curves.  */

  if (scoef != SISL_NULL) freearray(scoef);
  if (qcurve != SISL_NULL)
    {
      for (ki=0; ki<3*iedge; ki+=3)
	{
	  if (qcurve[ki+1] != SISL_NULL) freeCurve(qcurve[ki+1]);
	  if (qcurve[ki+2] != SISL_NULL) freeCurve(qcurve[ki+2]);
	}
      freearray(qcurve);
    }
  
  return;
}

  



   
#if defined(SISLNEEDPROTOTYPES)
   static
      void
	    sh1263_s9makec0(SISLCurve *vbound[],int iedge,int *jstat)
#else	       
static void sh1263_s9makec0(vbound,iedge,jstat)
   SISLCurve *vbound[];
   int iedge;
   int *jstat;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : The produced cross tangent curves are in some corners not
*              consistent with the tangents of the adjacent curves.
*              Change the cross tangent curves in order to achieve this
*              consistency.
*
*
*
* INPUT      : vbound     - Produced boundary curves to vertex blend.
*              iedge      - Number of edges in vertex region.
*                       
*
* OUTPUT     : vbound     - Boundary curves after modification.
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
* CALLS      : s1221, s6lenght
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 11.92.
*
*********************************************************************
*/
{
   int kstat = 0;                /* Local status variable.       */
   int ki,kj,kk,kl;              /* Counters.                    */
   int kn;                       /* Number of vertices of cross tangent curve.*/
   int kn2;                      /* Half the number of vertices. */
   int k1;                       /* Counter.                     */
   int kdim = vbound[0] -> idim; /* Dimension of geometry space. */
   int kleft = 0;                /* Parameter used in s1221.     */
   double tpar1,tpar2;           /* Parameter values of endpoints of curves. */
   double tref;                  /* Referance value used in comparisement.   */
   double tmult;                 /* Multiplicator used to smoothen curve.    */
   double sdiff[3];              /* Difference between wanted and actual
				    cross tangent vector.                    */
   double sder[9];               /* Result of curve evalutation.             */
   SISLCurve *qcross;            /* Cross tangent curve.                     */
   SISLCurve *qc;                /* Adjacent position curve.                 */
   
   /* Traverse the edges.  */
   
   for (ki=0; ki<iedge; ki++)
   {
      qcross = vbound[2*ki+1];
      
      /* Check if the cross tangent curve is consistent with
	 the adjacent tangent in the startpoint of the curve. */
      /* First fetch adjacent curve.   */
      
      kj = (ki > 0) ? ki-1 : iedge-1;
      qc = vbound[2*kj];
      
      /* Fetch parameter values at the current corner.  */
      
      tpar1 = *(qc->et + qc->in);
      tpar2 = *(qcross->et + qcross->ik - 1);
      
      /* Evaluate curves in corner.  */
      
      s1221(qc,1,tpar1,&kleft,sder,&kstat);
      if (kstat < 0) goto error;
      
      s1221(qcross,0,tpar2,&kleft,sder+2*kdim,&kstat);
      if (kstat < 0) goto error;
      
      tref = s6length(sder+kdim,kdim,&kstat);
      if (kstat < 0) goto error;
      
      /* Check equality in each dimension. */
      
      for (kk=0; kk<kdim; kk++)
	 if (DNEQUAL(tref+sder[kdim+kk],tref-sder[2*kdim+kk])) break;
      
      if (kk < kdim)
      {
	 /* Inconsistence. Change the first coefficient of the
	    cross tangent curve to be equal to the tangent (the cross 
	    tangent curve is k-regular). Then smooth this change
	    out along the cross tangent curve.  */
	 
	 /* The sign of the cross tangent curve and the sign of the
	    tangent is opposite.   */
	 
	 for (kl=0; kl<kdim; kl++) sder[kdim+kl] *= -(double)1.0;
	 
	 /* Find difference between coefficient of cross tangent
	    curve and wanted coefficient.  */
	 
	 s6diff(sder+kdim,qcross->ecoef,kdim,sdiff);
	 
	 /* Change cross tangent curve, most at the corner, then
	    smoothen out the effect of the change. */
	    
	 for (kn=qcross->in, kn2=kn/2, k1=0, tmult=(double)1.0;
	  k1<kn2; k1++, tmult -= (double)1.0/(double)kn2)
	    for (kl=0; kl<kdim; kl++)
	       qcross->ecoef[k1*kdim+kl] += tmult*sdiff[kl];
      }
	    
      /* Check if the cross tangent curve is consistent with
	 the adjacent tangent in the endpoint of the curve. */
      /* First fetch adjacent curve.   */
      
      kj = (ki < iedge-1) ? ki+1 : 0;
      qc = vbound[2*kj];
      
      /* Fetch parameter values at the current corner.  */
      
      tpar1 = *(qc->et + qc->ik - 1);
      tpar2 = *(qcross->et + qcross->in);
      
      /* Evaluate curves in corner.  */
      
      s1221(qc,1,tpar1,&kleft,sder,&kstat);
      if (kstat < 0) goto error;
      
      s1221(qcross,0,tpar2,&kleft,sder+2*kdim,&kstat);
      if (kstat < 0) goto error;
      
      tref = s6length(sder+kdim,kdim,&kstat);
      if (kstat < 0) goto error;
      
      /* Check equality in each dimension. */
      
      for (kk=0; kk<kdim; kk++)
	 if (DNEQUAL(sder[kdim+kk]+tref,sder[2*kdim+kk]+tref)) break;
      
      if (kk < kdim)
      {
	 /* Inconsistence. Change the last coefficient of the
	    cross tangent curve to be equal to the tangent (the cross 
	    tangent curve is k-regular). Then smooth this change
	    out along the cross tangent curve.  */
	 
	 /* Find difference between coefficient of cross tangent
	    curve and wanted coefficient.  */
	 
	 s6diff(sder+kdim,qcross->ecoef+(qcross->in-1)*kdim,kdim,sdiff);
	 
	 /* Change cross tangent curve, most at the corner, then
	    smoothen out the effect of the change. */
	    
	 for (kn=qcross->in, kn2=kn/2, k1=kn-1, tmult=(double)1.0;
	  k1>kn2; k1--, tmult -= (double)1.0/(double)kn2)
	    for (kl=0; kl<kdim; kl++)
	       qcross->ecoef[k1*kdim+kl] += tmult*sdiff[kl];
      }
      
   }

  *jstat = 0;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  goto out;
  
  out:
     return;
}     

   
#if defined(SISLNEEDPROTOTYPES)
   static
      void
	    sh1263_s9checkcrtan(SISLCurve *pcrtan,int *jstat)
#else	       
static void sh1263_s9checkcrtan(pcrtan,jstat)
     int *jstat;
     SISLCurve *pcrtan;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Check that the current cross tangent curve will not give
*              rise to a surface with intersecting constant parameter
*              curves in one parameter direction.
*
*
*
* INPUT/OUTPUT : pcrtan   - Cross tangent curve.
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
* CALLS      : 
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{

  *jstat = 0;
  
  return;
}     

   
   





