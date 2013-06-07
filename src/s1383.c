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
 * $Id: s1383.c,v 1.5 2009-03-18 13:30:55 vsk Exp $
 *
 */


#define S1383

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1383(SISLSurf *psurf,SISLCurve *pcurv,double aepsge,double amax,int ider,
	   SISLCurve **rcpos,SISLCurve **rcder1,SISLCurve **rcder2,int *jstat)
#else
void s1383(psurf,pcurv,aepsge,amax,ider,rcpos,rcder1,rcder2,jstat)
     SISLSurf   *psurf;
     SISLCurve  *pcurv;
     double aepsge;
     double amax;
     int    ider;
     SISLCurve  **rcpos;
     SISLCurve  **rcder1;
     SISLCurve  **rcder2;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To create a 3-D B-spline approximating the curve traced out
*              by a curve in the parameter plane.             
*
* INPUT      : psurf  - The surface object
*              pcurv  - The input B-spline curve in the parameter plane   
*              aepsge - Maximal deviation allowed between true 3-D curve
*                       and the approximated 3-D curve.
*              amax   - Maximal stepping length. Is negleceted if amax<=aepsge
*                       If amax <= 0.0 the 3-D SISLbox of the surface us used
*                       for estimating max step length
*              ider   - Derivativ indicator
*                        0 - Calculate only psositional curve
*                        1 - Calculate positional + derivative curves
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              rcpos  - Pointer the approximated position curve
*              rcder1 - Pointer the approximated position curve
*              rcder2 - Pointer the approximated position curve
*
* METHOD     : 
*
* EXAMPLE OF USE:
*              SISLCurve *qr;         
*              int    kstat;
*              .
*              .
*
* REFERENCES :
*
*-                                                 
*              
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. November 1988
* CHANGED BY : Ulf J. Krystad, Oslo, Norway. April 1992
*              call to s1312 changed to call to s1359.
*
*********************************************************************
*/
{
  int kmaxinf;        /* Number of vertices space is allocated for       */
  int knbinf=0;       /* Number of points stored so far                  */
  int kstat;          /* Status variable                                 */
  int kleftc=0;       /* Left indicator for point calculation            */
  int klefts1=0;      /* Left indicator for point calculation            */
  int klefts2=0;      /* Left indicator for point calculation            */
  int kn;             /* The number of B-splines, i.e., the dimension of
			 the spline space associated with the knot
			 vector.                                         */
  int kk;             /* The polynomial order of the curve.              */
  int kdimc;          /* The dimension of the space in which the curve
			 lies. Equivalently, the number of components
			 of each B-spline coefficient.                   */
  int kdims;          /* Dimension of space where the surface lies       */
  int notaccepted;    /* Loop control variable                           */
  int kcont;          /* Loop control variable                           */
  int kdiv;           /* Divergence indicator                            */
  int knbit;          /* Number of iterations                            */
  int kpos=0;         /* Position of error                               */
  int kpar = 1;       /* Indicate that parametrization array exist       */
  int kmult;           /* Multiplicity of knot                           */
  double smidd[3];    /* Middle point of current Bezier segement         */
  double smtang[3];   /* Tangent at smidd                                */
  double sdiff[3];    /* Difference of vectors                           */
  double tproj1;      /* Projection of vector                            */
  double tproj2;      /* Projection of vector                            */
  double tlast;       /* Value in last itertion                          */
  double tfak;        /* Necessary reduction of interval length          */
  double *s3dinf=SISL_NULL;/* Pointer to storage for point info (10 dobules pr
			 point when idim=3, 7 when idim=3)              */
  double *sudpos=SISL_NULL;/* Pointer to storage of u-derivatives             */
  double *sudder=SISL_NULL;/* Pointer to storage of ut-derivatives            */
  double *svdpos=SISL_NULL;/* Pointer to storage of u-derivatives             */
  double *svdder=SISL_NULL;/* Pointer to storage of vt-derivatives            */
  double *spar=SISL_NULL;  /* Pointer to array for storage of knots           */
  double *st;         /* Pointer to the first element of the knot vector
			 of the curve. The knot vector has [kn+kk]
			 elements.                                       */
  double sder[9];     /* Position, first and second derivative on curve  */
  double sderu[6];    /* Position and derivative of u-derivative         */
  double sderv[6];    /* Position and derivative of v-derivative         */
  double sdern[6];    /* Position and derivative of normal vector        */
  double tx1,tx2,txm; /* Parameter value */
  double tstep;       /* Step length     */
  double tlengthend;  /* Length of 1st derivative at end of segment */
  double tincre;      /* Parameter value increment */
  double *start;      /* Pointer to start of current segment */
  double tmax;        /* Local maximal step length                       */
  double tdist;       /* Distance */
  double tang;        /* Angle */
  double tnew;        /* New increment */
  double tlength;     /* Estimated length of current curve piece */
  
  
  /* Make maximal step length based on box-size of surface */
  
   sh1992su(psurf,0,aepsge,&kstat);
   if (kstat < 0) goto error;
  
  tmax = MAX(psurf->pbox->e2max[0][0] - psurf->pbox->e2min[0][0],
	     psurf->pbox->e2max[0][1] - psurf->pbox->e2min[0][1]);
  tmax = MAX(tmax,psurf->pbox->e2max[0][2] - psurf->pbox->e2min[0][2]);
  
  if (amax>DZERO) tmax = MIN(tmax,amax);
  
  /* Copy curve attributes to local parameters.  */
  
  kdimc = pcurv -> idim;
  kk    = pcurv -> ik;
  kn    = pcurv -> in;
  st    = pcurv -> et;
  kdims = psurf -> idim;
  
  if (kdimc != 2 || kdims != 3) goto err105;
  
  kmaxinf = 100;
  
  /* Allocate space for storage of points,tangents, curvature and radius of
     curvature */
  
  s3dinf = newarray((3*kdims+1)*kmaxinf,DOUBLE);                               
  if (s3dinf == SISL_NULL) goto err101;

  if (ider >= 1)
    {
      sudpos = newarray(kdims*kmaxinf,DOUBLE);                               
      if (sudpos == SISL_NULL) goto err101;
      sudder = newarray(kdims*kmaxinf,DOUBLE);                               
      if (sudder == SISL_NULL) goto err101;
      svdpos = newarray(kdims*kmaxinf,DOUBLE);                               
      if (svdpos == SISL_NULL) goto err101;
      svdder = newarray(kdims*kmaxinf,DOUBLE);                               
      if (svdder == SISL_NULL) goto err101;
    }
  
  /* Allocate space for parametrization array */
  
  spar = newarray(kmaxinf,DOUBLE);
  if (spar == SISL_NULL) goto err101;
  
  /* Store knot values at start of curve */
  
  tx1 = st[kk-1];
  spar[0] = tx1;
  
  
  /* Make start point and intital step length */
  
  s1384(pcurv,psurf,kdims,1,tx1,&kleftc,&klefts1,&klefts2,
	sder,sderu,sderv,sdern,&kstat);
  if (kstat<0) goto error;
  
  /* Calculate unit tangent and radius of curvature */
  
  s1307(sder,kdims,s3dinf,&kstat);
  if (kstat<0) goto error;
  knbinf = 1;                   
  
  /* Store the other calculated information */
  
  if (ider>=1)
    {
      memcopy(sudpos,sderu,kdims,DOUBLE);
      memcopy(sudder,sderu+kdims,kdims,DOUBLE);
      memcopy(svdpos,sderv,kdims,DOUBLE);
      memcopy(svdder,sderv+kdims,kdims,DOUBLE);
    }
  
  /* Calculate step length based on curvature */
  
  tstep = s1311(s3dinf[3*kdims],aepsge,tmax,&kstat);
  if (kstat<0) goto error;
  
  /* Remember length of start tangent, end of zero segment */
  
  tlengthend = s6length(sder+kdims,kdims,&kstat);
  if (kstat<0) goto error;                                            
  
  /* While end not reached */
  
  
  while (tx1 < st[kn])
    {
      
      /* Find candidate end point, make sure that no breaks in tangent or
	 curvature exists between start and endpoints of the segment      */
      
      /* Make step length equal to aepsge if the length is zero */
      
      /* Find parameter value of candidate end point of segment */
      
      if (DEQUAL(tlengthend,DZERO))
        { 
	  /* Step equal to parameter resolution */
	  tincre = max(tx1*((double)1.0+REL_PAR_RES),REL_PAR_RES);
        }
      else
        tincre = tstep/tlengthend;
      
      /* Make sure that we don't pass any knots */
      
      tx2 = MIN(tx1 + tincre,st[kleftc+1]);
      
      /* While segement not accepted */
      
      notaccepted = 1;
      
      while(notaccepted==1)
        {
	  
	  /* Make end point of segment, and store it */
	  
	  if (knbinf>=kmaxinf)
            {
	      kmaxinf = kmaxinf + 100;
	      s3dinf = increasearray(s3dinf,(3*kdims+1)*kmaxinf,DOUBLE);
	      if (s3dinf == SISL_NULL) goto err101;
	      spar   = increasearray(spar,kmaxinf,DOUBLE);
	      if (spar == SISL_NULL) goto err101;

	      if (ider >= 1)
		{
		  sudpos = increasearray(sudpos,kdims*kmaxinf,DOUBLE);
		  if (sudpos == SISL_NULL) goto err101;
		  sudder = increasearray(sudder,kdims*kmaxinf,DOUBLE);
		  if (sudder == SISL_NULL) goto err101;
		  svdpos = increasearray(svdpos,kdims*kmaxinf,DOUBLE);
		  if (svdpos == SISL_NULL) goto err101;
		  svdder = increasearray(svdder,kdims*kmaxinf,DOUBLE);
		  if (svdder == SISL_NULL) goto err101;
		}
            }
	  
	  
	  s1384(pcurv,psurf,kdims,-1,tx2,&kleftc,&klefts1,&klefts2,
		sder,sderu,sderv,sdern,&kstat);
	  if (kstat<0) goto error;
	  
	  
	  /* Remember length of start tangent, end of zero segment */
	  
	  tlengthend = s6length(sder+kdims,kdims,&kstat);
	  if (kstat<0) goto error;
	  
	  
	  
	  /* Calculate unit tangent and radius of curvature */
	  
	  s1307(sder,kdims,s3dinf+(3*kdims+1)*knbinf,&kstat);
	  if (kstat<0) goto error;
	  
	  /* Store the other calculated information */
	  
	  if (ider>=1)
            {
	      memcopy(&sudpos[kdims*knbinf],sderu,kdims,DOUBLE);
	      memcopy(&sudder[kdims*knbinf],sderu+kdims,kdims,DOUBLE);
	      memcopy(&svdpos[kdims*knbinf],sderv,kdims,DOUBLE);
	      memcopy(&svdder[kdims*knbinf],sderv+kdims,kdims,DOUBLE);
            }
	  
	  
	  /* Decide if Hermit shape acceptable and find position and tangent
	     at midpoint of segment */
	  
	  start = s3dinf + (3*kdims+1)*(knbinf-1);
	  
	  s1361(start,start+(3*kdims+1),kdims,smidd,smtang,&kstat);
	  if (kstat<0) goto error;
	  
	  /* Iterate to midpoint of segment, start from middle of [tx1,tx2].
	     The iteration is performed to find the intersection between the
	     plane described by smidd and smtang. */
	  
	  txm = (tx1+tx2)/(double)2.0;
	  
	  kcont = 1;
	  kdiv    = 0;
	  
	  knbit = 0;
	  while (kcont==1)
            {
	      
	      /*  Calculate position and tangent at txm */
	      
	      
	      
	      s1384(pcurv,psurf,kdims,-1,txm,&kleftc,&klefts1,&klefts2,
		    sder,sderu,sderv,sdern,&kstat);
	      if (kstat<0) goto error;
	      
	      /* Make difference of calculated point and smidd, project this
		 onto the normal of the plane. */
	      
	      s6diff(sder,smidd,kdims,sdiff);
	      tproj1 = s6scpr(sdiff,smtang,kdims);
	      tproj2 = s6scpr(&sder[kdims],smtang,kdims);
	      
	      /* If tproj2==0 then curve tangent normal to plane, half step
		 length */
	      
	      if (DEQUAL(tproj2,DZERO))
                {
		  kdiv = 1;
		  kcont = 0;
                }
	      else if (knbit==0)
                {
		  /* First iteration */
		  knbit = 1;
		  txm -= tproj1/tproj2;
		  tlast = fabs(tproj1);
                }
	      
	      else if (fabs(tproj1)>=tlast)
                {
		  /* Not convergence any longer */
		  kcont = 0;
                }
	      else
                {
		  /* Still convergence */
		  txm -= tproj1/tproj2;
		  tlast = fabs(tproj1);
		  knbit += 1;
		  if (txm <=tx1 || tx2 <= txm)
                    {
		      kdiv = 1;
		      kcont = 0;
                    }
                }
            }
	  /* Find how close the midpoint position and tangent of the segement
	     is to the true curve */
	  
	  tdist = s6dist(sder,smidd,kdims);
	  
	  tang  = s6ang(&sder[kdims],smtang,kdims);

	  tlength = s6dist(start,smidd,kdims) + 
	    s6dist(start+3*kdims+1,smidd,kdims);
	  
	  /* If the point is not within the resolution treat it as divergence
	   */
	  if (fabs(tdist) > aepsge || fabs(tang) > ANGULAR_TOLERANCE)
            {
	      kdiv = 1;
            }
	  
	  /* Dependent on previous conditions decide if the segment 
	     is acceptable or not */
	  
	  if (kdiv==0 || tlength < (double)2*aepsge)
            {
	      /* Segement acceptable */
	      notaccepted = 0;
            }
	  else
            {
	      /* Segment unacceptable. Remember that the error of the Hermit
		 interpolation is O(h**4). Thus taking this into consideration
		 we can determin the new parameter interval. */
	      
	      tfak = MAX(tdist/aepsge,(double)1.0);
	      tfak = (double)2.0*pow(tfak,ONE_FOURTH);
	      
	      tnew = MIN(tincre/(double)2.0,(tx2-tx1)/tfak);
	      if (DEQUAL(tmax+tnew,tmax+tincre)) goto err179;
	      tincre = tnew;
	      tx2 = tx1 + tincre;
            }
        }
      
      /* Store segment information */
      
      /* Make knots */
      spar[knbinf] = tx2;
      
      if (kstat<0) goto error;
      knbinf += 1;
      
      /* Make new step length */ 
      
      /* Calculate step length based on curvature */
      
      tstep = s1311(s3dinf[(3*kdims+1)*knbinf-1],aepsge,tmax,&kstat);
      if (kstat<0) goto error;
      
      /* Update start parameter value of segment, and calculate right
	 hand derivative */
      
      tx1 = tx2;
      
      s1384(pcurv,psurf,kdims,1,tx1,&kleftc,&klefts1,&klefts2,
	    sder,sderu,sderv,sdern,&kstat);
      if (kstat<0) goto error;
      
      kmult = s6knotmult(st,kk,kn,&kleftc,tx1,&kstat);
      if (kstat<0) goto error;
      
      if (kmult>=kk-1)
        {
	  if (knbinf>=kmaxinf)
            {
	      kmaxinf = kmaxinf + 100;
	      s3dinf = increasearray(s3dinf,(3*kdims+1)*kmaxinf,DOUBLE);
	      if (s3dinf == SISL_NULL) goto err101;
	      spar   = increasearray(spar,kmaxinf,DOUBLE);
	      if (spar == SISL_NULL) goto err101;

	      if (ider >= 1)
		{
		  sudpos = increasearray(sudpos,kdims*kmaxinf,DOUBLE);
		  if (sudpos == SISL_NULL) goto err101;
		  sudder = increasearray(sudder,kdims*kmaxinf,DOUBLE);
		  if (sudder == SISL_NULL) goto err101;
		  svdpos = increasearray(svdpos,kdims*kmaxinf,DOUBLE);
		  if (svdpos == SISL_NULL) goto err101;
		  svdder = increasearray(svdder,kdims*kmaxinf,DOUBLE);
		  if (svdder == SISL_NULL) goto err101;
		}
            }
	  
	  /* Remember length of start tangent, end of zero segment */
	  
	  tlengthend = s6length(sder+kdims,kdims,&kstat);
	  if (kstat<0) goto error;
	  
	  /* Calculate unit tangent and radius of curvature */
	  
	  s1307(sder,kdims,s3dinf+(3*kdims+1)*knbinf,&kstat);
	  if (kstat<0) goto error;
	  
	  /* Store the other calculated information */
	  
	  if (ider>=1)
            {
	      memcopy(&sudpos[kdims*knbinf],sderu,kdims,DOUBLE);
	      memcopy(&sudder[kdims*knbinf],sderu+kdims,kdims,DOUBLE);
	      memcopy(&svdpos[kdims*knbinf],sderv,kdims,DOUBLE);
	      memcopy(&svdder[kdims*knbinf],sderv+kdims,kdims,DOUBLE);
            }
	  
	  spar[knbinf] = spar[knbinf-1];
	  
	  knbinf += 1;
	  
        }
    }
  
  
  /*  Interpolate trace curve */
  
  /* UJK, 92.04. : s1312 and s1359 shadow functions
     s1312(s3dinf,kdims,knbinf,kpar,spar,rcpos,&kstat); */
  s1359(s3dinf,aepsge,kdims,knbinf,kpar,spar,rcpos,&kstat);
  if (kstat < 0) goto error;
  
  if (ider>=1)
    {
      /*   Interpolate u-derivative curve */
      s1379(sudpos,sudder,spar,knbinf,kdims,rcder1,&kstat);
      if (kstat<0) goto error;
      
      /*   Interpolate v-derivative curve */
      s1379(svdpos,svdder,spar,knbinf,kdims,rcder2,&kstat);
      if (kstat<0) goto error;
    }                
  
  
  *jstat = 0;
  goto out;
  
  /* Error in memory allocation */
  
 err101: *jstat = -101;
  s6err("s1383",*jstat,kpos);
  goto out;
  
  /* Error in input, dimension not equal to 2 or 3 */
  
 err105: *jstat = -105;
  s6err("s1383",*jstat,kpos);
  goto out;
    
  
 err179: *jstat = -179;
  s6err("s1383",*jstat,kpos);
  goto out;
  
  /* Error in lower level function */
  
 error:  *jstat = kstat;
  s6err("s1383",*jstat,kpos);
  goto out;
  
  
 out:
  
  if (s3dinf != SISL_NULL) freearray(s3dinf);
  if (sudpos != SISL_NULL) freearray(sudpos);
  if (sudder != SISL_NULL) freearray(sudder);
  if (svdpos != SISL_NULL) freearray(svdpos);
  if (svdder != SISL_NULL) freearray(svdder);
  if (spar   != SISL_NULL) freearray(spar);
  
  return;
}          
                                      
