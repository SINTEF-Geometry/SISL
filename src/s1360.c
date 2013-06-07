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
 * $Id: s1360.c,v 1.4 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1360

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1360(SISLCurve *pc,double aoffset,double aepsge,double enorm[],
	   double amax,int idim,SISLCurve **rc,int *jstat)
#else
void s1360(pc,aoffset,aepsge,enorm,amax,idim,rc,jstat)
     SISLCurve  *pc;
     double aoffset;
     double aepsge;
     double enorm[];
     double amax;
     int    idim;
     SISLCurve  **rc;
     int    *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To create a B-spline approximating the offset curve of
*              a NURBS curve. The curve must have continous first
*              derivative if the offset is different from zero. If the
*              offset is zero then the position must be continous.
*                           *
*
* INPUT      : pc     - The input NURBS curve.
*              aoffset- The offset distance.
*                       If idim=2 a positive offset value will place the
*                       offset curve on the positive side of the normal vector,
*                       and a negative value places the offset curve on the
*                       negative side of the normal vector.
*                       If idim=3 the offset is determined by the cross
*                       product of the tangent vector and the anorm vector.
*                       The offset distance is multiplied by this vector.
*              enorm  - normal vector
*              aepsge - Maximal deviation allowed between true offset curve
*                       and the approximated offset curve.
*              amax   - Maximal stepping length. Is negleceted if amax<=aepsge
*                       If amax==0 then a maximal step length of the longest
*                       SISLbox side is used.
*              idim   - The dimension of the space (2 or 3).
*
* OUTPUT     :
*              jstat  - status messages
*                            > 0      : warning
*                            = 0      : ok
*                            = -1     : curve is degenrate.
*				        nothing returned.
*                            < 0      : error
*              rc     - Pointer the approximated offset curve
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
* CALLS      : s1307, s1311, s1361, s1362, s6norm, s6length,
*              s6dist,s6ang,s6scpr,s6diff,s6err
*
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 24. May 1988
* REVISED BY : Michael Floater, SI, Oslo, 6/2/92.
*                   Give error if pc is degenerate.
* CHANGED BY : Ulf J. Krystad, Oslo, Norway. April 1992
*              call to s1312 changed to call to s1359.
*********************************************************************
*/
{
  int kmaxinf;        /* Number of vertices space is allocated for       */
  int knbinf=0;       /* Number of points stored so far                  */
  int kder=2;         /* Derivative indicator                            */
  int kstat,kstat1;   /* Status variable                                 */
  int kstat2;
  int kleft=0;        /* Left indicator for point calculation            */
  int kleftend=0;     /* Pointer telling left knot of end of last segment*/
  int kmult;          /* Multiplicity of knot at parameter value         */
  int kn;             /* The number of B-splines, i.e., the dimension of
			 the spline space associated with the knot
			 vector.                                         */
  int kk;             /* The polynomial order of the curve.              */
  int kdim;           /* The dimension of the space in which the curve
			 lies. Equivalently, the number of components
			 of each B-spline coefficient.                   */
  int notaccepted;    /* Loop control variable                           */
  int kcont;          /* Loop control variable                           */
  int kdiv;           /* Divergence indicator                            */
  int knbit;          /* Number of iterations                            */
  int kpos=0;         /* Position of error                               */
  int kpar = 1;       /* Indicate that parametrization array exist       */
  double sder1[9];    /* Derivatives from right                          */
  double sder2[9];    /* Derivatives from left                           */
  double smidd[3];    /* Middle point of current Bezier segement         */
  double smtang[3];   /* Tangent at smidd                                */
  double sdiff[3];    /* Difference of vectors                           */
  double tproj1;      /* Projection of vector                            */
  double tproj2;      /* Projection of vector                            */
  double tlast;       /* Value in last itertion                          */
  double tfak;        /* Necessary reduction of interval length          */
  double *s3dinf=SISL_NULL;/* Pointer to storage for point info (10 dobules pr
			 point when idim=3, 7 when idim=3)               */
  double *spar=SISL_NULL;  /* Pointer to array for storage of knots           */
  double *st;         /* Pointer to the first element of the knot vector
			 of the curve. The knot vector has [kn+kk]
			 elements.                                       */
  double *scoef;      /* Pointer to the first element of the curve's
			 B-spline coefficients. This is assumed to be an
			 array with [kn*kdim] elements stored in the
			 following order:
			 First the kdim components of the first B-spline
			 coefficient, then the kdim components of the
			 second B-spline coefficient and so on.          */
  double sder[9];     /* Poisition, first and second derivative on curve */
  double tx1,tx2,txm; /* Parameter value */
  double tstep;       /* Step length     */
  double tmax;        /* Local maximal step length                 */
  double tlengthend;  /* Length of 1st derivative at end of segment */
  double tincre;      /* Parameter value increment */
  double *start;      /* Pointer to start of current segment */
  double tdist;       /* Distance */
  double tang;        /* Angle */
  double tnew;        /* New increment */
  double tdum;        /* Temporary variable */
  double tl1,tl2;     /* Vector lengths */
  double tscal;       /* Scalar product */

  /* Test that dimension is 2 or 3 */
  if (idim !=2 && idim !=3) goto err105;

  if (aepsge <= DZERO) goto err184;



  /* Make maximal step length based on box-size of surface */

  sh1992cu(pc,0,aepsge,&kstat);
  if (kstat < 0) goto error;

  tmax = MAX(pc->pbox->e2max[0][0] - pc->pbox->e2min[0][0],
	     pc->pbox->e2max[0][1] - pc->pbox->e2min[0][1]);
  tmax = MAX(tmax,pc->pbox->e2max[0][2] - pc->pbox->e2min[0][2]);

  /* Check for degenerate curve. */

  if(tmax == DZERO)
  {
      /* Curve is degenerate -- pc is just a point. */
      *jstat = -1;
      goto out;
  }

  if (amax>DZERO) tmax = MIN(tmax,amax);
  /* Copy curve attributes to local parameters.  */

  kn    = pc -> in;
  kk    = pc -> ik;
  st    = pc -> et;
  scoef = pc -> ecoef;
  kdim  = pc -> idim;

  kmaxinf = 100;

  /* Allocate space for storage of points,tangents, curvature and radius of
     curvature */

  s3dinf = newarray((3*kdim+1)*kmaxinf,DOUBLE);
  if (s3dinf == SISL_NULL) goto err101;

  /* Allocate space for parametrization array */

  spar = newarray(kmaxinf,DOUBLE);
  if (spar == SISL_NULL) goto err101;

  /* Store knot values at start of curve */

  tx1     = st[kk-1];
  spar[0] = tx1;


  /* Make start point and intital step length */

  kder = 2;

  s1362(pc,aoffset,enorm,idim,kder,tx1,&kleftend,sder,&kstat);
  if (kstat<0) goto error;

  /* Calculate unit tangent and radius of curvature */

  s1307(sder,kdim,s3dinf,&kstat);
  if (kstat<0) goto error;
  knbinf = 1;

  /* Calculate step length based on curvature */

  tstep = s1311(s3dinf[3*kdim],aepsge,tmax,&kstat);
  if (kstat<0) goto error;

  /* Remember length of start tangent, end of zero segment */

  tlengthend = s6length(sder+kdim,kdim,&kstat);
  if (kstat<0) goto error;

  /* While end not reached */


  while (tx1 < st[kn])
    {

      /*  Find candidate end point, make sure that no breaks in tangent or
	  curvature exists between start and endpoints of the segment      */

      /*  Make step length equal to aepsge if the length is zero */

      /*  Find parameter value of candidate end point of segment */

      if (DEQUAL(tlengthend,DZERO))
        {
	  /* Step equal to computer resolution */
	  tincre = tx1*((double)1.0+REL_COMP_RES);
        }
      else
        tincre = tstep/tlengthend;

      /*  Make sure that we don't pass any knots */

      tx2 = MIN(tx1 + tincre,st[kleftend+1]);


      /*  While segement not accepted */

      notaccepted = 1;

      while(notaccepted==1)
        {

	  /* Make end point of segment, and store it */

	  if (knbinf+2>=kmaxinf)
            {
	      kmaxinf = kmaxinf + 100;
	      s3dinf = increasearray(s3dinf,(3*kdim+1)*kmaxinf,DOUBLE);
	      spar   = increasearray(spar,kmaxinf,DOUBLE);
            }

	  /* Check multiplicity of end point of segment make sure that if
	     the multiplicity is greater than kk that we calculate derivatives
	     from the left. */

	  kmult = s6knotmult(st,kk,kn,&kleft,tx2,&kstat);
	  if (tx2 >= st[kn]) kmult = 0;
	  if (kmult > kk-2)
	    {
	      /* Check right and left hand derivatives of curve */

	      s1227(pc,2,tx2,&kleft,sder2,&kstat);
	      if (kstat<0) goto error;

	      s1221(pc,2,tx2,&kleft,sder1,&kstat);
	      if (kstat<0) goto error;

	      tl1 = s6length(sder1+kdim,kdim,&kstat1);
	      tl2 = s6length(sder2+kdim,kdim,&kstat2);
	      tscal = s6scpr(sder1+kdim,sder2+kdim,kdim);
	      tdum = MAX(tl1,tl2);
	      if (kstat1<0 || kstat2<0) goto error;

	      if(DEQUAL(tl1*tl2+tdum,tscal+tdum)) kmult = 0; /* G1 continuity*/
	    }

	  /* Make correct value of kleftend */

	  s1219(st,kk,kn,&kleftend,tx2,&kstat);
	  if (kstat<0) goto error;


	  /* Real offset, test multiplicity before calculation */

          if (kmult > kk-2) tdum = tx2 * ((double)1.0-REL_COMP_RES);
          else tdum = tx2;

          /* We have the right value of kleftend, if no offset and
             break then s1362 will decrease kleft */

          if (aoffset != DZERO)
	    {
              s1362(pc,aoffset,enorm,idim,kder,tdum,&kleft,sder,&kstat);
	      if (kstat<0) goto error;
	    }
          else if (kmult > kk-2)
	    memcopy(sder,sder2,kdim*3,DOUBLE); /* Copy left derivatives */
          else
	    {
	      s1221(pc,2,tx2,&kleft,sder,&kstat);
	      if (kstat<0) goto error;
	    }

	  /* Remember length of start tangent, end of zero segment */

	  tlengthend = s6length(sder+kdim,kdim,&kstat);
	  if (kstat<0) goto error;



	  /* Calculate unit tangent and radius of curvature */

	  s1307(sder,kdim,s3dinf+(3*kdim+1)*knbinf,&kstat);
	  if (kstat<0) goto error;

	  /* Decide if Hermit shape acceptable and find position and tangent
	     at midpoint of segment */

	  start = s3dinf + (3*kdim+1)*(knbinf-1);

	  s1361(start,start+(3*kdim+1),kdim,smidd,smtang,&kstat);
	  if (kstat<0) goto error;

	  /* Iterate to midpoint of segment, start from middle of [tx1,tx2].
	     The iteration is performed to find the intersection between the
	     plane described by smidd and smtang. */

	  txm = (tx1+tx2)/(double)2.0;

	  kcont = 1;
	  kdiv  = 0;

	  knbit = 0;
	  while (kcont==1)
            {

	      /* Calculate position and tangent at txm */

	      s1362(pc,aoffset,enorm,idim,kder,txm,&kleft,sder,&kstat);
	      if (kstat<0) goto error;

	      /* Make difference of calculated point and smidd, project this
		 onto the normal of the plane. */

	      s6diff(sder,smidd,kdim,sdiff);
	      tproj1 = s6scpr(sdiff,smtang,kdim);
	      tproj2 = s6scpr(&sder[kdim],smtang,kdim);

	      /* If tproj2==0 then curve tangent normal to plane, half step
		 length */

	      if (DEQUAL(tproj2,DZERO))
                {
		  kdiv  = 1;
		  kcont = 0;
                }
	      else if (knbit==0)
                {
		  /* First iteration */
		  knbit = 1;
		  txm -= tproj1/tproj2;
		  txm = MAX(tx1,txm);
		  txm = MIN(tx2,txm);
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
		  txm   -= tproj1/tproj2;
		  tlast  = fabs(tproj1);
		  knbit += 1;
		  if (txm <=tx1 || tx2 <= txm)
                    {
		      kdiv  = 1;
		      kcont = 0;
                    }
                }
            }
	  /* Find how close the midpoint position and tangent of the segement
	     is to the true curve */

	  tdist = s6dist(sder,smidd,kdim);

	  tang  = s6ang(&sder[kdim],smtang,kdim);

	  /* If the point is not within the resolution treat it as divergence, except
	     when the segment length is less than aepsge */

	  if ((s6dist(start,smidd,kdim) <= aepsge) &&
	      (s6dist(start+(3*kdim+1),smidd,kdim) <= aepsge))
	    {
	      kdiv = 0;
	    }
	  else if (fabs(tdist) > aepsge || fabs(tang) > ANGULAR_TOLERANCE)
            {
	      kdiv = 1;
            }

	  /* Dependent on previous conditions decide if the segment is acceptable
	     or not */

	  if (kdiv==0)
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

      /*  Store segment information */

      if (kmult>kk-2)
	{
          /* End of segment is possible break, add exstra version of point */

          if (aoffset != DZERO)
	    {
	      spar[knbinf] = tx2-(double)0.1*(tx2-spar[knbinf-1]);
	      s1362(pc,aoffset,enorm,idim,kder,tx2,&kleft,sder,&kstat);
	      if (kstat<0) goto error;
	    }
          else
	    {
	      spar[knbinf] = tx2;
	      memcopy(sder,sder1,3*kdim,DOUBLE);/*Copy right hand derivatives*/
	    }

          knbinf++;

	  /* Remember length of start tangent, end of zero segment */

	  tlengthend = s6length(sder+kdim,kdim,&kstat);
	  if (kstat<0) goto error;

	  /* Calculate unit tangent and radius of curvature */

	  s1307(sder,kdim,s3dinf+(3*kdim+1)*knbinf,&kstat);
	  if (kstat<0) goto error;
	}

      /*  Make knots */
      spar[knbinf] = tx2;

      if (kstat<0) goto error;
      knbinf += 1;

      /*  Make new step length */

      /*  Calculate step length based on curvature */

      tstep = s1311(s3dinf[(3*kdim+1)*knbinf-1],aepsge,tmax,&kstat);
      if (kstat<0) goto error;

      /*  Update start parameter value of segment */

      tx1 = tx2;
    }


  /*  Interpolate offset curve */

  /* UJK, 92.04. : s1312 and s1359 shadow functions
     s1312(s3dinf,kdim,knbinf,kpar,spar,rc,&kstat); */
  s1359(s3dinf,aepsge,kdim,knbinf,kpar,spar,rc,&kstat);
  if (kstat < 0) goto error;

  *jstat = 0;
  goto out;

  /* Error in memory allocation */

 err101:
  *jstat = -101;
  s6err("s1360",*jstat,kpos);
  goto out;

  /* Error in input, dimension not equal to 2 or 3 */

 err105:
  *jstat = -105;
  s6err("s1360",*jstat,kpos);
  goto out;


 err179:
  *jstat = -179;
  s6err("s1360",*jstat,kpos);
  goto out;

  /* Negative aepsge */

 err184:
  *jstat = -184;
  s6err("s1360",*jstat,kpos);
  goto out;

  /* Error in lower level function */

 error:
  *jstat = kstat;
  s6err("s1360",*jstat,kpos);
  goto out;

 out:
  if (s3dinf != SISL_NULL) freearray(s3dinf);
  if (spar   != SISL_NULL) freearray(spar);
  return;
}
