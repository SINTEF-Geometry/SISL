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
 * $Id: sh1850.c,v 1.2 2001-03-19 15:59:06 afr Exp $
 *
 */


#define SH1850

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
sh1850(SISLCurve *pc1,double epoint[],double enorm[],int idim,double aepsco,
       double aepsge,int trackflag, int *jtrack, SISLTrack *** wtrack,
       int *jpt,double **gpar,int **pretop,int *jcrv,
       SISLIntcurve ***wcurve,int *jstat)
#else
void sh1850(pc1,epoint,enorm,idim,aepsco,aepsge,
	    trackflag,jtrack,wtrack,jpt,gpar,pretop,jcrv,wcurve,jstat)
     SISLCurve    *pc1;
     double   epoint[];
     double   enorm[];
     int      idim;
     double   aepsco;
     double   aepsge;
     int       trackflag;
     int       *jtrack;
     SISLTrack ***wtrack;
     int      *jpt;
     double   **gpar;
     int      **pretop;
     int      *jcrv;
     SISLIntcurve ***wcurve;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find all intersections between a curve and a plane if
*              the dimension is equal to three and a line if the 
*              dimension is two.
*
*
*
* INPUT      : pc1    - Pointer to the curve.
*              epoint - SISLPoint in the plane/line.
*              enorm  - Normal to the plane/line.
*              idim   - Dimension of the space in which the plane/line
*                       lies. idim should be equal to two or three.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*              trackflag - If true, create tracks.
*
*
*
* OUTPUT     : jtrack - Number of tracks created
*              wtrack - Array of pointers to tracks
*              jpt    - Number of single intersection points.
*              gpar   - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       interval of the curve. The points lie continuous. 
*                       Intersection curves are stored in wcurve.
*              pretop - Topology info. for single intersection points.
*              jcrv   - Number of intersection curves.
*              wcurve - Array containing descrjptions of the intersection
*                       curves. The curves are only described by points
*                       in the parameter interval. The curve-pointers points
*                       to nothing. (See descrjption of Intcurve
*                       in intcurve.dcl).
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The vertices of the curve is put into the equation of the
*              plane/line achieving a curve in the one-dimentional space.
*              The zeroes of this curve is found.
*
*
* REFERENCES :
*
*- 
* CALLS      : sh1761       - Find point/curve intersections.
*              hp_s1880       - Put intersections into output format.
*              s6diff      - Difference between two vectors.
*              s6scpr      - Scalar product of two vectors.
*              newCurve    - Create new curve.                        
*              newPoint    - Create new point.
*              newObject   - Create new object.
*              freeObject  - Free the space occupied by a given object.
*              freeIntdat  - Free space occupied by an intersection data.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
* REWRITTEN BY : Bjoern Olav Hoset, SI, 89-06.
* REVISED BY : Mike Floater, SI, 91-04 for rational curves.
* REVISED BY : UJK, SI, 92-04 periodicity.
* 
*********************************************************************
*/                                                               
{                                                                     
  double *nullp = SISL_NULL;
  int ikind;               /* kind of curve pc1 is                     */
  int kstat = 0;           /* Local status variable.                     */
  int kpos = 0;            /* Position of error.                         */
  int kn;                  /* Number of vertices of curve.               */
  int kk;                  /* Order of curve.                            */
  int ki;                  /* Counter                                    */
  int kdim = 1;            /* Dimension of curve in curve/point intersection.*/
  double *st;              /* Pointer to knotvector of curve.              */
  double *scoef;           /* Pointer to vertices of curve.                */
  double *sc = SISL_NULL;       /* Pointer to vertices of curve in curve/point
			      intersection.                                */
  double spoint[1];        /* SISLPoint in curve/point intersection.           */
  double *spar = SISL_NULL;     /* Values of intersections in the parameter 
			      area of the second object. Empty in this case. */
  double *sdiff = SISL_NULL;    /* Difference between point and coefficient.    */
  double *s1,*s2;          /* Pointers used to traverse coefficient arrays.*/
  SISLCurve *qc = SISL_NULL;        /* Pointer to curve in 
			      curve/point intersection.  */
  SISLPoint *qp = SISL_NULL;        /* Pointer to point in 
			      curve/point intersection.  */
  SISLObject *qo1 = SISL_NULL;      /* Pointer to curve in 
			      object/point intersection. */
  SISLObject *qo2 = SISL_NULL;      /* Pointer to point in 
			      object/point intersection    */
  SISLIntdat *qintdat = SISL_NULL;  /* Intersection result */
  double *rscoef;          /* Scaled coefficients if pc1 is rational        */
  double *s3;              /* pointer to weight if curve rational           */
  double wmin,wmax;        /* min and max values of the weights if rational */
  double scale;            /* factor for scaling weights if rational        */
  int i;                   /* loop variable                                 */
  int idimp1;              /* idim+1                                        */
  int      ksurf=0;         /* Dummy number of Intsurfs. */
  SISLIntsurf **wsurf=SISL_NULL;    /* Dummy array of Intsurfs. */
  int      kdeg=2000;       /* input to int_join_per. */
  SISLObject *track_obj=SISL_NULL;
  
  /*
  * Create new object and connect curve to object.
  * ------------------------------------------------
  */
  
  if (!(track_obj = newObject (SISLCURVE)))
    goto err101;
  track_obj->c1 = pc1;
  

  /* 
   * Check dimension.  
   * ----------------
   */

  *jpt  = 0;
  *jcrv = 0;
  *jtrack = 0;

  if ( !(idim == 2 || idim == 3)) goto err104;
  if (pc1 -> idim != idim) goto err103;

  /* 
   * Copy curve to local variables.  
   * ------------------------------
   */

  kn    = pc1 -> in;
  kk    = pc1 -> ik;
  st    = pc1 -> et;
  ikind = pc1 -> ikind;
  
  /* rational curves are a special case */
  if(ikind == 2 || ikind == 4)
  {
      /* scale the coeffs so that min. weight * max. weight = 1  */
      idimp1=idim+1;
      rscoef = pc1 -> rcoef;
      wmin=rscoef[idim];
      wmax=rscoef[idim];
      for(i=idim; i< kn*idimp1; i+=idimp1)
      {
          if(rscoef[i] < wmin) wmin=rscoef[i];
          if(rscoef[i] > wmax) wmax=rscoef[i];
      } 
      scale=1.0/sqrt(wmin*wmax);
      scoef=newarray(kn*idimp1,DOUBLE);
      for(i=0; i< kn*idimp1; i++)
      {
          scoef[i]=rscoef[i]*scale;
      } 
  }
  else
  {
      scoef = pc1 -> ecoef;
  }


  /* 
   * Allocate space for coeffecients of new curve and help array. 
   * ------------------------------------------------------------
   */
                         
  if ((sdiff = newarray(idim,double)) == SISL_NULL) goto err101;
  if ((sc = newarray(kn,double)) == SISL_NULL) goto err101;

  /* 
   * Put vertices into the equation of the plane/line. 
   * Compute new vertices.
   * -------------------------------------------------
   */

  for (s1=scoef,s2=sc,ki=0; ki<kn; ki++,s1+=idim,s2++)
    {
      if(ikind == 2 || ikind == 4)
      {
      /* curve is rational so we're using idim+1 - d homogeneous coords */
          s3=s1+idim;
          for(i=0; i<idim; i++)
          {
              sdiff[i]=(*(s1+i)) - (*s3) * epoint[i];
          }
          *s2 = s6scpr(sdiff,enorm,idim);
	  s1++;
      }
      else
      {
      /* curve is not rational so we're using ordinary idim - d coords */
          s6diff(s1,epoint,idim,sdiff);
          *s2 = s6scpr(sdiff,enorm,idim);
      }
    }                                 

  if(ikind == 2 || ikind == 4) freearray(scoef);

  /* 
   * Create new curve.  
   * -----------------
   */

  if(ikind == 2 || ikind == 4) ikind--;
  qc = newCurve(kn,kk,st,sc,ikind,kdim,0);
  if (qc == SISL_NULL) goto err101;
  qc->cuopen = pc1->cuopen;		  
      
  /* 
   * Create new object and connect curve to object.  
   * ----------------------------------------------
   */

  if (!(qo1 = newObject(SISLCURVE))) goto err101;
  qo1 -> c1 = qc;
  qo1 -> o1 = qo1;
  
  /*
   * Create new object and connect point to object.
   * ----------------------------------------------
   */

  if(!(qo2 = newObject(SISLPOINT))) goto err101;
  spoint[0] = DZERO;
  if(!(qp = newPoint(spoint,kdim,1))) goto err101;
  qo2 -> p1 = qp;

  /* 
   * Find intersections.  
   * -------------------
   */

  sh1761(qo1,qo2,aepsge,&qintdat,&kstat);
  if (kstat < 0) goto error;
  
  /* Join periodic curves */
  int_join_per( &qintdat,track_obj,track_obj,nullp,kdeg,aepsge,&kstat);
  if (kstat < 0)
    goto error;

  /* Create tracks */
  if (trackflag && qintdat)
    {
      make_tracks (qo1, qo2, 0, nullp,
		   qintdat->ilist, qintdat->vlist, 
		   jtrack, wtrack, aepsge, &kstat);
      if (kstat < 0)
	goto error;
    }

  /* 
   * Express intersections on output format.  
   * ---------------------------------------
   */

  if (qintdat)/* Only if there were intersections found */
    {
      hp_s1880(track_obj, track_obj, kdeg,
	       1,0,qintdat,jpt,gpar,&spar,pretop,jcrv,wcurve,&ksurf,&wsurf,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* 
   * Intersections found.  
   * --------------------
   */

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

 err101: *jstat = -101;
        s6err("sh1850",*jstat,kpos);
        goto out;

  /* Dimensions conflicting.  */

 err103: *jstat = -103;
        s6err("sh1850",*jstat,kpos);
        goto out;

  /* Dimension not equal to two or three.  */

 err104: *jstat = -104;                          
        s6err("sh1850",*jstat,kpos);
        goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
        s6err("sh1850",*jstat,kpos);
        goto out;

 out:

  /* Free allocated space.  */

  if (sc)      freearray(sc);
  if (sdiff)   freearray(sdiff);       
  if (spar)    freearray(spar);
  if (qo1)     freeObject(qo1);
  if (qo2)     freeObject(qo2);
  if (qintdat) freeIntdat(qintdat);
  if (track_obj)
    {
       track_obj->c1 = SISL_NULL;
       freeObject(track_obj);
    }

return;
}                                               
                                           
                       

