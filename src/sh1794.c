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

#define SH1794

#include "sislP.h"


/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static int sh1794_s9findedg(SISLSurf *ps1, SISLSurf *ps2, SISLIntpt **up,
			    int nmb_pt, SISLIntpt **pt1, SISLIntpt **pt2,
			    int *jstat);
#else
static int sh1794_s9findedg();
#endif

#if defined(SISLNEEDPROTOTYPES)
void sh1794(SISLObject *po1, SISLObject *po2, SISLIntpt **up, int nmb_pt,
	    double aepsge, int *jstat)
#else
  void sh1794(po1, po2, up, nmb_pt, aepsge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt **up;
     int nmb_pt; 
     double aepsge; 
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if an intersection curve following a constant
*              parameter direction between two specified intersection
*              points represent a tangential intersection curve
*              between two corners in the intersection problem
*
*
*
* INPUT      : po1    - First surface in intersection.
*              po2    - Second surface in intersection
*              pt1    - Intersection point in the start of the curve
*              pt2    - Intersection point in the end of the curve
*              idir   - Constant parameter direction (0 and 1 corresponds
*	       aepsge - Geometry resolution.
*                                                                     
*
* OUTPUT     : jstat  - status messages  
*                            = 1      : Tangential curve found
*			     = 0      : The curve is not tangential or
*			                do not end in corners
*                            < 0      : error
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
* WRITTEN BY : Vibeke Skytt, SINTEF, 2018-02
*
*********************************************************************
*/                                     
{
  int kstat = 0;   /* Local status variable                        */
  int ki, kj, kr;  /* Counters                                     */
  int kmax = 7;    /* Maximum number of iterations to find par. val.*/
  int kdim;        /* Dimension of geometry space                  */ 
  SISLSurf *qs1;   /* Surface with constant parameter intersection */
  SISLSurf *qs2;   /* The other surface                            */
  int nsample;     /* Number of sampling points                    */
  int kleft1=0, kleft2=0, kleft3=0, kleft4=0;  /* Indices in knot vector  */
  double spar1[2], spar2[2];  /* Parameter value of point on intersection
				 curve with respect to surface            */
  double tstart, tend;  /* Endparameters of intersection curve in surface */
  double tdel, tdel2;   /* Step length                                  */
  int kdir;        /* Parameter direction of constant parameter in surface */
  int kd;          /* Parameter direction                          */
  double scratch[60];  /* Storage for results of surface evaluation */
  double *sder1, *sder2, *snorm1, *snorm2;  /* Pointers to position and 
					       surface normals      */
  double *cvder1, *cvder2;
  double tang;    /* Angle between surface normals                  */
  double sstart[2], send[2];  /* Parameter limits of surface        */
  int kk, kn;     /* Order and number of coefficients of boundary 
		     intersection curve                             */
  double *st;     /* Knot vector along boundary intersection curve  */
  double angtol = 2.0*ANGULAR_TOLERANCE;  /* Needs tuning           */
  double sdist[2];  /* Specified distances between surfaces         */
  double dtol = 0.1*aepsge;/* Tolerance in search for zone par. val.*/
  double *zonepar = NULL;  /* Parameter values corresponding to distances */
  double *samplepar = NULL;  /* Parameter values at samples         */
  double *seg = NULL;   /* Segmentation parameters along intersection curve */
  int nseg = 0;         /* Number of segmentation parameters        */
  double ptol = 100*REL_PAR_RES;
  double tpar;    /* Parameter value at the end of the zone         */
  double tprev;   /* Previous estimate for parameter value          */
  double tdistprev;  /* Distance at previous parameter value        */
  double initdist;   /* Distance at initial parameter value        */
  double tdum;    /* Estimate for distance at given parameter value */

  double par1[2], par2[2], pos1[3], pos2[3];
  double tstart2, tend2, dist2;
  int sgn;
  double tmin, tmax;  /* Allowed endparameters of tangential zone   */
  double td1;         /* Distance between intersection points       */
  int failure;
  int idir;
  SISLIntpt *pt1, *pt2;

  sdist[0] = 1.2*aepsge;
  sdist[1] = 2.0*aepsge;

  /* Test input */
  if (po1->iobj != SISLSURFACE || po2->iobj != SISLSURFACE)
    goto warn0;  /* Not a surface-surface intersection */

  /* Identify edge intersection curve */
  idir = sh1794_s9findedg(po1->s1, po2->s1, up, nmb_pt, &pt1, &pt2, &kstat);
  if (kstat < 0) 
    goto error;
  if (idir < 0 || idir > 3)
    goto warn0;  /* No constant parameter direction specified */

  /* Initiate */
  *jstat = 0;

  /* Test for corner configuration */
  sh6isinside(po1, po2, pt1, &kstat);
  if (kstat < 0)
    goto error;
  if (kstat < 3)
    goto warn0;  /* Not a corner */

  sh6isinside(po1, po2, pt2, &kstat);
  if (kstat < 0)
    goto error;
  if (kstat < 3)
    goto warn0;  /* Not a corner */

  /* Evaluate the intersection curve in a number of sampling points
     to check whether it is tangential */
  qs1 = (idir <= 1) ? po1->s1 : po2->s1;
  qs2 = (idir <= 1) ? po2->s1 : po1->s1;
  kdim = qs1->idim;
  if (kdim != 3)
    goto err104;

  /* Prepare for checking */
  kdir = (idir <= 1) ? idir : idir-2;
  kd = (idir <= 1) ? 1-idir : 3-kdir;
  spar1[kdir] = pt1->epar[idir];
  tstart = pt1->epar[kd];
  tend =  pt2->epar[kd];
  if (tstart > tend)
    {
      tdum = tstart;
      tstart = tend;
      tend = tdum;
    }

  sstart[0] = qs2->et1[qs2->ik1-1];
  send[0] = qs2->et1[qs2->in1];
  sstart[1] = qs2->et2[qs2->ik2-1];
  send[1] = qs2->et2[qs2->in2];
  spar2[0] = (idir <= 1) ? pt1->epar[2] : pt1->epar[0];
  spar2[1] = (idir <= 1) ? pt1->epar[3] : pt1->epar[1];

  if (kdir == 0)
    {
      tstart2 = qs1->et1[qs1->ik1-1];
      tend2 = qs1->et1[qs1->in1];
    }
  else
    {
      tstart2 = qs1->et2[qs1->ik2-1];
      tend2 = qs1->et1[qs2->in2];
    }
  sgn = (spar1[kdir]-tstart < tend-spar1[kdir]) ? 1 : -1;

  /* Set pointers to evaluation results */
  sder1 = scratch;
  snorm1 = sder1 + 18;
  sder2 = snorm1 + 3;
  snorm2 = sder2 + 18;
  cvder1 = snorm2 + 3;
  cvder2 = cvder1 + 9;

  /* Compute number of sampling points */
  kn = (kdir == 0) ? qs1->in2 : qs1->in1;
  kk = (kdir == 0) ? qs1->ik2 : qs1->ik1;
  st = (kdir == 0) ? qs1->et2 : qs1->et1;
  s1219(st, kk, kn, &kleft1, tstart, &kstat);
  if (kstat < 0)
    goto error;
  s1219(st, kk, kn, &kleft2, tend, &kstat);
  if (kstat < 0)
    goto error;
  nsample = 3*(kleft2 - kleft1)*kk;
  nsample = max(3*kk, min(nsample, 20*kk));
  tdel = (tend - tstart)/(double)(nsample-1);

  /* Initial check to avoid overemphasising tiny curves */
  s1421(qs1, 0, pt1->epar+2*(idir>0), &kleft1, &kleft2, sder1, snorm1, &kstat);
  if (kstat < 0)
    goto error;
  s1421(qs1, 0, pt2->epar+2*(idir>0), &kleft1, &kleft2, sder2, snorm2, &kstat);
  if (kstat < 0)
    goto error;
  td1 = s6dist(sder1, sder2, kdim);
  if (td1 < aepsge && (DNEQUAL(tstart,st[kk-1]) || DNEQUAL(tend,st[kn])))
    goto warn0;

  nsample = max(2, min(nsample, (int)(td1/aepsge)));

  zonepar = new0array(2*nsample, DOUBLE);
  samplepar = newarray(nsample, DOUBLE);
  seg = newarray(2+nsample, DOUBLE);
  if (zonepar == NULL || samplepar == NULL || seg == NULL) 
    goto err101;

  for (ki=0, spar1[1-kdir]=tstart; ki<nsample; ++ki, spar1[1-kdir]+=tdel)
    {
      samplepar[ki] = spar1[1-kdir];

      /* Test if the two surfaces intersect tangentially in the
	 current sampling point.
	 First evaluate the surface in which there is an intersection
	 curve along the boundary                                     */
      s1421(qs1, 2, spar1, &kleft1, &kleft2, sder1, snorm1, &kstat);
      if (kstat < 0)
	goto error;

      /* Find the closest point in the other surface. */
      s1775(qs2, sder1, kdim, aepsge, sstart, send, spar2, spar2, &kstat);
      if (kstat < 0)
	goto error;

      /* Evaluate the second surface. */
      s1421(qs2, 2, spar2, &kleft3, &kleft4, sder2, snorm2, &kstat);
      if (kstat < 0)
	goto error;

      // Test if the current point is a touch point
      tang = s6ang(snorm1, snorm2, kdim);
      if (tang > angtol)
	break;

      /* Compute 0-2. derivative of corresponding curves in surfaces */
      /* Surface with boundary intersection */
      memmove(cvder1, sder1, kdim*sizeof(double));
      memmove(cvder1+kdim, sder1+(1+kdir)*kdim, kdim*sizeof(double));
      memmove(cvder1+2*kdim, sder1+(3+2*kdir)*kdim, kdim*sizeof(double));

      /* The other surface */
      s1291(cvder1, sder2, kdim, cvder2, &kstat);
      if (kstat < 0)
	goto error;

      /* Compute the parameter value where two curves have the
	 specified distances (approximation) */
      tprev = 0.0;
      initdist = tdistprev = s6dist(cvder1, cvder2, kdim);
      tdum = s6dist(cvder1+2*kdim, cvder2+2*kdim, kdim);

      for (kj=0; kj<2; ++kj)
	{
	  if (tdum > REL_PAR_RES)
	    {
	      tpar = 2*sdist[kj]/tdum;
	      tpar = sqrt(tpar);
	    }
	  else
	    tpar = 0.0;

	  /* Iterate to a better position */
	  /* Could also take the parameter value of the previous
	     sample point as input  */
	  failure = 0;
	  for (kr=0; kr<kmax; ++kr)
	    {
	      /* Check accuracy of current estimate */
	      par1[1-kdir] = spar1[1-kdir];
	      par1[kdir] = spar1[kdir] + sgn*tpar;
	      par2[0] = spar2[0];
	      par2[1] = spar2[1];
	      s1424(qs1, 0, 0, par1, &kleft1, &kleft2, pos1, &kstat);
	      if (kstat < 0)
		goto error;
	      s1775(qs2, pos1, kdim, aepsge, sstart, send, par2, par2, &kstat);
	      if (kstat < 0)
		goto error;
	      s1424(qs2, 0, 0, par2, &kleft3, &kleft4, pos2, &kstat);
	      if (kstat < 0)
		goto error;
	      dist2 = s6dist(pos1, pos2, kdim);
	      /*printf("Dist: %7.13f, expected: %7.13f \n",dist2, sdist[kj]); */
	      if (fabs(dist2 - sdist[kj]) < dtol)
		break;   /* The estimate for the tangential zone limit
			    is close enough                            */
	      if (DEQUAL(dist2, tdistprev) || dist2 < initdist)
		{
		  failure = 1;
		  break;
		}
	      if (fabs(dist2 - sdist[kj]) > fabs(tdistprev - sdist[kj]))
		{
		  tpar = 0.5*(tprev+tpar);  // We are loosing accuracy
		  continue;
		}

	      /* Find new position */
	      tdel2 = (sdist[kj]-dist2)*(tpar-tprev)/(dist2-tdistprev);
	      tprev = tpar;
	      tdistprev = dist2;
	      tpar += tdel2;
	    }
	  if (failure)
	    break;
	  zonepar[2*ki+kj] = tprev;
	}

      if (failure)
	break;

      if (zonepar[2*ki] > zonepar[2*ki+1])
	{
	  tdum = zonepar[2*ki];
	  zonepar[2*ki] = zonepar[2*ki+1];
	  zonepar[2*ki+1] = tdum;
	}
    }
  
  if (ki < nsample)
    {
      /* Not a tangential intersection curve */
      *jstat = 0;
      goto out;
    }

  /* Identify non-overlapping pairs of tangential zone end parameters */
  /* Check also whether the tangential intersection curve follows the
     entire surface boundary  */
  if (tstart > st[kk-1] + ptol)
    seg[nseg++] = tstart;

  for (ki=0; ki<nsample; ki=kj)
    {
      for (kj=ki+1; kj<nsample; ++kj)
	{
	  if (zonepar[2*ki] > zonepar[2*kj+1] ||
	      zonepar[2*ki+1] < zonepar[2*kj])
	    {
	      /* No overlap. Define segmentation parameter */
	      /* First look for a suitable knot value */
	      s1219(st, kk, kn, &kleft1, samplepar[ki], &kstat);
	      if (kstat < 0)
		goto error;
	      s1219(st, kk, kn, &kleft2, samplepar[kj], &kstat);
	      if (kstat < 0)
		goto error;
	      if (samplepar[kj] > st[kleft2]+ptol && kleft2 > kleft1)
		{
		  seg[nseg++] = st[kleft2];
		  for (kr=kj-1; kr>ki && samplepar[kr]>st[kleft2]; --kr);
		  kj = max(kr, ki+1);
		}
	      else if (kj > ki+1)
		{
		  seg[nseg++] = samplepar[kj-1];
		  --kj;
		}
	      else
		seg[nseg++] = 0.5*(samplepar[ki]+samplepar[kj]);

	      break;
	    }
	}
    }

  if (tend < st[kn] - ptol)
    seg[nseg++] = tend;
  
  if (nseg > 0)
    {
      /* Define segmentation parameters along the tangential 
	 intersection curve   */
      sh6setseg(qs1, 1-kdir, seg, nseg, LIMITING_SEG, &kstat);
      if (kstat < 0)
	goto error;
    }
  else
    {
      /* Set width of tangential zone */
      tmin = zonepar[0];
      tmax = zonepar[1];
      for (ki=1; ki<nsample; ++ki)
	{
	  tmin = max(tmin, zonepar[2*ki]);
	  tmax = min(tmax, zonepar[2*ki+1]);
	}
      tpar = spar1[kdir] + 0.5*sgn*(tmin + tmax);

      /* Define segmentation parameter perpendicular to the tangential
	 intersection curve */
      seg[nseg++] = tpar;
      sh6setseg(qs1, kdir, seg, nseg, (sgn == 1) ? TANGENTIAL_BELT_LEFT :
		TANGENTIAL_BELT_RIGHT, &kstat);
      if (kstat < 0)
	goto error;
    }

  *jstat = 1;
  goto out;

 warn0:
  /* Not a tangential curve */
  *jstat = 0;
  goto out;

error:
  /* Error in lower order function */
  *jstat = kstat;
  goto out;

 err101:
  /* Error in scratch allocation */
  *jstat = -101;
  goto out;

 err104:
  /* Dimension of geometry space not equal to 3 */
  *jstat = -104;
  goto out;

 out:
  if (zonepar != NULL) freearray(zonepar);
  if (samplepar != NULL) freearray(samplepar);
  if (seg != NULL) freearray(seg);

  return;
}


#if defined(SISLNEEDPROTOTYPES)
static int sh1794_s9findedg(SISLSurf *ps1, SISLSurf *ps2, SISLIntpt **up,
			    int nmb_pt, SISLIntpt **pt1, SISLIntpt **pt2,
			    int *jstat)
#else
static int sh1794_s9findedg(ps1, ps2, up, nmb_pt, pt1, pt2, jstat)
SISLSurf *ps1; 
SISLSurf *ps2; 
SISLIntpt **up;
int nmb_pt; 
SISLIntpt **pt1; 
SISLIntpt **pt2;
int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a chain of linked intersection points
*              contain an intersection curve along some surface edge
*              and select the most signicant candidate
*
*
*
* INPUT      : ps1    - First surface in intersection.
*              ps2    - Second surface in intersection
*              up     - Chain of intersection points
*	       nmb_pt - Number of intersection points
*                                                                     
*
* OUTPUT     : return value - Parameter direction of detected edge (both sfs)
*              pt1    - Intersection point in the start of the curve
*              pt2    - Intersection point in the end of the curve
*              jstat  - status messages  
*			     = 0      : OK
*                            < 0      : error
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
* WRITTEN BY : Vibeke Skytt, SINTEF, 2018-02
*
*********************************************************************
*/                                     
{
  int kstat = 0;  /* Local status variable                             */
  int kdir;  /* Parameter direction of possible tangential intersection curve */
  int kdir2;
  int pdir[4];    /* Maximum possible parameter directions             */
  int ix[8];      /* Maximum possible end parameters of constant chain */
  int ki, kj;     /* Counters                                          */
  int kidx = 0;   /* Current parameter direction                       */
  int klist1=0, klist2=0;
  double parlen1, parlen2;  /* Distance between intersection points in
			       parameter domain along constant direction */
  int p1, p2;               /* Parameter directions                      */
  double del1, del2;        /* Edge length of surface in parameter domain */

  pdir[0] = pdir[1] = pdir[2] = pdir[3] = -1;
  ix[2*kidx] = 0;
  for (kj=0; kj<nmb_pt; kj=max(kj+1,ki-1))
    {
      kdir = -1;
      for (ki=kj+1; ki<nmb_pt; ++ki)
	{
	  sh6getlist(up[ki-1], up[ki], &klist1, &klist2, &kstat);
	  if (kstat < 0)
	    goto error;
		       
	  for (kdir2=0; kdir2<up[ki-1]->ipar; ++kdir2)
	    {
	      if (up[ki-1]->curve_dir[klist1] & (1 << (kdir2+1)))
		break;
	    }
	  ix[2*kidx] = kj;
	  if (kdir2 == up[ki-1]->ipar)
	    {
	      if (kdir >= 0)
		{
		  pdir[kidx] = kdir;
		  ix[2*kidx+1] = ki-1;
		  kidx++;
		  kdir = 0;
		}
	      kdir2 = -1;
	      break;
	    }
	  else if (kdir >= 0 && kdir2 != kdir)
	    {
	      pdir[kidx] = kdir;
	      ix[2*kidx+1] = ki;
	      kidx++;
	      kdir = kdir2;
	      break;
	    }
	  else
	    {
	      kdir = kdir2;
	    }
	}
    }
  if (kdir2 >= 0 && (kidx == 0 || kdir2 != pdir[kidx-1]))
    {
      pdir[kidx] = kdir2;
      ix[2*kidx+1] = ki;
      kidx++;
    }
      
  kdir = -1;
  if (kidx > 0)
    {
      kdir = pdir[0];
      p1 = (kdir <= 1) ? 1 - kdir : 5 - kdir;
      parlen1 = fabs(up[ix[1]-1]->epar[p1]-up[ix[0]]->epar[p1]);
      if (kdir <= 0)
	del1 = (p1 == 0) ? ps1->et1[ps1->in1] - ps1->et1[ps1->ik1-1] :
	  ps1->et2[ps1->in2] - ps1->et2[ps1->ik2-1];
      else
	del1 = (p1 == 0) ? ps2->et1[ps2->in1] - ps2->et1[ps2->ik1-1] :
	  ps2->et2[ps2->in2] - ps2->et2[ps2->ik2-1];

      (*pt1) = up[ix[0]];
      (*pt2) = up[ix[1]-1];
      for (ki=1; ki<kidx; ++ki)
	{
	  /* Check if other edge intersections are more significant */
	  kdir2 = pdir[ki];
	  p2 = (kdir2 <= 1) ? 1 - kdir2 : 5 - kdir2;
	  parlen2 = fabs(up[ix[2*ki+1]-1]->epar[p2]-up[ix[2*ki]]->epar[p2]);
	  if (kdir2 <= 0)
	    del2 = (p2 == 0) ? ps1->et1[ps1->in1] - ps1->et1[ps1->ik1-1] :
	      ps1->et2[ps1->in2] - ps1->et2[ps1->ik2-1];
	  else
	    del2 = (p2 == 0) ? ps2->et1[ps2->in1] - ps2->et1[ps2->ik1-1] :
	      ps2->et2[ps2->in2] - ps2->et2[ps2->ik2-1];
	  if (parlen2/del2 > parlen1/del1)
	    {
	      kdir = kdir2;
	      (*pt1) = up[ix[2*ki]];
	      (*pt2) = up[ix[2*ki+1]-1];
	      parlen1 = parlen2;
	      del1 = del2;
	    }
	}
    }

  *jstat = 0;
  goto out;

 error:
  *jstat = kstat;
  goto out;

 out:
  return kdir;
}

