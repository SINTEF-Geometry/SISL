/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl-copyright.h"

#define SH1795

#include "sislP.h"


/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static int sh1795_s9iterate(SISLSurf *ps1, SISLSurf *ps2, double aepsge,
			    double eps, int dir, int sgn, 
			    double epar1[2], double epar2[2],
			    double initdist, double der2dist, double lim1,
			    double lim2, double *zonepar, int *jstat);
#else
  static int sh1795_s9iterate();
#endif

#if defined(SISLNEEDPROTOTYPES)
void sh1795(SISLObject *po1, SISLObject *po2, SISLIntpt *pt,
	    double aepsge, int *jstat)
#else
  void sh1795(po1, po2, pt, aepsge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt *pt;
     double aepsge; 
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a singular intersection point in a corner is 
*              isolated. In that case, define segmentation parameters 
*              bounding the singularity.
*
*
*
* INPUT      : po1    - First surface in intersection.
*              po2    - Second surface in intersection
*              pt     - Intersection point
*	       aepsge - Geometry resolution.
*                                                                     
*
* OUTPUT     : jstat  - status messages  
*                            = 2      : Entire surface withing tangential zone
*                            = 1      : Tangential zone found
*			     = 0      : Not an isolated singularity or tangential
*                                       not found
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
* WRITTEN BY : Vibeke Skytt, SINTEF, 2023-03
*
*********************************************************************
*/                                     
{
  int kstat = 0;   /* Local status variable                        */
  int kdim;
  int ki, kj, kr;
  double sdist[2];  /* Specified distances between surfaces         */
  double zonepar[12];  /* Parameter values corresponding to distances */
  double samplepar[6];  /* Parameter values at samples         */
  double tsize1, tsize2;
  int iobj;
  double *st1 = po1->s1->et1;
  double *st2 = po1->s1->et2;
  double *st3 = po2->s1->et1;
  double *st4 = po2->s1->et2;
  int kn1 = po1->s1->in1;
  int kn2 = po1->s1->in2;
  int kn3 = po2->s1->in1;
  int kn4 = po2->s1->in2;
  int kk1 = po1->s1->ik1;
  int kk2 = po1->s1->ik2;
  int kk3 = po2->s1->ik1;
  int kk4 = po2->s1->ik2;
  SISLSurf *qs1;   /* Surface with singular corner intersection */
  SISLSurf *qs2;   /* The other surface                            */
  double tlim1[2], tlim2[2];
  int kleft1=0, kleft2=0, kleft3=0, kleft4=0;  /* Indices in knot vector  */
  double spar1[2], spar2[2];  /* Parameter value of point on intersection */
  double scratch[69];  /* Storage for results of surface evaluation */
  double *sder1, *sder2, *snorm1, *snorm2;  /* Pointers to position and 
					       surface normals      */
  double *cvder1, *cvder2;
  double tang;    /* Angle between surface normals                  */
  int sgn[2];
  int converged = 0;
  double initdist;
  double tdum;
  double sstart[2], send[2];
  double ptol = 100*REL_PAR_RES;
  double tmin, tmax;
  double tpar;
  SISLCurve *qc = NULL;
  SISLObject *qoc = NULL;
  SISLIntdat *qintdat = NULL;
  double seg[2];
 
  sdist[0] = 1.2*aepsge;
  sdist[1] = 2.0*aepsge;

  /* Test configuration */
  if (po1->s1->idim != 3)
    goto err104;
  kdim = po1->s1->idim;
  
  sh6isinside(po1, po2, pt, &kstat);
  if (kstat < 0)
    goto error;
  if (kstat != 3 && kstat != 4)
    goto out;   /* Not a corner */

  /* Set pointers to evaluation results */
  sder1 = scratch;
  snorm1 = sder1 + 18;
  sder2 = snorm1 + 3;
  snorm2 = sder2 + 18;
  cvder1 = snorm2 + 3;
  cvder2 = cvder1 + 9;

  /* Identify corner and governing surface */
  tsize1 = (st1[kn1]-st1[kk1-1])*(st2[kn2]-st2[kk2-1]);
  tsize2 = (st3[kn3]-st3[kk3-1])*(st4[kn4]-st4[kk4-1]);
  iobj = (tsize1 >= tsize2) ? 0 : 1;
  if (iobj == 0 &&
      (!(DEQUAL(pt->epar[0],st1[kk1-1]) || DEQUAL(pt->epar[0],st1[kn1]))) &&
      (!(DEQUAL(pt->epar[1],st2[kk2-1]) || DEQUAL(pt->epar[1],st2[kn2]))))
    iobj = 1;
  if (iobj == 1 &&
      (!(DEQUAL(pt->epar[2],st3[kk3-1]) || DEQUAL(pt->epar[2],st3[kn3]))) &&
      (!(DEQUAL(pt->epar[3],st4[kk4-1]) || DEQUAL(pt->epar[3],st4[kn4]))))
    iobj = 0;
  qs1 = (iobj == 0) ? po1->s1 : po2->s1;
  qs2 = (iobj == 0) ? po2->s1 : po1->s1;
  samplepar[0] = pt->epar[2*iobj];
  samplepar[3] = pt->epar[2*iobj+1];
  tlim1[0] = qs1->et2[qs1->ik2-1];
  tlim2[0] = qs1->et2[qs1->in2];
  tlim1[1] = qs1->et1[qs1->ik1-1];
  tlim2[1] = qs1->et1[qs1->in1];
  sgn[0] = (DEQUAL(samplepar[0], qs1->et1[qs1->ik1-1])) ? 1 : -1;
  sgn[1] = (DEQUAL(samplepar[3], qs1->et2[qs1->ik2-1])) ? 1 : -1;
  spar1[0] = pt->epar[2*iobj];
  spar1[1] = pt->epar[2*iobj+1];
  spar2[0] = pt->epar[2*(1-iobj)];
  spar2[1] = pt->epar[2*(1-iobj)+1];
  sstart[0] = qs2->et1[qs2->ik1-1];
  sstart[1] = qs2->et2[qs2->ik2-1];
  send[0] = qs2->et1[qs2->in1];
  send[1] = qs2->et2[qs2->in2];

 
  s1421(qs1, 2, spar1, &kleft1, &kleft2, sder1, snorm1, &kstat);
  if (kstat < 0)
    goto error;
  s1421(qs2, 2, spar2, &kleft3, &kleft4, sder2, snorm2, &kstat);
  if (kstat < 0)
    goto error;
  tang = s6ang(snorm1, snorm2, 3);
  if (tang > ANGULAR_TOLERANCE && M_PI-tang > ANGULAR_TOLERANCE)
    goto out; /* Not a singularity */

  /* Corner */
  for (ki=0; ki<2; ++ki)
    {
      memmove(cvder1, sder1, kdim*sizeof(double));
      memmove(cvder1+kdim, sder1+(2-ki)*kdim, kdim*sizeof(double));
      memmove(cvder1+2*kdim, sder1+(5-2*ki)*kdim, kdim*sizeof(double));
      s1291(cvder1, sder2, kdim, cvder2, &kstat);
      if (kstat < 0)
	goto error;
      
      /* Compute the parameter value where two curves have the
	 specified distances (approximation) */
      initdist = s6dist(cvder1, cvder2, kdim);
      tdum = s6dist(cvder1+2*kdim, cvder2+2*kdim, kdim);

      for (kj=0; kj<2; ++kj)
	{
	  converged = sh1795_s9iterate(qs1, qs2, aepsge, sdist[kj], 1-ki, sgn[1-ki],
				       spar1, spar2, initdist, 
				       tdum, tlim1[1-ki], tlim2[1-ki],
				       zonepar+6*ki+kj, &kstat);
	  if (kstat < 0)
	    goto error;
	  if (!converged)
	    goto out;
	}

    }

  /* Additional samples */
  samplepar[2] = samplepar[0] + 0.9*sgn[0]*zonepar[6];
  samplepar[1] = 0.5*(samplepar[0] + samplepar[2]);
  samplepar[5] = samplepar[3] + 0.9*sgn[1]*zonepar[0];
  samplepar[4] = 0.5*(samplepar[3] + samplepar[5]);
  for (ki=0; ki<2; ++ki)
    {
      spar1[1-ki] = pt->epar[2*iobj+1-ki];
      spar2[0] = pt->epar[2*(1-iobj)];
      spar2[1] = pt->epar[2*(1-iobj)+1];
      for (kr=1; kr<=2; ++kr)
	{
	  spar1[ki] = samplepar[3*ki+kr];

	  /* Evaluate first surface in the start point of the search */
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

	  /* Compute 0-2. derivative of corresponding curves in surfaces */
	  /* Surface with boundary intersection */
	  memmove(cvder1, sder1, kdim*sizeof(double));
	  memmove(cvder1+kdim, sder1+(2-ki)*kdim, kdim*sizeof(double));
	  memmove(cvder1+2*kdim, sder1+(5-2*ki)*kdim, kdim*sizeof(double));

	  /* The other surface */
	  s1291(cvder1, sder2, kdim, cvder2, &kstat);
	  if (kstat < 0)
	    goto error;

	  /* Compute the parameter value where two curves have the
	     specified distances (approximation) */
	  initdist = s6dist(cvder1, cvder2, kdim);
	  tdum = s6dist(cvder1+2*kdim, cvder2+2*kdim, kdim);

	  for (kj=0; kj<2; ++kj)
	    {
	      if (initdist > sdist[kj])
		zonepar[6*ki+2*kr+kj] = 0.0;
	      else
		converged = sh1795_s9iterate(qs1, qs2, aepsge, sdist[kj],
					     1-ki, sgn[1-ki],
					     spar1, spar2, initdist, 
					     tdum, tlim1[1-ki], tlim2[1-ki],
					     zonepar+6*ki+2*kr+kj, &kstat);
	      if (kstat < 0)
		goto error;
	      if (!converged)
		goto out;
	    }
	}
    }

  for (ki=0; ki<2; ++ki)
    {
      /* Define test curve */
      tmin = zonepar[6*ki];
      tmax = zonepar[6*ki+1];
      for (kr=1; kr<3; ++kr)
	{
	  tmin = max(tmin, zonepar[6*ki+2*kr]);
	  tmax = min(tmax, zonepar[6*ki+2*kr+1]);
	}
      tpar = pt->epar[2*iobj+1-ki] + sgn[1-ki]*max(0.5*(tmin+tmax),zonepar[6*ki]);
      
      /* Check with possible help points */
      for (kj=0; kj<pt->no_of_curves; ++kj)
	{
	  if (sh6ishelp(pt->pnext[kj]))
	    {
	      sh6isinside(po1, po2, pt->pnext[kj], &kstat);
	      if (kstat < 0)
		goto error;
	      if (kstat == 0)
		continue;
	      if (fabs(pt->pnext[kj]->epar[2*iobj+1-ki]-pt->epar[2*iobj+1-ki]) >
		  fabs(tpar - pt->epar[2*iobj+1-ki]))
		tpar = pt->pnext[kj]->epar[2*iobj+1-ki] + sgn[1-ki]*ptol;
	    }
	}
	  
      /* Check for intersections along limiting curve */
      if (ki == 1)
	s1437(qs1, tpar, &qc, &kstat);
      else
	s1436(qs1, tpar, &qc, &kstat);
      if (kstat < 0)
	goto error;

      /* Intersect with the other surface and count intersections */
      if ((qoc = newObject(SISLCURVE)) == NULL)
	goto err101;
      qoc->c1 = qc;

      sh1761((iobj == 0) ? po2 : po1, qoc, aepsge, &qintdat, &kstat);
      if (kstat < 0)
	goto error;
      if (qintdat)
	goto out;  /* Not an isolated singularity */
      seg[ki] = tpar;
    }

  /* Define segmentation parameters */
  for (ki=0; ki<2; ++ki)
    {
      sh6setseg(qs1, 1-ki, seg+ki, 1, (sgn[1-ki] == 1) ? TANGENTIAL_BELT_LEFT :
		TANGENTIAL_BELT_RIGHT, &kstat);
      if (kstat < 0)
	goto error;
   }

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
  if (qoc != NULL) 
    {
      freeObject(qoc);
      qc = NULL;
    }
  if (qc != NULL) freeCurve(qc);
  if (qintdat != NULL) freeIntdat(qintdat);
  return;
}


#if defined(SISLNEEDPROTOTYPES)
static int sh1795_s9iterate(SISLSurf *ps1, SISLSurf *ps2, double aepsge,
			    double eps, int dir, int sgn, 
			    double epar1[2], double epar2[2],
			    double initdist, double der2dist, double lim1,
			    double lim2, double *zonepar, int *jstat)
#else
  static int sh1795_s9iterate(ps1, ps2, aepsge, eps, dir, sgn, epar1,
			      epar2, initdist, der2dist, lim1, lim2,
			      zonepar, jstat)
     SISLSurf *ps1; 
     SISLSurf *ps2; 
double aepsge;
double eps;
int dir;
int sgn; 
double epar1[2];
double epar2[2];
double initdist;
double der2dist;
double lim1;
double lim2;
double *zonepar;
int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : 
*
*
*
* INPUT      : ps1    - First surface in intersection.
*              ps2    - Second surface in intersection
*              
*                                                                     
*
* OUTPUT     : 
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
* WRITTEN BY : Vibeke Skytt, SINTEF, 2023-03
*
*********************************************************************
*/                                     
{
  int kstat = 0;
  int kdim = ps1->idim;
  int kmax = 7;    /* Maximum number of iterations to find par. val.*/
  double tprev = 0.0;
  double tdistprev = HUGE;
  double dist;
  double tdel;
  int failure = 0;
  int kr;
  double par1[2], par2[2], pos1[3], pos2[3];
  int kleft1 = 0, kleft2 = 0, kleft3 = 0, kleft4 = 0;
  double sstart[2], send[2];
  double dtol = 0.1*aepsge;/* Tolerance in search for zone par. val.*/
  double tdist0;
  double tdel2;
  int converge = 0;

  sstart[0] = ps2->et1[ps2->ik1-1];
  sstart[1] = ps2->et2[ps2->ik2-1];
  send[0] = ps2->et1[ps2->in1];
  send[1] = ps2->et2[ps2->in2];

  if (der2dist > REL_PAR_RES)
    {
      tdel = 2*eps/der2dist;
      tdel = sqrt(tdel);
      if (tdel > lim2 -lim1)
	tdel = lim2 - lim1;
    }
  else
    tdel = lim2 - lim1;
  
  /* Iterate to a better position */
  for (kr=0; kr<kmax; ++kr)
    {
      /* Check accuracy of current estimate */
      par1[1-dir] = epar1[1-dir];
      par1[dir] = epar1[dir] + sgn*tdel;
      par2[0] = epar2[0];
      par2[1] = epar2[1];
      s1424(ps1, 0, 0, par1, &kleft1, &kleft2, pos1, &kstat);
      if (kstat < 0)
	goto error;
      s1775(ps2, pos1, kdim, aepsge, sstart, send, par2, par2, &kstat);
      if (kstat < 0)
	goto error;
      s1424(ps2, 0, 0, par2, &kleft3, &kleft4, pos2, &kstat);
      if (kstat < 0)
	goto error;
      dist = s6dist(pos1, pos2, kdim);
      if (fabs(dist - eps) < dtol ||
	  dist < eps && tdel >= lim2 - lim1 - REL_PAR_RES)
	{
	  tprev = tdel;
	  break;   /* The estimate for the tangential zone limit
		      is close enough                            */
	}

      if (fabs(dist - eps) > fabs(tdistprev - eps))
	{
	  tdel = 0.5*(tprev+tdel);  // We are loosing accuracy
	  continue;
	}
      
      /* Find new position */
      tdist0 = (kr == 0) ? initdist : tdistprev;
      tdel2 = (eps-dist)*(tdel-tprev)/(dist-tdist0);
      if (DEQUAL(tdel2, 0.0))
	break;
      tprev = tdel;
      tdistprev = dist;
      tdel += tdel2;
    }
  
  if (fabs(dist - eps) < dtol ||
      dist < eps && tdel >= lim2 - lim1 - REL_PAR_RES)
    {
      *zonepar = tprev;
      converge = 1;
    }
  else
      converge = 0;

  goto out;

error:
  /* Error in lower order function */
  *jstat = kstat;
  goto out;

 out:
  return converge;
}
