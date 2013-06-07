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
 * $Id: s1990.c,v 1.5 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1990

#include "sislP.h"
/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void s1990_s9edg(double [],double [],double [],double,double *,
			int,int *);
/*
static void s1990_s9smooth(double [],int,int,int,double,double [],int *);
*/
#else
static void s1990_s9edg();
/* static void s1990_s9smooth(); */
#endif

#if defined(SISLNEEDPROTOTYPES)
void
     s1990(SISLSurf *ps,double aepsge,int *jstat)
#else
void s1990(ps,aepsge,jstat)
     SISLSurf   *ps;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To make the orientation surface on the unit sphere to
*	       a b-spline surface, the surface is representated with
*	       a surrounding cone piced from the unit sphere.
*
*
*
* INPUT      : ps     - The original B-spline surface.
*              aepsge - Geometry resolution.
*
*
* OUTPUT     : jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : We are making a cone surrounding the orientating surface
*	       on the unit sphere. The cone is representated with senter
*	       coordinates and an angle. The orientation is computed
*	       from aproximation of the normal to the surface.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-01.
*              UJK, Changed to accept equality between vertices.
* REWISED BY : Vibeke Skytt, SI, 91-02.
*              UJK, SI, 91-10.Accepting surf's degenarated to a curve lying in
*                   a plane. Necessary for 2D !
* Revised by : Christophe Rene Birkeland, SINTEF OSLO, June 1993.
*
*********************************************************************
*/
{
  int kpos = 0;     /* Position of the error.                             */
  int kstat = 0;    /* Local status variable.                             */
  int kfirst = 1;   /* Flag to mark if the first patch is treating.       */
  int kcount;       /* Counts number of vanishing normals.                */
  int kn1;          /* Number of vertices of surface in 1. par. direction.*/
  int kn2;          /* Number of vertices of surface in 2. par. direction.*/
  int kdim;	   /* Dimension of the space in which the objects lie.   */
  int kdim4;	   /* Help variable to contain  4*kdim.			 */
  int kver,khor;    /* The index to the vertice in the upper left corner 
		       to the patch to treat.				 */
  int k1,k2,k3,k4;  /* Control variables in loop. 			 */
  int ki;           /* Control variable in loop.  			 */
  int lcone[4];     /* Flag telling if the cone has been generated.       */
  double *t=SISL_NULL;   /* Allocating t[5][kdim]. Five tangents around the
		       patch, the first and the last is the same.         */
  double *tn;       /* Allocating tn[4][kdim]. Four normals in the corner
		       of the patch.					 */
  double *tsen;     /* Allocating tsen[4][kdim] for senter in edge cones. */
  double *ttan;     /* Allocating ttan[kdim] for tangent on edges.        */
  double tmax,tmin; /* Maximum and minimum coordinates to the narmals in
		       the first patch.					 */
  double tlen;      /* The length of a vector.				 */
  double tnlen;     /* The length of a normal vector.	   	         */
  double tang;	   /* An angle between two vectors.			 */
  double t1,t2;     /* Help variables.					 */
  double sang[4];   /* Angel to the cones to edges.                       */
  double svec1[3];  /* Vectors used to determin degeneration.             */
  double svec2[3];  /* Vectors used to determin degeneration.             */
  double *scoef;    /* Pointer to smoothed coefficient vector.            */
  double slen[5];   /* Distances between coefficients.                    */
  double scorn[4];  /* Angle between derivatives in corner of patch.      */
  
  /* Initiate output status */

  *jstat = 0;
  
  /* Test if the surfaces already have been treated.  */
  
  if (ps->pdir != SISL_NULL) goto out;
  
  /* Initialate dimentions. */
  
  kdim = ps -> idim;
  kn1  = ps -> in1;
  kn2  = ps -> in2;
  kdim4 = 4*kdim;
  
  lcone[0] = 1;
  lcone[1] = 1;
  lcone[2] = 1;
  lcone[3] = 1;
    
  /*Make a new direction cone. */
  
  if ((ps->pdir = newdir(kdim)) == SISL_NULL) goto err101;
  
  ps->pdir->aang = DZERO;
  for (k1=0;k1<kdim;k1++) ps->pdir->ecoef[k1] = DZERO;
  
  /* Allocate scratch for smoothed coefficients.  */
  
  if ((ps->pdir->esmooth = newarray(kn1*kn2*kdim,DOUBLE)) == SISL_NULL) goto err101;
  scoef = ps->pdir->esmooth;
  
  /* Compute coefficients of smoothed curve.  */
  
  /* s1990_s9smooth(ps->ecoef,kn1,kn2,kdim,aepsge,scoef,&kstat);
  if (kstat < 0) goto error; */
  
  memcopy(scoef,ps->ecoef,kn1*kn2*kdim,DOUBLE); 
  
  /* Allocate local used matrices, t[5][kdim] and tn[4][kdim]. */
  
  if ((t = newarray(14*kdim,double)) == SISL_NULL) goto err101;
  tn   = t + 5*kdim;
  tsen = tn + 4*kdim;
  ttan = tsen + 4*kdim;
  
  /* Here we are treating each patch in the control polygon separately.*/
  
  for (kver=0; kver < (kn2-1); kver++)
    for (khor=0; khor < (kn1-1); khor++)
      {
	slen[0] = slen[1] = slen[2] = slen[3] = DZERO;
	scorn[0] = scorn[1] = scorn[2] = scorn[3] = DZERO;
	
	/* Here we make the tangents in each corner of the patch,
           and in direction with the clock. The first and the last
	   vector contains both the first tangent. */
	
	k2 = (kver*kn1+khor)*kdim;
	
	for (k1=0; k1 < kdim; k1++,k2++)
	  {
	    t[kdim+k1]   = scoef[k2+kdim] - scoef[k2];
	    t[2*kdim+k1] = scoef[k2+(kn1+1)*kdim]-scoef[k2+kdim];
	    t[3*kdim+k1] = scoef[k2+kn1*kdim]-scoef[k2+(kn1+1)*kdim];
	    t[kdim4+k1] = t[k1] = scoef[k2]-scoef[k2+kn1*kdim];
	    
	    slen[0] += t[k1]*t[k1];
	    slen[1] += t[k1+kdim]*t[k1+kdim];
	    slen[2] += t[k1+2*kdim]*t[k1+2*kdim];
	    slen[3] += t[k1+3*kdim]*t[k1+3*kdim];
	  }
	slen[4] = slen[0] = sqrt(slen[0]);
	slen[1] = sqrt(slen[1]);
	slen[2] = sqrt(slen[2]);
	slen[3] = sqrt(slen[3]);
	
	scorn[0] = s6ang(t,t+kdim,kdim);
	scorn[1] = s6ang(t+kdim,t+2*kdim,kdim);
	scorn[2] = s6ang(t+2*kdim,t+3*kdim,kdim);
	scorn[3] = s6ang(t+3*kdim,t,kdim);
	
	/* If problems on edges is found we jump to the surface. */
	
	if (ps->pdir->igtpi > 0) goto next;
	
	/* Computing cones of edges in ends of parameter two. */
	
	if (kver == 0)
	  {
	    if (lcone[0])
	      {
		/* First time to generate cone. */
		 
		 memcopy(tsen,t+kdim,kdim,DOUBLE);
		 tlen = slen[1];
		
		if (tlen > aepsge)
		  {
		    for (k1=0; k1 < kdim; k1++) tsen[k1] /= tlen;
		    lcone[0] = 0;
		    sang[0] = (double)0;
		  }
	      }
	    else
	      {
		/* Modify existing cone. */
		 s1990_s9edg(t+(kdim),ttan,tsen,aepsge,sang,kdim,&kstat);
		
		if (kstat)   ps->pdir->igtpi = 10;
	      }
	  } 
	if (kver == kn2-2)
	  {
	    if (lcone[1])
	      {
		/* First time to generate cone. */
		 
		 memcopy(tsen+kdim,t+3*kdim,kdim,DOUBLE);
		 tlen = slen[3];
		
		if (tlen > aepsge)
		  {
		    for (k1=0; k1 < kdim; k1++) tsen[kdim+k1] /= tlen;
		    lcone[1] = 0;
		    sang[1] = (double)0;
		  }
	      }
	    else
	      {
		 s1990_s9edg(t+(3*kdim),ttan,tsen+kdim,aepsge,sang+1,kdim,&kstat);
		if (kstat) ps->pdir->igtpi = 10;
	      }
	  }
	
	/* Computing cones of edges in ends of parameter one. */
	
	if (khor == 0)
	  {
	    if (lcone[2])
	      /* First time to generate cone. */
	      {
		 memcopy(tsen+2*kdim,t,kdim,DOUBLE);
		 tlen = slen[0];
		
		if (tlen > aepsge)
		  {
		    for (k1=0; k1 < kdim; k1++) tsen[2*kdim+k1] /= tlen;
		    lcone[2] = 0;
		    sang[2] = (double)0;
		  }
	      }
	    else
	      {
		 s1990_s9edg(t,ttan,tsen+(2*kdim),aepsge,sang+2,kdim,&kstat);
		if (kstat) ps->pdir->igtpi = 10;
	      }
	  } 
	if (khor == kn1-2)
	  {
	    if (lcone[3])
	      {
		 memcopy(tsen+3*kdim,t+2*kdim,kdim,DOUBLE);
		 tlen = slen[2];
		
		if (tlen > aepsge)
		  {
		    for (k1=0; k1 < kdim; k1++) tsen[3*kdim+k1] /= tlen;
		    lcone[3] = 0;
		    sang[3] = (double)0;
		  }
	      }
	    else
	      {
		 s1990_s9edg(t+(2*kdim),ttan,tsen+(3*kdim),aepsge,sang+3,kdim,&kstat);
		if (kstat)  ps->pdir->igtpi = 10;
	      }
	  }
	
      next:
	
	/* Here we makes the normales in each corner of the patch.
	   We are using a cross product between two tangents.
	   The normals is also normalized by deviding with its
	   own length. */
	
	for (kcount=0, ki=0, k1=0; k1 < kdim4; k1+=kdim, ki++)
	  {
	    for (tlen=DZERO,k2=0,k3=1,k4=2; k2 < kdim; k2++,k3++,k4++)
	      {
		if(k3 == kdim) k3 = 0;
		if(k4 == kdim) k4 = 0;
		tn[k1+k2] = t[k1+k3]*t[k1+kdim+k4]-t[k1+k4]*t[k1+kdim+k3];
		
		tlen += tn[k1+k2]*tn[k1+k2];
	      }
	    tlen = sqrt(tlen);
	    /* KYS 070494 : multiplied ANGULAR_TOLERANCE by 1.0e-2 */
	    if (slen[ki]>aepsge && slen[ki+1]>aepsge &&
		scorn[ki] > 1.0e-2*ANGULAR_TOLERANCE)
	      for (k2=0; k2 < kdim; k2++) tn[k1+k2] /= tlen;
	    else 
	      {
	      for (k2=0; k2 < kdim; k2++) tn[k1+k2] = ps->pdir->ecoef[k2];
	      kcount++;
	      }
	  }
	
	if (kcount == 4) continue;   /* Degenerate control polygon patch */
	
	/* We are treating the first patch. */
	
	if (kfirst)
	  {
	    /* Computing the center coordinates of the cone.*/
	    
	    for (tlen=DZERO,k1=0; k1 < kdim; k1++)
	      {
		tmin = (double)1.0;
		tmax = - tmin;
		for (k2=0; k2 < kdim4; k2+=kdim)
		  {
		    tmax = max(tn[k2+k1],tmax);
		    tmin = min(tn[k2+k1],tmin);
		  }
		ps->pdir->ecoef[k1]=(tmax+tmin)/(double)2;
		
		tlen += ps->pdir->ecoef[k1]*ps->pdir->ecoef[k1];
	      }
	    tlen = sqrt(tlen);
	    if (tlen > DZERO)
	      for (k1=0; k1 < kdim; k1++) ps->pdir->ecoef[k1] /= tlen;
	    else
	      /* KYS 070494 : 'continue' replaced by the following block {} */
	      /* There are nonzero normals pointing in
		 opposite directions, i.e. not simple case */
	      {
		if (khor <= kver)
		  ps->pdir->igtpi = 1;
		else
		  ps->pdir->igtpi = 2;
		ps->pdir->aang = PI;
		goto out;
	      }
	    
	    /* Computing the angle of the cone. */
	    
	    for (ps->pdir->aang=DZERO,k1=0; k1<kdim4; k1+=kdim)
	      {
		 for (tnlen=DZERO,tlen=DZERO,k2=0;k2<kdim;k2++)
		   {
		      tlen += ps->pdir->ecoef[k2]*tn[k1+k2];
		      tnlen += tn[k1+k2]*tn[k1+k2];
		   }
		
		if (tlen >= DZERO) tlen = min((double)1.0,tlen);
		else               tlen = max((double)-1.0,tlen);
		
		tlen = acos(tlen);
		if (sqrt(tnlen) < aepsge) tlen = DZERO;
		
		ps->pdir->aang = max(ps->pdir->aang,tlen);
	      }
	    
	    kfirst = 0;   /* The first patch have been treated.*/
	  } 
	else
	  for (k1=0; k1<kdim4; k1+=kdim)
	    {
	      /* Computing the angle beetween the senter of the cone
		 and the normal. */
	      
	      for (tnlen=DZERO,tang=DZERO,k2=0;k2<kdim;k2++)
		{
		   tang += ps->pdir->ecoef[k2]*tn[k1+k2];
		   tnlen += tn[k1+k2]*tn[k1+k2];
		}
	      
	      if (tang >= DZERO) tang = MIN((double)1.0,tang);
	      else               tang = MAX((double)-1.0,tang);
	      
	      tang = acos(tang);
	      if (sqrt(tnlen) < aepsge) tang = DZERO;
	      
	      if (tang + ps->pdir->aang >= PI)
		{
		  /* The angle is to great, give a meesage
		     how to subdivied and exit this function. */
		  
		  if (khor <= kver)
		    ps->pdir->igtpi = 1;
		  else	
		    ps->pdir->igtpi = 2;
		  goto out;
		}
	      else if (tang > ps->pdir->aang)
		{
		  /* The normal is not inside the cone, than we
		     have to compute a new cone. */
		  
		  /* Computing the center coordinates.*/
		  
	          double sin_tang = sin(tang);                     /*@  hke  */
	          double delta    = (tang - ps->pdir->aang)/2.0;   /*@  hke  */

	          t1 = sin(delta)/sin_tang;                        /*@  hke  */
	          t2 = sin(tang - delta)/sin_tang;                 /*@  hke  */

		  /*
		  t1 = (tang - ps->pdir->aang)/((double)2*tang);
		  t2 = (double)1 - t1;
		  */
		  
		  for (tlen=DZERO,k2=0; k2<kdim; k2++)
		    {
		      ps->pdir->ecoef[k2] = 
			ps->pdir->ecoef[k2]*t2 + tn[k1+k2]*t1;
		      tlen += ps->pdir->ecoef[k2]*ps->pdir->ecoef[k2];
		    }
		  tlen = sqrt(tlen);
		  
		  for (k2=0; k2 < kdim; k2++)  ps->pdir->ecoef[k2] /= tlen;
		  
		  /* Computing the angle of the cone. */
		  
		  ps->pdir->aang = (tang + ps->pdir->aang)/(double)2;
		}
	    }
	
	if (ps->pdir->aang >= SIMPLECASE)
	  {
	    /* The angle is to great, give a meesage
	       how to subdivied and exit this function. */
	    
	    if (khor <= kver)
	      ps->pdir->igtpi = 10;
	    else	
	      ps->pdir->igtpi = 20;
	  }
      }			
  
  /* A final check if we have made a cone. */
  /* UJK, SI, 91-10, when 2D, return values from edge case */
  if (kfirst && kdim != 2)
    {
      /* No cone has been generated. We must examin if the surface is 
	 degenerated to a point or line. */
      for (k1 = 1; k1 < kn1*kn2; k1++)
	if (s6dist(scoef,scoef + (k1*kdim),kdim) >aepsge) break;
      
      if (k1 == kn1*kn2)
	{
	  /* Degenerated to a point. */
	  ps->pdir->igtpi = 0;
	  ps->pdir->aang  = DZERO;
	  ps->pdir->ecoef[0] = (double) 1.0;
	  for (k1 = 1; k1 < kdim; k1++) ps->pdir->ecoef[k1] = DZERO;
	}
      else
	{
	  s6diff(scoef,scoef + (k1*kdim),kdim,svec1);
	  
	  for (k2 = k1 + 1; k2 < kn1*kn2; k2++)
	    if (s6dist(scoef,scoef + (k2*kdim),kdim) >aepsge)
	      {
		s6diff(scoef,scoef + (k2*kdim),kdim,svec2);
		if (s6ang(svec1,svec2,kdim) > 1.0e-2*ANGULAR_TOLERANCE) break;
	      }
	  
	  if (k2 == kn1*kn2)
	    {
	      /* Degenerated to a line. */
	      ps->pdir->igtpi = 0;
	      ps->pdir->aang  = DZERO;
	      ps->pdir->ecoef[0] = (double) 1.0;
	      for (k1 = 1; k1 < kdim; k1++) ps->pdir->ecoef[k1] = DZERO;
	    }
	  else
	    {
	       /* Three points describing a plane found, continue subdividing. */
	       if (ps->et1[kn1] - ps->et1[ps->ik1-1] >=
		   ps->et2[kn2] - ps->et2[ps->ik2-1])
		  ps->pdir->igtpi = 1;
	       else
	       ps->pdir->igtpi = 2; 
	    }
	}
    }
  
  /* success */
  
  goto out;
  
  /* Error in space allacation.  */
  
  err101: 
    *jstat = -101;
    s6err("s1990",*jstat,kpos);
    goto out;
  
  /* Error in lower level routine.  */
  
  /* error : 
    *jstat = kstat;
    goto out; 
  */
  
  /* Free local used memory. */
  
  out:    
    if (t != SISL_NULL) freearray(t);
}

#if defined(SISLNEEDPROTOTYPES)
static  void
  s1990_s9edg(double et[],double etan[],double esen[],double aepsge,
	      double *cang,int idim,int *jstat)
#else
static void s1990_s9edg(et,etan,esen,aepsge,cang,idim,jstat)
     double et[];
     double etan[];
     double esen[];
     double aepsge;
     double *cang;
     int    idim;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To make the orientation surface on the unit sphere to
*	       the edges of a b-spline surface, the surface is
*	       representated with a surrounding cone piced from
*              the unit sphere.
*
*
*
* INPUT      : et[]    - The tangent vector to the edge.
*              etan[]  - The normalized tangent vector to the edge.
*              idim    - The dimention of the surface.
*
*
* INPUT/OUTPUT:esen[]  - The senter vector of the cone.
*              cang    - The angel of the cone..
*
*
*
* OUTPUT     : jstat  - status messages  
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
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-05.
*
*********************************************************************
*/
{
  int ki;
  double tlen;
  double tang;
  double t1,t2;
  
  
  /* Normalizing the tangent. */
  
  for (tlen = DZERO,ki=0; ki < idim; ki++)
    {
      etan[ki] = et[ki];
      tlen += etan[ki]*etan[ki];
    }
  tlen = sqrt(tlen);
  
  if (tlen > aepsge)
    for (ki=0; ki < idim; ki++) etan[ki] /= tlen;
  else
    {
      *jstat = 0;
      goto out;
    }
  
  
  /* Computing the angle beetween the senter of the cone
     and the tangent. */
  
  for (tang=DZERO,ki=0;ki<idim;ki++)
    tang += esen[ki]*etan[ki];
  
  if (tang >= DZERO) tang = min((double)1.0,tang);
  else               tang = max((double)-1.0,tang);
  
  tang = acos(tang);
  
  
  if (tang + *cang >= PI)
    {
      /* The angle is to great, give a meesage
	 to subdivied and exit this function. */
      
      *jstat = 1;
      goto out;
    }
  else if (tang > *cang)
    {
      /* The tangent is not inside the cone, and we
	 have to compute a new cone. */
      
      /* Computing the center coordinates.*/
      
      t1 = (tang - *cang)/((double)2*tang);
      t2 = (double)1 - t1;
      
      for (tlen=DZERO,ki=0; ki<idim; ki++)
        {
	  esen[ki] = esen[ki]*t2 + etan[ki]*t1;
	  tlen += esen[ki]*esen[ki];
        }
      tlen = sqrt(tlen);
      
      if (tlen > DZERO)
	for (ki=0; ki < idim; ki++) esen[ki] /= tlen;
      else
	{
	  /* Vi have to be aware of colapsed polygon. */
	  
	  *jstat = 1;
	  goto out;
	}
      
      /* Computing the angle of the cone. */
      
      *cang = (tang + *cang)/(double)2;
    }
  
  
  if (*cang >= SIMPLECASE)
    {
      /* The angle is to large, give a meesage
	 to subdivied and exit this function. */
      
      *jstat = 1;
      goto out;
    }
  
  
  *jstat = 0;
  
 out: ;
}

#if 0   
#if defined(SISLNEEDPROTOTYPES)
static  void
  s1990_s9smooth(double ecoef1[],int in1,int in2,int idim,
		 double aepsge,double ecoef2[],int *jstat)
#else
static void s1990_s9smooth(ecoef1,in1,in2,idim,aepsge,ecoef2,jstat)
   double ecoef1[];
   int    in1;
   int    in2;
   int    idim;
   double aepsge;
   double ecoef2[];
   int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Perform noise filthering at the corners of the control
*              polygon of a B-spline surface.
*
*
*
* INPUT      : ecoef1 - Original coefficients of surface
*              in1    - Number of coefficients in 1. par dir.
*              in2    - Number of coefficients in 2. par dir.
*              idim   - Dimension of geometry space.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : ecoef2 - New coefficients after smoothing.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*              
*
*
* REFERENCES :
*
*-
* CALLS      : s6dplane  -  Distance to given plane.
*              s6dline   -  Distance to given line.
*              s6dist    -  Distance between two points.
*
* WRITTEN BY : Vibeke Skytt, SI, 91-02.
*
*********************************************************************
*/
{
   int kstat = 0;     /* Local status variable.        */
   int kn = MIN(in1/2,in2/2)+1;  /* Maximum numbers of 
				  coefficients to smooth. */
   int ki,kj,kh,kl;   /* Counters.                     */
   int kc;            /* Index of current corner.      */
   int k1;            /* Sign of change in 1. par dir  */
   int k2;            /* Sign of change in 2. par dir  */
   int lcorn[4];      /* Indexes of corners.           */
   int lsgn1[4];      /* Sign of changes in 1. par dir */
   int lsgn2[4];      /* Sign of changes in 2. par dir */
   double tdist;      /* Distance to closest point in plane. */
   
   /* Set contents of arrays.  */
   
   lcorn[0] = 0;
   lcorn[1] = (in1-1)*idim;
   lcorn[2] = (in1*in2-1)*idim;
   lcorn[3] = in1*(in2-1)*idim;
   
   lsgn1[0] = 1;
   lsgn1[1] = -1;
   lsgn1[2] = -1;
   lsgn1[3] = 1;
   
   lsgn2[0] = 1;
   lsgn2[1] = 1;
   lsgn2[2] = -1;
   lsgn2[3] = -1;
   
   /* Copy coefficients to output array.  */
   
   memcopy(ecoef2,ecoef1,in1*in2*idim,DOUBLE);

   /* For each corner, try to smooth the coefficients in the
      neighbourhood of the corner.  */
   
   for (ki=0; ki<4; ki++)
   {
      kc = lcorn[ki];   /* Index of current corner.   */
      k1 = lsgn1[ki];   /* Sign change in 1. par dir. */
      k2 = lsgn2[ki];   /* Sign change in 2. par dir. */
      
      /* Try to smooth coefficients on center line.  */
	 
      for (kj=2; kj<kn; kj++)
      {
	 if (s6dist(ecoef2+kc,ecoef2+kc+(k2*kj*in1+k1*kj)*idim,
		    idim) < aepsge) continue;
	 
	 for (kh=1; kh<kj; kh++)
	 {
	    tdist = s6dline(ecoef2+kc,ecoef2+kc+(k2*kj*in1+k1*kj)*idim,
			    ecoef2+kc+(k2*kh*in1+k1*kh)*idim,idim,&kstat);
	    if (kstat < 0) goto error;
	    if (kstat || tdist >= aepsge) break;
	 }
	 if (kh < kj) break;
      }
      
      /* Perform smoothing.  */
      
      kj--;
      for (kh=1; kh<kj; kh++)
	 memcopy(ecoef2+kc+(k2*kh*in1+k1*kh)*idim,ecoef2+kc,
		 idim,DOUBLE);
      
      /* Try to smooth coefficients on lower triangle.  */
      
      for (kj=2; kj<kn; kj++)
      {
	 for (kh=1; kh<kj; kh++)
	 {
	    for (kl=0; kl<kh; kl++)
	    {
	       tdist = s6dplane(ecoef2+kc,ecoef2+kc+k1*kj*idim,
				ecoef2+kc+(k2*kj*in1+k1*kj)*idim,
			        ecoef2+kc+(k2*kl*in1+k1*kh)*idim,
				idim,&kstat);
	       if (tdist >= aepsge) break;
	    }
	    if (tdist >= aepsge) break;
	 }
	 if (kh < kj) break;
      }
      
      /* Perform smoothing.  */
      
      kj--;
      for (kh=1; kh<kj; kh++)
	 for (kl=0; kl<kh; kl++)
	    memcopy(ecoef2+kc+(k2*kl*in1+k1*kh)*idim,ecoef2+kc,
		    idim,DOUBLE);
      
      /* Try to smooth coefficients on upper triangle.  */
      
      for (kj=2; kj<kn; kj++)
      {
	 for (kh=0; kh<kj; kh++)
	 {
	    for (kl=kh+1; kl<kj; kl++)
	    {
	       tdist = s6dplane(ecoef2+kc,ecoef2+kc+k2*kj*in1*idim,
				ecoef2+kc+(k2*kj*in1+k1*kj)*idim,
			        ecoef2+kc+(k2*kl*in1+k1*kh)*idim,
				idim,&kstat);
	       if (tdist >= aepsge) break;
	    }
	    if (tdist >= aepsge) break;
	 }
	 if (kh < kj) break;
      }
      
      /* Perform smoothing.  */
      
      kj--;
      for (kh=0; kh<kj; kh++)
	 for (kl=kh+1; kl<kj; kl++)
	    memcopy(ecoef2+kc+(k2*kl*in1+k1*kh)*idim,ecoef2+kc,
		    idim,DOUBLE);
   }
   
   /* Smoothing performed. */
   *jstat = 0;
   goto out;
   
   /* Error in lower level routine.  */
   
   error : *jstat = kstat;
   goto out;
   
   out :
      
   return;
}
 
#endif /* if 0 */
