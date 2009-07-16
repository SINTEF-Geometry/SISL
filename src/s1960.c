/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s1960.c,v 1.2 2001-03-19 15:58:58 afr Exp $
 *
 */
#define S1960

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1960(SISLPoint *ppoint, SISLSurf *psurf, double gpos[], int *jstat)
#else
void s1960(ppoint,psurf,gpos,jstat)
     SISLPoint        *ppoint;
     SISLSurf         *psurf;
     double       gpos[];
     int          *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Estimate parameter-pair of guess-point (used by closest point
*              calculation).
*
*
* INPUT      : ppoint   - Pointer to the point.
*              psurf    - Pointer to the surface.
*
* OUTPUT     : gpos    - Parameter values of the found guess-point.
*              jstat   - status messages  
*                                = 0   : Guess-point found.
*                                = 1   : Closest point as guess-point.
*                                < 0   : error.
*
*
* METHOD     : Quadrant analysis and Schoenberger equation.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Per Evensen, SI, August 1991
* REVISED BY : Michael Floater, SI, December 1991
*
*********************************************************************
*/                       
{  
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  int i,j;                  /* Running indexes                             */
  int iind,jind;            /* Index variables                             */
  int k1;                   /* The polynomial order in 1. parameter 
                               direction of the surface (psurf)            */
  int k2;                   /* The polynomial order in 2. parameter 
                               direction of the surface (psurf)            */
  int n1;                   /* The number of vertices in 1. parameter 
                               direction of the surface (psurf)            */
  int n2;                   /* The number of vertices in 2. parameter 
                               direction of the surface (psurf)            */
  int dim;                  /* Dimension of space the surface lies in      */
  double *et1;              /* Knots in 1. parameter direction of the 
                               surface (psurf)                             */
  double *et2;              /* Knots in 2. parameter direction of the 
                               surface (psurf)                             */
  double *lpoint;           /* Coefficients of the point (ppoint)          */
  double *scoef;            /* Coefficients of the surface (psurf)         */
  double tdist;             /* Distance variable                           */
  double tmin;              /* Minimum distance variable                   */
  double vec1[3],vec2[3],vec3[3],vec4[3];
                            /* Vectors defining the quadrants surrounding 
                               a vertice                                   */
  double vecp[3];           /* Relative point vector                       */
  double lqua1=DZERO;
  double lqua2=DZERO;
  double lqua3=DZERO;
  double lqua4=DZERO;
                            /* Length of vec1,....,vec4                    */
  double lprj1=DZERO;
  double lprj2=DZERO;
  double lprj3=DZERO;
  double lprj4=DZERO;
                            /* Length of projection of vecp on 
                               vec1,....,vec4                              */
  double svals,svale,tvals,tvale;
                            /* Local start and end parameter values in s 
                               and t direction                             */
  /* --------------------------------------------------------------------- */
  
  /* Test input.  */
  if (ppoint->idim != psurf->idim || psurf->idim <= 1) goto err106;

  /* Initiate local variables.  */
  k1  = psurf->ik1;
  k2  = psurf->ik2;
  n1  = psurf->in1;
  n2  = psurf->in2;
  et1 = psurf->et1;
  et2 = psurf->et2;
  scoef = psurf->ecoef;
  dim = psurf->idim;
  lpoint = ppoint->ecoef;
   
  /* Find vertice closest to point.  */
  tdist=s6dist(scoef,lpoint,dim);
  tmin=tdist;
  iind = 0;
  jind = 0;
  for (i=0; i<n1; i++)
  {
     for (j=0; j<n2; j++)
     {
        tdist=s6dist(scoef,lpoint,dim);
        if (tdist<tmin)
        {
           tmin=tdist;
           iind = i;
           jind = j;
        }
        scoef+=3;
     }
  }
  
  /* 
  * Generate the "vectors of the quadrant".
  
                     I vec2
                     l
           vec3      l       vec 1
             <-------X------->
closest            / l   x
vertice (iind,jind)  l    \ point
                     I vec4
  */
  scoef = psurf->ecoef;
  
  if (iind < (n1-1))
     s6diff(&(scoef[(iind+1)*dim+jind*dim*n1]),
            &(scoef[iind*dim+jind*dim*n1]),dim,vec1);
  if (jind < (n2-1))
     s6diff(&(scoef[iind*dim+(jind+1)*dim*n1]),
            &(scoef[iind*dim+jind*dim*n1]),dim,vec2);
  if (iind > 0)
     s6diff(&(scoef[(iind-1)*dim+jind*dim*n1]),
            &(scoef[iind*dim+jind*dim*n1]),dim,vec3);
  if (jind > 0)
     s6diff(&(scoef[iind*dim+(jind-1)*dim*n1]),
            &(scoef[iind*dim+jind*dim*n1]),dim,vec4);
  
  /* Generate the point - closest vertice vector. */
  s6diff(lpoint,&(scoef[iind*dim+jind*dim*n1]),dim,vecp);
     
  /* Calculate the length of the quadrant vectors. */
  if (iind < (n1-1)) lqua1 = s6length(vec1,dim,&kstat);
  if (jind < (n2-1)) lqua2 = s6length(vec2,dim,&kstat);
  if (iind > 0) lqua3 = s6length(vec3,dim,&kstat);
  if (jind > 0) lqua4 = s6length(vec4,dim,&kstat);

  /* Calculate the length of the projection of 'vecp' on the quadrant 
     vectors. */
  if (iind < (n1-1)) lprj1 = s6lprj(vecp,vec1,dim);
  if (jind < (n2-1)) lprj2 = s6lprj(vecp,vec2,dim);
  if (iind > 0) lprj3 = s6lprj(vecp,vec3,dim);
  if (jind > 0) lprj4 = s6lprj(vecp,vec4,dim);

  /* Decide in which quadrant the point lies. */
  if (iind == 0 && jind == 0)
  {
     /* Point lies in 1. quadrant */
     
     /* Calculate knot values of vertices */
     svals = s6schoen(et1,k1,iind);
     svale = s6schoen(et1,k1,iind+1);
     tvals = s6schoen(et2,k2,jind);
     tvale = s6schoen(et2,k2,jind+1);
     
     /* Calculate estimated parameter values of point */
     if (lqua1 != DZERO) gpos[0] = svals + (lprj1/lqua1)*(svale-svals);
     else gpos[0] = svals;
     if (lqua2 != DZERO) gpos[1] = tvals + (lprj2/lqua2)*(tvale-tvals);
     else gpos[1] = tvals;
  }
  else if (iind == (n1-1) && jind == 0)
  {
     /* Point lies in 2. quadrant */
     
     /* Calculate knot values of vertices */
     svals = s6schoen(et1,k1,iind-1);
     svale = s6schoen(et1,k1,iind);
     tvals = s6schoen(et2,k2,jind);
     tvale = s6schoen(et2,k2,jind+1);
     
     /* Calculate estimated parameter values of point */
     if (lqua3 != DZERO) gpos[0] = svals + ((lqua3-lprj3)/lqua3)*(svale-svals);
     else gpos[0] = svals;
     if (lqua2 != DZERO) gpos[1] = tvals + (lprj2/lqua2)*(tvale-tvals);
     else gpos[1] = tvals;
  }
  else if (iind == (n1-1) && jind == (n2-1))
  {
     /* Point lies in 3. quadrant */
     
     /* Calculate knot values of vertices */
     svals = s6schoen(et1,k1,iind-1);
     svale = s6schoen(et1,k1,iind);
     tvals = s6schoen(et2,k2,jind-1);
     tvale = s6schoen(et2,k2,jind);
     
     /* Calculate estimated parameter values of point */
     if (lqua3 != DZERO) gpos[0] = svals + ((lqua3-lprj3)/lqua3)*(svale-svals);
     else gpos[0] = svals;
     if (lqua4 != DZERO) gpos[1] = tvals + ((lqua4-lprj4)/lqua4)*(tvale-tvals);
     else gpos[1] = tvals;
  }
  else if (iind == 0 && jind == (n2-1))
  {
     /* Point lies in 4. quadrant */
     
     /* Calculate knot values of vertices */
     svals = s6schoen(et1,k1,iind);
     svale = s6schoen(et1,k1,iind+1);
     tvals = s6schoen(et2,k2,jind-1);
     tvale = s6schoen(et2,k2,jind);
     
     /* Calculate estimated parameter values of point */
     if (lqua1 != DZERO) gpos[0] = svals + (lprj1/lqua1)*(svale-svals);
     else gpos[0] = svals;
     if (lqua4 != DZERO) gpos[1] = tvals + ((lqua4-lprj4)/lqua4)*(tvale-tvals);
     else gpos[1] = tvals;
  }
  else if (iind > 0 && iind < (n1-1) && jind == 0)
  {
     /* Evaluate 1. and 2. quadrant */
     
     if (lprj1 > lprj3)
        {
           /* Point lies in 1. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind);
           svale = s6schoen(et1,k1,iind+1);
           tvals = s6schoen(et2,k2,jind);
           tvale = s6schoen(et2,k2,jind+1);
           
           /* Calculate estimated parameter values of point */
           if (lqua1 != DZERO) gpos[0] = svals + (lprj1/lqua1)*(svale-svals);
           else gpos[0] = svals;
           if (lqua2 != DZERO) gpos[1] = tvals + (lprj2/lqua2)*(tvale-tvals);
           else gpos[1] = tvals;
        }
        else if (lprj3 > lprj1)
        {
           /* Point lies in 2. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind-1);
           svale = s6schoen(et1,k1,iind);
           tvals = s6schoen(et2,k2,jind);
           tvale = s6schoen(et2,k2,jind+1);
           
           /* Calculate estimated parameter values of point */
           if (lqua3 != DZERO) gpos[0] = svals + ((lqua3-lprj3)/lqua3)*(svale-svals);
           else gpos[0] = svals;
           if (lqua2 != DZERO) gpos[1] = tvals + (lprj2/lqua2)*(tvale-tvals);
           else gpos[1] = tvals;
        }
        else
        {
           /* lprj1 and lprj3 are equal. */
           /* Choose original control point. */
           goto usvert;
        }
  }
  else if (iind == (n1-1) && jind > 0 && jind < (n2-1))
  {
     /*  Evaluate 2. and 3. quadrant */
     
     if (lprj2 > lprj4)
        {
           /* Point lies in 2. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind-1);
           svale = s6schoen(et1,k1,iind);
           tvals = s6schoen(et2,k2,jind);
           tvale = s6schoen(et2,k2,jind+1);
           
           /* Calculate estimated parameter values of point */
           if (lqua3 != DZERO) gpos[0] = svals + ((lqua3-lprj3)/lqua3)*(svale-svals);
           else gpos[0] = svals;
           if (lqua2 != DZERO) gpos[1] = tvals + (lprj2/lqua2)*(tvale-tvals);
           else gpos[1] = tvals;
        }
        else if (lprj4 > lprj2)
        {
           /* Point lies in 3. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind-1);
           svale = s6schoen(et1,k1,iind);
           tvals = s6schoen(et2,k2,jind-1);
           tvale = s6schoen(et2,k2,jind);
           
           /* Calculate estimated parameter values of point */
           if (lqua3 != DZERO) gpos[0] = svals + ((lqua3-lprj3)/lqua3)*(svale-svals);
           else gpos[0] = svals;
           if (lqua4 != DZERO) gpos[1] = tvals + ((lqua4-lprj4)/lqua4)*(tvale-tvals);
           else gpos[1] = tvals;
        }
        else
        {
           /* lprj2 and lprj4 are equal. */
           /* Choose original control point. */
           goto usvert;
        }
  }
  else if (iind > 0 && iind < (n1-1) && jind == (n2-1))
  {
     /*  Evaluate 3. and 4. quadrant */
     
     if (lprj3 > lprj1)
        {
           /* Point lies in 3. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind-1);
           svale = s6schoen(et1,k1,iind);
           tvals = s6schoen(et2,k2,jind-1);
           tvale = s6schoen(et2,k2,jind);
           
           /* Calculate estimated parameter values of point */
           if (lqua3 != DZERO) gpos[0] = svals + ((lqua3-lprj3)/lqua3)*(svale-svals);
           else gpos[0] = svals;
           if (lqua4 != DZERO) gpos[1] = tvals + ((lqua4-lprj4)/lqua4)*(tvale-tvals);
           else gpos[1] = tvals;
        }
        else if (lprj1 > lprj3)
        {
           /* Point lies in 4. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind);
           svale = s6schoen(et1,k1,iind+1);
           tvals = s6schoen(et2,k2,jind-1);
           tvale = s6schoen(et2,k2,jind);
           
           /* Calculate estimated parameter values of point */
           if (lqua1 != DZERO) gpos[0] = svals + (lprj1/lqua1)*(svale-svals);
           else gpos[0] = svals;
           if (lqua4 != DZERO) gpos[1] = tvals + ((lqua4-lprj4)/lqua4)*(tvale-tvals);
           else gpos[1] = tvals;
        }
        else
        {
           /* lprj1 and lprj3 are equal. */
           /* Choose original control point. */
           goto usvert;
        }
  }
  else if (iind == 0 && jind > 0 && jind < (n2-1))
  {
     /*  Evaluate 4. and 1. quadrant */
     
     if (lprj2 > lprj4)
        {
           /* Point lies in 1. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind);
           svale = s6schoen(et1,k1,iind+1);
           tvals = s6schoen(et2,k2,jind);
           tvale = s6schoen(et2,k2,jind+1);
           
           /* Calculate estimated parameter values of point */
           if (lqua1 != DZERO) gpos[0] = svals + (lprj1/lqua1)*(svale-svals);
           else gpos[0] = svals;
           if (lqua2 != DZERO) gpos[1] = tvals + (lprj2/lqua2)*(tvale-tvals);
           else gpos[1] = tvals;
        }
        else if (lprj4 > lprj2)
        {
           /* Point lies in 4. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind);
           svale = s6schoen(et1,k1,iind+1);
           tvals = s6schoen(et2,k2,jind-1);
           tvale = s6schoen(et2,k2,jind);
           
           /* Calculate estimated parameter values of point */
           if (lqua1 != DZERO) gpos[0] = svals + (lprj1/lqua1)*(svale-svals);
           else gpos[0] = svals;
           if (lqua4 != DZERO) gpos[1] = tvals + ((lqua4-lprj4)/lqua4)*(tvale-tvals);
           else gpos[1] = tvals;
        }
        else
        {
           /* lprj2 and lprj4 are equal. */
           /* Choose original control point. */
           goto usvert;
        }
  }
  else if (iind > 0 && iind < (n1-1) && jind > 0 && jind < (n2-1))
  {
     /* Evaluate all four quadrants */
     
     if (lprj1 > lprj3)
        {
           /* Point lies in 1. or 4. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind);
           svale = s6schoen(et1,k1,iind+1);
           
           /* Calculate estimated parameter values of point */
           if (lqua1 != DZERO) gpos[0] = svals + (lprj1/lqua1)*(svale-svals);
           else gpos[0] = svals;
        }
        else if (lprj3 > lprj1)
        {
           /* Point lies in 2. or 3. quadrant */
           
           /* Calculate knot values of vertices */
           svals = s6schoen(et1,k1,iind-1);
           svale = s6schoen(et1,k1,iind);
           
           /* Calculate estimated parameter values of point */
           if (lqua3 != DZERO) gpos[0] = svals + ((lqua3-lprj3)/lqua3)*(svale-svals);
           else gpos[0] = svals;
        }
        else
        {
           /* lprj1 and lprj3 are equal. */
           /* Choose original control point. */
           gpos[0] = s6schoen(et1,k1,iind);
        }


        if (lprj2 > lprj4)
        {
           /* Point lies in 1. or 2. quadrant */
           
           /* Calculate knot values of vertices */
           tvals = s6schoen(et2,k2,jind);
           tvale = s6schoen(et2,k2,jind+1);
           
           if (lqua2 != DZERO) gpos[1] = tvals + (lprj2/lqua2)*(tvale-tvals);
           else gpos[1] = tvals;
        }


        else if (lprj4 > lprj2)
        {
           /* Point lies in 3. or 4. quadrant */
           
           /* Calculate knot values of vertices */
           tvals = s6schoen(et2,k2,jind-1);
           tvale = s6schoen(et2,k2,jind);
           
           if (lqua4 != DZERO) gpos[1] = tvals + ((lqua4-lprj4)/lqua4)*(tvale-tvals);
           else gpos[1] = tvals;
        }
        else
        {
           /* lprj2 and lprj4 are equal. */
           /* Choose original control point. */
           gpos[1] = s6schoen(et2,k2,jind);
        }
  }
  else
  {
     /* Error */
     goto usvert;
  }
  
  /* Check that values are within parameter plane of surface. */
  if (gpos[0]<et1[k1-1]) gpos[0]=et1[k1-1];
  else if (gpos[0]>et1[n1]) gpos[0]=et1[n1];
  if (gpos[1]<et2[k2-1]) gpos[1]=et2[k2-1];
  else if (gpos[1]>et2[n2]) gpos[1]=et2[n2];
  *jstat = 0;
  
  /* Calculation completed.  */
  goto out;
  
  /* No intermediate parameter values found,
     use parameter values of closest vertice */
  usvert: *jstat = 1;
           
  /* Calculate knot values of closest vertice */
  gpos[0] = s6schoen(et1,k1,iind);
  gpos[1] = s6schoen(et2,k2,jind);

  /* Check that values are within parameter plane of surface. */
  if (gpos[0]<et1[k1-1]) gpos[0]=et1[k1-1];
  else if (gpos[0]>et1[n1]) gpos[0]=et1[n1];
  if (gpos[1]<et2[k2-1]) gpos[1]=et2[k2-1];
  else if (gpos[1]>et2[n2]) gpos[1]=et2[n2];
  goto out;                  
  
 /* --------------------------------------------------------------------- */ 
  /* Error in input. Dimension not equal to 1 */
 err106: *jstat = -106;
  s6err("s1960",*jstat,kpos);
  goto out;   
  
 out:
    return;
}

