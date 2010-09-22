/* SISL - SINTEF Spline Library version 4.4.                              */
/* Definition and interrogation of NURBS curves and surface.              */
/* Copyright (C) 1978-2005, SINTEF ICT, Applied Mathematics, Norway.      */

/* This program is free software; you can redistribute it and/or          */
/* modify it under the terms of the GNU General Public License            */
/* as published by the Free Software Foundation version 2 of the License. */

/* This program is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of         */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          */
/* GNU General Public License for more details.                           */

/* You should have received a copy of the GNU General Public License      */
/* along with this program; if not, write to the Free Software            */
/* Foundation, Inc.,                                                      */
/* 59 Temple Place - Suite 330,                                           */
/* Boston, MA  02111-1307, USA.                                           */

/* Contact information: e-mail: tor.dokken@sintef.no                      */
/* SINTEF ICT, Department of Applied Mathematics,                         */
/* P.O. Box 124 Blindern,                                                 */
/* 0314 Oslo, Norway.                                                     */

/* SISL commercial licenses are available for:                            */
/* - Building commercial software.                                        */
/* - Building software whose source code you wish to keep private.        */


#define DBGqqq
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#include "sisl_aux.h"
#include "aux2.h"






//----------------------------------------------------------------------
//
// 000817: Discretizing a surface for plotting by lowering the degree
//         to 1 in both direction, so that we can plot the control
//         polygon.
// 000829: Oopps! s1968 doesn't do general degree lowering!
//
//----------------------------------------------------------------------

void lower_degree_to_linear(SISLSurf **srf, double e)
{
  CRIT_ERR(puts("Not implemented."));
}






//----------------------------------------------------------------------
//
// 000817: Discretizing a surface for plotting by lowering the degree
//         to 3 in both directions and inserting new knots, so that we
//         can plot the control polygon with increased probability of
//         plotting something that actually looks like the surface.
//
//----------------------------------------------------------------------

void lower_degree_and_subdivide(SISLSurf **srf,
				int new_knots_per_interval,
				int max_number_of_knots)
{


//   printf("\n\n\n ############# %d %d \n\n\n\n\n", new_knots_per_interval,
// 	 max_number_of_knots);

  SISLSurf *tmp=*srf;
  const double tol=1e-7;

  if (tmp->idim!=3)
    CRIT_ERR(puts("Dimension is supposed to be three."));
  
  double eeps[3];
  eeps[0]=eeps[1]=eeps[2]=tol;

  int nend[4];
  nend[0]=nend[1]=nend[2]=nend[3]=2;

  double edgeps[4*3];
  {
    int i;
    
    for (i=0; i<12; i++)
      edgeps[i]=tol;
  }

  if ((tmp->ik1>4) || (tmp->ik2>4))
    {
      int stat;
      
      s1968(tmp, eeps, nend, tmp->cuopen_1, tmp->cuopen_2, edgeps, 3, 2,
	    srf, &stat);
      if (stat!=0)
	CRIT_ERR(printf("s1968 returned %d.\n", stat));
      freeSurf(tmp);
    }

  if (new_knots_per_interval==0)
    return;

  //
  // Now, we insert new knots.
  // 1. pick out distinct knots. 2. generate new knots. 3. insert.
  //
  tmp=*srf;

  //
  // 001101: From now on we use a different refinement value for each
  //         parameter direction, and we also limit the number of knots.
  //         NB! This could be done much more precise, but we'll go for the
  //         quick and dirty solution.
  int new_knots_per_interval1 = new_knots_per_interval;

  if ((tmp->in1 - tmp->ik1 + 1) * new_knots_per_interval1 > max_number_of_knots) {
      new_knots_per_interval1 = (int)(max_number_of_knots / (tmp->in1 - tmp->ik1 + 1));
  }

  int new_knots_per_interval2 = new_knots_per_interval;
  
  if ((tmp->in2 - tmp->ik2 + 1) * new_knots_per_interval2 > max_number_of_knots) {
      new_knots_per_interval2 = (int)(max_number_of_knots / (tmp->in2 - tmp->ik2 + 1));
  }

  //
  // Allocating space for *max* number of new knots. Note that all of
  // it will be used only if we have no knots of multiplicity greater
  // than one.
  // @@@ Note! Looks like this code assumes k-multiple knots at the ends,
  //     at least it looks like no new knots are inserted between the
  //     k end-knots at either end...?!?! Should be fixed...
  //
  double *new_u=new double[(tmp->in1-tmp->ik1+1)*new_knots_per_interval1];
  if (new_u==NULL)
    CRIT_ERR(puts("Couldn't allocate memory for new knots."));
  double *new_v=new double[(tmp->in2-tmp->ik2+1)*new_knots_per_interval2];
  if (new_v==NULL)
    CRIT_ERR(puts("Couldn't allocate memory for new knots."));

  int new_u_knots, new_v_knots;
  
  // printf("New knots in direction 1:");
  {
    int i, j, p;

    for (i=0, p=0; i<tmp->in1-tmp->ik1+1; i++)
      if (fabs(tmp->et1[tmp->ik1+i]-tmp->et1[tmp->ik1-1+i]) > 1e-12)
	for (j=0; j<new_knots_per_interval1; j++, p++)
	  {
	    new_u[p]=(tmp->et1[tmp->ik1-1+i]+
		      (tmp->et1[tmp->ik1+i]-tmp->et1[tmp->ik1-1+i])/
		      (new_knots_per_interval1+1.0)*(j+1.0));
	    // printf(" %f", new_u[p]);
	  }
    // printf("\n");
    new_u_knots=p;
  }

  // printf("New knots in direction 2:");
  {
    int i, j, p;

    for (i=0, p=0; i<tmp->in2-tmp->ik2+1; i++)
      if (fabs(tmp->et2[tmp->ik2+i]-tmp->et2[tmp->ik2-1+i]) > 1e-12)
	for (j=0; j<new_knots_per_interval2; j++, p++)
	  {
	    new_v[p]=(tmp->et2[tmp->ik2-1+i]+
		      (tmp->et2[tmp->ik2+i]-tmp->et2[tmp->ik2-1+i])/
		      (new_knots_per_interval2+1.0)*(j+1.0));
	    // printf(" %f", new_v[p]);
	  }
    // printf("\n");
    new_v_knots=p;
  }

  {
    int stat;

    s1025(tmp, new_u, new_u_knots, new_v, new_v_knots, srf, &stat);
    if (stat!=0)
      CRIT_ERR(printf("s1025 returned %d.\n", stat));
  }

  freeSurf(tmp);
  delete new_u;
  delete new_v;
}






//----------------------------------------------------------------------
//
// 000828: Evaluating normals in all knot-pairs. If degree is linear in
//         both directions, what will it produce?
//
// ngrid: (Pointer to) array of 3*n1*n2 doubles containing normals. Layout:
//        n1 normals for smallest second parameter value, then n1 for
//        next second parameter value and so on. n1 and n2 are the numbers
//        of control vertices in first and second parameter direction,
//        respectively. The normals are evaluated in the parameters
//
//          ts_i = (t_{i+1} + ... + t_{i+d})/d. 
//
//----------------------------------------------------------------------

void compute_surface_normals(SISLSurf *srf, double **ngrid)
{
  if (srf->idim!=3)
    CRIT_ERR(puts("This should be a 3d surface."));

  *ngrid=new double[3*(srf->in1)*(srf->in2)];

  int i, j, left1=0, left2=0;
  
  for (i=0; i<srf->in2; i++)
    for (j=0; j<srf->in1; j++)
      {
	double par[2];
	{
	  int k;
	  
	  par[0]=0.0;
	  for (k=1; k<srf->ik1; k++)
	    par[0]+=srf->et1[j+k];
	  par[0]/=(srf->ik1-1.0);

	  par[1]=0.0;
	  for (k=1; k<srf->ik2; k++)
	    par[1]+=srf->et2[i+k];
	  par[1]/=(srf->ik2-1.0);
	}

	double der[9];
	double norm[3];
	
	int stat;
	s1421(srf, 1, par, &left1, &left2, der, norm, &stat);
	if (stat!=0)
	  {
	    printf("s1421 returned %d.\n"
		   " (""Surface is degenerate at the point, normal has "
		   "zero length"")\n", stat);
	    printf("Trying to continue anyway.\n");
	    printf("Location: %s:%d\n", __FILE__, __LINE__);
	    norm[0]=norm[1]=norm[2]=1.0;
	  }

	double l=s6length(norm, 3, &stat);
	if (stat==0)
	  CRIT_ERR(printf("s6length-vector had length 0.\n"));

	int k;
	for (k=0; k<3; k++)
	  (*ngrid)[3*(i*(srf->in1)+j) + k]=norm[k]/l;
      }
}

