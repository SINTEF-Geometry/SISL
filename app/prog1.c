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

#include "sisl.h"
#include <stdlib.h>
#include <stdio.h>
int main()
{
  SISLCurve *pc=NULL;
  double aepsco,aepsge,top[3],axispt[3],conept[3];
  double st[100],stcoef[100],*spar=NULL;
  int kstat;
  int cone_exists=0;
  int kk,kn,kdim,ki;
  int kpt,kcrv;
  SISLIntcurve **qrcrv=NULL;
  char ksvar[100];
  kdim=3;
  aepsge=0.001; /* geometric tolerance */
  aepsco=0.000001; /* computational tolerance This parameter is included from historical reasons
                      and no longer used  */

  ksvar[0] = '0';  /* arbitrary character */
  while (ksvar[0] != 'q')
    {
      printf("\n cu - define a new B-spline curve");
      printf("\n co - define a new cone");
      printf("\n i - intersect the B-spline curve with the cone");
      printf("\n q - quit");
      printf("\n> ");
      scanf("%s",ksvar);

      if (ksvar[0] == 'c' && ksvar[1] == 'u')
	{
	  printf("\n Give number of vertices, order of curve: ");
	  scanf("%d %d", &kn, &kk);
	  printf("Give knots values in ascending order: \n");
	  for (ki=0;ki<kn+kk;ki++)
	    {
	      scanf("%lf",&st[ki]);
	    }
	  printf("Give vertices \n");
	  for (ki=0;ki<kn*kdim;ki++)
	    {
	      scanf("%lf",&stcoef[ki]);
	    }
	  if(pc) freeCurve(pc);
	  pc = newCurve(kn,kk,st,stcoef,1,kdim,1);
	}
      else if (ksvar[0] == 'c' && ksvar[1] == 'o')
	{
	  printf("\n Give top point: ");
	  scanf("%lf %lf %lf",&top[0],&top[1],&top[2]);
	  printf("\n Give a point on the axis: ");
	  scanf("%lf %lf %lf",&axispt[0],&axispt[1],&axispt[2]);
	  printf("\n Give a point on the cone surface: ");
	  scanf("%lf %lf %lf",&conept[0],&conept[1],&conept[2]);
	  cone_exists=1;
	}
      else if (ksvar[0] == 'i' && cone_exists && pc)
	{
	  s1373(pc,top,axispt,conept,kdim,aepsco,aepsge,
		&kpt,&spar,&kcrv,&qrcrv,&kstat);
	  printf("\n kstat %d",kstat);
	  printf("\n kpt %d",kpt);
	  printf("\n kcrv %d",kcrv);
	  for (ki=0;ki<kpt;ki++)
	    {
	      printf("\nIntersection point %lf",spar[ki]);
	    }
	  if (spar)
	    {
	      free (spar);
	      spar=NULL;
	    }
	  if (qrcrv)
	    {
	      freeIntcrvlist(qrcrv,kcrv);
	      qrcrv=NULL;
	    }
	}
    }
  return 0;
}
