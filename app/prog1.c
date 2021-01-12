#include "sisl.h"
#include <stdlib.h>
#include <stdio.h>
int main()
{
  SISLCurve *pc=0;
  double aepsco,aepsge,top[3],axispt[3],conept[3];
  double st[100],stcoef[100],*spar;
  int kstat;
  int cone_exists=0;
  int kk,kn,kdim,ki;
  int kpt,kcrv;
  SISLIntcurve **qrcrv;
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
	      spar=0;
	    }
	  if (qrcrv)
	    {
	      freeIntcrvlist(qrcrv,kcrv);
	      qrcrv=0;
	    }
	}
    }
  return 0;
}
