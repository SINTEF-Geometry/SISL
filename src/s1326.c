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


#define S1326

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
 s1326(SISLSurf *ps, int power, double *ecimp, int inarr, double **et1,
       double **et2, double **ecoef, int *ik1, int *ik2, int *in1,
       int *in2, int *numprd, int *jstat)
#else
void
    s1326(ps, power, ecimp, inarr, et1, et2, ecoef, ik1, ik2,
	  in1, in2, numprd, jstat)
       SISLSurf *ps;
       int power;
       double *ecimp;
       int inarr;
       double **et1;
       double **et2;
       double **ecoef;
       int *ik1;
       int *ik2;
       int *in1;
       int *in2;
       int *numprd;
       int *jstat;
#endif       
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute x(s,t)^i y(s,t)^j z(s,t)^k h(s,t)^(power-i-j-k)
*                 for
*              i=0,...,power. j=0,...,power-i. k=0,power-i-j.
*              If implicit coefficients are known, these are multiplied
*              with these factors to produce the coefficients of a 1D surface 
*              where the input surface are put into an implicit function.
*
* INPUT      : ps     - Pointer to the surface 
*              power  - The maximum power to be produced
*              ecimp  - Coefficients of implicit function, if
*                       no such function is known ecimp = NULL.
*              inarr  - Number of parallell arrays of implicit coeff.
*
* OUTPUT     : et1    - Knot vector in first direction
*                       belonging to the powers
*                       of the curve. The knots produced
*                       are the Bernstein knot vector, e.g.
*                       all knots has maximal multiplicity.
*              et2    - Knot vector in second direction
*                       belonging to the powers
*                       of the curve. The knots produced
*                       are the Bernstein knot vector, e.g.
*                       all knots has maximal multiplicity.
*              ecoef  - The coefficients produced stored
*                       in a long array. For each Bernstein
*                       segment the ordering is:
*                       vertices of h^power, h^(power-1)y, h^(power-2)y^2,...,
*                       x h^(power-1), x y h^(power-2),..., x^(power).
*                       If ecimp != NULL, the implicit coefficients are
*                       multiplied by the powers of the Bezier surfaces to
*                       produce the coefficients of a 1D surface where the
*                       input surface are put into an implicit function.
*              ik1     - Order in first parameter direction of the
*                        produced products
*              ik2     - Order in first parameter direction of the
*                        produced products
*              in1     - Number of vertices in first parameter direction
*                        in the products
*              in2     - Number of vertices in second parameter direction
*                        in the products
*              numprd - The number of products.
*
* METHOD     : 
* REFERENCES :
*
*-
* CALLS      : s1731    - Split B-spline surface into Bezier segments.
*              s1733    - Pick Bezier segment.
*              s6multsfs  - Multiply two 1D Bezier surfaces.
*              s6bezpowsf - Compute power of Bezier surface.
*              s6err    - Error handling routine 
*
*
* WRITTEN BY : Tor Dokken, SINTEF SI, November 1993
* REWISED BY : Vibeke Skytt, SINTEF, 06.94.
* 
* 
*********************************************************************/  
{
 int i,i1,i2,j,k,l;           /* Loop variables */
 int r, k1, k2;               /* Counters.      */
 int knvar;
  int knumb1;             /* Number of Bezier segments in first direction*/
  int knumb2;             /* Number of Bezier segments in second direction*/
  int kdum1,kdum2;        /* Dummy varaible */
  int kstat;
  int kdim;              /* Spatial dimension */
  int krat;              /* Indicator rational/nonrational curve */
  int kpos=0;
  int kkps1=ps->ik1;       /* The order of the input curve */
  int kgradps1=kkps1-1;    /* The degree of the input curve */
  int kkps2=ps->ik2;       /* The order of the input curve */
  int kkps12 = kkps1*kkps2;
  int kgradps2=kkps2-1;    /* The degree of the input curve */
  int order1=(power*(kgradps1))+1;  /* The polynomial order of the products
                            in first parameter direction*/
  int order2=(power*(kgradps2))+1;  /* The polynomial order of the products
                            in second parameter direction*/
  int order12 = order1*order2;
  int kk12;
  int *poffset=NULL;       /* Array describing offsets into the arrays
			      for powers of the x, y,, z and h coordinates */
  int order=MAX(order1,order2);/* Maximum of order1 and order 2*/
  int klength;           /* Length of arrays allocated for powers
                            of the components */
  int *variables=NULL;   /* An array containing the indicies of
                            all homogeneous combinations of the
                            curve of degree power */
  int *kpek;             /* Pointer into variables */
  int numbvar=0;           /* The number of tupples in variables */
  int kpow11,kpow12;      /* Polynomial degree of products */
  int kpow21,kpow22;      /* Polynomial degree of products */
  double sstart1,sstart2; /* Start of current Bezier segment */
  double send1,send2;     /* End of current Bezier segment */
  double *scoef1=NULL;    /* Array for storage of Bezier coefficients in 
			    mixed  x,y,z sequence*/
  double *scoef2=NULL;    /* Array for storage of Bezier coefficients in 
			    first all xsequence*/
  double *qs1,*qs2;      /* Pointers */
  double *qsx,*qsy,*qsz,*qsh;
  double *temp1=NULL;    /* Temporary storage of product of Bezier curves */
  double *temp2=NULL;    /* Temporary storage of product of Bezier curves */
  double *temp3=NULL;    /* Temporary storage of product of Bezier curves */
  double *temppek;       /* Pointer to product of Bezeir curves */
  double *Pascal=NULL;   /* Pointer to the binomial coefficients */
  double *psl1=NULL;     /* Pointer used in Pascals triangle */
  double *psl2=NULL;     /* Pointer used in Pascals triangle */
  double *powcomp=NULL;  /* Array for storage of pwoers of components */
  double *products=NULL; /* Array for storage of inal products */
  double *xyzh;          /* Pointer into products */
  double *knots1=NULL;   /* Pointer to output knots */
  double *knots2=NULL;   /* Pointer to output knots */
  SISLSurf *qsbez=NULL; /* Input curve represetned in Bernstein knot vector */
  int kimpl = (ecimp == NULL) ? 0 : 1;  /* Flage indicating whether an
					   input array of implicit coefficients
					   are given.                        */

  /* Test input.  */
  
  if (inarr <1 || inarr > 3 || (inarr > 1 && !kimpl)) goto err172;
   

/* Convert the input surface to a surface with the Bernstein knot vector */

  s1731(ps,&qsbez,&kstat);
  if(kstat<0) goto error;



  /* Calculate number of Bezier segments */

  knumb1 = qsbez->in1/qsbez->ik1;
  knumb2 = qsbez->in2/qsbez->ik2;
  kk12 = qsbez->ik1*qsbez->ik2;


/* Make number components in homogeneous polynomial of degree power in 2 or 3 dimensions */

  if(ps->idim==3)
    numbvar= (power+1)*(power+2)*(power+3)/6;

 *numprd=numbvar;


/* Allocate space for temporary arrays */
 knvar = numbvar*order12;
  temp1 = newarray(order12,DOUBLE);
  if(temp1==NULL) goto err101;
  temp2 = newarray(order12,DOUBLE);
  if(temp2==NULL) goto err101;
  temp3 = newarray(knvar,DOUBLE);
  if(temp3==NULL) goto err101;

/* Allocate space for binomial coefficients */
  
  Pascal = newarray((order+1)*(order+2)/2, DOUBLE);
  if (Pascal == NULL) goto err101;

  /*Allocate space for offsets into the arrays representing
    the powers of the x,y, z and h component */
  poffset = newarray(power+2,INT);
  if (poffset==NULL) goto error;

/* Make the actual offsets */

  poffset[0]=0;

  for(i=0;i<=power;i++)
    poffset[i+1] = poffset[i] + (i*kgradps1+1)*(i*kgradps2+1);


/* Find size of arrays for storage of powers of the components */

  klength = poffset[power+1];


/* Allocate space for powers of the components */
  powcomp = newarray((ps->idim+1)*klength,DOUBLE);
  if (powcomp==NULL) goto error;

/* Allocate space for final product. Number of Bezier segments x
   order of products x number of products, make at the same time
   the number of vertices and the order */
   *in1 = knumb1*order1;
   *in2 = knumb2*order2;
   *ik1 = order1;
   *ik2 = order2;
   if (kimpl)
      products = new0array((*in1)*(*in2)*inarr, DOUBLE);
   else
      products = newarray((*in1)*(*in2)*numbvar,DOUBLE);
  if(products==NULL) goto err101;
  knots1 = newarray((*in1)+(*ik1),DOUBLE);
  if (knots1==NULL) goto err101;
  knots2 = newarray((*in2)+(*ik2),DOUBLE);
  if (knots2==NULL) goto err101;

 for(i=0,psl2=Pascal; i<=order ; i++,psl1=psl2,psl2+=i)
     {  
       psl2[0] = (double)1;
       
       for(j=1;j<i;j++)
	 psl2[j] = psl1[j-1] + psl1[j];

       psl2[i] =(double)1;
     }
 

/* Detect if rational surface */

  if(ps->ikind==1 || ps->ikind==3) 
    krat=0;
  else
    krat=1;
  kdim = ps->idim;
             
  scoef1=newarray((qsbez->idim+1)*kk12,DOUBLE);
  if(scoef1==NULL)  goto err101;
  scoef2=newarray((qsbez->idim+1)*kk12,DOUBLE);
  if(scoef2==NULL)  goto err101;


/* Make the actual combinations */

  variables = newarray(numbvar*(ps->idim+1),INT);
  if(variables==NULL) goto err101;
  
  for(i=0,kpek=variables;i<=power;i++)
    for(j=0;j<=power-i;j++)
      for(k=0;k<=power-i-j;k++,kpek+=4)
    {
      kpek[0] = i;
      kpek[1] = j;
      kpek[2] = k;
      kpek[3] = power-i-j-k;
    }

  /* Travers all Bezier segments and compute the B-spline surface
     put into the implicit function.   */
  
  for(i2=0,xyzh=products;i2<knumb2;i2++)
    for(i1=0;i1<knumb1;i1++,xyzh+=(1-kimpl)*knvar)
     {
       /* Pick Bezier segment of the surface */

       s1733(qsbez,i1,i2,&sstart1,&send1,&sstart2,&send2,scoef1,&kstat);
       if(kstat<0) goto error;
      
       /* Store knots. */
       
       if (i1 == 0)
       {
	  if (i2 == 0)
	  {
	     for (j=0; j<order2; j++) knots2[j] = sstart2;
	  }
	  for (j=0; j<order2; j++) knots2[(i2+1)*order2+j] = send2;
       }
       if (i2 == 0)
       {
	  if (i1 == 0)
	  {
	     for (j=0; j<order1; j++) knots1[j] = sstart1;
	  }
	  for (j=0; j<order1; j++) knots1[(i1+1)*order1+j] = send1;
       }
       
       /* Order the coefficients in separate arrays */   
       
       for(j=0,qs1=scoef1,qs2=scoef2;j<kk12;j++,qs2++)
	 for(l=0;l<kdim+krat;l++,qs1++)
	   {
	     qs2[l*kk12] = *qs1;
	   }
       if(krat==0)
	 {
	   /* Fill in "1"s inn the numerator */
	   for(j=0,qs2=scoef2+kk12*qsbez->idim;
	       j<kk12;j++,qs2++)
	     *qs2 = (double)1.0;
	 }

       /* Make all homogeneous polynomial from the ps->idim coordinates
          pluss the numerator */

       for(j=0;j<ps->idim+1;j++)
	 {

	    s6bezpowsf(&scoef2[j*kkps12],kkps1,kkps2,power,Pascal,
                       &powcomp[j*klength]);

	 } 

       /* Make the homogeneous combinations of the x,y,...,h components */

       for(j=0,kpek=variables,qs1=temp3;
	   j<numbvar;j++,kpek+=(ps->idim)+1,qs1+=(1-kimpl)*order12) 
	 {
	   /* Prepare multiplication of x and y powers */


	   qsx = &powcomp[poffset[kpek[0]]];
	   qsy = &powcomp[klength+poffset[kpek[1]]];

	   /* Multiply x and y powers */

	   s6multsfs(qsx,MAX(1,(kgradps1)*kpek[0]+1),
		      MAX(1,(kgradps2)*kpek[0]+1),qsy,
		      MAX(1,(kgradps1)*kpek[1]+1),
		      MAX(1,(kgradps2)*kpek[1]+1),
		      Pascal,temp1,&kpow11,&kpow12);

	   /* Prepare multiplication of z and h powers */

	   qsz = &powcomp[2*klength+poffset[kpek[2]]];
	   qsh = &powcomp[3*klength+poffset[kpek[3]]];

	   /* Multiply z and homogenous coordinate */

	   s6multsfs(qsz,
		      MAX(1,(kgradps1)*kpek[2]+1),
		      MAX(1,(kgradps2)*kpek[2]+1),
		      qsh,
		      MAX(1,(kgradps1)*kpek[3]+1),
		      MAX(1,(kgradps2)*kpek[3]+1),
		      Pascal,temp2,&kpow21,&kpow22);
	   temppek = temp2;

	   s6multsfs(temp1, kpow11,kpow12, temppek,kpow21,kpow22,
		      Pascal, qs1,&kdum1,&kdum2);
	   
/*           fprintf(stdout,"\n %d, %d, %d, %d",kpek[0],kpek[1],kpek[2],kpek[3]);*/

	   if (kimpl)
	   {
	      /* Multiply by implicit coefficint and put the product into
		 output array.  */
	      
	      for (r=0; r<inarr; r++)
		 for (k1=0; k1<order1; k1++)
		    for (k2=0; k2<order2; k2++)
		       products[((i2*order2+k2)*(*in1) + i1*order1+k1)*inarr+r]
			  += qs1[k2*order1+k1]*ecimp[r*numbvar+j];
	   }
	 }


       if (!kimpl)
       {
       /* Copy transposed version of the result onto the output
	  array */
 
       for(j=0,qs2=xyzh,qs1=temp3;j<order12;j++,qs1=temp3+j)
	 for(l=0;l<numbvar;l++,qs2++,qs1+=order12)
	   *qs2=(*qs1);
       }
     }
  
  *jstat = 0;
  goto out;
  
  /* Error. Allocation error, not enough memory.  */
  
 err101: *jstat = -101;
  s6err("s1326",*jstat,kpos);
  goto out;
  
  /* Wrong dimension of inarr. */
  
  err172 :
     *jstat = -172;
  s6err("s13267",*jstat,kpos);
  goto out;
  
/* Error in lower level routine.  */

 error:  *jstat = kstat;
         s6err("s1326",*jstat,kpos);
         goto out;
  
  
  /* Free local used memory. */
  
 out: 
  if (temp1)  freearray(temp1);
  if (temp2)  freearray(temp2);
  if (temp3)  freearray(temp3);
  if (scoef1) freearray(scoef1);
  if (scoef2) freearray(scoef2);
  if (Pascal) freearray(Pascal);
  if (variables) freearray(variables);
  if (powcomp) freearray(powcomp);
  if (poffset) freearray(poffset);
  if (qsbez) freeSurf(qsbez);

  if(*jstat<0)
    {
      if (products) freearray(products);
      if (knots1)    freearray(knots1);
      if (knots2)    freearray(knots2);
    }
  else
    {
      *ecoef = products;
      *et1    = knots1;
      *et2    = knots2;
    }

}

