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
 * $Id: s1308.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1308

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1308(double ep[],int idim,double eimpli[],int ideg,double enorm[],int *jstat)
#else
void s1308(ep,idim,eimpli,ideg,enorm,jstat)
     double ep[];
     int    idim;
     double eimpli[];
     int    ideg;
     double enorm[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the normal vector at point in an implicit
*              function.
*
* INPUT      : ep     - The coordinates of the point
*              idim   - The dimension of the space
*              eimpli - Description of the implicit surface
*              ideg   - The degree of the implicit surface
*                        ideg=1: Plane
*                        ideg=2; Quadric surface
*              
*
* OUTPUT:      enorm  - The normal vector
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : For degree=1 the 3 first component of eimpli is the
*              normal vector. For ideg=2:                     T
*              The surface can be represented as  (X ,1)A(X,1) where A is
*              a matrix (idim+1xidim+1) representing eimpli, and X a point.
*              Assume that X is parameterized with respect to some parameters
*              X = X(s,...). If we take the derivative of X with respect
*              to s, we get a tangent in the surface.
*
*              dX/ds A X = 0. Thus the expression AX is a normal to
*              the implicit function. The idim first components of 
*              AX is the acutal normal
*
*
* REFERENCES :
*
*-
* CALLS      : 
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 3. July 1988
*
*********************************************************************
*/
{            
  int ki,kj,kl;       /* Variables in loop                           */
  int kdimp1=idim+1;  /* Dimension + 1                               */
  int kpos=0;         /* Position of error                           */
  int kstat=0;        /* Local error                                 */
  double tsum;        /* Dummy variable                              */
  
  if (ideg != 1 && ideg !=2 && ideg != 1001) goto err175; 
  
  if (ideg == 1)
    {
      /*  First degree implicit surface normal vector is eimpli[0:idim-1] */
      memcopy(enorm,eimpli,idim,DOUBLE);
    }
  else if (ideg==2)
    {
      
      /* Calculate the matrix product */
      
      for (ki=0;ki<idim;ki++)
        {
	  tsum = eimpli[idim*kdimp1+ki];
	  for (kj=0,kl=ki ; kj<idim ; kj++,kl+=kdimp1)
            {
	      tsum +=(eimpli[kl]*ep[kj]);
            }
	  enorm[ki] = tsum;
        }
    }
  else if (ideg==1001)
    {  
      /*  Torus surface */
      
      double *scentr;  /* The center of the torus */
      double *snorm;   /* The normal of the torus symmetry plane */
      double tbigr;    /* The big radius of the torus */ 
      double tsmalr;   /* The small radius of the torus */
      double sdum1[3]; /* Temporary storage for point */
      double sdum2[3]; /* Temporary storage for point */
      double tproj;    /* Projection of vector onto snorm */
      
      
      scentr = eimpli;
      snorm  = eimpli+3;
      tbigr  = *(eimpli+6);
      tsmalr = *(eimpli+7);
      
      /*  Find projection of vector from torus center on to torus axis */
      s6diff(ep,scentr,3,sdum1);
      tproj = s6scpr(sdum1,snorm,3);
      
      /*  Project vector from torus center to ep onto torus plane */
      for (ki=0;ki<3;ki++)
        sdum2[ki] = sdum1[ki] - tproj*snorm[ki];
      (void)s6norm(sdum2,3,sdum2,&kstat);
      if (kstat<0) goto error;
      
      /*  Find vector from torus circle to ep */
      for (ki=0;ki<3;ki++)
        sdum1[ki] = sdum1[ki] - tbigr*sdum2[ki];
      
      /*  Normalize this vector */
      (void)s6norm(sdum1,3,enorm,&kstat);
      if (kstat<0) goto error;
    }
  
  *jstat = 0;
  goto out;
  
  /* IDEG NOT 1 OR 2 */
 err175:
  *jstat = -175;
  s6err("s1308",*jstat,kpos);
  goto out;
  
  
  /* Error in lower leve function */
 error:
  *jstat = kstat;
  s6err("s1308",*jstat,kpos);
  goto out;
  
 out:
  return;
}
                                                                              
