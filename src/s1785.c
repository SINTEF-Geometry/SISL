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
 * $Id: s1785.c,v 1.3 2005-02-28 09:04:49 afr Exp $
 *
 */
#define S1785

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
/* static void s1785_s9eval(double [], double [], double [],double [], int, int *); */
#else
/* static void s1785_s9eval(); */
#endif

#if defined(SISLNEEDPROTOTYPES)
void
        s1785(SISLCurve *pcurve,SISLSurf *psurf,double aepsge,
	   double epar1[],double epar2[],int icur,int *jstat)
#else
   void s1785(pcurve,psurf,aepsge,epar1,epar2,icur,jstat)
      SISLCurve  *pcurve;
      SISLSurf   *psurf;
      double     aepsge;
      double     epar1[];
      double     epar2[];
      int        icur;
      int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Test if a curve and a surface coincide beetween
*              two intersection points.
*              This function is used when the first derivative
*              for the curve and the surface is matching in both
*              intersection points.
*
*
* INPUT      : pcurve   - Pointer to the curve.
*              psurf    - Pointer to the surface.
*              aepsge   - Geometry resolution.
*              epar1[3] - Parameter values for the first intersection point.
*              epar2[3] - Parameter values for the second intersection point.
*              icur     -         = 1 , epar[0] is the parameter value
*                                       for the curve and epar[1],epar[2] is
*                                       the parameter values for the surface.
*                                 = 0 , epar[2] is the parameter value
*                                       for the curve and epar[0],epar[1] is
*                                       the parameter values for the surface.
*
*
*
* OUTPUT     :  jstat   - status messages  
*                                = 1   : Coincidence of curve and surface
*                                = 0   : No coincidence.
*                                < 0   : Error.
*
*
* METHOD     : March along the curve and iterate down to the surface for
*              each midpoint and endpoint of each step. The steplenght
*              is computed from the curvature of the curve and the surface
*              in the direction of the curve, and from the knot vector
*              of the curve. If any knot lines of the surface are crossed,
*              the marching are drawn back to the knot line, and a new
*              step is initiated. NB. To speed up the marching, the step
*              length estimate from the surface is commented out.
*
*
* REFERENCES :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. July 1989
* REWISED BY : Vibeke Skytt, SI, Oslo, Norway. Oct. 1990
*
*********************************************************************
*/
{
	 int kstat;          /* Status variable                                 */
      int ki;             /* Counter.                                        */
      int kleftc=0;       /* Left indicator for point calculation            */
      int kleft1=0;       /* Left indicator for point calculation in 1. par.
			     direction of surface.                           */
      int kleft2=0;       /* Left indicator for point calculation in 2. par dir.*/
      int kleft1prev,kleft2prev;  /* Previous left indicators of surface.    */
      int khelp;          /* Help index of knot vector.                      */
      int kn;             /* The number of B-splines, i.e., the dimension of
			     the spline space associated with the knot
			     vector.                                         */
      int kk;             /* The polynomial order of the curve.              */
      int kk1,kk2,kn1,kn2;/* Orders and nu,ber of vertices of surface        */
      int kdimc;          /* The dimension of the space in which the curve
			     lies. Equivalently, the number of components
			     of each B-spline coefficient.                   */
      int kdims;          /* Dimension of space where the surface lies       */
      int kpos=0;         /* Position of error                               */
      int kderc=2;        /* Number of derivatives to be claculated on curve */
      int kders=1;        /* Number of derivatives to be calculated on surface
			     If step lenght is to be generated from surface,
			     kders must be equal to 2.                       */
      int kdum;           /* Temporary variable                              */
      int kpar;           /* Parameter value of constant parameter curve.    */
      double tclose1,tclose2;  /* Parameter values of closest point between curves. */
      double snorm[3];    /* Normal vector of surface                        */
      double s3dinf1[10]; /* Pointer to storage for point info of curve
			     (10 dobules prpoint when idim=3, 7 when idim=3) */
      double *st;         /* Pointer to the first element of the knot vector
			     of the curve. The knot vector has [kn+kk]
			     elements.                                       */
      double *st1;        /* First knot direction of surface                 */
      double *st2;        /* Second knot direction of surface                */
      double sfirst[2];   /* Start parameter par in surface                  */
      double slast[2];    /* End parameter par in surface                    */
      double tfirst;      /* Fist parameter on curve                         */
      double tend;        /* Last parameter on curve                         */
      double sderc[9];    /* Position, first and second derivative of curve  */
      double sders[18];   /* Position, first and second derivatives of surface */
      double tx,tx1,tx2;  /* Parameter value */
      double tcstep;      /* Step length based on curvature of objects.   */
      double tstep;       /* Final step length     */
      double tmaxinc;     /* Maximal increment in parameter value along curve*/
      double tlengthend;  /* Length of 1st derivative at end of segment */
      double tincre;      /* Parameter value increment */
      double tsmax,tcmax; /* Local maximal step length based of boxsizes of objects */
      double tdist;       /* Distance */
      double tref;        /* Referance value in equality test.               */
      double sstart[2];   /* Lower boundary of parameter intervals */
      double send[2];     /* Upper bounadry of parameter intervals */
      double snext[3];    /* Existing iteration point on  surface            */
      double spos[3];     /* New iteration  point on surface                 */
      double snext2[2];   /* Help parameter values.                          */
      SISLPoint *qpoint=SISL_NULL;
      SISLCurve *qc = SISL_NULL;   /* Constant parameter curve.                       */
      
      *jstat = 0;
      
      /* Make maximal step length based on box-size of surface */
      
      sh1992su(psurf,0,aepsge,&kstat);
      if (kstat < 0) goto error;
      
      tsmax = MAX(psurf->pbox->e2max[0][0] - psurf->pbox->e2min[0][0],
		  psurf->pbox->e2max[0][1] - psurf->pbox->e2min[0][1]);
      tsmax = MAX(tsmax,psurf->pbox->e2max[0][2] - psurf->pbox->e2min[0][2]);
      
      /* Make maximal step length based on box-size of curve */
      
      sh1992cu(pcurve,0,aepsge,&kstat);
      if (kstat < 0) goto error;
      
      tcmax = MAX(pcurve->pbox->e2max[0][0] - pcurve->pbox->e2min[0][0],
		  pcurve->pbox->e2max[0][1] - pcurve->pbox->e2min[0][1]);
      tcmax = MAX(tcmax,pcurve->pbox->e2max[0][2] - pcurve->pbox->e2min[0][2]);
      
      /* Copy curve attributes to local parameters.  */
      
      kdimc = pcurve -> idim;
      kk    = pcurve -> ik;
      kn    = pcurve -> in;
      st    = pcurve -> et;
      
      /* Copy surface attributes to local parameters.  */
      
      kdims = psurf -> idim;
      kk1   = psurf -> ik1;
      kk2   = psurf -> ik2;
      kn1   = psurf -> in1;
      kn2   = psurf -> in2;
      st1   = psurf -> et1;
      st2   = psurf -> et2;
      
      /* Set reference value.  */
      
      tref = MAX(st[kn]-st[kk-1],MAX(st1[kn1]-st1[kk1-1],st2[kn2]-st2[kk2-1]));

      /* Check that dimensions are 3 */
      
      if (kdimc != 3 || kdims != 3) goto err105;
      
      sstart[0] = st1[kk1-1];
      sstart[1] = st2[kk2-1];
      send[0] = st1[kn1];
      send[1] = st2[kn2];
      
      /* Copy interval description into local variables */
      
      if (icur ==1)
	 if ( epar1[0]<epar2[0] )
	 {
	    sfirst[0] = epar1[1];
	    sfirst[1] = epar1[2];
	    slast[0]  = epar2[1];
	    slast[1]  = epar2[2];
	    tfirst    = epar1[0];
	    tend      = epar2[0];
	 }
	 else
	 {
	    sfirst[0] = epar2[1];
	    sfirst[1] = epar2[2];
	    slast[0]  = epar1[1];
	    slast[1]  = epar1[2];
	    tfirst    = epar2[0];
	    tend      = epar1[0];
	 }
      else
	 if ( epar1[2]<epar2[2] )
	 {
	    sfirst[0] = epar1[0];
	    sfirst[1] = epar1[1];
	    slast[0]  = epar2[0];
	    slast[1]  = epar2[1];
	    tfirst    = epar1[2];
	    tend      = epar2[2];
	 }
	 else
	 {
	    sfirst[0] = epar2[0];
	    sfirst[1] = epar2[1];
	    slast[0]  = epar1[0];
	    slast[1]  = epar1[1];
	    tfirst    = epar2[2];
	    tend      = epar1[2];
	 }
      
      /* To make sure we do not start outside or end outside the curve we
	 truncate tstart and tend to the knot interval of the curve */
      
      tfirst = MAX(tfirst,st[kk-1]);
      tend   = MIN(tend,st[kn]);
      if (DEQUAL(tfirst,tend)) goto out;
      
      /* Set start point of iteration on surface */
      
      spos[0] = sfirst[0];
      spos[1] = sfirst[1];
      
      /* Store knot values at start of curve */
      
      tx2 = tfirst;
      kdum = MAX(kk1,kk2);
      kdum = MAX(kdum,kk);
      tmaxinc = (tend-tfirst)/(kdum*kdum);      
      
      /* Make start point of curve  */
      
      s1221(pcurve,kderc,tx2,&kleftc,sderc,&kstat);
      if (kstat<0) goto error;
      
      /* Make start point of surface.  */
      
      s1421(psurf,kders,spos,&kleft1,&kleft2,sders,snorm,&kstat);
      if (kstat < 0) goto error;
      
      /* While end not reached */
      
      while (tx2 < tend)
      {
	 /* Save parameters of previous step.   */
	 
	 tx1 = tx2;
	 snext[0] = spos[0];
	 snext[1] = spos[1];
	 kleft1prev = kleft1;
	 kleft2prev = kleft2;
	 
	 /* Calculate unit tangent and radius of curvature of curve. */
	 
	 s1307(sderc,kdimc,s3dinf1,&kstat);
	 if (kstat<0) goto error;
	 
	 /* Calculate step length based on curvature */
	 
	 tcstep = s1311(s3dinf1[3*kdimc],aepsge,tsmax,&kstat);
	 if (kstat<0) goto error;
	 
	 /* Remember length of start tangent, end of zero segment */
	 
	 tlengthend = s6length(sderc+kdimc,kdimc,&kstat);
	 if (kstat<0) goto error;     
	 
	 /* Compute position, first and second derivative of the curve in the
	    surface going through the evaluated point in this point. 
	 
		      s1785_s9eval(sders,snorm,sderc+kdimc,sder2,kdims,&kstat);
	 if (kstat < 0) goto error;
	 
	  Calculate unit tangent and radius of curvature of curve in surface. 
	 
	 s1307(sder2,kdims,s3dinf2,&kstat);
	 if (kstat<0) goto error;
	 
	  Calculate step length based on curvature 
	 
	 tsstep = s1311(s3dinf2[3*kdims],aepsge,tcmax,&kstat);
	 if (kstat<0) goto error;
	 
	  Compute minimum step length.  
	 
	 tstep = MIN(tcstep,tsstep);  */

	 tstep = tcstep;
	 
	 /* Find candidate end point, make sure that no breaks in tangent or
	    curvature exists between start and endpoints of the segment      */
	 
	 /* Make step length equal to resolution if the length is zero */
	 
	 /* Find parameter value of candidate end point of segment */
	 
	 if (DEQUAL(tlengthend,DZERO))
	    tincre = REL_PAR_RES;
	 else
	    tincre = tstep/tlengthend;
	 
	 /* Make sure that we don't pass any knots of the curve. */
	 
	 tincre = MIN(tincre,tmaxinc);
	 tx2 = MIN(tx1 + tincre,st[kleftc+1]);
	 
	 for (ki=0, tx=(tx1+tx2)/(double)2.0; ki<2; ki++, tx=tx2)
	 {
	    if (tx >= tend) break;
	    
	    /* Make point sderc at curve at tx */
	    
	    s1221(pcurve,kderc,tx,&kleftc,sderc,&kstat);
	    if (kstat<0) goto error;
	    
	    /* Find closest point on surface to sderc */
	    
	    qpoint = newPoint(sderc,kdimc,0);
	    if (qpoint==SISL_NULL) goto err101;
	    
	    snext2[0] = snext[0];
	    snext2[1] = snext[1];
	    s1773(qpoint,psurf,aepsge,sstart,send,snext2,spos,&kstat);
	    if(kstat<0) goto error;
	    
	    freePoint(qpoint);  qpoint = SISL_NULL;
	    
	    /* Calculate point and derivatives in surface */
	    
	    s1421(psurf,kders,spos,&kleft1,&kleft2,sders,snorm,&kstat);
	    if (kstat<0) goto error;
	    
	    /* Check if point on curve and surface are within positional and
	       angular tolerances */
	    
	    tdist = s6dist(sderc,sders,kdimc);
	    
	    if (tdist>aepsge)
	    {
	       /* Points not within tolerances, curve and surface do not
		  coincide */
	       goto war01;
	    }
	    
	    /* Check if any parameter lines of the surface is crossed in the 1. 
	       parameter direction.  */
	    
	    /* changed by Michael Metzger, Feb 1993 */
	    /* for (khelp=kleft1prev-1; DEQUAL(st1[khelp],st1[kleft1prev]); khelp--); */
	    for (khelp=kleft1prev-1; khelp >= 0 && DEQUAL(st1[khelp],st1[kleft1prev]); khelp--);
	    if (kleft1 != kleft1prev && 
		((DNEQUAL(spos[0]+tref,st1[kleft1]+tref) &&
		 DNEQUAL(snext[0]+tref,st1[kleft1]+tref)) || 
		 kleft1 != kleft1prev+1) &&
		((DNEQUAL(snext[0]+tref,st1[kleft1prev]+tref) &&
		 DNEQUAL(spos[0]+tref,st1[kleft1prev]+tref)) || kleft1 != khelp))
	    {
	       /* At least one parameter line is crossed. Fetch the constant parameter
		  curve at the closest parameter line in the direction of the marching. */
	       
	       if (kleft1 > kleft1prev) kpar = kleft1prev + 1;
	       else if (snext[0] != st1[kleft1prev]) kpar = kleft1prev;
	       else kpar = khelp;
	       
	       /* Pick constant parameter curve.   */
	       
	       s1437(psurf,st1[kpar],&qc,&kstat);
	       if (kstat < 0) goto error;
	       
	       /* Find the closest point between the input curve and the constant
		  parameter curve.    */
	       
	       s1770(pcurve,qc,aepsge,tx1,st2[kk2-1],tx,st2[kn2],(tx1+tx)/(double)2.0,
		     st2[kleft2],&tclose1,&tclose2,&kstat);
	       if (kstat < 0) goto error;
	       
	       /* Set new parameter values to the iteration.  */
	       
	       spos[0] = st1[kpar];
	       spos[1] = tclose2;
	       tx2 = tclose1;
	       
              /* Test midpoint of reduced step. First evaluate curve in midpoint. */
	       
	       tx = (tx1 + tx2)/(double)2.0;
	       s1221(pcurve,kderc,tx,&kleftc,sderc,&kstat);
	       if (kstat<0) goto error;
	       
	       /* Find closest point on surface to sderc */
	       
	       qpoint = newPoint(sderc,kdimc,0);
	       if (qpoint==SISL_NULL) goto err101;
	       
	       snext2[0] = snext[0];
	       snext2[1] = snext[1];
	       s1773(qpoint,psurf,aepsge,sstart,send,snext2,snext2,&kstat);
	       if(kstat<0) goto error;
	       
	       freePoint(qpoint);  qpoint = SISL_NULL;
	       
	       /* Calculate point and derivatives in surface */
	       
	       s1421(psurf,kders,snext2,&kleft1,&kleft2,sders,snorm,&kstat);
	       if (kstat<0) goto error;
	       
	       /* Check if point on curve and surface are within positional and
		  angular tolerances */
	       
	       tdist = s6dist(sderc,sders,kdimc);
	       
	       if (tdist>aepsge)
	       {
		  /* Points not within tolerances, curve and surface do not
		     coincide */
		  goto war01;
	       }

	       /* Calculate point and derivatives in the curve in the endpoint of the step. */
	       
	       s1221(pcurve,kderc,tx2,&kleftc,sderc,&kstat);
	       if (kstat<0) goto error;
	       
	       /* Calculate point and derivatives in the surface.  */
	       
	       s1421(psurf,kders,spos,&kleft1,&kleft2,sders,snorm,&kstat);
	       if (kstat<0) goto error;
	       
	       /* Check if point on curve and surface are within positional and
		  angular tolerances */
	       
	       tdist = s6dist(sderc,sders,kdimc);
	       
	       if (tdist>aepsge)
	       {
		  /* Points not within tolerances, curve and surface do not
		     coincide */
		  goto war01;
	       }
	       
	       /* Mark that a new step is to be initiated.  */
	       
	       ki = 2;
	       
	       /* Free constant parameter curve.  */
	       
	       if (qc != SISL_NULL) freeCurve(qc);  qc = SISL_NULL;
	    }
	    
	    /* Check if any parameter lines of the surface is crossed in the 2. 
	       parameter direction.  */
	    
	    /* changed by Michael Metzger, Feb 1993 */
	    /* for (khelp=kleft2prev-1; DEQUAL(st2[khelp],st2[kleft2prev]); khelp--); */
	    for (khelp=kleft2prev-1; khelp >= 0 && DEQUAL(st2[khelp],st2[kleft2prev]); khelp--);
	    if (kleft2 != kleft2prev && 
		((DNEQUAL(spos[1]+tref,st2[kleft2]+tref) &&
		 DNEQUAL(snext[1]+tref,st2[kleft2]+tref)) || 
		 kleft2 != kleft2prev+1) &&
		((DNEQUAL(snext[1]+tref,st2[kleft2prev]+tref) &&
		 DNEQUAL(spos[1]+tref,st2[kleft2prev]+tref)) ||
		 kleft2 != khelp))
	    {
	       /* At least one parameter line is crossed. Fetch the constant parameter
		  curve at the closest parameter line in the direction of the marching. */
	       
	       if (kleft2 > kleft2prev) kpar = kleft2prev + 1;
	       else if (snext[1] != st2[kleft2prev]) kpar = kleft2prev;
	       else kpar = khelp;
	       
	       /* Pick constant parameter curve.   */
	       
	       s1436(psurf,st2[kpar],&qc,&kstat);
	       if (kstat < 0) goto error;
	       
	       /* Find the closest point between the input curve and the constant
		  parameter curve.    */
	       
	       s1770(pcurve,qc,aepsge,tx1,st1[kk1-1],tx,st1[kn1],(tx1+tx)/(double)2.0,
		     st1[kleft1],&tclose1,&tclose2,&kstat);
	       if (kstat < 0) goto error;
	       
	       /* Set new parameter values to the iteration.  */
	       
	       spos[0] = tclose2;
	       spos[1] = st2[kpar];
	       tx2 = tclose1;
	       
	       /* Test midpoint of reduced step. First evaluate curve in midpoint. */
	       
	       tx = (tx1 + tx2)/(double)2.0;
	       s1221(pcurve,kderc,tx,&kleftc,sderc,&kstat);
	       if (kstat<0) goto error;
	       
	       /* Find closest point on surface to sderc */
	       
	       qpoint = newPoint(sderc,kdimc,0);
	       if (qpoint==SISL_NULL) goto err101;
	       
	       snext2[0] = snext[0];
	       snext2[1] = snext[1];
	       s1773(qpoint,psurf,aepsge,sstart,send,snext2,snext2,&kstat);
	       if(kstat<0) goto error;
	       
	       freePoint(qpoint);  qpoint = SISL_NULL;
	       
	       /* Calculate point and derivatives in surface */
	       
	       s1421(psurf,kders,snext2,&kleft1,&kleft2,sders,snorm,&kstat);
	       if (kstat<0) goto error;
	       
	       /* Check if point on curve and surface are within positional and
		  angular tolerances */
	       
	       tdist = s6dist(sderc,sders,kdimc);
	       
	       if (tdist>aepsge)
	       {
		  /* Points not within tolerances, curve and surface do not
		     coincide */
		  goto war01;
	       }

	       /* Calculate point and derivatives in the curve.    */
	       
	       s1221(pcurve,kderc,tx2,&kleftc,sderc,&kstat);
	       if (kstat<0) goto error;
	       
	       /* Calculate point and derivatives in the surface.  */
	       
	       s1421(psurf,kders,spos,&kleft1,&kleft2,sders,snorm,&kstat);
	       if (kstat<0) goto error;
	       
	       /* Check if point on curve and surface are within positional and
		  angular tolerances */
	       
	       tdist = s6dist(sderc,sders,kdimc);
	       
	       if (tdist>aepsge)
	       {
		  /* Points not within tolerances, curve and surface do not
		     coincide */
		  goto war01;
	       }
	       
	       /* Mark that a new step is to be initiated.  */
	       
	       ki = 2;
	       
	       /* Free constant parameter curve.  */
	       
	       if (qc != SISL_NULL) freeCurve(qc);  qc = SISL_NULL;
	    }
	 }
      }
      
      /* Curves within tolerance. Test on whether the start- and
	 endpoint of the curve are identical.                      */
      
      *jstat = (DEQUAL(tfirst,tend)) ? 0 : 1;
      goto out;
      
      /* Curve and surface not within tolerance */
      war01: *jstat = 0;
      goto out;
      
      /* Error in memory allocation */
      
      err101: *jstat = -101;
      s6err("S1785",*jstat,kpos);
      goto out;
      
      /* Error in input, dimension not equal to 2 or 3 */
      
      err105: *jstat = -105;
      s6err("S1785",*jstat,kpos);
      goto out;
      
      /* Error in lower level function */
      
      error:  *jstat = kstat;
      s6err("S1785",*jstat,kpos);
      goto out;
      
      
      out:
	 
	 return;
      }          
	 
#if 0
	 
#if defined(SISLNEEDPROTOTYPES)
   static
      void
            s1785_s9eval(double eders[],double enorms[],double etanc[], 
		   double ederc[],int idim, int *jstat)
#else
static void s1785_s9eval(eders,enorms,etanc,ederc,idim,jstat)
	double eders[];
	double enorms[];
	double etanc[];
	double ederc[];
	int idim;
	int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Compute the position, first and second derivative of 
*              a curve going through a given point of the surface when
*              the 0-2'th derivatives of the surface is given. The
*              tangent of the wanted curve is parallel to the projection
*              of a given vector into the tangent plane of the surface.
*
*
* INPUT      : eders    - 0-2'th derivatives of the surface. Dimension
*                         is 6*idim.
*              enorms   - Normal vector of the surface. Dimension is idim.
*              etanc    - Vector to be projected into the tangent plane
*                         of the surface. Dimension is idim.
*              idim     - Dimension of geometry space.
*
*
*
* OUTPUT     : ederc   - 0-2'th derivative of the curve in the surface.
*              jstat   - status messages  
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SI, Oslo, Norway. Oct. 1990
*
*********************************************************************
*/
{
   int kstat = 0;         /* Status variable.  */
   int ki;                /* Counter.          */
   int ksign = 1;         /* Parameter used in s6findfac.     */
   double tfac1,tfac2,tfac3;  /* Factors found by s6findfac.  */
   
   /* Copy position of surface to output array.   */
   
   memcopy(ederc,eders,idim,DOUBLE);
   
   /* Compute the factors used to express etanc by the derivatives and normal
      of the surface.  */
   
   s6findfac(eders+idim,eders+2*idim,enorms,etanc,idim,ksign,&tfac1,&tfac2,
	     &tfac3,&kstat);
   if (kstat < 0) goto error;
   
   /* Compute first and second derivative of the curve in the surface.  */
   
   for (ki=0; ki<idim; ki++)
   {
      ederc[idim+ki] = tfac1*eders[idim+ki] + tfac2*eders[2*idim+ki];
      ederc[2*idim+ki] = tfac1*tfac1*eders[3*idim+ki] 
	 + (double)2.0*tfac1*tfac2*eders[4*idim+ki] + tfac2*tfac2*eders[5*idim+ki];
   }
   
   *jstat = 0;
   goto out;
   
   /* Error in lower level routine.  */
   
   error:
      *jstat = kstat;
   goto out;
   
   out:
      return;
}

#endif /* if 0 */
