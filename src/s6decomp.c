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
 * $Id: s6decomp.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6DECOMP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6decomp(double ea[],double gx[],double eb1[],double eb2[],double eb3[],int *jstat)
#else
void s6decomp(ea,gx,eb1,eb2,eb3,jstat)
     double ea[];
     double gx[];
     double eb1[];
     double eb2[];
     double eb3[];
     int    *jstat;
#endif
/*
***********************************************************************
*
************************************************************************
*
*   PURPOSE : In dimention tree we change basis from a euclides bases
*             to eb1, eb2, eb3.
*
*
*   INPUT   : ea   - The orginale vector.
*             eb1  - the first bases vector.
*             eb2  - the first bases vector.
*             eb3  - the first bases vector.
*
*
*   
*   OUTPUT  : gx    - The new vector
*             jstat - Status variable.
*                       < 0 : error
*                       = 0 : ok
*                       > 0 : warning
*
*
*   METHOD  : 
*
*
*   REFERENCES : 
*-
*   CALLS      :
*
*   WRITTEN BY : Arne Laksaa, SI, 89-07.
*
************************************************************************
*/
{
  int kstat =0;       /* Local status variable.    */
  int ki;             /* Counter.                  */
  int n1[3];          /* Array for use in lufac.   */
  double sc[9],se[3]; /* Matrix and help vector.   */
  
  
  /* Copy new bases into local matrix.  */
  
  memcopy(sc,eb1,3,double);
  memcopy(sc+3,eb2,3,double);
  memcopy(sc+6,eb3,3,double);
  
  
  s6lufacp(sc,n1,3,&kstat);
  if (kstat < 0) goto error;
  else if (kstat > 0) goto warn1;
  
  for (ki=0; ki<3; ki++)
    {                      
      se[0] = se[1] = se[2] = DZERO;
      se[ki] = (double)1;
      
      s6lusolp(sc,se,n1,3,&kstat);
      if (kstat < 0) goto error;
      else if (kstat > 0) goto warn1;
      
      gx[ki] = s6scpr(ea,se,3);
    }
  
  /* Change of bases performed.  */
  
  *jstat = 0;
  goto out;

/* Singular equation system.  */

warn1 : *jstat = 1;
        goto out;

/* Error in subrutines.  */

error: *jstat = kstat;
        s6err("s6decomp",*jstat,0);
        goto out;

 out: ;
}

                                    

                                        
