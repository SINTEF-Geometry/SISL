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

#include <stdio.h>
#include <math.h>

// #ifdef _MSC_VER
// #  include <vector>
// #else
// #  include <vector.h>
// #endif
#include <vector>

#ifdef _MSC_VER
#  define max(a, b) ((a)>(b)? (a) : (b))
#endif

#if defined(_MSC_VER) || defined(LINUX)
#  include <sys/timeb.h>
#endif

#ifndef _MSC_VER
#  include <strings.h>
#  include <unistd.h>
#else
#  include <stdlib.h>
#endif

#include <time.h>




double jon_timer;


static double jon_sec(void)
{
#ifndef _MSC_VER
#ifdef LINUX
  struct timeb jon_ts;
#else
  struct timespec jon_ts;
#endif
#else
  struct _timeb jon_ts;
#endif
  
#ifndef _MSC_VER
#ifdef LINUX
  ftime(&jon_ts);
  return (jon_ts.time + 0.001*jon_ts.millitm);
#else
  clock_gettime(CLOCK_REALTIME, &jon_ts);
  return (int)(jon_ts.tv_sec) + ((double)(jon_ts.tv_nsec)/1e9);
#endif
#else
  _ftime(&jon_ts);
  return (jon_ts.time + 0.001*jon_ts.millitm);
#endif
}


void tic(void)
{
  jon_timer=jon_sec();
}


void toc(void)
{
  printf("Time elapsed: %f\n", jon_sec()-jon_timer);
}
