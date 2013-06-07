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

#ifndef AUX2_H_INCLUDED

#include <vector>
using std::vector;
using std::pair;

#include <iostream>
using std::cout;
using std::endl;

#include <stdio.h>

#include "jonvec.h"


#ifdef MICROSOFT
#  define CRIT_ERR(stmnt) \
    printf("\nIn file %s, line %d:\n  ", __FILE__, __LINE__), \
    (stmnt), printf("\n%d\n", getchar()), exit(0)
#else
#  define CRIT_ERR(stmnt) \
    printf("\nIn file %s, line %d:\n  ", __FILE__, __LINE__), \
    (stmnt), exit(0)
#endif

#define ASSERT2(a, b) if (!(a)) CRIT_ERR(b)

#ifndef _ERRORMACROS_H
#define _ERRORMACROS_H

#include <exception>
#include <iostream>

/// Usage: REPORT;
/// Usage: MESSAGE("Message string.");
#ifdef NVERBOSE // Not verbose mode
#  define REPORT 0
#  define MESSAGE(x) 0
#else // Verbose mode
#  define REPORT cout << "\nIn file " << __FILE__ << ", line " << __LINE__ << endl
#  define MESSAGE(x) cout << "\nIn file " << __FILE__ << ", line " << __LINE__ << ": " << x << endl
#endif

/// Usage: THROW("Error message string.");
#define THROW(x) MESSAGE(x), throw std::exception()

/// Usage: ASSERT(condition)
/// Usage: ASSERT2(condition, "Error message string.")
/// Usage: ERROR_IF(condition, "Error message string.");
#ifdef NDEBUG // Not in debug mode
#  ifndef QTMODE
//   030916: There is a conflict. Hope I don't use ASSERT for anything...
#    define ASSERT(x)
#  endif
#  define ASSERT3(cond, x)
#  define ERROR_IF(cond, x)
#else // Debug mode
#  ifndef QTMODE
//   030916: There is a conflict. Hope I don't use ASSERT for anything...
#    define ASSERT(cond) if (!(cond)) THROW("Assertation \'" #cond "\' failed.")
#  endif
#  define ASSERT3(cond, x) if (!(cond)) THROW(x)
#  define ERROR_IF(cond, x) if (cond) THROW(x)
#endif


#endif // _ERRORMACROS_H






//
// 030102: There are some routines for this in 'aux1.C', but these macros
//         are more in line with what is used in the CoCreate project,
//         and in sisl. (Just hope there won't be any conflicts...)
// 030116: Ajajaj... Horrible mistake, forgot the max(..., 1)...
// 040411: Legger inn en egen ikke-sisl-relatert variant med spesifisert
//         toleranse.
//

#define DEQUALX(a, b, tol) (fabs((a)-(b))<=(tol) * std::max(std::max(fabs(a), fabs(b)), 1.0))

#define DEQUAL(a, b) (fabs((a)-(b))<=1e-12 * std::max(std::max(fabs(a), fabs(b)), 1.0))
#define NDEQUAL(a, b) (!DEQUAL(a, b))
#define DNEQUAL(a, b) (NDEQUAL(a, b))
#define DLESS(a, b) (((a)<(b)) && (NDEQUAL((a), (b))))
#define DGREATER(a, b) (((a)>(b)) && (NDEQUAL((a), (b))))
#define DLESSEQ(a, b) (((a)<(b)) || (DEQUAL((a), (b))))
#define DGREATEREQ(a, b) (((a)>(b)) || (DEQUAL((a), (b))))

//
// 030105: For some reason that I'm not completely sure of, some relatively
//         simple computations (angles in triangles, 2nd degree equations,
//         not very fancy) seem to be riddled with extreme loss of precision.
//         (Geos. project, see 'building2.C', 'split_triangle_corner_edge'.
//
#define DEQUAL2(a, b) (fabs((a)-(b))<=1e-8 * std::max(std::max(fabs(a), fabs(b)), 1.0))
// 030711: The 3-version for floats...
#define DEQUAL3(a, b) (fabs((a)-(b))<=1e-6 * std::max(std::max(fabs(a), fabs(b)), 1.0))
#define NDEQUAL2(a, b) (!DEQUAL2(a, b))
#define DNEQUAL2(a, b) (NDEQUAL2(a, b))
#define DLESSEQ2(a, b) (((a)<(b)) || (DEQUAL2((a), (b))))
#define DGREATEREQ2(a, b) (((a)>(b)) || (DEQUAL2((a), (b))))


#define DEQUAL4(a, b) (fabs((a)-(b))<=1e-14 * std::max(std::max(fabs(a), fabs(b)), 1.0))





#define SQR(a) ((a)*(a))

#ifndef PI
#  define PI 3.1415926535
#endif

#ifndef MIN
#  define MIN(a,b) ((a)<(b)? (a):(b))
#endif
#ifndef MAX
#  define MAX(a,b) ((a)>(b)? (a):(b))
#endif



  


int eps_equal(const double &a, const double &b);
int eps_less(const double &a, const double &b);
int eps_greater(const double &a, const double &b);
int eps_less_eq(const double &a, const double &b);
int eps_greater_eq(const double &a, const double &b);

void tic(void);
void toc(void);






#define AUX2_H_INCLUDED
#endif
