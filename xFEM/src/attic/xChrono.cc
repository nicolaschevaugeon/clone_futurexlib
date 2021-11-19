/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef CHRONO_CC
#define CHRONO_CC

// ----------------------------------------------------------------------------
// HEADERS
// ----------------------------------------------------------------------------
#include "xChrono.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>

//#include <sys/times.h>
#include <unistd.h>

using namespace std;


namespace xfem
{
// ----------------------------------------------------------------------------
// Chrono::Chrono()
// ----------------------------------------------------------------------------
Chrono::Chrono()
{
    struct timezone tz;
    gettimeofday( &start, &tz );
}

// ----------------------------------------------------------------------------
// Chrono::top()
// ----------------------------------------------------------------------------
double Chrono::top() const
{
    struct timezone tz;
    timeval current;
    gettimeofday( &current, &tz );

    double scale = 1.0e-6;
    double val = static_cast < double >( current.tv_sec - start.tv_sec ) +
                 scale * static_cast < double >( current.tv_usec - start.tv_usec );
    return val;
}



double Chrono::top(std::string stage) const
{
    struct timezone tz;
    timeval current;
    gettimeofday( &current, &tz );

    double scale = 1.0e-6;
    double val = static_cast < double >( current.tv_sec - start.tv_sec ) +
                 scale * static_cast < double >( current.tv_usec - start.tv_usec );


    return val-times[stage];
}

void Chrono::setStage(std::string stage)
{
    struct timezone tz;
    timeval current;
    gettimeofday( &current, &tz );

    double scale = 1.0e-6;
    double val = static_cast < double >( current.tv_sec - start.tv_sec ) +
                 scale * static_cast < double >( current.tv_usec - start.tv_usec );

    times[stage] = val;
}

} // end of namespace
#endif
// == END OF FILE =============================================================
