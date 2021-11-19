/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef CHRONO_HH
#define CHRONO_HH

// ----------------------------------------------------------------------------
// HEADERS
// ----------------------------------------------------------------------------
#include <sys/times.h>
#include <sys/time.h>
#include <string>
#include <map>
#include <vector>
#include "mpi.h"

//#define USE_FINE_THREAD 1

#ifdef USE_FINE_THREAD
extern "C" {
#include <time.h>
}
#endif

namespace xfem
{
// ----------------------------------------------------------------------------
// Chrono
// ----------------------------------------------------------------------------
/*! \class Chrono
    \ingroup Utilities
    \brief Process timing measurment.

    Chrono is a small class that stores the current time when an instance is
    created and provides a function returning the number of seconds elapsed
    since the instance creation.

  
*/
class Chrono
{
  public:
    Chrono();
    ~Chrono()= default;;

    double top() const; //!< Return elapsed time since object creation.
    double top(std::string stage) const;//!< Return elapsed time since mark has been set.
    void setStage(std::string stage);//!< set mark.
      
  private:
    timeval start;      //!< Creation time of the object.
    mutable std::map<std::string, double> times;
};


} // end of namespace
#endif
// == END OF FILE =============================================================
