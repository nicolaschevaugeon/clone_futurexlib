/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef XDELTAMEMORY_H
#define XDELTAMEMORY_H

#include <map>
#include <string>
#include <vector>

#include "mpi.h"

#ifdef __GNUC__
extern "C"
{
#include <malloc.h>
}
#endif

namespace xtool
{
class xDeltaMemory
{
  public:
   xDeltaMemory(MPI_Comm world_ = MPI_COMM_WORLD);
   ~xDeltaMemory();

   int initAccu(std::string stage, bool local = false);
   inline void startAccu(int id)
   {
      switch_hook_on();
      set(m[id]);
   }
   inline void endAccu(int id)
   {
      set();
      dm[id] += (m_cur - m[id]);
      switch_hook_off();
   }

   void print();
   double get(int id);
   double get(int id, double &mi, double &mx, double &me, double &su);

  private:
   MPI_Comm world;
   int n, nb_proc, proc_id;
   std::map<std::string, int> strtoind;
   std::vector<long long int> dm;
   std::vector<long long int> m;
   long long int zero;
   long long int m_cur;
   void set();
   void set(long long int &m);
   void reduce(double &aloc, double *min, double *max, double *mean, double *su) const;
   void setOldPointer(void);
   void switch_hook_on(void);
   void switch_hook_off(void);
};

class xNoDeltaMemory
{
  public:
   xNoDeltaMemory(MPI_Comm world_ = MPI_COMM_WORLD) {}
   ~xNoDeltaMemory() = default;
   inline int initAccu(std::string stage, bool local = false) { return 0; }
   inline void startAccu(int id) {}
   inline void endAccu(int id) {}
   void print() {}
   inline double get(int id) { return 0; }
   inline double get(int id, double &mi, double &mx, double &me, double &su) { return 0.; }
};

}  // namespace xtool
#endif
