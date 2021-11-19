/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

// ----------------------------------------------------------------------------
// HEADERS
// ----------------------------------------------------------------------------
#include "xDeltaMemory.h"

#include <sys/resource.h>
#include <sys/time.h>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

// =================================================================
// malloc,calloc,free hooks
//
// Declaration part
static void *(*old_malloc_hook)(size_t, const void *);
static void *(*old_realloc_hook)(void *, size_t, const void *);
static void (*old_free_hook)(void *, const void *);
static long long int sum = 0;
static size_t hook_counter = 0;
static void *xDeltaMemory_malloc_hook(size_t size, const void *caller);
static void *xDeltaMemory_realloc_hook(void *ptr, size_t size, const void *caller);
static void xDeltaMemory_free_hook(void *ptr, const void *caller);
static void switch_hooks_on()
{
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
   __malloc_hook = xDeltaMemory_malloc_hook;
   __realloc_hook = xDeltaMemory_realloc_hook;
   __free_hook = xDeltaMemory_free_hook;
#pragma GCC diagnostic pop
#endif
}
static void switch_hooks_off()
{
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
   __malloc_hook = old_malloc_hook;
   __realloc_hook = old_realloc_hook;
   __free_hook = old_free_hook;
#pragma GCC diagnostic pop
#endif
}

// Implementation part
static void *xDeltaMemory_malloc_hook(size_t size, const void *caller)
{
   void *result = nullptr;
   // switch from specific hook to std hook
   switch_hooks_off();

   // std call
   result = malloc(size);

   //   cout<<"mh "<<sum<<" "<<malloc_usable_size(result)<<endl;

   // counter update (take into account padding but not extra bytes assoiatted to block.
   // Thus it is a lower estimate of what is consumed but a upper estimate of what is
   // asked)
   sum += malloc_usable_size(result);

   // not clear why those 3 lines are required !
   // commented for now
   // old_malloc_hook = __malloc_hook;
   // old_realloc_hook = __realloc_hook;
   // old_free_hook = __free_hook;

   // switch back to specific hook
   switch_hooks_on();
   return result;
}
static void *xDeltaMemory_realloc_hook(void *ptr, size_t size, const void *caller)
{
   void *result = nullptr;
   // switch from specific hook to std hook
   switch_hooks_off();

   // counter update before reallocation: remove already allocated stuff to avoid counting it
   // twice. If ptr is null realloc is malloc and nothing have to be remove to counting
   size_t delta = 0;
   if (ptr) delta = malloc_usable_size(ptr);

   // std call
   result = realloc(ptr, size);

   //   cout<<"rh "<<sum<<" "<<delta<<" "<<malloc_usable_size(result)<<endl;

   // counter update: new size counted for non null size
   if (size) sum = sum - delta + malloc_usable_size(result);
   // counter update: freed size counted for null size
   else
      sum = sum - delta;

   // not clear why those 3 lines are required !
   // commented for now
   // old_malloc_hook = __malloc_hook;
   // old_realloc_hook = __realloc_hook;
   // old_free_hook = __free_hook;

   // switch back to specific hook
   switch_hooks_on();
   return result;
}
static void xDeltaMemory_free_hook(void *ptr, const void *caller)
{
   // switch from specific hook to std hook
   switch_hooks_off();

   //   cout<<"dh "<<sum<<" "<<malloc_usable_size(ptr)<<endl;

   // counter update (see remark in malloc hook)
   sum -= malloc_usable_size(ptr);

   // std call
   free(ptr);

   // not clear why those 3 lines are required !
   // commented for now
   // old_malloc_hook = __malloc_hook;
   // old_realloc_hook = __realloc_hook;
   // old_free_hook = __free_hook;

   // switch back to specific hook
   switch_hooks_on();
}
// ===============================================================================================
namespace xtool
{
std::pair<double, std::string> scaleUnit(double val)
{
   const double B = 1.;
   const double KB = 1024. * B;
   const double MB = 1024. * KB;
   const double GB = 1024. * MB;
   const double TB = 1024. * GB;

   if (val < 0) val = -val;

   if (val > TB)
   {
      // TB
      return std::make_pair(1. / TB, "TB");
   }
   else if (val > GB)
   {
      // GB
      return std::make_pair(1. / GB, "GB");
   }
   else if (val > MB)
   {
      // MB
      return std::make_pair(1. / MB, "MB");
   }
   else if (val > KB)
   {
      // KB
      return std::make_pair(1. / KB, "KB");
   }
   else
   {
      // B
      return std::make_pair(1., "B ");
   }
   throw -1234;
}

xDeltaMemory::xDeltaMemory(MPI_Comm world_) : world(world_), n(0), nb_proc(1), proc_id(0)
{
   MPI_Comm_size(world, &nb_proc);
   MPI_Comm_rank(world, &proc_id);

   // init zero
   zero = 0;

   // set hooks saving pointers (old)
   setOldPointer();
}
xDeltaMemory::~xDeltaMemory() { switch_hooks_off(); }

void xDeltaMemory::set() { m_cur = sum; }
void xDeltaMemory::set(long long int &m) { m = sum; }

int xDeltaMemory::initAccu(std::string stage, bool local)
{
   switch_hooks_off();  // remove xDeltaMemory from counts
   int i = dm.size();
   if (local)
      strtoind[stage] = -i;
   else
      strtoind[stage] = i;
   dm.push_back(zero);
   m.push_back(zero);
   switch_hooks_on();  // counts again
   return i;
}

void xDeltaMemory::print()
{
   switch_hook_off();  // remove xDeltaMemory from counts if not in a middle of a count otherwise thing will be counted
   struct rusage rusage;
   getrusage(RUSAGE_SELF, &rusage);
   double salloc = rusage.ru_maxrss * 1024.;
   assert(salloc != 0.);
   double mi, mx, me, su;
   reduce(salloc, &mi, &mx, &me, &su);
   string Heapu("Heap asked");
   string Peack("Peack  memory usage up to this point (");

   auto scale = scaleUnit(salloc);
   auto scale2 = scale;
   cout << " " << endl;
   cout
       << "==xDeltaMemory output ================================================================================================"
       << fixed << setprecision(5) << endl;
   cout << Peack + scale.second << ") : " << salloc * scale.first;
   if (!proc_id) cout << " min/max/average : " << mi * scale.first << "/" << mx * scale.first << "/" << me * scale.first;
   cout
       << endl
       << "======================================================================================================================"
       << endl;

   bool do_more_print = false;
   int j;
   std::map<std::string, int>::iterator it = strtoind.begin();
   std::map<std::string, int>::iterator itend = strtoind.end();
   if (it != itend) do_more_print = true;
   if (do_more_print)
   {
      bool got_local = false;
      if (proc_id)
      {
         cout << setfill('=') << setw(82) << " " << endl;
         cout << setfill(' ') << "| " << setw(32) << Heapu << " | " << setw(42) << "for"
              << " |" << endl;
         cout << setfill('=') << setw(82) << " " << endl;
         cout << setfill(' ') << "| " << setw(2) << "U."
              << " | " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(44) << " |" << endl;
         for (; it != itend; ++it)
         {
            j = it->second;
            if (j > -1)
            {
               double uallocj = dm[j];
               reduce(uallocj, &mi, &mx, &me, &su);
               scale = scaleUnit(uallocj);
               cout << "| " << setw(2) << scale.second << " | " << setprecision(2) << setw(12) << uallocj * scale.first << " | "
                    << setprecision(1) << setw(12) << 100 * uallocj / salloc << " | " << setw(42) << it->first << " |" << endl;
            }
            else
               got_local = true;
         }
         cout << setfill('=') << setw(82) << " " << setfill(' ') << endl;
      }
      else
      {
         cout << setfill('=') << setw(130) << " " << endl;
         cout << setfill(' ') << "| " << setw(80) << Heapu << " | " << setw(42) << "for"
              << " |" << endl;
         cout << setfill('=') << setw(130) << " " << endl;
         cout << setfill(' ') << "| " << setw(9) << "val sum"
              << " | " << setw(2) << "U."
              << " | " << setw(8) << "val min"
              << " | " << setw(8) << "val max"
              << " | " << setw(8) << "val"
              << " | " << setw(8) << "% total"
              << " | " << setw(8) << "val avg"
              << " | " << setw(8) << "% t. avg"
              << " | " << setw(44) << " |" << endl;
         for (; it != itend; ++it)
         {
            j = it->second;
            if (j > -1)
            {
               double uallocj = dm[j];
               reduce(uallocj, &mi, &mx, &me, &su);
               scale = scaleUnit(me);
               scale2 = scaleUnit(su);
               cout << "| " << fixed << setprecision(2) << setw(8) << su * scale2.first << setw(2) << scale2.second << "| "
                    << setw(2) << scale.second << " | " << setw(8) << mi * scale.first << " | " << setw(8) << mx * scale.first
                    << " | " << setw(8) << uallocj * scale.first << " | " << setprecision(1) << setw(8) << 100 * uallocj / salloc
                    << " | " << setprecision(2) << setw(8) << me * scale.first << " | " << setprecision(1) << setw(8)
                    << 100 * me / salloc;
               cout << " | " << setw(42) << it->first << " |" << endl;
            }
            else
               got_local = true;
         }
         cout << setfill('=') << setw(130) << " " << setfill(' ') << endl;
      }
      if (got_local)
      {
         cout << setfill(' ') << "| " << setw(32) << "local"
              << " |" << endl;
         cout << setfill('=') << setw(82) << " " << endl;
         cout << setfill(' ') << "| " << setw(32) << Heapu << " | " << setw(42) << "for"
              << " |" << endl;
         cout << setfill('=') << setw(82) << " " << endl;
         cout << setfill(' ') << "| " << setw(2) << "U."
              << " | " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(44) << " |" << endl;
         for (it = strtoind.begin(); it != itend; ++it)
         {
            j = it->second;
            if (j < 0)
            {
               j = -j;
               double uallocj = dm[j];
               scale = scaleUnit(uallocj);
               cout << "| " << setw(2) << scale.second << " | " << setprecision(2) << setw(12) << uallocj * scale.first << " | "
                    << setprecision(1) << setw(12) << 100 * uallocj / salloc << " | " << setw(42) << it->first << " |" << endl;
            }
         }
         cout << setfill('=') << setw(82) << " " << setfill(' ') << endl;
      }
   }

   cout
       << "==End xDeltaMemory output ============================================================================================"
       << endl;
   switch_hook_on();  // counts again
}
double xDeltaMemory::get(int id)
{
   if (id < 0) id = -id;
   return dm[id];
}
double xDeltaMemory::get(int id, double &mi, double &mx, double &me, double &su)
{
   if (id < 0) id = -id;
   double val = dm[id];
   reduce(val, &mi, &mx, &me, &su);
   return val;
}

void xDeltaMemory::reduce(double &val, double *min, double *max, double *mean, double *su) const
{
   double *data = &val;
   MPI_Reduce(data, max, 1, MPI_DOUBLE, MPI_MAX, 0, world);
   MPI_Reduce(data, min, 1, MPI_DOUBLE, MPI_MIN, 0, world);
   MPI_Reduce(data, su, 1, MPI_DOUBLE, MPI_SUM, 0, world);
   *mean = (*su) / nb_proc;
   return;
}

void xDeltaMemory::setOldPointer()
{
   if (hook_counter)
   {
      cout << "hook_counter when invoking setOldPointer should be null " << endl;
      cout << "Only one xDeltaMemory possible" << endl;
      throw -78;
   }
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
   old_malloc_hook = __malloc_hook;
   old_realloc_hook = __realloc_hook;
   old_free_hook = __free_hook;
#pragma GCC diagnostic pop
#endif
}
void xDeltaMemory::switch_hook_on()
{
#ifdef __GNUC__
   if (++hook_counter < 2)
   {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      __malloc_hook = xDeltaMemory_malloc_hook;
      __realloc_hook = xDeltaMemory_realloc_hook;
      __free_hook = xDeltaMemory_free_hook;
#pragma GCC diagnostic pop
      sum = 0;
   }
#endif
}
void xDeltaMemory::switch_hook_off()
{
#ifdef __GNUC__
   if (--hook_counter < 1)
   {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      __malloc_hook = old_malloc_hook;
      __realloc_hook = old_realloc_hook;
      __free_hook = old_free_hook;
#pragma GCC diagnostic pop
      sum = 0;
   }
#endif
}

}  // namespace xtool
