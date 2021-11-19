/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

// ----------------------------------------------------------------------------
// HEADERS
// ----------------------------------------------------------------------------
#include "xDeltaTime.h"

#include <unistd.h>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

namespace xtool
{
void xDeltaTimeData::operator=(const xDeltaTimeData &rhs)
{
#ifdef USE_FINE_THREAD
   this->ncpu = rhs.ncpu;
#endif
   this->cpu = rhs.cpu;
   this->sys = rhs.sys;
   this->elapse = rhs.elapse;
   this->elapsem = rhs.elapsem;
}

xDeltaTimeData xDeltaTimeData::operator+(const xDeltaTimeData &rhs)
{
   xDeltaTimeData res;
   res.sys = this->sys + rhs.sys;
   res.elapse = this->elapse + rhs.elapse;
   res.elapsem = this->elapsem + rhs.elapsem;
   res.cpu = this->cpu + rhs.cpu;
#ifdef USE_FINE_THREAD
   res.ncpu = this->ncpu + rhs.ncpu;
#endif
   return res;
}

xDeltaTimeData xDeltaTimeData::operator-(const xDeltaTimeData &rhs)
{
   xDeltaTimeData res;
   res.sys = this->sys - rhs.sys;
   res.elapse = this->elapse - rhs.elapse;
   res.elapsem = this->elapsem - rhs.elapsem;
   res.cpu = this->cpu - rhs.cpu;
#ifdef USE_FINE_THREAD
   res.ncpu = this->ncpu - rhs.ncpu;
#endif
   return res;
}

void xDeltaTimeData::operator-=(const xDeltaTimeData &rhs)
{
   this->sys -= rhs.sys;
   this->elapse -= rhs.elapse;
   this->elapsem -= rhs.elapsem;
   this->cpu -= rhs.cpu;
#ifdef USE_FINE_THREAD
   this->ncpu -= rhs.ncpu;
#endif
}

void xDeltaTimeData::operator/=(const xDeltaTimeData &rhs)
{
   this->sys = rhs.sys - this->sys;
   this->elapse = rhs.elapse - this->elapse;
   this->elapsem = rhs.elapsem - this->elapsem;
   this->cpu = rhs.cpu - this->cpu;
#ifdef USE_FINE_THREAD
   this->ncpu = rhs.ncpu - this->ncpu;
#endif
}

xDeltaTime::xDeltaTime(MPI_Comm world_) : world(world_), n(0), nb_proc(1), proc_id(0), scale(1.0e-6), scalen(1.0e-9)
{
   MPI_Comm_size(world, &nb_proc);
   MPI_Comm_rank(world, &proc_id);

   // initializing periode for time calculation
   periode = (double)1. / ((double)sysconf(_SC_CLK_TCK));

   // init zero
   zero.cpu = 0;
#ifdef USE_FINE_THREAD
   zero.ncpu = 0;
#endif
   zero.sys = 0;
   zero.elapse = 0;
   zero.elapsem = 0;

   // set initial reference
   set(init);
}

void xDeltaTime::get(double *cpu, double *sys, double *elapse)
{
   set();
#ifdef USE_FINE_THREAD
   *cpu = static_cast<double>(t_cur.cpu) + scalen * static_cast<double>(t_cur.ncpu);
#else
   *cpu = static_cast<double>(t_cur.cpu);
#endif
   *sys = static_cast<double>(t_cur.sys * periode);
   *elapse = static_cast<double>(t_cur.elapse) + scale * static_cast<double>(t_cur.elapsem);
}
double xDeltaTime::getElapseId(int id)
{
   return static_cast<double>(dt[id].elapse) + scale * static_cast<double>(dt[id].elapsem);
}
void xDeltaTime::set()
{
   struct tms buff;
   timeval current;
#ifdef USE_FINE_THREAD
   struct timespec tp;
   if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tp) < 0) return;
#endif
   if (times(&buff) < 0) return;
   gettimeofday(&current, nullptr);
#ifdef USE_FINE_THREAD
   t_cur.cpu = tp.tv_sec;
   t_cur.ncpu = tp.tv_nsec;
#else
   t_cur.cpu = buff.tms_utime;
#endif
   t_cur.sys = buff.tms_stime;
   t_cur.elapse = current.tv_sec;
   t_cur.elapsem = current.tv_usec;
}
void xDeltaTime::set(xDeltaTimeData &t)
{
   struct tms buff;
   struct timezone tz;
   timeval current;
#ifdef USE_FINE_THREAD
   struct timespec tp;
   if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tp) < 0) return;
#endif
   if (times(&buff) < 0) return;
   gettimeofday(&current, &tz);
#ifdef USE_FINE_THREAD
   t.cpu = tp.tv_sec;
   t.ncpu = tp.tv_nsec;
#else
   t.cpu = buff.tms_utime;
#endif
   t.sys = buff.tms_stime;
   t.elapse = current.tv_sec;
   t.elapsem = current.tv_usec;
}
int xDeltaTime::start(std::string stage, bool local)
{
   set();
   int i = dt.size();
   if (local)
      strtoind[stage] = -i;
   else
      strtoind[stage] = i;
   dt.push_back(t_cur);
   t.push_back(zero);
   return i;
}
void xDeltaTime::end(int id)
{
   set();
   dt[id] /= t_cur;
}
void xDeltaTime::end(std::string stage)
{
   set();
   int id = strtoind[stage];
   dt[(id > 0) ? id : -id] /= t_cur;
}

int xDeltaTime::initAccu(std::string stage, bool local)
{
   int i = dt.size();
   if (local)
      strtoind[stage] = -i;
   else
      strtoind[stage] = i;
   dt.push_back(zero);
   t.push_back(zero);
   return i;
}

void xDeltaTime::print(double tcpu, double tsys, double telapse, double min[], double max[], double mean[])
{
   bool do_more_print = false;
   double mi[3], mx[3], me[3];
   int j;
   std::map<std::string, int>::iterator it = strtoind.begin();
   std::map<std::string, int>::iterator itend = strtoind.end();
   if (it != itend) do_more_print = true;
   if (do_more_print)
   {
      bool got_local = false;
      tcpu = (tcpu == 0.0) ? 100.0 : tcpu;
      tsys = (tsys == 0.0) ? 100.0 : tsys;
      telapse = (telapse == 0.0) ? 100.0 : telapse;
      if (proc_id)
      {
         cout << setfill('=') << setw(137) << " " << endl;
         cout << setfill(' ') << "| " << setw(27) << "elapse time"
              << " | " << setw(27) << "system time"
              << " | " << setw(27) << "cpu time"
              << " | " << setw(42) << "for"
              << " |" << endl;
         cout << setfill('=') << setw(137) << " " << endl;
         cout << setfill(' ') << "| " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(44) << " |" << endl;
         for (; it != itend; ++it)
         {
            j = it->second;
            if (j > -1)
            {
               double cp, sy, el;
               computeDeltaToDouble(j, cp, sy, el);
               reduce(cp, sy, el, mi, mx, me);
               cout << "| " << fixed << setprecision(5) << setw(12) << el << " | " << setprecision(1) << setw(12)
                    << 100 * el / telapse;
               cout << " | " << setprecision(5) << setw(12) << sy << " | " << setprecision(1) << setw(12) << 100 * sy / tsys;
               cout << " | " << setprecision(5) << setw(12) << cp << " | " << setprecision(1) << setw(12) << 100 * cp / tcpu
                    << " | " << setw(42) << it->first << " |" << endl;
            }
            else
               got_local = true;
         }
         cout << setfill('=') << setw(137) << " " << setfill(' ') << endl;
      }
      else
      {
         double tcpu_mean = (mean[0] == 0.0) ? 100.0 : mean[0];
         double tsys_mean = (mean[1] == 0.0) ? 100.0 : mean[1];
         double tela_mean = (mean[2] == 0.0) ? 100.0 : mean[2];
         cout << setfill('=') << setw(137) << " " << endl;
         cout << setfill(' ') << "| " << setw(87) << "cpu time"
              << " | " << setw(42) << "for"
              << " |" << endl;
         cout << setfill('=') << setw(137) << " " << endl;
         cout << setfill(' ') << "| " << setw(12) << "val min"
              << " | " << setw(12) << "val max"
              << " | " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(12) << "val avg"
              << " | " << setw(12) << "% total avg"
              << " | " << setw(44) << " |" << endl;
         ostringstream ossys;
         ostringstream osela;
         for (; it != itend; ++it)
         {
            j = it->second;
            if (j > -1)
            {
               double cp, sy, el;
               computeDeltaToDouble(j, cp, sy, el);
               reduce(cp, sy, el, mi, mx, me);
               cout << "| " << fixed << setprecision(5) << setw(12) << mi[0] << " | " << setw(12) << mx[0] << " | " << setw(12)
                    << cp << " | " << setprecision(1) << setw(12) << 100 * cp / tcpu << " | " << setprecision(5) << setw(12)
                    << me[0] << " | " << setprecision(1) << setw(12) << 100 * me[0] / tcpu_mean;
               cout << " | " << setw(42) << it->first << " |" << endl;

               ossys << "| " << fixed << setprecision(5) << setw(12) << mi[1] << " | " << setw(12) << mx[1] << " | " << setw(12)
                     << sy << " | " << setprecision(1) << setw(12) << 100 * sy / tsys << " | " << setprecision(5) << setw(12)
                     << me[1] << " | " << setprecision(1) << setw(12) << 100 * me[1] / tsys_mean;
               ossys << " | " << setw(42) << it->first << " |" << endl;

               osela << "| " << fixed << setprecision(5) << setw(12) << mi[2] << " | " << setw(12) << mx[2] << " | " << setw(12)
                     << el << " | " << setprecision(1) << setw(12) << 100 * el / telapse << " | " << setprecision(5) << setw(12)
                     << me[2] << " | " << setprecision(1) << setw(12) << 100 * me[2] / tela_mean;
               osela << " | " << setw(42) << it->first << " |" << endl;
            }
            else
               got_local = true;
         }
         cout << setfill('=') << setw(137) << " " << setfill(' ') << endl;
         cout << setfill(' ') << "| " << setw(87) << "system time"
              << " | " << setw(42) << "for"
              << " |" << endl;
         cout << setfill('=') << setw(137) << " " << endl;
         cout << setfill(' ') << "| " << setw(12) << "val min"
              << " | " << setw(12) << "val max"
              << " | " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(12) << "val avg"
              << " | " << setw(12) << "% total avg"
              << " | " << setw(44) << " |" << endl;
         cout << ossys.str();
         cout << setfill('=') << setw(137) << " " << setfill(' ') << endl;
         cout << setfill(' ') << "| " << setw(87) << "elapse time"
              << " | " << setw(42) << "for"
              << " |" << endl;
         cout << setfill('=') << setw(137) << " " << endl;
         cout << setfill(' ') << "| " << setw(12) << "val min"
              << " | " << setw(12) << "val max"
              << " | " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(12) << "val avg"
              << " | " << setw(12) << "% total avg"
              << " | " << setw(44) << " |" << endl;
         cout << osela.str();
         cout << setfill('=') << setw(137) << " " << setfill(' ') << endl;
      }
      if (got_local)
      {
         cout << setfill(' ') << "| " << setw(27) << "local"
              << " |" << endl;
         cout << setfill('=') << setw(137) << " " << endl;
         cout << setfill(' ') << "| " << setw(27) << "elapse time"
              << " | " << setw(27) << "system time"
              << " | " << setw(27) << "cpu time"
              << " | " << setw(42) << "for"
              << " |" << endl;
         cout << setfill('=') << setw(137) << " " << endl;
         cout << setfill(' ') << "| " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(12) << "val"
              << " | " << setw(12) << "% total"
              << " | " << setw(44) << " |" << endl;
         for (it = strtoind.begin(); it != itend; ++it)
         {
            j = it->second;
            if (j < 0)
            {
               j = -j;
               double cp, sy, el;
               computeDeltaToDouble(j, cp, sy, el);
               cout << "| " << fixed << setprecision(5) << setw(12) << el << " | " << setprecision(1) << setw(12)
                    << 100 * el / telapse;
               cout << " | " << setprecision(5) << setw(12) << sy << " | " << setprecision(1) << setw(12) << 100 * sy / tsys;
               cout << " | " << setprecision(5) << setw(12) << cp << " | " << setprecision(1) << setw(12) << 100 * cp / tcpu
                    << " | " << setw(42) << it->first << " |" << endl;
            }
         }
         cout << setfill('=') << setw(137) << " " << setfill(' ') << endl;
      }
   }
}

void xDeltaTime::print(std::string stage)
{
   int j = strtoind[stage];
#ifdef USE_FINE_THREAD
   double tcpu = static_cast<double>(dt[j].cpu) + scalen * static_cast<double>(dt[j].ncpu);
#else
   double tcpu = static_cast<double>(dt[j].cpu * periode);
#endif
   double tsys = static_cast<double>(dt[j].sys * periode);
   double telapse = static_cast<double>(dt[j].elapse) + scale * static_cast<double>(dt[j].elapsem);
   double min[3], max[3], mean[3];
   reduce(tcpu, tsys, telapse, min, max, mean);
   cout << " " << endl;
   cout << "==xDeltaTime output ================================================================================================"
        << fixed << setprecision(5) << endl;
   cout << "CPU    times since creation : " << tcpu;
   if (!proc_id) cout << " min/max/average : " << min[0] << "/" << max[0] << "/" << mean[0];
   cout << endl << "SYSTEM times since creation : " << tsys;
   if (!proc_id) cout << " min/max/average : " << min[1] << "/" << max[1] << "/" << mean[1];
   cout << endl << "ELAPSE times since creation : " << telapse;
   if (!proc_id) cout << " min/max/average : " << min[2] << "/" << max[2] << "/" << mean[2];
   cout << endl
        << "===================================================================================================================="
        << endl;

   print(tcpu, tsys, telapse, min, max, mean);

   cout << "==End xDeltaTime output ============================================================================================"
        << endl;
}

void xDeltaTime::print()
{
   set();
   t_cur -= init;
#ifdef USE_FINE_THREAD
   double tcpu = static_cast<double>(t_cur.cpu) + scalen * static_cast<double>(t_cur.ncpu);
#else
   double tcpu = static_cast<double>(t_cur.cpu * periode);
#endif
   double tsys = static_cast<double>(t_cur.sys * periode);
   double telapse = static_cast<double>(t_cur.elapse) + scale * static_cast<double>(t_cur.elapsem);
   double min[3], max[3], mean[3];
   reduce(tcpu, tsys, telapse, min, max, mean);
   cout << " " << endl;
   cout << "==xDeltaTime output ================================================================================================"
        << fixed << setprecision(5) << endl;
   cout << "CPU    times since creation : " << tcpu;
   if (!proc_id) cout << " min/max/average : " << min[0] << "/" << max[0] << "/" << mean[0];
   cout << endl << "SYSTEM times since creation : " << tsys;
   if (!proc_id) cout << " min/max/average : " << min[1] << "/" << max[1] << "/" << mean[1];
   cout << endl << "ELAPSE times since creation : " << telapse;
   if (!proc_id) cout << " min/max/average : " << min[2] << "/" << max[2] << "/" << mean[2];
   cout << endl
        << "===================================================================================================================="
        << endl;

   print(tcpu, tsys, telapse, min, max, mean);

   cout << "==End xDeltaTime output ============================================================================================"
        << endl;
}

void xDeltaTime::reduce(double &tcpu, double &tsys, double &telapse, double min[], double max[], double mean[]) const
{
   double data[3];
   data[0] = tcpu;
   data[1] = tsys;
   data[2] = telapse;
   // mean nead to be initialized since MPI_SUM add to the current value.
   mean[0] = 0;
   mean[1] = 0;
   mean[2] = 0;
   MPI_Reduce(data, max, 3, MPI_DOUBLE, MPI_MAX, 0, world);
   MPI_Reduce(data, min, 3, MPI_DOUBLE, MPI_MIN, 0, world);
   MPI_Reduce(data, mean, 3, MPI_DOUBLE, MPI_SUM, 0, world);
   mean[0] /= nb_proc;
   mean[1] /= nb_proc;
   mean[2] /= nb_proc;
   return;
}

void xDeltaTime::get(std::string stage, double val[], double mi[], double mx[], double me[]) const
{
   const int &j = strtoind.at(stage);
   int i = (j > -1) ? j : -j;
   computeDeltaToDouble(i, val[0], val[1], val[2]);
   if (j > -1)
      reduce(val[0], val[1], val[2], mi, mx, me);
   else
   {
      for (int l = 0; l < 3; ++l) mi[l] = me[l] = mx[l] = val[l];
   }
}
void xDeltaTime::computeDeltaToDouble(int i, double &cp, double &sy, double &el) const
{
   assert(i > -1);
   el = static_cast<double>(dt[i].elapse) + scale * static_cast<double>(dt[i].elapsem);
   sy = static_cast<double>(dt[i].sys * periode);
#ifdef USE_FINE_THREAD
   cp = static_cast<double>(dt[i].cpu) + scalen * static_cast<double>(dt[i].ncpu);
#else
   cp = static_cast<double>(dt[i].cpu * periode);
#endif
}

}  // namespace xtool
