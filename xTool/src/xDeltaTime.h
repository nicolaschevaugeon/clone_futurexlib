/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef XDELTATIME_H
#define XDELTATIME_H

// ----------------------------------------------------------------------------
// HEADERS
// ----------------------------------------------------------------------------
#include <sys/time.h>
#include <sys/times.h>

#include <map>
#include <string>
#include <vector>

#include "mpi.h"

//#define USE_FINE_THREAD 1

#ifdef USE_FINE_THREAD
extern "C"
{
#include <time.h>
}
#endif

namespace xtool
{
// ----------------------------------------------------------------------------
// xDeltaTime
// ----------------------------------------------------------------------------
/*! \class xDeltaTime
    \ingroup Utilities
    \brief Process timing measurment.

<br/>    xDeltaTime is a small class that give the ability to instrument in a lightweight maner
         any part of a code.
<br/>    It's not a magic implementation.
<br/>    Instrumentation based on this class will have a intrinsic cost which is not negligible when
         used in small fine grain functions.

<br/>    To catch more fine grain functions and to deal with threaded application the macro USE_FINE_THREAD, if set,
         give better result with the burden of linking your application with -lrt. Ok on linux but not tested else where
         As this is optional it doesn't integrate optimally the current implementation : for one measure 3 functions are
         called !

<br/>    principal :
<br/>    xDeltaTime is storing in a container elapse, cpu and system time measured from 2 points.
<br/>    Indexation of this container is done with a string chosen by the user.
<br/>    * The starting and ending point may be for a unique usage in the life of the object instance of this class.
<br/>    In this case method "start" is used to initiate the first points of timing measurement and to index, with
         the given string, this measure in the container. This method return the index of this measure as a int.
         The second point is then set with "end" method with correct index of the measure given. Its the responsibility
         of the user to check consistence of starting an ending points.
<br/>    * The starting and ending point may be for a multiple usage in the life of the object instance of this class.
         Typically in a loop some one may want the contribution of part of this loop to the overall execution of it.
<br/>    In this case method "initAccu" is used to initiate and return the index of the measure in the container corresponding
         to the string given. This operation is done only once before the loop to instrument. In the loop itself "startAccu"
         and "endAccu" are used with index id to define starting and ending point of the measure.
<br/>    * Starting method have an extra boolean parameter call local defaulted to false. When in parallel context one may
         want to monitor a task only in some process. In this case when calling "start" or "initAccu" you must turn local to
         true. It avoid when reducing information from all process on zero one to try to collect information not available on
         all process. Those information are just outputted locally where they are present.
<br/>    Data and statistic are given by a call to "print". This method print in stdout (cout) all measure created since the
         creation of the  object instance of this class. Is given the elapse time, the system time and the cpu time. A percentage
         for all these data is given by using some reference measure. If called with no argument the reference is the time from
         the creation of the object instance to the call of the print. If called with  a string (corresponding to a particular
         measure) the reference is that particular measure.

<br/>    If user want to managed is measure himself the "get" method give him directly a timing point without interfering with the
         container mechanism.

<br/>    nota :  "startAccu","endAccu","end" methode existe with string argument. They should not be used as they intrinsically
               cost more then their equivalent with index argument. For small part of code instrumented their use may become
               preponderant in term of cpu time !!!! They are here only if user really didn't find any solution to bring index
               id at the level of theire call.

<br/>    Example :

<pre>
void a_methode_somewhere
{

  // creation of the object
  xDeltaTime dte;

  // somme measure instrumented in a loop of a_methode_somewhere
  int iddt1=dte.initAccu("ApplyNeuman");
  int iddt2=dte.initAccu("Residual calc");
  int iddt3=dte.initAccu("Assemble Mat");
  int iddt4=dte.initAccu("Solve");


  int iddt0=dte.start("some stuff");
  // some stuff done only once in a_methode_somewhere
  xIntegrationRule toto(2);
  ....
  ....
  ....
  ....
  // End of some stuff
  dte.end(iddt0);

  // Befor a big loop with many things in it
  iddt0=dte.start("The loop");
  while ( some_condition )
  {
      xlinalg::xCSRVector solnl(numdofs);
      xlinalg::xCSRVector bnl(numdofs);
      ......
      ......
      ......

     // --------------------- neumann - Fext
     dte.startAccu(iddt1);
     applyNeumannBC(dis................
     dte.endAccu(iddt1);
      ......
      ......
      ......
      ......
     // --------------------- Residual calc
     dte.startAccu(iddt2);
     assembl.....................
      ......
      ......
      ......
     dte.endAccu(iddt2);
    // --------------------- Matrix Assemble
    dte.startAccu(iddt3);
    xFormBilinearWithLaw..............
      ......
      ......
      ......
    dte.endAccu(iddt3);
    // ---------- SOLVING ---------------
    xlinalg::xLinearSystemSolverTaucs solver;
    dte.startAccu(iddt4);
    solver.solve(Anl, bnl, solnl);
    dte.endAccu(iddt4);
      ......
      ......
      ......
  }
  // after end of the loop
  dte.end(iddt0);

  // print times statistic since creation of dte = time reference is the times
  // passed in a_methode_somewhere
  dte.print();

  // print times statistic with measure of "The Loop" as reference
  // => give percentage of ApplyNeuman for example against "The Loop" and not "a_methode_somewhere"
  dte.print("The Loop");

  //end of a_methode_somewhere
}

</pre>
*/
struct xDeltaTimeData
{
  public:
#ifdef USE_FINE_THREAD
   time_t cpu;
   long ncpu;
#else
   clock_t cpu;
#endif
   clock_t sys;
   time_t elapse;
   suseconds_t elapsem;
   void operator=(const xDeltaTimeData &rhs);
   xDeltaTimeData operator+(const xDeltaTimeData &rhs);
   xDeltaTimeData operator-(const xDeltaTimeData &rhs);
   void operator-=(const xDeltaTimeData &rhs);
   void operator/=(const xDeltaTimeData &rhs);
};

class xDeltaTime
{
  public:
   xDeltaTime(MPI_Comm world_ = MPI_COMM_WORLD);
   ~xDeltaTime() = default;

   int initAccu(std::string stage, bool local = false);
   inline void startAccu(int id) { set(t[id]); }
   inline void endAccu(int id)
   {
      set();
      cum(dt[id], t_cur, t[id]);
   }

   int start(std::string stage, bool local = false);
   void end(int id);

   void get(double *cpu, double *sys, double *elapse);
   void get(std::string stage, double val[], double mi[], double mx[], double me[]) const;
   double getElapseId(int id);

   void print();
   void print(std::string stage);

   // unefficient methode used in desaspare case (a int can't comme from initAccu)
   inline void startAccu(std::string stage)
   {
      int i = strtoind[stage];
      set(t[(i > 0) ? i : -i]);
   }
   inline void endAccu(std::string stage)
   {
      set();
      int j = strtoind[stage];
      j = (j > 0) ? j : -j;
      cum(dt[j], t_cur, t[j]);
   }
   void end(std::string stage);

  private:
   MPI_Comm world;
   // inverse of clock ticks per second to transforme clock_t value of times in second
   int n, nb_proc, proc_id;
   double periode;
   const double scale;
   const double scalen;
   xDeltaTimeData init;
   xDeltaTimeData zero;
   xDeltaTimeData t_cur;
   std::map<std::string, int> strtoind;
   std::vector<xDeltaTimeData> dt;
   std::vector<xDeltaTimeData> t;
   void set();
   void set(xDeltaTimeData &it);
   void print(double tcpu, double tsys, double telapse);
   void print(double tcpu, double tsys, double telapse, double min[], double max[], double mean[]);
   void reduce(double &tcpu, double &tsys, double &telapse, double min[], double max[], double mean[]) const;
   inline void cum(xDeltaTimeData &c, xDeltaTimeData &u, xDeltaTimeData &r)
   {
#ifdef USE_FINE_THREAD
      c.ncpu += (u.ncpu - r.ncpu);
      r.ncpu = u.ncpu;
#endif
      c.cpu += (u.cpu - r.cpu);
      c.sys += (u.sys - r.sys);
      c.elapse += (u.elapse - r.elapse);
      c.elapsem += (u.elapsem - r.elapsem);
      r.cpu = u.cpu;
      r.sys = u.sys;
      r.elapse = u.elapse;
      r.elapsem = u.elapsem;
   }
   void computeDeltaToDouble(int i, double &cp, double &sy, double &el) const;
};

/*! \class xNoDeltaTime
    \ingroup Utilities
    \brief No timing measurement class but providing same interface as xDeltaTime.

<br/> xNoDeltaTime give a way to simplify monitoring by writing DeltaTime sequence definitively
      in the code and use xDeltaTime at creation time if monitoring is asked and xNoDeltaTime if no
      monitoring is asked. xNoDeltaTime provides a fake object that cost nothing if well used (if
      encapsulated in a std::function it start costing !!!). The compiler optimize out all the calls.
      This fake object having the same API as xDeltaTime it permit to leave monitoring instructions
      (Start,....) in your code without use of macro to switch from monitoring or not (except eventually
      at creation time).

<br/> If monitoring object are passed as function argument do use a template argument to let this object be
      a xDeltaTime or a xNoDeltaTime object. DON'T use a std::function to encapsulate the object because if
      your function is called with a xNoDeltaTime instance it will cost some time as compiler will not have the
      opportunity to optimize out the monitoring calls and std::functions has an intrinsic cost when used.

<br/> Example :
With the following function
<pre>
size_t foo(xtool::xNoDeltaTime &ndt,xtool::xDeltaTime &dt,const char *msg)
{
   int iddti0 = dt.initAccu(msg);
   int iddti1 = ndt.initAccu(msg);
   size_t j=0,k=0;
   for (size_t i = 0; i < 1000000; ++i)
   {
      ndt.startAccu(iddti1);
      dt.startAccu(iddti0);
      k=i-j;j+=k;// some-thing;
      dt.endAccu(iddti0);
      ndt.endAccu(iddti1);
   }
   return k+j;
}
</pre>
Callgrind shows that no instruction are used with calls  ndt.startAccu(iddti1); ndt.endAccu(iddti1); and 66% of the function
cost is related to dt.startAccu(iddti0);  dt.endAccu(iddti0); calls. This is natural as "some-thing" do not cost at all and
the only real work done in this function is monitoring (call to times).

 */
class xNoDeltaTime
{
  public:
   xNoDeltaTime(MPI_Comm world_ = MPI_COMM_WORLD) {}
   ~xNoDeltaTime() = default;

   int initAccu(std::string stage, bool local = false) { return 0; }
   inline void startAccu(int id) {}
   inline void endAccu(int id) {}

   int start(std::string stage, bool local = false) { return 0; }
   void end(int id) {}

   void get(double *cpu, double *sys, double *elapse) {}
   void get(std::string stage, double val[], double mi[], double mx[], double me[]) const {}
   double getElapseId(int id) { return 0.; }

   void print() {}
   void print(std::string stage) {}

   // unefficient methode used in desaspare case (a int can't comme from initAccu)
   inline void startAccu(std::string stage) {}
   inline void endAccu(std::string stage) {}
   void end(std::string stage) {}
};

}  // namespace xtool
#endif
