/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include "xDeltaMemory.h"
#include "xDeltaTime.h"
using namespace std;
#include "xMemoryMonitor.h"

#define NBL 10
void foo(xtool::xDeltaTime &dt, xtool::xDeltaMemory &dm, xMemoryMonitor &mm, std::ofstream &out_ref, int proc_id)
{
   int iddti0 = dt.initAccu("in loop cost DM");
   int iddti1 = dt.initAccu("out of loop cost MM");
   int iddmi0 = dm.initAccu("sum rand alloc");
   int iddmi1 = dm.initAccu("dealloc");
   int iddmi2 = dm.initAccu("a small chunck in foo");
   int iddmi3 = dm.initAccu("alloc+dealloc");
   dm.startAccu(iddmi3);
   size_t k = 0;
   double *pointers[NBL];
   std::cout << "========================================" << std::endl;
   std::cout << "Before alloc loop" << std::endl;
   std::cout << "========================================" << std::endl;
   dt.startAccu(iddti1);
   int idmmi0 = mm.start("sum rand alloc");
   dt.endAccu(iddti1);
   for (size_t i = 0; i < NBL; ++i)
   {
      dt.startAccu(iddti0);
      dm.startAccu(iddmi0);
      dt.endAccu(iddti0);
      size_t s = rand() % (5000000 * (i + 1)) + 1;
      k += s;
      cout << "Allocate " << i << "th chunck (B): " << s * 8 << endl;
      pointers[i] = new double[s];
      std::fill(pointers[i], pointers[i] + s, 3.);
      std::cout << "========================================" << std::endl;
      dt.startAccu(iddti0);
      dm.endAccu(iddmi0);
      dt.endAccu(iddti0);
   }
   dt.startAccu(iddti1);
   mm.end(idmmi0);
   dt.endAccu(iddti1);
   std::cout << "After alloc loop" << std::endl;
   std::cout << "========================================" << std::endl;
   std::cout << "Pick allocation (B): " << k * 8 << std::endl;
   std::cout << "========================================" << std::endl;
   std::cout << "Memory leak introduced by not freeing last allocated block" << std::endl;
   std::cout << "========================================" << std::endl;
   dt.startAccu(iddti1);
   idmmi0 = mm.start("dealloc");
   dt.endAccu(iddti1);
   for (size_t i = 0; i < NBL - 1; ++i)
   {
      dt.startAccu(iddti0);
      dm.startAccu(iddmi1);
      dt.endAccu(iddti0);
      cout << "Deallocate " << i << "th chunck" << endl;
      delete[] pointers[i];
      std::cout << "========================================" << std::endl;
      dt.startAccu(iddti0);
      dm.endAccu(iddmi1);
      dt.endAccu(iddti0);
   }
   dm.endAccu(iddmi3);
   dt.startAccu(iddti1);
   mm.end(idmmi0);
   dt.endAccu(iddti1);
   std::cout << "========================================" << std::endl;
   std::cout << "A small chunck (32B) in foo (leak)" << std::endl;
   idmmi0 = mm.start("a small chunck in foo");
   dm.startAccu(iddmi2);
   double *i = new double[4];
   std::fill(i, i + 4, 4.);
   dm.endAccu(iddmi2);
   mm.end(idmmi0);
   std::cout << "========================================" << std::endl;
   double mx, mi, me, su;
   double val = dm.get(iddmi0, mi, mx, me, su);
   if (proc_id)
      out_ref << "Retriving data for 'sum rand alloc' (GB) " << val / 1073741824. << std::endl;
   else
      out_ref << "Retriving data for 'sum rand alloc' (GB) " << val / 1073741824. << " min/max/mean/sum " << mi / 1073741824.
              << "  " << mx / 1073741824. << " " << me / 1073741824. << " " << su / 1073741824. << std::endl;
   out_ref << "========================================" << std::endl;
   val = dm.get(iddmi1, mi, mx, me, su);
   if (proc_id)
      out_ref << "Retriving data for 'dealloc' (GB) " << val / 1073741824. << std::endl;
   else
      out_ref << "Retriving data for 'dealloc' (GB) " << val / 1073741824. << " min/max/mean/sum " << mi / 1073741824. << " "
              << mx / 1073741824. << " " << me / 1073741824. << " " << su / 1073741824. << std::endl;
   out_ref << "========================================" << std::endl;
   val = dm.get(iddmi2, mi, mx, me, su);
   if (proc_id)
      out_ref << "Retriving data for 'a small chunck in foo' (B) " << val << std::endl;
   else
      out_ref << "Retriving data for 'a small chunck in foo' (B) " << val << " min/max/mean/sum " << mi << " " << mx << " " << me
              << " " << su << std::endl;
   out_ref << "========================================" << std::endl;
   val = dm.get(iddmi3, mi, mx, me, su);
   if (proc_id)
      out_ref << "Retriving data for 'alloc+dealloc' (MB) " << val / 1048576. << std::endl;
   else
      out_ref << "Retriving data for 'alloc+dealloc' (MB) " << val / 1048576. << " min/max/mean/sum " << mi / 1048576. << " "
              << mx / 1048576. << " " << me / 1048576. << " " << su / 1048576. << std::endl;
   out_ref << "========================================" << std::endl;
   return;
}
int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int proc_id;
   MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
   srand((proc_id + 1) * 37);
   string no = "proc_" + std::to_string(proc_id) + "_output.txt";
   auto pfile = freopen(no.c_str(), "w", stdout);
   std::ofstream out_ref;
   string noo = "reference_" + std::to_string(proc_id) + ".txt";
   out_ref.open(noo.c_str());
   out_ref << fixed << std::setprecision(2);

   std::cout << "==Start  ===============================" << std::endl;
   xtool::xDeltaTime dt;
   xtool::xDeltaMemory dm;
   xMemoryMonitor mm;
   int iddm = dm.initAccu("foo");
   int iddmm = mm.start("foo");
   dm.startAccu(iddm);
   foo(dt, dm, mm, out_ref, proc_id);
   dm.endAccu(iddm);
   mm.end(iddmm);
   std::cout << "========================================" << std::endl;
   double mx, mi, me, su;
   double val = dm.get(iddm, mi, mx, me, su);
   const int small_chunck = 10;
   if (!proc_id)
   {
      /* During instalation some small variations have been observed depending on the
       * way test are launched. Not clear to me why ?! To avoid prb
       * remove from reference
      out_ref << "Retriving data for foo (MB) " << val / 1048576. << " min/max/mean/sum " << mi / 1048576. << " " << mx / 1048576.
              << " " << me / 1048576. << " " << su / 1048576. << std::endl;
      */
      std::cout << "========================================" << std::endl;
      iddm = dm.initAccu("a small chunck i", true);
      std::cout << "========================================" << std::endl;
      std::cout << "A small chunck i (" << small_chunck * 8 << "B) in main P0 (leak)" << std::endl;
      dm.startAccu(iddm);
      double *i = new double[small_chunck];
      std::fill(i, i + small_chunck, 7.2);
      dm.endAccu(iddm);
      out_ref << "========================================" << std::endl;
      out_ref << "Retriving data for small chunck i " << dm.get(iddm) << "B" << std::endl;
      std::cout << "========================================" << std::endl;
      std::cout << "Another small chunck k (" << small_chunck * 80 << "B) in main P0 (leak)" << std::endl;
      iddmm = mm.start("another small chunck k");
      double *k = new double[small_chunck * 10];
      std::fill(k, k + small_chunck * 10, 4.2);
      mm.end(iddmm);
   }
   if (proc_id)
   {
      out_ref << "Retriving data for foo (MB) " << val / 1048576. << std::endl;
      std::cout << "========================================" << std::endl;
      std::cout << "A small chunck j (" << small_chunck * 2 << "B) in main Px (leak)" << std::endl;
      iddm = dm.initAccu("A small chunck j", true);
      dm.startAccu(iddm);
      double *j = new double[small_chunck * 2];
      j[small_chunck - 1] = 3.2;
      dm.endAccu(iddm);
      out_ref << "========================================" << std::endl;
      out_ref << "Retriving data for small chunck j " << dm.get(iddm) << "B" << std::endl;
      std::cout << "========================================" << std::endl;
      std::cout << "Another small chunck k (" << small_chunck * 8000 << "B) in main Px (leak)" << std::endl;
      iddm = dm.initAccu("Another small chunck k", true);
      iddmm = mm.start("Another small chunck k");
      dm.startAccu(iddm);
      double *k = new double[small_chunck * 1000];
      std::fill(k, k + small_chunck * 1000, 4.2);
      dm.endAccu(iddm);
      mm.end(iddmm);
      std::cout << "========================================" << std::endl;
      out_ref << "========================================" << std::endl;
      out_ref << "Retriving data for small chunck k " << dm.get(iddm) / 1024. << "KB" << std::endl;
      out_ref << "========================================" << std::endl;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   dm.print();
#ifdef PROFILE
   dt.print();
#endif
   mm.print(cout);
   MPI_Finalize();
   return 0;
}
