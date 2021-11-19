#include <fstream>

#include "mpi.h"
#include "xDeltaTime.h"
#include "xMPIEnv.h"

using namespace std;

// a function to test monitoring
template <typename DT>
size_t foo(DT &dt, const char *msg)
{
   int iddti0 = dt.initAccu(msg);
   size_t j = 0, k = 0;
   for (size_t i = 0; i < 1000000; ++i)
   {
      dt.startAccu(iddti0);
      k = i - j;
      j += k;  // some-thing;
      dt.endAccu(iddti0);
   }
   return k + j;
}

int main(int argc, char *argv[])
{
   xtool::xMPIEnv::init(argc, argv);

   // local variable
   MPI_Comm world = MPI_COMM_WORLD;
   int proc_id, iddt0;

   // get rank for world
   MPI_Comm_rank(world, &proc_id);

   // testing class
   xtool::xDeltaTime dt;
   xtool::xNoDeltaTime ndt;

   // test with true monitoring
   iddt0 = dt.start("foo dt");
   foo(dt, "foo dt oper");
   dt.end(iddt0);

   // test with false monitoring
   iddt0 = dt.start("foo ndt");
   foo(ndt, "foo ndt oper");
   dt.end(iddt0);

   // output stats
   dt.print();
   ndt.print();

   // open reference
   ofstream ofs("proc_" + std::to_string(proc_id) + "_ref.txt", ios::out | ios::trunc);

   // without monitoring time must be null
   double val[3];
   double mi[3];
   double mx[3];
   double me[3];
#ifndef NDEBUG
   double eps = 1.E-1;
#else
   double eps = 1.E-4;
#endif
   dt.get("foo ndt", val, mi, mx, me);
   ofs << "foo ndt should be null" << std::endl;
   ofs << "mi " << ((mi[0] < eps) ? 0. : mi[0]) << " " << ((mi[1] < eps) ? 0. : mi[1]) << " " << ((mi[2] < eps) ? 0. : mi[2])
       << std::endl;
   ofs << "mx " << ((mx[0] < eps) ? 0. : mx[0]) << " " << ((mx[1] < eps) ? 0. : mx[1]) << " " << ((mx[2] < eps) ? 0. : mx[2])
       << std::endl;
   ofs << "me " << ((me[0] < eps) ? 0. : me[0]) << " " << ((me[1] < eps) ? 0. : me[1]) << " " << ((me[2] < eps) ? 0. : me[2])
       << std::endl;

   // ending prog
   ofs.close();
   xtool::xMPIEnv::finalize();
   return 0;
}
