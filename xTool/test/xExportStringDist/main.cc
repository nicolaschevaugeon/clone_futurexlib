/*
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
 */
#include <iostream>
#include <string>
#include "xExportStringDist.h"
#include "xMPIEnv.h"

using namespace xtool;
using namespace std;

int main(int argc, char *argv[])
{
   // initialize mpi universe
   xMPIEnv::init(argc, argv);

   // local variable
   MPI_Comm world = MPI_COMM_WORLD;
   MPI_Comm colworld;
   int proc_id;

   // get rank for world
   MPI_Comm_rank(world, &proc_id);

   // output (distributed)
   std::ostringstream os;
   os << "message from world proc " << proc_id << endl;
   string msg("message from world proc ");
   msg += to_string(proc_id) + "\n";

   // test on global comunicator
   xExportStringDist("world_os.txt", os, world);
   xExportStringDist("world_string.txt", msg, world);

   // test on splitted comunicator: One file per color (io optimisation for high number of proc)
   if (proc_id < 2)
   {
      MPI_Comm_split(world, 0, 0, &colworld);
      xExportStringDist("colworld1_os.txt", os, colworld);
   }
   else
   {
      MPI_Comm_split(world, 1, 0, &colworld);
      xExportStringDist("colworld2_os.txt", os, colworld);
   }

   return xMPIEnv::finalize();
}
