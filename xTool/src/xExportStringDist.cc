/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#include "xExportStringDist.h"
namespace xtool
{
void xExportStringDist(const std::string& f_name, const std::string& s, MPI_Comm world)
{
   // local
   MPI_File file;
   MPI_Status status;

   // open file an clean it
   MPI_File_open(world, f_name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
   MPI_File_set_size(file, 0);

   // nb_info store nb of char in string (char considerer here to be 1 Byte long)
   const MPI_Offset nb_info = s.size();

   // compute offset
   MPI_Offset offsets = 0;
   MPI_Scan(&nb_info, &offsets, 1, MPI_OFFSET, MPI_SUM, world);
   MPI_Offset dist;
   dist = offsets - nb_info;

   // write info
   if (nb_info) MPI_File_write_at(file, dist, s.c_str(), nb_info, MPI_CHAR, &status);

   // close file
   MPI_File_close(&file);
}
void xExportStringDist(const std::string& f_name, const std::ostringstream& os, MPI_Comm world)
{
   xExportStringDist(f_name, os.str(), world);
}
}  // end namespace
