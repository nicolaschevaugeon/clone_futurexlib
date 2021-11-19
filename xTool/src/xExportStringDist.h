/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef xEXPORTSTRINGDIST_H
#define xEXPORTSTRINGDIST_H
#include <sstream>
#include <string>
#include "mpi.h"
namespace xtool
{
/// A function to export a set of distributed string into a common file.
//! The only assertion is about string order in the file. It follows rank
//! order in the given communicator. String in proc 0 will be written first
//! then string of proc 1 written second, etc
void xExportStringDist(const std::string& f_name, const std::string& s, MPI_Comm world);
/// A function to export a set of distributed string generated in ostringstream into a common file.
//! It use the above string version
void xExportStringDist(const std::string& f_name, const std::ostringstream& os, MPI_Comm world);
}  // end namespace
#endif
