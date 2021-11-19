/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _NODEANDCONNECTEDEDGEDIST_H
#error "This header must not be included directly"
#endif
// xTool
#include "xExportStringDist.h"

namespace xgraph
{
namespace internal
{
template <typename NODE, class DATAMANAGER, typename ID, typename L, typename W>
void exportDot(const DATAMANAGER &graph, std::string f_name, MPI_Comm world, const ID &fid, L *fl, W *fw)
{
   // local
   int proc_id, nb_proc;
   std::ostringstream os;

   MPI_Comm_rank(world, &proc_id);
   MPI_Comm_size(world, &nb_proc);
   auto node_range = xtool::make_range(graph.beginKey(), graph.endKey());
   size_t nb_node = node_range.size();

   if (!proc_id) os << "digraph X {" << std::endl;
   if (nb_proc > 1 && nb_node)
   {
      os << "  subgraph cluster_" << proc_id << " {" << std::endl;
      os << "     label=\"p" << proc_id << "\";" << std::endl;
   }

   if (fl)
   {
      for (auto &key : node_range)
      {
         unsigned int parent = fid(key);
         os << parent << "[label=\"" << (*fl)(key) << "\"];" << std::endl;
      }
   }
   else
   {
      for (auto &key : node_range)
      {
         unsigned int parent = fid(key);
         os << parent << ";" << std::endl;
      }
   }
   if (nb_proc > 1 && nb_node) os << "  }" << std::endl;
   if (fw)
   {
      for (auto &key : node_range)
      {
         const NODE *node = graph.getData(*key);
         assert(node);
         unsigned int parent = fid(key);
         for (auto pair : node->rangeChild())
         {
            os << parent << " -> " << fid(pair.first) << " [label=" << (*fw)(pair.second) << "];" << std::endl;
         }
      }
   }
   else
   {
      for (auto &key : node_range)
      {
         const NODE *node = graph.getData(*key);
         assert(node);
         unsigned int parent = fid(key);
         for (auto pair : node->rangeChild())
         {
            os << parent << " -> " << fid(pair.first) << ";" << std::endl;
         }
      }
   }
   if (proc_id == nb_proc - 1) os << "}" << std::endl;

   // export string to file
   xtool::xExportStringDist(f_name.c_str(), os, world);
}
}  // end namespace
}  // end namespace
