/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _NODEANDCONNECTEDEDGEDIST_H
#define _NODEANDCONNECTEDEDGEDIST_H
#include <numeric>
#include "nodeAndConnectedEdge.h"
#include "nodeAndConnectedEdgeDistInternal.h"
#include "xMPIDataType.h"
namespace xgraph
{
//==============================================================================================
/// A simple function that output a distributed nodeTo graph in dot format
/*!
 * NODE the type of nodes of the graph (nodeTo<KEYTYPE,ET,NT>)
 * DATAMANAGER the data manager type that hold NODE nodes associated to KEYTYPE keys.
 * ID a function that associate a unique ID across process to the key holding a node
 */
template <typename NODE, class DATAMANAGER, typename ID>
void exportDot(const DATAMANAGER &graph, std::string f_name, MPI_Comm world, const ID &fid)
{
   auto *fw_null = static_cast<unsigned int (*)(const double)>(nullptr);
   auto *fl_null = static_cast<char (*)(const typename NODE::key_t *)>(nullptr);
   internal::exportDot<NODE, DATAMANAGER, ID>(graph, f_name, world, fid, fl_null, fw_null);
}
/// A simple function that output a distributed nodeTo graph in dot format
/*!
 * NODE the type of nodes of the graph (nodeTo<KEYTYPE,ET,NT>)
 * DATAMANAGER the data manager type that hold NODE nodes associated to KEYTYPE keys.
 * ID a function that associate a unique ID across process to the key holding a node
 * L a function that associate a label(anything that iostream handle corectelly) to the node related
 * to the key given as argument
 */
template <typename NODE, class DATAMANAGER, typename ID, typename L>
void exportDot(const DATAMANAGER &graph, std::string f_name, MPI_Comm world, const ID &fid, const L &fl)
{
   auto *fw_null = static_cast<unsigned int (*)(const double)>(nullptr);
   internal::exportDot<NODE, DATAMANAGER, ID>(graph, f_name, world, fid, &fl, fw_null);
}
/// A simple function that output a distributed nodeTo graph in dot format
/*!
 * NODE the type of nodes of the graph (nodeTo<KEYTYPE,ET,NT>)
 * DATAMANAGER the data manager type that hold NODE nodes associated to KEYTYPE keys.
 * ID a function that associate a unique ID across process to the key holding a node
 * L a function that associate a label(anything that iostream handle corectelly) to the node related
 * to the key given as argument
 * W a function that associate a unsigned int to the weight of each edges.
 */
template <typename NODE, class DATAMANAGER, typename ID, typename L, typename W>
void exportDot(const DATAMANAGER &graph, std::string f_name, MPI_Comm world, const ID &fid, const L &fl, const W &fw)
{
   internal::exportDot<NODE, DATAMANAGER, ID>(graph, f_name, world, fid, &fl, &fw);
}

/// A breath first search algorithm dedicated to distributed graph traversal via edges following nodeTo API.
/*!
 * Graph distribution follow the principal that edges are not linking node across process. Node of the graph
 * may be present in more then one process. But locally all its child edges point to local nodes. A Boolean
 * given by hasRemoteChild method indicate that a node got a remote counterpart with child. A node and one of
 * its children being duplicate on same remote process  must have only one edge connecting them across those process.
 *
 * Tree traversal respect breath first search algorithm by passing in review nodes that are
 * topologically closer from source first. All children of a node are traversed via its outgoing
 * edges before passing to children of the children. Algorithm insure that edges are traversed only once.
 * The distributed nature of the graph is taken into account by using distance to source criteria. Locally
 * each proc stop its job when all local nodes at a specific distance are visited. Then communication occurs for
 * nodes having remote counterpart. Proc receiving information add child to the next generation locally. If new
 * generation is not empty in at least one process then all process do their job on this generation. And so on
 *
 * NODE the type of nodes of the graph (Following nodeTo API NODE must provide at least key_t,weight_t and child_t types)
 * DATAMANAGER the data manager type that hold NODE nodes associated to key_t keys.
 * EX A type representing an exchanger tuned by user to act properly when communication is involved. It should derive from
 * breadthFirstSearchDistExchanger and potentially integrate F work (i.e. F is not called for incoming edge traversal
 * in remote process arriving in this process).
 * F a function that is run for every edges traversed. Its argument are the 2 nodes it connect (as key_t).
 *   As the edge is oriented the source is the last argument and the target is in the first pair member of
 *   the first argument. The second pair member is the weight associated to the edge.
 */
template <typename NODE, template <typename> class DATAMANAGER, typename EX, typename F>
void breadthFirstSearchDist(const DATAMANAGER<NODE> &graph, const typename NODE::key_t *key_source, EX &exchanger, F func_node)
{
   // local type
   typedef typename NODE::key_t node_key_t;
   static_assert(std::is_same<node_key_t, typename EX::node_key_t>::value, "Exchanger node_key_t must be of NODE::key_t type");
   typedef std::pair<const node_key_t *, typename NODE::weight_t> data_pair_t;
   // pair hash function
   struct hashFunc
   {
      typedef size_t result_type;
      typedef data_pair_t argument_type;
      result_type operator()(const argument_type &p) const noexcept { return hashKey(p.first); }
      std::hash<const node_key_t *> hashKey;
   };
   typedef std::unordered_set<data_pair_t, hashFunc> data_t;

   // local
   DATAMANAGER<data_t> bnd_data;
   DATAMANAGER<size_t> lenght;
   size_t local_lenght = 1;
   size_t next_local_lenght = 2;
   // q to store edge visited for curent generation
   std::deque<std::pair<typename NODE::child_t, const typename NODE::key_t *>> q;
   // qn to store edge visited for next generation
   std::deque<std::pair<typename NODE::child_t, const typename NODE::key_t *>> qn;

   class BFSIM
   {
     public:
      // traits
      typedef typename EX::information_key_t information_key_t;
      typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
      typedef xtool::send_only_keys_communication_trait communication_trait;
      // methods
      void setInfo(const xtool::xMpiOutputBuffer &buff, int receivedfrom)
      {
         // unpack node key
         const node_key_t *nk;
         buff.unPack(&nk, 1, MPI_AINT);
         // unpack length
         size_t nl;
         buff.unPack(&nl, 1, xtool::xMPIDataType<size_t>());
         if (nl > local_lenght)
         {
            std::cout << " received length " << nl << " from " << receivedfrom << " is greater then current local_lenght !!"
                      << std::endl;
            throw - 632;
         }
         // unpack data size
         size_t s;
         buff.unPack(&s, 1, xtool::xMPIDataType<size_t>());
         // unpack data
         for (size_t i = 0; i < s; ++i)
         {
            exchanger.unPack(buff, nk, nl, receivedfrom);
            // retrieve current node
            const NODE *node = graph.getData(*nk);
            if (node)
            {
               // functor normally applied by unpack method of echange because only client
               // know what to exchange with remotes
               if (!node->isVisited())
               {
                  node->resetVisited(true);
                  // store child if any in next generation queue
                  for (auto pair : node->rangeChild())
                  {
                     auto &l = lenght.setData(*pair.first);
                     if (l > next_local_lenght) lenght.setData(*pair.first) = next_local_lenght;
                     qn.push_back(std::make_pair(pair, nk));
                  }
               }
            }
            else
            {
               std::cout << "A child is not present in the graph !?" << std::endl;
               throw - 345;
            }
         }
      }
      void getInfo(information_key_t key, xtool::xMpiInputBuffer &buff, int sendto)
      {
         // pack remote node key
         const node_key_t *ro = exchanger.getRemoteNodeKey(key, sendto);
         buff.pack(&ro, 1, MPI_AINT);
         // pack length
         const node_key_t *nk = exchanger.toNodeKeyType(key);
         size_t *pnl = lenght.getData(*nk);
         assert(pnl);
         buff.pack(pnl, 1, xtool::xMPIDataType<size_t>());
         // pack user data
         const data_t *pdata = bnd_data.getData(*nk);
         assert(pnl);
         size_t s = pdata->size();
         buff.pack(&s, 1, xtool::xMPIDataType<size_t>());
         for (auto &data : *pdata)
         {
            exchanger.pack(buff, nk, data.first, data.second, *pnl, sendto);
         }
      }
      size_t getApproxDataSize(void) { return sizeof(typename NODE::child_t) + sizeof(size_t) + exchanger.getApproxDataSize(); }
      BFSIM(EX &exchanger_, DATAMANAGER<size_t> &lenght_, const DATAMANAGER<data_t> &bnd_data_, const DATAMANAGER<NODE> &graph_,
            std::deque<std::pair<typename NODE::child_t, const node_key_t *>> &qn_, const size_t &local_lenght_,
            const size_t &next_local_lenght_)
          : exchanger(exchanger_),
            lenght(lenght_),
            bnd_data(bnd_data_),
            graph(graph_),
            qn(qn_),
            local_lenght(local_lenght_),
            next_local_lenght(next_local_lenght_)
      {
      }

     private:
      EX &exchanger;
      DATAMANAGER<size_t> &lenght;
      const DATAMANAGER<data_t> &bnd_data;
      const DATAMANAGER<NODE> &graph;
      std::deque<std::pair<typename NODE::child_t, const typename NODE::key_t *>> &qn;
      const size_t &local_lenght;
      const size_t &next_local_lenght;

      data_style_trait dummy1;     // remove -Wunused-local-typedefs warnning
      communication_trait dummy2;  // remove -Wunused-local-typedefs warnning
   };
   BFSIM info_man(exchanger, lenght, bnd_data, graph, qn, local_lenght, next_local_lenght);

   // set length flag and visited
   for (auto &key : xtool::make_range(graph.beginKey(), graph.endKey()))
   {
      const NODE *node = graph.getData(*key);
      assert(node);
      node->resetVisited();
      lenght.setData(*key) = std::numeric_limits<size_t>::max();
   }

   // only if key given for source
   if (key_source)
   {
      // apply functor on key_source for 0 weight artificial entering edge
      func_node(std::make_pair(key_source, xtool::xDataType<typename NODE::weight_t>::zero()), nullptr, local_lenght);

      // init queue with child edge of source
      const NODE *source = graph.getData(*key_source);
      assert(source);
      for (auto pair : source->rangeChild())
      {
         qn.push_back(std::make_pair(pair, key_source));
         lenght.setData(*pair.first) = local_lenght;
      }
   }

   unsigned char active = 1;

   // loop on generation
   while (active)
   {
      // swap queue
      q = qn;
      qn.clear();

      exchanger.clearInfoKey();

      // loop till queue is empty
      while (!q.empty())
      {
         // retrieve current edge from parent
         auto &curent = q.front();

         // node key
         auto &nk = *curent.first.first;

         // current lenght
         size_t *pnl = lenght.getData(nk);
         assert(pnl);
         size_t nl = *pnl;

         if (nl <= local_lenght)
         {
            // retrieve current node
            const NODE *node = graph.getData(nk);
            if (node)
            {
               // apply functor to current
               bool ok = func_node(curent.first, curent.second, nl);
               if (ok)
               {
                  if (!node->isVisited())
                  {
                     node->resetVisited(true);
                     // store child if any in next generation queue
                     for (auto pair : node->rangeChild())
                     {
                        auto &l = lenght.setData(*pair.first);
                        if (l > next_local_lenght) lenght.setData(*pair.first) = next_local_lenght;
                        qn.push_back(std::make_pair(pair, &nk));
                     }
                     if (node->hasRemoteChild())
                     {
                        data_t &data = bnd_data.setData(nk);
                        data.insert(std::make_pair(curent.second, curent.first.second));
                        exchanger.toExchange(nk);
                     }
                  }
               }
               else
               {
                  // postponed by user function decision
                  q.push_back(curent);
               }
            }
            else
            {
               std::cout << "A child is not present in the graph !?" << std::endl;
               throw - 345;
            }
         }
         else
         {
            // postpone
            q.push_back(curent);
         }
         // delete current
         q.pop_front();
      }

      // exchanging
      exchanger.exchange(info_man);
      exchanger.exchange(qn.empty(), active);

      local_lenght = next_local_lenght;
      ++next_local_lenght;
   }
}

/// Base class concept of exchanger for breadthFirstSearchDist function. Abstract class
template <typename PM, typename KM, typename NODE>
class breadthFirstSearchDistExchanger
{
  public:
   typedef typename KM::information_key_t information_key_t;
   typedef typename NODE::key_t node_key_t;
   breadthFirstSearchDistExchanger(PM &pm_) : pm(pm_), world(pm.getComm()), km(pm), kc(world) {}
   void clearInfoKey(void) { kc.clearKeys(); }
   template <typename IM>
   void exchange(IM &im)
   {
      xtool::exchangeInformation(kc, im);
   }
   void exchange(bool empty, unsigned char &active)
   {
      active = (empty) ? 0 : 1;
      MPI_Allreduce(MPI_IN_PLACE, &active, 1, MPI_UNSIGNED_CHAR, MPI_MAX, world);
   }
   virtual const node_key_t *toNodeKeyType(information_key_t k) { return static_cast<const node_key_t *>(k); }
   virtual const node_key_t *getRemoteNodeKey(information_key_t key, int sendto)
   {
      return static_cast<const node_key_t *>(pm.getConstPartitionObject(*key).getRemoteObjectOn(sendto));
   }
   void toExchange(const node_key_t &nk) { kc.accumulateKeys(nk, km); }
   virtual size_t getApproxDataSize() = 0;
   virtual void pack(xtool::xMpiInputBuffer &buff, const node_key_t *nk, const node_key_t *nks, typename NODE::weight_t w,
                     size_t &l, int sendto) = 0;
   virtual void unPack(const xtool::xMpiOutputBuffer &buff, const node_key_t *nk, size_t &l, int receivedfrom) = 0;

  private:
   PM &pm;
   MPI_Comm world;
   KM km;
   xtool::xKeyContainerSendOrRecv<information_key_t> kc;
};

}  // end namespace
#endif
