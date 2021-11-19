/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _NODEANDCONNECTEDEDGE_H
#define _NODEANDCONNECTEDEDGE_H
#include <deque>
#include <unordered_map>
#include "xDataType.h"
#include "xIteratorTools.h"
namespace xgraph
{
//==============================================================================================
/*!
 * nodeFrom class
 * A simple class that represent nodes of a directed graph.  It store incoming weighted edges.
 * Having a node of this graph you can only access to its parent nodes. Number of parent is
 * limited by NB. NB is supposed to be unique for all graph. This way no dynamic allocation
 * holds but there will be memory waste for the nodes  with less then NB parents.
 * The wall graph is supposed to be stored in a exlibris DATAMANAGER container and as such
 * KEYTYPE correspond to the key type used with this container.
 * ET corresponds to the type of weigh associated to incoming edges (i.e. the ones connecting
 * this node to its parents)
 * NT corresponds to the type of data associated to node represented by nodeFrom instance
 * Note that from a structural point of view edges are incoming edges only because associated algorithm
 * treat them as such.
 */
template <typename KEYTYPE, typename ET, typename NT, size_t NB>
class nodeFrom
{
  public:
   typedef KEYTYPE key_t;
   typedef ET weight_t;
   typedef NT node_data_t;
   typedef std::pair<const KEYTYPE *, ET> edge_t;
   typedef std::integral_constant<size_t, NB> NB_t;
   typedef typename std::array<edge_t, NB> data_container_t;
   typedef typename data_container_t::const_iterator const_iterator;
   nodeFrom() : it_end(const_cast<const data_container_t &>(from_data).begin())  // empty init
   {
      std::fill(from_data.begin(), from_data.end(), std::make_pair(nullptr, xtool::xDataType<ET>::zero()));
   }
   nodeFrom(const nodeFrom &other)
       : from_data(other.from_data),
         it_end(const_cast<const data_container_t &>(from_data).begin() + std::distance(other.from_data.begin(), other.it_end))
   {
      NT *od = other.attached_data.get();
      if (od) attached_data.reset(new NT(*od));
   }
   nodeFrom &operator=(const nodeFrom &other)
   {
      from_data = other.from_data;
      it_end = const_cast<const data_container_t &>(from_data).begin() + std::distance(other.from_data.begin(), other.it_end);
      NT *od = other.attached_data.get();
      if (od) attached_data.reset(new NT(*od));

      return *this;
   }
   template <typename K>
   void set(K &keys, ET *vals)
   {
      size_t i = 0;
      for (auto key : keys)
      {
         assert(i < NB);
         auto &p = from_data[i];
         p.first = key;
         p.second = vals[i];
         ++i;
      }
      it_end = (const_cast<const data_container_t &>(from_data)).begin() + i;
   }
   void set(const KEYTYPE *&key, ET val)
   {
      from_data[0].first = key;
      from_data[0].second = val;
      it_end = (const_cast<const data_container_t &>(from_data)).begin() + 1;
   }
   xtool::xRange<const_iterator> getParentEdges() const { return xtool::make_range(from_data.begin(), it_end); }
   void setAttachedData(const node_data_t &data)
   {
      if (attached_data)
      {
         *attached_data = data;
      }
      else
      {
         attached_data.reset(new NT(data));
      }
   }
   const node_data_t *getAttachedData() const { return attached_data.get(); }
   void delAttachedData() { attached_data.reset(nullptr); }

  private:
   data_container_t from_data;
   const_iterator it_end;
   std::unique_ptr<NT> attached_data;
   template <class DATAMANAGER, typename NODE>
   friend std::array<const NODE *, NODE::NB_t::value> &&getParent(const DATAMANAGER &graph, const NODE &node);
};

// A function that gives a array of pointer to parent nodes instance.
template <class DATAMANAGER, typename NODE>
std::array<const NODE *, NODE::NB_t::value> &&getParent(const DATAMANAGER &graph, const NODE &node)
{
   // static_assert(std::is_same< typename DATAMANAGER::data_t, NODE >::value, " must be derived from AOMD::mEntity") ;
   std::array<const NODE *, NODE::NB_t::value> parents;
   auto it = parents.begin();
   for_each(node.from_data.begin(), node.it_end, [&it, &graph, &node](const typename NODE::edge_t e) {
      const NODE *pnode = graph.getData(*(e.first));
      if (pnode)
      {
         *(it++) = pnode;
      }
   });
   return std::move(parents);
}
//==============================================================================================
/*!
 * nodeTo class
 * A simple class that represent nodes of a directed graph.  It store outgoing weighted edges.
 * Having a node of this graph you can only access to its child nodes. Number of child is
 * unlimited.
 * The wall graph is supposed to be stored in a exlibris DATAMANAGER container and as such
 * KEYTYPE correspond to the key type used with this container.
 * ET corresponds to the type of weigh associated to outgoing edges (i.e. the ones connecting
 * this node to its child)
 * NT corresponds to the type of data associated to node represented by nodeTo instance
 * Note that from a structural point of view edges are outgoing edges only because associated algorithm
 * treat them as such.
 */
template <typename KEYTYPE, typename ET, typename NT>
class nodeTo
{
  public:
   typedef KEYTYPE key_t;
   typedef ET weight_t;
   typedef NT node_data_t;
   typedef std::pair<const KEYTYPE *, ET> child_t;
   typedef std::unordered_map<const KEYTYPE *, ET> child_container_t;
   void insertChild(const KEYTYPE *k, const ET &w) { to_childs.insert(std::make_pair(k, w)); }
   size_t removeChild(const KEYTYPE *k)
   {
      to_childs.erase(k);
      // return number of child after removal if any
      return (to_childs.size());
   }
   size_t nbChild() const { return to_childs.size(); }
   auto rangeChild() const -> xtool::xRange<typename child_container_t::const_iterator>
   {
      return xtool::make_range(to_childs.begin(), to_childs.end());
   }
   void resetVisited(bool val = false) { visited = val; }
   void resetVisited(bool val = false) const { visited = val; }
   bool isVisited() const { return visited; }

   void setAttachedData(const node_data_t &data)
   {
      if (attached_data)
      {
         *attached_data = data;
      }
      else
      {
         attached_data.reset(new NT(data));
      }
   }
   const node_data_t *getAttachedData() const { return attached_data.get(); }
   void delAttachedData() { attached_data.reset(nullptr); }
   void setHasRemoteChild(bool setting) { has_remote_child = setting; }
   bool hasRemoteChild() const { return has_remote_child; }

  private:
   std::unique_ptr<NT> attached_data;
   child_container_t to_childs;
   bool has_remote_child = false;
   mutable bool visited;
};

// run a function for each node with no child
template <typename NODE, typename DATAMANAGER, typename F>
void apllyToNoChild(const DATAMANAGER &graph, F func_node)
{
   for (auto &key : xtool::make_range(graph.beginKey(), graph.endKey()))
   {
      const NODE *node = graph.getData(*key);
      assert(node);
      if (!node->nbChild()) func_node(key, node);
   }
}
// run a function for each node fulfilling condition
template <typename NODE, typename DATAMANAGER, typename C, typename F>
void apllyToConditionalNode(const DATAMANAGER &graph, F func_node)
{
   for (auto &key : xtool::make_range(graph.beginKey(), graph.endKey()))
   {
      const NODE *node = graph.getData(*key);
      assert(node);
      if (!C(key, node)) func_node(key, node);
   }
}
/// A breath first search algorithm dedicated to nodeTo graph traversal via edges
/*!
 * Tree traversal respect breath first search algorithm by passing in review nodes that are
 * topologically closer from source first. All children of a node are traversed via its outgoing
 * edges before passing to children of the children. Algorithm insure that edges are traversed only once.
 *
 * NODE the type of nodes of the graph (nodeTo<KEYTYPE,ET,NT>)
 * DATAMANAGER the data manager type that hold NODE nodes associated to KEYTYPE keys.
 * F a function that is run for every edges traversed. Its argument are the 2 nodes it connect (as KEYTYPE).
 *   As the edge is oriented the source is the last argument and the target is in the first pair member of
 *   the first argument. The second pair member is the weight associated to the edge.
 *
 */
template <typename NODE, class DATAMANAGER, typename F>
void breadthFirstSearch(const DATAMANAGER &graph, const typename NODE::key_t *key_source, F func_node)
{
   // reset visited flag
   for (auto &key : xtool::make_range(graph.beginKey(), graph.endKey()))
   {
      const NODE *node = graph.getData(*key);
      assert(node);
      node->resetVisited();
   }

   // q to store edge visited
   std::deque<std::pair<typename NODE::child_t, const typename NODE::key_t *>> q;

   // apply functor on key_source for 0 weight artificial entering edge
   func_node(std::make_pair(key_source, xtool::xDataType<typename NODE::weight_t>::zero()), nullptr);

   // init queue with child edge of source
   const NODE *source = graph.getData(*key_source);
   assert(source);
   for (auto pair : source->rangeChild())
   {
      q.push_back(std::make_pair(pair, key_source));
   }

   // loop till queue is empty
   while (!q.empty())
   {
      // retrieve current edge from parent
      auto &curent = q.front();
      // retrieve current node
      const NODE *node = graph.getData(*curent.first.first);
      if (node)
      {
         // apply functor to current
         func_node(curent.first, curent.second);
         if (!node->isVisited())
         {
            node->resetVisited(true);
            // store child if any in queue
            for (auto pair : node->rangeChild())
            {
               q.push_back(std::make_pair(pair, curent.first.first));
            }
         }
      }
      else
      {
         std::cout << "A child is not present in the graph !?" << std::endl;
         throw - 345;
      }
      // delete treated node
      q.pop_front();
   }
}
/// A breath first search algorithm dedicated to nodeTo graph traversal via edges under condition
/*!
 * Tree traversal respect breath first search algorithm by passing in review nodes that are
 * topologically closer form source first. All children of a node are traversed via its outgoing
 * edges before passing to children of the children. Algorithm insure that edges are traversed only once.
 * A specific condition is imposed by F callback (see bellow) that accept or reject traversal through an edge.
 * If F return always true this function behave exactly like breadthFirstSearch function.
 * If F return false edge re-enter the traversal queue. It is user responsibility to insure that all edges
 * traversed are validated by F once and in correct order of the traversal.
 *
 * NODE the type of nodes of the graph (nodeTo<KEYTYPE,ET,NT>)
 * DATAMANAGER the data manager type that hold NODE nodes associated to KEYTYPE keys.
 * F a function that is run for every edges traversed. Its argument are the 2 nodes it connect (as KEYTYPE).
 *   As the edge is oriented the source is the last argument and the target is in the first pair member of
 *   the first argument. The second pair member is the weight associated to the edge.
 *   This function return true or false to accept or postpone edge treatments.
 *
 */
template <typename NODE, class DATAMANAGER, typename F>
void conditionalBreadthFirstSearch(const DATAMANAGER &graph, const typename NODE::key_t *key_source, F func_node)
{
   // reset visited flag
   for (auto &key : xtool::make_range(graph.beginKey(), graph.endKey()))
   {
      const NODE *node = graph.getData(*key);
      assert(node);
      node->resetVisited();
   }

   // q to store edge visited
   std::deque<std::pair<typename NODE::child_t, const typename NODE::key_t *>> q;

   // apply functor on key_source for 0 weight artificial entering edge
   func_node(std::make_pair(key_source, xtool::xDataType<typename NODE::weight_t>::zero()), nullptr);

   // init queue with child edge of source
   const NODE *source = graph.getData(*key_source);
   assert(source);
   for (auto pair : source->rangeChild())
   {
      q.push_back(std::make_pair(pair, key_source));
   }

   // loop till queue is empty
   while (!q.empty())
   {
      // retrieve current edge from parent
      auto &curent = q.front();
      // retrieve current node
      const NODE *node = graph.getData(*curent.first.first);
      if (node)
      {
         // apply functor to current
         bool ok = func_node(curent.first, curent.second);
         if (ok)
         {
            if (!node->isVisited())
            {
               node->resetVisited(true);
               // store child if any in queue
               for (auto pair : node->rangeChild())
               {
                  q.push_back(std::make_pair(pair, curent.first.first));
               }
            }
         }
         else
            q.push_back(curent);
      }
      else
      {
         std::cout << "A child is not present in the graph !?" << std::endl;
         throw - 345;
      }
      // delete treated node
      q.pop_front();
   }
}
}  // end namespace
#endif
