/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef __ENTITY_TO_ENTITY__H
#define __ENTITY_TO_ENTITY__H

#include <iostream>
#include <string>
#include "mEntity.h"
#include "xEntityFilter.h"

namespace xfem
{
class xMesh;
typedef std::function<AOMD::mEntity*(AOMD::mEntity*)> xEntityToEntity;

struct xUpperAdjacency
{
   AOMD::mEntity* operator()(AOMD::mEntity* e);
};

// Same as xUpperAdjacency but
// if there is a partition in the upper adjacency (usually volume element)
// it returns the element in that partition which is connected to entity e instead of
// the upper adjacency itself.
// This is needed in the implementation of Neumann BC when elements close
// to the boundary are cut by an interface.
class xUpperAdjacencyPartition
{
  public:
   AOMD::mEntity* operator()(AOMD::mEntity* e);
};

struct xCreator
{
   xCreator();
   AOMD::mEntity* operator()(AOMD::mEntity* e);
};

struct xCreatorRecursive
{
   xCreatorRecursive();
   AOMD::mEntity* operator()(AOMD::mEntity* e);
};

struct xUpperCreator
{
   xUpperCreator();
   AOMD::mEntity* operator()(AOMD::mEntity* e);
};

struct xPartitionCreatorRecursive
{
   xPartitionCreatorRecursive();
   AOMD::mEntity* operator()(AOMD::mEntity* e);
};

struct xUpperCreatorFiltered : public xUpperCreator
{
   xUpperCreatorFiltered(const xEntityFilter f) : xUpperCreator(), filter(f) {}
   AOMD::mEntity* operator()(AOMD::mEntity* e);

  protected:
   xEntityFilter filter;
};

/// return the entity of level dim that created the one given to operator(mEntity *)
/* dim as a default value of 3, change it for 2d problem  */
struct xUpperCreatorRecursive
{
   xUpperCreatorRecursive(int dim = 3);
   AOMD::mEntity* operator()(AOMD::mEntity* e);

  private:
   int dim;
};

class xClassifyOn
{
  public:
   xClassifyOn(const std::string& z);
   AOMD::mEntity* operator()(AOMD::mEntity* e);

  private:
   std::string zone_name;
   int zone_id;
};

/*class xClassifyOnCreatorIfNotClassified {
public:
   xClassifyOnCreatorIfNotClassified();
   AOMD::mEntity* operator()(AOMD::mEntity* e);
private:
   std::string zone_name;
   unsigned int zone_tag;
   unsigned int was_created_by_tag;
};*/

/// this  more general then xCreator. It gives any attached entity to e with tag
struct xAttached
{
   xAttached(unsigned int tag);
   AOMD::mEntity* operator()(AOMD::mEntity* e);

  protected:
   unsigned int tag;
};

}  // namespace xfem

#endif
