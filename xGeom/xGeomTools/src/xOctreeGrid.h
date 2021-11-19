/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XOCTREEGRID_H
#define _XOCTREEGRID_H

// std
#include <forward_list>
#include <memory>

// xtensor
#include "xPoint.h"
#include "xVector.h"
// xgeom
#include "xBoundingBox.h"

namespace xgeom
{
/// This  is the template class to represent the octree search data structure.
//! The template parameters are:
//! maxElem is the maximum number of element centeroid that can be found in an octant
//! maxLevel is the maxmimum number of Level of subdivision
//! numDiv is the number of subdivision per axis. ( for a proper octree nuMDiv =2)
//! Each octant bucket of the octree represent a bounding box and entities that covers it are
//! accessible from the buckets. The octree refinment strategie is based on the number of entity
//! whose centroid are inside the bucket.
//! If this number is larger than The template parameter maxElem and the level of the bucket is
//! lower than maxLevel, the bucket is divided and the entities are moved to higher level buckets,
//! recursively, until the maxElem or maxLevel constrains are reached.
template <size_t maxElem = 5, size_t maxLevel = 20, size_t numDiv = 2>
class octantBucket
{
  public:
   struct elem_data
   {
      void *elem = nullptr;                    /* the pointer to a mesh Db region */
      xtensor::xPoint centroid = {0., 0., 0.}; /* centroid of element bounding box inside of the octant */
      xgeom::xBoundingBox bb = {{0., 0., 0.}, {0., 0., 0.}};
      // \note following constructors are not need with recent compiler (gcc9) ... but it need to be spelled out at least for
      // gcc4.8
      elem_data(void *_elem, const xtensor::xPoint &_centroid, const xgeom::xBoundingBox &_bb)
          : elem(_elem), centroid(_centroid), bb(_bb)
      {
      }
      elem_data() = default;
   };
   octantBucket(const xgeom::xBoundingBox &_bb, size_t _level = 0, const octantBucket *_parent = nullptr)
       : bb(_bb), level(_level), parent(_parent)
   {
   }
   octantBucket() = default;
   /// Return the level of the bucket (Number of division since root bucket)
   size_t getLevel() const { return level; }
   /// return the bucket Bounding box
   const xgeom::xBoundingBox &getBoundingBox() const { return bb; }
   /// return the list of data attached to the buckets.
   const std::forward_list<elem_data> &getDatas() const { return covered_elements; }
   /// call fun on all leaves of the octree, (depth first)
   void visitLeaves(std::function<void(const octantBucket &buck)> fun) const;
   /// return the leaf bucket containing p if exist it exist, nullptr otherwise
   const octantBucket<maxElem, maxLevel, numDiv> *findBucket(const xtensor::xPoint &p) const;
   /// Add element to leaves octant bucket's list if edata.bb cover the leaf boundingbox.
   //! if a leaf has more than maxElem element centeroid, the leaf is subdivided and all data are moved to these new leaves
   //! recursivelly
   void addElement(const elem_data &edata);

  private:
   xgeom::xBoundingBox bb;                // bounding box.
   size_t level = 0;                      // level in the octree of the bucket
   const octantBucket *parent = nullptr;  // link to the parent bucket */
   // number of element centeroid contains in the bucket.
   size_t contained_centroids_size = 0;
   // list of element whose BB intersect octant.
   std::forward_list<elem_data> covered_elements;
   static constexpr size_t numChilds = numDiv * numDiv * numDiv;
   // Pointer to all the childs (numBucks ... ) if any.
   std::unique_ptr<std::array<octantBucket, numChilds>> childs;
   /// try to refine
   //! prerequist : bucket has no child
   //! Returns true for success, false for failure (no memory left).
   bool subdivideOctantBucket();
};

using xOctreeGrid = octantBucket<5, 20, 2>;

template <size_t maxElem, size_t maxLevel, size_t numDiv>
void octantBucket<maxElem, maxLevel, numDiv>::visitLeaves(
    std::function<void(const octantBucket<maxElem, maxLevel, numDiv> &buck)> fun) const
{
   if (childs)
   {
      for (const auto &child : *childs) child.visitLeaves(fun);
      return;
   }
   fun(*this);
}

template <size_t maxElem, size_t maxLevel, size_t numDiv>
const octantBucket<maxElem, maxLevel, numDiv> *octantBucket<maxElem, maxLevel, numDiv>::findBucket(const xtensor::xPoint &p) const
{
   if (!bb.contains(p)) return nullptr;
   if (!childs) return this;
   for (const octantBucket &child : *childs)
   {
      const octantBucket *found = child.findBucket(p);
      if (found) return found;
   }
   return nullptr;
}

template <size_t maxElem, size_t maxLevel, size_t numDiv>
void octantBucket<maxElem, maxLevel, numDiv>::addElement(const elem_data &edata)
{
   if (!edata.bb.cover(bb)) return;
   if (childs)
   {
      for (octantBucket &bucket : *childs) bucket.addElement(edata);
      return;
   }
   // At this point, I'am on a leaf that cover the bounding box of elem.
   // checking if edata already in the list. Do nothing if its the case.
   for (const auto &data : covered_elements)
   {
      if (data.elem == edata.elem)
      {
         if (!(data.bb == edata.bb && data.centroid == edata.centroid))
            std::cout << "Warning !!! added same element with diff centeroid and bb" << std::endl;
         return;
      }
   }
   covered_elements.push_front(edata);
   if (bb.contains(edata.centroid)) ++contained_centroids_size;
   if (contained_centroids_size <= maxElem) return;
   if (level < maxLevel)
      subdivideOctantBucket();
   else
      std::cout << "Warning : maxLevel Reach in xOctreeGrid, more than  " << maxElem << " In currrent octant " << __FILE__
                << " line " << __LINE__ << std::endl;
   return;
}

template <size_t maxElem, size_t maxLevel, size_t numDiv>
bool octantBucket<maxElem, maxLevel, numDiv>::subdivideOctantBucket()
{
   if (childs) throw;
   // std::make_unique is c++14
   // childs = std::make_unique<std::array<octantBucket, numChilds>>();
   childs = std::unique_ptr<std::array<octantBucket, numChilds>>(new std::array<octantBucket, numChilds>);
   if (!childs)
   {
      std::cerr << "Error, subdivideOctantBucket could not allocate enough space" << std::endl;
      return false;
   }
   const xtensor::xPoint tmp = (bb.max - bb.min) / numDiv;
   for (size_t k = 0; k < numDiv; ++k)
   {
      size_t shift = k * numDiv * numDiv;
      for (size_t j = 0; j < numDiv; ++j)
      {
         for (size_t i = 0; i < numDiv; ++i)
         {
            octantBucket &child = (*childs)[shift + i];
            child.parent = this;
            child.level = level + 1;
            child.bb.min = bb.min + xtensor::xPoint(tmp[0] * i, tmp[1] * j, tmp[2] * k);
            child.bb.max = bb.min + xtensor::xPoint(tmp[0] * (i + 1), tmp[1] * (j + 1), tmp[2] * (k + 1));
         }
         shift += numDiv;
      }
   }
   contained_centroids_size = 0;
   for (auto &child : *childs)
      for (const elem_data &data : covered_elements) child.addElement(data);
   covered_elements.clear();
   return true;
}

}  // namespace xgeom

#endif
