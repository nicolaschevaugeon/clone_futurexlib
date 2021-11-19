/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _DISTANCE_NEAREST_POINT_GENERATOR_IMP_H
#define _DISTANCE_NEAREST_POINT_GENERATOR_IMP_H
#include <vector>

#include "mEntity.h"
#include "xMapping.h"
namespace xgeom
{
#ifdef HAVE_AOMD
template <typename GEN>
void genAndTagFromCut(GEN &generator, AOMD::mEntity *e, typename GEN::Container_t &container)
{
   // local variable
   int id;
   double min, max;
   const double zero = 0.0;
   std::vector<AOMD::mEntity *> subentities_vect;

   // getting constitutve node of entity e
   xinterface::aomd::getSubentities(e, 0, subentities_vect);

   // getting value of the level set ls for all nodes of entity  e
   std::vector<double> vals = generator.ls.getVals(subentities_vect);

   // min max
   min = *min_element(vals.begin(), vals.end());
   max = *max_element(vals.begin(), vals.end());

   // if the element is related to iso-zero it have to be inspected
   // -------------------------------------------------------------
   if (min * max <= zero)
   {
      // cut reference element of e and check if something was done
      id = generator.cutter(vals);

      // if element is strictly cut
      if (id > 0)
      {
         // define a mapping to declare new nodes in physical coordinate
         xmapping::xMapping *mapping = xfem::xMappingBuilderHolderSingleton::instance().buildMapping(*e);
         if (mapping)
         {
            generator.modifyCoord(mapping);
         }
         else
            throw -1;

         // the element is strictly cut
         // tag the node of the entity
         for_each(subentities_vect.begin(), subentities_vect.end(), generator.tagger);

         // generate geometric entity from points
         generator.addList(container);

         // if 2 geometric entity
         if (id > 1)
         {
            generator.modifyCoord2(mapping);
            generator.addList2(container);
         }

         delete mapping;
      }
      // element e is not cut but touch iso-zero
      // else if (id <0)
      // must be negative as min*max <= zero and id <=0 : id <= zero => min*max = 0 => id <0
      else
      {
         id = -id;
         // if a boundary is entirelly in iso zero generate geometric entity
         if (generator.tagTouchingIso(e, id))
         {
            // define a mapping to declare new nodes in physical coordinate
            xmapping::xMapping *mapping = xfem::xMappingBuilderHolderSingleton::instance().buildMapping(*e);
            if (mapping)
            {
               generator.modifyCoord(mapping);

               // generate geometric entity from points
               generator.addList(container);

               delete mapping;
            }
            else
            {
               throw -1;
            }
         }
      }
   }
   // e is not related to iso-zero
   // ----------------------------
   else
   {
      ;  // nothing to do for now
   }

   return;
}
#endif
}  // namespace xgeom

#endif
