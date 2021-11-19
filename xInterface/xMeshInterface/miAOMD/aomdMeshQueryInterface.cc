#include "aomdMeshQueryInterface.h"

#include <stdio.h>

#include <memory>
#include <vector>

#include "AOMDfwd.h"
#include "GEntity.h"
#include "aomdAttachablePointer.h"
#include "mAOMD.h"
#include "mEdge.h"
#include "mEntity.h"
#include "mFace.h"
#include "mHex.h"
#include "mMesh.h"
#include "mRegion.h"
#include "mTet.h"
#include "mVertex.h"
#include "modeler.h"
#include "mpi.h"
#include "xEntityIterators.h"
#include "xPoint.h"

namespace xinterface
{
namespace xmeshinterface
{
// aomdMeshQueryInterface::aomdMeshQueryInterface(): xMeshQueryInterface() { cout<< "aomdMeshQueryInterface() constuctor"<<endl; }

int aomdMeshQueryInterface::dim() const { return getMeshSize(3) ? 3 : (getMeshSize(2) ? 2 : (getMeshSize(1) ? 1 : 0)); }

eType aomdMeshQueryInterface::getType(const any& entity_identifier) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(entity_identifier);
   switch (pe->getType())
   {
      case AOMD::mEntity::mType::VERTEX:
         return eType::VERTEX;
         break;
      case AOMD::mEntity::mType::EDGE:
         return eType::EDGE;
         break;
      case AOMD::mEntity::mType::TRI:
         return eType::TRI;
         break;
      case AOMD::mEntity::mType::QUAD:
         return eType::QUAD;
         break;
      case AOMD::mEntity::mType::TET:
         return eType::TET;
         break;
      case AOMD::mEntity::mType::HEX:
         return eType::HEX;
         break;
      case AOMD::mEntity::mType::PRISM:
         return eType::PRISM;
         break;
      case AOMD::mEntity::mType::PYRAMID:
         return eType::PYRAMID;
         break;
      default:
         throw;
   }
};

int aomdMeshQueryInterface::getId(const any& entity_identifier) const
{
   return (any_cast<AOMD::mEntity*>(entity_identifier))->getId();
};

xtensor::xPoint aomdMeshQueryInterface::getCoordinates(const any& entity_identifier) const
{
   xtensor::xPoint _pt = (any_cast<AOMD::mVertex*>(entity_identifier))->point();
   return _pt;
};

int aomdMeshQueryInterface::size(const any& entity_identifier, int i) const
{
   return (any_cast<AOMD::mEntity*>(entity_identifier))->size(i);
};

xVertex aomdMeshQueryInterface::getVertex(const any& entity_identifier, int j) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(entity_identifier);
   if (pe->size(0) == 0)
   {
      cout << "Unable to evaluate getVertex(" << j << ") of the entity. Please execute :" << endl;
      cout << "     mesh.modifyState(0," << pe->getLevel() << ",true) " << endl;
      cout << " and mesh.modifyState(" << pe->getLevel() << ",0,true) " << endl;
   }
   return xVertex(pe->get(0, j), *this);
};

xEdge aomdMeshQueryInterface::getEdge(const any& entity_identifier, int j) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(entity_identifier);
   if (pe->size(1) == 0)
   {
      cout << "Unable to evaluate getVertex(" << j << ") of the entity. Please execute :" << endl;
      cout << "     mesh.modifyState(1," << pe->getLevel() << ",true) " << endl;
      cout << " and mesh.modifyState(" << pe->getLevel() << ",1,true) " << endl;
   }
   return xEdge(pe->get(1, j), *this);
};
;

xFace aomdMeshQueryInterface::getFace(const any& entity_identifier, int j) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(entity_identifier);
   if (pe->size(2) == 0)
   {
      cout << "Unable to evaluate getFace(" << j << ") of the entity. Please execute :" << endl;
      cout << "     mesh.modifyState(2," << pe->getLevel() << ",true) " << endl;
      cout << " and mesh.modifyState(" << pe->getLevel() << ",2,true) " << endl;
   }
   return xFace(pe->get(2, j), *this);
};

xSolid aomdMeshQueryInterface::getSolid(const any& entity_identifier, int j) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(entity_identifier);
   if (pe->size(3) == 0)
   {
      cout << "Unable to evaluate getVertex(" << j << ") of the entity. Please execute :" << endl;
      cout << "     mesh.modifyState(3," << pe->getLevel() << ",true) " << endl;
      cout << " and mesh.modifyState(" << pe->getLevel() << ",3,true) " << endl;
   }
   return xSolid(pe->get(3, j), *this);
};

xClassification aomdMeshQueryInterface::getClassification(const any& entity_identifier) const
{
   AOMD::mEntity* pent = any_cast<AOMD::mEntity*>(entity_identifier);
   AOMD::pGEntity pge = pent->getClassification();
   return xClassification(static_cast<void*>(pge), pge->tag(), pge->dim());
};

int aomdMeshQueryInterface::isAdjacencyCreated(const any& entity_identifier, int level) const
{
   AOMD::mEntity* pent = any_cast<AOMD::mEntity*>(entity_identifier);
   return pent->isAdjacencyCreated(level);
};

bool aomdMeshQueryInterface::equalEntityIdentifier(const any& entity_identifier, const any& other_identifier) const
{
   return any_cast<AOMD::mEntity*>(entity_identifier) == any_cast<AOMD::mEntity*>(other_identifier);
};

bool aomdMeshQueryInterface::lessthanEntityIdentifier(const any& entity_identifier, const any& other_identifier) const
{
   return any_cast<AOMD::mEntity*>(entity_identifier) < any_cast<AOMD::mEntity*>(other_identifier);
};

void* aomdMeshQueryInterface::getAttachedPointer(const any& entity_identifier, const unsigned int& tag) const
{
   AOMD::mEntity* pent = any_cast<AOMD::mEntity*>(entity_identifier);
   auto data = pent->getData(tag);
   if (data)
   {
      aomdAttachablePointer* attachedPt = static_cast<aomdAttachablePointer*>(data);
      return static_cast<void*>(attachedPt->getPointer());
   }
   return nullptr;
}

void aomdMeshQueryInterface::attachPointer(const any& entity_identifier, const unsigned int& tag, void* data) const
{
   AOMD::mEntity* pent = any_cast<AOMD::mEntity*>(entity_identifier);
   aomdAttachablePointer* attachedPt = (aomdAttachablePointer*)pent->getData(tag);
   if (!attachedPt)
   {
      attachedPt = new aomdAttachablePointer(data);
      pent->attachData(tag, attachedPt);
   }
   else
   {
      attachedPt->resetPointer(data);
   }
}

void aomdMeshQueryInterface::deleteAttachment(const any& entity_identifier, const unsigned int& tag) const
{
   AOMD::mEntity* pent = any_cast<AOMD::mEntity*>(entity_identifier);
   pent->deleteData(tag);
}

int aomdMeshQueryInterface::getMeshSize(int i) const { return mesh.size(i); }

std::vector<xVertex> aomdMeshQueryInterface::getVertices(const any& entity_identifier) const
{
   std::vector<xVertex> vertices;
   return getVertices(entity_identifier, vertices);
};
std::vector<xVertex>& aomdMeshQueryInterface::getVertices(const any& entity_identifier, std::vector<xVertex>& vertices) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(entity_identifier);
   const int size = pe->size(0);
   const size_t needed_size = vertices.size() + size;
   if (vertices.capacity() < needed_size) vertices.reserve(needed_size);
   for (auto i = 0; i < size; ++i) vertices.emplace_back(pe->get(0, i), *this);
   return vertices;
};
std::vector<xEdge> aomdMeshQueryInterface::getEdges(const any& entity_identifier) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(entity_identifier);
   const int size = pe->size(1);
   std::vector<xEdge> edges(size);
   for (auto i = 0; i < size; ++i)
   {
      edges[i] = xEdge(pe->get(1, i), *this);
   }
   return edges;
};
std::vector<xFace> aomdMeshQueryInterface::getFaces(const any& entity_identifier) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(entity_identifier);
   const int size = pe->size(2);
   std::vector<xFace> faces(size);
   for (auto i = 0; i < size; ++i)
   {
      faces[i] = xFace(pe->get(2, i), *this);
   }
   return faces;
};
std::vector<xSolid> aomdMeshQueryInterface::getSolids(const any& entity_identifier) const
{
   std::vector<xSolid> solids;
   return getSolids(entity_identifier, solids);
};

std::vector<xSolid>& aomdMeshQueryInterface::getSolids(const any& entity_identifier, std::vector<xSolid>& solids) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(entity_identifier);
   const int size = pe->size(3);
   const size_t needed_size = solids.size() + size;
   if (solids.capacity() < needed_size) solids.reserve(needed_size);
   for (auto i = 0; i < size; ++i) solids.emplace_back(pe->get(3, i), *this);
   return solids;
};

// ITERATORS

xVertex aomdMeshQueryInterface::fromVertexIteratorIdentifier(const any& iterator_identifier) const
{
   AOMD::mMesh::iter it = any_cast<AOMD::mMesh::iter>(iterator_identifier);
   return xVertex((*it), *this);
};
any aomdMeshQueryInterface::nextVertexIteratorIdentifier(const any& iterator_identifier) const
{
   AOMD::mMesh::iter it = any_cast<AOMD::mMesh::iter>(iterator_identifier);
   ++it;
   return it;
};
bool aomdMeshQueryInterface::equalVertexIteratorIdentifier(const any& first, const any& second) const
{
   return ((any_cast<AOMD::mMesh::iter>(first)) == (any_cast<AOMD::mMesh::iter>(second)));
};

xEntity aomdMeshQueryInterface::fromEdgeIteratorIdentifier(const any& iterator_identifier) const
{
   AOMD::mMesh::iter it = any_cast<AOMD::mMesh::iter>(iterator_identifier);
   return xEntity((*it), *this);
};
any aomdMeshQueryInterface::nextEdgeIteratorIdentifier(const any& iterator_identifier) const
{
   AOMD::mMesh::iter it = any_cast<AOMD::mMesh::iter>(iterator_identifier);
   ++it;
   return it;
};
bool aomdMeshQueryInterface::equalEdgeIteratorIdentifier(const any& first, const any& second) const
{
   return ((any_cast<AOMD::mMesh::iter>(first)) == (any_cast<AOMD::mMesh::iter>(second)));
};

xEntity aomdMeshQueryInterface::fromFaceIteratorIdentifier(const any& iterator_identifier) const
{
   AOMD::mMesh::iter it = any_cast<AOMD::mMesh::iter>(iterator_identifier);
   return xEntity((*it), *this);
};
any aomdMeshQueryInterface::nextFaceIteratorIdentifier(const any& iterator_identifier) const
{
   AOMD::mMesh::iter it = any_cast<AOMD::mMesh::iter>(iterator_identifier);
   ++it;
   return it;
};
bool aomdMeshQueryInterface::equalFaceIteratorIdentifier(const any& first, const any& second) const
{
   return ((any_cast<AOMD::mMesh::iter>(first)) == (any_cast<AOMD::mMesh::iter>(second)));
};

xEntity aomdMeshQueryInterface::fromSolidIteratorIdentifier(const any& iterator_identifier) const
{
   AOMD::mMesh::iter it = any_cast<AOMD::mMesh::iter>(iterator_identifier);
   return xEntity((*it), *this);
};
any aomdMeshQueryInterface::nextSolidIteratorIdentifier(const any& iterator_identifier) const
{
   AOMD::mMesh::iter it = any_cast<AOMD::mMesh::iter>(iterator_identifier);
   ++it;
   return it;
};
bool aomdMeshQueryInterface::equalSolidIteratorIdentifier(const any& first, const any& second) const
{
   return ((any_cast<AOMD::mMesh::iter>(first)) == (any_cast<AOMD::mMesh::iter>(second)));
};

xVertexIterator aomdMeshQueryInterface::beginVertex() const { return xVertexIterator(*this, mesh.begin(0)); };
xEdgeIterator aomdMeshQueryInterface::beginEdge() const { return xEdgeIterator(*this, mesh.begin(1)); };
xFaceIterator aomdMeshQueryInterface::beginFace() const { return xFaceIterator(*this, mesh.begin(2)); };
xSolidIterator aomdMeshQueryInterface::beginSolid() const { return xSolidIterator(*this, mesh.begin(3)); };
xVertexIterator aomdMeshQueryInterface::endVertex() const { return xVertexIterator(*this, mesh.end(0)); };
xEdgeIterator aomdMeshQueryInterface::endEdge() const { return xEdgeIterator(*this, mesh.end(1)); };
xFaceIterator aomdMeshQueryInterface::endFace() const { return xFaceIterator(*this, mesh.end(2)); };
xSolidIterator aomdMeshQueryInterface::endSolid() const { return xSolidIterator(*this, mesh.end(3)); };

xEntity aomdMeshQueryInterface::find(xEntity e) const
{
   AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(e.getEntityIdentifier());
   AOMD::mEntity* found = mesh.find(pe);
   if (!found) return xEntity();
   return xEntity(found, *this);
};

xtensor::xPoint aomdMeshQueryInterface::getCentroid(const any& entity_identifier) const
{
   AOMD::mEntity* pent = any_cast<AOMD::mEntity*>(entity_identifier);
   Trellis_Util::mPoint pt = pent->getCentroid();
   return xtensor::xPoint(pt(0), pt(1), pt(2));
};

xVertex aomdMeshQueryInterface::getMeshVertexFromId(int ID) const { return xVertex(mesh.getVertex(ID), *this); };

// ---------------  nested translateToPartitionManagerXentity --------------------

void aomdMeshQueryInterface::translateToPartitionManagerXentity(xPartitionManagerXentity& part_man_XENT,
                                                                const partmanAOMD_t& part_man_AOMD)
{
   // MPI_Comm _comm = part_man_AOMD.getComm();

   // the new partition manager (to be return)
   // xPartitionManagerXentity part_man_XENT(*this, _comm);

   // for all element of the previous partition manager
   for (auto it = part_man_AOMD.beginObject(); it != part_man_AOMD.endObject(); ++it)
   {
      // the key of the current element
      AOMD::mEntity* penti_AOMD = *it;

      // the key of the element in the new partition manager
      xEntity enti_XENT(penti_AOMD, *this);

      // takes the partition object (i.e. here it is the collection of remote copies of the entity)
      auto po_XENT = part_man_XENT.getPartitionObject(enti_XENT);

      // takes/create the remote copies container (=xPartitionObject) in the new partition manager
      auto po_AOMD = part_man_AOMD.getConstPartitionObject(*penti_AOMD);

      // for all the remote objects (remote copies of the object)
      for (auto ro : po_AOMD.getRemoteObjectsCollectionRange())
      {
         // insert corresponding remote object in the new partition manager

         //!!!!! ici insert s'attend Ã  recevoir un xEntity* car dÃ©pend d'un xPartitionManager< DataManager<Xentity> >
         //   po_XENT.insertUntyped(   ro.getProcId(),   ro.getUntypedObjectAddress()   );
         po_XENT.insert(ro.getProcId(), static_cast<const xEntity*>(static_cast<const void*>(ro.getObjectAddress())));
      }
   }
   // return part_man_XENT;
}

}  // namespace xmeshinterface
}  // namespace xinterface
