#ifndef _MI_AOMD_QUERY_H_
#define _MI_AOMD_QUERY_H_

#include <boost/function.hpp>

#include "AOMDfwd.h"
#include "GEntity.h"
#include "mAOMD.h"
#include "mAttachableDataContainer.h"
#include "mEdge.h"
#include "mEntity.h"
#include "mFace.h"
#include "mHex.h"
#include "mMesh.h"
#include "mRegion.h"
#include "mTet.h"
#include "mVertex.h"
#include "modeler.h"
#include "xAttachedDataManagerAOMD.h"
#include "xEntity.h"
#include "xMeshQueryInterface.h"
#include "xPartitionManager.h"
#include "xPartitionManagerXentity.h"
#include "xVertex.h"

namespace xinterface
{
namespace xmeshinterface
{
using namespace xinterface::xmeshinterface;

using any = small_any;

// only to translate xAttachedDataManagerAOMD into xPartitionManagerXentity with funtion translateToPartitionManagerXentity after
// mesh distribution by xDistMesh
using partmanAOMD_t = xtool::xPartitionManager<xinterface::aomd::xAttachedDataManagerAOMD>;

class aomdMeshQueryInterface : public xMeshQueryInterface
{
  private:
   AOMD::mMesh& mesh;

   //    protected:
   //      aomdMeshQueryInterface() ;
  public:
   aomdMeshQueryInterface(AOMD::mMesh& _m) : mesh(_m){};
   //    aomdMeshQueryInterface( AOMD::mMesh* _pm ) :  mesh(_pm) { /*cout<< "aomdMeshQueryInterface() constructor from a
   //    mMesh*"<<endl;*/ };
   ~aomdMeshQueryInterface(){/* cout<<" ~aomdMeshQueryInterface()"<<endl;*/};

   inline AOMD::mMesh& getMesh() { return mesh; };
   int dim() const;
   eType getType(const any& entity_identifier) const;
   int getId(const any& entity_identifier) const;
   xtensor::xPoint getCoordinates(const any& entity_identifier) const;
   int size(const any& entity_identifier, int i) const;
   xVertex getVertex(const any& entity_identifier, int j) const;
   xEdge getEdge(const any& entity_identifier, int j) const;
   xFace getFace(const any& entity_identifier, int j) const;
   xSolid getSolid(const any& entity_identifier, int j) const;
   xClassification getClassification(const any& entity_identifier) const;
   int isAdjacencyCreated(const any& entity_identifier, int level) const;
   // for comparison operators for xValKey
   bool equalEntityIdentifier(const any& entity_identifier, const any& other_identifier) const;
   bool lessthanEntityIdentifier(const any& entity_identifier, const any& other_identifier) const;

   void* getAttachedPointer(const any& entity_identifier, const unsigned int& tag) const;
   void attachPointer(const any& entity_identifier, const unsigned int& tag, void* data) const;
   void deleteAttachment(const any& entity_identifier, const unsigned int& tag) const;

   /// for partition manager and used to compare entity
   inline void* getUniqueAddressFromXentity(const xEntity e) const override { return any_cast<void*>(e.getEntityIdentifier()); };

   inline xEntity getXentityFromUniqueAddress(const void* add) const override
   {
      AOMD::mEntity* pent = static_cast<AOMD::mEntity*>(const_cast<void*>(add));
      return xEntity(pent, *this);
   };

   int getMeshSize(int i) const;

   std::vector<xVertex> getVertices(const any& entity_identifier) const override;
   std::vector<xVertex>& getVertices(const any& entity_identifier, std::vector<xVertex>& vertices) const override;
   std::vector<xEdge> getEdges(const any& entity_identifier) const override;
   std::vector<xFace> getFaces(const any& entity_identifier) const override;
   std::vector<xSolid> getSolids(const any& entity_identifier) const override;
   std::vector<xSolid>& getSolids(const any& entity_identifier, std::vector<xSolid>& solids) const override;

   void translateToPartitionManagerXentity(xPartitionManagerXentity& _pmx, const partmanAOMD_t& pm);

   // for iterator

   xVertex fromVertexIteratorIdentifier(const any& iterator_identifier) const;
   any nextVertexIteratorIdentifier(const any& iterator_identifier) const;
   bool equalVertexIteratorIdentifier(const any& first, const any& second) const;

   xEntity fromEdgeIteratorIdentifier(const any& iterator_identifier) const;
   any nextEdgeIteratorIdentifier(const any& iterator_identifier) const;
   bool equalEdgeIteratorIdentifier(const any& first, const any& second) const;

   xEntity fromFaceIteratorIdentifier(const any& iterator_identifier) const;
   any nextFaceIteratorIdentifier(const any& iterator_identifier) const;
   bool equalFaceIteratorIdentifier(const any& first, const any& second) const;

   xEntity fromSolidIteratorIdentifier(const any& iterator_identifier) const;
   any nextSolidIteratorIdentifier(const any& iterator_identifier) const;
   bool equalSolidIteratorIdentifier(const any& first, const any& second) const;

   xVertexIterator beginVertex() const;
   xEdgeIterator beginEdge() const;
   xFaceIterator beginFace() const;
   xSolidIterator beginSolid() const;
   xVertexIterator endVertex() const;
   xEdgeIterator endEdge() const;
   xFaceIterator endFace() const;
   xSolidIterator endSolid() const;

   xEntity find(xEntity e) const;
   xtensor::xPoint getCentroid(const any& entity_identifier) const;
   xVertex getMeshVertexFromId(int id) const;

   /* ATTIC

   //______________ nested class aomdMeshQueryInterface :: miMeshEntityIterator ______________
   public:
   class miMeshVertexIterator {
   public:
   miMeshVertexIterator( const aomdMeshQueryInterface  &_pmi , const AOMD::mMesh::iter &_it);
   xVertex operator*();
   miMeshVertexIterator  operator++(int);
   miMeshVertexIterator &operator++();
   bool operator!=(const miMeshVertexIterator &in);
   bool operator==(const miMeshVertexIterator &in);
   private:
   AOMD::mMesh::iter it;
   const aomdMeshQueryInterface* pQuery;
   };

   class miMeshEntityIterator {
   public:
   miMeshEntityIterator( const aomdMeshQueryInterface  &_pmi , const AOMD::mMesh::iter &_it);
   xEntity operator*();
   miMeshEntityIterator  operator++(int);
   miMeshEntityIterator &operator++();
   bool operator!=(const miMeshEntityIterator &in);
   bool operator==(const miMeshEntityIterator &in);
   private:
   AOMD::mMesh::iter it;
   const aomdMeshQueryInterface* pQuery;
   };
   // specific to AOMD:
   using miMeshEdgeIterator = miMeshEntityIterator;
   using miMeshFaceIterator = miMeshEntityIterator;
   using miMeshSolidIterator = miMeshEntityIterator;


   xPartitionManagerXentity  translateToPartitionManagerXentity(const partmanAOMD_t & pm) ;


   //____________________________________________________________________________________

   */

};  // end of aomdMeshQueryInterface

//____________________________________________________________________________________

// ATTIC
//******************************************************************************************************
//---------------- specialization of template functions for aomdMeshQueryInterface -------

/*

  template<>
  inline xVertexIterator const meshVertexIteratorBegin(const aomdMeshQueryInterface& mi)
  {
  return  xVertexIterator( mi, (mi.mesh)->begin(0) );
  }

  template<>
  inline xVertexIterator const meshVertexIteratorEnd(const aomdMeshQueryInterface& mi)
  {
  return xVertexIterator( mi , (mi.mesh)->end(0) );
  }

  template<>
  inline xEdgeIterator const meshEdgeIteratorBegin(const aomdMeshQueryInterface& mi)
  {
  return  xEdgeIterator( mi, (mi.mesh)->begin(1) );
  }

  template<>
  inline xEdgeIterator const meshEdgeIteratorEnd(const aomdMeshQueryInterface& mi)
  {
  return xEdgeIterator( mi , (mi.mesh)->end(1) );
  }

  template<>
  inline xFaceIterator const meshFaceIteratorBegin(const aomdMeshQueryInterface& mi)
  {
  return  xFaceIterator( mi, (mi.mesh)->begin(2) );
  }

  template<>
  inline xFaceIterator const meshFaceIteratorEnd(const aomdMeshQueryInterface& mi)
  {
  return xFaceIterator( mi , (mi.mesh)->end(2) );
  }

  template<>
  inline xSolidIterator const meshSolidIteratorBegin(const aomdMeshQueryInterface& mi)
  {
  cout<<"meshSolidIteratorBegin(const aomdMeshQueryInterface& mi)"<<endl;
  return  xSolidIterator( mi, (mi.mesh)->begin(3) );
  }

  template<>
  inline xSolidIterator const meshSolidIteratorEnd(const aomdMeshQueryInterface& mi)
  {
  cout<<"meshSolidIteratorEnd(const aomdMeshQueryInterface& mi)"<<endl;
  return xSolidIterator( mi , (mi.mesh)->end(3) );
  }


*/

//******************************************************************************************************

}  // namespace xmeshinterface
}  // namespace xinterface

#endif
