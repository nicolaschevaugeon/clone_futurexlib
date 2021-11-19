#ifndef _MIT_XENTITY_
#define _MIT_XENTITY_

#include <vector>
//#include <type_traits>
#include <string>

#include "small_any.h"
#include "xAttachableData.h"
#include "xClassification.h"
#include "xDataType.h"
#include "xEntityType.h"
#include "xMeshQueryInterface.h"
#include "xPoint.h"

namespace xinterface
{
namespace xmeshinterface
{
using any = small_any;

class xVertex;
class xEdge;
class xFace;
class xSolid;

class xEntity
{
  protected:
   // here, a pointer to a query because it can be nullptr: use fonction isValid() to test. Use isGhost() to test if the query is
   // know through its name (tag)
   const xMeshQueryInterface* pQuery = nullptr;
   any entity_identifier = 0;
   mutable unsigned int tagOfQuery = 0;
   mutable void* uniqueAddress = nullptr;

  public:
   xEntity(const void* a, const unsigned int t);

   template <typename T>
   xEntity(const T& thing, const xMeshQueryInterface& _mi)
       : pQuery(&_mi), entity_identifier(thing), tagOfQuery(0), uniqueAddress(0)
   {
   }

  public:
   xEntity() = default;
   xEntity& operator=(xEntity&& rhs) = default;
   xEntity(xEntity&& rhs) = default;
   xEntity(const xEntity& rhs) = default;
   xEntity& operator=(const xEntity& rhs) = default;

   virtual ~xEntity() {}

   bool isGhost() const;
   bool isValid() const;
   const any& getEntityIdentifier() const;
   eType getType() const;
   int getId() const;
   int size(int i) const;
   xVertex getVertex(int j) const;
   xEdge getEdge(int j) const;
   xFace getFace(int j) const;
   xSolid getSolid(int j) const;
   xEntity get(int i, int j) const;
   xClassification getClassification() const;
   int isAdjacencyCreated(int level) const;
   bool operator==(const xEntity& other) const;
   bool operator!=(const xEntity& other) const;
   bool operator<(const xEntity& other) const;
   bool operator>(const xEntity& other) const;

  public:
   const xMeshQueryInterface* getQuery() const;
   int getLevel() const;
   std::string printType() const;
   void print() const;
   void printGhost() const;
   std::vector<xVertex> getVertices() const;
   std::vector<xVertex>& getVertices(std::vector<xVertex>& vertices) const;

   std::vector<xEdge> getEdges() const;
   std::vector<xFace> getFaces() const;
   std::vector<xSolid> getSolids() const;
   std::vector<xSolid>& getSolids(std::vector<xSolid>& solids) const;

   void* getUniqueAddress() const;
   unsigned int getTagOfQuery() const;
   xtensor::xPoint getCentroid() const;

  public:
   // ATTACH Methods------------------------------------------
   template <typename T>
   void attachData(const unsigned int& tag, const T& objet) const
   {
      void* stored = pQuery->getAttachedPointer(entity_identifier, tag);
      if (stored)
      {
         (static_cast<xAttachableData<T>*>(stored))->setData(objet);
      }
      xAttachableData<T>* to_store = new xAttachableData<T>(objet);
      const_cast<xMeshQueryInterface*>(pQuery)->attachPointer(const_cast<small_any&>(entity_identifier), tag,
                                                              static_cast<void*>(to_store));
   }

   template <typename T>
   T& getData(const unsigned int& tag) const
   {
      void* stored = pQuery->getAttachedPointer(entity_identifier, tag);
      if (!stored)
      {
         std::cout << "ERROR: using xEntity::getData(tag=" << tag << ") while no data is attached." << std::endl;
         std::cout << "       Use function xEntity::hasData(tag) to check if data exists" << std::endl;
         throw;
      }
      return (static_cast<xAttachableData<T>*>(stored))->getData();
   }

   template <typename T>
   void deleteData(const unsigned int& tag) const  // deleteData
   {
      void* stored = pQuery->getAttachedPointer(entity_identifier, tag);
      if (stored) delete static_cast<xAttachableData<T>*>(stored);
      const_cast<xMeshQueryInterface*>(pQuery)->deleteAttachment(const_cast<small_any&>(entity_identifier), tag);
   }

   bool hasData(const unsigned int& tag) const  // hasData
   {
      void* stored = pQuery->getAttachedPointer(entity_identifier, tag);
      if (!stored) return false;
      return true;
   }
};

}  // namespace xmeshinterface
}  // namespace xinterface

namespace xtool
{
/// Specialization of xDatatype for xEntity
template <>
class xDataType<xinterface::xmeshinterface::xEntity>
{
  public:
   static xinterface::xmeshinterface::xEntity zero() { return xinterface::xmeshinterface::xEntity(); }
   static std::string stype() { return "xEntity"; }
};
}  // namespace xtool

namespace std
{
template <>
inline bool std::equal_to<xinterface::xmeshinterface::xEntity>::operator()(const xinterface::xmeshinterface::xEntity& e,
                                                                           const xinterface::xmeshinterface::xEntity& o) const
{
   return (e == o);
}

template <>
struct hash<xinterface::xmeshinterface::xEntity>
{
   int operator()(const xinterface::xmeshinterface::xEntity& e) const { return hash<int>()(e.getId()); }
};

}  // namespace std

#endif
