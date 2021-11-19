#ifndef _MI_MESH_QUERY_INTERFACE_H_
#define _MI_MESH_QUERY_INTERFACE_H_
#include <boost/any.hpp>
#include <memory>
#include <fstream>
#include <iostream>
#include <vector>

#include "small_any.h"

//#include "xEntity.h"
//#include "xEdge.h"
//#include "xFace.h"
//include "xSolid.h"
#include "xPoint.h"
#include "xClassification.h"
#include "xAttachableData.h"
#include "xEntityType.h"

#include "xStringManager.h"


namespace xinterface{

  namespace xmeshinterface{

  class xEntity;
  class xVertex;
  class xFace;
  class xSolid;
  class xEdge;
  class xEntityIterator;
  class xVertexIterator;
  class xFaceIterator;
  class xSolidIterator;
  class xEdgeIterator;


  void exportGmshEdge (std::ostream&,const xEdge&, int&, int );

    void exportGmshFace (std::ostream&,const xFace&, int&, int );

    void exportGmshSolid(std::ostream&,const xSolid&, int&, int );


    class xMeshQueryInterface  { 
    private:
      xtool::xStringManager<unsigned int,xtool::hashString<unsigned int> > QueryStringManager ;
      unsigned int tag=0;
      std::string name;
    public:
      xMeshQueryInterface()  ;  
    public:
      void setName(std::string _name) ;
      std::string getName() const ;
      inline unsigned int getTag() const   { return tag ;  }
      xEntity getGhostXentityWithRemoteAddress(  const void* add  ) const;
      unsigned int lookupQueryDataId(const std::string& s) ;
      xEntityIterator begin(int what) const;
      xEntityIterator end(int what)   const;
      void printMesh() const;
      void exportGMSH(std::string _filename) const ;  
      void exportGMSH(std::ofstream & _f) const ; 

    public:
      // vitual functions below
      virtual int dim() const =0;
      virtual eType getType( const small_any& entity_identifier ) const =0;
      virtual int     getId( const small_any& entity_identifier ) const =0;
      virtual xtensor::xPoint  getCoordinates( const small_any& entity_identifier ) const =0;
      /// return the number of adjacencies of dimension i
      virtual int     size( const small_any& entity_identifier , int i ) const =0;
      virtual xVertex getVertex( const small_any& entity_identifier , int j) const =0;
      virtual xEdge getEdge  ( const small_any& entity_identifier , int j) const =0;
      virtual xFace getFace  ( const small_any& entity_identifier , int j) const =0;
      virtual xSolid getSolid ( const small_any& entity_identifier , int j) const =0;
      virtual xClassification getClassification( const small_any& entity_identifier ) const =0;
      virtual int isAdjacencyCreated( const small_any& entity_identifier , int level) const =0;

      // comparison operators for xValKey
      virtual bool equalEntityIdentifier(    const small_any& entity_identifier, const  small_any& other_identifier) const =0;
      virtual bool lessthanEntityIdentifier( const small_any& entity_identifier, const  small_any& other_identifier) const =0;

      // --
      virtual void* getAttachedPointer( const small_any& entity_identifier , const unsigned int& tag )              const =0;
      virtual void  attachPointer(      const small_any& entity_identifier , const unsigned int& tag , void* data ) const =0;
      virtual void  deleteAttachment(         const small_any& entity_identifier , const unsigned int& tag )              const =0;

      virtual void*   getUniqueAddressFromXentity(  const xEntity e )  const =0; 
      virtual xEntity getXentityFromUniqueAddress(  const void*   a )  const =0;   
      virtual int   getMeshSize( int i )    const=0 ;
  
      virtual xVertex fromVertexIteratorIdentifier  (const small_any& iterator_identifier) const =0;
      virtual small_any     nextVertexIteratorIdentifier  (const small_any& iterator_identifier) const =0;
      virtual bool    equalVertexIteratorIdentifier (const small_any& first,const small_any& second) const =0;

      virtual xEntity fromEdgeIteratorIdentifier  (const small_any& iterator_identifier) const =0;
      virtual small_any     nextEdgeIteratorIdentifier  (const small_any& iterator_identifier) const =0;
      virtual bool    equalEdgeIteratorIdentifier (const small_any& first,const small_any& second) const =0;

      virtual xEntity fromFaceIteratorIdentifier  (const small_any& iterator_identifier) const =0;
      virtual small_any     nextFaceIteratorIdentifier  (const small_any& iterator_identifier) const =0;
      virtual bool    equalFaceIteratorIdentifier (const small_any& first,const small_any& second) const =0;

      virtual xEntity fromSolidIteratorIdentifier  (const small_any& iterator_identifier) const =0;
      virtual small_any     nextSolidIteratorIdentifier  (const small_any& iterator_identifier) const =0;
      virtual bool    equalSolidIteratorIdentifier (const small_any& first,const small_any& second) const =0;

      virtual xVertexIterator beginVertex() const =0;
      virtual   xEdgeIterator beginEdge()   const =0;
      virtual   xFaceIterator beginFace()   const =0;
      virtual  xSolidIterator beginSolid()  const =0;
      virtual xVertexIterator endVertex()   const =0;
      virtual   xEdgeIterator endEdge()     const =0;
      virtual   xFaceIterator endFace()     const =0;
      virtual  xSolidIterator endSolid()    const =0;

      virtual std::vector<xVertex> getVertices(const small_any& entity_identifier ) const =0;
      virtual std::vector<xVertex> &getVertices(const small_any& entity_identifier, std::vector<xVertex >& vertices ) const = 0;

      virtual std::vector<xEdge>   getEdges   (const small_any& entity_identifier ) const =0;
      virtual std::vector<xFace>   getFaces   (const small_any& entity_identifier ) const =0;
      virtual std::vector<xSolid>  getSolids  (const small_any& entity_identifier ) const =0;
      virtual std::vector<xSolid>  &getSolids  (const small_any& entity_identifier, std::vector<xSolid> &solids ) const =0;


      virtual xEntity find(xEntity e) const =0;
      virtual xtensor::xPoint getCentroid(const small_any& entity_identifier) const =0;
  
      // returns the first ID of the first Vertex for xPostPro.cc
      virtual xVertex getMeshVertexFromId(int ID) const=0; 
 
    };

  } // namepsace xmeshinterface
} // namespace xinterface


#endif




 
