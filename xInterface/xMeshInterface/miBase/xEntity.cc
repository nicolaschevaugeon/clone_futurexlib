#include "xEntity.h" 
#include "xVertex.h" 
#include "xMeshQueryInterface.h"
#include "xEdge.h"
#include "xFace.h"
#include "xSolid.h"


namespace xinterface{

  namespace xmeshinterface{



    eType   xEntity::getType() const { return pQuery->getType( entity_identifier) ; };

    const any&     xEntity::getEntityIdentifier() const { return  entity_identifier ; };

    int     xEntity::getId()     const { return pQuery->getId( entity_identifier) ;   };

    int     xEntity::size(int i) const { return pQuery->size ( entity_identifier , i ) ;   };

    xVertex xEntity::getVertex( int j) const { return pQuery->getVertex(entity_identifier ,j ) ; }

    xEdge xEntity::getEdge(   int j) const { return pQuery->getEdge(  entity_identifier ,j ) ; }

    xFace xEntity::getFace(   int j) const { return pQuery->getFace(  entity_identifier ,j ) ; }

    xSolid xEntity::getSolid(  int j) const { return pQuery->getSolid(  entity_identifier ,j ) ; }

    xEntity xEntity::get(int i, int j) const { 
      switch(i) {
      case 0 : return getVertex(j); // retrun a xVertex
      case 1 : return getEdge(j); 
      case 2 : return getFace(j);
      case 3 : return getSolid(j);
      default:  
	try                     { throw std::logic_error( "ERROR: function xEntity::get(i,j) is impossible. The first parameter must be between 0 and 3.  " );  } 
	catch ( std::exception & e ) { std::cerr << e.what()<<std::endl; exit(0) ;} 
      }
    }

    xClassification  xEntity::getClassification() const { return pQuery->getClassification( entity_identifier ) ; }  

    int xEntity::isAdjacencyCreated(int level) const { return pQuery->isAdjacencyCreated( entity_identifier , level ) ;  } 


    bool xEntity::operator==(const xEntity& other) const {
      if ( getTagOfQuery()   != other.getTagOfQuery())   return false;
      return ( getUniqueAddress() == other.getUniqueAddress() ) ;
    };

    bool xEntity::operator!=(const xEntity& other) const {
      if (getTagOfQuery()   != other.getTagOfQuery())  return true;
      return   ( getUniqueAddress() != other.getUniqueAddress() );
    };

    bool xEntity::operator<(const xEntity& other) const {
      const unsigned int t1= getTagOfQuery();
      const unsigned int t2=  other.getTagOfQuery();
      if (t1 <  t2 )   return true;
      if (t1>t2)       return false;
      return ( getUniqueAddress() < other.getUniqueAddress() ) ;
    };

    bool xEntity::operator>(const xEntity& other) const {
      const unsigned int t1= getTagOfQuery();
      const unsigned int t2=  other.getTagOfQuery();
      if (t1 >  t2 )   return true;
      else if (t1<t2) return false;
      return ( getUniqueAddress() > other.getUniqueAddress() ) ;
    };
/*
    void* xEntity::getAttachedPointer(const unsigned int& tag ) const      { return pQuery->getAttachedPointer( entity_identifier , tag ) ;  }

    void  xEntity::attachPointer(     const unsigned int& tag, void* data) const { return pQuery->attachPointer(      entity_identifier , tag , data ) ;  }

    void  xEntity::deleteData(        const unsigned int& tag )     const       { return pQuery->deleteData(         entity_identifier , tag ) ;  }
*/



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // MI-independent functions
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const xMeshQueryInterface* xEntity::getQuery()      const { return pQuery ;            }

    int xEntity::getLevel() const
    {
      switch ( getType() ) {
      case eType::VERTEX : return 0  ; break;
      case eType::EDGE   : return 1  ; break;
      case eType::TRI    : return 2  ; break;
      case eType::QUAD   : return 2  ; break;
      case eType::HEX    : return 3  ; break;
      case eType::PRISM  : return 3  ; break;
      case eType::TET    : return 3  ; break;
      case eType::MIRROR : return 0  ; break;
      default : abort();
      }
      abort();
    };


    std::string xEntity::printType() const { 
      switch ( getType() )    {
      case eType::VERTEX : return "VERTEX"; break;
      case eType::EDGE   : return "EDGE"  ; break;
      case eType::TRI    : return "TRI"   ; break;
      case eType::QUAD   : return "QUAD"  ; break;
      case eType::HEX    : return "HEX"   ; break;
      case eType::PRISM  : return "PRISM" ; break;
      case eType::TET    : return "TET"   ; break;
      default : abort();
      }
      abort();
    } 


    void xEntity::print() const
    {      
      std::cout << "[print]  " << printType() << "  ID " << getId() << "  Addr. "<< getUniqueAddress() <<std::endl;
      std::cout << "         Query: \""<<getQuery()->getName() <<"\" = tag: "<< getTagOfQuery() <<std::endl;
      std::cout << "         ";
      if ( getType() == eType::VERTEX ) {
	std::cout<< "with XYZ: "  << (static_cast< const xVertex *>(this))->point()  <<  std::endl; 
      }
      else if ( getType() == eType::EDGE ) {
	std::cout<< "with VERTICES " ;
	for(int j=0 ; j < size(0) ; ++j){ std::cout<< ( getVertex(j) ).getId() << " "; };
    std::cout <<std::endl;
      }
      else if ( getType() == eType::TRI || getType() == eType::QUAD  ) {
	std::cout<< "with EDGE " ;
	for(int j=0 ; j < size(1) ; ++j){ std::cout << ( getEdge(j) ).getId() << " "; }; 
	std::cout<< "with VERTICES " ;
	for(int j=0 ; j < size(0) ; ++j){ std::cout << ( getVertex(j) ).getId() << " "; };
	std::cout <<std::endl;
      }
      else { 
	std::cout<< "with FACES " ;
	for(int j=0 ; j < size(2) ; ++j){ std::cout << ( getFace(j) ).getId() << " "; }
	std::cout << "with VERTICES " ;
	for(int j=0 ; j < size(0) ; ++j){ std::cout << ( getVertex(j) ).getId() << " "; };   
	std::cout <<std::endl;
    
      }
    }

    void xEntity::printGhost() const
    {      
      std::cout << "[print]  "<< " Addr. "<< getUniqueAddress()<< ", pQuery: "<< getQuery();
      std::cout<< ", tag_of_query: "<< getTagOfQuery() <<std::endl;
    }


    std::vector<xVertex> xEntity::getVertices() const {return pQuery->getVertices( entity_identifier) ;}
    std::vector<xVertex> &xEntity::getVertices(std::vector <xVertex > &vertices ) const {return pQuery->getVertices( entity_identifier, vertices) ;}

    std::vector<xEdge>   xEntity::getEdges()   const {return pQuery->getEdges   ( entity_identifier) ;};
    std::vector<xFace>   xEntity::getFaces()   const {return pQuery->getFaces   ( entity_identifier) ;};
    std::vector<xSolid>  xEntity::getSolids()  const {return pQuery->getSolids  ( entity_identifier) ;};
    std::vector<xSolid>&  xEntity::getSolids(std::vector <xSolid > &solids)  const {return pQuery->getSolids  ( entity_identifier, solids) ;};




    void* xEntity::getUniqueAddress() const {
      /*if((!uniqueAddress)&&pQuery) return pQuery->getUniqueAddressFromXentity(*this);
      if (uniqueAddress) return uniqueAddress;
      return nullptr;*/
      if ( uniqueAddress ) return uniqueAddress;
      if(!pQuery) return nullptr;
      return pQuery->getUniqueAddressFromXentity(*this);
    }

    unsigned int xEntity::getTagOfQuery() const { 
     /* if ((tagOfQuery == 0)&& (pQuery)) return pQuery->getTag();
      if( tagOfQuery  ) return tagOfQuery;
      return 0;*/
     if ( tagOfQuery != 0 ) return tagOfQuery;
     if( !pQuery ) return 0;
     return pQuery->getTag();
    }

    // ghost entity
    xEntity::xEntity(const void* a , const unsigned int t ){
      pQuery=nullptr;
      entity_identifier =0;
      uniqueAddress=const_cast<void*>(a);
      tagOfQuery=t;
      pQuery=nullptr;
    };


    bool xEntity::isGhost() const 
    {
      if ( isValid()       ) return false;
      if ( tagOfQuery==0   ) return false;
      return !uniqueAddress;
      //if ( uniqueAddress==0) return false;
      //return true;
    };

    bool xEntity::isValid() const 
    {
      if( pQuery==nullptr     ) return false;
      if( entity_identifier==0) return false;
      return true;
    };

    xtensor::xPoint xEntity::getCentroid() const 
    {
      return pQuery->getCentroid( entity_identifier );
    };
  } // namepsace xmeshinterface
} // namespace xinterface
