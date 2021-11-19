#ifndef _MIT_XVERTEX_
#define _MIT_XVERTEX_

#include "xEntity.h"
#include "xPoint.h"
#include "xMeshQueryInterface.h"

 
namespace xinterface{

  namespace xmeshinterface{


    class xVertex : public xEntity {
    public:
      template< typename T > 
	xVertex( const T& thing , const xMeshQueryInterface& _mi ): xEntity(thing , _mi) {};
    xVertex( const xEntity& ent): xEntity(ent) {};
    xVertex(): xEntity() {};
      xVertex(const xVertex&)=default;
      // xVertex(xVertex&&)=default;
   
 
      xtensor::xPoint  point() const ;
      void   print() const ; 


    };



  } // namespace xmeshinterface
} // namespace xinterface

#endif

