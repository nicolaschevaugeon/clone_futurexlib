#ifndef _MIT_XEDGE_
#define _MIT_XEDGE_

#include <assert.h>
#include "xEntity.h"

namespace xinterface
{

  namespace xmeshinterface
  {
 
    class xEdge: public xEntity
    {
    public:
      friend class xMeshQueryInterface;
      ~xEdge(){};
      xEdge(const xEntity& e);
      xEdge()=default;
      xEdge(const xEdge&)=default;
      //  xEdge(xEdge&&)=default;
   
   
      eType getType() const ;
    public:
      template< typename T > 
    xEdge( const T& thing , const xMeshQueryInterface& _mi ): xEntity(thing, _mi)  { };

    } ;                                      


  }
}                                       


#endif


