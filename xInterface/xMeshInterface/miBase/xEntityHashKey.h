#ifndef _MESH_GENERIC_ENTITYHASHKEY_H_
#define _MESH_GENERIC_ENTITYHASHKEY_H_

#include "xEntity.h"
#include "xMeshQueryInterface.h"

namespace xinterface{

  namespace xmeshinterface{


    class xEntityHashKey {
    public:
      int operator()(xEntity e) const   {  
	return  e.getId()  ;  }
    };


    class xEntityEqualKey 
    {
    public:
      bool operator()(xEntity e, xEntity o) const {
	return (e==o);
      }
    };



    class xEntityLessThanKey 
    {
    public:
      bool operator()(xEntity e, xEntity o) const {
	return (e<o);
      }
    };

  } // namepsace xmeshinterface
} // namespace xinterface




#endif
