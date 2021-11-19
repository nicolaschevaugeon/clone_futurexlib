#ifndef _M_ATTACHABLEPOINTER_H_
#define _M_ATTACHABLEPOINTER_H_


#include "mEntity.h"
#include "mAttachableDataContainer.h"
 
namespace xinterface{

  namespace xmeshinterface{

    class aomdAttachablePointer : public AOMD::mAttachableData
    {
    public :
    aomdAttachablePointer(void* _p) : data(_p) {}
      virtual ~aomdAttachablePointer (){};
      void  resetPointer(void* _p) { data=_p;}
      void* getPointer() { return data;}
    private:
      void* data;
    };

  } // namespace xmeshinterface
} // namespace xinterface

 
#endif
