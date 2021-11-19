#include <assert.h>
#include "xEdge.h"



namespace xinterface
{

  namespace xmeshinterface
  {

    xEdge::xEdge(const xEntity& e) : xEntity(e) 
    {
      if(isValid())      assert ( e.getType() == eType::EDGE );
    };
   
    eType xEdge::getType() const
    {
      return eType::EDGE;
    };
    
  }
}
