#include <assert.h>
#include "xFace.h"


namespace xinterface
{

  namespace xmeshinterface
  {

    xFace::xFace(const xEntity& e) : xEntity(e) 
    { 
      if(isValid())  assert ( e.getType() == eType::TRI || e.getType() == eType::QUAD );
    }



    xTri::xTri(const xFace& e) : xFace(e) 
    { 
     if(isValid())   assert ( e.getType() == eType::TRI  );
    }

    eType xTri::getType()
    { 
      return eType::TRI;
    }



    xQuad::xQuad(const xFace& e) : xFace(e) 
    { 
     if(isValid())   assert ( e.getType() == eType::QUAD  );
    }

    eType xQuad::getType()
    { 
      return eType::QUAD;
    }

  }
}
