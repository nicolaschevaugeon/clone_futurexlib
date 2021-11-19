#include <assert.h>
#include "xSolid.h"


namespace xinterface
{

  namespace xmeshinterface
  {

    xSolid::xSolid(const xEntity& e) : xEntity(e) 
    {	
     if(isValid())   assert ( e.getType() == eType::TET || e.getType() == eType::HEX || e.getType() == eType::PYRAMID ||e.getType() == eType::PRISM );
    };
 
  
    xTet::xTet(const xSolid& e) : xSolid(e) 
    { 
     if(isValid())   assert ( e.getType() == eType::TET );
    }

    eType xTet::getType()
    {
      return eType::TET;
    }
 
    xHex::xHex(const xSolid& e) : xSolid(e) 
    { 
     if(isValid())   assert ( e.getType() == eType::HEX );
    }

    eType xHex::getType()
    {
      return eType::HEX;
    }
 
 
    xPyramid::xPyramid(const xSolid& e) : xSolid(e) 
    { 
     if(isValid())   assert ( e.getType() == eType::PYRAMID );
    }

    eType xPyramid::getType()
    {
      return eType::PYRAMID;
    }
 
    xPrism::xPrism(const xSolid& e) : xSolid(e) 
    { 
     if(isValid())   assert ( e.getType() == eType::PRISM );
    }

    eType xPrism::getType()
    {
      return eType::PRISM;
    }


  }                 
}
