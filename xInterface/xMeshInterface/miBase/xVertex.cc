#include "xMeshQueryInterface.h" 
#include "xVertex.h" 



namespace xinterface{

  namespace xmeshinterface{


    xtensor::xPoint xVertex::point() const { return pQuery->getCoordinates( entity_identifier ) ; }; 

    void xVertex::print() const {
      std::cout << "[print]  " << printType() ;
      std::cout<< " ID " << getId() ;
      if ( getType() == eType::VERTEX ) { std::cout<< " XYZ: "<< point() <<std::endl;
      };
    };


  } // namespace xmeshinterface
} // namespace xinterface



