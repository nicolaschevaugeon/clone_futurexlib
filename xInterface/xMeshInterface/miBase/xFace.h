#ifndef _MIT_XFACE_
#define _MIT_XFACE_

#include "xEntity.h"

namespace xinterface{

  namespace xmeshinterface{

 
    class xFace: public xEntity
    {
    public:
    friend class xMeshInterface;
        ~xFace() {};
      xFace(const xEntity& e) ;
      xFace()=default;
      xFace(const xFace&)=default;
      //  xFace(xFace&&)=default;
   
     public:
      template< typename T > 
    xFace( const T& thing , const xMeshQueryInterface& _mi ): xEntity(thing, _mi)  { };
    };


    class xTri: public xFace
    {
    public:
      ~xTri() {};
      xTri(const xFace& e);
      eType getType();
       };

    class xQuad: public xFace
    {
      ~xQuad() {};
      xQuad(const xFace& e) ;
      eType getType();
     };

  }                                       

}                                       

#endif


