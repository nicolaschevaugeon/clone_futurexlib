#ifndef _MIT_XSOLID_
#define _MIT_XSOLID_


#include "xEntity.h"

namespace xinterface
{

  namespace xmeshinterface
  {


 
    class xSolid: public xEntity
    {
    public:
   friend class xMeshInterface;
         ~xSolid(){};
      xSolid(const xEntity& e) ;
      xSolid()=default;
      xSolid(const xSolid&)=default;
      // xSolid(xSolid&&)=default;
   
    public:
      template< typename T > 
    xSolid( const T& thing , const xMeshQueryInterface& _mi ): xEntity(thing, _mi)  { };
    };


  

    class xTet : public xSolid
    {
    public:
      ~xTet(){};
      xTet(const xSolid& e) ;
      eType getType();
      };
 

    class xHex : public xSolid
    {
    public:
      ~xHex(){};
      xHex(const xSolid& e) ;
      eType getType() ;
      };
 

    class xPyramid : public xSolid
    {
    public:
      ~xPyramid(){};
      xPyramid(const xSolid& e);
      eType getType() ;
       };

 
    class xPrism : public xSolid
    {
    public:
      ~xPrism(){};
      xPrism(const xSolid& e);
      eType getType();
     };


  } // end namespace xmeshinterface

}// end namespace xinterface


namespace std {
  template<>
    inline bool std::equal_to< xinterface::xmeshinterface::xSolid >::operator() ( const xinterface::xmeshinterface::xSolid & e,
                                                                                 const xinterface::xmeshinterface::xSolid & o ) const
    {
      return (e==o);
    }



 template <> struct hash<xinterface::xmeshinterface::xSolid>
  {
    int operator()(const xinterface::xmeshinterface::xSolid & e) const
    {
     return hash<int>()( e.getId() );
    }
  };


}

                                       

#endif





