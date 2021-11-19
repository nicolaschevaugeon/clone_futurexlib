/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/
#ifndef _MESHINTERFACEAOMDMMESH_
#define _MESHINTERFACEAOMDMMESH_

#include "linearalgebra3d.h"
#include "mAOMD.h"
#include "mMesh.h"
#include "mVertex.h"
#include "mFace.h"
#include "mEdge.h"
#include "mRegion.h"

/// Define the interface between AOMD mMesh end the Fast Marching Library
/*!
  All the function defined here need to be defined for the MESHINTERFACE class passed to the fast marching interface.
*/

namespace xfastmarching
{

  template <class IT>
    class entitytovertexiteratorconvertor{
  public:
    typedef const AOMD::mVertex * value_type;
  entitytovertexiteratorconvertor(IT _it):it(_it){}
    entitytovertexiteratorconvertor<IT> operator++(int){
      entitytovertexiteratorconvertor<IT> rit(it);
      ++(*this);
      return rit;
    }
    entitytovertexiteratorconvertor<IT>& operator++(){
      ++it;
      return *this;
    }
    bool operator!=(const  entitytovertexiteratorconvertor& other){
      return it != other.it;
    }
    void operator=(value_type v){
      it=const_cast<AOMD::mVertex*>(v);
    }
    const AOMD::mVertex * operator*(){
      return static_cast<const AOMD::mVertex *>(*it);
    }
  private:
    IT it;
  };

  class meshinterfaceAOMD{
  public:
  meshinterfaceAOMD(const AOMD::mMesh& _m):m(_m){}
    typedef AOMD::mVertex vertex;
    typedef AOMD::mEdge   edge;
    typedef AOMD::mFace   face;
    typedef AOMD::mRegion region;
    typedef AOMD::mEntity entity;
    const AOMD::mMesh  &m;
  };

  inline std::size_t getTag(const meshinterfaceAOMD &mi, const std::string &tagname){
    return AOMD::AOMD_Util::Instance()->lookupMeshDataId (tagname.c_str());
  };

  inline std::size_t getNewTag(const meshinterfaceAOMD &mi, const std::string &tagname){
    return AOMD::AOMD_Util::Instance()->newMeshDataId (tagname.c_str());
  }

  inline void releaseTag(const meshinterfaceAOMD &mi, const std::string &tagname){
    int tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId (tagname.c_str());
    AOMD::AOMD_Util::Instance()->deleteMeshDataId (tag);
  }

  inline void releaseTag(const meshinterfaceAOMD &mi, const std::size_t &tag){
    AOMD::AOMD_Util::Instance()->deleteMeshDataId (tag);
  }

  class AttachableDataPointer : public AOMD::mAttachableData{
  public:
  AttachableDataPointer(void * pnt_):pnt(pnt_){}
    ~AttachableDataPointer() override = default;
    void *pnt;
  };

  inline void * getAttachedDataPointer(const meshinterfaceAOMD &mi, AOMD::mEntity &en, size_t tag){
    AttachableDataPointer * dataatt = static_cast< AttachableDataPointer *> (en.getData(tag));
    if (dataatt) return dataatt->pnt;
    else return nullptr;
  }

  inline void deleteData(const meshinterfaceAOMD &mi, AOMD::mEntity &en, size_t tag){
    en.deleteData(tag);
  }
  inline void attachDataPointer(const meshinterfaceAOMD &mi, AOMD::mEntity &en, size_t tag, void * data){
    AttachableDataPointer * dataatt = static_cast< AttachableDataPointer *> ( en.getData(tag));
    if (dataatt){ dataatt->pnt = data; return;}
    dataatt = new  AttachableDataPointer(data);
    en.attachData(tag, dataatt);
  }

  inline void * getAttachedDataPointer(const meshinterfaceAOMD &mi, AOMD::mVertex &v, size_t tag){
    AttachableDataPointer * dataatt = static_cast< AttachableDataPointer *> (v.getData(tag));
    if (dataatt) return dataatt->pnt;
    else return nullptr;
  }

  inline void attachDataPointer(const meshinterfaceAOMD &mi, AOMD::mVertex &v, size_t tag, void * data){
    AttachableDataPointer * dataatt = static_cast< AttachableDataPointer *> ( v.getData(tag));
    if (dataatt){ dataatt->pnt = data; return;}
    dataatt = new  AttachableDataPointer(data);
    v.attachData(tag, dataatt);
  }

  inline void deleteData(const meshinterfaceAOMD &mi, AOMD::mVertex &v, size_t tag){
    v.deleteData(tag);
  }

  inline void * getAttachedDataPointer(const meshinterfaceAOMD &mi, AOMD::mFace &v, size_t tag){
    AttachableDataPointer * dataatt = static_cast< AttachableDataPointer *> (v.getData(tag));
    if (dataatt) return dataatt->pnt;
    else return nullptr; 
  }

  inline void attachDataPointer(const meshinterfaceAOMD &mi, AOMD::mFace &v, size_t tag, void * data){
    AttachableDataPointer * dataatt = static_cast< AttachableDataPointer *> ( v.getData(tag));
    if (dataatt){ dataatt->pnt = data; return;}
    dataatt = new  AttachableDataPointer(data);
    v.attachData(tag, dataatt);
  }

  inline void deleteData(const meshinterfaceAOMD &mi, AOMD::mFace &v, size_t tag){
    v.deleteData(tag);
  }

  template <class VERTEXINPUTITERATOR> 
    void getVerticesNeighbors(const meshinterfaceAOMD &mi, const AOMD::mVertex &v, VERTEXINPUTITERATOR it ){
    size_t nbedges = v.size(1);
    for (size_t i = 0; i < nbedges; ++i){
      AOMD::mEntity * e = v.get( 1, i );
      for (size_t k = 0; k < 2; ++k){
        const AOMD::mVertex *vik = static_cast<const AOMD::mVertex *  > (e->get(0, k));
        if (vik != &v) *it++ = vik;
      }
    }
  }


  template <class EDGEINPUTITERATOR>
    void getEdges(const meshinterfaceAOMD &mi, const AOMD::mVertex &v, EDGEINPUTITERATOR it){
    size_t n = v.size(1);
    for (int i=0; i < n; ++i){
      auto* e=v.get(1,i);
      *it=static_cast<AOMD::mEdge*>(e); 
      ++it;
    }
  }

  template <class FACEINPUTITERATOR>
    void getFaces(const meshinterfaceAOMD &mi, const AOMD::mVertex &v, FACEINPUTITERATOR it){
    size_t n = v.size(2);
    for (int i=0; i < n; ++i){
      auto* e=v.get(2,i);
      *it=static_cast<AOMD::mFace*>(e); 
      ++it;
    }
  }

  template <class REGIONINPUTITERATOR>
    void getRegions(const meshinterfaceAOMD &mi, const AOMD::mVertex &v, REGIONINPUTITERATOR it){
    size_t n = v.size(3);
    for (int i=0; i < n; ++i){
      auto* e=v.get(3,i);
      *it=static_cast<AOMD::mRegion*>(e); 
      ++it;
    }
  }


  template <class FACEINPUTITERATOR>
    void getFaces (const meshinterfaceAOMD &mi, const AOMD::mEdge &v, FACEINPUTITERATOR it){
    size_t n = v.size(2);
    for (int i=0; i < n; ++i){
      auto* e=v.get(2,i);
      *it=static_cast<AOMD::mFace*>(e); 
      ++it;
    }
  }

  template <class REGIONINPUTITERATOR>
    void getRegions (const meshinterfaceAOMD &mi, const AOMD::mEdge &v, REGIONINPUTITERATOR it){
    size_t n = v.size(3);
    for (int i=0; i < n; ++i){
      auto* e=v.get(3,i);
      *it=static_cast<AOMD::mRegion*>(e); 
      ++it;
    }
  }


  template <class VERTEXINPUTITERATOR>
    void getVertices (const meshinterfaceAOMD &mi, const AOMD::mEdge &v, VERTEXINPUTITERATOR it){
    size_t n = v.size(0);
    for (int i=0; i < n; ++i){
      auto* e=v.get(0,i);
      *it=static_cast<AOMD::mVertex*>(e); 
      ++it;
    }
  }

  template <class VERTEXINPUTITERATOR>
    void getVertices (const meshinterfaceAOMD &mi, const AOMD::mFace &v, VERTEXINPUTITERATOR it){
    size_t n = v.size(0);
    for (int i=0; i < n; ++i){
      auto* e=v.get(0,i);
      *it=static_cast<AOMD::mVertex*>(e); 
      ++it;
    }
  }

  template <class VERTEXINPUTITERATOR>
    void getVertices (const meshinterfaceAOMD &mi, const AOMD::mRegion &v, VERTEXINPUTITERATOR it){
    size_t n = v.size(0);
    for (int i=0; i < n; ++i){
      auto* e=v.get(0,i);
      *it=static_cast<AOMD::mVertex*>(e); 
      ++it;
    }
  }


  inline double getX(const meshinterfaceAOMD &mi, const AOMD::mVertex &v){
    return v.point()(0);
  }

  inline double getY(const meshinterfaceAOMD &mi, const AOMD::mVertex &v){
    return v.point()(1);
  }

  inline double getZ(const meshinterfaceAOMD &mi, const AOMD::mVertex &v){
    return v.point()(2);
  }

  inline int getId(const meshinterfaceAOMD &mi, const AOMD::mVertex &v){
    return v.getId();
  }

  template <class GEOMVECT>
    double getCoord(const meshinterfaceAOMD &mi, const AOMD::mVertex &v, GEOMVECT& x){
    std::cout << "getcoord not implemented for the used type " << std::endl;
    assert(0);
  }

  template <class T>
    double getCoord(const meshinterfaceAOMD &mi, const AOMD::mVertex &v, vector2d<T>& x){
    x={v.point()(0), v.point()(1)};
  }

  template <class T>
    double getCoord(const meshinterfaceAOMD &mi, const AOMD::mVertex &v, vector3d<T>& x){
    x={v.point()(0), v.point()(1), v.point()(2)};
  }

} // end namespace
#endif
