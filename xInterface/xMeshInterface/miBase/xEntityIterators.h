#ifndef _MI_XENTITYITERATOR_HEADER_
#define _MI_XENTITYITERATOR_HEADER_ 


#include "xEntity.h" 
#include "xVertex.h" 
#include "xEdge.h" 
#include "xFace.h" 
#include "xSolid.h" 
#include "small_any.h"

namespace xinterface{

  namespace xmeshinterface{

    using any = small_any;
    class xVertexIterator;

    class xEntityIteratorBasic {
    public:
      friend class xVertexIterator;
      friend class xEdgeIterator;
      friend class xFaceIterator;
      friend class xSolidIterator;
      friend class xEntityIterator;
     
    xEntityIteratorBasic( ) : iterator_identifier(0), pQuery(0)   { };

      template<typename ITER>
    xEntityIteratorBasic(const xMeshQueryInterface& _mi , const ITER &_it ) : iterator_identifier(_it), pQuery(&_mi)   { };

      xEntityIteratorBasic(const xEntityIteratorBasic& )         = default ;
      xEntityIteratorBasic( xEntityIteratorBasic& )              = default ;
      xEntityIteratorBasic( xEntityIteratorBasic&&rhs)           = default;  
      xEntityIteratorBasic& operator=(xEntityIteratorBasic&&rhs) = default;
      xEntityIteratorBasic& operator=(const xEntityIteratorBasic &rhs)  = default;
      typedef std::forward_iterator_tag iterator_category;
      typedef std::ptrdiff_t difference_type;
      typedef  xEntity value_type; 
      typedef  value_type * pointer;
      typedef  value_type  reference;
     protected:
      any iterator_identifier;
      const xMeshQueryInterface* pQuery;
    };



    class xVertexIterator : public xEntityIteratorBasic {
    public:
    xVertexIterator(): xEntityIteratorBasic() {};   
    xVertexIterator(const xEntityIteratorBasic& iter ) : xEntityIteratorBasic(iter){}  ;
      template<typename ITER>
    xVertexIterator(const xMeshQueryInterface& _mi , const ITER &_it) : xEntityIteratorBasic(_mi,_it)   {};
      xVertex operator*() const;
      typedef  xVertex value_type; 
      typedef  value_type * pointer;
      typedef  value_type  reference;
 
      xVertexIterator  operator++(int);
      xVertexIterator &operator++() ;
      bool operator!=(const xVertexIterator &in) const;
      bool operator==(const xVertexIterator &in) const;
    };


    class xEdgeIterator : public xEntityIteratorBasic {
    public:
    xEdgeIterator(): xEntityIteratorBasic() {};   
    xEdgeIterator(const xEntityIteratorBasic& iter ) : xEntityIteratorBasic(iter) {}
      template<typename ITER>
    xEdgeIterator(const xMeshQueryInterface& _mi , const ITER &_it )  : xEntityIteratorBasic(_mi,_it) {};
      xEdge operator*() const;
      typedef  xEntity value_type; 
      typedef  value_type * pointer;
      typedef  value_type  reference;
      xEdgeIterator  operator++(int) ;
      xEdgeIterator &operator++() ;
      bool operator!=(const xEdgeIterator &in) const;
      bool operator==(const xEdgeIterator &in) const;
    };



    class xFaceIterator : public xEntityIteratorBasic {
    public:
    xFaceIterator(): xEntityIteratorBasic() {}
    xFaceIterator(const xEntityIteratorBasic& iter ) : xEntityIteratorBasic(iter) {}
      template<typename ITER>
    xFaceIterator(const xMeshQueryInterface& _mi , const ITER &_it ) : xEntityIteratorBasic(_mi,_it) {};
      xFace operator*() const;
      typedef  xEntity value_type;
      typedef  value_type * pointer;
      typedef  value_type  reference;
      xFaceIterator  operator++(int);
      xFaceIterator &operator++();
      bool operator!=(const xFaceIterator &in) const;
      bool operator==(const xFaceIterator &in) const;
    };



    class xSolidIterator : public xEntityIteratorBasic {
    public:
    xSolidIterator(): xEntityIteratorBasic() {};   
    xSolidIterator(const xEntityIteratorBasic& iter ) : xEntityIteratorBasic(iter) {} ;
      template<typename ITER>
    xSolidIterator(const xMeshQueryInterface& _mi , const ITER &_it ) : xEntityIteratorBasic(_mi,_it) {};
      xSolid operator*() const ;
      typedef  xEntity value_type;
      typedef  value_type * pointer;
      typedef  value_type  reference;
      xSolidIterator  operator++(int) ;
      xSolidIterator &operator++() ;
      bool operator!=(const xSolidIterator &in) const;
      bool operator==(const xSolidIterator &in) const;
    };





    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////// general iterator with a type = xIter /////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   class xEntityIterator {
    public:
      enum iType {VERTEXITERATOR,EDGEITERATOR,FACEITERATOR,SOLIDITERATOR};
      xEntityIterator();
      xEntityIterator(const xEntityIterator& )         = default ;
      xEntityIterator( xEntityIterator& )              = default ;
      xEntityIterator( xEntityIterator&&rhs)           = default;  
      xEntityIterator& operator=(xEntityIterator&&rhs) = default;
      xEntityIterator& operator=(const xEntityIterator &rhs)  = default;
      typedef std::forward_iterator_tag iterator_category;
      typedef std::ptrdiff_t difference_type;
      typedef xEntity value_type; 
      typedef value_type * pointer;
      typedef value_type   reference;
    xEntityIterator(const xVertexIterator& _it) : iter(_it),type(iType::VERTEXITERATOR) {};
    xEntityIterator(const xEdgeIterator&   _it) : iter(_it),type(iType::EDGEITERATOR) {} ;
    xEntityIterator(const xFaceIterator&   _it) : iter(_it),type(iType::FACEITERATOR) {} ;
    xEntityIterator(const xSolidIterator&  _it) : iter(_it),type(iType::SOLIDITERATOR) {} ;

    public:
      xEntity operator*() const;
      xEntityIterator  operator++(int) ;
      xEntityIterator &operator++() ;
      bool operator!=(const xEntityIterator &in) const;
      bool operator==(const xEntityIterator &in) const;

    private:
      xEntityIteratorBasic iter;
      iType type;
    };

  } // namespace xmeshinterface
} // namespace xinterface





#endif

