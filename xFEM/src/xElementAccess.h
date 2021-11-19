/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

  
#ifndef ELEMENT_ACCESSOR_HH
#define ELEMENT_ACCESSOR_HH

// ----------------------------------------------------------------------------
// HEADERS
// ----------------------------------------------------------------------------

// -- AOMD --------------------------------------------------------------------
#include "mEntity.h"

// -- XFEM --------------------------------------------------------------------
#include "xField.h"
#include "xValueManager.h"
#include "xStateOfValue.h"

//mat
#include "ParUtil.h"
#include "xTag.h"


namespace xfem
{

// ----------------------------------------------------------------------------
// CLASS ElementAccessor
// ----------------------------------------------------------------------------
/*! \ingroup Parallel
    \brief Access class for data exchange between partitions.

    This class provides access to vector entries and is designed to work with
    the Exchange or ExchangeSum classes. The vector entries are assumed to
    be based on the numbering of degree of freedom of a given field (which
    is used to reference value).

    It is also assumed that the VECTOR type is an indexed container and that
    this container provides the value_type type definition and access
    operation by means of brackets [ ].

    
*/
template< typename VECTOR, typename VT >
class VectorElementByEntity
{
  public:
    // -- Traits --------------------------------------------------------------
    typedef VECTOR                                              storage_type;
    typedef size_t                                              indexator;
    typedef typename storage_type::value_type                   value_type;
    typedef Tag< int >                                          tag_info;
    
    // -- Constructors & destructor -------------------------------------------
    //! Constructor requires the reference field and the storage structure.
    VectorElementByEntity( xField<VT>& f,
                           storage_type& v ) : field( f ), data( v ) {}
    
    ~VectorElementByEntity() {}

    // -- Functionality -------------------------------------------------------
    //! Fills pvals with the indexator related to e and returns pvals size.
    size_t operator()( AOMD::mEntity* e, std::vector< indexator >& pvals ) const;

    //! Returns value referenced by indexator i.
    value_type get( const indexator& i ) const { return data[i]; }

    //! Assigns value referenced by indexator i with received value val.
    void       set( const indexator& i, value_type& val ) { data[i] = val; }

  private:
    xField<VT>&       field;
    storage_type&  data;
};
    
    
    
// ----------------------------------------------------------------------------
// CLASS MatrixBlocByEntity
// ----------------------------------------------------------------------------
/*! \ingroup Parallel
    \brief Access class for data exchange between partitions.
    
    This class provides access to matrix entries and is designed to work with
    the Exchange or ExchangeSum classes. The matrix entries are assumed to
    be based on the numbering of degree of freedom of a given field (which
    is used to reference value).

    It is also assumed that the MATRIX type is an indexed container and that
    this container provides the value_type type definition and access
    operation by means of brackets [ ].

    Rem: here, indexator is a pair of integer (row and column indexes).
*/
template< typename MATRIX, typename VT >
class MatrixBlocByEntity
{
  public:
    // -- Traits --------------------------------------------------------------
    typedef MATRIX                                                matrix_type;
    typedef std::pair< size_t, size_t >                           indexator;
    typedef typename matrix_type::value_type                      value_type;
    typedef Tag< int >                                            tag_info;

    // -- Constructors & destructor -------------------------------------------
    //! Constructor requires the reference field and the storage structure.
    MatrixBlocByEntity( xField<VT>& f, matrix_type& m )
      : field( f ), data_send( m ), data_recv( m ) {}

    //! Constructor requires the reference field and storage structures
    //! (target differs from source).
    MatrixBlocByEntity( xField<VT>& f, matrix_type& source, matrix_type& target )
      : field( f ), data_send( source ), data_recv( target ) {}

    ~MatrixBlocByEntity(){}
    
    // -- Functionality -------------------------------------------------------
    //! Fills pvals with the indexator related to e and return pvals size.
    size_t operator()( AOMD::mEntity* e, std::vector< indexator >& pvals ) const;

    //! Returns value referenced by indexator i.
    value_type get( const indexator& i ) const 
      { return data_send[ i.first][ i.second ]; }

    //! Assigns value referenced by indexator i with received value val.
    void       set( const indexator& i, value_type val ) 
      { data_recv[ i.first ][ i.second ] = val; } 

  private:

    xField<VT>&     field;
    matrix_type& data_send;
    matrix_type& data_recv;
};


//ajout mathieu 
// ----------------------------------------------------------------------------
// CLASS  VectorBoolByEntity 
// ----------------------------------------------------------------------------

//template< typename BOOL_TYPE>
template< typename X >
class AccessByEntity
{
  public:
    // -- Traits --------------------------------------------------------------
    //typedef BOOL_TYPE                                         storage_type;
    typedef std::map<AOMD::mEntity*,X>                                storage_type;
    
    //typedef size_t                                              indexator;
    typedef typename storage_type::iterator                       indexator;
    //typedef typename storage_type::value_type                   value_type;
    typedef X                                                     value_type;
    typedef Tag< int >                                          tag_info;
 
    // -- Constructors & destructor -------------------------------------------
    AccessByEntity(storage_type& v ) : data( v ) {}
    
    ~AccessByEntity() {}

    // -- Functionality -------------------------------------------------------
    size_t operator()( AOMD::mEntity* e, std::vector< indexator >& pvals ) const;

    //! Returns value referenced by indexator i.
    value_type get( const indexator& i ) const { return (*i).second; }

    //! Assigns value referenced by indexator i with received value val.
    void       set( const indexator& i, value_type& val ) { (*i).second = val; }

  private:
    storage_type&  data;
};




} // end of namespace

// ----------------------------------------------------------------------------
// IMPLEMENTATION
// ----------------------------------------------------------------------------
#include "xElementAccess_imp.h"


#endif
// == END OF FILE =============================================================


