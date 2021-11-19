/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

  
#ifndef EXCHANGE_HH
#define EXCHANGE_HH

// ----------------------------------------------------------------------------
// HEADERS
// ----------------------------------------------------------------------------

// -- AOMD --------------------------------------------------------------------
//#include "pmExchangeData.h"
#include "mEntity.h"
#include "mExchangeData.h"


// -- XFEM --------------------------------------------------------------------
#include "xStateOfValue.h"
#include "xFiniteElement.h"
// ----------------------------------------------------------------------------
// ENUM Tag
// ----------------------------------------------------------------------------
#include "xTag.h"


namespace xfem
{
// ----------------------------------------------------------------------------
// CLASS Exchange
// ----------------------------------------------------------------------------
// template type ACCESS requirements :
//  - access provides value_type
//  - provides  value_type* operator()( AOMD::mEntity* e )
//
/*! \ingroup Parallel
    \brief Exchange method for mesh based indexed container.

    The Exchange class provides a concrete realization of the abstract
    AOMD_DataExchanger class. This class is templatized on an ACCESS class,
    providing access to the data to be exchanged. This class is used for
    a simple data exchange.

    This access class must have the following internal types defined:
    - indexator : a reference type to data storage
    - value_type: the type of data associated with one indexator

    In addition to these requirements:
    - the indexator class must support the const and non-const dereferencing
      operator
    - the access class must have the public method <BR>
      <CENTER>
      size_t operator()( AOMD::mEntity*, std::vector< indexator >& ) 
      </CENTER>
      filling the vector with all the indexators that the given AOMD::mEntity object
      involves.
    - the access class must have the public method <BR>
      <CENTER>
      value_type get( const indexator& ) const
      </CENTER>
      returning the value associated with the given indexator.
    - the access class must have the public method <BR>
      <CENTER>
      void set( const indexator&, value_type& ) 
      </CENTER>
      assigning the given value to the value associated with the given
      indexator.
    
    Examples of access classes: MatrixBlocByEntity, VectorElementByEntity
*/
template< typename ACCESS >
class Exchange : public AOMD::AOMD_DataExchanger
{
  public:
    // -- Traits --------------------------------------------------------------
    typedef ACCESS                                                 access_type;
    typedef typename access_type::tag_info                         tag_info;
    typedef typename access_type::indexator                        indexator;
    typedef typename access_type::value_type                       value_type;

    // -- Constructors & destructor -------------------------------------------
    //! Constructor requires a tag and the access object.
    Exchange( tag_info t,  access_type& a );
    virtual ~Exchange();
    
    // -- Abstract methods ----------------------------------------------------
    //! Sending data to other partitions
     void* AP_alloc_and_fill_buffer( AOMD::mEntity* e, AOMD::AOMD_SharedInfo& si, int t1 );
    //! Receiving data from other partitions
    virtual void receiveData( int from, void* buf );
    
    // -- Infos methods -------------------------------------------------------
    //! returns the tag used
    int                tag() const { return extag; }
        
  protected:
    int              extag;  
    access_type&     access;  

    int              nbytes;
    int              dataSize;
};

// ----------------------------------------------------------------------------
// CLASS ExchangeSum
// ----------------------------------------------------------------------------
/*! \ingroup Parallel
    \brief Additive exchange method for mesh based indexed container.

    Refinement of Exchange class. This class simply redefines the receiveData()
    method in order to sum the received data to the local values.

    See Exchange for more informations
*/
template< typename ACCESS >
class ExchangeSum : public Exchange< ACCESS >
{
  public:

    typedef ACCESS                                                 access_type;
    typedef typename access_type::tag_info                         tag_info;
    typedef typename access_type::indexator                        indexator;
    typedef typename access_type::value_type                       value_type;

    // -- Constructors & destructor -------------------------------------------
    //! Constructor requires a tag and the access object.
    ExchangeSum( tag_info t,  access_type& a ) : 
      Exchange< ACCESS >( t, a ), info(false) {}
    virtual ~ExchangeSum() {}
    
    // -- Abstract methods ----------------------------------------------------
    //! Addingdata from other partition.
    void receiveData( int from, void* buf );
    void setDebugInfo(bool in){info=in;};
  private:
    bool info;
};

// ----------------------------------------------------------------------------
// CLASS ExchangeMin
// ----------------------------------------------------------------------------
// NOT NEEDED NOW

// ----------------------------------------------------------------------------
// CLASS ExchangeMax
// ----------------------------------------------------------------------------
// NOT NEEDED NOW



template <class FIELD>
class xExchangeStateOfValue : public AOMD::AOMD_DataExchanger
{
public:
  xExchangeStateOfValue(const FIELD &_fct ):fct(_fct),extag(33){};
  //~xExchangeStateOfValue(); 
  //! Sending data to other partitions
  void* AP_alloc_and_fill_buffer( AOMD::mEntity* e, AOMD::AOMD_SharedInfo& si, int t1 ){
    typename FIELD::ValueManager* val_manager = fct.getValueManager();
    xFiniteElement FEM;
    FEM.setKeys(e, fct.begin(), fct.end());
    std::vector<xValue<double>*> vals;
    val_manager->getValPtr(FEM.beginKey(), FEM.endKey(), vals);
    int n = vals.size();
    std::vector<int> fixed(n);
    std::vector<double> fixedv(n);
    for (int i=0; i< n; ++i){
      xStateOfValueFixed  *isfixed=  dynamic_cast <xStateOfValueFixed  *> (vals[i]->getState());
      fixed[i] = isfixed?1:0;
      fixedv[i] = isfixed?vals[i]->getVal():0.;
    }
    int dataSize= sizeof(int);
    int buffer_size = sizeof( AOMD::mEntity* ) + n * sizeof(int) + n*sizeof(double);
    void * buf = 0;

    
    //void* buf = this->AP_alloc( si.pid(), extag, buffer_size );
    AOMD::mEntity** ebuf = reinterpret_cast< AOMD::mEntity ** >( buf );
    *(ebuf++) = si.getRemotePointer();
    int* ibuf = reinterpret_cast< int* >( ebuf );
    for (int i=0; i<n; ++i, ++ibuf )
      *ibuf = fixed[i];
    
    double* dbuf = reinterpret_cast< double* >( ibuf );
    for (int i=0; i<n; ++i, ++dbuf )
      *dbuf = fixedv[i];
    
    return buf;
  };
  //! Receiving data from other partitions
  void receiveData( int from, void* buf ){
      AOMD::mEntity** ebuf = reinterpret_cast< AOMD::mEntity** >( buf );
      AOMD::mEntity*  e    = *(ebuf++); 
    typename FIELD::ValueManager* val_manager = fct.getValueManager();
    xFiniteElement FEM;
    FEM.setKeys(e, fct.begin(), fct.end());
    std::vector<xValue<double>*> vals;
    val_manager->getValPtr(FEM.beginKey(), FEM.endKey(), vals);
    int n = vals.size();
    
    std::vector<int> fixed(n);
    std::vector<double> fixedv(n);
    if( n == 0 ) return;
    
    int *ibuf = reinterpret_cast< int * >( ebuf );
    
    for( int i =0; i <n; ++i, ++ibuf)
      {
	fixed[i] = *ibuf;
      }
    
    double *dbuf = reinterpret_cast< double * >( ibuf );
    
    for( int i =0; i <n; ++i, ++dbuf)
      {
	fixedv[i] = *dbuf;
      }
    
    for (int i=0; i<n; ++i)
      if (fixed[i]) {
	vals[i]->setState(new xStateOfValueFixed(vals[i]));
	vals[i]->setVal(fixedv[i]);
      }
    
    
  };
  int  tag() const { return extag; }
private:
  const FIELD & fct;
  const int extag;
};
}
// end of namespace
// ----------------------------------------------------------------------------
// IMPLEMENTATION
// ----------------------------------------------------------------------------

#include "xExchange_imp.h"


#endif
// == END OF FILE =============================================================
