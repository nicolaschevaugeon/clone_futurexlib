/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

  
#ifndef EXCHANGE_HXX
#define EXCHANGE_HXX

namespace xfem
{
// ----------------------------------------------------------------------------
// Exchange::Exchange()
// ----------------------------------------------------------------------------
template< typename ACCESS >
Exchange< ACCESS >::Exchange( tag_info t, access_type& a )
  : extag( t() ), access( a )
{
  dataSize = sizeof( value_type );
  nbytes   = sizeof( AOMD::mEntity* ) + dataSize;
}

// ----------------------------------------------------------------------------
// Exchange::~Exchange()
// ----------------------------------------------------------------------------
template< typename ACCESS >
Exchange< ACCESS >::~Exchange() 
{
}

// ----------------------------------------------------------------------------
// Exchange::AP_alloc_and_fill_buffer()
// ----------------------------------------------------------------------------
template< typename ACCESS >
void*
Exchange< ACCESS >::AP_alloc_and_fill_buffer( AOMD::mEntity* e,
                                              AOMD::AOMD_SharedInfo& si, int t1 )
{
  // -- Find if a valid references exists for that entity ---------------------
  std::vector< indexator > pvals;
  size_t nbVar = access( e, pvals ); 

  if( nbVar == 0 ) return 0;
    
  // -- Create buffer ---------------------------------------------------------
  int buffer_size = sizeof( AOMD::mEntity* ) + nbVar * dataSize;

  void* buf = reinterpret_cast<void*>( new char[ buffer_size ] );//Obsolete?


  AOMD::mEntity** ebuf = reinterpret_cast< AOMD::mEntity ** >( buf );
  *(ebuf++) = si.getRemotePointer(); // Rem: increment after dereference !!

  value_type* fbuf = reinterpret_cast< value_type* >( ebuf );

  for( typename std::vector< indexator >::const_iterator it = pvals.begin();
       it != pvals.end(); ++fbuf, ++it )
    *fbuf = access.get( *it );

  return buf;
}

// ----------------------------------------------------------------------------
// Exchange::receiveData()
// ----------------------------------------------------------------------------
template< typename ACCESS >
void Exchange< ACCESS >::receiveData( int from, void* buf )
{
    AOMD::mEntity** ebuf = reinterpret_cast< AOMD::mEntity** >( buf );
    AOMD::mEntity*  e    = *(ebuf++); // Rem: increment after dereference !!

  std::vector< indexator > pvals;
  size_t nbVar = access( e, pvals );

  if( nbVar == 0 ) return;

  value_type* dSent = reinterpret_cast< value_type* >( ebuf );
 
  for( typename std::vector< indexator >::const_iterator it = pvals.begin();
       it != pvals.end(); ++dSent, ++it )
  {
    value_type val = (*dSent);
    access.set( *it, val );
  }
}
  
// ----------------------------------------------------------------------------
// ExchangeSum::receiveData()
// ----------------------------------------------------------------------------
template< typename ACCESS >
void ExchangeSum< ACCESS >::receiveData( int from, void* buf ) 
{
  //  const bool info = true;
    AOMD::mEntity** ebuf = reinterpret_cast< AOMD::mEntity** >( buf );
    AOMD::mEntity*  e    = *(ebuf++); // Rem: increment after dereference !!

  std::vector< indexator > pvals;
  size_t nbVar = this->access( e, pvals );

  if( nbVar == 0 ) return;

  value_type* dSent = reinterpret_cast< value_type* >( ebuf );
 
  for( typename std::vector< indexator >::const_iterator it = pvals.begin();
       it != pvals.end(); ++dSent, ++it )
  {
    value_type val = this->access.get( *it );
    if (info) 
      {
	if (val != (*dSent)) 
	  {
	    std::cout << "com in exchangesum for entity:" << std::endl;
	    e->print();
	  }
      }
    val += (*dSent);

    this->access.set( *it, val );
  }
}


} // end of namespace
#endif
// == END OF FILE =============================================================

