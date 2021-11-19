/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef ELEMENT_ACCESSOR_HXX
#define ELEMENT_ACCESSOR_HXX


namespace xfem
{

// ----------------------------------------------------------------------------
// VectorElementByEntity< STORAGE >::operator()( AOMD::mEntity* e ) const
// ----------------------------------------------------------------------------
template< typename VECTOR >
size_t
VectorElementByEntity< VECTOR >::operator()( AOMD::mEntity* e,
                                             std::vector<indexator>& pvals ) 
const
{ 
  pvals.clear();
  xFiniteElement fem;
  fem.setKeys( e, field.begin(), field.end(), "shape" );
  vector< xValKey >* keys = fem.getKeys( "shape" );
  for( std::vector< xValKey >::const_iterator itk = keys->begin();
       itk != keys->end(); ++itk ) 
  {
    if( itk->getEnti() == e ) 
    { // Avoid key for sub component entitites
      xValue<double>* p = field.getValueManager()->find( *itk );
      xStateOfValueDof* s = dynamic_cast< xStateOfValueDof* >( p->getState() );
      if( s != 0 )
	{
	  pvals.push_back( s->Numdof-1 );
	  // STORAGE is assumed 0-based indexed storage
	  // TOCHANGE : introduce an indextag
	  // coming from INDEXED_STORAGE, or index_tag      
	} 
    }
  }
  
  return pvals.size();
}

// ----------------------------------------------------------------------------
// MatrixBlocByEntity< MATRIX >::operator()
// ----------------------------------------------------------------------------
template< typename MATRIX >
size_t 
MatrixBlocByEntity< MATRIX >::operator()( AOMD::mEntity* e, 
                                          std::vector< indexator >& pvals )
const
{
  pvals.clear();
  vector< bool > mask;

  xFiniteElement fem;  
  // Simplification : assumes that the same interpolation is used
  // for test and trial 
  fem.setKeys( e, field.begin(), field.end(), "shape" ); 
  vector< xValKey >* keys = fem.getKeys( "shape" );
  size_t e_size = keys->size();
  mask.resize( e_size * e_size, true );
  int dim = e->getLevel();

  // -- Set the matrix mask --------------------------------------------------
  if( dim != 0 ) 
  {
    for( int k=0; k < e->size( dim-1 ); ++k )
    {
        AOMD::mEntity* sub_e = e->get( dim-1, k );
      xFiniteElement sub_fem;
      sub_fem.setKeys( sub_e, field.begin(), field.end(), "subshape" );
      vector< xValKey >* sub_keys = sub_fem.getKeys( "subshape" );
      for( vector< xValKey >::iterator i = sub_keys->begin();
           i != sub_keys->end(); ++i )
      {
        for( vector< xValKey >::iterator j = i; j != sub_keys->end(); ++j )
        {
          vector< xValKey >::iterator k1 = find( keys->begin(), keys->end(),
                                                  *i );
          vector< xValKey >::iterator k2 = find( keys->begin(), keys->end(),
                                                  *j );
          if( k1 != keys->end() && k2 != keys->end() ) 
          {
            size_t idx1 = distance( keys->begin(), k1 );
            size_t idx2 = distance( keys->begin(), k2 );
            mask[ idx1 * e_size + idx2 ] = false;
            mask[ idx2 * e_size + idx1 ] = false;
          }
          // else ... nothing                                  
        }
      }

    }
  }

  // -- Loop over matrix entries ----------------------------------------------
  for( size_t i=0; i < e_size; ++i )
  {
    for( size_t j=0; j < e_size; ++j )
    {
      if( mask[ i * e_size + j ] )
      {
        xValue< double >* p = 0;
        p = field.getValueManager()->find( (*keys)[i] );
        xStateOfValueDof* state1 = dynamic_cast< xStateOfValueDof* >( p->getState() );
        p = field.getValueManager()->find( (*keys)[j] );
        xStateOfValueDof* state2 = dynamic_cast< xStateOfValueDof* >( p->getState() );

        if( state1 != 0 && state2 != 0 )
        {
          pvals.push_back( make_pair( state1->Numdof-1, state2->Numdof-1 ) );
          // STORAGE is assumed 0-based indexed storage
          // TOCHANGE : introduce an indextag
          // coming from INDEXED_STORAGE, or index_tag 
        }
      }
    }
  }

  return pvals.size();
}



//ajout mathieu
// ----------------------------------------------------------------------------
// AccessByEntity< STORAGE >::operator()( AOMD::mEntity* e ) const
// ----------------------------------------------------------------------------
template< typename X >
size_t AccessByEntity< X >::operator()( AOMD::mEntity* e, std::vector<indexator>& pvals ) const
{
  pvals.clear();
  pvals.push_back(data.find(e));
  return pvals.size();
}


} // end of namespace
#endif
// == END OF FILE =============================================================


