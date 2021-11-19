#ifndef _RENUMBERING_HH
#Error Do NOT include xRenumbering_imp.h alone !
#endif




//------------------------------------------------------------------------------------
template<typename VT>
pair<size_t, size_t> xReverseCutHillMcKeeNumbering( xField<VT>& field,
                                                   xIter start, xIter last )
{
  size_t dofsize = field.getValueManager()->size( "dofs" );

  vector< set< size_t > >   graph( dofsize );
  vector< size_t >          numbering( dofsize );
  vector< size_t >          inv_numbering( dofsize );
  vector< bool >            numbered( dofsize, false );

  // -- 1. Graph construction -------------------------------------------------
  for( ; start != last; ++start )
  {
    xFiniteElement fem;
    fem.setKeys( *start, field.begin(), field.end() );

    vector< xValKey >* keys = fem.getKeys();

    for( vector< xValKey >::const_iterator itk1 = keys->begin();
         itk1 != keys->end(); ++itk1 )
    {
      xValue< double >* v1 = field.getValueManager()->find( *itk1 );
      vector< xValKey >::const_iterator itk2 = itk1;

      // if( itk2 != keys->end() ) ++itk2;

      for( ; itk2 != keys->end(); ++itk2 )
      {
        xValue< double >* v2 = field.getValueManager()->find( *itk2 );

        if( v1 != nullptr && v2 != nullptr )
        {
          xStateOfValueDof* s1 = dynamic_cast< xStateOfValueDof* >( v1->getState() );
          xStateOfValueDof* s2 = dynamic_cast< xStateOfValueDof* >( v2->getState() );

          if( s1 != nullptr && s2 != nullptr )
          {
            graph[ s1->Numdof -1 ].insert( s2->Numdof -1 );
            graph[ s2->Numdof -1 ].insert( s1->Numdof -1 );
          }

        }
      }
    }
  }

  // -- 2. Compute original bandwidth -----------------------------------------
  size_t old_bw = 0;
  std::vector<int> band(dofsize);
  for( size_t i=0; i < dofsize; ++i ){
    int iband = 0;
    std::set<size_t>::const_iterator it   = graph[i].begin();
    std::set<size_t>::const_iterator itend = graph[i].end();
    while(it!=itend){
      int dist = i - (*it);
      iband = std::max(iband, std::abs(dist));
      ++it;
    }
    band[i] = iband;
  }
  old_bw =  *max_element( band.begin(), band.end() ) ;

  /*
     std::cout<< " nb dof = " << graph.size() << "orbw " << old_bw << std::endl;
     ofstream testbw1("testbw1.txt");
     for( size_t i=0; i < dofsize; ++i )
     testbw1 << band[i] << std::endl;
  */
  /*
    testbw1 << std::endl;
  for( size_t i=0; i < dofsize; ++i ){
    std::set<size_t>::const_iterator it   = graph[i].begin();
    std::set<size_t>::const_iterator itend = graph[i].end();
    while(it!=itend){
      testbw1 << (*it) << " " ;
      ++it;
    }
    testbw1 << std::endl;
  }


  testbw1.close();
  */
  // -- 3. Cuthill-McKee-Numbering --------------------------------------------
  size_t current = dofsize / 2; // Start node
  numbering[ 0 ] = current;
  inv_numbering[ current ] = 0;
  numbered[ current ] = true;
  size_t node_progress = 0;
  size_t index_progress = 0;

  for( size_t i=1; i  < dofsize; ++i )
  {

    multimap< size_t, size_t > neighbors;
    for( set< size_t >::const_iterator j = graph[ current ].begin();
         j != graph[ current ].end(); ++j )
    {
      if( !numbered[ *j ] )
        neighbors.insert( std::make_pair( graph[*j].size(), *j ) );
    }

    while( !neighbors.empty() )
    {
      size_t candidate = neighbors.begin()->second;
      numbering[ ++index_progress ] = candidate;
      inv_numbering[ candidate ] = index_progress;
      numbered[ candidate ] = true;
      neighbors.erase( neighbors.begin() );
    }
    current = numbering[ ++node_progress ];
  }

  // -- 4. Reverse Numbering --------------------------------------------------
  vector< size_t > old_inv_numbering( inv_numbering.begin(),
                                      inv_numbering.end());

  for( size_t i = 0 ; i < dofsize ; ++i )
  {
    inv_numbering[ i ] = dofsize-1 - old_inv_numbering[ i ];
    numbering[ inv_numbering[ i ] ] = i;
  }

  // -- 5. Update graph -------------------------------------------------------
  vector< set< size_t > > newgraph( dofsize );

  for( size_t i=0; i < dofsize; ++i )
  {
    size_t oldref = numbering[ i ];
    for( set< size_t >::const_iterator j = graph[ oldref ].begin();
         j != graph[ oldref ].end(); ++j )
      newgraph[ i ].insert( inv_numbering[ *j ] );
  }

  // -- 6. Compute new bandwidth ----------------------------------------------
  size_t new_bw = 0;
  std::vector<int> new_band(dofsize);
  for( size_t i=0; i < dofsize; ++i ){
    int iband = 0;
    std::set<size_t>::const_iterator it   = graph[i].begin();
    std::set<size_t>::const_iterator itend = graph[i].end();
    while(it!=itend){
      int dist = i - (*it);
      iband = std::max(iband, std::abs(dist));
      ++it;
    }
    new_band[i] = iband;
  }
  new_bw =  *max_element( new_band.begin(),new_band.end() ) ;


  for( size_t i=1; i < dofsize; ++i )
    new_bw = std::max( new_bw,
                  std::max( i, *( max_element( newgraph[i].begin(),
                                          newgraph[i].end() ) ) ) -
                  std::min( i, *( min_element( newgraph[i].begin(),
                                          newgraph[i].end() ) ) ) );

  // -- 7. Apply Renumbering --------------------------------------------------
  typedef vector< size_t >::const_iterator const_iterator;
  std::transform( field.getValueManager()->begin( "dofs" ),  field.getValueManager()->end( "dofs" ) , field.getValueManager()->begin( "dofs" ), xNumDofRenumberingVisitor< const_iterator >( inv_numbering.begin() ) );


  /* Visit( xNumDofRenumberingVisitor< const_iterator >( inv_numbering.begin() ),
         field.getValueManager()->begin( "dofs" ),
         field.getValueManager()->end( "dofs" ) );
  */
  return std::make_pair( old_bw, new_bw );
}

// ----------------------------------------------------------------------------
// xReverseCutHillMcKeeNumberingBoost
// ----------------------------------------------------------------------------
template<typename VT>
pair<int,int> xReverseCutHillMcKeeNumberingBoost( xField<VT>& field,
                                                  xIter start, xIter last )
{
  typedef boost::adjacency_list<
    boost::setS, boost::vecS, boost::undirectedS,
    boost::property< boost::vertex_color_t, boost::default_color_type,
    boost::property< boost::vertex_degree_t, short > > > Graph;

  Graph g;

  // -- 1. Graph construction -------------------------------------------------
  for( ; start != last; ++start )
  {
    xFiniteElement fem;
    fem.setKeys( *start, field.begin(), field.end() );

    vector< xValKey >* keys = fem.getKeys();

    for( vector< xValKey >::const_iterator itk1 = keys->begin();
         itk1 != keys->end(); ++itk1 )
    {
      xValue< double >* v1 = field.getValueManager()->find( *itk1 );
      vector< xValKey >::const_iterator itk2 = itk1;

      //for( ++itk2; itk2 != keys->end(); ++itk2 ) (self in the graph or not ??)
      for( ; itk2 != keys->end(); ++itk2 )
      {
        xValue< double >* v2 = field.getValueManager()->find( *itk2 );

        if( v1 != nullptr && v2 != nullptr )
        {
          xStateOfValueDof* s1 = dynamic_cast< xStateOfValueDof* >( v1->getState() );
          xStateOfValueDof* s2 = dynamic_cast< xStateOfValueDof* >( v2->getState() );

          if( s1 != nullptr && s2 != nullptr ){
            boost::add_edge( s1->Numdof - 1, s2->Numdof - 1, g );
      }
        }
      }
    }
  }

  int old_bw = boost::bandwidth( g ); //there is a bug in boost::bandwidth
  //std::vector<int> band(num_vertices(g));
  //for (int i=0; i< num_vertices(g); ++i)
  //  band[i] = ith_bandwidth(i, g);
  //ofstream bandboost("bandboost.txt");
  //for( size_t i=0; i < num_vertices(g); ++i )
  //  bandboost << band[i] << std::endl;
  //bandboost << std::endl;
  //g::vertex_iterator it, itend;

  /*
  Graph::vertex_iterator it;
  Graph::vertex_iterator itend;
  boost::tie(it, itend) =  vertices(g);
  while(it!=itend){
    boost::graph_traits<Graph>::adjacency_iterator itad;
    boost::graph_traits<Graph>::adjacency_iterator itadend;
    boost::tie(itad, itadend) =   adjacent_vertices(*it,g);
    while (itad !=itadend){
      bandboost << (*itad) << " -  " ;
      ++itad;
    }
    bandboost << std::endl;
    ++it;
  }
  */




  //ndboost.close();




  // -- 2. New order computation ----------------------------------------------
  typedef boost::graph_traits< Graph >::vertex_descriptor    gVertex;
  typedef boost::graph_traits< Graph >::vertices_size_type   g_size_type;

  boost::property_map< Graph, boost::vertex_index_t >::type index_map =
    boost::get( boost::vertex_index, g );

  vector< gVertex >     inv_perm( boost::num_vertices( g ) );
  vector< g_size_type > perm( boost::num_vertices( g ) );

  boost::cuthill_mckee_ordering( g, inv_perm.rbegin(),
                                 boost::get( boost::vertex_color, g ),
                                 boost::make_degree_map( g ) );

  for( g_size_type c = 0; c != inv_perm.size(); ++c )
    perm[ index_map[ inv_perm[c] ] ] = c;

  int new_bw = boost::bandwidth( g,
                                 boost::make_iterator_property_map( &perm[0],
                                                                    index_map,
                                                                    perm[0] ));

  // -- 3. Renumbering --------------------------------------------------------
  typedef vector< g_size_type >::const_iterator const_index_iterator;
 std::transform( field.getValueManager()->begin( "dofs" ),  field.getValueManager()->end( "dofs" ) , field.getValueManager()->begin( "dofs" ), xNumDofRenumberingVisitor< const_index_iterator >( perm.begin() ) );
 /*
  Visit( xNumDofRenumberingVisitor< const_index_iterator >( perm.begin() ),
         field.getValueManager()->begin( "dofs" ),
         field.getValueManager()->end( "dofs" ) );
 */


  return std::make_pair( old_bw, new_bw );
}
//------------------------------------------------------------------------------------

