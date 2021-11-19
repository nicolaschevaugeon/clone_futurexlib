
#include "xPartitionManagerXentity.h"

////////////////////////////////////////////////////// 
// ---------- xPartitionManagerXentity ---------------
//////////////////////////////////////////////////////

namespace xinterface{

  namespace xmeshinterface{

    using std::cout;
    using std::cerr;
    using std::endl;

    void xPartitionManagerXentity::print() const {
      cout<<"[print] xPartitionManagerXentity::print() " <<query.getName()<<endl ;
      int count=0;
      for (auto it = beginObject() ; it != endObject(); ++it ) {
	xEntity enti = *it;
	count++;
	cout<<"  - local entity #"<<count<<endl;
	enti.print();
	auto po   = getConstPartitionObject(enti);
	for (auto ro : po.getRemoteObjectsCollectionRange( ) ){
	  cout<<"           -> remote object "<<ro.getObjectAddress() <<" on  proc " << ro.getProcId()  <<endl ;
	}
      }
    };


    xEntity xPartitionManagerXentity::getXentityFromUniqueAddress( const void* address  ) const{
      return query.getXentityFromUniqueAddress( address );
    };



    xEntity xPartitionManagerXentity::getXentityFromUniqueAddressAndQueryTag( const void* address , unsigned int tag ) const{
      if (tag != query.getTag()) 
	{
	  try  {     
	    throw std::logic_error( "ERROR in  xPartitionManagerXentity::getXentityFromUniqueAddressAndQueryTag(.,.): the specified tag is not the tag of the current query interface. Use xPartitionManagerXentity::getXentityFromUniqueAddress(.) instead, or verify why it is not the good one." ); 
	  } 
	  catch ( std::exception & e ) { cerr << e.what()<<endl; exit(0) ;} 
	}

      return  getXentityFromUniqueAddress(  address  );
    };



 

    ////////////////////////////////////////////////////////// 
    // ---------- xPartitionManagerXentityUnion --------------
    ////////////////////////////////////////////////////////// 


    xPartitionManagerXentityUnion::~xPartitionManagerXentityUnion() {
      //    partman_container.clear();
    }


    void xPartitionManagerXentityUnion::add(const xPartitionManagerXentity& new_element) { 
      const xPartitionManagerXentity* tmp=&new_element;
      partman_container.push_back(tmp);
    }

 
    const xtool::xConstPartitionObject < xEntity > xPartitionManagerXentityUnion::getConstPartitionObject( xEntity& enti ) const{
      // scan all partition of the union
      for (auto it = partman_container.begin(); it != partman_container.end(); ++it) {
	const xPartitionManagerXentity* partman=*it;
	cout<<"xPartitionManagerXentityUnion::getConstPartitionObject: query name: "<<partman->getQuery().getName()<<",   tag: "<<partman->getQuery().getTag()<<endl;
	//  find if the enti is in the partman
	if( partman->getConstPartitionObject( enti ).hasRemoteObject())   return partman->getConstPartitionObject(enti);
      } 
      // if here enti has no remote object, but may have a empty PartitionObject. Any partman of the partman_container can return a empty PatitionObject.
      return (*(partman_container.begin()))->getConstPartitionObject( enti );
    };

    /*    xEntity xPartitionManagerXentityUnion::getXentityFromUniqueAddress( const void* address  ) const{
      try  {  throw logic_error( "ERROR in  xPartitionManagerXentityUnion::getXentityFromUniqueAddress(.): this function is unavailable for xPartitionManagerXentityUnion because it has no unique query. Use xPartitionManagerXentityUnion::getXentityFromUniqueAddressAndQueryTag(.,.) instead and specify the tag of the query you want to use." );       } 
      catch ( exception & e ) { cerr << e.what()<<endl; exit(0) ;} 
    };
    */

    xEntity xPartitionManagerXentityUnion::getXentityFromUniqueAddressAndQueryTag( const  void* address , unsigned int tag) const {
      // scan all partition of the union
      for (auto it = partman_container.begin(); it != partman_container.end(); ++it) {
	const xPartitionManagerXentity* partman=*it;
	//    cout<<"partman query "<<  partman->getQuery()->getName() <<endl;
	//    cout<<"partman query "<<  partman->getQuery()->getTag() <<endl;
	//    cout<<" to compare with "<<tag<<endl;
	if ( partman->getQuery().getTag() == tag ) {
	  return partman->xPartitionManagerXentity::getXentityFromUniqueAddress(address);
	}
      } 
      // if here, the enti does not belong to any partman : should not be possible
      cout<<"Warning in xPartitionManagerXentityUnion::getXentityFromUniqueAddressAndQueryTag(.,.):  the tag specified is: '"<<tag<<"'"<<endl;
      try  { throw std::logic_error( "ERROR: this tag does not belong to any partition manager. The entity cannot be found." );  } 
      catch ( std::exception & e ) { cerr << e.what()<<endl; exit(0) ;}    
    };




    void xPartitionManagerXentityUnion::print() const {
      cout<<"print xPartitionManagerXentityUnion:" <<endl ;
      for (auto it = partman_container.begin(); it != partman_container.end(); ++it) {
	const xPartitionManagerXentity* partman=*it;
	partman->print();
      }
    };

    /*
    const xMeshQueryInterface* xPartitionManagerXentityUnion::getQuery() const  { 
      try  { throw std::logic_error( "ERROR in xPartitionManagerXentityUnion::getQuery(): never use getQuery() function for xPartitionManagerXentityUnion. the latter contains several xPartitionManagerXentity with one query for each xPartitionManagerXentity. That means that xPartitionManagerXentityUnion have not a unique query interface. This function is only available for xPartitionManagerXentity. " );  } 
      catch ( std::exception & e ) { cerr << e.what()<<endl; exit(0) ;}    
    };
    */

  } // namepsace xmeshinterface
} // namespace xinterface
