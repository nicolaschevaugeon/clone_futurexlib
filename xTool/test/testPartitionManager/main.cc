#include "xPartitionManager.h"
#include "xUnorderedMapDataManager.h"
#include "mpi.h"
#include <iostream>
#include <vector>


class entity {
public: 
  entity( ):n{'a'}{}
  entity(  char _n):n{_n}{}
  char n;
};

//! A small test of xPartition Manager on 4 mpi procs.
// a group of 6 entity  a, b, c, d, e, f are partitionned between 4 proc :
// proc 0 has a, b, c, d
// proc 1 has b, d, f
// proc 2 has c, d, e, f
// proc 3 has a, c, e.
// each proc r create each entity x it has and then send to the procs which also have x the address of x on r. 
// upon receive each proc update their partition manager.
// owned are counted by each processor, then checkPartition is called.
// Then entity c is removed from proc 2
// checkPartition is called -> we expect that it return false for proc 1 and 2
// since the partition manager is not udated yet.
// Then the partition manager is updated (testing xPartitionManager::remove )
// check partition is then called on a copy of the initial partition manager ( the intent is to
//  test the partition manager copy constructor)


//xPartitionManager function tested
//  directly:
//    constructor
//    destructor
//    insert
//    range
//    remove( object, rank)
//    remove( object )
//    copy constructor
//    clear
//    getOwner !!!
//    beginObject,endObject
//  from  countOwned :
//    isOwner
//  from checkPartition :
//    getComm
//    getRemoteObjectOn
// xDataManager is tested also throught the calls on the partition manager :
//   tested :
//   constructor, destructor, copy constructor
//   getData( const KEY & ) const
//   getData( KEY & ) 
//   setData( KEY &)
//   deleteDat(KEY &)
//   clear()
//   beginKey,endKey
//   constIterKey class


 template <class T>
  using datman =  xtool::xUnorderedMapDataManager< entity, T > ;

typedef xtool::xPartitionManager< datman > partman_t;

int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank);
  MPI_Comm_size( MPI_COMM_WORLD, &size);
  std::string fo = "proc_"+std::to_string(rank) +"_output.txt";
  auto ok = freopen (fo.c_str(),"w",stdout);
  if(!ok){std::cout << "Can't reopen stdout on file "<< fo << __FILE__<< __LINE__ << std::endl; throw;}
  
  if (size != 4){
    if (!rank) std::cout << argv[0] << " Need at 4 mpi process to run "<< std::endl ; 
    MPI_Finalize();  
    return 0;
  }
  
  
  partman_t partman(MPI_COMM_WORLD);
  std::vector <entity *> entities;

  if (rank == 0){
    entities.push_back(new entity('a'));
    entities.push_back(new entity('b'));
    entities.push_back(new entity('c'));
    entities.push_back(new entity('d'));
  }
  if (rank == 1){
    entities.push_back(new entity('b'));
    entities.push_back(new entity('d'));
    entities.push_back(new entity('f'));
  }
  if (rank == 2){
    entities.push_back(new entity('c'));
    entities.push_back(new entity('d'));
    entities.push_back(new entity('e'));
    entities.push_back(new entity('f'));
  }
  if (rank == 3){
    entities.push_back(new entity('a'));
    entities.push_back(new entity('c'));
    entities.push_back(new entity('e'));
  }
 
  
  MPI_Status status;
  if (rank == 0){
    MPI_Send(&entities[0], sizeof(entity *) , MPI_BYTE, 3, 1 ,MPI_COMM_WORLD);
    MPI_Send(&entities[1], sizeof(entity *), MPI_BYTE, 1, 1 ,MPI_COMM_WORLD); 
    MPI_Send(&entities[2], sizeof(entity *), MPI_BYTE, 2, 1 ,MPI_COMM_WORLD); 
    MPI_Send(&entities[2], sizeof(entity *), MPI_BYTE, 3, 2 ,MPI_COMM_WORLD);
    MPI_Send(&entities[3], sizeof(entity *), MPI_BYTE, 1, 2 ,MPI_COMM_WORLD);
    MPI_Send(&entities[3], sizeof(entity *), MPI_BYTE, 2, 2 ,MPI_COMM_WORLD); 
    auto partentity0 = partman.getPartitionObject(*entities[0]);
    auto partentity1 = partman.getPartitionObject(*entities[1]);
    auto partentity2 = partman.getPartitionObject(*entities[2]);
    auto partentity3 = partman.getPartitionObject(*entities[3]);
   
    entity *recv;
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 1, 1 ,MPI_COMM_WORLD, &status);
    partentity1.insert( 1, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 1, 2 ,MPI_COMM_WORLD, &status);  
    partentity3.insert( 1, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 2, 1 ,MPI_COMM_WORLD, &status);
    partentity2.insert( 2, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 2, 2 ,MPI_COMM_WORLD, &status);
    partentity3.insert( 2, recv); 
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 3, 1 ,MPI_COMM_WORLD, &status);
    partentity0.insert( 3, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 3, 2 ,MPI_COMM_WORLD, &status);
    partentity2.insert(3, recv);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "PROC" << rank << " :" << std::endl;
    for(auto e : entities ){
      std::cout << " has " << e->n << " at " << e;
      auto partobj = partman.getConstPartitionObject(*e);
      auto rmotes = partobj.getRemoteObjectsCollectionRange( );
      std::cout << " remote copies : " <<  rmotes.size() << " : ";
      for (auto r : rmotes ){
	std::cout << " " << r.getProcId() << " " << r.getObjectAddress();
      }
      auto owner = partobj.getOwner ( );
      std::cout << " Owner is " <<  owner.getProcId()  << " " <<  owner.getObjectAddress() << std::endl;
      std::cout << " Am I Owner ? "  << partobj.isOwner () << std::endl;
    }
    std::cout << std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (rank == 1){
    auto partentity0= partman.getPartitionObject(*entities[0]);
    auto partentity1= partman.getPartitionObject(*entities[1]);
    auto partentity2= partman.getPartitionObject(*entities[2]);
    
    entity *recv;
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 0, 1 ,MPI_COMM_WORLD, &status);
    partentity0.insert( 0, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 0, 2 ,MPI_COMM_WORLD, &status);
    partentity1.insert( 0, recv);
    
    MPI_Send(&entities[0], sizeof(entity *), MPI_BYTE, 0, 1 ,MPI_COMM_WORLD); 
    MPI_Send(&entities[1], sizeof(entity *), MPI_BYTE, 0, 2 ,MPI_COMM_WORLD); 
    MPI_Send(&entities[1], sizeof(entity *), MPI_BYTE, 2, 1 ,MPI_COMM_WORLD); 
    MPI_Send(&entities[2], sizeof(entity *), MPI_BYTE, 2, 2 ,MPI_COMM_WORLD);

    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 2, 1 ,MPI_COMM_WORLD, &status);
    partentity1.insert(2, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 2, 2 ,MPI_COMM_WORLD, &status);
    partentity2.insert( 2, recv);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "PROC" << rank << " :" << std::endl;
    for(auto e : entities ){
      std::cout << " has " << e->n << " at " << e;
      auto partobj = partman.getConstPartitionObject(*e);
      auto rmotes = partobj.getRemoteObjectsCollectionRange( );
      std::cout << " remote copies : " <<  rmotes.size() << " : ";
      for (auto r : rmotes){
	std::cout << " " << r.getProcId() << " " << r.getObjectAddress();
      }
      auto owner =  partobj.getOwner ( );
      std::cout << " Owner is " <<  owner.getProcId()  << " " <<  owner.getObjectAddress() << std::endl;
      std::cout << " Am I Owner ? "  << partobj.isOwner () << std::endl;
    }
    std::cout << std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
  } 
  if (rank == 2){
    auto partentity0= partman.getPartitionObject(*entities[0]);
    auto partentity1= partman.getPartitionObject(*entities[1]);
    auto partentity2= partman.getPartitionObject(*entities[2]);
    auto partentity3= partman.getPartitionObject(*entities[3]);

    entity *recv;
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 0, 1 ,MPI_COMM_WORLD, &status);
    partentity0.insert( 0, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 0, 2 ,MPI_COMM_WORLD, &status);
    partentity1.insert( 0, recv);
   
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 1, 1 ,MPI_COMM_WORLD, &status);
    partentity1.insert( 1, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 1, 2 ,MPI_COMM_WORLD, &status);
    partentity3.insert( 1, recv);
    
    MPI_Send(&entities[0], sizeof(entity *), MPI_BYTE, 0, 1 ,MPI_COMM_WORLD); 
    MPI_Send(&entities[0], sizeof(entity *), MPI_BYTE, 3, 1 ,MPI_COMM_WORLD); 
    MPI_Send(&entities[1], sizeof(entity *), MPI_BYTE, 0, 2 ,MPI_COMM_WORLD);
    MPI_Send(&entities[1], sizeof(entity *), MPI_BYTE, 1, 1 ,MPI_COMM_WORLD);
    MPI_Send(&entities[2], sizeof(entity *), MPI_BYTE, 3, 2 ,MPI_COMM_WORLD);
    MPI_Send(&entities[3], sizeof(entity *), MPI_BYTE, 1, 2 ,MPI_COMM_WORLD);
    
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 3, 1 ,MPI_COMM_WORLD, &status);
    partentity0.insert( 3, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 3, 2 ,MPI_COMM_WORLD, &status);
    partentity2.insert( 3, recv);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "PROC" << rank << " :" << std::endl;
    
    for(auto e : entities ){
      std::cout << " has " << e->n << " at " << e;
      auto partobj = partman.getConstPartitionObject(*e);
      auto rmotes = partobj.getRemoteObjectsCollectionRange( );
      std::cout << " remote copies : " <<  rmotes.size() << " : ";
      for (auto r : rmotes){
	std::cout << " " << r.getProcId() << " " << r.getObjectAddress();
      } 
      auto owner = partobj.getOwner ();
      std::cout << " Owner is " <<  owner.getProcId()  << " " <<  owner.getObjectAddress() << std::endl;
      std::cout << " Am I Owner ? "  << partobj.isOwner () << std::endl;
    }
    std::cout << std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
    
  }
  if (rank == 3){
    auto partentity0 = partman.getPartitionObject(*entities[0]);
    auto partentity1 = partman.getPartitionObject(*entities[1]);
    auto partentity2 = partman.getPartitionObject(*entities[2]);
    
    entity *recv;
    MPI_Recv(&recv, sizeof( entity *), MPI_BYTE, 0, 1 ,MPI_COMM_WORLD, &status);
    partentity0.insert( 0, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 0, 2 ,MPI_COMM_WORLD, &status);
    partentity1.insert( 0, recv);
    
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 2, 1 ,MPI_COMM_WORLD, &status);
    partentity1.insert( 2, recv);
    MPI_Recv(&recv, sizeof(entity *), MPI_BYTE, 2, 2 ,MPI_COMM_WORLD, &status);
    partentity2.insert( 2, recv);

    MPI_Send(&entities[0], sizeof(entity *), MPI_BYTE, 0, 1 ,MPI_COMM_WORLD); 
    MPI_Send(&entities[1], sizeof(entity *), MPI_BYTE, 0, 2 ,MPI_COMM_WORLD); 
    MPI_Send(&entities[1], sizeof(entity *), MPI_BYTE, 2, 1 ,MPI_COMM_WORLD);
    MPI_Send(&entities[2], sizeof(entity *), MPI_BYTE, 2, 2 ,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "PROC" << rank << " :" << std::endl;
    for(auto e : entities ){
      std::cout <<  " has " << e->n << " at " << e;
      auto partobj = partman.getConstPartitionObject(*e);
      auto rmotes = partobj.getRemoteObjectsCollectionRange( );
      std::cout << " remote copies : " <<  rmotes.size() << " : ";
      for (auto r : rmotes ){
	std::cout << " " << r.getProcId() << " " << r.getObjectAddress();
      }
      auto owner = partobj.getOwner ();
      std::cout << " Owner is " <<  owner.getProcId()  << " " <<  owner.getObjectAddress() << std::endl;
      std::cout << " Am I Owner ? "  << partobj.isOwner () << std::endl;
    }
    std::cout << std::flush;
  }
  
  int owned = countOwned(  entities.begin(), entities.end(), partman  );
  std::vector<int>   owned_rcv(size);
  MPI_Gather(&owned, 1, MPI_INT,  owned_rcv.data(), 1, MPI_INT, 0, MPI_COMM_WORLD );
  if(!rank)    for (int i =0; i < size; ++i){
      std::cout<< "PROC"<< i << " own " << owned_rcv[i] << " entities" << std::endl;
    }
  
  
  int locpartok = checkPartition(  entities.begin(), entities.end(),  partman);
  
  int globpartok;
  MPI_Reduce(&locpartok, &globpartok, 1, MPI_INT,  MPI_MIN, 0, MPI_COMM_WORLD);
  if(!rank){
    if (globpartok) {
      std::cout <<  " OK The partition is correct " << std::endl<< std::flush ;
    }
    else{
      std::cout <<  " ERROR The partition is wrong " << std::endl<< std::flush ;
    }
  }
  
  // now we remove entity c from d, and update the partitions manager
  if(rank == 2){
    partman.remove( *entities[0]);
    entities.erase( entities.begin()); 
  }
  
  locpartok = checkPartition(  entities.begin(), entities.end(),  partman);
  MPI_Reduce(&locpartok, &globpartok, 1, MPI_INT,  MPI_MIN, 0, MPI_COMM_WORLD); 
  if(!rank){
    if (globpartok){
      std::cout <<  " ERROR The partition is correct : it's expected to be WRONG !!" << std::endl<< std::flush ;
    }
    else{
      std::cout <<  " OK The partition is wrong, as expected" << std::endl<< std::flush ;
    }
  }
  
  if(rank == 0){
    partman.getPartitionObject(*entities[2]).remove(2);
  }
  if(rank == 3){
    partman.getPartitionObject(*entities[1]).remove(2);
  } 
  
  locpartok = checkPartition(  entities.begin(), entities.end(),  partman);
  MPI_Reduce(&locpartok, &globpartok, 1, MPI_INT,  MPI_MIN, 0, MPI_COMM_WORLD);
  if(!rank) {
    if (globpartok){
      std::cout <<  " OK The partition is correct" << std::endl<< std::flush ;
    }
    else{
      std::cout <<  " ERROR The partition is wrong !!" << std::endl<< std::flush ;
    }
  }
  
  // Testing Copy Constructor
  partman_t partman2(partman);
  if(!rank){
    std::cout << "partman is copied to partman 2 ! Check by counting owned :" << std::endl;  
  }
  
  
  owned = countOwned(  entities.begin(), entities.end(), partman2  );
  //owned_rcv;
  MPI_Gather(&owned, 1, MPI_INT,  owned_rcv.data(), 1, MPI_INT, 0, MPI_COMM_WORLD );
  if(!rank)    for (int i =0; i < size; ++i){
      std::cout<< "PROC"<< i << " own " << owned_rcv[i] << " entities" << std::endl;
    }
  
  // Testing iterator on entity having remote object on other proc
  std::cout<<"For PROC "<<rank<<" entities having remote one on other proc are :";
  for (auto it=partman.beginObject(),ite=partman.endObject(); it!=ite; ++it)
  {
      entity *e=*it;
      std::cout<<" "<<e->n;
  }
  std::cout<<std::endl<<std::flush;

  if(!rank)
  {
      //testing iterator operator itself
      // allready tested above prefix ++, !=, =, *
      auto it_test=partman.beginObject();
      //copy constructor
      partman_t::c_iter_object_t cp_it_test(it_test);
      std::cout<<"copy constructor + = oper test "<<(it_test==cp_it_test)<<std::endl;
      std::cout<<"post fix ++ oper test "<<(it_test==(cp_it_test++)); std::cout<<(*cp_it_test)->n<<std::endl;

  }

  
  partman.clear();
  if(!rank){
    std::cout << "partman is cleared ! :  entity do not have remoteOject from it's point of so that all entities are owned " << std::endl;  
  }
  owned = countOwned(  entities.begin(), entities.end(), partman  );
  MPI_Gather(&owned, 1, MPI_INT,  owned_rcv.data(), 1, MPI_INT, 0, MPI_COMM_WORLD );
  if(!rank)    for (int i =0; i < size; ++i){
      std::cout<< "PROC"<< i << " own " << owned_rcv[i] << " entities" << std::endl;
    }
  
  
  
  MPI_Finalize();  
}

