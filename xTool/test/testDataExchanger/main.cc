#include "xPartitionManager.h"
#include "xUnorderedMapDataManager.h"
#include "xDataExchanger.h"
#include "mpi.h"
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <functional>


class entity {
public: 

  entity( ):n{'a'}{}
  entity(  char _n):n{_n}{}
  char n;
};
template <class T>
  using datman =  xtool::xUnorderedMapDataManager< entity, T > ;

typedef xtool::xPartitionManager< datman > partman_t;

template< template < class  >  class DATAMANAGER >
class entityKeyManagerSendAndRecv  {
 public:
  typedef const entity * information_key_t;
  entityKeyManagerSendAndRecv( const  xtool::xPartitionManager<DATAMANAGER > &_partman ):partman(_partman) {};
  information_key_t localObjectKey(const entity & o){
    return  &o;
  }
  information_key_t remoteObjectKey(const xtool::xRemoteObject<entity > & ro, const entity  &lo ){
    return ro.getObjectAddress();
  }
  xtool::xConstPartitionObject< entity > getConstPartitionObject( const entity & o){
    return partman.getConstPartitionObject(o);
  }
 private:
  const xtool::xPartitionManager<DATAMANAGER > &partman;
};

class entityKeyManagerSendOrRecv  {
 public:
  typedef const entity * information_key_t;
  entityKeyManagerSendOrRecv( const  MPI_Comm &comm){
    MPI_Comm_rank(comm, &mpi_rank);
    MPI_Comm_size(comm, &mpi_size);
  };
  information_key_t localObjectKey(const entity & o){
    return  &o;
  }
  std::set<int> getMessageRanks( const information_key_t & lk){
    const int nbtarget = ((4*(lk->n-'a'+1)))%(mpi_size-1);
    std::set<int> res;
    for (int i = 0; i < nbtarget; ++i){
      int target = (((i+1)*(lk->n-'a'+1)))%mpi_size;
      int k=mpi_rank;
      while (target == mpi_rank) target = (++k)%mpi_size;
      res.insert(target);
    }
    if (!res.empty()){
      std::cout << "Proc " << mpi_rank << " Will send/receive data associated to object "  << lk->n <<  " to/from procs ";
      for(auto ires : res) std::cout << " " <<ires;
      std::cout << std::endl;
    }
    
    return res;
  }
 private:
  int mpi_rank, mpi_size;
};

class msgManager 
{
 public:
  void printMsg()
  {
      std::function<void(const std::string&)> out =[](const std::string& s) -> void
              {
                  std::cout<<s<<std::endl;
              };
      for_each(get_msg.begin(),get_msg.end(),out);
      for_each(set_msg.begin(),set_msg.end(),out);
  }
  void clearMsg()
  {
    get_msg.clear();
    set_msg.clear();
  }
 protected:
  std::set<std::string> set_msg;
  std::set<std::string> get_msg;
};

class randomvectorofdoubleinfomanager : public msgManager {
 public:

  typedef xtool::nonhomogeneous_data_style_trait  data_style_trait;
  typedef xtool::send_and_recv_keys_communication_trait communication_trait ;
  
  randomvectorofdoubleinfomanager( const MPI_Comm &comm  ){
    MPI_Comm_rank(comm, &mpi_rank);
  };

  typedef const entity * information_key_t;
  
  
  void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff , int sendto){
    int nb =  (2*(key->n-'a'+1))%5;
    buff.pack(&nb, 1, MPI_INT );
    std::vector< double >  v(nb);
    double r=45.;
    std::generate(v.begin(), v.end(), [&key,&sendto,&r]()-> double {
            r*=(r+34.)/23.;
            return ((sendto+1.)*(key->n-'a'+8.))*54985./(31.+r);
            });
    buff.pack(v.data(), nb, MPI_DOUBLE);
    std::string s("P"+std::to_string(mpi_rank)+ " Send for "+key->n+" info "+std::to_string(nb));
    for(auto d : v) s+=" "+std::to_string(d);
    s+=" to "+std::to_string(sendto);
    get_msg.insert(s);
  }
  
  size_t getApproxDataSize(){
    return sizeof(int) + 3*sizeof(double);
  }
  
  void setInfo(information_key_t key, const xtool::xMpiOutputBuffer & buff , int receivedfrom){
    int nb;
    buff.unPack(&nb, 1, MPI_INT);
    std::vector<double> v(nb);
    buff.unPack(v.data(), nb, MPI_DOUBLE);
    std::string s("P"+std::to_string(mpi_rank)+ " Receive for "+key->n+" info "+std::to_string(nb));
    for(auto d : v) s+=" "+std::to_string(d);
    s+=" from "+std::to_string(receivedfrom);
    set_msg.insert(s);
  }  
 private:
  int mpi_rank;
};


class randomintinfomanager : public msgManager {
 public:
  
  typedef xtool::homogeneous_data_style_trait  data_style_trait;
  typedef xtool::send_and_recv_keys_communication_trait communication_trait ;
  
  randomintinfomanager( const MPI_Comm &comm  ){
    MPI_Comm_rank(comm, &mpi_rank);
  };

  typedef const entity * information_key_t;
  typedef int            information_t;
  
  information_t getInfo(information_key_t key, int sendto){
    int info = 100*(mpi_rank+1) + 10*sendto + (key->n-'a'+1)%10;
    get_msg.insert("P"+std::to_string(mpi_rank)+ " Send for "+key->n+" info "+std::to_string(info)+" to "+std::to_string(sendto));
    return info;
  }
  void setInfo(information_key_t key, const information_t &info, int receivedfrom){
    set_msg.insert("P"+std::to_string(mpi_rank)+ " Receive for "+key->n+" info "+std::to_string(info)+" from "+std::to_string(receivedfrom));
  }  
 private:
  int mpi_rank;
};
 


class randomintinfomanager_sendonly : public msgManager  {
 public:
  typedef xtool::homogeneous_data_style_trait data_style_trait;
  typedef xtool::send_only_keys_communication_trait  communication_trait;
  
  randomintinfomanager_sendonly( const MPI_Comm &comm  ) {
    MPI_Comm_rank(comm, &mpi_rank);
    MPI_Comm_size(comm, &mpi_size);
  };
  typedef const entity * information_key_t;
  typedef int            information_t;
  
  information_t getInfo(information_key_t key, int sendto){
    int info = 100*(mpi_rank+1) + 10*sendto + (2*(key->n-'a'+1))%10;
    get_msg.insert("P"+std::to_string(mpi_rank)+ " Send for "+key->n+" info "+std::to_string(info)+" to "+std::to_string(sendto));
    return info;
  }
  
  void setInfo( const std::vector<information_t> &infos, int receivedfrom){ 
    std::string s("P"+std::to_string(mpi_rank)+ " Receive  "+std::to_string(infos.size())+" infos from P "+std::to_string(receivedfrom)+":");
    for (const auto& info : infos ) s+=" "+std::to_string(info);
    
    for ( int i =0; i < static_cast<int>(infos.size()); ++i ) {
      int modified_data = infos[i]-mpi_rank-receivedfrom;
      new_data_for_each_setted_info[receivedfrom].push_back(modified_data);
    }
    
    s+=" For each recieved info a int is computed with rank an will be sent back:";
    for (const auto& newinfo : new_data_for_each_setted_info.at(receivedfrom) ) s+=" "+std::to_string(newinfo);

    set_msg.insert(s);
    
  } 
  
  // the following Member are not uses by data exchanger ... They are particular for their local usage in this main.
  const std::map<int, std::vector<int > > & getNewDataCreatedFromEachReceivedInfoFrom( ){
    return new_data_for_each_setted_info;
  };
  
  
 private:
  int mpi_rank, mpi_size;
  std::map<int, std::vector<int > > new_data_for_each_setted_info; 
};


class randomintinfomanager_recvonly : public msgManager  {
 public:
  typedef xtool::homogeneous_data_style_trait data_style_trait;
  typedef xtool::recv_only_keys_communication_trait communication_trait;
  typedef const entity * information_key_t;
  typedef int            information_t;
  randomintinfomanager_recvonly( const MPI_Comm & comm, const std::map<int,  std::vector< int  > > &_send_to): send_to(_send_to)  { 
    MPI_Comm_rank(comm, &mpi_rank);
  };
  
  std::vector<information_t> getInfo( int sendto){ 
    auto it = send_to.find(sendto);
    if (it != send_to.end()) return it->second;
    return emptyinfo;
  }
  
  void setInfo( const information_key_t & key, information_t &info, int recv_from){ 
    set_msg.insert("P"+std::to_string(mpi_rank)+ " Receive for "+key->n+" info "+std::to_string(info)+" from "+std::to_string(recv_from));
  } 
 
  
 private:
  int mpi_rank ;
  std::map< int, std::vector< int  > > send_to;
  std::vector<information_t> emptyinfo;
}; 

class randomvectorofdoubleinfomanager_sendonly : public msgManager {
 public:

  typedef xtool::nonhomogeneous_data_style_trait  data_style_trait;
  typedef xtool::send_only_keys_communication_trait  communication_trait;
  
  randomvectorofdoubleinfomanager_sendonly( const MPI_Comm &comm  ){
    MPI_Comm_rank(comm, &mpi_rank);
  };

  typedef const entity * information_key_t;
  
  
  void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff , int sendto){
    int nb =  (3*(key->n-'a'+1))%5;
    buff.pack(&nb, 1, MPI_INT );
    std::vector< double >  v(nb);
    double r=33.;
    std::generate(v.begin(), v.end(), [&key,&sendto,&r]()-> double {
            r*=(r+64.)/83.;
            return ((sendto+1.)*(key->n-'a'+7.))*78925./(29.+r);
            });
    buff.pack(v.data(), nb, MPI_DOUBLE);
    std::string s("P"+std::to_string(mpi_rank)+ " Send for "+key->n+" info "+std::to_string(nb));
    for(auto d : v) s+=" "+std::to_string(d);
    s+=" to "+std::to_string(sendto);
    get_msg.insert(s);
  }
  
  size_t getApproxDataSize(){
    return sizeof(int) + 3*sizeof(double);
  }
  
  void setInfo(const xtool::xMpiOutputBuffer & buff , int receivedfrom){
    int nb,tot=0;
    std::string s("P"+std::to_string(mpi_rank)+" Receive from P "+std::to_string(receivedfrom)+" ");
    while (!buff.exhausted())
    {
        buff.unPack(&nb, 1, MPI_INT);
        std::vector<double> v(nb);
        buff.unPack(v.data(), nb, MPI_DOUBLE);
        s+= " a block of "+std::to_string(nb)+" infos :";
        for(auto d : v) s+=" "+std::to_string(d);
        tot+=nb;
        for ( int i =0; i < nb; ++i ) {
            double modified_data = v[i]-mpi_rank*10.-receivedfrom*1000.;
            new_data_for_each_setted_info[receivedfrom].push_back(modified_data);
        } 
    }
    if (tot)
    {
        s+=" For each vector values a double is computed with rank and will be sent back when implemented:";
        for (const auto& newinfo : new_data_for_each_setted_info.at(receivedfrom) ) s+=" "+std::to_string(newinfo);
    }
    set_msg.insert(s);
  }  

  // the following Member are not uses by data exchanger ... They are particular for their local usage in this main.
  const std::map<int, std::vector<double > > & getNewDataCreatedFromEachReceivedInfoFrom( ){
    return new_data_for_each_setted_info;
  };
  
  
 private:
  int mpi_rank;
  std::map<int, std::vector<double > > new_data_for_each_setted_info; 
};

//! A small test of xPartition Manager on 4 mpi procs.
// a group of 6 entity  a, b, c, d, e, f are partitionned between 4 proc :
// proc 0 has a, b, c, d
// proc 1 has b, d, f
// proc 2 has c, d, e, f
// proc 3 has a, c, e , g   //note that g is an object which is not on partition boundary
// each proc r create each entity x it has and then send to the procs which also have x the address of x on r. 
// upon receive each proc update their partition manager.
// owned are counted by each processor, then checkPartition is called.
// Then entity c is removed from proc 2
// checkPartition is called -> we expect that it return false for proc 1 and 2
// since the partition manager is not udated yet.
// Then the partition manager is updated (testing xPartitionManager::remove )
// check partition is then called on a copy of the initial partition manager ( the intent is to
//  test the partition manager copy constructor)



int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank);
  MPI_Comm_size( MPI_COMM_WORLD, &size);
  std::string fo = "proc_"+std::to_string(rank) +"_output.txt";
  auto ok = freopen (fo.c_str(),"w",stdout);
  if(!ok){std::cout << "Can't reopen stdout on file "<< fo << __FILE__<< __LINE__ << std::endl; throw;}
  
  if (size < 4){
    if (!rank) std::cout << argv[0] << " Need at least 4 mpi process to run "<< std::endl ; 
    MPI_Finalize();  
    return 0;
  }

  /*std::string fo = "proc_"+std::to_string(rank) +"_output.txt";
  auto ok = freopen (fo.c_str(),"w",stdout);
  if(!ok){std::cout << "Can't reopen stdout on file "<< fo << __FILE__<< __LINE__ << std::endl; throw;}*/
  partman_t partman(MPI_COMM_WORLD);
  
  std::function < void ( const entity & ) > print_entity = [&partman]( const entity & e){
    std::cout << " has " << e.n;
    auto partobj = partman.getConstPartitionObject(e);
    auto rmotes =  partobj.getRemoteObjectsCollectionRange( );
    std::cout << " number of remote copies : " <<  rmotes.size() << " on procs : ";
    for (auto r : rmotes ){
      std::cout << " " << r.getProcId();
    }
    auto owner = partobj.getOwner ( );
    std::cout << " Owner is " <<  owner.getProcId()  << std::endl;
    std::cout << " Am I Owner ? "  << partobj.isOwner () << std::endl;
  };
  
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
    entities.push_back(new entity('g'));
    
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
    for(auto e : entities ) print_entity(*e);
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
    for(auto e : entities )  print_entity(*e);
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
    
    for(auto e : entities ) print_entity(*e);
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
    for(auto e : entities ) print_entity(*e);
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

  if(!rank)    std::cout <<"Partman is build" << std::endl;
  
  // send_and_recv_keys_association_trait ====================================================================================================
  entityKeyManagerSendAndRecv<datman > keymansar(partman);
  xtool::xKeyContainerSendAndRecv<  entityKeyManagerSendAndRecv< datman  >::information_key_t > key_container_int(partman.getComm());

  randomintinfomanager infomanint(partman.getComm());

  // send_and_recv_keys_communication_trait
  // homogeneous_data_style_trait
  std::cout <<"Testing exchange OwnerGather for int data" << std::endl;
  key_container_int.accumulateKeysOwnerGather( entities.begin(), entities.end(), keymansar );
  xtool::exchangeInformation(key_container_int,infomanint);
  infomanint.printMsg();
  key_container_int.clearKeys();
  infomanint.clearMsg();
  std::cout <<"OwnerGather Done" << std::endl;
  
  std::cout <<"Testing exchange OwnerScatter for int data" << std::endl; 

  key_container_int.accumulateKeysOwnerScatter( entities.begin(), entities.end(), keymansar );
  xtool::exchangeInformation(key_container_int,infomanint);
  infomanint.printMsg();
  key_container_int.clearKeys();
  infomanint.clearMsg();
  std::cout <<"OwnerScatter Done" << std::endl;

  std::cout <<"Testing exchange AllGather for int data" << std::endl; 

  key_container_int.accumulateKeysAllGather( entities.begin(), entities.end(), keymansar );
  xtool::exchangeInformation(key_container_int,infomanint);
  infomanint.printMsg();
  std::cout <<"AllGather Done" << std::endl;

  // send_and_recv_keys_communication_trait
  // nonhomogeneous_data_style_trait
  randomvectorofdoubleinfomanager infomanvecdouble(partman.getComm());

  xtool::xKeyContainerSendAndRecv<entityKeyManagerSendAndRecv<datman >::information_key_t> key_container_vec_double(partman.getComm());
  std::cout <<"Testing exchange OwnerGather for pseudo-random size vector of pseudo-random double data" << std::endl;
  key_container_vec_double.accumulateKeysOwnerGather( entities.begin(), entities.end(), keymansar );
  xtool::exchangeInformation(key_container_vec_double,infomanvecdouble);
  infomanvecdouble.printMsg();
  key_container_vec_double.clearKeys();
  infomanvecdouble.clearMsg();
  std::cout <<"OwnerGather Done" << std::endl;
  
  std::cout <<"Testing exchange OwnerScatter for pseudo-random size vector of pseudo-random double data" << std::endl; 
  key_container_vec_double.accumulateKeysOwnerScatter( entities.begin(), entities.end(),keymansar );
  xtool::exchangeInformation(key_container_vec_double,infomanvecdouble);
  infomanvecdouble.printMsg();
  infomanvecdouble.clearMsg();
  std::cout <<"OwnerScatter Done" << std::endl;
 
  // send_or_recv_keys_association_trait====================================================================================================
  entityKeyManagerSendOrRecv keymansor(partman.getComm());
  xtool::xKeyContainerSendOrRecv<  entityKeyManagerSendOrRecv::information_key_t > key_container_oneside_int(partman.getComm());

  // send_only_keys_communication_trait
  // homogeneous_data_style_trait
  std::cout <<"Testing exchange Sender Only Keys for int to pseudo-random targets" << std::endl;
  randomintinfomanager_sendonly infomanrandomintrandomtargets(partman.getComm() );
  key_container_oneside_int.accumulateKeys( entities.begin(), entities.end(), keymansor );
  xtool::exchangeInformation(key_container_oneside_int,infomanrandomintrandomtargets);
  infomanrandomintrandomtargets.printMsg();
  infomanrandomintrandomtargets.clearMsg();
  std::cout <<" exchange Sender Only Keys for int to pseudo-random targets Done" << std::endl;
  
  // recv_only_keys_communication_trait
  // homogeneous_data_style_trait
  std::cout <<"Testing exchange Receive Only Keys for int from sources defined above" << std::endl;
  randomintinfomanager_recvonly infomanrandomintrandomsources(partman.getComm(),infomanrandomintrandomtargets.getNewDataCreatedFromEachReceivedInfoFrom( ) );
  xtool::exchangeInformation(key_container_oneside_int,infomanrandomintrandomsources);
  infomanrandomintrandomsources.printMsg();
  infomanrandomintrandomsources.clearMsg();
  std::cout <<" exchange Receive Only Keys for int from pseudo-random sources Done" << std::endl;
  
  // send_only_keys_communication_trait
  // nonhomogeneous_data_style_trait
  std::cout <<"Testing exchange Sender Only Keys for  pseudo-random size vector of pseudo-random double data to pseudo-random targets" << std::endl;
  randomvectorofdoubleinfomanager_sendonly infomanrandomvectdoublerandomsources(partman.getComm());
  xtool::exchangeInformation(key_container_oneside_int,infomanrandomvectdoublerandomsources);
  infomanrandomvectdoublerandomsources.printMsg();
  infomanrandomvectdoublerandomsources.clearMsg();
  std::cout <<" exchange Sender Only Keys for double vect to pseudo-random targets Done" << std::endl;
  MPI_Finalize();  
}
