#include "xMPIEnv.h"
//#include "mpi.h"
#include "xLoadBalanceTools.h"
#include "mMesh.h"
#include "mAOMD.h"

#include <fstream>
void printMeshStat(const AOMD::mMesh &m){
  std::cout << "Mesh statistc for current partition : " << std::endl;
  std::cout << " V "<< m.size(0)
	    << " E "<< m.size(1)
	    << " F "<< m.size(2)
	    << " R "<< m.size(3) 
	    << std::endl;
}

int main(int argc, char *argv[]){
  xtool::xMPIEnv::init(argc,argv);
  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank);
  MPI_Comm_size( MPI_COMM_WORLD, &size);
  //std::cout << "### I'M PROC "<< rank  << std::endl; 
  
  std::string fo = "proc_"+std::to_string(rank) +"_output.txt";
  auto ok = freopen (fo.c_str(),"w",stdout);
  if(!ok){std::cout << "Can't reopen stdout on file "<< fo << __FILE__<< __LINE__ << std::endl; throw;}
  std::cout << "### I'M PROC "<< rank  << std::endl; 
  
 
  std::string meshfilename("data/test.msh");
  if (argc > 1){ 
    meshfilename = argv[1]; 
  }
  
  AOMD::mMesh m;
  xmeshtool::partitionManager part_man(MPI_COMM_WORLD);
    
  if (!rank){
    std::ifstream meshfile(meshfilename.c_str());
    if (!meshfile.is_open()) {
      std::cout << "can't open meshfile "<< meshfilename << " in xMesh Constructor " << __FILE__ << " " << __LINE__ << std::endl;
      throw;
    }
    meshfile.close();   
    AOMD::AOMD_Util::Instance()->import(meshfilename.c_str(), &m);
    m.modifyState(3,2, true);
    m.modifyState(2,3, true);
    m.modifyState(2,1, true);
    m.modifyState(1,2, true);
    m.modifyState(1,0, true);
    m.modifyState(0,1, true); 
    m.modifyState(3,1, false);
    m.modifyState(3,0, false);
    m.modifyState(2,0, false);
    m.modifyState(1,3, false);
    m.modifyState(0,2, false);
    m.modifyState(0,3, false);
    std::cout << "Mesh Loaded on proc 0" << std::endl;
  }
  
  
  printMeshStat(m);
  const size_t start = meshfilename.find_last_of("/")+1;
  const std::string outfilename_base = meshfilename.substr(start, meshfilename.size()-4-start );
  
  int nb_part = std::max(size/2,1);
  std::cout <<"## FIRST PASS RANDOM PARTITIONING ON "<< nb_part<< " PARTITIONS" << std::endl;
  srand(1);  // to make sure that we reporduce alwys the same partition !
  xinterface::aomd::xAttachedDataManagerAOMD< int > target_proc = xmeshtool::setRandomPartition(m,  nb_part );
  xmeshtool::migrateEntities(m, part_man, target_proc);
  printMeshStat(m);
  std::string outfilename = outfilename_base + "_FIRST_PATH_RANDOM";
  xmeshtool::exportGMSH_V2_Dist(m, part_man, outfilename.c_str());
  
  std::cout <<"## SECOND PASS METIS PARTITIONNING ON "<< nb_part << " PARTITIONS"  << std::endl;
  target_proc=xmeshtool::setParmetisPartition(m, part_man,  nb_part  );
  xmeshtool::migrateEntities(m, part_man, target_proc);
  printMeshStat(m);
  outfilename = outfilename_base + "_SECOND_PATH_PARMETIS";
  xmeshtool::exportGMSH_V2_Dist(m, part_man, outfilename.c_str());
  
  nb_part =1;
  std::cout <<"## THIRD PASS METIS PARTITIONNING -> BACK TO "<< nb_part << " PARTITION"  << std::endl;
  target_proc=xmeshtool::setParmetisPartition(m, part_man,  nb_part  );
  xmeshtool::migrateEntities(m, part_man, target_proc);
  printMeshStat(m);
  outfilename = outfilename_base + "_THIRD_PATH_PARMETIS";
  xmeshtool::exportGMSH_V2_Dist(m, part_man, outfilename.c_str());
  
  nb_part =size;
  std::cout <<"## FOURTH PASS METIS PARTITIONNING " << nb_part << " PARTITIONS" << std::endl;
  target_proc=xmeshtool::setParmetisPartition(m, part_man,  nb_part  );
  xmeshtool::migrateEntities(m, part_man, target_proc);
  printMeshStat(m);
  outfilename = outfilename_base + "_FOURTH_PATH_PARMETIS";
  xmeshtool::exportGMSH_V2_Dist(m, part_man, outfilename.c_str());
  
  // int pok = checkPartition (m.begin(0), m.end(0),  part_man);
  //std::cout <<" Check part " << pok << std::endl;
  
  return xtool::xMPIEnv::finalize();
}
