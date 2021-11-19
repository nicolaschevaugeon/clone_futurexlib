/* 
    This file is a part of eXlibris C++ Library
    under the GNU General Public License:
    See the LICENSE.md files for terms and 
    conditions.
*/
#include <fstream>
#include <iostream>
#include <vector>
#include <memory>

// eXlibris_tools
#include "xMPIEnv.h"

// Trzellis
#include "mAOMD.h"

// xMeshQueryInterface
#include "aomdMeshQueryInterface.h"
#include "aomdMeshModifierInterface.h"
#include "xEntity.h"
#include "xEntityIterators.h"

using namespace xtool;
using namespace xinterface::xmeshinterface;
using namespace AOMD;


int main(int argc, char *argv[])  
{  

    xMPIEnv::init(argc,argv);
  
  //=========================================================
 

  AOMD::mMesh                                           mesh;
  xinterface::xmeshinterface::aomdMeshQueryInterface    query(mesh);

  fprintf(stderr, "Starting reading the master data file\n");
  AOMD::AOMD_Util::Instance()->import("data/cube.msh", &mesh);
  fprintf(stderr, "Done with reading the master data file\n");

  cout<<" begin modifyAllState()"<<endl;
  mesh.modifyState(2,1,true);
  mesh.modifyState(2,0,true);
  mesh.modifyState(1,0,true);
  mesh.modifyState(0,1,true);
  mesh.modifyState(0,2,true);
  mesh.modifyState(1,2,true);
  cout<<" end modifyAllState()"<<endl;
    

  enum tags {  entity_tag,  int_tag , double_tag , pointer_tag , string_tag , unique_ptr_tag, double_ptr_tag};

  auto itf1 = query.beginFace() ;
  auto itf2 = itf1;
  itf2++;

  xEntity e1=*(itf1);
  xEntity e2=*(itf2);

  std::cout<<std::endl;


  while (itf1 !=  itf2 ) {
      // INT

      cout<<" 1. attach<int> = 999"<<std::endl;
      e1.attachData<int>(int_tag, 999)  ;

      std::cout<<std::endl<<" 2. attach<xEntity>" <<std::endl;
      e2.print();
      e1.attachData<xEntity>(entity_tag, e2)  ;

      // POINTER
       cout<<std::endl<<" 3. attach<std::string*> = \"tentative\": " ;
      std::string* pt=  new string("tentative");
      e1.attachData<string*>(pointer_tag, pt);
      std::cout<<std::endl;

      std::cout<<std::endl<<" 4. attach<double> = 1.e-3"<<std::endl<<std::endl;
      e1.attachData<double>(double_tag, 1.e-3)  ;
      std::cout<<std::endl;

      // string
      cout<<std::endl<<" 5. attach<std::string> \"tentative\": " ;
      std::string st="tentative2";
      e1.attachData<std::string>(string_tag, st);
      std::cout<<std::endl;

/*
      // UNIQUE POINTER
      cout<<std::endl<<" 6. attach<unique_ptr<double>> = 123.456: " ;
      unique_ptr<double> uptr ( new  double(123.456) );
      std::cout<< *uptr<<std::endl;
      e1.attachData< unique_ptr<double> >(unique_ptr_tag, std::move(uptr) );
      std::cout<<std::endl;
*/

      cout<<std::endl<<" 7. attach<double*> = 546.321: " ;
      double* dptr ( new  double(546.321) );
      std::cout<< *dptr<<std::endl;
      e1.attachData< double* >(double_ptr_tag, dptr );
      std::cout<<std::endl;
      itf1++;
  };

  // INT
  std::cout<<std::endl<<"----------- out of scope ------------"<<std::endl;
  std::cout<<std::endl;
  cout<<" 1. getData<int> = "<< e1.getData<int>(int_tag) <<std::endl ;

  //xEntity
  std::cout<<std::endl;
  std::cout << " 2. getData<xEntity> = " <<std::endl;
  (e1.getData<xEntity>(entity_tag)).print() ;

  // pointer
  std::cout<<std::endl;
  cout<<" 3. getData<string*> = \""<< *(e1.getData<std::string*>(pointer_tag)) <<"\""<<std::endl ;

  // pointer
  std::cout<<std::endl;
  cout<<" 4. getData<double> = "<< e1.getData<double>(double_tag) <<std::endl ;

  // string
  std::cout<<std::endl;
  cout<<" 5. getData<string> = \""<< e1.getData<string>(string_tag) <<"\""<<std::endl ;

  // pointer
  std::cout<<std::endl;
  cout<<" 6.0 getData<unique_ptr<double>> (once)             = "<< *(e1.getData<unique_ptr<double>>(unique_ptr_tag)) <<std::endl ;
  cout<<" 6.1 getData<unique_ptr<double>> (twice)            = "<< *(e1.getData<unique_ptr<double>>(unique_ptr_tag)) <<std::endl ;
  double* tmp=(e1.getData<unique_ptr<double>>(unique_ptr_tag)).get();
  cout<<" 6.2 double* ptr =getData<unique_ptr<double>>.get() = "<< *tmp <<std::endl ;
//  AOMD::mEntity* pe = any_cast<AOMD::mEntity*>(  e1.getEntityIdentifier());
  e1.deleteData<unique_ptr<double>>(unique_ptr_tag);
  cout<<" 6.3 ptr after delete                                 = *("<< tmp << ") = " << *tmp <<std::endl ;

  // pointer
  std::cout<<std::endl;
  cout<<" 7.0 getAtached<double*>            = "<< *(e1.getData< double*>(double_ptr_tag)) <<std::endl ;
  double* tmp2=e1.getData<double*>(double_ptr_tag);
  cout<<" 7.2 double* ptr =getData<double*> = "<< *tmp2 <<std::endl ;
  e1.deleteData<double*>(unique_ptr_tag);
  cout<<" 7.3 ptr after delete                                 = *("<< tmp2 << ") = " << *tmp2 <<std::endl ;





  //=========================================================
  
  // wait every destruction done
  MPI_Barrier(MPI_COMM_WORLD);
  return xMPIEnv::finalize();
  
}
