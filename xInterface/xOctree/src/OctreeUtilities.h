/* 
   xfem : C++ Finite Element Library
   developed under the GNU Lesser General Public License
   See the NOTICE & LICENSE files for conditions.
*/
  
#ifndef OCTREE_UTILITIES_H
#define OCTREE_UTILITIES_H

#include "oOctree.h"
#include "xIntegrationRule.h"
#include "xFemMatrix.h"
#include "xCommandOnGeomElem.h"
#include "xValueManager.h"
#include "xFiniteElement.h"
#include "xLevelSet.h"




namespace xinterface {
  namespace xoctree {

    // struct scalings{
    const int scaling_pow1[16]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} ;
    const int scaling_pow2[16] ={1, 2, 4, 8, 16, 32, 64, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864, 268435456, 1073741824} ;
    const int scaling_pow4[16] ={1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864, 268435456, 1073741824} ;
    // };

    template <class BILINEARFORM, class FIELD, class ITER, class ASSEMBLER, template <class > class DATAMANAGER =xfem::xMesh::datamanager_t >
      void OctreeAssemble(BILINEARFORM& bilin,
			  ASSEMBLER& assembler,
              const xfem::xIntegrationRule& integration_rule,
			  const FIELD& test,
			  const FIELD& trial,
			  ITER it, ITER end,
			  const int * scaling2, const int * scaling3,
              xfem::xEntityToEntity integ2appro_test,
              xfem::xEntityToEntity integ2appro_trial,
              const DATAMANAGER<int > &octreeLevel = xfem::xMesh::get_const_octree_level()
                          )
    {
      bool first=true;
      xfem::xFemMatrix<double> matrix;
      int levelSaved=-10; //set arbitrary to dummy value -10 for -Wmaybe-uninitialized

      const bool debug =  false;   //xDebugSingleton::instance().flag;
      xfem::xIntegrateFormCommand command(&bilin);
      if (&trial == &test) {
          if(1) std::cout<<"Assemble : TEST==TRIAL\n";
          xfem::xValueManagerDist<double>* double_manager = trial.getValueManager();
          for(; it != end; ++it){
              AOMD::mEntity*e_integ  = *it;
              AOMD::mEntity*e_appro  = integ2appro_test(e_integ);
              xfem::xFiniteElement FEM;
              FEM.setKeysAndFcts(e_appro, trial.begin(), trial.end(), "trial");
              if (FEM.sizeKey("trial"))
              {
                  //        const bool isSym=bilin.isSym(FEM.getFcts("trial"), FEM.getFcts("trial"));
                  if(debug) std::cout<<"assembler.ass-1\n";
                  xfem::xGeomElem geo_appro(e_appro);

                  if(first){
                      bilin.init(&geo_appro, FEM.getFcts("trial"));
                      integration_rule.accept(command, e_integ);
                      bilin.finalize();
                      //if(isSym) bilin.getMatrix().symmetrize();
                      matrix=bilin.getMatrix();
                      const int *pLevelSaved = octreeLevel.getData(*e_integ);
                      levelSaved=(pLevelSaved)?*pLevelSaved:0;
                      first=false;
                  }

                  std::vector<xfem::xValue<double>*> vals;
                  double_manager->getValPtr(FEM.beginKey("trial"), FEM.endKey("trial"), vals);
                  if(debug) std::cout<<"assembler.assemble\n";

                  const int *poLevel =  octreeLevel.getData(*e_integ);
                  const int olevel = (poLevel)?(*poLevel):0;
                  const int levelDiff=olevel -levelSaved;
                  double coeff=1.;

                  if(e_integ->getLevel()==2){
                      if(levelDiff<0) coeff=scaling2[abs(levelDiff)];
                      else coeff=1./scaling2[abs(levelDiff)];
                      assembler.setCoeff(coeff);
                  }else{
                      if(levelDiff<0) coeff=scaling3[abs(levelDiff)];
                      else coeff=1./scaling3[abs(levelDiff)];
                      assembler.setCoeff(coeff);
                  }


                  // 	  if(e_integ->getLevel()==3){
                  //        const int *poLevel =  octreeLevel.getData(*e_integ);
                  //        const int olevel = (poLevel)?(*poLevel):0;
                  //        const int levelDiff=olevel -levelSaved;
                  // 	    double coeff=1.;
                  // 	    if(levelDiff<0) coeff=scaling3[abs(levelDiff)];
                  // 	    else coeff=1./pow_base2[abs(levelDiff)];
                  // 	    assembler.setCoeff(coeff);
                  // 	  }
                  assembler.assemble_matrix(vals.begin(), vals.end(), matrix);

                  // cout<<"matrix "<<bilin.getMatrix()<<endl;
                  // cout<<"level saved="<<levelSaved<<endl;
                  // const int *pcurrent_level = octreeLevel.getData(*e_integ);
                  // const int current_level = (pcurrent_level)?(*pcurrent_level):0;
                  // cout<<"current level="<< current_level <<endl;
                  // cout<<"diff="<<levelDiff<<endl;
                  // cout<<"areaCoeff="<<areaCoeff<<endl;
                  // cout<<"matrix "<<matrix<<endl;

              }
          }
          first=true;
          assembler.setCoeff(1.0);
      }
      else {
          if(1) std::cout<<"Assemble : TEST # TRIAL\n";
          xfem::xValueManagerDist<double>* double_manager_trial = trial.getValueManager();
          xfem::xValueManagerDist<double>* double_manager_test  =  test.getValueManager();
          for(; it != end; ++it){
              AOMD::mEntity*e_integ  = *it;
              AOMD::mEntity*e_appro_test   = integ2appro_test(e_integ);
              AOMD::mEntity*e_appro_trial  = integ2appro_trial(e_integ);
              if (debug)
              {
                  std::cout << "e_integ is " << std::endl;
                  e_integ->print();
                  std::cout << "e_appro_test is " << std::endl;
                  e_appro_test->print();
                  std::cout << "e_appro_trial is " << std::endl;
                  e_appro_trial->print();
              }
              xfem::xFiniteElement FEM_trial;
              FEM_trial.setKeysAndFcts(e_appro_trial, trial.begin(), trial.end(), "trial");
              if (debug) std::cout << "size trial " << FEM_trial.sizeKey("trial") << std::endl;
              xfem::xFiniteElement FEM_test;
              FEM_test.setKeysAndFcts(e_appro_test, test.begin(),   test.end(), "test");
              if (debug) std::cout << "size test  " << FEM_test.sizeKey("test") << std::endl;
              if (FEM_trial.sizeKey("trial") && FEM_test.sizeKey("test"))
              {
                  xfem::xGeomElem geo_appro_test(e_appro_test);
                  xfem::xGeomElem geo_appro_trial(e_appro_trial);

                  if(first){
                      bilin.init(&geo_appro_test, FEM_test.getFcts("test"), &geo_appro_trial, FEM_trial.getFcts("trial"));
                      integration_rule.accept(command, e_integ);
                      bilin.finalize();
                      matrix=bilin.getMatrix();
                      const int *pLevelSaved = octreeLevel.getData(*e_integ);
                      levelSaved=(pLevelSaved)?*pLevelSaved:0;
                      first=false;
                  }

                  std::vector<xfem::xValue<double>*> valstrial, valstest;
                  double_manager_trial->getValPtr(FEM_trial.beginKey("trial"), FEM_trial.endKey("trial"), valstrial);
                  double_manager_test->getValPtr(FEM_test.beginKey("test"),  FEM_test.endKey("test"),  valstest);
                  if (debug) std::cout << "size val trial " << valstrial.size() <<
                                          "size val test  " << valstest.size() << std::endl;

                  const int *poLevel =  octreeLevel.getData(*e_integ);
                  const int olevel = (poLevel)?(*poLevel):0;
                  const int levelDiff=olevel -levelSaved;
                  double coeff=1.;

                  if(e_integ->getLevel()==2){
                      if(levelDiff<0) coeff=scaling2[abs(levelDiff)];
                      else coeff=1./scaling2[abs(levelDiff)];
                      assembler.setCoeff(coeff);
                  }else{
                      if(levelDiff<0) coeff=scaling3[abs(levelDiff)];
                      else coeff=1./scaling3[abs(levelDiff)];
                      assembler.setCoeff(coeff);
                  }


                  // 	  if(e_integ->getLevel()==3){
                  // 	  int levelDiff=e_integ->getAttachedInt(octreeLevel_tag)-levelSaved;
                  // 	  double coeff=1.;
                  // 	  if(levelDiff<0) coeff=pow_base2[abs(levelDiff)];
                  // 	  else coeff=1./pow_base2[abs(levelDiff)];
                  // 	  assembler.setCoeff(coeff);
                  // 	  }


                  assembler.assemble_matrix(valstest.begin(),  valstest.end(),
                                            valstrial.begin(), valstrial.end(),
                                            matrix);

              }
          }
          first=true;
          assembler.setCoeff(1.0);
      }

      return;
    }

    template <class BILINEARFORM, class FIELD, class ITER, class ASSEMBLER> 
      void OctreeAssemble(BILINEARFORM& bilin,
			  ASSEMBLER & assembler,
              const xfem::xIntegrationRule& integration_rule,
			  const FIELD& test,
			  const FIELD& trial,
			  ITER it, ITER end,
			  const int * scaling2, const int * scaling3,
              xfem::xEntityToEntity integ2appro  = xtool::xIdentity<AOMD::mEntity*>())
    {
      OctreeAssemble(bilin, assembler, integration_rule, test, trial, it, end, scaling2, scaling3,
		     integ2appro, integ2appro);
    } 


    //Caution : this one is VERY inefficient...
    class xLevelSetTo_ls_xyz_t_Adaptator{       
    public:
        xLevelSetTo_ls_xyz_t_Adaptator(xfem::xLevelSet &ls_);
        double operator()(const double&x, const double&y, const double&z);

    private:
        xfem::xLevelSet &ls;
        xfem::xMesh *mesh;
    };


    class xPointToDoubleTo_ls_xyz_t_Adaptator{
    public:
    xPointToDoubleTo_ls_xyz_t_Adaptator(xfem::xPointToDouble &xpd_);
      double operator()(const double&x, const double&y, const double&z);
    private:
      xfem::xPointToDouble &xpd;
    };
  } //end namespace xoctree
} // end namespace xinterface



#endif
