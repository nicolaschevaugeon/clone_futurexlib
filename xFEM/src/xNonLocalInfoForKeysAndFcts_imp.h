/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
:    and conditions.
*/

//#ifndef _XNONLOCALINFO_IMP_H_
//#define _XNONLOCALINFO_IMP_H_

#include "xAssembler.h"
#include "xCSRVector.h"
#include "xDenseMatrix.h"
#include "xFiniteElement.h"
#include "xForm.h"
#include "xGenericSparseMatrix.h"
#include "xGraphMatrix.h"
#include "xIntegrationRule.h"
#include "xNonLocalInfoForKeysAndFcts.h"
#include "xOperators.h"

//#define TIMING_MONITORING 1
//#define CALLGRIND_MONITORING 1
#ifdef TIMING_MONITORING
#include "xDeltaTime.h"
#endif
#ifdef CALLGRIND_MONITORING
#include <valgrind/callgrind.h>
#endif

namespace xfem
{
// to avoid to many inclusion xAlgorithm.h is not use ans a direct declaration is taken
// for setDofNum. TODO : make this thing clean (i.e. move setDofNum in a correct location xValue ?? xState ??)
extern int setDofNum(xValue<double> *v);

template <typename S, template <class> class DATAMANAGER>
xNonLocalInfoGeneratorForKeysAndFctsForHangingNode<S, DATAMANAGER>::xNonLocalInfoGeneratorForKeysAndFctsForHangingNode(
    const DATAMANAGER<AOMD::mEntity *> &isHangingOn_, const DATAMANAGER<std::vector<AOMD::mEntity *>> &downGroup_,
    const DATAMANAGER<std::vector<AOMD::mEntity *>> &bndGroup_, double _thresold)
    : xNonLocalInfoGeneratorForKeysAndFcts(),
      isHangingOn(isHangingOn_),
      downGroup(downGroup_),
      bndGroup(bndGroup_),
      thresold(_thresold),
      pnli_data(new DATAMANAGER<xNonLocalInfoForKeysAndFcts *>),
      nli_data(*pnli_data)
{
}

template <typename S, template <class> class DATAMANAGER>
xNonLocalInfoGeneratorForKeysAndFctsForHangingNode<S, DATAMANAGER>::xNonLocalInfoGeneratorForKeysAndFctsForHangingNode(
    xMesh &_mesh, const DATAMANAGER<AOMD::mEntity *> &isHangingOn_, const DATAMANAGER<std::vector<AOMD::mEntity *>> &downGroup_,
    const DATAMANAGER<std::vector<AOMD::mEntity *>> &bndGroup_, double _thresold)
    : xNonLocalInfoGeneratorForKeysAndFcts(_mesh),
      isHangingOn(isHangingOn_),
      downGroup(downGroup_),
      bndGroup(bndGroup_),
      thresold(_thresold),
      pnli_data(new DATAMANAGER<xNonLocalInfoForKeysAndFcts *>),
      nli_data(*pnli_data)
{
}

template <typename S, template <class> class DATAMANAGER>
xNonLocalInfoGeneratorForKeysAndFctsForHangingNode<S, DATAMANAGER>::xNonLocalInfoGeneratorForKeysAndFctsForHangingNode(
    xMesh &_mesh, DATAMANAGER<xNonLocalInfoForKeysAndFcts *> &_nli_data, const DATAMANAGER<AOMD::mEntity *> &isHangingOn_,
    const DATAMANAGER<std::vector<AOMD::mEntity *>> &downGroup_, const DATAMANAGER<std::vector<AOMD::mEntity *>> &bndGroup_,
    double _thresold)
    : xNonLocalInfoGeneratorForKeysAndFcts(_mesh),
      isHangingOn(isHangingOn_),
      downGroup(downGroup_),
      bndGroup(bndGroup_),
      thresold(_thresold),
      pnli_data(nullptr),
      nli_data(_nli_data)
{
}

template <typename S, template <class> class DATAMANAGER>
const xNonLocalInfoForKeysAndFcts *xNonLocalInfoGeneratorForKeysAndFctsForHangingNode<S, DATAMANAGER>::getNonLocalInfo(
    const AOMD::mEntity &e) const
{
   xNonLocalInfoForKeysAndFcts *const *pnli = nli_data.getData(e);
   if (pnli) return *pnli;
   return nullptr;
}

template <typename S, template <class> class DATAMANAGER>
void xNonLocalInfoGeneratorForKeysAndFctsForHangingNode<S, DATAMANAGER>::clearNonLocalInfoForKeysAndFcts()
{
   if (generated)
   {
      clearNonLocalInfoContainer();
      int dimp1 = dim + 1;
      for (int i = 1; i < dimp1; ++i)
      {
         for (AOMD::mEntity *pe : mesh->range(i))
         {
            // nullifyAttachedNonLocalInfoForKeysAndFcts(*it,nli_tag);
            nli_data.deleteData(*pe);
         }
      }
      generated = false;
   }
}
template <typename S, template <class> class DATAMANAGER>
xNonLocalInfoGeneratorForKeysAndFctsForHangingNode<S, DATAMANAGER>::~xNonLocalInfoGeneratorForKeysAndFctsForHangingNode()
{
   clearNonLocalInfoForKeysAndFcts();
}

template <typename S, template <class> class DATAMANAGER>
void xNonLocalInfoGeneratorForKeysAndFctsForHangingNode<S, DATAMANAGER>::generateNonLocalInfoForKeysAndFcts(
    spacePtr_t space, const int polynomial_order)
{
   // This is comming from DeclareInterpolationHanging
   // It is more restrictive as no recursion mecanism is implemented
   // DeclareInterpolationHanging autorize recursive nested relation : when there is more then
   // one jump level for neighbour element in octree.
   // Here we consider that the rule "no more then on jump level is possible" is strictly folowed by "Octree to AOMD" interface

   //
   if (generated) return;

   // local
   int numdof_master, numdof_slave, numdof;
   int dimp1 = dim + 1;
#ifdef CALLGRIND_MONITORING
   CALLGRIND_START_INSTRUMENTATION
#endif
#ifdef TIMING_MONITORING
   xtool::xDeltaTime dte;
#endif

   // tag from Octree to AOMD
   // const unsigned int id_is_hanging_on = xMesh::get_is_hanging_on_tag();
   // const unsigned int id_down_groupe = xMesh::get_down_groupe_tag();
   // const unsigned int id_bnd_groupe = xMesh::get_bnd_groupe_tag();

   // matrix container for L2 projection resolution
   typedef xlinalg::xGenericSparseMatrix<double,
                                         xTraitMatrixDefinitePositive,  // matrix is defined positive
                                         xTraitMatrixUnSym,  // storing only lower part of sym matrix as Taucs Pastix and Mumps REFACTOR: unsym for superlu
                                                                // wan't it
                                         typename S::matrix_storage, typename S::matrix_indexing>
       l2_matrix_t;
   // assembler
   xAssemblerBasic<l2_matrix_t, xlinalg::xDenseMatrix, double> assembler_mass;
   xAssemblerBasic<xlinalg::xDenseMatrix, xlinalg::xCSRVector, double> assembler_rhs;

   // formulation for L2 projection problem
   xFormBilinearWithoutLaw<xValOperator<xtool::xIdentity<double>>, xValOperator<xtool::xIdentity<double>>> l2_bilin;

   // Integration
   xIntegrationRuleBasic l2_ir(2 * polynomial_order);
   xIntegrateFormCommand command(&l2_bilin);

   // space iterator
   // nota : this look rather strange but here is what all this is about :
   //    - by creating a vector it give a simple way to obtain iterators on spacePtr
   //    - having this iterator we can bypass field creation for setKeysAndFcts methode of xFiniteElement :
   //      field in this context yould have been only created (at every hanging entity) to give space to setKeysAndFcts.
   std::vector<xFiniteElement::spacePtr> space_container(1, space);
   std::vector<xFiniteElement::spacePtr>::iterator space_begin = space_container.begin();
   std::vector<xFiniteElement::spacePtr>::iterator space_end = space_container.end();

   // means to get upper/master element  from slave entity
   // xAttached slave2master(id_is_hanging_on);
   xEntityToEntity slave2master = [this](AOMD::mEntity *es) {
      AOMD::mEntity *const *ppmaster = isHangingOn.getData(*es);
      if (ppmaster) return const_cast<AOMD::mEntity *>(*ppmaster);
      AOMD::mEntity *pe = nullptr;
      return pe;
   };
   // adding here a filter via  xRegion is straitforward but does it have much sens ??
#ifdef TIMING_MONITORING
   int iddt1 = dte.initAccu("get attached def");
   int iddt2 = dte.initAccu("slave/master dof generation");
   int iddt3 = dte.initAccu("matrix graph generation");
   int iddt4 = dte.initAccu("mass and RHS matrix generation");
   int iddt5 = dte.initAccu("solve");
   int iddt6 = dte.initAccu("slave<->mater relations generation");
   int iddt7 = dte.initAccu("master->slave storing");
   int iddt8 = dte.initAccu("attaching loop");
#endif
   // node iterator
   xIter it = mesh->begin(0);
   xIter itend = mesh->end(0);
   // loop on nodes to generate NLI definition
   for (; it != itend; ++it)
   {
      // set as a pointer
      AOMD::mEntity *nh = (*it);

      // if hanging
      if (slave2master(nh))
      {
#ifdef TIMING_MONITORING
         dte.startAccu(iddt1);
#endif
         // local
         std::vector<xValKey> slave_keys;
         std::vector<xValKey> master_keys;

         // get attached vector
         std::vector<AOMD::mEntity *> &slave = const_cast<std::vector<AOMD::mEntity *> &>(*downGroup.getData(*nh));
         std::vector<AOMD::mEntity *> &exclude = const_cast<std::vector<AOMD::mEntity *> &>(*bndGroup.getData(*nh));

         // set iterator
         std::vector<AOMD::mEntity *>::iterator it_slave_elem = slave.begin();
         std::vector<AOMD::mEntity *>::iterator it_slave_elem_beg = it_slave_elem;
         std::vector<AOMD::mEntity *>::iterator it_slave_elem_end = slave.end();
         std::vector<AOMD::mEntity *>::iterator it_exclude_elem = exclude.begin();
         std::vector<AOMD::mEntity *>::iterator it_exclude_elem_beg = it_exclude_elem;
         std::vector<AOMD::mEntity *>::iterator it_exclude_elem_end = exclude.end();

         // volatil double manager
         xValueManagerDist<double> local_master_val_manager;
         xValueManagerDist<double> local_slave_val_manager;

         // FEM and vals maps
         std::map<std::vector<AOMD::mEntity *>::iterator, std::vector<xValue<double> *>> slave_vals;
         std::map<std::vector<AOMD::mEntity *>::iterator, xFiniteElement *> slave_FEM;
         std::map<AOMD::mEntity *, std::vector<xValue<double> *>> master_vals;
         std::map<AOMD::mEntity *, xFiniteElement *> master_FEM;

         // initiate a empty NLI definition
         xNonLocalInfoForKeysAndFcts *nli = new xNonLocalInfoForKeysAndFcts();

#ifdef TIMING_MONITORING
         dte.endAccu(iddt1);
         dte.startAccu(iddt2);
#endif
         // loop on slave to generate slave and master dofs in local double manager
         // and fill FEM/vals map
         for (it_slave_elem = it_slave_elem_beg; it_slave_elem != it_slave_elem_end; ++it_slave_elem)
         {
            // local to loop
            std::vector<xValue<double> *> vals;
            xFiniteElement *FEM = new xFiniteElement;

            // finding keys and function in space
            FEM->setKeysAndFcts((*it_slave_elem), space_begin, space_end, "slave");

            AOMD::mEntity *master_elem = slave2master(*it_slave_elem);

            if (!master_elem)
            {
               cout << "File " << __FILE__ << " line " << __LINE__ << " Master element is not attached to slave element " << endl;
               throw;
            }
            // generate master key if not yet done
            else if (master_FEM.find(master_elem) == master_FEM.end())
            {
               // local to if
               std::vector<xValue<double> *> valsm;
               xFiniteElement *FEMm = new xFiniteElement;

               // finding keys and function for master element "master"
               FEMm->setKeysAndFcts(master_elem, space_begin, space_end, "master");

               // loop on keys
               xFiniteElement::iterKey itkem = FEMm->endKey("master");
               for (xFiniteElement::iterKey itkm = FEMm->beginKey("master"); itkm != itkem; ++itkm)
               {
                  xValue<double> *vm;
                  xValKey &keym = (*itkm);

                  // generate std value in local val manager
                  if (!(vm = (xValue<double> *)local_master_val_manager.insert(keym, creator_reg)))
                  {
                     cout << "File " << __FILE__ << " line " << __LINE__
                          << " It appends !!! generation of a xValue<T> failed for master key in local val manager " << (*itkm)
                          << endl;
                     throw;
                  }

                  // set state of std value as dof using numdof
                  numdof = local_master_val_manager.size("master") + 1;
                  if (!vm->getState())
                  {
                     local_master_val_manager.add(vm, "master");
                     vm->setState(new xStateOfValueDof(numdof));
                     master_keys.push_back(keym);
                  }
                  // store value
                  valsm.push_back(vm);
               }

               // store for loop below
               master_FEM.insert(make_pair(master_elem, FEMm));
               master_vals.insert(make_pair(master_elem, valsm));
            }

            // create slave values in ValManager
            xFiniteElement::iterKey itke = FEM->endKey("slave");
            for (xFiniteElement::iterKey itk = FEM->beginKey("slave"); itk != itke; ++itk)
            {
               xValue<double> *v;
               xValKey &key = (*itk);

               // if not already created in local slave val manager insert std value
               if (!(v = (xValue<double> *)local_slave_val_manager.insert(key, creator_reg)))
               {
                  cout << "File " << __FILE__ << " line " << __LINE__
                       << " It appends !!! generation of a xValue<T> failed for slave key in val manager" << key << endl;
                  throw;
               }
               numdof = local_slave_val_manager.size("slave") + 1;
               if (!v->getState())
               {
                  local_slave_val_manager.add(v, "slave");
                  v->setState(new xStateOfValueDof(numdof));
                  slave_keys.push_back(key);
               }
               vals.push_back(v);
            }

            // store for loop below
            slave_FEM.insert(make_pair(it_slave_elem, FEM));
            slave_vals.insert(make_pair(it_slave_elem, vals));
         }

#ifdef TIMING_MONITORING
         dte.endAccu(iddt2);
         dte.startAccu(iddt3);
#endif

         // number of dofs
         numdof_slave = local_slave_val_manager.size("slave");
         numdof_master = local_master_val_manager.size("master");

         // setting RHS container
         xlinalg::xDenseMatrix B(numdof_slave, numdof_master);
         assembler_rhs.setTarget(B);

         // setting solution container
         xlinalg::xDenseMatrix SOL(numdof_slave, numdof_master);

         // loop on slave to generate matrix graph
         xlinalg::xGraphMatrix graph(numdof_slave, numdof_slave, 0);//REFACTOR: Sym should be taken from solver
         for (it_slave_elem = it_slave_elem_beg; it_slave_elem != it_slave_elem_end; ++it_slave_elem)
         {
            std::vector<xValue<double> *> &vals_slave = slave_vals.find(it_slave_elem)->second;
            std::vector<int> dof_num(vals_slave.size());
            std::transform(vals_slave.begin(), vals_slave.end(), dof_num.begin(), setDofNum);

            int i, nb_dofs = dof_num.size();

            // treating dofs if any
            if (nb_dofs > 1)
            {
               // sort numdof in ascending order
               std::vector<int>::iterator itbeg = dof_num.begin();
               std::vector<int>::iterator itend = dof_num.end();
               sort(itbeg, itend);

               // loop on pseudo colone of the elementary matrix
               for (i = 0; i < nb_dofs; ++i)
               {
                  graph.addLinesSymBlock(i, nb_dofs, dof_num[i], &dof_num[0]);
               }
            }
            else
            {
               i = dof_num[0];
               graph.add(i, i);
            }
         }

#ifdef TIMING_MONITORING
         dte.endAccu(iddt3);
         dte.startAccu(iddt4);
#endif

         // setting nnz
         graph.countNNZ();

         // setting matrix container
         l2_matrix_t A(graph);
         assembler_mass.setTarget(A);

         // loop on slave to generate mass matrix and RHS
         for (it_slave_elem = it_slave_elem_beg; it_slave_elem != it_slave_elem_end; ++it_slave_elem)
         {
            // finding keys and values for slave element
            xFiniteElement *FEM_slave = slave_FEM.find(it_slave_elem)->second;
            std::vector<xValue<double> *> &vals_slave = slave_vals.find(it_slave_elem)->second;

            // integration of slave for mass matrix
            xGeomElem geo_slave((*it_slave_elem));
            l2_bilin.init(&geo_slave, FEM_slave->getFcts("slave"));
            l2_ir.accept(command, (*it_slave_elem));
            assembler_mass.assemble_matrix(vals_slave.begin(), vals_slave.end(), l2_bilin.getMatrix());

            // get master element attached to this slave. No check as it have already been tested above
            AOMD::mEntity *master_elem = slave2master(*it_slave_elem);

            // finding keys and values for master element
            xFiniteElement *FEM_master = master_FEM.find(master_elem)->second;
            std::vector<xValue<double> *> &vals_master = master_vals.find(master_elem)->second;

            // integration of slave on master for RHS matrix
            xGeomElem geo_master(master_elem);
            l2_bilin.init(&geo_slave, FEM_slave->getFcts("slave"), &geo_master, FEM_master->getFcts("master"));
            l2_ir.accept(command, (*it_slave_elem));
            assembler_rhs.assemble_matrix(vals_slave.begin(), vals_slave.end(), vals_master.begin(), vals_master.end(),
                                          l2_bilin.getMatrix());

            // free slave memory
            delete FEM_slave;
         }

#ifdef TIMING_MONITORING
         dte.endAccu(iddt4);
         dte.startAccu(iddt5);
#endif

         // free master memory
         std::map<AOMD::mEntity *, xFiniteElement *>::iterator itMFe = master_FEM.end();
         for (std::map<AOMD::mEntity *, xFiniteElement *>::iterator itMF = master_FEM.begin(); itMF != itMFe; ++itMF)
            delete ((*itMF).second);

         // solve problem
         l2_solver.connectMatrix(A);
         l2_solver.solve(B, SOL);

#ifdef TIMING_MONITORING
         dte.endAccu(iddt5);
         dte.startAccu(iddt6);
#endif

         // master -> slave relation container
         std::vector<xNonLocalInfoForKeysAndFcts::femRelKeys> tmp_slaves(numdof_master);

         // slave entity from space definition
         std::set<AOMD::mEntity *> slave_from_key;

         // loop on slave dofs to generate in xNonLocalInfoForKeysAndFcts :
         //          slave->master relations
         //          master->slave relations
         for (numdof = 0; numdof < numdof_slave; ++numdof)
         {
            xValKey &key = slave_keys[numdof];

            AOMD::mEntity *slave_entity = key.getEnti();
            slave_from_key.insert(slave_entity);

            // filter excluded entities from relation creation
            if (std::find(it_exclude_elem_beg, it_exclude_elem_end, slave_entity) == it_exclude_elem_end)
            {
               // master key to be stored in xNonLocalInfoForKeysAndFcts for this slave key
               typename xSpace::femKeys tmp_masters;

               // reserve memory for maximum non nul coefficient eaven if  some (most) are null
               tmp_masters.reserve(numdof_master);

               for (int j = 0; j < numdof_master; ++j)
               {
                  const double r = SOL(numdof, j);

                  if (r > thresold || r < -thresold)
                  {
                     tmp_masters.push_back(master_keys[j]);
                     tmp_slaves[j].insert(std::make_pair(key, r));
                  }
               }

               // add slave->master relations in xNonLocalInfoForKeysAndFcts
               nli->addSlave(key, tmp_masters);
            }
         }

#ifdef TIMING_MONITORING
         dte.endAccu(iddt6);
         dte.startAccu(iddt7);
#endif

         nli->reserveMaster(numdof_master);

         // store  master -> slave relation
         for (int j = 0; j < numdof_master; ++j)
         {
            if (tmp_slaves[j].size()) nli->addMaster(master_keys[j], tmp_slaves[j]);
         }

#ifdef TIMING_MONITORING
         dte.endAccu(iddt7);
         dte.startAccu(iddt8);
#endif

         // container to minimise merge usage
         std::set<AOMD::mEntity *> treated;
         std::pair<set<AOMD::mEntity *>::iterator, bool> ret;
         std::map<xNonLocalInfoForKeysAndFcts *, xNonLocalInfoForKeysAndFcts *> new_non_local;
         std::map<xNonLocalInfoForKeysAndFcts *, xNonLocalInfoForKeysAndFcts *>::iterator new_non_local_end = new_non_local.end();
         std::map<xNonLocalInfoForKeysAndFcts *, xNonLocalInfoForKeysAndFcts *>::iterator itfind;

         // loop on slave entity to tag them and attache xNonLocalInfoForKeysAndFcts or merge/create
         std::set<AOMD::mEntity *>::iterator itss;
         std::set<AOMD::mEntity *>::iterator itssend = slave_from_key.end();
         for (itss = slave_from_key.begin(); itss != itssend; ++itss)
         {
            AOMD::mEntity *slave_entity = *itss;

            // filter excluded entities from relation creation
            if (std::find(it_exclude_elem_beg, it_exclude_elem_end, slave_entity) == it_exclude_elem_end)
            {
               // attache nli definition to current slave
               // attachNonLocalInfoForKeysAndFcts(slave_entity, nli_tag, nli);
               nli_data.setData(*slave_entity) = nli;
               // and it's surounding entity
               // nota : by using treated merging is donne only one time per entity
               for (int i_dim = slave_entity->getType() + 1; i_dim < dimp1; ++i_dim)
               {
                  int i_e_tot = slave_entity->size(i_dim);
                  for (int i_e = 0; i_e < i_e_tot; ++i_e)
                  {
                     AOMD::mEntity *adj = slave_entity->get(i_dim, i_e);

                     ret = treated.insert(adj);

                     if (ret.second)
                     {
                        // xNonLocalInfoForKeysAndFcts *old_nli = getAttachedNonLocalInfoForKeysAndFcts(adj, nli_tag);
                        xNonLocalInfoForKeysAndFcts **pold_nli = nli_data.getData(*adj);
                        if (pold_nli && *pold_nli != nli)
                        {
                           xNonLocalInfoForKeysAndFcts *old_nli = *pold_nli;
                           if ((itfind = new_non_local.find(old_nli)) != new_non_local_end)
                           {
                              // attachNonLocalInfoForKeysAndFcts(adj, nli_tag, itfind->second);
                              nli_data.setData(*adj) = itfind->second;
                           }
                           else
                           {
                              xNonLocalInfoForKeysAndFcts *new_nli = new xNonLocalInfoForKeysAndFcts(*old_nli);
                              new_nli->merge(nli);
                              nli_container.push_back(new_nli);
                              new_non_local.insert(make_pair(old_nli, new_nli));
                              new_non_local_end = new_non_local.end();
                              // attachNonLocalInfoForKeysAndFcts(adj, nli_tag, new_nli);
                              nli_data.setData(*adj) = new_nli;
                           }
                        }
                        else if (!pold_nli)
                        {
                           // attachNonLocalInfoForKeysAndFcts(adj, nli_tag, nli);
                           nli_data.setData(*adj) = nli;
                        }
                     }
                  }
               }
            }
         }

#ifdef TIMING_MONITORING
         dte.endAccu(iddt8);
#endif

         // add nli definition to container
         nli_container.push_back(nli);
      }
   }

#ifdef TIMING_MONITORING
   dte.print();
#endif
#ifdef CALLGRIND_MONITORING
   CALLGRIND_STOP_INSTRUMENTATION;
#endif

   generated = true;
}

}  // namespace xfem

//#endif
