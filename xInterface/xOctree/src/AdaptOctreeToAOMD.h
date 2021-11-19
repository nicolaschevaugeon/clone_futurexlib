/*
   xfem : C++ Finite Element Library
   developed under the GNU Lesser General Public License
   See the NOTICE & LICENSE files for conditions.
*/

#ifndef ADAPT_OCTREEDB_TO_AOMD__H
#define ADAPT_OCTREEDB_TO_AOMD__H

#include <string>

#include "oKeyManager.h"
#include "oLevelSet.h"
#include "oOctree.h"
#include "xAlgorithm.h"
#include "xIntegrationRule.h"
#include "xLevelSet.h"
#include "xMesh.h"
#include "xSubMesh.h"
#include "xSubMeshManager.h"

// Pour exporter les niveaux du maillage :
class evalLevel : public xfem::xEval<double>
{
  public:
   evalLevel(const xinterface::aomd::xAttachedDataManagerAOMD<int> &octreeLevelContainer_ = xfem::xMesh::get_const_octree_level())
       : octreeLevelContainer{octreeLevelContainer_}
   {
   }

   void operator()(const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ, result_type &result) const override
   {
      const int *poctreeLevel = octreeLevelContainer.getData(*geo_appro->getEntity());
      result = (poctreeLevel) ? (*poctreeLevel) : 0;
   }

  private:
   const xinterface::aomd::xAttachedDataManagerAOMD<int> &octreeLevelContainer;
};

// InterfaceOctreeToAOMD2G
// attachCell
// oCoarseningCriteriaOnEvaluator
// xAttachableCell
// refinementCriteriaCommand
// oRefinementCriteriaOnEvaluator
// ReadLsetFromITK2D
// ReadLsetFromITK3D

namespace xinterface
{
namespace xoctree
{
using ::xoctree::oField;
using ::xoctree::oKeyManager;
using ::xoctree::oLevelSet;
using ::xoctree::oMapping;
using ::xoctree::oOctree;
using ::xoctree::oRefinementCriteria;

/// \brief Interface for conversion from the octree DB to AOMD mesh. Focus on mesh adaptation
/// \param[in] fullModify : if false only the minimal AOMD database will be constructed. WARNING: when set to false, the
/// boundaries are not tagged ! \param[in] simplex : if true, Quad/Hex are replaced by Tri/Tet \param[in] hhg : if false mesh
/// modification to remove hanging nodes
template <template <class> class DATAMANAGER = xfem::xMesh::datamanager_t>
void InterfaceOctreeToAOMDG(const oOctree &octree, const oField &field, oLevelSet &o_ls, xfem::xMesh &mesh, xfem::xLevelSet &ls,
                            DATAMANAGER<int> &octreeLevel = xfem::xMesh::get_octree_level(),
                            DATAMANAGER<AOMD::mEntity *> &isHangingOn = xfem::xMesh::get_is_hanging_on(),
                            const bool fullModify = true, const bool simplex = true, bool hgg = true);

/// Attach cell to an entity
inline void attachCell(AOMD::mEntity *Entity, unsigned int tagNum, const oOctree::cell_type *cell);

/// \brief Read ITK levelset from a file (2D ONLY !)
/// \param[out] lsVec : Discrete LS values vector (row major ?)
/// \param[out] sizeLS : number of pixels along x and y
/// \param[out] stepPx : pixel size along x and y
/// \param[out] origin : origin of the domain
void ReadLsetFromITK2D(string filename, std::vector<double> &lsVec, std::vector<int> &sizeLs, std::vector<double> &stepPx,
                       std::vector<double> &origin);

/// \brief Read ITK levelset from a file (3D ONLY !) NEW FORMAT ! (.lso)
/// \param[out] lsVec : Discrete LS values vector (row major ?)
/// \param[out] sizeLS : number of pixels along x, y and z
/// \param[out] stepPx : pixel size along x, y and z
/// \param[out] origin : origin of the domain
/// \param[in]  invertSide : invert in and out
void ReadLsetFromITK3D(string filename, std::vector<double> &lsVec, std::vector<int> &sizeLs, std::vector<double> &stepPx,
                       std::vector<double> &origin, bool invertSide = false);

/// Attachable cells (cf aomd)
class xAttachableCell : public AOMD::mAttachableData
{
  public:
   const oOctree::cell_type *cell;
};

/// \brief Command used to refine the mesh.
/// \param[in] mesh, error evaluator
/// \param[out] cells to be refined

template <template <class> class DATAMANAGER = xfem::xMesh::datamanager_t>
class refinementCriteriaCommand : public xfem::xCommandOnGeomElem
{
  public:
   refinementCriteriaCommand(xfem::xMesh *mesh_, xfem::xEval<int> &eval_,
                             std::set<std::pair<const oOctree::cell_type *, int>> *cellsContext_,
                             const DATAMANAGER<int> &octreeLevelContainer_ = xfem::xMesh::get_const_octree_level())
       : mesh(mesh_), eval(eval_), cellsContext(cellsContext_), octreeLevelContainer(octreeLevelContainer_)
   {
      // submesh_manager.allocateSubsetEntities("toRefine");
      mesh->createSubMesh("toRefine");
      creatorCell_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("creator_cell");
   }

   void execute(xfem::xGeomElem *geo_integ) override
   {
      const bool debug = false;
      const int nb = geo_integ->GetNbIntegrationPoints();
      xfem::xSubMesh &toRefine = mesh->getSubMesh("toRefine");
      for (int k = 0; k < nb; k++)
      {
         geo_integ->setUVW(k);
         if (geom_appro != geo_integ)
            geom_appro->setUVWForXYZ(geo_integ->getXYZ());
         else
            geom_appro->setUVW(k);

         // if(debug) geom_appro->getEntity()->print();
         int res = -1;
         eval(geom_appro, geo_integ, res);

         if (res == 1)
         { /*Cellule a raffiner*/
            AOMD::mEntity *e_appro = geom_appro->getEntity();
            if (debug) e_appro->print();
            if (debug) cout << "Element to Refine\n";
            // 	submesh_manager.add_sub(geom_appro->getEntity(),"toRefine");
            toRefine.add(e_appro);

            for (int i = 0; i < e_appro->size(0); ++i)  // Only 2D !!!
            {
               AOMD::mEntity *vert = e_appro->get(0, i);
               toRefine.add(vert);
               // 		submesh_manager.add_sub(vert,"toRefine");
            }

            xAttachableCell *ac = (xAttachableCell *)(geom_appro->getEntity())->getData(creatorCell_tag);
            if (ac)
            {
               const int *plevel = octreeLevelContainer.getData(*e_appro);
               int level = (plevel) ? (*plevel) : 0;

               std::pair<const oOctree::cell_type *, int> context(ac->cell, level);
               cellsContext->insert(context);
               //			cout<<"aLevel="<<level<<endl;
               //			cout<<"add cell "<< ac->cell << " no be refined\n";
            }
            else
            {
               cout << "The mesh was not Tagged !\n";
               assert(1 == 0);
            }

            break; /*si un pdg depasse le critere, pas la peine d'aller voir les autres...*/
         }
      }
   }

  private:
   xfem::xMesh *mesh;
   // xSubMeshManager submesh_manager;
   xfem::xEval<int> &eval;
   std::set<std::pair<const oOctree::cell_type *, int>> *cellsContext;
   unsigned int creatorCell_tag;
   // unsigned int octreeLevel_tag;
   const DATAMANAGER<int> &octreeLevelContainer;
};

/// Coarsenind with evaluator : To be improved
class oCoarseningCriteriaOnEvaluator : public oRefinementCriteria
{
  public:
   oCoarseningCriteriaOnEvaluator(xfem::xEval<int> &eval_) : eval(eval_) {}

   template <class ITER>
   int findCellsToFreeze(xfem::xMesh *mesh, xfem::xIntegrationRule &integration_rule, ITER it, ITER end)
   {
      refinementCriteriaCommand<> commandDummy(mesh, eval, &cellsContext);
      xfem::ApplyCommandOnIntegrationRule(commandDummy, integration_rule, it, end);
      (mesh->getSubMesh("toRefine")).exportGmsh("toFreeze.msh");
      foundCells = cellsContext.size();
      return foundCells;
   }

   bool operator()(oOctree::const_iterator cell, int level, const int *ijk, oOctree::const_iterator children_beg,
                   oOctree::const_iterator children_end) const override
   {
      std::pair<const oOctree::cell_type *, int> pair(cell, level);
      return (!(cellsContext.find(pair) != cellsContext.end()));
   }

  private:
   std::set<std::pair<const oOctree::cell_type *, int>> cellsContext;
   int foundCells;
   xfem::xEval<int> &eval;
};

// Old implementation : still usefull ?
#if 0
    class oRefinementCriteriaOnEvaluator : public oRefinementCriteria
    {
    public:
    oRefinementCriteriaOnEvaluator(xEval<int> &eval_): eval(eval_) {}
    oRefinementCriteriaOnEvaluator(): eval(xfem::xEvalConstant<int>(1)) {}

      template <class ITER>
	int findCellsToRefine(xfem::xMesh *mesh,xIntegrationRule &integration_rule,ITER it, ITER end)
	{
	  refinementCriteriaCommand commandDummy(mesh,eval,&cellsContext);
	  xfem::ApplyCommandOnIntegrationRule(commandDummy,integration_rule, it,end);
	  mesh->export_sub("toRefine.msh", "toRefine");
	  foundCells=cellsContext.size();
	  return foundCells;
	}

      bool operator()(oOctree::const_iterator cell, int level, 
		      const int * ijk,
		      oOctree::const_iterator children_beg,
		      oOctree::const_iterator children_end) const
      {
	std::pair<const oOctree::cell_type*, int> pair(cell,level);
	return (cellsContext.find(pair)!=cellsContext.end());
      }

      bool inContext(std::pair<const oOctree::cell_type*, int> &paire){
	return (cellsContext.find(paire) != cellsContext.end());
      };

    private:
      std::set<std::pair<const oOctree::cell_type*, int>  > cellsContext;
      int foundCells;
      xfem::xEval<int> &eval;
    };
#endif

/// Refinement class used by the octree to refine the database.
/// @todo Virer le code du header et mettre dans le .cc : non, car templatise...
class oRefinementCriteriaOnEvaluator : public oRefinementCriteria
{
  public:
   oRefinementCriteriaOnEvaluator(xfem::xEval<int> *eval_ = nullptr) : eval(eval_) {}

   template <class ITER>
   int findCellsToRefine(xfem::xMesh *mesh, xfem::xIntegrationRule &integration_rule, ITER beg, ITER end)
   {
      if (!eval)
      {
         cout << "No evaluator !!!\n";
         return 0;
      }
      refinementCriteriaCommand<> commandDummy(mesh, *eval, &cellsContext);
      //      xRegion all(mesh);
      //       ApplyCommandOnIntegrationRule(commandDummy,integration_rule, all.begin(),all.end());
      xfem::ApplyCommandOnIntegrationRule(commandDummy, integration_rule, beg, end);
      mesh->getSubMesh("toRefine").exportGmsh("toRefine.msh");

      foundCells = cellsContext.size();
      return foundCells;
   }

   bool operator()(oOctree::const_iterator cell, int level, const int *ijk, oOctree::const_iterator children_beg,
                   oOctree::const_iterator children_end) const override
   {
      std::pair<const oOctree::cell_type *, int> pair(cell, level);
      //      cout<<"tLevel="<<level<<endl;
      //      cout<<"Trouve ???="<<(cellsContext.find(pair)!=cellsContext.end())<<endl;
      return (cellsContext.find(pair) != cellsContext.end());
   }

   bool inContext(std::pair<const oOctree::cell_type *, int> &paire) { return (cellsContext.find(paire) != cellsContext.end()); };

   void setEvaluator(xfem::xEval<int> *eval_) { eval = eval_; }

   void clearCellsToRefine() { cellsContext.clear(); }

   int translateToParents(oOctree &octree)
   {
      std::set<std::pair<const oOctree::cell_type *, int>> newContext;
      std::set<std::pair<const oOctree::cell_type *, int>>::iterator it = cellsContext.begin();
      std::set<std::pair<const oOctree::cell_type *, int>>::iterator itE = cellsContext.end();

      for (; it != itE;)
      {
         int level = it->second;
         const oOctree::cell_type *cellSon = it->first;
         cout << "level=" << level << endl;
         const oOctree::cell_type *cellFather = octree.getAncestor(1, cellSon, level);
         std::pair<const oOctree::cell_type *, int> context(cellFather, level - 1);
         newContext.insert(context);
         it = cellsContext.erase(it);
         ++it;
      }

      it = newContext.begin();
      itE = newContext.end();
      for (; it != itE; ++it)
      {
         cellsContext.insert(*it);
      }

      return newContext.size();
   }

   /// Experimental (?) add one layer of cells
   //    int addLayer(oOctree &octree){

   //      const bool debug=false;
   //      int Dim=octree.getDim();
   //      int on_what = Dim-1;
   //      bool flag;

   //      std::set<std::pair<const oOctree::cell_type*, int>  > layerContext;
   //      std::set<std::pair<const oOctree::cell_type*, int>  >::iterator it=cellsContext.begin();
   //      std::set<std::pair<const oOctree::cell_type*, int>  >::iterator ite=cellsContext.end();

   //      for(;it!=ite;++it){//boucle cells a raffiner...
   //        //         cout<<"cell num="<<num++<<endl;

   //        std::pair<const oOctree::cell_type*, int> ctx=*it;

   //        int ijk[3];
   //        int ijk_up[3];
   //        fill(ijk, ijk+3, 0);
   //        int l=ctx.second;
   //        octree.octree2cartesian(ctx.first,l,ijk);

   //        // on_what = Dim-1 => Dim=3 on_what=2=face, Dim=2 on_what=1=edge; in both case nextNeighbor dont use id_neighbor =>
   //        id_neighbor set to zero
   //        // nota : A.S. 6/01/12 this was introduced in a blunt maner. No check to see if loop on neigbors is requested or not
   //        in 3D for this methode
   //        // to be check
   //        int id_neighbor=0;
   //        oTopo topo=octree.getTopo();
   //        int count = 0;
   //        int ijkn[3];
   //        //         int nei=0;
   //        while(topo.nextNeighbor(ijk, l, ijkn, on_what, flag, Dim, octree.getPeriodicity(), count,id_neighbor))
   //        {
   //          //c'est un set (elements uniques): pas la peine de se faire suer, on ajoute tous les voisins sans se poser de
   //          question !

   //          oOctree::cell_type* neighbor = octree.cartesian2octree(ijkn, l);
   //          //           cout<<"neighbour="<<neighbor<<endl;

   //          layerContext.insert(make_pair(neighbor,l));
   //          count++;
   //        }
   //      }

   //      it=layerContext.begin();
   //      ite=layerContext.end();
   //      for(;it!=ite;++it){
   //        cellsContext.insert(*it);
   //      }

   //      return cellsContext.size();
   //    }

  private:
   std::set<std::pair<const oOctree::cell_type *, int>> cellsContext;
   int foundCells;
   xfem::xEval<int> *eval;
};

/// Refines the mesh one time
template <template <class> class DATAMANAGER = xfem::xMesh::datamanager_t>
void buildOneLevelFinerMesh(oOctree &octree, const oField &lsoct_field, oLevelSet &o_lset, oKeyManager &key_manager,
                            xfem::xMesh &meshFine, xfem::xLevelSet &lsFine,
                            DATAMANAGER<int> &octreeLevel = xfem::xMesh::get_octree_level(),
                            DATAMANAGER<AOMD::mEntity *> &isHangingOn = xfem::xMesh::get_is_hanging_on(),
                            const bool fullModify = true, const bool simplex = true, bool hgg = true, int nb = 1);

}  // namespace xoctree

}  // namespace xinterface
#include "AdaptOctreeToAOMD_imp.h"

#endif
