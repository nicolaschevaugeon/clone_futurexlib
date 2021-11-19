/*
  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE & LICENSE files for conditions.
*/

#include <algorithm>
#include <functional>

#include "InterfaceOctreeToAOMD.h"
#include "mEdge.h"
#include "mFace.h"
#include "mHex.h"
#include "mTet.h"
#include "oExport.h"
#include "oKeyManager.h"
#include "oOctree.h"
#include "xLevelSet.h"
#include "xMesh.h"

namespace xinterface
{
namespace xoctree
{
// xmesh
// octree_level
// is_hanging_on
// is_hanging_by
// down_group
// bnd_groupe
template <template <class> class DATAMANAGER>
void InterfaceOctreeToAOMD(const oOctree& octree, const oField& lsoct, xfem::xMesh& mesh, xfem::xLevelSet& ls, const bool simplex,
                           iClassifyCriteria* ptr_clfCriteria, bool alternate, DATAMANAGER<int>& octreeLevel,
                           DATAMANAGER<AOMD::mEntity*>& isHangingOn, DATAMANAGER<AOMD::mEntity*>& isHangingBy,
                           DATAMANAGER<std::vector<AOMD::mEntity*>>& downGroup,
                           DATAMANAGER<std::vector<AOMD::mEntity*>>& bndGroup)
{
   const bool debug = false;
   const oMapping& mapping = octree.getMapping();
   const oTopo& topo(octree.getTopo());
   const int dim = octree.getDim();
   const oKeyManager& key_manager = lsoct.getKeyManager();
   AOMD::mMesh& mmesh = mesh.getMesh();
   for (const oKey* key : key_manager)
   {
      double xyz[3];
      mapping.ijk2xyz(key->getIJK(), octree.getLevelMax(), xyz);
      if (debug)
      {
         const int* ijk = key->getIJK();
         cout << " vertex created with id " << key->getId() << " at ijk_fine " << endl;
         std::copy(ijk, ijk + 3, std::ostream_iterator<int>(std::cout, " "));
         cout << endl;
         cout << " at location " << endl;
         std::copy(xyz, xyz + 3, std::ostream_iterator<double>(std::cout, " "));
         cout << endl;
      }
      mmesh.createVertex(key->getId(), xyz[0], xyz[1], xyz[2], nullptr);
   }

   for (int l = 0; l <= octree.getLevelMax(); l++)
   {
      oOctree::const_iterator it = octree.begin(l);
      oOctree::const_iterator ite = octree.end(l);
      int ids[8];
      for (; it != ite; ++it)
      {
         if (*it == 1)
         {
            int ijk[3];
            octree.octree2cartesian(it, l, ijk);
            key_manager.getNodeIdsOnElt(ijk, l, ids);
            oOctree::const_iterator children;
            if (l == octree.getLevelMax())
               children = octree.end();
            else
               children = octree.getChildren(it, l);
            int physical_index = init_physical_index;
            if (nullptr != ptr_clfCriteria)
               physical_index = (*ptr_clfCriteria)(it, l, ijk, children, children + octree.getNbChildren());
            auto pGeom = mmesh.getGEntity(physical_index, octree.getDim());
            if (octree.getDim() == 2)
            {
               if (!simplex)
               {
                  // quad case
                  AOMD::mFace* e1 = mmesh.createFaceWithVertices(ids[0], ids[1], ids[2], ids[3], pGeom);
                  octreeLevel.setData(*e1) = l;
               }
               else
               {
                  // simplex case
                  if ((ijk[0] + ijk[1]) % 2 == 0)
                  {
                     AOMD::mFace* e1 = mmesh.createFaceWithVertices(ids[0], ids[1], ids[2], pGeom);
                     AOMD::mFace* e2 = mmesh.createFaceWithVertices(ids[0], ids[2], ids[3], pGeom);
                     octreeLevel.setData(*e1) = l;
                     octreeLevel.setData(*e2) = l;
                  }
                  else
                  {
                     AOMD::mFace* e1 = mmesh.createFaceWithVertices(ids[0], ids[1], ids[3], pGeom);
                     AOMD::mFace* e2 = mmesh.createFaceWithVertices(ids[1], ids[2], ids[3], pGeom);
                     octreeLevel.setData(*e1) = l;
                     octreeLevel.setData(*e2) = l;
                  }
               }
            }
            else
            {
               if (!simplex)
               {
                  // hex
                  AOMD::mHex* hex =
                      mmesh.createHexWithVertices(ids[0], ids[1], ids[2], ids[3], ids[4], ids[5], ids[6], ids[7], pGeom);
                  octreeLevel.setData(*hex) = l;
               }
               else
               {
                  // Maillage non alterne
                  AOMD::mTet* tet1 = mmesh.createTetWithVertices(ids[0], ids[1], ids[5], ids[3], pGeom);
                  AOMD::mTet* tet2 = mmesh.createTetWithVertices(ids[1], ids[2], ids[3], ids[5], pGeom);
                  AOMD::mTet* tet3 = mmesh.createTetWithVertices(ids[2], ids[3], ids[5], ids[6], pGeom);
                  AOMD::mTet* tet4 = mmesh.createTetWithVertices(ids[0], ids[3], ids[4], ids[5], pGeom);
                  AOMD::mTet* tet5 = mmesh.createTetWithVertices(ids[3], ids[4], ids[5], ids[7], pGeom);
                  AOMD::mTet* tet6 = mmesh.createTetWithVertices(ids[3], ids[5], ids[6], ids[7], pGeom);
                  octreeLevel.setData(*tet1) = l;
                  octreeLevel.setData(*tet2) = l;
                  octreeLevel.setData(*tet3) = l;
                  octreeLevel.setData(*tet4) = l;
                  octreeLevel.setData(*tet5) = l;
                  octreeLevel.setData(*tet6) = l;
               }
            }
         }
      }
   }

   if (simplex)  // note : why only for simplex ??
   {
      mmesh.modifyState(3, 1, true);
   }

   // classify corners
   //   int ijk_ori[3] = {0,0,0};
   //   int ijk_corner[3];
   //   int offset = 1;
   //   int count = 0;
   //   while(topo.next_ijk_node_on_element(ijk_corner, 0, ijk_ori, offset, dim, count))
   //     {
   //       int node_id = key_manager.getNodeId(ijk_corner, 0);
   //       mVertex* v = mesh.getVertex(node_id);
   //       if (debug) cout << " classifying a corner node on geometry id " << count << " for node " << node_id << endl;
   //       v->classify(mesh.getGEntity(++count,0));
   //     }

   // New...
   int ijk_corner[3];
   int count = 0;
   for (int iC = 0; iC < topo.pow_base2[dim]; ++iC)
   {
      ijk_corner[0] = topo.ijk_corner_node_level0[iC][0];
      ijk_corner[1] = topo.ijk_corner_node_level0[iC][1];
      ijk_corner[2] = topo.ijk_corner_node_level0[iC][2];
      int node_id = key_manager.getNodeId(ijk_corner, 0);
      if (node_id != -1)
      {
         AOMD::mVertex* v = mmesh.getVertex(node_id);
         if (debug) cout << " classifying a corner node on geometry id " << count << " for node " << node_id << endl;
         v->classify(mmesh.getGEntity(++count, 0));
      }
      else
      {
         ++count;
         cout << "Not found\n";
      }
   }

   int level_max = octree.getLevelMax();
   for (int l = 0; l <= level_max; ++l)
   {
      oOctree::const_iterator it = octree.begin(l);
      oOctree::const_iterator ite = octree.end(l);
      for (; it != ite; ++it)
      {
         if (*it == 1)
         {
            int ijk[3];
            octree.octree2cartesian(it, l, ijk);
            for (int bnd = 0; bnd < topo.nb_entities[dim][dim - 1]; ++bnd)
            {
               if (octree.is_on_boundary(bnd, ijk, l))
               {
                  auto pGeom = mmesh.getGEntity(bnd + 1, dim - 1);
                  if (dim == 2)
                  {
                     int iloc1 = topo.edge_connect[bnd][0];
                     int iloc2 = topo.edge_connect[bnd][1];
                     int unit_offset = 1;
                     int ijk_first[3], ijk_second[3];
                     topo.next_ijk_node_on_element(ijk_first, 0, ijk, unit_offset, dim, iloc1);
                     topo.next_ijk_node_on_element(ijk_second, 0, ijk, unit_offset, dim, iloc2);
                     int node_id_first = key_manager.getNodeId(ijk_first, l);
                     int node_id_second = key_manager.getNodeId(ijk_second, l);
                     mmesh.createEdge(node_id_first, node_id_second, pGeom);
                  }
                  else
                  {  // dim ==3
                     int iloc1 = topo.face_connect[bnd][0];
                     int iloc2 = topo.face_connect[bnd][1];
                     int iloc3 = topo.face_connect[bnd][2];
                     int iloc4 = topo.face_connect[bnd][3];
                     int unit_offset = 1;
                     int ijk_first[3], ijk_second[3], ijk_third[3], ijk_fourth[3];
                     topo.next_ijk_node_on_element(ijk_first, 0, ijk, unit_offset, dim, iloc1);
                     topo.next_ijk_node_on_element(ijk_second, 0, ijk, unit_offset, dim, iloc2);
                     topo.next_ijk_node_on_element(ijk_third, 0, ijk, unit_offset, dim, iloc3);
                     topo.next_ijk_node_on_element(ijk_fourth, 0, ijk, unit_offset, dim, iloc4);
                     int node_id_first = key_manager.getNodeId(ijk_first, l);
                     int node_id_second = key_manager.getNodeId(ijk_second, l);
                     int node_id_third = key_manager.getNodeId(ijk_third, l);
                     int node_id_fourth = key_manager.getNodeId(ijk_fourth, l);

                     if (!simplex)
                     {
                        mmesh.createFaceWithVertices(node_id_first, node_id_second, node_id_third, node_id_fourth, pGeom);
                        if (debug) cout << "-----Face Created------\n";
                     }
                     if (simplex)
                     {
                        if (!alternate)
                        {
                           std::map<int, int> mapIlocId;
                           mapIlocId[iloc1] = node_id_first;
                           mapIlocId[iloc2] = node_id_second;
                           mapIlocId[iloc3] = node_id_third;
                           mapIlocId[iloc4] = node_id_fourth;

                           // 6 cas qui peuvent etre discrimines a partir du prod et de la somme des iloc
                           int somme = iloc1 + iloc2 + iloc3 + iloc4;
                           int produit = iloc1 * iloc2 * iloc3 * iloc4;

                           if (somme == 10)
                           {
                              mmesh.createFaceWithVertices(mapIlocId[0], mapIlocId[1], mapIlocId[5], pGeom);
                              mmesh.createFaceWithVertices(mapIlocId[0], mapIlocId[4], mapIlocId[5], pGeom);
                           }

                           if (somme == 14 && produit == 60)
                           {
                              mmesh.createFaceWithVertices(mapIlocId[1], mapIlocId[2], mapIlocId[5], pGeom);
                              mmesh.createFaceWithVertices(mapIlocId[2], mapIlocId[5], mapIlocId[6], pGeom);
                           }

                           if (somme == 18)
                           {
                              mmesh.createFaceWithVertices(mapIlocId[2], mapIlocId[3], mapIlocId[6], pGeom);
                              mmesh.createFaceWithVertices(mapIlocId[3], mapIlocId[7], mapIlocId[6], pGeom);
                           }

                           if (somme == 14 && produit == 0)
                           {
                              mmesh.createFaceWithVertices(mapIlocId[0], mapIlocId[3], mapIlocId[4], pGeom);
                              mmesh.createFaceWithVertices(mapIlocId[3], mapIlocId[4], mapIlocId[7], pGeom);
                           }

                           if (somme == 6)
                           {
                              mmesh.createFaceWithVertices(mapIlocId[0], mapIlocId[1], mapIlocId[3], pGeom);
                              mmesh.createFaceWithVertices(mapIlocId[1], mapIlocId[2], mapIlocId[3], pGeom);
                           }

                           if (somme == 22)
                           {
                              mmesh.createFaceWithVertices(mapIlocId[4], mapIlocId[5], mapIlocId[7], pGeom);
                              mmesh.createFaceWithVertices(mapIlocId[5], mapIlocId[6], mapIlocId[7], pGeom);
                           }
                        }
                        else
                        {  // alternate == true
                           AOMD::mVertex* v1 = mmesh.getVertex(node_id_first);
                           AOMD::mVertex* v2 = mmesh.getVertex(node_id_second);
                           AOMD::mVertex* v3 = mmesh.getVertex(node_id_third);
                           AOMD::mVertex* v4 = mmesh.getVertex(node_id_fourth);
                           AOMD::mEntity* ed1 = mmesh.getEdge(v1, v3);
                           AOMD::mEntity* ed2 = mmesh.getEdge(v2, v4);

                           if (ed1)
                           {
                              //              cout<<"cas 1\n";
                              mmesh.createFaceWithVertices(node_id_first, node_id_second, node_id_third, pGeom);
                              mmesh.createFaceWithVertices(node_id_first, node_id_fourth, node_id_third, pGeom);
                           }

                           if (ed2)
                           {
                              //              cout<<"cas 2\n";
                              mmesh.createFaceWithVertices(node_id_first, node_id_second, node_id_fourth, pGeom);
                              mmesh.createFaceWithVertices(node_id_second, node_id_fourth, node_id_third, pGeom);
                           }

                        }  // endif alternate
                     }     // quad

                  }  // dim ==3
               }     // is_on_boundary bnd
            }        // bnd
         }
      }
   }  // it

   classifyUnclassifiedVerices(&mmesh);
   // mesh.modifyAllState();
   mmesh.modifyState(3, 2, true);
   mmesh.modifyState(3, 1, true);
   mmesh.modifyState(3, 0, true);
   mmesh.modifyState(2, 1, true);
   mmesh.modifyState(2, 0, true);
   mmesh.modifyState(1, 0, true);
   mmesh.modifyState(0, 1, true);
   mmesh.modifyState(0, 2, true);
   mmesh.modifyState(0, 3, true);
   mmesh.modifyState(1, 2, true);
   mmesh.modifyState(1, 3, true);
   mmesh.modifyState(2, 3, true);

   // telling the mesh what is the level max.
   for (int l = level_max; l >= 0; --l)
   {
      if (octree.size(l))
      {
         mesh.setOctreeLevelMax(l);
         break;
      }
   }

   // setting tag to the hanging nodes.
   for (oKeyManager::const_iterator it = key_manager.begin(); it != key_manager.end(); ++it)
   {
      const oKey* key = *it;
      if (!key->isRegular())
      {
         AOMD::mVertex* v = mmesh.getVertex(key->getId());
         const std::vector<oKey*>& ties = key->getTies();
         if ((dim == 2) || ((dim == 3) && (ties.size() == 2)))
         {
            // init groupes
            std::vector<AOMD::mEntity*> groupe_down;
            groupe_down.reserve(2);
            std::vector<AOMD::mEntity*> groupe_border;
            groupe_border.reserve(2);

            // first upper node
            AOMD::mVertex* v1 = mmesh.getVertex(ties[0]->getId());
            groupe_border.push_back(v1);

            // second upper node
            AOMD::mVertex* v2 = mmesh.getVertex(ties[1]->getId());
            groupe_border.push_back(v2);

            // upper edge
            AOMD::mEntity* edgu = mmesh.getEdge(v1, v2);
            if (!edgu)
            {
               cout << "bugggg wrong edge definition";
               assert(0);
            }
            // attache
            isHangingOn.setData(*v) = edgu;
            // v->attachEntity(id_is_hanging_on, edgu);
            bndGroup.setData(*v) = groupe_border;
            // attachEntitiesVector(v,id_bnd_groupe, groupe_border);

            // first down edge
            AOMD::mEntity* edg = mmesh.getEdge(v1, v);
            if (!edg) edg = mmesh.getEdge(v, v1);
            groupe_down.push_back(edg);
            // attache
            isHangingOn.setData(*edg) = edgu;
            // edg->attachEntity(id_is_hanging_on, edgu);

            // second down edge
            edg = mmesh.getEdge(v, v2);
            if (!edg) edg = mmesh.getEdge(v2, v);
            groupe_down.push_back(edg);
            // attache
            isHangingOn.setData(*edg) = edgu;
            // edg->attachEntity(id_is_hanging_on, edgu);

            // attache
            downGroup.setData(*v) = groupe_down;
            // attachEntitiesVector(v,id_down_groupe, groupe_down);

            // only for 2D attache the hanging entity to the upper entity
            // in 3D this is not used so it is not generated
            // This information is attached to the upper entity to minimise
            // attachment, i.e it could have been attached to both down edges
            if (dim == 2)
            {
               isHangingBy.setData(*v) = edgu;
               // edgu->attachEntity(id_is_hanging_by, v);
            }

            if (debug)
            {
               cout << "creating a hanging node on edge info" << endl;
               cout << " node is " << endl;
               v->print();
               cout << " edge is " << endl;
               edgu->print();
               cout << "Last  hanging edge info" << endl;
               cout << " small edge is " << endl;
               edg->print();
               cout << " big   edge is " << endl;
               edgu->print();
            }
         }
         else if (((dim == 3) && (ties.size() == 4)))
         {
            // init groupes
            std::vector<AOMD::mEntity*> groupe_down;
            groupe_down.reserve(8);
            std::vector<AOMD::mEntity*> groupe_border;
            groupe_border.reserve(16);
            std::vector<AOMD::mEntity*> groupe_up;
            groupe_up.reserve(2);
            std::vector<std::vector<AOMD::mEntity*>> nodes_per_elem(2);

            // set of potentiel face of groupe down
            std::set<AOMD::mEntity*> potential_face;

            // loop creation of corner nodes
            for (size_t j = 0; j < 4; ++j)
            {
               AOMD::mVertex* vj = mmesh.getVertex(ties[j]->getId());
               groupe_border.push_back(vj);
               potential_face.insert(vj->begin(2), vj->end(2));
            }

            // loop creation of border hanging nodes
            for (size_t j = 0; j < 4; ++j)
            {
               // intermediate node hi in betwen 2 corner nodes
               int ijk_hi[3];
               const int* ijk_v1 = ties[j]->getIJK();
               const int* ijk_v2 = ties[(j + 1) % 4]->getIJK();
               for (int k = 0; k < 3; ++k)
               {
                  int dk = ijk_v2[k] - ijk_v1[k];
                  if (dk > 0)
                     ijk_hi[k] = ijk_v1[k] + dk / 2;
                  else if (dk < 0)
                     ijk_hi[k] = ijk_v2[k] - dk / 2;
                  else
                     ijk_hi[k] = ijk_v1[k];
               }
               oKey k_p(ijk_hi);
               AOMD::mEntity* v_hi = mmesh.getVertex(k_p.getId());
               groupe_border.push_back(v_hi);
               potential_face.insert(v_hi->begin(2), v_hi->end(2));
            }

#if QUAD
            AOMD::mEntity* fa =
                mesh.getQuad(static_cast<AOMD::mVertex*>(groupe_border[0]), static_cast<mVertex*>(groupe_border[1]),
                             static_cast<mVertex*>(groupe_border[2]), static_cast<AOMD::mVertex*>(groupe_border[3]));
            groupe_up.push_back(fa);
            potential_face.erase(potential_face.find(fa));
#else
            AOMD::mEntity* fa =
                mmesh.getTri(static_cast<AOMD::mVertex*>(groupe_border[0]), static_cast<AOMD::mVertex*>(groupe_border[1]),
                             static_cast<AOMD::mVertex*>(groupe_border[2]));

            if (!fa)
            {
               fa = mmesh.getTri(static_cast<AOMD::mVertex*>(groupe_border[0]), static_cast<AOMD::mVertex*>(groupe_border[1]),
                                 static_cast<AOMD::mVertex*>(groupe_border[3]));
               groupe_up.push_back(fa);
               potential_face.erase(potential_face.find(fa));
               nodes_per_elem[0].push_back(groupe_border[0]);
               nodes_per_elem[0].push_back(groupe_border[1]);
               nodes_per_elem[0].push_back(groupe_border[3]);
               nodes_per_elem[0].push_back(groupe_border[4]);
               nodes_per_elem[0].push_back(groupe_border[7]);
               nodes_per_elem[0].push_back(v);

               fa = mmesh.getTri(static_cast<AOMD::mVertex*>(groupe_border[3]), static_cast<AOMD::mVertex*>(groupe_border[1]),
                                 static_cast<AOMD::mVertex*>(groupe_border[2]));
               groupe_up.push_back(fa);
               potential_face.erase(potential_face.find(fa));
               nodes_per_elem[1].push_back(groupe_border[3]);
               nodes_per_elem[1].push_back(groupe_border[1]);
               nodes_per_elem[1].push_back(groupe_border[2]);
               nodes_per_elem[1].push_back(groupe_border[5]);
               nodes_per_elem[1].push_back(groupe_border[6]);
               nodes_per_elem[1].push_back(v);
            }
            else
            {
               groupe_up.push_back(fa);
               potential_face.erase(potential_face.find(fa));
               nodes_per_elem[0].push_back(groupe_border[0]);
               nodes_per_elem[0].push_back(groupe_border[1]);
               nodes_per_elem[0].push_back(groupe_border[2]);
               nodes_per_elem[0].push_back(groupe_border[4]);
               nodes_per_elem[0].push_back(groupe_border[5]);
               nodes_per_elem[0].push_back(v);

               fa = mmesh.getTri(static_cast<AOMD::mVertex*>(groupe_border[0]), static_cast<AOMD::mVertex*>(groupe_border[2]),
                                 static_cast<AOMD::mVertex*>(groupe_border[3]));
               groupe_up.push_back(fa);
               potential_face.erase(potential_face.find(fa));
               nodes_per_elem[1].push_back(groupe_border[0]);
               nodes_per_elem[1].push_back(groupe_border[2]);
               nodes_per_elem[1].push_back(groupe_border[3]);
               nodes_per_elem[1].push_back(groupe_border[6]);
               nodes_per_elem[1].push_back(groupe_border[7]);
               nodes_per_elem[1].push_back(v);
            }

            if (!fa)
            {
               cout << "no face found\n";
               assert(0);
            }

#endif
            // attache
            isHangingOn.setData(*v) = fa;
            // v->attachEntity(id_is_hanging_on, fa);

            if (debug)
            {
               cout << "creating a hanging node on face info" << endl;
               cout << " node is " << endl;
               v->print();
               cout << " face is " << endl;
               fa->print();
            }

            // potential faces down iterator
            std::set<AOMD::mEntity*>::iterator itf = potential_face.begin();
            std::set<AOMD::mEntity*>::iterator itfe = potential_face.end();

            // nodes of down level
            std::vector<AOMD::mEntity*> down_nodes(groupe_border);
            down_nodes.push_back(v);
            std::vector<AOMD::mEntity*>::iterator itb = down_nodes.begin();
            std::vector<AOMD::mEntity*>::iterator itbe = down_nodes.end();
            std::vector<AOMD::mEntity*>::iterator itbf;
#ifndef QUAD
            // nodes of down level for first upper element
            std::vector<AOMD::mEntity*>::iterator itb1 = nodes_per_elem[0].begin();
            std::vector<AOMD::mEntity*>::iterator itb1e = nodes_per_elem[0].end();
#endif

            // loop on all potential face to find down faces and down edge border
            for (; itf != itfe; ++itf)
            {
               // face pointer
               AOMD::mEntity* fa = (*itf);

               // node of face
               AOMD::mAdjacencyContainer::iter itn = fa->begin(0);
               AOMD::mAdjacencyContainer::iter itnb = itn;
               AOMD::mAdjacencyContainer::iter itne = fa->end(0);
               short state[3] = {0, 0, 0};
               bool skeep_face = false;

               // does all these nodes are in node border
               for (; itn != itne; ++itn)
               {
                  if ((itbf = std::find(itb, itbe, (*itn))) != itbe)
                  {
                     short d = short(itbf - itb);
                     if (d < 4)
                        state[itn - itnb] = 1;
                     else if (d < 8)
                        state[itn - itnb] = 2;
                  }
                  else
                  {
                     skeep_face = true;
                     break;
                  }
               }

               if (skeep_face) continue;

               groupe_down.push_back(fa);

               // loop on edges to find out if on border
               for (int i = 0; i < 3; ++i)
               {
                  int i1 = (i + 1) % 3;
                  if (state[i] + state[i1] == 3) groupe_border.push_back(fa->get(1, i));
               }

#if QUAD
               isHangingOn.setData(*fa) = groupe_up[0];
               // fa->attachEntity(id_is_hanging_on, groupe_up[0]);
#else
               int upper_is = 0;
               // find upper element
               for (itn = itnb; itn != itne; ++itn)
               {
                  if ((itbf = std::find(itb1, itb1e, (*itn))) != itb1e)
                     ;
                  else
                  {
                     upper_is = 1;
                     break;
                  }
               }
               isHangingOn.setData(*fa) = groupe_up[size_t(upper_is)];
               // fa->attachEntity(id_is_hanging_on, groupe_up[upper_is]);
#endif
            }

            // attache
            downGroup.setData(*v) = groupe_down;
            bndGroup.setData(*v) = groupe_border;
            for (auto pe : groupe_up) isHangingBy.setData(*pe) = v;
            /*
            attachEntitiesVector(v,id_down_groupe, groupe_down);
            attachEntitiesVector(v,id_bnd_groupe, groupe_border);
                    std::vector < mEntity * >::iterator itgu = groupe_up.begin();
            // I think there is a bug !!!
            // original :
            // std::vector < mEntity * >::iterator itgue = groupe_up.begin();
            // corrected to :
            std::vector < mEntity * >::iterator itgue = groupe_up.end();
                    for (; itgu != itgue; ++itgu) ( *itgu )->attachEntity(id_is_hanging_by, v);
            */
         }
         else
            assert(0);
      }
   }

   if (debug) cout << " EXPORTING THE MESH  " << endl;
   if (debug) AOMD::AOMD_Util::Instance()->ex_port("mesh.msh", &mmesh);
   if (debug) cout << " done EXPORTING THE MESH  " << endl;

   ls.setSupport(&mesh, 0.0);
   for (oKeyManager::const_iterator it = key_manager.begin(); it != key_manager.end(); ++it)
   {
      const oKey* key = *it;
      AOMD::mVertex* v = mmesh.getVertex(key->getId());
      ls(v) = lsoct.getVal(key);
   }
}

template <template <class> class DATAMANAGER>
void InterfaceOctreeToAOMDNoHanging(const oOctree& octree, const oField& lsoct, xfem::xMesh& mesh, xfem::xLevelSet& ls,
                                    const bool coarser, const bool simplex, iClassifyCriteria* ptr_clfCriteria,
                                    DATAMANAGER<int>& octreeLevel)
{
   const bool debug = false;
   const oMapping& mapping = octree.getMapping();
   const oTopo& topo(octree.getTopo());
   const int dim = octree.getDim();
   const int* Periodicity = octree.getPeriodicity();
   const oKeyManager& key_manager = lsoct.getKeyManager();
   const int level_max = octree.getLevelMax();
   //  const unsigned int vertex_corner_tag_face =
   //     AOMD::AOMD_Util::Instance()->lookupMeshDataId("OctreeToAOMDNoHanging_vertex_corner_tag_face");
   //  const unsigned int vertex_corner_tag_vol =
   //      AOMD::AOMD_Util::Instance()->lookupMeshDataId("OctreeToAOMDNoHanging_vertex_corner_tag_vol");
   xinterface::aomd::xAttachedDataManagerAOMD<std::vector<AOMD::mEntity*>> vertex_corner_face_tagger;
   xinterface::aomd::xAttachedDataManagerAOMD<std::vector<AOMD::mEntity*>> vertex_corner_vol_tagger;
   int ijk[3];
   int ijk_n[3];
   AOMD::mMesh& mmesh = mesh.getMesh();
   for (oKeyManager::const_iterator it = key_manager.begin(); it != key_manager.end(); ++it)
   {
      const oKey* key = *it;
      double xyz[3];
      mapping.ijk2xyz(key->getIJK(), level_max, xyz);
      if (debug)
      {
         const int* ijk = key->getIJK();
         cout << " vertex created with id " << key->getId() << " at ijk_fine " << endl;
         std::copy(ijk, ijk + 3, std::ostream_iterator<int>(std::cout, " "));
         cout << endl;
         cout << " at location " << endl;
         std::copy(xyz, xyz + 3, std::ostream_iterator<double>(std::cout, " "));
         cout << endl;
      }

      int id = key->getId();
      mmesh.createVertex(id, xyz[0], xyz[1], xyz[2], nullptr);
   }

   // init level set with 0 value for curent mesh made only of node
   ls.setSupport(&mesh, 0.0);

   for (int l = 0; l <= level_max; l++)
   {
      oOctree::const_iterator it = octree.begin(l);
      oOctree::const_iterator ite = octree.end(l);
      const int offset = topo.pow_base2[level_max - l];
      int ids[8];
      for (; it != ite; ++it)
      {
         if (*it == 1)
         {
            octree.octree2cartesian(it, l, ijk);
            key_manager.getNodeIdsOnElt(ijk, l, ids);
            oOctree::const_iterator children;
            if (l == octree.getLevelMax())
               children = octree.end();
            else
               children = octree.getChildren(it, l);
            int physical_index = init_physical_index;
            if (nullptr != ptr_clfCriteria)
               physical_index = (*ptr_clfCriteria)(it, l, ijk, children, children + octree.getNbChildren());
            auto pGeom = mmesh.getGEntity(physical_index, octree.getDim());
            // 2D
            if (dim == 2)
            {
               int hanging_on_neighbor_edge[4];
               // look  if hanging is mandatory by inspecting neighbor cell of curent one
               int count = 0;
               int is_hanging = 0;
               int id_hanging = 0;  // to avoid -Wmaybe-uninitialized. Value not important as it is modified when is_hanging is
                                    // incremented and used only if is_hanging>0
               int on_what = 1;
               int neighbor = 0;
               bool exist;
               while (topo.nextNeighbor(ijk, l, ijk_n, on_what, exist, dim, Periodicity, count, neighbor))
               {
                  hanging_on_neighbor_edge[count] = 0;
                  if (exist)
                  {
                     int a = (int)(*octree.cartesian2octree(ijk_n, l));
                     if (a < 1)
                     {
                        is_hanging++;
                        id_hanging = count;
                        hanging_on_neighbor_edge[count] = 1;
                     }
                  }
                  // if boundary generate edge to classify
                  else
                  {
                     int id1 = topo.edge_connect[count][0];
                     int id2 = topo.edge_connect[count][1];
                     mmesh.createEdge(ids[id1], ids[id2], mmesh.getGEntity(count + 1, 1));
                  }

                  count++;
               }

               if ((coarser && is_hanging > 1) || (!coarser && is_hanging > 0))
               {
                  int ijk_c[3];
                  int ijk_o[3];
                  topo.cartesian2finest_cartesian(ijk, l, ijk_o, level_max);
                  topo.ijk_center_on_element(ijk_c, ijk_o, offset, dim);
                  oKey center(ijk_c);
                  double xyz[3];
                  mapping.ijk2xyz(ijk_c, level_max, xyz);
                  int cid = center.getId();
                  AOMD::mVertex* v = mmesh.createVertex(cid, xyz[0], xyz[1], xyz[2], nullptr);
                  AOMD::mEntity* central_vertex = static_cast<AOMD::mEntity*>(v);
                  count = 0;
                  double mean_ls = 0.;
                  while (topo.next_ijk_node_on_element(ijk_n, 0, ijk_o, offset, dim, count))
                  {
                     const oKey* key = key_manager.find(ijk_n);
                     mean_ls += lsoct.getVal(key);
                     count++;
                  }
                  ls(v) = mean_ls / count;

                  count = 0;
                  while (topo.next_ijk_node_on_element(ijk_n, 1, ijk_o, offset, dim, count))
                  {
                     int id1 = ids[count];
                     int id2 = ids[(count + 1) % 4];
                     // if it is hanging edge
                     if (hanging_on_neighbor_edge[count])
                     {
                        oKey hanging(ijk_n);
                        const int hid = hanging.getId();
                        AOMD::mFace* e1 = mmesh.createFaceWithVertices(id1, hid, cid, pGeom);
                        AOMD::mFace* e2 = mmesh.createFaceWithVertices(hid, id2, cid, pGeom);
                        octreeLevel.setData(*e1) = l;
                        octreeLevel.setData(*e2) = l;
                     }
                     // it is regular edge
                     else
                     {
                        AOMD::mFace* e1 = mmesh.createFaceWithVertices(id1, id2, cid, pGeom);
                        octreeLevel.setData(*e1) = l;
                     }

                     // push vertex in entity vector of id1
                     std::vector<AOMD::mEntity*>* vertex_vector;
                     AOMD::mEntity* vcorner = (mmesh.getVertex(id1));
                     //   if ((vertex_vector = xinterface::aomd::getAttachedEntitiesVector(vcorner, vertex_corner_tag_face)))
                     if (vertex_vector = vertex_corner_face_tagger.getData(*vcorner))
                     {
                        vertex_vector->push_back(central_vertex);
                     }
                     else
                     {
                        std::vector<AOMD::mEntity*> new_vertex_vector(1, central_vertex);
                        // xinterface::aomd::attachEntitiesVector(vcorner, vertex_corner_tag_face, new_vertex_vector);
                        vertex_corner_face_tagger.setData(*vcorner) = new_vertex_vector;
                     }
                     count++;
                  }
               }
               else if (coarser && is_hanging > 0)
               {
                  int ijk_o[3];
                  topo.cartesian2finest_cartesian(ijk, l, ijk_o, level_max);
                  topo.next_ijk_node_on_element(ijk_n, 1, ijk_o, offset, dim, id_hanging);
                  oKey hanging(ijk_n);
                  const int hid = hanging.getId();
                  int id0 = ids[id_hanging];
                  int id1 = ids[(id_hanging + 1) % 4];
                  int id2 = ids[(id_hanging + 2) % 4];
                  int id3 = ids[(id_hanging + 3) % 4];
                  AOMD::mFace* e1 = mmesh.createFaceWithVertices(id0, hid, id3, pGeom);
                  AOMD::mFace* e2 = mmesh.createFaceWithVertices(hid, id1, id2, pGeom);
                  AOMD::mFace* e3 = mmesh.createFaceWithVertices(hid, id2, id3, pGeom);
                  octreeLevel.setData(*e1) = l;
                  octreeLevel.setData(*e2) = l;
                  octreeLevel.setData(*e3) = l;
               }
               else
               {
                  if ((ijk[0] + ijk[1]) % 2 == 0)
                  {
                     AOMD::mFace* e1 = mmesh.createFaceWithVertices(ids[0], ids[1], ids[2], pGeom);
                     AOMD::mFace* e2 = mmesh.createFaceWithVertices(ids[0], ids[2], ids[3], pGeom);

                     octreeLevel.setData(*e1) = l;
                     octreeLevel.setData(*e2) = l;
                  }
                  else
                  {
                     AOMD::mFace* e1 = mmesh.createFaceWithVertices(ids[0], ids[1], ids[3], pGeom);
                     AOMD::mFace* e2 = mmesh.createFaceWithVertices(ids[1], ids[2], ids[3], pGeom);
                     octreeLevel.setData(*e1) = l;
                     octreeLevel.setData(*e2) = l;
                  }
               }
            }
            // 3D
            else
            {
               assert(dim == 3);

               // zizag mesh selector
               int ii = ijk[0] + ijk[1] + ijk[2];
               int hanging_on_neighbor_edge[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
               int hanging_on_neighbor_face[6] = {0, 0, 0, 0, 0, 0};
               int id_hanging_edge[6];
               int neighbor = 0;
               bool exist;
               int is_hanging_edge = 0, is_hanging_face = 0, nb_hanging, iedge_cell;
               int iface, iedge;
               AOMD::mTet *tet1, *tet2, *tet3, *tet4, *tet5, *tet6, *tet7, *tet8;

               // loop on edge to set hangings by inspecting neighbor cell of curent one
               iedge = 0;
               while (topo.nextNeighbor(ijk, l, ijk_n, 1, exist, dim, Periodicity, iedge, neighbor))
               {
                  if (exist)
                  {
                     if (*octree.cartesian2octree(ijk_n, l) < 1)
                     {
                        ++is_hanging_edge;
                        hanging_on_neighbor_edge[iedge] = 1;
                        neighbor = 2;
                     }
                  }
                  if (neighbor > 1)
                  {
                     neighbor = 0;
                     iedge++;
                  }
                  else
                     ++neighbor;
               }

               // loop on face to set hangings by inspecting neighbor cell of curent one
               iface = 0;
               while (topo.nextNeighbor(ijk, l, ijk_n, 2, exist, dim, Periodicity, iface, neighbor))
               {
                  if (exist && *octree.cartesian2octree(ijk_n, l) < 1)
                  {
                     ++is_hanging_face;
                     hanging_on_neighbor_face[iface] = -4;
                  }
                  else
                  {
                     nb_hanging = 0;
                     for (iedge = 0; iedge < 4; ++iedge)
                     {
                        if (hanging_on_neighbor_edge[topo.face_edge_connect[iface][iedge]])
                        {
                           id_hanging_edge[iface] = iedge;
                           ++nb_hanging;
                        }
                     }
                     hanging_on_neighbor_face[iface] = nb_hanging;
                  }

                  iface++;
               }

               // First implementation step : no 3D coarser approche (only one hanging face for example => no need of central
               // node). Only use of coarser for face. Maybe in future something will/might be done if hanging edge (hanging_face
               // imply hanging edge so there is no need to check is_hanging_face. this last one is computed for 3D coarser
               // approche in future) => central point
               if (is_hanging_edge > 0)
               {
                  // center node creation
                  int ijk_c[3];
                  int ijk_o[3];
                  topo.cartesian2finest_cartesian(ijk, l, ijk_o, level_max);
                  topo.ijk_center_on_element(ijk_c, ijk_o, offset, dim);
                  oKey center(ijk_c);
                  double xyz[3];
                  mapping.ijk2xyz(ijk_c, level_max, xyz);
                  int cid = center.getId();
                  AOMD::mVertex* v = mmesh.createVertex(cid, xyz[0], xyz[1], xyz[2], nullptr);
                  AOMD::mEntity* central_vertex = static_cast<AOMD::mEntity*>(v);

                  // ls center node evaluation and tag for future projection
                  int count = 0;
                  double mean_ls = 0.;
                  std::vector<AOMD::mEntity*>* vertex_vector;
                  while (topo.next_ijk_node_on_element(ijk_n, 0, ijk_o, offset, dim, count))
                  {
                     // find corner key
                     const oKey* key = key_manager.find(ijk_n);
                     assert(key);
                     int id_corner = key->getId();

                     // accumulate value
                     mean_ls += lsoct.getVal(key);

                     // push central vertex in entity vector of corner's
                     AOMD::mEntity* vcorner = mmesh.getVertex(id_corner);
                     // if ((vertex_vector = xinterface::aomd::getAttachedEntitiesVector(vcorner, vertex_corner_tag_vol)))
                     if (vertex_vector = vertex_corner_vol_tagger.getData(*vcorner))
                     {
                        vertex_vector->push_back(central_vertex);
                     }
                     else
                     {
                        std::vector<AOMD::mEntity*> new_vertex_vector(1, central_vertex);
                        //   xinterface::aomd::attachEntitiesVector(vcorner, vertex_corner_tag_vol, new_vertex_vector);
                        vertex_corner_vol_tagger.setData(*vcorner) = new_vertex_vector;
                     }

                     count++;
                  }
                  ls(v) = mean_ls / count;

                  // loop on face
                  int i0, i1, i2, i3;
                  for (int iface = 0; iface < 6; ++iface)
                  {
                     int nbhf = hanging_on_neighbor_face[iface];
                     bool do_bnd = (octree.is_on_boundary(iface, ijk, l));

                     // if face got hanging pattern
                     if (nbhf > 0)
                     {
                        // if not coarser or coarser with more then one hanging use center point
                        if ((coarser && nbhf > 1) || (!coarser))
                        {
                           // create/retrive vertex on center of face
                           if (!(topo.next_ijk_node_on_element(ijk_n, 2, ijk_o, offset, dim, iface))) throw -8961;
                           oKey fc(ijk_n);
                           int ifc = fc.getId();
                           double xyz[3];
                           mapping.ijk2xyz(ijk_n, level_max, xyz);
                           AOMD::mVertex* v = mmesh.getVertex(ifc);
                           if (!v)
                           {
                              v = mmesh.createVertex(ifc, xyz[0], xyz[1], xyz[2], nullptr);
                              central_vertex = static_cast<AOMD::mEntity*>(v);
                              double mean_ls = 0.;
                              for (int inode = 0; inode < 4; ++inode)
                              {
                                 // find face corner key
                                 if (!(topo.next_ijk_node_on_element(ijk_n, 0, ijk_o, offset, dim,
                                                                     topo.face_connect[iface][inode])))
                                    throw -9657;
                                 const oKey* key = key_manager.find(ijk_n);
                                 assert(key);
                                 int id_corner = key->getId();

                                 // accumulate value
                                 mean_ls += lsoct.getVal(key);

                                 // push central vertex in entity vector of corner's
                                 AOMD::mEntity* vcorner = mmesh.getVertex(id_corner);
                                 //   if ((vertex_vector = xinterface::aomd::getAttachedEntitiesVector(vcorner,
                                 //   vertex_corner_tag_face)))
                                 if (vertex_vector = vertex_corner_face_tagger.getData(*vcorner))
                                 {
                                    vertex_vector->push_back(central_vertex);
                                 }
                                 else
                                 {
                                    std::vector<AOMD::mEntity*> new_vertex_vector(1, central_vertex);
                                    //    xinterface::aomd::attachEntitiesVector(vcorner, vertex_corner_tag_face,
                                    //    new_vertex_vector);
                                    vertex_corner_face_tagger.setData(*vcorner) = new_vertex_vector;
                                 }
                              }
                              ls(v) = mean_ls / 4.;
                           }

                           // genarate element by looping on edge of face and do 2D algo for base face construction
                           for (iedge = 0; iedge < 4; ++iedge)
                           {
                              iedge_cell = topo.face_edge_connect[iface][iedge];

                              // this use the fact that face_edge_connect discribe edges of face of vol in same order as
                              // face_connect discribe node of face. edge_connect is not used as its definition doesn't
                              // give way to respect normal inward volume of genarated tet
                              i0 = ids[topo.face_connect[iface][(0 + iedge) % 4]];
                              i1 = ids[topo.face_connect[iface][(1 + iedge) % 4]];

                              if (hanging_on_neighbor_edge[iedge_cell])
                              {
                                 if (!(topo.next_ijk_node_on_element(ijk_n, 1, ijk_o, offset, dim, iedge_cell))) throw -7698;
                                 oKey ec(ijk_n);
                                 int ie = ec.getId();
                                 tet1 = mmesh.createTetWithVertices(i0, ifc, ie, cid, mmesh.getGEntity(physical_index, 3));
                                 tet2 = mmesh.createTetWithVertices(ie, ifc, i1, cid, mmesh.getGEntity(physical_index, 3));
                                 octreeLevel.setData(*tet2) = l;
                                 if (do_bnd)
                                 {
                                    mmesh.createFaceWithVertices(i0, ie, ifc, mmesh.getGEntity(iface + 1, 2));
                                    mmesh.createFaceWithVertices(ie, i1, ifc, mmesh.getGEntity(iface + 1, 2));
                                 }
                              }
                              else
                              {
                                 tet1 = mmesh.createTetWithVertices(ifc, i1, i0, cid, mmesh.getGEntity(physical_index, 3));
                                 if (do_bnd) mmesh.createFaceWithVertices(ifc, i0, i1, mmesh.getGEntity(iface + 1, 2));
                              }
                              octreeLevel.setData(*tet1) = l;
                           }
                        }
                        // if coarser with only one hanging no extra center point
                        else
                        {
                           int ief = id_hanging_edge[iface];
                           iedge_cell = topo.face_edge_connect[iface][ief];
                           if (!(topo.next_ijk_node_on_element(ijk_n, 1, ijk_o, offset, dim, iedge_cell))) throw -11618;
                           oKey ec(ijk_n);
                           int ie = ec.getId();
                           i0 = ids[topo.face_connect[iface][(0 + ief) % 4]];
                           i1 = ids[topo.face_connect[iface][(1 + ief) % 4]];
                           i2 = ids[topo.face_connect[iface][(2 + ief) % 4]];
                           i3 = ids[topo.face_connect[iface][(3 + ief) % 4]];
                           tet1 = mmesh.createTetWithVertices(i0, i3, ie, cid, mmesh.getGEntity(physical_index, 3));
                           tet2 = mmesh.createTetWithVertices(ie, i2, i1, cid, mmesh.getGEntity(physical_index, 3));
                           tet3 = mmesh.createTetWithVertices(i3, i2, ie, cid, mmesh.getGEntity(physical_index, 3));
                           octreeLevel.setData(*tet1) = l;
                           octreeLevel.setData(*tet2) = l;
                           octreeLevel.setData(*tet3) = l;
                           if (do_bnd)
                           {
                              mmesh.createFaceWithVertices(i0, ie, i3, mmesh.getGEntity(iface + 1, 2));
                              mmesh.createFaceWithVertices(ie, i1, i2, mmesh.getGEntity(iface + 1, 2));
                              mmesh.createFaceWithVertices(i3, ie, i2, mmesh.getGEntity(iface + 1, 2));
                           }
                        }
                     }
                     // if face is a hanging face
                     else if (nbhf < 0)
                     {
                        assert(!do_bnd);
                        if (!(topo.next_ijk_node_on_element(ijk_n, 2, ijk_o, offset, dim, iface))) throw -8745;
                        oKey fc(ijk_n);
                        int ifc = fc.getId();
                        i0 = ids[topo.face_connect[iface][0]];
                        i1 = ids[topo.face_connect[iface][1]];
                        i2 = ids[topo.face_connect[iface][2]];
                        i3 = ids[topo.face_connect[iface][3]];
                        if (!(topo.next_ijk_node_on_element(ijk_n, 1, ijk_o, offset, dim, topo.face_edge_connect[iface][0])))
                           throw -8548;
                        oKey ec1(ijk_n);
                        int ie1 = ec1.getId();
                        if (!(topo.next_ijk_node_on_element(ijk_n, 1, ijk_o, offset, dim, topo.face_edge_connect[iface][1])))
                           throw -8548;
                        oKey ec2(ijk_n);
                        int ie2 = ec2.getId();
                        if (!(topo.next_ijk_node_on_element(ijk_n, 1, ijk_o, offset, dim, topo.face_edge_connect[iface][2])))
                           throw -8548;
                        oKey ec3(ijk_n);
                        int ie3 = ec3.getId();
                        if (!(topo.next_ijk_node_on_element(ijk_n, 1, ijk_o, offset, dim, topo.face_edge_connect[iface][3])))
                           throw -8548;
                        oKey ec4(ijk_n);
                        int ie4 = ec4.getId();
                        // zigzag mesh for coarser mesh
                        if (coarser)
                        {
                           tet1 = mmesh.createTetWithVertices(i0, ie4, ie1, cid, mmesh.getGEntity(physical_index, 3));
                           tet2 = mmesh.createTetWithVertices(i1, ie1, ie2, cid, mmesh.getGEntity(physical_index, 3));
                           tet3 = mmesh.createTetWithVertices(i2, ie2, ie3, cid, mmesh.getGEntity(physical_index, 3));
                           tet4 = mmesh.createTetWithVertices(i3, ie3, ie4, cid, mmesh.getGEntity(physical_index, 3));
                           tet5 = mmesh.createTetWithVertices(ifc, ie1, ie4, cid, mmesh.getGEntity(physical_index, 3));
                           tet6 = mmesh.createTetWithVertices(ifc, ie2, ie1, cid, mmesh.getGEntity(physical_index, 3));
                           tet7 = mmesh.createTetWithVertices(ifc, ie3, ie2, cid, mmesh.getGEntity(physical_index, 3));
                           tet8 = mmesh.createTetWithVertices(ifc, ie4, ie3, cid, mmesh.getGEntity(physical_index, 3));
                        }
                        // regular mesh
                        else
                        {
                           switch (iface)
                           {
                              case 0:
                              case 3:
                              case 4:
                              {
                                 tet1 = mmesh.createTetWithVertices(ifc, ie1, i0, cid, mmesh.getGEntity(physical_index, 3));
                                 tet2 = mmesh.createTetWithVertices(ifc, ie2, ie1, cid, mmesh.getGEntity(physical_index, 3));
                                 tet3 = mmesh.createTetWithVertices(ie1, ie2, i1, cid, mmesh.getGEntity(physical_index, 3));
                                 tet4 = mmesh.createTetWithVertices(i2, ie2, ifc, cid, mmesh.getGEntity(physical_index, 3));
                                 tet5 = mmesh.createTetWithVertices(i2, ifc, ie3, cid, mmesh.getGEntity(physical_index, 3));
                                 tet6 = mmesh.createTetWithVertices(i3, ie3, ie4, cid, mmesh.getGEntity(physical_index, 3));
                                 tet7 = mmesh.createTetWithVertices(ie3, ifc, ie4, cid, mmesh.getGEntity(physical_index, 3));
                                 tet8 = mmesh.createTetWithVertices(i0, ie4, ifc, cid, mmesh.getGEntity(physical_index, 3));
                                 break;
                              }
                              default:
                              {
                                 tet1 = mmesh.createTetWithVertices(i0, ie4, ie1, cid, mmesh.getGEntity(physical_index, 3));
                                 tet2 = mmesh.createTetWithVertices(i1, ie1, ifc, cid, mmesh.getGEntity(physical_index, 3));
                                 tet3 = mmesh.createTetWithVertices(i1, ifc, ie2, cid, mmesh.getGEntity(physical_index, 3));
                                 tet4 = mmesh.createTetWithVertices(i2, ie2, ie3, cid, mmesh.getGEntity(physical_index, 3));
                                 tet5 = mmesh.createTetWithVertices(ie2, ifc, ie3, cid, mmesh.getGEntity(physical_index, 3));
                                 tet6 = mmesh.createTetWithVertices(i3, ie3, ifc, cid, mmesh.getGEntity(physical_index, 3));
                                 tet7 = mmesh.createTetWithVertices(i3, ifc, ie4, cid, mmesh.getGEntity(physical_index, 3));
                                 tet8 = mmesh.createTetWithVertices(ie1, ie4, ifc, cid, mmesh.getGEntity(physical_index, 3));
                                 break;
                              }
                           }
                        }

                        octreeLevel.setData(*tet1) = l;
                        octreeLevel.setData(*tet2) = l;
                        octreeLevel.setData(*tet3) = l;
                        octreeLevel.setData(*tet4) = l;
                        octreeLevel.setData(*tet5) = l;
                        octreeLevel.setData(*tet6) = l;
                        octreeLevel.setData(*tet7) = l;
                        octreeLevel.setData(*tet8) = l;
                     }
                     // if face is normal (no hanging, current standard mesh for this level)
                     else
                     {
                        // zigzag mesh for coarser mesh
                        if (coarser)
                        {
                           if (ii % 2 == 0)
                           {
                              // Cas A
                              i0 = 0;
                              i1 = 1;
                              i2 = 2;
                              i3 = 3;
                           }
                           else
                           {
                              // Cas B
                              i0 = 1;
                              i1 = 2;
                              i2 = 3;
                              i3 = 0;
                           }
                           if (iface == 4)
                           {
                              ++i0;
                              ++i1;
                              ++i2;
                              ++i3;
                           }
                        }
                        // regular mesh
                        else
                        {
                           if (iface < 2)
                           {
                              i0 = 1;
                              i1 = 2;
                              i2 = 3;
                              i3 = 0;
                           }
                           else if (iface != 4)
                           {
                              i0 = 0;
                              i1 = 1;
                              i2 = 2;
                              i3 = 3;
                           }
                           else
                           {
                              i0 = 1;
                              i1 = 2;
                              i2 = 3;
                              i3 = 4;
                           }
                        }
                        i0 = ids[topo.face_connect[iface][(i0 + iface) % 4]];
                        i1 = ids[topo.face_connect[iface][(i1 + iface) % 4]];
                        i2 = ids[topo.face_connect[iface][(i2 + iface) % 4]];
                        i3 = ids[topo.face_connect[iface][(i3 + iface) % 4]];
                        tet1 = mmesh.createTetWithVertices(i0, i3, i1, cid, mmesh.getGEntity(physical_index, 3));
                        tet2 = mmesh.createTetWithVertices(i3, i2, i1, cid, mmesh.getGEntity(physical_index, 3));
                        octreeLevel.setData(*tet1) = l;
                        octreeLevel.setData(*tet2) = l;

                        if (do_bnd)
                        {
                           mmesh.createFaceWithVertices(i0, i1, i3, mmesh.getGEntity(iface + 1, 2));
                           mmesh.createFaceWithVertices(i3, i1, i2, mmesh.getGEntity(iface + 1, 2));
                        }
                     }
                  }
               }
               // no hanging => std mesh
               else
               {
                  int i0, i1, i2, i3;
                  int ii0, ii1, ii2, ii3;

                  // zigzag mesh for coarser mesh
                  if (coarser)
                  {
                     // ZigZag mesh
                     if (ii % 2 == 0)
                     {
                        // Cas A
                        i0 = 0;
                        i1 = 1;
                        i2 = 2;
                        i3 = 3;
                        tet1 = mmesh.createTetWithVertices(ids[0], ids[1], ids[3], ids[4], mmesh.getGEntity(physical_index, 3));
                        tet2 = mmesh.createTetWithVertices(ids[5], ids[1], ids[4], ids[6], mmesh.getGEntity(physical_index, 3));
                        tet3 = mmesh.createTetWithVertices(ids[2], ids[3], ids[1], ids[6], mmesh.getGEntity(physical_index, 3));
                        tet4 = mmesh.createTetWithVertices(ids[4], ids[6], ids[1], ids[3], mmesh.getGEntity(physical_index, 3));
                        tet5 = mmesh.createTetWithVertices(ids[7], ids[6], ids[4], ids[3], mmesh.getGEntity(physical_index, 3));
                     }
                     else
                     {
                        // Cas B
                        i0 = 1;
                        i1 = 2;
                        i2 = 3;
                        i3 = 0;
                        tet1 = mmesh.createTetWithVertices(ids[1], ids[2], ids[0], ids[5], mmesh.getGEntity(physical_index, 3));
                        tet2 = mmesh.createTetWithVertices(ids[4], ids[5], ids[0], ids[7], mmesh.getGEntity(physical_index, 3));
                        tet3 = mmesh.createTetWithVertices(ids[0], ids[5], ids[2], ids[7], mmesh.getGEntity(physical_index, 3));
                        tet4 = mmesh.createTetWithVertices(ids[6], ids[2], ids[5], ids[7], mmesh.getGEntity(physical_index, 3));
                        tet5 = mmesh.createTetWithVertices(ids[3], ids[7], ids[0], ids[2], mmesh.getGEntity(physical_index, 3));
                     }
                  }
                  // regular mesh
                  else
                  {
                     // regular mesh
                     tet1 = mmesh.createTetWithVertices(ids[1], ids[3], ids[0], ids[5], mmesh.getGEntity(physical_index, 3));
                     tet2 = mmesh.createTetWithVertices(ids[1], ids[2], ids[3], ids[5], mmesh.getGEntity(physical_index, 3));
                     tet3 = mmesh.createTetWithVertices(ids[6], ids[3], ids[2], ids[5], mmesh.getGEntity(physical_index, 3));
                     tet4 = mmesh.createTetWithVertices(ids[4], ids[0], ids[3], ids[5], mmesh.getGEntity(physical_index, 3));
                     tet5 = mmesh.createTetWithVertices(ids[4], ids[3], ids[7], ids[5], mmesh.getGEntity(physical_index, 3));
                     tet6 = mmesh.createTetWithVertices(ids[7], ids[6], ids[5], ids[3], mmesh.getGEntity(physical_index, 3));
                     octreeLevel.setData(*tet6) = l;
                  }
                  octreeLevel.setData(*tet1) = l;
                  octreeLevel.setData(*tet2) = l;
                  octreeLevel.setData(*tet3) = l;
                  octreeLevel.setData(*tet4) = l;
                  octreeLevel.setData(*tet5) = l;

                  // boundary
                  for (int iface = 0; iface < 6; ++iface)
                  {
                     if (octree.is_on_boundary(iface, ijk, l))
                     {
                        if (coarser)
                        {
                           if (iface == 4)
                           {
                              ii0 = i0 + 1;
                              ii1 = i1 + 1;
                              ii2 = i2 + 1;
                              ii3 = i3 + 1;
                           }
                           else
                           {
                              ii0 = i0;
                              ii1 = i1;
                              ii2 = i2;
                              ii3 = i3;
                           }
                        }
                        else
                        {
                           if (iface < 2)
                           {
                              ii0 = 1;
                              ii1 = 2;
                              ii2 = 3;
                              ii3 = 0;
                           }
                           else if (iface != 4)
                           {
                              ii0 = 0;
                              ii1 = 1;
                              ii2 = 2;
                              ii3 = 3;
                           }
                           else
                           {
                              ii0 = 1;
                              ii1 = 2;
                              ii2 = 3;
                              ii3 = 4;
                           }
                        }
                        ii0 = ids[topo.face_connect[iface][(ii0 + iface) % 4]];
                        ii1 = ids[topo.face_connect[iface][(ii1 + iface) % 4]];
                        ii2 = ids[topo.face_connect[iface][(ii2 + iface) % 4]];
                        ii3 = ids[topo.face_connect[iface][(ii3 + iface) % 4]];
                        mmesh.createFaceWithVertices(ii0, ii1, ii3, mmesh.getGEntity(iface + 1, 2));
                        mmesh.createFaceWithVertices(ii3, ii1, ii2, mmesh.getGEntity(iface + 1, 2));
                     }
                  }
               }
            }
         }
      }
   }

   // classify corners
   int ijk_corner[3];
   int count = 0;
   for (int iC = 0; iC < topo.pow_base2[dim]; ++iC)
   {
      ijk_corner[0] = topo.ijk_corner_node_level0[iC][0];
      ijk_corner[1] = topo.ijk_corner_node_level0[iC][1];
      ijk_corner[2] = topo.ijk_corner_node_level0[iC][2];
      int node_id = key_manager.getNodeId(ijk_corner, 0);
      if (node_id != -1)
      {
         AOMD::mVertex* v = mmesh.getVertex(node_id);
         if (debug) cout << " classifying a corner node on geometry id " << count << " for node " << node_id << endl;
         v->classify(mmesh.getGEntity(++count, 0));
      }
      else
      {
         ++count;
         cout << "Not found\n";
      }
   }

   classifyUnclassifiedVerices(&mmesh);
   // mesh.modifyAllState();
   mmesh.modifyState(3, 2, true);
   mmesh.modifyState(3, 1, true);
   mmesh.modifyState(3, 0, true);
   mmesh.modifyState(2, 1, true);
   mmesh.modifyState(2, 0, true);
   mmesh.modifyState(1, 0, true);
   mmesh.modifyState(0, 1, true);
   mmesh.modifyState(0, 2, true);
   mmesh.modifyState(0, 3, true);
   mmesh.modifyState(1, 2, true);
   mmesh.modifyState(1, 3, true);
   mmesh.modifyState(2, 3, true);

   // telling the mesh what is the level max.
   for (int l = level_max; l >= 0; --l)
   {
      if (octree.size(l))
      {
         mesh.setOctreeLevelMax(l);
         break;
      }
   }

   // change support of the ls to take into account all elements newly created
   // this doesn't change allready seted value on nodes and on extra central nodes
   xfem::xRegion mesh_region(&mesh);
   ls.reduceSupport(mesh_region);

   if (debug) cout << " EXPORTING THE MESH  " << endl;
   if (debug) AOMD::AOMD_Util::Instance()->ex_port("mesh.msh", &mmesh);
   if (debug) cout << " done EXPORTING THE MESH  " << endl;

   for (const oKey* key : key_manager)
   {
      AOMD::mVertex* v = mmesh.getVertex(key->getId());
      ls(v) = lsoct.getVal(key);
   }
}

}  // namespace xoctree

}  // end namespace xinterface
