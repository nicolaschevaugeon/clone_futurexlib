/*
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
// std
#include <algorithm>
// xinterface/xoctree
#include "AdaptOctreeToAOMD.h"
// aomd
#include "mEdge.h"
#include "mFace.h"
#include "mHex.h"
#include "mTet.h"
// xoctree
#include "oExport.h"
#include "oKeyManager.h"
#include "oMapping.h"
#include "oOctree.h"
// xfem
#include "xLevelSet.h"
#include "xMesh.h"

// main

// xz/normale -y -> surf 1
// xz/normale +y -> surf 3
// xy/normale -z -> surf 5
// xy/normale +z -> surf 6
// yz/normale -x -> surf 4
// yz/normale +x -> surf 2

namespace xinterface
{
namespace xoctree
{
inline void attachCell(AOMD::mEntity* Entity, unsigned int tagNum, const oOctree::cell_type* cell)
{
   xAttachableCell* ac = (xAttachableCell*)Entity->getData(tagNum);
   if (!ac)
   {
      ac = new xAttachableCell;
      Entity->attachData(tagNum, ac);
   }
   ac->cell = cell;
}

template <template <class> class DATAMANAGER>
void InterfaceOctreeToAOMDG(const oOctree& octree, const oField& lsoct, oLevelSet& o_ls, xfem::xMesh& mesh, xfem::xLevelSet& ls,
                            DATAMANAGER<int>& octreeLevel, DATAMANAGER<AOMD::mEntity*>& isHangingOn, const bool fullModify,
                            const bool simplex, bool hgg)
{
   const bool debug = false;
   AOMD::mMesh& mmesh = mesh.getMesh();
   cout << "Creation du maillage";
   if (hgg)
   {
      cout << " avec noeuds hangings.\n";
   }
   else
   {
      cout << " sans noeud hanging.\n";
   }

   const ::xoctree::oMapping& mapping = octree.getMapping();
   const ::xoctree::oTopo& topo(octree.getTopo());
   const oKeyManager& key_manager = lsoct.getKeyManager();
   const int dim = octree.getDim();

   const unsigned int creatorCell_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("creator_cell");

   cout << "ccel tag=" << creatorCell_tag << endl;

   double xyz[3];
   for (const ::xoctree::oKey* key : key_manager)
   {
      mapping.ijk2xyz(key->getIJK(), octree.getLevelMax(), xyz);
      if (debug)
      {
         const int* ijk = key->getIJK();
         cout << " vertex created with id " << key->getId() << " at ijk_fine " << endl;
         std::copy(ijk, ijk + 3, std::ostream_iterator<int>(std::cout, " "));
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

            if (dim == 2)
            {
               // quad case
               if (!simplex)
               {
                  AOMD::mFace* Quad =
                      mmesh.createFaceWithVertices(ids[0], ids[1], ids[2], ids[3], mmesh.getGEntity(100, octree.getDim()));
                  octreeLevel.setData(*Quad) = l;
                  attachCell(Quad, creatorCell_tag, it);
               }
               else
               {
                  // simplex case
                  if ((ijk[0] + ijk[1]) % 2 == 0)
                  {
                     AOMD::mFace* Tri1 = mmesh.createFaceWithVertices(ids[0], ids[1], ids[2], mmesh.getGEntity(100, dim));
                     AOMD::mFace* Tri2 = mmesh.createFaceWithVertices(ids[0], ids[2], ids[3], mmesh.getGEntity(100, dim));
                     // On associe la cellule et son niveau aux elements...
                     octreeLevel.setData(*Tri1) = l;
                     octreeLevel.setData(*Tri2) = l;
                     attachCell(Tri1, creatorCell_tag, it);
                     attachCell(Tri2, creatorCell_tag, it);
                  }
                  else
                  {
                     AOMD::mFace* Tri1 = mmesh.createFaceWithVertices(ids[0], ids[1], ids[3], mmesh.getGEntity(100, dim));
                     AOMD::mFace* Tri2 = mmesh.createFaceWithVertices(ids[1], ids[2], ids[3], mmesh.getGEntity(100, dim));
                     // On associe la cellule et son niveau aux elements...
                     octreeLevel.setData(*Tri1) = l;
                     octreeLevel.setData(*Tri2) = l;
                     attachCell(Tri1, creatorCell_tag, it);
                     attachCell(Tri2, creatorCell_tag, it);
                  }
               }
            }
            else
            {  // dim3
               if (!simplex)
               {
                  // hex
                  AOMD::mHex* Hex = mmesh.createHexWithVertices(ids[0], ids[1], ids[2], ids[3], ids[4], ids[5], ids[6], ids[7],
                                                                mmesh.getGEntity(100, dim));
                  octreeLevel.setData(*Hex) = l;
                  attachCell(Hex, creatorCell_tag, it);
               }
               else
               {
                  // Maillage non alterne
                  AOMD::mTet* tet = nullptr;
                  tet = mmesh.createTetWithVertices(ids[0], ids[1], ids[5], ids[3], mmesh.getGEntity(100, 3));
                  attachCell(tet, creatorCell_tag, it);
                  octreeLevel.setData(*tet) = l;
                  tet = mmesh.createTetWithVertices(ids[1], ids[2], ids[3], ids[5], mmesh.getGEntity(100, 3));
                  attachCell(tet, creatorCell_tag, it);
                  octreeLevel.setData(*tet) = l;
                  tet = mmesh.createTetWithVertices(ids[2], ids[3], ids[5], ids[6], mmesh.getGEntity(100, 3));
                  attachCell(tet, creatorCell_tag, it);
                  octreeLevel.setData(*tet) = l;
                  tet = mmesh.createTetWithVertices(ids[0], ids[3], ids[4], ids[5], mmesh.getGEntity(100, 3));
                  attachCell(tet, creatorCell_tag, it);
                  octreeLevel.setData(*tet) = l;
                  tet = mmesh.createTetWithVertices(ids[3], ids[4], ids[5], ids[7], mmesh.getGEntity(100, 3));
                  attachCell(tet, creatorCell_tag, it);
                  octreeLevel.setData(*tet) = l;
                  tet = mmesh.createTetWithVertices(ids[3], ids[5], ids[6], ids[7], mmesh.getGEntity(100, 3));
                  attachCell(tet, creatorCell_tag, it);
                  octreeLevel.setData(*tet) = l;
               }
            }
         }
      }
   }

   // classify corners
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
         v->classify(mmesh.getGEntity(++count, 0));
         if (debug) std::cout << " classifying a corner node on geometry id " << count << " for node " << node_id << endl;
      }
      else
      {
         ++count;
         cout << "Not found\n";
      }
   }

   // Otherwise, cannot detect edges and create boundary faces
   if (dim == 3)
   {
      mmesh.modifyState(1, 3);
      mmesh.modifyState(3, 1);
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
                     mmesh.createEdge(node_id_first, node_id_second, mmesh.getGEntity(bnd + 1, dim - 1));
                  }
                  else
                  {
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
                     //			  mesh.createFaceWithVertices(node_id_first, node_id_second,
                     //						      node_id_third, node_id_fourth,
                     //						      mesh.getGEntity(bnd+1,dim-1));

                     //                mesh.modifyState(3,1,true);//Si pas commente, performance d'un escargot !
                     AOMD::mVertex* v1 = mmesh.getVertex(node_id_first);
                     AOMD::mVertex* v2 = mmesh.getVertex(node_id_second);
                     AOMD::mVertex* v3 = mmesh.getVertex(node_id_third);
                     AOMD::mVertex* v4 = mmesh.getVertex(node_id_fourth);

                     if (simplex)
                     {
                        AOMD::mEdge* ed1 = mmesh.getEdge(v1, v3);
                        AOMD::mEdge* ed2 = mmesh.getEdge(v2, v4);
                        if (ed1)
                        {
                           mmesh.createFaceWithVertices(node_id_first, node_id_second, node_id_third,
                                                        mmesh.getGEntity(bnd + 1, dim - 1));
                           mmesh.createFaceWithVertices(node_id_first, node_id_fourth, node_id_third,
                                                        mmesh.getGEntity(bnd + 1, dim - 1));
                        }
                        if (ed2)
                        {
                           mmesh.createFaceWithVertices(node_id_first, node_id_second, node_id_fourth,
                                                        mmesh.getGEntity(bnd + 1, dim - 1));
                           mmesh.createFaceWithVertices(node_id_second, node_id_fourth, node_id_third,
                                                        mmesh.getGEntity(bnd + 1, dim - 1));
                        }
                     }
                     else
                     {
                        mmesh.createFaceWithVertices(node_id_first, node_id_second, node_id_third, node_id_fourth,
                                                     mmesh.getGEntity(bnd + 1, dim - 1));
                        if (debug)
                           std::cout << " classifying a face node on geometry id " << bnd + 1 << " for node " << node_id_first
                                     << " " << node_id_second << " " << node_id_third << " " << node_id_fourth << endl;
                     }
                  }
               }
            }
         }
      }
   }

   classifyUnclassifiedVerices(&mmesh);
   if (fullModify)
   {
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
   }
   else
   {
      cout << "Modify 3->1\n";
      mmesh.modifyState(3, 1);
      cout << "Modify 3->2\n";
      mmesh.modifyState(3, 2);
      cout << "Modify 1->0\n";
      mmesh.modifyState(1, 0);
      cout << "Modify 2->0\n";
      mmesh.modifyState(2, 0);
      cout << "Modify Done\n";
   }

   // No Hanging :
   // liste des cellules octrees qui ont un noeud orphelin au milieu d'une des aretes :
   map<::xoctree::oOctree::const_iterator, vector<AOMD::mVertex*>> crs_cells_hgg;
   // Niveau octree des cellules actives
   map<::xoctree::oOctree::const_iterator, int> active_cells_lvl;

   // setting tag to the hanging nodes.
   for (::xoctree::oKey* key : key_manager)
   {
      if (key->isHanging())
      {
         AOMD::mEntity* v = mmesh.getVertex(key->getId());
         if (dim == 2)
         {
            const std::vector<::xoctree::oKey*>& ties = key->getTies();
            AOMD::mVertex* v1 = mmesh.getVertex(ties[0]->getId());
            AOMD::mVertex* v2 = mmesh.getVertex(ties[1]->getId());
            AOMD::mEdge* ed = mmesh.getEdge(v1, v2);
            isHangingOn.setData(*v) = ed;
            // v->attachEntity(xMesh::get_is_hanging_on_tag(), ed);
            if (debug)
            {
               cout << "creating a hanging node info" << endl;
               cout << " node is " << endl;
               v->print();
               cout << " edge is " << endl;
               ed->print();
            }

            // listage des noeuds hgg dans la map cree plus haut
            if (!hgg)  // Si on ne veut pas qu'il y ait de noeud hanging
            {
               AOMD::mEntity* f1 = ed->get(2, 0);  // on recupere l'element (peut etre pas utile)
               xAttachableCell* ac = (xAttachableCell*)f1->getData(creatorCell_tag);  // On recupere l'attachable cell

               oOctree::const_iterator cell =
                   (oOctree::const_iterator)ac->cell;  // On recupere la cellule a laquelle l'element est lie

               if (crs_cells_hgg.find(cell) == crs_cells_hgg.end())  // Si pas dans la map
               {
                  vector<AOMD::mVertex*> vec;
                  crs_cells_hgg[cell] = vec;  // on ajoute la cle
               }
               crs_cells_hgg[cell].push_back(static_cast<AOMD::mVertex*>(v));  // On rajoute le hangnig au vecteur
               // active_cells_lvl[cell]= f1->getAttachedInt(octreeLevel_tag);//On recupere le niveau de la cellule
               auto plevel = octreeLevel.getData(*f1);
               active_cells_lvl[cell] = (plevel) ? *plevel : 0;
            }
         }
         else
         {
            const std::vector<::xoctree::oKey*>& ties = key->getTies();
            int nbTies = int(ties.size());
            AOMD::mEntity* ed = nullptr;

            AOMD::mVertex* v1 = mmesh.getVertex(ties[0]->getId());
            AOMD::mVertex* v2 = mmesh.getVertex(ties[1]->getId());
            if (nbTies == 2) ed = mmesh.getEdge(v1, v2);

            // cout << " node is " << endl; v->print();
            // cout<<"ID1="<<ties[0]->getId()<<" ID2="<<ties[1]->getId()<<endl;

            if (nbTies != 2)
            {
               AOMD::mVertex* v3 = mmesh.getVertex(ties[2]->getId());
               AOMD::mVertex* v4 = mmesh.getVertex(ties[3]->getId());

               if (simplex)
               {
                  ed = mmesh.getEdge(v1, v3);
                  if (!ed) ed = mmesh.getEdge(v2, v4);
                  if (!ed) throw;
               }
               else
               {
                  //                        v1->print();
                  //                        v2->print();
                  //                        v3->print();
                  //                        v4->print();
                  ed = mmesh.getQuad(v1, v2, v3, v4);
               }
            }

            isHangingOn.setData(*v) = ed;
            // v->attachEntity(xMesh::get_is_hanging_on_tag(), ed);
            if (debug)
            {
               cout << "creating a hanging node info" << endl;
               cout << " node is " << endl;
               v->print();
               cout << " entity is " << endl;
               ed->print();
            }
         }
      }
   }

   ls.setSupport(&mesh, 0.0);

   // maintenant qu'on a fait la liste, on traite les cellules concernees
   if (!hgg)  // Si on ne veut pas qu'il y ait de noeud hanging
   {
      for (map<oOctree::const_iterator, std::vector<AOMD::mVertex*>>::iterator iter = crs_cells_hgg.begin();
           iter != crs_cells_hgg.end(); iter++)  // pour chaque cellule qui a au moins 1 hgg
      {
         // Recuperation des 4 noeuds de la cellule *******************
         oOctree::const_iterator cell = iter->first;  // On recupere la cellule concernee
         int level = active_cells_lvl[cell];          // On recupere le niveau de la cellule
         int ijk[3], ids[8];                          // ijk = coordonne ijk de la cellule et ids=id des noeuds de la cellule
         octree.octree2cartesian(cell, level, ijk);   // on recupere ijk
         /*active_nodes*/ key_manager.getNodeIdsOnElt(ijk, level, ids);  // on recupere les ids des noeuds
         AOMD::mVertex *v0 = mmesh.getVertex(ids[0]), *v1 = mmesh.getVertex(ids[1]), *v2 = mmesh.getVertex(ids[2]),
                       *v3 = mmesh.getVertex(ids[3]);  // Les noeuds (du maillage)

         // Suppression des faces du maillage *************************
         pGEntity classification_face;    // Classifiaction des elements
         if ((ijk[0] + ijk[1]) % 2 == 0)  // On recupere les faces 1 et 2 qu'il faut supprimer en fonction de leur orientation
         {
            classification_face = (mmesh.getTri(v0, v1, v2))->getClassification();  // On recupere la classification de la face 1
            mmesh.DEL_updateAdj(mmesh.getTri(v0, v1, v2));                          // Suppression face 1
            mmesh.DEL_updateAdj(mmesh.getTri(v0, v2, v3));                          // Suppression face 2
            mmesh.DEL_updateAdj(mmesh.getEdge(v0, v2));                             // Suppression de l'edge diagonale
            // les autres edges inutiles seront enlevees plus tard
         }
         else
         {
            classification_face = (mmesh.getTri(v0, v1, v3))->getClassification();  // On recupere la classification de la face 1
            mmesh.DEL_updateAdj(mmesh.getTri(v0, v1, v3));                          // Supression face 1
            mmesh.DEL_updateAdj(mmesh.getTri(v1, v2, v3));                          // Suppresion face 2
            mmesh.DEL_updateAdj(mmesh.getEdge(v1, v3));                             // Suppression de l'edge diagonale
            // les autres edges inutiles seront enlevees plus tard
         }

         // Rajout du point du milieu ************************
         Trellis_Util::mPoint Pm = (v0->point() + v2->point()) * 0.5;                   // coordonnees du point milieu
         AOMD::mVertex* p_milieu = mmesh.createVertex(Pm, nullptr /*classification*/);  // il faut creer le vertex du milieu...
         // id_AN++;

         // Affectation de la valeur de la LS au noeud milieu *****************
         int ijk_m_lsup[3];    // tableau position octree noeud milieu niveau level+1
         int ijk_m_finest[3];  // idem au plus grand niveau (le plus fin)
         int un[3] = {1, 1, 0};
         transform(ijk, ijk + 3, ijk_m_lsup, bind2nd(std::multiplies<int>(), 2));
         transform(ijk_m_lsup, ijk_m_lsup + 3, un, ijk_m_lsup, std::plus<int>());  // Position noeud milieu level+1
         // octree.cartesian2finest_cartesian(ijk_m_lsup,level+1,ijk_m_finest);//Position noeud milieu level_max ////existe pas...
         // octree.cartesian2octree( ijk_m_lsup,level+1)
         octree.octree2finest_cartesian(octree.cartesian2octree(ijk_m_lsup, level + 1), level + 1, ijk_m_finest);

         // ls(p_milieu) = lsoct.getLevelSet(ijk_m_finest);//On effecte la valeur de la levelset
         //      oKey* key_milieu=(oKey*) key_manager.find(ijk_m_finest);//Clef du noeud du milieu
         //      ls(p_milieu) = lsoct.getVal(key_milieu);//On effecte la valeur de la levelset
         ls(p_milieu) = o_ls.getLevelSet(ijk_m_finest);

         // Reperage des cote avec noeuds hgg et sans noeud hgg ************************************
         map<AOMD::mEntity*, AOMD::mVertex*> edge_hgg;  // map des cotes qui ont un noeud hanging.
         for (int i = 0; i <= static_cast<int>(crs_cells_hgg[cell].size()) - 1;
              i++)  // On tague les cote qui ont des noeuds hanging
         {
            AOMD::mEntity* n_hgg = (crs_cells_hgg[cell][i]);
            // mEntity *edge=n_hgg->getAttachedEntity(xMesh::get_is_hanging_on_tag());
            AOMD::mEntity** ppedge = isHangingOn.getData(*n_hgg);
            AOMD::mEntity* edge = (ppedge) ? *ppedge : nullptr;
            edge_hgg[edge] = crs_cells_hgg[cell][i];  // on ajoute le noeud a la face
         }

         // Ajout des nouvelles edges (en fonction de noeuds hgg ou pas) *************************
         int Idm = p_milieu->getId();  // ID noeud milieu
         for (int i = 0; i <= 3; i++)  // Pour chaque cote de la cellule
         {
            int Id1 = ids[i];            // ID noeud 1
            int Id2 = ids[(i + 1) % 4];  // ID noeud 2

            AOMD::mEntity* edge = mmesh.getEdge(mmesh.getVertex(Id1), mmesh.getVertex(Id2));  // On recupere l'edge

            if (edge_hgg.find(edge) == edge_hgg.end())  // si le cote n'a pas un noeud hanging (cf au dessus)
            {
               AOMD::mEntity* Tri = mmesh.createFaceWithVertices(
                   Id1, Id2, Idm, classification_face /*mesh.getGEntity(100, octree.getDim())*/);  // On cree un element tout seul
               attachCell(Tri, creatorCell_tag, cell);                                             // Ca marche ??
            }
            else  // si le cote a un noeud hanging
            {
               int Idh = edge_hgg[edge]->getId();  // ID du noeud hhg de l'edge
               AOMD::mEntity* Tri1 = mmesh.createFaceWithVertices(
                   Id1, Idh, Idm,
                   classification_face /*mesh.getGEntity(100, octree.getDim())*/);  // On cree les 2 petits elements
               AOMD::mEntity* Tri2 =
                   mmesh.createFaceWithVertices(Id2, Idm, Idh, classification_face /*mesh.getGEntity(100, octree.getDim())*/);
               attachCell(Tri1, creatorCell_tag, cell);  // On y attache la cellule octree d'origine (ca peut etre utile)
               attachCell(Tri2, creatorCell_tag, cell);
               isHangingOn.deleteData(*mmesh.getVertex(Idh));
               //(mesh.getVertex(Idh))->deleteData(xMesh::get_is_hanging_on_tag());//efface le tag hanging
               mmesh.DEL_updateAdj(edge);  // On efface l'edge qui ne sert plus a rien
            }
         }

      }  // fin du for

      classifyUnclassifiedVerices(&mmesh);
      //      mesh.modifyAllState();

   }  // fin du if(!hgg)

   if (debug) cout << " EXPORTING THE MESH  " << endl;
   AOMD::AOMD_Util::Instance()->ex_port("mesh.msh", &mmesh);
   if (debug) cout << " done EXPORTING THE MESH  " << endl;

   for (oKeyManager::const_iterator it = key_manager.begin(); it != key_manager.end(); ++it)
   {
      ::xoctree::oKey* key = *it;
      AOMD::mVertex* v = mmesh.getVertex(key->getId());
      ls(v) = lsoct.getVal(key);
      //      cout<<"ijk="<<key->getIJK()[0]<<" "<<key->getIJK()[1]<<endl;
   }
}

// Grégory a recoder dans le nouveau cadre, désolé ...

template <template <class> class DATAMANAGER>
void buildOneLevelFinerMesh(oOctree& octree, const oField& lsoct_field, oLevelSet& o_lset, oKeyManager& key_manager,
                            xfem::xMesh& meshFine, xfem::xLevelSet& lsFine, DATAMANAGER<int>& octreeLevel,
                            DATAMANAGER<AOMD::mEntity*>& isHangingOn, const bool fullModify, const bool simplex, const bool hgg,
                            int nb)
{
   cout << "Build Finer Mesh\n";
   const bool debug = false;
   ::xoctree::oRefinementCriteriaTrue oneLevel;

   ::xoctree::oObserverKeyManager obsKeyM(key_manager);
   ::xoctree::oObserverRecorder obsRec;
   octree.attach(obsKeyM);
   octree.attach(obsRec);
   for (int ii = 0; ii < nb; ++ii)
   {
      octree.refine(oneLevel);
   }
   InterfaceOctreeToAOMDG(octree, lsoct_field, o_lset, meshFine, lsFine, octreeLevel, isHangingOn, fullModify, simplex, hgg);
   if (debug) AOMD::AOMD_Util::Instance()->ex_port("meshFine.msh", &meshFine.getMesh());
   // Maintenant qu'on a fait ce qu'on avait a faire, on revient a l'etat initial
   octree.modifyOctree(obsRec);
   octree.detach(obsRec);
   octree.detach(obsKeyM);
}

}  // namespace xoctree

}  // namespace xinterface

// ---------------------------------------
