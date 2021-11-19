/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include "oKeyManager.h"

#include <fstream>
#include <functional>

#include "oTopo.h"

using std::cout;
using std::endl;
using std::fill;

namespace xoctree
{
int oKey::i_field = 0;
double oKey::step_[3] = {0., 0., 0.};
int oKey::coeff_[3] = {0, 0, 0};

void oKeyManager::getKeysOnElt(const int* ijk_ori, int level, std::vector<oKey*>& elt_keys)
{
   int ijk_ori_fine[3], ijk[3];
   topo.cartesian2finest_cartesian(ijk_ori, level, ijk_ori_fine, level_max);
   const int offset = topo.pow_base2[level_max - level];
   int count = 0;
   while (topo.next_ijk_node_on_element(ijk, 0, ijk_ori_fine, offset, dim, count))
   {
      oKey* key = find(ijk);
      elt_keys.push_back(key);
      count++;
   }

   return;
}
void oKeyManager::getKeysOnElt(const int* ijk_ori, int level, std::vector<const oKey*>& elt_keys) const
{
   int ijk_ori_fine[3], ijk[3];
   topo.cartesian2finest_cartesian(ijk_ori, level, ijk_ori_fine, level_max);
   const int offset = topo.pow_base2[level_max - level];
   int count = 0;
   while (topo.next_ijk_node_on_element(ijk, 0, ijk_ori_fine, offset, dim, count))
   {
      const oKey* key = find(ijk);
      elt_keys.push_back(key);
      count++;
   }

   return;
}

int oKeyManager::getNodeId(int* ijk, int level) const
{
   int ijk_fine[3];
   topo.cartesian2finest_cartesian(ijk, level, ijk_fine, level_max);
   const oKey* key = find(ijk_fine);
   // assert(key);
   if (nullptr == key)
      return -1;
   else
      return key->getId();
}

void oKeyManager::getNodeIdsOnElt(int* ijk_ori, int level, int* node_ids) const
{
   int ijk_ori_fine[3], ijk[3];
   topo.cartesian2finest_cartesian(ijk_ori, level, ijk_ori_fine, level_max);
   const int offset = topo.pow_base2[level_max - level];
   int count = 0;
   while (topo.next_ijk_node_on_element(ijk, 0, ijk_ori_fine, offset, dim, count))
   {
      const oKey* key = find(ijk);
      node_ids[count] = key->getId();
      count++;
   }

   return;
}

void oKeyManager::createStencils()
{
   // we first reserve the proper memory and initialize with no neighbors
   // four lines below could be written as one.
   std::for_each(begin(), end(), bind2nd(std::mem_fun(&oKey::resizeStencil), dim * 2));
   // iterator it =  begin();
   // iterator ite = end();
   // const int d = octree.dim * 2;
   // for (;it != ite;++it) (*it)->resizeStencil(d);

   // we now assemble the stencils
   for (int l = level_max; l >= 0; --l)
   {
      int ijk[3];
      fill(ijk, ijk + 3, 0);
      do
      {
         oOctree::iterator cell = octree.cartesian2octree(ijk, l);
         if (*cell == 1)  // we consider only active cells
         {
            std::vector<oKey*> elt_keys;
            getKeysOnElt(ijk, l, elt_keys);
            for (int i = 0; i < topo.nb_entities[dim][0]; ++i)
            {
               oKey* key = elt_keys[i];
               if (key->isRegular())
               {
                  for (int j = 0; j < topo.nb_nodes_connected_to_node[dim]; ++j)
                  {
                     key->setStencil(topo.stencil_location[dim][i][j], elt_keys[topo.connected_nodes[dim][i][j]]);
                  }
               }
               if (key->isPeriodic())
               {
                  // on assemble sur le copain
                  oKey* keyo = key->getMirror();
                  for (int j = 0; j < topo.nb_nodes_connected_to_node[dim]; ++j)
                  {
                     keyo->setStencil(topo.stencil_location[dim][i][j], elt_keys[topo.connected_nodes[dim][i][j]]);
                  }
               }
            }
         }
      } while (octree.next_ijk(ijk, l));
   }
}

void oKeyManager::createKeys()
{
   const bool debug = false;
   keys.clear();

   // declare active nodes = union of the nodes of the active elements.
   for (int l = 0; l <= level_max; l++)
   {
      oOctree::const_iterator it = octree.begin(l);
      oOctree::const_iterator ite = octree.end(l);
      const int offset = topo.pow_base2[level_max - l];
      int ijk_ori[3];
      int ijk[3];
      for (; it != ite; ++it)
      {
         if (*it == 1)
         {
            octree.octree2finest_cartesian(it, l, ijk_ori);
            int count = 0;
            while (topo.next_ijk_node_on_element(ijk, 0, ijk_ori, offset, dim, count++))
            {
               insert(ijk);
            }
         }
      }
   }
   if (debug) cout << " size  active_nodes after decl " << size() << endl;

   // we now separate the active nodes into two categories
   //   free nodes and
   //   tied nodes
   //   both free nodes and tied nodes have their own id.
   //   this id is stored in id[0]
   //   other ids may be stored after id[0] for tied nodes.
   //   - in case of periodicity, id[1] gives the related node
   //   - in 2D, a hanging node is related to 2 nodes given by
   //            id[0], id[1]
   //   - in 3D, a node hanging on an edge  is related to 2 nodes given by
   //            id[0], id[1],
   //   - in 3D, a node hanging on a face is related to 4 nodes given by
   //            id[0], id[1], id[2], id[3]
   //   The size of the id vector is given by nb_ids

   for (int l = 0; l <= level_max; l++)
   {
      oOctree::const_iterator it = octree.begin(l);
      oOctree::const_iterator ite = octree.end(l);
      const int offset = topo.pow_base2[level_max - l];
      int ijk_ori[3];
      int ijk[3];
      for (; it != ite; ++it)
      {
         if (*it == 1)
         {
            octree.octree2finest_cartesian(it, l, ijk_ori);
            // loop over the nodes
            int count = 0;
            if (debug) cout << " starting loop on nodes " << endl;
            std::vector<oKey*> elt_keys;
            while (topo.next_ijk_node_on_element(ijk, 0, ijk_ori, offset, dim, count))
            {
               oKey* key = find(ijk);  // create above with same algo so it exist => no check
               elt_keys.push_back(key);
               count++;
            }

            // loop over the edges
            if (l != level_max)
            {
               int count = 0;
               if (debug) cout << " starting loop on edges " << endl;
               while (topo.next_ijk_node_on_element(ijk, 1, ijk_ori, offset, dim, count))
               {
                  if (oKey* key = find(ijk))
                  {
                     key->setHanging(elt_keys[topo.edge_connect[count][0]], elt_keys[topo.edge_connect[count][1]]);

                     // periodique a coder
                  }
                  count++;
               }
            }

            // loop over the faces
            if (l != level_max && dim == 3)
            {
               int count = 0;
               if (debug) cout << " starting loop on faces " << endl;
               while (topo.next_ijk_node_on_element(ijk, 2, ijk_ori, offset, dim, count))
               {
                  if (oKey* key = find(ijk))
                  {
                     key->setHanging(elt_keys[topo.face_connect[count][0]], elt_keys[topo.face_connect[count][1]],
                                     elt_keys[topo.face_connect[count][2]], elt_keys[topo.face_connect[count][3]]);

                     // periodique a coder
                  }
                  count++;
               }
            }
         }
      }
   }

   return;
}

oObserverKeyManager::oObserverKeyManager(oKeyManager& k) : key_manager(k) {}
void oObserverKeyManager::refine(const oOctree& octree, oOctree::const_iterator cell, int level, const int* ijk)
{
   key_manager.refine(cell, level, ijk);
}
void oObserverKeyManager::derefine(const oOctree& octree, oOctree::const_iterator cell, int level, const int* ijk)
{
   key_manager.derefine(cell, level, ijk);
}

void oKeyManager::derefine(oOctree::const_iterator cell, int level, const int* ijk)
{
   const bool debug = false;
   int ijk_fine[3];
   int ijk_ori_fine[3];
   topo.cartesian2finest_cartesian(ijk, level, ijk_ori_fine, level_max);

   if (debug)
   {
      cout << "in  oLevelKeyManager::derefine ijk_ori_fine is ";
      std::copy(ijk_ori_fine, ijk_ori_fine + 3, std::ostream_iterator<int>(std::cout, " "));
      cout << " level " << level << endl;
   }

   const int offset = topo.pow_base2[level_max - level];
   topo.ijk_center_on_element(ijk_fine, ijk_ori_fine, offset, dim);
   erase(ijk_fine);

   int exist_on_neighbor_edge[12];
   int exist_on_neighbor_face[6];
   std::vector<oKey*> keys_elt;
   int edge_iterator{0};

   if (dim == 2)
   {
      octree.lookForNodeOnNeighbors2D(cell, level, ijk, exist_on_neighbor_edge);
      edge_iterator = 4;
   }
   else if (dim == 3)
   {
      octree.lookForNodeOnNeighbors3D(cell, level, ijk, exist_on_neighbor_edge, exist_on_neighbor_face);
      edge_iterator = 12;
   }

   // this is available for 2D and partially in 3D.
   // if in list on edge exist we will need to transform the value as hanging
   if (*std::max_element(exist_on_neighbor_edge, exist_on_neighbor_edge + edge_iterator) == 1)
   {
      getKeysOnElt(ijk, level, keys_elt);
   }

   // loop on edges
   int count = 0;
   while (topo.next_ijk_node_on_element(ijk_fine, 1, ijk_ori_fine, offset, dim, count))
   {
      oKey* key = find(ijk_fine);
      if (key)
      {
         // if (!exist_on_neighbor_edge[count] || isOnBBoxEdges(key))  erase(ijk_fine);
         // isOnBBoxEdges no more needed as now lookForNodeOnNeighbors.. are doing the job correctely
         // leaved here for history
         if (!exist_on_neighbor_edge[count])
            erase(ijk_fine);
         else  // the values become a hanging
         {
            key->setHanging(keys_elt[topo.edge_connect[count][0]], keys_elt[topo.edge_connect[count][1]]);
         }
      }
      count++;
   }

   if (dim == 3)
   {
      // check
      if (keys_elt.size() == 0 && *std::max_element(exist_on_neighbor_face, exist_on_neighbor_face + 6) == 1)
      {
         getKeysOnElt(ijk, level, keys_elt);
      }

      // loop on face
      int count = 0;
      while (topo.next_ijk_node_on_element(ijk_fine, 2, ijk_ori_fine, offset, dim, count))
      {
         oKey* key = find(ijk_fine);
         if (key)
         {
            if (!exist_on_neighbor_face[count])
               erase(ijk_fine);
            else  // the values become a hanging
            {
               key->setHanging(keys_elt[topo.face_connect[count][0]], keys_elt[topo.face_connect[count][1]],
                               keys_elt[topo.face_connect[count][2]], keys_elt[topo.face_connect[count][3]]);
            }
         }
         count++;
      }
   }
}

void oKeyManager::refine(oOctree::const_iterator cell, int level, const int* ijk)
{
   const bool debug = false;

   std::vector<oKey*> keys_elt;
   getKeysOnElt(ijk, level, keys_elt);

   int ijk_ori_fine[3];
   topo.cartesian2finest_cartesian(ijk, level, ijk_ori_fine, level_max);
   if (debug)
   {
      cout << "in  oLevelKeyManager::refine ijk_ori_fine is ";
      std::copy(ijk_ori_fine, ijk_ori_fine + 3, std::ostream_iterator<int>(std::cout, " "));
      cout << " level " << level << endl;
   }

   // double vc = std::accumulate(keys_elt.begin(), keys_elt.end(), 0.)/(double) keys_elt.size();
   // a titre de debug vérifier que la valeur est bine insérée par le retour bool
   int ijk_fine[3];
   const int offset = topo.pow_base2[level_max - level];
   topo.ijk_center_on_element(ijk_fine, ijk_ori_fine, offset, dim);  // Centre du carre (2D) ou du CUBE (3D)

   std::pair<oKey*, bool> check = insert(ijk_fine);
   oKey* keyc = check.first;
   assert(check.second);

   std::vector<oKey*>::const_iterator it = keys_elt.begin();
   std::vector<oKey*>::const_iterator ite = keys_elt.end();
   std::vector<double> valskey(nb_field);
   for (; it != ite; ++it)
   {
      const oKey* key = *it;
      key->getVals(valskey);
      std::transform(valskey.begin(), valskey.end(), keyc->beginVal(), keyc->beginVal(), std::plus<double>());
   }
   std::transform(keyc->beginVal(), keyc->endVal(), keyc->beginVal(), bind2nd(std::divides<double>(), (double)keys_elt.size()));

   if (debug)
   {
      cout << " creating  center value at ijk_fine ";
      std::copy(ijk_fine, ijk_fine + 3, std::ostream_iterator<int>(std::cout, " "));
      cout << " of value " << endl;
      std::copy(keyc->beginVal(), keyc->endVal(), std::ostream_iterator<double>(std::cout, " "));
      cout << endl;
   }

   if (dim == 2)
   {  // 2D CASE
#if 0
        /* this have been modified on december 2011 /January 2012 by A.S.
         * this work fine but to be consistant with 3D go back to use of octree instead of 
         * keymanager
         * see #else below
         */
      int count = 0;
      while(topo.next_ijk_node_on_element(ijk_fine, 1, ijk_ori_fine, offset, dim, count))
      {
        //je pense que l'on a plus besoin des exist_on_neighbor
        //on cherche simplement si la clef ijk_fine existe déjà ou pas.
        //  si elle existe déjà on la met classique
        //  si elle n'existe pas on l'introduit comme hanging
        //on crée une valeur hanging
        std::pair<oKey*, bool> check =  insert(ijk_fine);
        oKey* key = check.first;


        if (check.second)
        {
            key->setHanging(keys_elt[topo.edge_connect[count][0]],
            keys_elt[topo.edge_connect[count][1]]);
        }

        if(!check.second || (isOnBoundary(key) && !(key->isPeriodic())))
        {
            key->setRegular();
        }
            count++;
      }
#else
      int hanging_on_neighbor_edge[4];
      // look  if hanging is mandatory by inspecting neighbor cell of curent one
      octree.lookForNodeHangingOnNeighbors2D(cell, level, ijk, hanging_on_neighbor_edge);
      // loop on edges
      int count = 0;
      while (topo.next_ijk_node_on_element(ijk_fine, 1, ijk_ori_fine, offset, dim, count))
      {
         // find/create key
         std::pair<oKey*, bool> check = insert(ijk_fine);
         oKey* key = check.first;

         // if new key inserted
         // it is regular by default
         if (check.second)
         {
            // if it is hanging
            if (hanging_on_neighbor_edge[count])
            {
               key->setHanging(keys_elt[topo.edge_connect[count][0]], keys_elt[topo.edge_connect[count][1]]);
            }
            // it is regular
            // values have to be set as it's a new key
            else
            {
               for (int ik = 0; ik < 2; ++ik)
               {
                  oKey* key_ik = keys_elt[topo.edge_connect[count][ik]];
                  key_ik->getVals(valskey);
                  std::transform(valskey.begin(), valskey.end(), key->beginVal(), key->beginVal(), std::plus<double>());
               }
               std::transform(key->beginVal(), key->endVal(), key->beginVal(), bind2nd(std::divides<double>(), (double)2));
            }
         }
         // if already created key
         else
         {
            // if it is not hanging anymore set as regular
            if (!hanging_on_neighbor_edge[count]) key->setRegular();
         }
         if (debug) cout << "edge count " << count << key;
         count++;
      }
#endif
   }
   else
   {  // 3D CASE
#if 0
        /* this have been modified on december 2011 /January 2012 by A.S.
         * after adding the creation of values in new bounding keys there is still
         * a problem 
         * Go back to lookForNodeOnNeighbors3D and exist_on_neighbor wich where debuged with derefine methode
         * see #else below
         */
        // ATTENTION  : experimental !!! (surtout pour les hanging edges...)
        // modifier pour prendre en compte la valeur pour les nouvelle clé
        // semble y avoir un bug qcq par quand même. Re

        //Create new nodes on edges
        int count = 0;
        while(topo.next_ijk_node_on_element(ijk_fine, 1, ijk_ori_fine, offset, dim, count))
        {
            std::pair<oKey*, bool> check =  insert(ijk_fine);
            oKey* key = check.first;


            // if new key inserted
            // it is regular by default
            if (check.second )
            {
                // if on Bounding edge it is regular
                // but values have to be set as it's a new key
                if (isOnBBoxEdges(key))
                {
                    for (int ik=0; ik<2; ++ik)
                    {
                        oKey* key_ik = keys_elt[topo.edge_connect[count][ik]];
                        key_ik->getVals(valskey);
                        std::transform(valskey.begin(), valskey.end(), key->beginVal(),  key->beginVal(), std::plus<double>());
                    }
                    std::transform(key->beginVal(), key->endVal(), key->beginVal(), 
                            bind2nd(std::divides<double>(), (double)2 ));
                }
                // otherwise it is hanging
                else
                {
                    key->setHanging(keys_elt[topo.edge_connect[count][0]],
                                keys_elt[topo.edge_connect[count][1]]);
                }
            }
            // otherwise it's a already created key
            else
            {
                int nbreg=key->countSetRegular();
                // set key as regular and fix values from ties at the same time if :
                // * isOnBoundary is true as this key is then only connected to 2 cell, the one created here
                //   and the one which create the key
                // * if nbreg > 2 , i.e nbreg==3 , 4 cell are connected to this key which is the maximum possible
                //   this key is not hanging anymore
                if(nbreg >2 || isOnBoundary(key))
                { 
                    key->setRegular();
                }
            }

            count++;
        }

        //Create nodes on faces
        count = 0;
        while(topo.next_ijk_node_on_element(ijk_fine, 2, ijk_ori_fine, offset, dim, count))
        {
            std::pair<oKey*, bool> check =  insert(ijk_fine);
            oKey* key = check.first;


            if (check.second)
            {
                key->setHanging(keys_elt[topo.face_connect[count][0]],
                                keys_elt[topo.face_connect[count][1]],
                                keys_elt[topo.face_connect[count][2]],
                                keys_elt[topo.face_connect[count][3]]);
            }

            if(!check.second || (isOnBoundary(key) && !(key->isPeriodic())))
            {
                key->setRegular();
            }

//          cout<<"creating key Face"<<*key<<endl;

            count++;
        }

#else
      int hanging_on_neighbor_edge[12];
      int hanging_on_neighbor_face[6];
      // look  if hanging is mandatory by inspecting neighbor cell of curent one
      octree.lookForNodeHangingOnNeighbors3D(cell, level, ijk, hanging_on_neighbor_edge, hanging_on_neighbor_face);
      // loop on edges
      int count = 0;
      while (topo.next_ijk_node_on_element(ijk_fine, 1, ijk_ori_fine, offset, dim, count))
      {
         // find/create key
         std::pair<oKey*, bool> check = insert(ijk_fine);
         oKey* key = check.first;

         // if new key inserted
         // it is regular by default
         if (check.second)
         {
            // if it is hanging
            if (hanging_on_neighbor_edge[count])
            {
               key->setHanging(keys_elt[topo.edge_connect[count][0]], keys_elt[topo.edge_connect[count][1]]);
            }
            // it is regular
            // values have to be set as it's a new key
            else
            {
               for (int ik = 0; ik < 2; ++ik)
               {
                  oKey* key_ik = keys_elt[topo.edge_connect[count][ik]];
                  key_ik->getVals(valskey);
                  std::transform(valskey.begin(), valskey.end(), key->beginVal(), key->beginVal(), std::plus<double>());
               }
               std::transform(key->beginVal(), key->endVal(), key->beginVal(), bind2nd(std::divides<double>(), (double)2));
            }
         }
         // if already created key
         else
         {
            // if it is not hanging anymore set as regular
            if (!hanging_on_neighbor_edge[count]) key->setRegular();
         }
         if (debug) cout << "edge count " << count << key;
         count++;
      }
      // loop on face
      count = 0;
      while (topo.next_ijk_node_on_element(ijk_fine, 2, ijk_ori_fine, offset, dim, count))
      {
         // find/create key
         std::pair<oKey*, bool> check = insert(ijk_fine);
         oKey* key = check.first;

         // if new key inserted
         // it is regular by default
         if (check.second)
         {
            // if it is hanging
            if (hanging_on_neighbor_face[count])
            {
               key->setHanging(keys_elt[topo.face_connect[count][0]], keys_elt[topo.face_connect[count][1]],
                               keys_elt[topo.face_connect[count][2]], keys_elt[topo.face_connect[count][3]]);
            }
            // it is regular
            // values have to be set as it's a new key
            else
            {
               for (int ik = 0; ik < 4; ++ik)
               {
                  oKey* key_ik = keys_elt[topo.face_connect[count][ik]];
                  key_ik->getVals(valskey);
                  std::transform(valskey.begin(), valskey.end(), key->beginVal(), key->beginVal(), std::plus<double>());
               }
               std::transform(key->beginVal(), key->endVal(), key->beginVal(), bind2nd(std::divides<double>(), (double)4));
            }
         }
         // if already created key
         else
         {
            // if it is not hanging anymore set as regular
            if (!hanging_on_neighbor_face[count]) key->setRegular();
         }
         if (debug) cout << "face count " << count << key;
         count++;
      }
#endif
   }
}

/// Check if a node is on the boundary
bool oKeyManager::isOnBoundary(oKey* key) const
{
   const bool debug = false;
   const int* ijk = key->getIJK();                                                      // On recupere les coordonnees
   bool result_x = (ijk[0] == 0 || ijk[0] == topo.index_max_for_level[level_max] + 1);  // Sur le bord en x ?
   bool result_y =
       ((dim > 1 && ijk[1] == 0) || ijk[1] == topo.index_max_for_level[level_max] + 1);  // Sur les bords en y ? (en 2D ou 3D)
   bool result_z =
       ((dim > 2 && ijk[2] == 0) || ijk[2] == topo.index_max_for_level[level_max] + 1);  // sur les bords en z ? (en 3D)

   if (debug)
      cout << "Est-ce que le noeud (" << ijk[0] << "," << ijk[1] << "," << ijk[2] << ") est sur le bord ?   (X->" << result_x
           << ",y->" << result_y << ",z->" << result_z << ") TOTAL=" << (result_x || result_y || result_z)
           << "     (dim = " << dim << "  et  index_max = " << topo.index_max_for_level[level_max] + 1 << ")\n";

   return (result_x || result_y || result_z);
}

bool oKeyManager::isOnBBoxEdges(oKey* key) const
{
   const int* ijk = key->getIJK();
   bool maxMinI = ijk[0] == 0 || (ijk[0] == topo.idx_max_for_level(level_max, 0) + 1);
   bool maxMinJ = ijk[1] == 0 || (ijk[1] == topo.idx_max_for_level(level_max, 1) + 1);
   bool maxMinK = ijk[2] == 0 || (ijk[2] == topo.idx_max_for_level(level_max, 2) + 1);

   return ((maxMinI && maxMinJ) || (maxMinI && maxMinK) || (maxMinJ && maxMinK));
}

void oKeyManager::printForDebug(const std::string& filename)
{
   std::ofstream out(filename.c_str());
   printForDebug(out);
   out.close();
}

void oKeyManager::printForDebug(std::ostream& out)
{
   oKeyManager::iterator it = this->begin();
   oKeyManager::iterator ite = this->end();
   for (; it != ite; ++it)
   {
      out << *(*it);
   }
}

double distance(const oKey* k1, const oKey* k2, int i) { return fabs(((double)(k1->ijk[i] - k2->ijk[i])) * oKey::step(i)); }

void distance(const oKey* k, const std::vector<oKey*>& keys, std::vector<double>& hs)
{
   int i = 0;
   std::vector<oKey*>::const_iterator it = keys.begin();
   std::vector<oKey*>::const_iterator ite = keys.end();
   for (; it != ite; ++it, ++i)
      if (*it) hs[i] = distance(k, *it, i / 2);
}

}  // namespace xoctree
