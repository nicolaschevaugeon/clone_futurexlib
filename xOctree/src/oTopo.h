/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#ifndef _OTOPO_H__
#define _OTOPO_H__

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <iterator>
#include <vector>

namespace xoctree
{
class oTopo
{
  public:
   oTopo(int* motif_ = nullptr)
   {
      if (!motif_)
         std::fill(motif, motif + 3, 0);
      else
         std::copy(motif_, motif_ + 3, motif);

      assert(MAX_DEPTH + (*std::max_element(motif, motif + 3)) < MAX_ABSOLUTE_DEPTH);

      for (int l = 0; l < MAX_ABSOLUTE_DEPTH; ++l)
      {
         index_max_for_levels[0][l] = pow_base2[motif[0] + l] - 1;
         index_max_for_levels[1][l] = pow_base2[motif[1] + l] - 1;
         index_max_for_levels[2][l] = pow_base2[motif[2] + l] - 1;
      }

      for (int iC = 0; iC < 8; ++iC)
      {
         ijk_corner_node_level0[iC][0] = ijk_node_shift[0][iC][0] * pow_base2[motif[0]];
         ijk_corner_node_level0[iC][1] = ijk_node_shift[0][iC][1] * pow_base2[motif[1]];
         ijk_corner_node_level0[iC][2] = ijk_node_shift[0][iC][2] * pow_base2[motif[2]];
      }
   };

   static const int MAX_ABSOLUTE_DEPTH = 27;  //(29 en long int) max depth taking into account the initial motif
   static const int MAX_DEPTH = 15;           //(17 en long int) depth relative to initial motif

   static const int base2_ijk[8][3];
   static const int nb_entities[4][3];
   static const int ijk_node_shift[3][12][3];
   static const int ijk_center_shift[4][3];
   static const int edge_connect[12][2];
   static const int face_connect[6][4];
   static const int face_edge_connect[6][4];
   static const int ijk_edge_neighbor_shift_2D[4][3];
   static const int ijk_edge_neighbor_shift_3D[12][3][3];
   static const int ijk_face_neighbor_shift[6][3];
   static const int nb_nodes_connected_to_node[4];
   static const int connected_nodes[4][8][3];
   static const int stencil_location[4][8][3];

   //   static const int motif[3] ;
   int motif[3];

   // for (int l=0; l <= topo::MAX_DEPTH ; ++l) {  index_max_for_level[l] = powint(2,l)-1;}
   // modified as following :
   //     static const int index_max_for_level[topo::MAX_DEPTH+1] ;
   static const int generation_size[3][MAX_DEPTH + 1];
   static const int pow_base2[MAX_ABSOLUTE_DEPTH + 1];
   static const int index_max_for_level[MAX_ABSOLUTE_DEPTH + 1];
   int index_max_for_levels[3][MAX_ABSOLUTE_DEPTH + 1];
   int ijk_corner_node_level0[8][3];

   int idx_max_for_level(int l, int idir) const { return index_max_for_levels[idir][l]; }

   const int* getMotif() const { return motif; }

   inline static unsigned int powint(unsigned int nValue, int nPower)
   {
      unsigned int nReturn = 1;
      while (nPower > 0)
      {
         nReturn *= nValue;
         nPower -= 1;
      }
      return nReturn;
   }

   inline void cartesian2finest_cartesian(const int* ijk, int level, int* ijk_fine, int level_max) const
   {
      assert(level <= level_max);
      std::transform(ijk, ijk + 3, ijk_fine, bind2nd(std::multiplies<int>(), pow_base2[level_max - level]));
   }

   inline static void ijk_center_on_element(int* ijk, const int* ijk_ori, int offset, int dim)
   {
      offset /= 2;
      transform(ijk_center_shift[dim], ijk_center_shift[dim] + 3, ijk, bind2nd(std::multiplies<int>(), offset));
      transform(ijk_ori, ijk_ori + 3, ijk, ijk, std::plus<int>());
   }

   inline static bool next_ijk_node_on_element(int* ijk, int on_what, const int* ijk_ori, int offset, int dim, int count)
   {
      const bool debug = false;
      if (on_what > 0) offset /= 2;
      if (count >= nb_entities[dim][on_what]) return false;
      transform(ijk_node_shift[on_what][count], ijk_node_shift[on_what][count] + 3, ijk, bind2nd(std::multiplies<int>(), offset));

      transform(ijk_ori, ijk_ori + 3, ijk, ijk, std::plus<int>());
      if (debug)
      {
         std::cout << "offset" << offset << " count " << count << " ijk_ori ";
         std::copy(ijk_ori, ijk_ori + 3, std::ostream_iterator<int>(std::cout, " "));
         std::cout << " ijk found ";
         std::copy(ijk, ijk + 3, std::ostream_iterator<int>(std::cout, " "));
         std::cout << std::endl;
      }
      return true;
   }

   inline static void decompose(int val, int base, int* beg, int& size)
   {
      const bool debug = false;
      if (debug)
      {
         std::cout << " decomposition of " << val << " in base " << base << " is " << std::endl;
      }
      size = 0;
      do
      {
         beg[size++] = val % base;
         val /= base;
      } while (val);
      std::reverse(beg, beg + size);

      if (debug)
      {
         std::copy(beg, beg + size, std::ostream_iterator<int>(std::cout, " "));
         std::cout << std::endl << " and its size is " << size << std::endl;
      }
   }

   inline static void compose(int& val, int base, int* beg, int size)
   {
      val = *beg++;
      for (int i = 0; i < size - 1; ++i) val = base * val + *beg++;
   }

   struct periodic : public std::binary_function<int, int, int>
   {
      periodic(int imax_) : imax(imax_) {}
      int operator()(int i, int p)
      {
         if (i < 0) return (p) ? imax : -1;
         if (i > imax) return (p) ? 0 : -1;
         return i;
      }
      int imax;
   };

   inline bool nextNeighbor(const int* ijk, int l, int* ijkn, int on_what, bool& flag, int dim, const int* per, int count,
                            int id_neighbor) const
   {
      if (count < nb_entities[dim][on_what])
      {
         assert(dim > 1 && dim < 4);
         assert(id_neighbor > -1 && id_neighbor < 3);
         assert(on_what > 0 && on_what < 3);

         const int* shift =
             (dim > 2) ? ((on_what < 2) ? ijk_edge_neighbor_shift_3D[count][id_neighbor] : ijk_face_neighbor_shift[count])
                       : ((on_what < 2) ? ijk_edge_neighbor_shift_2D[count] : (int*)nullptr);

         if (!shift)
         {
            std::cout << "unknown case \nfile " << __FILE__ << " line " << __LINE__ << std::endl;
            throw;
         }

         transform(ijk, ijk + 3, shift, ijkn, std::plus<int>());

         periodic perio0(index_max_for_levels[0][l]);
         periodic perio1(index_max_for_levels[1][l]);
         periodic perio2(index_max_for_levels[2][l]);

         ijkn[0] = perio0(ijkn[0], per[0]);
         ijkn[1] = perio1(ijkn[1], per[1]);
         ijkn[2] = perio2(ijkn[2], per[2]);

         flag = (ijkn[0] < 0) ? false : ((ijkn[1] < 0) ? false : ((ijkn[2] < 0) ? false : true));
      }
      else
         return false;

      return true;
   }
};

}  // namespace xoctree

#endif
