/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

// -*- C++ -*-

#ifndef _OMAPPING_H__
#define _OMAPPING_H__

#include "oTopo.h"

namespace xoctree
{
class oMapping
{
  public:
   virtual ~oMapping() = default;
   ;
   oMapping(const oTopo& topo_, int d_) : Dim(d_), topo(topo_) {}
   virtual void getBox(const int* ijk, int level, double* box_inf, double* box_sup) const = 0;
   virtual void ijk2xyz(const int* ijk, int level, double* xyz) const = 0;
   // virtual double getStep(int level, int c) const = 0;
   virtual const double* getStep(int level) const = 0;
   int getDim() const { return Dim; }
   const oTopo& getTopo() const { return topo; }

  protected:
   const int Dim;
   const oTopo& topo;
};

// metric
class oMappingCartesian : public oMapping
{
  public:
   oMappingCartesian(const oTopo& topo_, int dim_, double* box_inf_, double* box_sup_) : oMapping(topo_, dim_)
   {
      const bool debug = false;
      if (box_inf_ == nullptr)
         std::fill(box_inf_level_0, box_inf_level_0 + 3, 0.);
      else
         std::copy(box_inf_, box_inf_ + 3, box_inf_level_0);
      if (box_sup_ == nullptr)
         std::fill(box_sup_level_0, box_sup_level_0 + 3, 1.);
      else
         std::copy(box_sup_, box_sup_ + 3, box_sup_level_0);
      if (debug)
      {
         std::cout << " box inf " << std::endl;
         std::copy(box_inf_level_0, box_inf_level_0 + 3, std::ostream_iterator<double>(std::cout, " "));
         std::cout << " box sup " << std::endl;
         std::copy(box_sup_level_0, box_sup_level_0 + 3, std::ostream_iterator<double>(std::cout, " "));
      }
      step = new double*[topo.MAX_DEPTH + 1];
      for (int l = 0; l <= topo.MAX_DEPTH; ++l)
      {
         step[l] = new double[3];
         //       double p =  (double) topo.pow_base2[l];
         for (int i = 0; i < 3; ++i)
         {
            double p = (double)topo.pow_base2[l + topo.motif[i]];
            step[l][i] = (box_sup_level_0[i] - box_inf_level_0[i]) / p;
            if (debug) std::cout << " step at level " << l << " in direction " << i << " is " << step[l][i] << std::endl;
         }
      }
   }

   void getBox(const int* ijk, int level, double* box_inf, double* box_sup) const override
   {
      const bool debug = false;
      ijk2xyz(ijk, level, box_inf);
      transform(box_inf, box_inf + Dim, step[level], box_sup, std::plus<double>());
      if (Dim == 2) box_sup[2] = box_inf[2] = 0;
      if (debug)
      {
         std::cout << " for ijk " << std::endl;
         std::copy(ijk, ijk + 3, std::ostream_iterator<int>(std::cout, " "));
         std::cout << std::endl << " at level " << level << std::endl;
         std::cout << std::endl << " box inf " << std::endl;
         std::copy(box_inf, box_inf + 3, std::ostream_iterator<double>(std::cout, " "));
         std::cout << std::endl << " box sup " << std::endl;
         std::copy(box_sup, box_sup + 3, std::ostream_iterator<double>(std::cout, " "));
         std::cout << std::endl;
      }
      return;
   }
   void ijk2xyz(const int* ijk, int level, double* xyz) const override
   {
      // printf("Coucouuuuuuuuuuuuuuuuuuu %d  %d  %d\n", ijk[0], ijk[1], ijk[2]);
      transform(ijk, ijk + 3, step[level], xyz, std::multiplies<double>());
      transform(box_inf_level_0, box_inf_level_0 + 3, xyz, xyz, std::plus<double>());
      // printf("XXXXXXXXXYYYYYYYYYYYYYZZZ %e  %e  %e \n", xyz[0], xyz[1], xyz[2]);
   }

   // double getStep(int level, int c) const {return step[level][c]; }
   const double* getStep(int level) const override { return step[level]; }

  private:
   double box_inf_level_0[3];
   double box_sup_level_0[3];
   //   double step[topa::MAX_DEPTH][3];
   double** step;
};

}  // namespace xoctree

#endif
