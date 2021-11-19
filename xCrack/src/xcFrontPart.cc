/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
// AOMD
#include "AOMD.h"
#include "mEdge.h"
// xcrack
#include "xcFrontPart.h"

using namespace AOMD;
using namespace xfem;

xcFrontPartBase::xcFrontPartBase(const xMesh &_front_mesh, const std::string &_frontpartname, const xMesh &_mesh,
                                 std::function<double(mVertex *)> _front_distance, const xParseData &_parameters)
    : front_mesh(_front_mesh),
      front_part_name(_frontpartname),
      front_region(&_front_mesh, _frontpartname),
      mesh(_mesh),
      front_distance(_front_distance),
      fronttype(None),
      parameters(_parameters){};

xcFrontPartBase::~xcFrontPartBase() = default;

// ------------------------------------------------------------------------------------------------------------------------------------

xcFrontPartPoint::xcFrontPartPoint(const xMesh &_front_mesh, const std::string &_frontname, const xMesh &_mesh,
                                   std::function<double(mVertex *)> _front_distance, const xParseData &_parameters)
    : xcFrontPartBase(_front_mesh, _frontname, _mesh, _front_distance, _parameters)
{
   fronttype = Point;
   totallenght = 0.;
}

void xcFrontPartPoint::ParametrizeFront(){};

void xcFrontPartPoint::restrict_ls1d_support(xRegion new_support, double init_val) {}

// ------------------------------------------------------------------------------------------------------------------------------------

xcFrontPartLineOpen::xcFrontPartLineOpen(const xMesh &_front_mesh, const std::string &_frontname, const xMesh &_mesh,
                                         std::function<double(mVertex *)> _front_distance, const xParseData &_parameters,
                                         mVertex *_start, mVertex *_end)
    : xcFrontPartBase(_front_mesh, _frontname, _mesh, _front_distance, _parameters), start(_start), end(_end)
{
   fronttype = LineOpen;
   lss1d.setSupport(front_region);
   ParametrizeFront();
};

void xcFrontPartLineOpen::ParametrizeFront()
{
   // Parametrize1dMeshBranch
   xinterface::aomd::xAttachedDataManagerAOMD<double> s;
   mVertex *vcurrent = start;
   s.setData(*vcurrent) = 0.;
   mEdge *ecurrent = static_cast<mEdge *>(vcurrent->get(1, 0));
   totallenght = 0.;
   local_vertices.clear();
   local_vertices.push_back(vcurrent);
   while (vcurrent != end)
   {
      mVertex *vnext = static_cast<mVertex *>(E_otherVertex(ecurrent, vcurrent));
      xtensor::xVector<> v01(vnext->point(), vcurrent->point());
      totallenght += v01.mag();
      s.setData(*vnext) = totallenght;
      mEdge *enext = static_cast<mEdge *>(vnext->get(1, 0));
      if (enext == ecurrent)
      {
         if (vnext->size(1) > 1) enext = static_cast<mEdge *>(vnext->get(1, 1));
      }
      ecurrent = enext;
      vcurrent = vnext;
      local_vertices.push_back(vcurrent);
   }
   for (mVertex *pv : local_vertices)
   {
      lss1d(pv) = s.at(*pv) * 2. / totallenght - 1.;
      s.deleteData(*pv);
   }
}

void xcFrontPartLineOpen::restrict_ls1d_support(xRegion new_support, double init_val) { lss1d.restrictTo(new_support, init_val); }

// ------------------------------------------------------------------------------------------------------------------------------------

xcFrontPartLineClose::xcFrontPartLineClose(const xMesh &_front_mesh, const std::string &_frontname, const xMesh &_mesh,
                                           std::function<double(mVertex *)> _front_distance, const xParseData &_parameters,
                                           mVertex *_start)
    : xcFrontPartBase(_front_mesh, _frontname, _mesh, _front_distance, _parameters), start(_start)
{
   fronttype = LineClose;
   lss1dCos.setSupport(front_region);
   lss1dSin.setSupport(front_region);
   ParametrizeFront();
}

void xcFrontPartLineClose::ParametrizeFront()
{
   // Parametrize1dMeshLoop
   xinterface::aomd::xAttachedDataManagerAOMD<double> s;
   mVertex *vcurrent = start;
   s.setData(*vcurrent) = 0.;
   mEdge *ecurrent = static_cast<mEdge *>(vcurrent->get(1, 0));
   totallenght = 0.;
   local_vertices.clear();
   local_vertices.push_back(vcurrent);
   while (1)
   {
      mVertex *vnext = static_cast<mVertex *>(E_otherVertex(ecurrent, vcurrent));
      xtensor::xVector<> v01(vnext->point(), vcurrent->point());
      totallenght += v01.mag();
      if (vnext == start) break;
      mEdge *enext = static_cast<mEdge *>(vnext->get(1, 0));
      if (enext == ecurrent)
      {
         enext = static_cast<mEdge *>(vnext->get(1, 1));
      }
      ecurrent = enext;
      vcurrent = vnext;
      s.setData(*vcurrent) = totallenght;
      local_vertices.push_back(vcurrent);
   }
   for (mVertex *v : local_vertices)
   {
      double theta = s.at(*v) * 2 * M_PI / totallenght;
      // cout << theta << endl;
      lss1dCos(v) = cos(theta);
      lss1dSin(v) = sin(theta);
      s.deleteData(*v);
   }
}

void xcFrontPartLineClose::restrict_ls1d_support(xRegion new_support, double init_val)
{
   lss1dCos.restrictTo(new_support, init_val);
   lss1dSin.restrictTo(new_support, init_val);
}
