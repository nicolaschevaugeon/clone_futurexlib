/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/

#include <iostream>
#include <sstream>

#include "xPhysSurfParameter.h"

namespace xcut
{
  using AOMD::mEntity;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// xPhysSurfByTagging class implementation ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Constructor ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  xPhysSurfParameter::xPhysSurfParameter(const xfem::xEntityToEntity  in, const xfem::xEntityToEntity  out)
    : classify_in(in), classify_out(out),fit(true),fittol(1.e-2),keep_old_partition(false),recursive(false), only_higher_dim_partition(false)
  {}
  xPhysSurfParameter::xPhysSurfParameter(const xfem::xEntityToEntity  in, const xfem::xEntityToEntity  out,std::vector<xPhysSurfByTagging *> promotors_)
    : classify_in(in), classify_out(out),fit(true),fittol(1.e-2),keep_old_partition(false),recursive(false), only_higher_dim_partition(false),promotors(promotors_)
  {}
  xPhysSurfParameter::xPhysSurfParameter(std::vector<xPhysSurfByTagging *> promotors_)
    : classify_in(xtool::xIdentity<mEntity*>()), classify_out(xtool::xIdentity<mEntity*>()),fit(true),fittol(1.e-2),keep_old_partition(false),recursive(false), only_higher_dim_partition(false),promotors(promotors_)
  {}
  xPhysSurfParameter::xPhysSurfParameter()
    : classify_in(xtool::xIdentity<mEntity*>()), classify_out(xtool::xIdentity<mEntity*>()),fit(true),fittol(1.e-2),keep_old_partition(false),recursive(false), only_higher_dim_partition(false)
  {}
  /////////////////////////////////////// End Constructor ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Destructor /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End Destructor /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Private methode ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End Private methode ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Public methode /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void xPhysSurfParameter::setFitting(bool fit_)
  {
    fit=fit_;
  }
  void xPhysSurfParameter::setFittol(double fittol_)
  {
    fittol=fittol_;
  }
  void xPhysSurfParameter::setRecursive(bool recursive_)
  {
    recursive=recursive_;
  }
  void xPhysSurfParameter::setKeepOldPartition(bool keep_old_partition_)
  {
    keep_old_partition=keep_old_partition_;
  }
  void xPhysSurfParameter::setOnlyHigherDimensionPartition(bool only_higher_dim_partition_)
  {
    only_higher_dim_partition=only_higher_dim_partition_;
  }

  /////////////////////////////////////// End Public methode /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End xPhysSurfByTagging class implementation ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


} // end of namespace
