/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/

#ifndef ___XPHYSSURFPARAMETER_H
#define ___XPHYSSURFPARAMETER_H
#include <vector>

#include "xEntityToEntity.h"

namespace xfem
{
  class xLevelSet;

}

namespace xcut
{
  class xPhysSurfByTagging;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// xPhysSurfParameter class ///////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // html for doxygen : please respect html tag will changing  comments
  /// This class is a utility class for  xPhysSurfByTagging.  It define a set of parameter for class xPhysSurfByTagging
  //! <br/> 
  //
  class xPhysSurfParameter
  {
  public:
    // constructor
    xPhysSurfParameter(const xfem::xEntityToEntity  in, const xfem::xEntityToEntity  out,std::vector<xPhysSurfByTagging *> promotors);
    xPhysSurfParameter(const xfem::xEntityToEntity  in, const xfem::xEntityToEntity  out);
    xPhysSurfParameter(std::vector<xPhysSurfByTagging *> promotors);
    xPhysSurfParameter();

    // public methodes ///////////////
    void setFitting(bool fit_);
    void setFittol(double fittol_);
    void setRecursive(bool recurcive_);
    void setKeepOldPartition(bool keep_old_partition_);
    void setOnlyHigherDimensionPartition(bool only_higher_dim_partition_);
    //
    /// Give classifyer ("in")
    xfem::xEntityToEntity &getClassifyerIn(){return classify_in; }
    /// Give classifyer ("out")
    xfem::xEntityToEntity &getClassifyerOut(){return classify_out; }
    /// Give Fitting parameter
    bool getFitting(){return fit; }
    /// Give Filter tolerance
    double getFittol(){return fittol; }
    /// Give keep_old_partition parameter
    bool getKeepOldPartition(){return keep_old_partition; }
    /// Give recursive parameter
    bool getRecursive(){return recursive; }
    /// Give promotors
    std::vector<xPhysSurfByTagging *>  getPromotors(){return promotors; }
    /// Give whether only high order partition has to be done
    bool getPartitionParameter(){return only_higher_dim_partition;}

  private:
    // private members ////////////////
        
    // classifyer
    xfem::xEntityToEntity classify_in;
    xfem::xEntityToEntity classify_out;

    // level set parameter
    bool fit;
    double fittol;

    // recurtion and partition
    bool keep_old_partition;
    bool recursive;
    bool only_higher_dim_partition;

    // promotors container
    std::vector<xPhysSurfByTagging *> promotors;

  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End xPhysSurfParameter class ///////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // end of namespace

#endif
