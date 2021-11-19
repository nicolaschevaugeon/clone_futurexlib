/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XPOSTPRO_H_
#define _XPOSTPRO_H_

#include "xEntityFilter.h"
#include "xField.h"
#include "xLevelSet.h"
#include "xMesh.h"
//#include "xPhysSurf.h"

#include <algorithm>
#include <fstream>
#include <vector>

namespace xexport
{
class xPostPro
{
  protected:
   std::fstream file;

   const xfem::xIter begin;
   const xfem::xIter end;
   const xfem::xEntityFilter filter;

   const xfem::xLevelSet *lsn;
   const xfem::xLevelSet *lst;

   std::vector<AOMD::mVertex *> nodesList;
   std::vector<std::vector<int>> connectivity;
   std::vector<AOMD::mEntity::mType> element_type;
   std::vector<AOMD::mEntity *> approNodeList;
   std::vector<AOMD::mEntity *> integNodeList;
   std::vector<AOMD::mEntity *> approElemList;
   std::vector<AOMD::mEntity *> integElemList;

  public:
   /// Export dedicated constructor
   xPostPro(const xfem::xIter &, const xfem::xIter &, const xfem::xEntityFilter, const xfem::xLevelSet *_lsn = nullptr,
            const xfem::xLevelSet *_lst = nullptr);
   /// Dump dedicated constructor
   xPostPro(xfem::xMesh *m);
   virtual ~xPostPro() = default;
   ;

   virtual void openFile(const char *filename) = 0;
   virtual void closeFile() = 0;

   virtual void exportLevelSet(const xfem::xLevelSet *ls, const char *lsname) = 0;

   virtual void exportScalarAtNode(const xfem::xEval<double> &eval, const char *fieldname) = 0;
   virtual void exportVectorAtNode(const xfem::xEval<xtensor::xVector<>> &eval, const char *fieldname) = 0;
   virtual void exportTensor2AtNode(const xfem::xEval<xtensor::xTensor2<>> &eval, const char *fieldname) = 0;

   virtual void exportScalarAtElem(const xfem::xEval<double> &eval, const char *fieldname) = 0;
   virtual void exportVectorAtElem(const xfem::xEval<xtensor::xVector<>> &eval, const char *fieldname) = 0;
   virtual void exportTensor2AtElem(const xfem::xEval<xtensor::xTensor2<>> &eval, const char *fieldname) = 0;

   const xfem::xIter &getBegin() const { return begin; }
   const xfem::xIter &getEnd() const { return end; }
   const xfem::xEntityFilter &getFilter() const { return filter; }
   const std::vector<AOMD::mVertex *> &getNodesList() const { return nodesList; }
   const std::vector<AOMD::mEntity *> &getApproNodeList() const { return approNodeList; }
   const std::vector<AOMD::mEntity *> &getIntegNodeList() const { return integNodeList; }
   const std::vector<AOMD::mEntity *> &getApproElemList() const { return approElemList; }
   const std::vector<AOMD::mEntity *> &getIntegElemList() const { return integElemList; }
   const std::vector<std::vector<int>> &getConnectivity() const { return connectivity; }

  private:
   void mergeMesh();
   void mergeMesh(const xfem::xLevelSet *);
   void mergeMesh(const xfem::xLevelSet *, const xfem::xLevelSet *);
};
}  // namespace xexport

#endif /* _XPOSTPRO_H_ */
