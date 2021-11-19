/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _XPOSTPROMSH_H_
#define _XPOSTPROMSH_H_

#include "xPostPro.h"
#include "xVector.h"
#include "xTensor2.h"


namespace xexport{

 
class xPostProMSH : virtual public xPostPro{
  
public:
  xPostProMSH(const xfem::xIter &, const xfem::xIter &, const xfem::xEntityFilter, const xfem::xLevelSet *lsn=nullptr, const xfem::xLevelSet *lst=nullptr);

  void openFile(const char* filename) override;
  void closeFile() override;
  
  void exportLevelSet(const xfem::xLevelSet *ls, const char* lsname) override;
  
  void exportScalarAtNode(const xfem::xEval<double> &eval, const char* fieldname) override;
  void exportVectorAtNode(const xfem::xEval<xtensor::xVector<> > &eval, const char* fieldname) override;
  void exportTensor2AtNode(const xfem::xEval<xtensor::xTensor2<> > &eval, const char* fieldname) override;

  void exportScalarAtElem(const xfem::xEval<double> &eval, const char* fieldname) override;
  void exportVectorAtElem(const xfem::xEval<xtensor::xVector<> > &eval, const char* fieldname) override;
  void exportTensor2AtElem(const xfem::xEval<xtensor::xTensor2<> > &eval, const char* fieldname) override;
  
};
}


#endif /* _XPOSTPROMSH_H_ */
