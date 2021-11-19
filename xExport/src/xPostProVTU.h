/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _XPOSTPROVTU_H_
#define _XPOSTPROVTU_H_

#include "xPostPro.h"
#include "xVector.h"
#include "xTensor2.h"


namespace xexport{
  
class xPostProVTU : virtual public xPostPro{

  template <class T>
  struct PostProField{
  PostProField(const std::string _name):name(_name){}
    std::vector<T> field;
    const std::string name;
  };

  typedef PostProField<double> ScalarField;
  typedef PostProField<xtensor::xVector<> > VectorField;
  typedef PostProField<xtensor::xTensor2<> > Tensor2Field;
  
  std::vector<ScalarField> scalarAtNode;
  std::vector<ScalarField> scalarAtElem;
  std::vector<VectorField> vectorAtNode;
  std::vector<VectorField> vectorAtElem;
  std::vector<Tensor2Field> tensor2AtNode;
  std::vector<Tensor2Field> tensor2AtElem;
  
public:
  xPostProVTU(const xfem::xIter &, const xfem::xIter &, const xfem::xEntityFilter, const xfem::xLevelSet *lsn=nullptr, const xfem::xLevelSet *lst=nullptr);
  xPostProVTU(xfem::xMesh *m);

  void openFile(const char* filename) override;
  void closeFile() override;
  
  void exportLevelSet(const xfem::xLevelSet *ls, const char * lsname) override;
  
  void exportScalarAtNode(const xfem::xEval<double> &eval, const char * fieldname) override;
  void exportVectorAtNode(const xfem::xEval<xtensor::xVector<> > &eval, const char * fieldname) override;
  void exportTensor2AtNode(const xfem::xEval<xtensor::xTensor2<> > &eval, const char * fieldname) override;

  void exportScalarAtElem(const xfem::xEval<double> &eval, const char * fieldname) override;
  void exportVectorAtElem(const xfem::xEval<xtensor::xVector<> > &eval, const char * fieldname) override;
  void exportTensor2AtElem(const xfem::xEval<xtensor::xTensor2<> > &eval, const char * fieldname) override;

private:
  void writeScalarField(const ScalarField &f);
  void writeVectorField(const VectorField &f);
  void writeTensor2Field(const Tensor2Field &f);

  const std::string Tab{"  "};
  std::string Tab2, Tab3, Tab4, Tab5, Tab6, Tab7;

};
}


#endif /* _XPOSTPROVTU_H_ */
