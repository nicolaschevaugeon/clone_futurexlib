/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xPostProMSH.h"

using namespace std;

namespace xexport
{
xPostProMSH::xPostProMSH(const xfem::xIter &_begin, const xfem::xIter &_end, const xfem::xEntityFilter _filter,
                         const xfem::xLevelSet *_lsn, const xfem::xLevelSet *_lst)
    : xPostPro(_begin, _end, _filter, _lsn, _lst)
{
}

void xPostProMSH::openFile(const char *filename)
{
   file.open(filename, fstream::out);
   file << "$MeshFormat" << endl;
   file << "2.2 0 8" << endl;
   file << "$EndMeshFormat" << endl;
   file << "$Nodes" << endl;
   file << nodesList.size() << endl;
   for (unsigned int i = 0; i < nodesList.size(); ++i)
   {
      file << i + 1 << ' ' << nodesList[i]->point()(0) << ' ' << nodesList[i]->point()(1) << ' ' << nodesList[i]->point()(2)
           << endl;
   }
   file << "$EndNodes" << endl;
   file << "$Elements " << endl;
   file << connectivity.size() << endl;
   for (unsigned int i = 0; i < connectivity.size(); ++i)
   {
      file << i + 1 << ' ';
      switch (element_type[i])
      {
         case AOMD::mEntity::TRI:  // triangle
            file << 2 << ' ';
            break;
         case AOMD::mEntity::TET:  // tetra
            file << 4 << ' ';
            break;
         default:
            break;
      }
      file << 0 << ' ';
      for (unsigned int j = 0; j < connectivity[i].size(); ++j)
      {
         file << connectivity[i][j] + 1 << ' ';
      }
      file << endl;
   }
   file << "$EndElements" << endl;
}

void xPostProMSH::closeFile() { file.close(); }

void xPostProMSH::exportLevelSet(const xfem::xLevelSet *ls, const char *lsname)
{
   xfem::xEvalLevelSet<xtool::xIdentity<double>> eval_ls(*ls);
   exportScalarAtNode(eval_ls, lsname);
}

void xPostProMSH::exportScalarAtNode(const xfem::xEval<double> &eval, const char *fieldname)
{
   double d;
   file << "$NodeData" << endl;
   file << 1 << endl;
   file << '"' << fieldname << '"' << endl;
   file << 1 << endl;
   file << 0.0 << endl;
   file << 3 << endl;
   file << 0 << endl;
   file << 1 << endl;
   file << nodesList.size() << endl;
   for (unsigned int i = 0; i < nodesList.size(); ++i)
   {
      xfem::xGeomElem geo_appro(approNodeList[i]);
      xfem::xGeomElem geo_integ(integNodeList[i]);
      geo_appro.setUVWForXYZ(nodesList[i]->point());
      geo_integ.setUVWForXYZ(nodesList[i]->point());
      eval(&geo_appro, &geo_integ, d);
      file << i + 1 << ' ' << d << endl;
   }
   file << "$EndNodeData" << endl;
}

void xPostProMSH::exportVectorAtNode(const xfem::xEval<xtensor::xVector<>> &eval, const char *fieldname)
{
   xtensor::xVector<> v;
   file << "$NodeData" << endl;
   file << 1 << endl;
   file << '"' << fieldname << '"' << endl;
   file << 1 << endl;
   file << 0.0 << endl;
   file << 3 << endl;
   file << 0 << endl;
   file << 3 << endl;
   file << nodesList.size() << endl;
   for (unsigned int i = 0; i < nodesList.size(); ++i)
   {
      xfem::xGeomElem geo_appro(approNodeList[i]);
      xfem::xGeomElem geo_integ(integNodeList[i]);
      geo_appro.setUVWForXYZ(nodesList[i]->point());
      geo_integ.setUVWForXYZ(nodesList[i]->point());
      eval(&geo_appro, &geo_integ, v);
      file << i + 1 << ' ' << v(0) << ' ' << v(1) << ' ' << v(2) << endl;
   }
   file << "$EndNodeData" << endl;
}

void xPostProMSH::exportTensor2AtNode(const xfem::xEval<xtensor::xTensor2<>> &eval, const char *fieldname)
{
   xtensor::xTensor2<> t;
   file << "$NodeData" << endl;
   file << 1 << endl;
   file << '"' << fieldname << '"' << endl;
   file << 1 << endl;
   file << 0.0 << endl;
   file << 3 << endl;
   file << 0 << endl;
   file << 9 << endl;
   file << nodesList.size() << endl;
   for (unsigned int i = 0; i < nodesList.size(); ++i)
   {
      xfem::xGeomElem geo_appro(approNodeList[i]);
      xfem::xGeomElem geo_integ(integNodeList[i]);
      geo_appro.setUVWForXYZ(nodesList[i]->point());
      geo_integ.setUVWForXYZ(nodesList[i]->point());
      eval(&geo_appro, &geo_integ, t);
      file << i + 1 << ' ';
      for (int j = 0; j < 3; ++j)
      {
         for (int k = 0; k < 3; ++k)
         {
            file << t(j, k) << ' ';
         }
      }
      file << endl;
   }
   file << "$EndNodeData" << endl;
}

void xPostProMSH::exportScalarAtElem(const xfem::xEval<double> &eval, const char *fieldname)
{
   double d;
   file << "$ElementData" << endl;
   file << 1 << endl;
   file << '"' << fieldname << '"' << endl;
   file << 1 << endl;
   file << 0.0 << endl;
   file << 3 << endl;
   file << 0 << endl;
   file << 1 << endl;
   file << connectivity.size() << endl;
   for (unsigned int i = 0; i < connectivity.size(); ++i)
   {
      xtensor::xPoint p = integElemList[i]->getCentroid();
      xfem::xGeomElem geo_appro(approElemList[i]);
      xfem::xGeomElem geo_integ(integElemList[i]);
      geo_appro.setUVWForXYZ(p);
      geo_integ.setUVWForXYZ(p);
      eval(&geo_appro, &geo_integ, d);
      file << i + 1 << ' ' << d << endl;
   }
   file << "$EndElementData" << endl;
}

void xPostProMSH::exportVectorAtElem(const xfem::xEval<xtensor::xVector<>> &eval, const char *fieldname)
{
   xtensor::xVector<> v;
   file << "$ElementData" << endl;
   file << 1 << endl;
   file << '"' << fieldname << '"' << endl;
   file << 1 << endl;
   file << 0.0 << endl;
   file << 3 << endl;
   file << 0 << endl;
   file << 3 << endl;
   file << connectivity.size() << endl;
   for (unsigned int i = 0; i < connectivity.size(); ++i)
   {
      xtensor::xPoint p = integElemList[i]->getCentroid();
      xfem::xGeomElem geo_appro(approElemList[i]);
      xfem::xGeomElem geo_integ(integElemList[i]);
      geo_appro.setUVWForXYZ(p);
      geo_integ.setUVWForXYZ(p);
      eval(&geo_appro, &geo_integ, v);
      file << i + 1 << ' ' << v(0) << ' ' << v(1) << ' ' << v(2) << endl;
   }
   file << "$EndElementData" << endl;
}

void xPostProMSH::exportTensor2AtElem(const xfem::xEval<xtensor::xTensor2<>> &eval, const char *fieldname)
{
   xtensor::xTensor2<> t;
   file << "$ElementData" << endl;
   file << 1 << endl;
   file << '"' << fieldname << '"' << endl;
   file << 1 << endl;
   file << 0.0 << endl;
   file << 3 << endl;
   file << 0 << endl;
   file << 9 << endl;
   file << connectivity.size() << endl;
   for (unsigned int i = 0; i < connectivity.size(); ++i)
   {
      xtensor::xPoint p = integElemList[i]->getCentroid();
      xfem::xGeomElem geo_appro(approElemList[i]);
      xfem::xGeomElem geo_integ(integElemList[i]);
      geo_appro.setUVWForXYZ(p);
      geo_integ.setUVWForXYZ(p);
      eval(&geo_appro, &geo_integ, t);
      file << i + 1 << ' ';
      for (int j = 0; j < 3; ++j)
      {
         for (int k = 0; k < 3; ++k)
         {
            file << t(j, k) << ' ';
         }
      }
      file << endl;
   }
   file << "$EndElementData" << endl;
}

}  // namespace xexport
