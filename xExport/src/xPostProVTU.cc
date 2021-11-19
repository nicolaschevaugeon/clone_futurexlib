/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xPostProVTU.h"
using namespace std;
using AOMD::mVertex;

namespace xexport
{
xPostProVTU::xPostProVTU(xfem::xMesh *m) : xPostPro(m)
{
   Tab2 = Tab + Tab;
   Tab3 = Tab + Tab + Tab;
   Tab4 = Tab + Tab + Tab + Tab;
   Tab5 = Tab + Tab + Tab + Tab + Tab;
   Tab6 = Tab + Tab + Tab + Tab + Tab + Tab;
   Tab7 = Tab + Tab + Tab + Tab + Tab + Tab + Tab;
}

xPostProVTU::xPostProVTU(const xfem::xIter &_begin, const xfem::xIter &_end, const xfem::xEntityFilter _filter,
                         const xfem::xLevelSet *_lsn, const xfem::xLevelSet *_lst)
    : xPostPro(_begin, _end, _filter, _lsn, _lst)
{
   Tab2 = Tab + Tab;
   Tab3 = Tab + Tab + Tab;
   Tab4 = Tab + Tab + Tab + Tab;
   Tab5 = Tab + Tab + Tab + Tab + Tab;
   Tab6 = Tab + Tab + Tab + Tab + Tab + Tab;
   Tab7 = Tab + Tab + Tab + Tab + Tab + Tab + Tab;
}

void xPostProVTU::openFile(const char *filename)
{
   file.open(filename, fstream::out);
   file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
   file << Tab + "<UnstructuredGrid>" << endl;
   file << Tab2 + "<Piece NumberOfPoints=\"" << nodesList.size() << "\" NumberOfCells=\"" << connectivity.size() << "\">" << endl;
   file << Tab3 + "<Points>" << endl;
   file << Tab4 + "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
   file << Tab5;
   for (std::vector<mVertex *>::iterator it = nodesList.begin(); it != nodesList.end(); ++it)
   {
      file << (*it)->point()(0) << ' ' << (*it)->point()(1) << ' ' << (*it)->point()(2) << ' ';
   }
   file << endl;
   file << Tab4 + "</DataArray>" << endl;
   file << Tab3 + "</Points>" << endl;
   file << Tab3 + "<Cells>" << endl;
   file << Tab4 + "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
   file << Tab5;
   for (std::vector<std::vector<int>>::iterator it = connectivity.begin(); it != connectivity.end(); ++it)
   {
      for (vector<int>::iterator itt = it->begin(); itt != it->end(); ++itt)
      {
         file << *itt << ' ';
      }
   }
   file << endl;
   file << Tab4 + "</DataArray>" << endl;
   file << Tab4 + "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
   file << Tab5;
   int offset = 0;
   for (unsigned int i = 0; i < connectivity.size(); ++i)
   {
      offset += connectivity[i].size();
      file << offset << ' ';
   }
   file << endl;
   file << Tab4 + "</DataArray>" << endl;
   file << Tab4 + "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
   for (unsigned int i = 0; i < connectivity.size(); ++i)
   {
      switch (element_type[i])
      {
         case AOMD::mEntity::TRI:
            file << 5 << ' ';
            break;
         case AOMD::mEntity::TET:
            file << 10 << ' ';
            break;
         case AOMD::mEntity::EDGE:
            file << 3 << ' ';
            break;
         default:
            std::cout << "ERROR type number " << element_type[i] << std::endl;
            exit(0);
      }
   }
   file << endl;
   file << Tab4 + "</DataArray>" << endl;
   file << Tab3 + "</Cells>" << endl;
}

void xPostProVTU::closeFile()
{
   file << Tab3 + "<PointData>" << endl;
   for (unsigned int i = 0; i < scalarAtNode.size(); ++i)
   {
      writeScalarField(scalarAtNode[i]);
   }
   for (unsigned int i = 0; i < vectorAtNode.size(); ++i)
   {
      writeVectorField(vectorAtNode[i]);
   }
   for (unsigned int i = 0; i < tensor2AtNode.size(); ++i)
   {
      writeTensor2Field(tensor2AtNode[i]);
   }
   file << Tab3 + "</PointData>" << endl;
   file << Tab3 + "<CellData>" << endl;
   for (unsigned int i = 0; i < scalarAtElem.size(); ++i)
   {
      writeScalarField(scalarAtElem[i]);
   }
   for (unsigned int i = 0; i < vectorAtElem.size(); ++i)
   {
      writeVectorField(vectorAtElem[i]);
   }
   for (unsigned int i = 0; i < tensor2AtElem.size(); ++i)
   {
      writeTensor2Field(tensor2AtElem[i]);
   }
   file << Tab3 + "</CellData>" << endl;
   file << Tab2 + "</Piece>\n";
   file << Tab + "</UnstructuredGrid>\n</VTKFile>" << endl;
   file.close();
}

void xPostProVTU::exportLevelSet(const xfem::xLevelSet *ls, const char *fieldname)
{
   xfem::xEvalLevelSet<xtool::xIdentity<double>> eval_ls(*ls);
   exportScalarAtNode(eval_ls, fieldname);
}

void xPostProVTU::exportScalarAtNode(const xfem::xEval<double> &eval, const char *fieldname)
{
   ScalarField f(fieldname);
   double v;
   for (unsigned int i = 0; i < nodesList.size(); ++i)
   {
      xfem::xGeomElem geo_appro(approNodeList[i]);
      xfem::xGeomElem geo_integ(integNodeList[i]);
      geo_appro.setUVWForXYZ(nodesList[i]->point());
      geo_integ.setUVWForXYZ(nodesList[i]->point());
      eval(&geo_appro, &geo_integ, v);
      f.field.push_back(v);
   }
   scalarAtNode.push_back(f);
}

void xPostProVTU::exportVectorAtNode(const xfem::xEval<xtensor::xVector<>> &eval, const char *fieldname)
{
   VectorField f(fieldname);
   xtensor::xVector<> v;
   for (unsigned int i = 0; i < nodesList.size(); ++i)
   {
      xfem::xGeomElem geo_appro(approNodeList[i]);
      xfem::xGeomElem geo_integ(integNodeList[i]);
      geo_appro.setUVWForXYZ(nodesList[i]->point());
      geo_integ.setUVWForXYZ(nodesList[i]->point());
      eval(&geo_appro, &geo_integ, v);
      f.field.push_back(v);
   }
   vectorAtNode.push_back(f);
}

void xPostProVTU::exportTensor2AtNode(const xfem::xEval<xtensor::xTensor2<>> &eval, const char *fieldname)
{
   Tensor2Field f(fieldname);
   xtensor::xTensor2<> v;
   for (unsigned int i = 0; i < nodesList.size(); ++i)
   {
      xfem::xGeomElem geo_appro(approNodeList[i]);
      xfem::xGeomElem geo_integ(integNodeList[i]);
      geo_appro.setUVWForXYZ(nodesList[i]->point());
      geo_integ.setUVWForXYZ(nodesList[i]->point());
      eval(&geo_appro, &geo_integ, v);
      f.field.push_back(v);
   }
   tensor2AtNode.push_back(f);
}

void xPostProVTU::exportScalarAtElem(const xfem::xEval<double> &eval, const char *fieldname)
{
   ScalarField f(fieldname);
   double v;
   for (unsigned int i = 0; i < connectivity.size(); ++i)
   {
      xtensor::xPoint p = integElemList[i]->getCentroid();
      xfem::xGeomElem geo_appro(approElemList[i]);
      xfem::xGeomElem geo_integ(integElemList[i]);
      geo_appro.setUVWForXYZ(p);
      geo_integ.setUVWForXYZ(p);
      eval(&geo_appro, &geo_integ, v);
      f.field.push_back(v);
   }
   scalarAtElem.push_back(f);
}

void xPostProVTU::exportVectorAtElem(const xfem::xEval<xtensor::xVector<>> &eval, const char *fieldname)
{
   VectorField f(fieldname);
   xtensor::xVector<> v;
   for (unsigned int i = 0; i < connectivity.size(); ++i)
   {
      xtensor::xPoint p = integElemList[i]->getCentroid();
      xfem::xGeomElem geo_appro(approElemList[i]);
      xfem::xGeomElem geo_integ(integElemList[i]);
      geo_appro.setUVWForXYZ(p);
      geo_integ.setUVWForXYZ(p);
      eval(&geo_appro, &geo_integ, v);
      f.field.push_back(v);
   }
   vectorAtElem.push_back(f);
}

void xPostProVTU::exportTensor2AtElem(const xfem::xEval<xtensor::xTensor2<>> &eval, const char *fieldname)
{
   Tensor2Field f(fieldname);
   xtensor::xTensor2<> v;
   for (unsigned int i = 0; i < connectivity.size(); ++i)
   {
      xtensor::xPoint p = integElemList[i]->getCentroid();
      xfem::xGeomElem geo_appro(approElemList[i]);
      xfem::xGeomElem geo_integ(integElemList[i]);
      geo_appro.setUVWForXYZ(p);
      geo_integ.setUVWForXYZ(p);
      eval(&geo_appro, &geo_integ, v);
      f.field.push_back(v);
   }
   tensor2AtElem.push_back(f);
}

void xPostProVTU::writeScalarField(const ScalarField &f)
{
   file << Tab4 + "<DataArray type=\"Float32\" Name=\"" << f.name << "\" format=\"ascii\">" << endl;
   file << Tab5;
   for (unsigned int i = 0; i < f.field.size(); ++i)
   {
      file << f.field[i] << ' ';
   }
   file << endl;
   file << Tab4 + "</DataArray>" << endl;
}

void xPostProVTU::writeVectorField(const VectorField &f)
{
   file << Tab4 + "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"" << f.name << "\" format=\"ascii\">" << endl;
   file << Tab5;
   for (unsigned int i = 0; i < f.field.size(); ++i)
   {
      file << f.field[i](0) << ' ' << f.field[i](1) << ' ' << f.field[i](2) << ' ';
   }
   file << endl;
   file << Tab4 + "</DataArray>" << endl;
}

void xPostProVTU::writeTensor2Field(const Tensor2Field &f)
{
   file << Tab4 + "<DataArray type=\"Float32\" NumberOfComponents=\"9\" Name=\"" << f.name << "\" format=\"ascii\">" << endl;
   file << Tab5;
   for (unsigned int i = 0; i < f.field.size(); ++i)
   {
      for (int j = 0; j < 3; ++j)
      {
         for (int k = 0; k < 3; ++k)
         {
            file << f.field[i](j, k) << ' ';
         }
      }
   }
   file << endl;
   file << Tab4 + "</DataArray>" << endl;
}
}  // namespace xexport
