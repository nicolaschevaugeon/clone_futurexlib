/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include <fstream>
#include <sstream>

#ifdef WITH_GZIPOUTPUT
#include <boost/iostreams/filter/gzip.hpp>
#endif


#include "workInProgress.h"
#include "xExportGmsh.h"

namespace xexport
{
using xtensor::xPoint;
using namespace std;

xExportGmsh::xExportGmsh(MPI_Comm world_) : xExport(world_) { filename_extension = ".pos"; }

void xExportGmsh::exportPoint(const xPoint &P1, const double &val1)
{
   double x[1] = {P1(0)};
   double y[1] = {P1(1)};
   double z[1] = {P1(2)};
   double v[1] = {val1};
   addDrawable(xGmshDrawable ::SP, 1, x, y, z, v);
}

void xExportGmsh::exportPoint(const xPoint &P1, const xtensor::xVector<> &val1)
{
   double x[1] = {P1(0)};
   double y[1] = {P1(1)};
   double z[1] = {P1(2)};
   double v[3] = {val1(0), val1(1), val1(2)};
   addDrawable(xGmshDrawable ::VP, 1, x, y, z, v);
}

void xExportGmsh::exportPoint(const xPoint &P1, const xtensor::xTensor2<> &val1)
{
   double x[1] = {P1(0)};
   double y[1] = {P1(1)};
   double z[1] = {P1(2)};
   double v[9] = {val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2)};
   addDrawable(xGmshDrawable ::TP, 1, x, y, z, v);
}

void xExportGmsh::exportLine(const xPoint &P1, const xPoint &P2, const double &val1, const double &val2)
{
   double x[2] = {P1(0), P2(0)};
   double y[2] = {P1(1), P2(1)};
   double z[2] = {P1(2), P2(2)};
   double v[2] = {val1, val2};
   addDrawable(xGmshDrawable ::SL, 1, x, y, z, v);
}

void xExportGmsh::exportLine(const xPoint &P1, const xPoint &P2, const xtensor::xVector<> &val1, const xtensor::xVector<> &val2)
{
   double x[2] = {P1(0), P2(0)};
   double y[2] = {P1(1), P2(1)};
   double z[2] = {P1(2), P2(2)};
   double v[6] = {val1(0), val1(1), val1(2), val2(0), val2(1), val2(2)};
   addDrawable(xGmshDrawable ::VL, 1, x, y, z, v);
}

void xExportGmsh::exportLine(const xPoint &P1, const xPoint &P2, const xtensor::xTensor2<> &val1, const xtensor::xTensor2<> &val2)
{
   double x[2] = {P1(0), P2(0)};
   double y[2] = {P1(1), P2(1)};
   double z[2] = {P1(2), P2(2)};
   double v[18] = {val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2),
                   val2(0, 0), val2(0, 1), val2(0, 2), val2(1, 0), val2(1, 1), val2(1, 2), val2(2, 0), val2(2, 1), val2(2, 2)};
   addDrawable(xGmshDrawable ::TL, 1, x, y, z, v);
}

void xExportGmsh::exportTriangle(const xPoint &P1, const xPoint &P2, const xPoint &P3, const double &val1, const double &val2,
                                 const double &val3)
{
   double x[3] = {P1(0), P2(0), P3(0)};
   double y[3] = {P1(1), P2(1), P3(1)};
   double z[3] = {P1(2), P2(2), P3(2)};
   double v[3] = {val1, val2, val3};
   addDrawable(xGmshDrawable ::ST, 1, x, y, z, v);
}

void xExportGmsh::exportTriangle(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xtensor::xVector<> &val1,
                                 const xtensor::xVector<> &val2, const xtensor::xVector<> &val3)
{
   double x[3] = {P1(0), P2(0), P3(0)};
   double y[3] = {P1(1), P2(1), P3(1)};
   double z[3] = {P1(2), P2(2), P3(2)};
   double v[9] = {val1(0), val1(1), val1(2), val2(0), val2(1), val2(2), val3(0), val3(1), val3(2)};
   addDrawable(xGmshDrawable ::VT, 1, x, y, z, v);
}

void xExportGmsh::exportTriangle(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xtensor::xTensor2<> &val1,
                                 const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3)
{
   double x[3] = {P1(0), P2(0), P3(0)};
   double y[3] = {P1(1), P2(1), P3(1)};
   double z[3] = {P1(2), P2(2), P3(2)};
   double v[27] = {val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2),
                   val2(0, 0), val2(0, 1), val2(0, 2), val2(1, 0), val2(1, 1), val2(1, 2), val2(2, 0), val2(2, 1), val2(2, 2),
                   val3(0, 0), val3(0, 1), val3(0, 2), val3(1, 0), val3(1, 1), val3(1, 2), val3(2, 0), val3(2, 1), val3(2, 2)};
   addDrawable(xGmshDrawable ::TT, 1, x, y, z, v);
}

void xExportGmsh::exportQuad(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const double &val1,
                             const double &val2, const double &val3, const double &val4)
{
   double x[4] = {P1(0), P2(0), P3(0), P4(0)};
   double y[4] = {P1(1), P2(1), P3(1), P4(1)};
   double z[4] = {P1(2), P2(2), P3(2), P4(2)};
   double v[4] = {val1, val2, val3, val4};
   addDrawable(xGmshDrawable ::SQ, 1, x, y, z, v);
}

void xExportGmsh::exportQuad(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4,
                             const xtensor::xVector<> &val1, const xtensor::xVector<> &val2, const xtensor::xVector<> &val3,
                             const xtensor::xVector<> &val4)
{
   double x[4] = {P1(0), P2(0), P3(0), P4(0)};
   double y[4] = {P1(1), P2(1), P3(1), P4(1)};
   double z[4] = {P1(2), P2(2), P3(2), P4(2)};
   double v[12] = {val1(0), val1(1), val1(2), val2(0), val2(1), val2(2), val3(0), val3(1), val3(2), val4(0), val4(1), val4(2)};
   addDrawable(xGmshDrawable ::VQ, 1, x, y, z, v);
}

void xExportGmsh::exportQuad(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4,
                             const xtensor::xTensor2<> &val1, const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3,
                             const xtensor::xTensor2<> &val4)
{
   double x[4] = {P1(0), P2(0), P3(0), P4(0)};
   double y[4] = {P1(1), P2(1), P3(1), P4(1)};
   double z[4] = {P1(2), P2(2), P3(2), P4(2)};
   double v[36] = {val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2),
                   val2(0, 0), val2(0, 1), val2(0, 2), val2(1, 0), val2(1, 1), val2(1, 2), val2(2, 0), val2(2, 1), val2(2, 2),
                   val3(0, 0), val3(0, 1), val3(0, 2), val3(1, 0), val3(1, 1), val3(1, 2), val3(2, 0), val3(2, 1), val3(2, 2),
                   val4(0, 0), val4(0, 1), val4(0, 2), val4(1, 0), val4(1, 1), val4(1, 2), val4(2, 0), val4(2, 1), val4(2, 2)};
   addDrawable(xGmshDrawable ::TQ, 1, x, y, z, v);
}

void xExportGmsh::exportTetra(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const double &val1,
                              const double &val2, const double &val3, const double &val4)
{
   double x[4] = {P1(0), P2(0), P3(0), P4(0)};
   double y[4] = {P1(1), P2(1), P3(1), P4(1)};
   double z[4] = {P1(2), P2(2), P3(2), P4(2)};
   double v[4] = {val1, val2, val3, val4};
   addDrawable(xGmshDrawable ::SS, 1, x, y, z, v);
}

void xExportGmsh::exportTetra(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4,
                              const xtensor::xVector<> &val1, const xtensor::xVector<> &val2, const xtensor::xVector<> &val3,
                              const xtensor::xVector<> &val4)
{
   double x[4] = {P1(0), P2(0), P3(0), P4(0)};
   double y[4] = {P1(1), P2(1), P3(1), P4(1)};
   double z[4] = {P1(2), P2(2), P3(2), P4(2)};
   double v[12] = {val1(0), val1(1), val1(2), val2(0), val2(1), val2(2), val3(0), val3(1), val3(2), val4(0), val4(1), val4(2)};
   addDrawable(xGmshDrawable ::VS, 1, x, y, z, v);
}

void xExportGmsh::exportTetra(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4,
                              const xtensor::xTensor2<> &val1, const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3,
                              const xtensor::xTensor2<> &val4)
{
   double x[4] = {P1(0), P2(0), P3(0), P4(0)};
   double y[4] = {P1(1), P2(1), P3(1), P4(1)};
   double z[4] = {P1(2), P2(2), P3(2), P4(2)};
   double v[36] = {val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2),
                   val2(0, 0), val2(0, 1), val2(0, 2), val2(1, 0), val2(1, 1), val2(1, 2), val2(2, 0), val2(2, 1), val2(2, 2),
                   val3(0, 0), val3(0, 1), val3(0, 2), val3(1, 0), val3(1, 1), val3(1, 2), val3(2, 0), val3(2, 1), val3(2, 2),
                   val4(0, 0), val4(0, 1), val4(0, 2), val4(1, 0), val4(1, 1), val4(1, 2), val4(2, 0), val4(2, 1), val4(2, 2)};
   addDrawable(xGmshDrawable ::TS, 1, x, y, z, v);
}

void xExportGmsh::exportHex(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xPoint &P5,
                            const xPoint &P6, const xPoint &P7, const xPoint &P8, const double &val1, const double &val2,
                            const double &val3, const double &val4, const double &val5, const double &val6, const double &val7,
                            const double &val8)
{
   double x[8] = {P1(0), P2(0), P3(0), P4(0), P5(0), P6(0), P7(0), P8(0)};
   double y[8] = {P1(1), P2(1), P3(1), P4(1), P5(1), P6(1), P7(1), P8(1)};
   double z[8] = {P1(2), P2(2), P3(2), P4(2), P5(2), P6(2), P7(2), P8(2)};
   double v[8] = {val1, val2, val3, val4, val5, val6, val7, val8};
   addDrawable(xGmshDrawable ::SH, 1, x, y, z, v);
}

void xExportGmsh::exportHex(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xPoint &P5,
                            const xPoint &P6, const xPoint &P7, const xPoint &P8, const xtensor::xVector<> &val1,
                            const xtensor::xVector<> &val2, const xtensor::xVector<> &val3, const xtensor::xVector<> &val4,
                            const xtensor::xVector<> &val5, const xtensor::xVector<> &val6, const xtensor::xVector<> &val7,
                            const xtensor::xVector<> &val8)
{
   double x[8] = {P1(0), P2(0), P3(0), P4(0), P5(0), P6(0), P7(0), P8(0)};
   double y[8] = {P1(1), P2(1), P3(1), P4(1), P5(1), P6(1), P7(1), P8(1)};
   double z[8] = {P1(2), P2(2), P3(2), P4(2), P5(2), P6(2), P7(2), P8(2)};
   double v[24] = {val1(0), val1(1), val1(2), val2(0), val2(1), val2(2), val3(0), val3(1), val3(2), val4(0), val4(1), val4(2),
                   val5(0), val5(1), val5(2), val6(0), val6(1), val6(2), val7(0), val7(1), val7(2), val8(0), val8(1), val8(2)};
   addDrawable(xGmshDrawable ::VH, 1, x, y, z, v);
}

void xExportGmsh::exportHex(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xPoint &P5,
                            const xPoint &P6, const xPoint &P7, const xPoint &P8, const xtensor::xTensor2<> &val1,
                            const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3, const xtensor::xTensor2<> &val4,
                            const xtensor::xTensor2<> &val5, const xtensor::xTensor2<> &val6, const xtensor::xTensor2<> &val7,
                            const xtensor::xTensor2<> &val8)
{
   double x[8] = {P1(0), P2(0), P3(0), P4(0), P5(0), P6(0), P7(0), P8(0)};
   double y[8] = {P1(1), P2(1), P3(1), P4(1), P5(1), P6(1), P7(1), P8(1)};
   double z[8] = {P1(2), P2(2), P3(2), P4(2), P5(2), P6(2), P7(2), P8(2)};
   double v[72] = {val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2),
                   val2(0, 0), val2(0, 1), val2(0, 2), val2(1, 0), val2(1, 1), val2(1, 2), val2(2, 0), val2(2, 1), val2(2, 2),
                   val3(0, 0), val3(0, 1), val3(0, 2), val3(1, 0), val3(1, 1), val3(1, 2), val3(2, 0), val3(2, 1), val3(2, 2),
                   val4(0, 0), val4(0, 1), val4(0, 2), val4(1, 0), val4(1, 1), val4(1, 2), val4(2, 0), val4(2, 1), val4(2, 2),
                   val5(0, 0), val5(0, 1), val5(0, 2), val5(1, 0), val5(1, 1), val5(1, 2), val5(2, 0), val5(2, 1), val5(2, 2),
                   val6(0, 0), val6(0, 1), val6(0, 2), val6(1, 0), val6(1, 1), val6(1, 2), val6(2, 0), val6(2, 1), val6(2, 2),
                   val7(0, 0), val7(0, 1), val7(0, 2), val7(1, 0), val7(1, 1), val7(1, 2), val7(2, 0), val7(2, 1), val7(2, 2),
                   val8(0, 0), val8(0, 1), val8(0, 2), val8(1, 0), val8(1, 1), val8(1, 2), val8(2, 0), val8(2, 1), val8(2, 2)};
   addDrawable(xGmshDrawable ::TH, 1, x, y, z, v);
}

xExportGmshAscii::xExportGmshAscii(MPI_Comm world) : xExportGmsh(world) {}
void xExportGmshAscii::openFile(const string &fName)
{
   string loc = getFileNamePrefix() + fName + getFileNameExtension();
   ofs = new std::ofstream(loc.c_str());
   process_started = true;
}

void xExportGmshAscii::startView(const string &commentar) { (*ofs) << "View \"" << commentar << "\" {" << std::endl; }

void xExportGmshAscii::addDrawable(xGmshDrawable::gType t, int n, double *x, double *y, double *z, double *v)
{
   xGmshDrawable d(t, n, x, y, z, v);
   (*ofs) << d;
}

void xExportGmshAscii::endView() { (*ofs) << "};" << std::endl; }

void xExportGmshAscii::closeFile()
{
   if (!process_started) return;
   ofs->close();
   delete ofs;
   process_started = false;
}

xExportGmshAscii::~xExportGmshAscii() { closeFile(); }

#ifdef WITH_GZIPOUTPUT
xExportGmshAsciigz::xExportGmshAsciigz(MPI_Comm world) : xExportGmsh(world) { filename_extension += ".gz"; }
void xExportGmshAsciigz::openFile(const string &fName)
{
   // add to the chain gzip compresor filter
   ofsf.push(boost::iostreams::gzip_compressor());

   // open a standard stream
   string loc = getFileNamePrefix() + fName + getFileNameExtension();
   ofs = new std::ofstream(loc.c_str(), ios::out | ios::binary);

   if (!ofs) throw;

   // add the standard stream to the chain
   ofsf.push(*ofs);

   process_started = true;
}

void xExportGmshAsciigz::startView(const string &commentar) { ofsf << "View \"" << commentar << "\" {" << std::endl; }

void xExportGmshAsciigz::addDrawable(xGmshDrawable::gType t, int n, double *x, double *y, double *z, double *v)
{
   xGmshDrawable d(t, n, x, y, z, v);
   ofsf << d;
}

void xExportGmshAsciigz::endView() { ofsf << "};" << std::endl; }

void xExportGmshAsciigz::closeFile()
{
   if (!process_started) return;
   // remove from the chain the standard stream
   ofsf.pop();
   // remove from the chain the gzip compressor
   ofsf.pop();
   // close standard stream
   ofs->close();
   // delete allocation
   delete ofs;
   process_started = false;
}

xExportGmshAsciigz::~xExportGmshAsciigz() { closeFile(); }
#endif

xExportGmshAsciiSort::xExportGmshAsciiSort(MPI_Comm world) : xExportGmsh(world) {}

void xExportGmshAsciiSort::openFile(const string &fName)
{
   string loc = getFileNamePrefix() + fName + getFileNameExtension();
   ofs = new std::ofstream(loc.c_str());
   process_started = true;
   d_container.clear();
}

void xExportGmshAsciiSort::startView(const string &commentar) { (*ofs) << "View \"" << commentar << "\" {" << std::endl; }

void xExportGmshAsciiSort::addDrawable(xGmshDrawable::gType t, int n, double *x, double *y, double *z, double *v)
{
   xGmshDrawable d(t, n, x, y, z, v);
   d_container.insert(d);
}

void xExportGmshAsciiSort::endView()
{
   xGmshDrawableSortContainer::iterator it = d_container.begin();
   xGmshDrawableSortContainer::iterator itend = d_container.end();
   for (; it != itend; ++it)
   {
      (*ofs) << (*it);
   }
   (*ofs) << "};" << std::endl;
   d_container.clear();
}

void xExportGmshAsciiSort::closeFile()
{
   if (!process_started) return;
   ofs->close();
   delete ofs;
   process_started = false;
}

xExportGmshAsciiSort::~xExportGmshAsciiSort() { closeFile(); }

// $PostFormat
// 1.2 file_type data_size
// $EndPostFormat
// $View
// view_name nb_time_steps
// nb_scalar_points nb_vector_points nb_tensor_points
// nb_scalar_lines nb_vector_lines nb_tensor_lines
// nb_scalar_triangles nb_vector_triangles nb_tensor_triangles
// nb_scalar_quadrangles nb_vector_quadrangles nb_tensor_quadrangles
// nb_scalar_tetrahedra nb_vector_tetrahedra nb_tensor_tetrahedra
// nb_scalar_hexahedra nb_vector_hexahedra nb_tensor_hexahedra
// nb_scalar_prisms nb_vector_prisms nb_tensor_prisms
// nb_scalar_pyramids nb_vector_pyramids nb_tensor_pyramids
// nb_text2d nb_text2d_chars nb_text3d nb_text3d_chars
// time_step_values
// < scalar_point_value > ...
// < vector_point_value > ...
// < tensor_point_value > ...
// < scalar_line_value > ...
// < vector_line_value > ...
// < tensor_line_value > ...
// < scalar_triangle_value > ...
// < vector_triangle_value > ...
// < tensor_triangle_value > ...
// < scalar_quadrangle_value > ...
// < vector_quadrangle_value > ...
// < tensor_quadrangle_value > ...
// < scalar_tetrahedron_value > ...
// < vector_tetrahedron_value > ...
// < tensor_tetrahedron_value > ...
// < scalar_hexahedron_value > ...
// < vector_hexahedron_value > ...
// < tensor_hexahedron_value > ...
// < scalar_prism_value > ...
// < vector_prism_value > ...
// < tensor_prism_value > ...
// < scalar_pyramid_value > ...
// < vector_pyramid_value > ...
// < tensor_pyramid_value > ...
// < text2d > ... < text2d_chars > ...
// < text3d > ... < text3d_chars > ...
// $EndView

void xExportGmshBinary::clear()
{
   view_name = "";
   nb_time_steps = 1;
   time_step_values.clear();
   time_step_values.push_back(1.0);

   nb_scalar_points = 0;
   nb_vector_points = 0;
   nb_tensor_points = 0;
   nb_scalar_lines = 0;
   nb_vector_lines = 0;
   nb_tensor_lines = 0;
   nb_scalar_triangles = 0;
   nb_vector_triangles = 0;
   nb_tensor_triangles = 0;
   nb_scalar_quadrangles = 0;
   nb_vector_quadrangles = 0;
   nb_tensor_quadrangles = 0;
   nb_scalar_tetrahedra = 0;
   nb_vector_tetrahedra = 0;
   nb_tensor_tetrahedra = 0;
   nb_scalar_hexahedra = 0;
   nb_vector_hexahedra = 0;
   nb_tensor_hexahedra = 0;
   nb_scalar_prisms = 0;
   nb_vector_prisms = 0;
   nb_tensor_prisms = 0;
   nb_scalar_pyramids = 0;
   nb_vector_pyramids = 0;
   nb_tensor_pyramids = 0;
   nb_text2d = 0;
   nb_text2d_chars = 0;
   nb_text3d = 0;
   nb_text3d_chars = 0;

   scalar_points.clear();
   vector_points.clear();
   tensor_points.clear();
   scalar_lines.clear();
   vector_lines.clear();
   tensor_lines.clear();
   scalar_triangles.clear();
   vector_triangles.clear();
   tensor_triangles.clear();
   scalar_quadrangles.clear();
   vector_quadrangles.clear();
   tensor_quadrangles.clear();
   scalar_tetrahedra.clear();
   vector_tetrahedra.clear();
   tensor_tetrahedra.clear();
   scalar_hexahedra.clear();
   vector_hexahedra.clear();
   tensor_hexahedra.clear();
   scalar_prisms.clear();
   vector_prisms.clear();
   tensor_prisms.clear();
   scalar_pyramids.clear();
   vector_pyramids.clear();
   tensor_pyramids.clear();
}

xExportGmshBinary::xExportGmshBinary(MPI_Comm world) : xExportGmsh(world) { clear(); }

void xExportGmshBinary::openFile(const string &fName)
{
   clear();
   string loc = getFileNamePrefix() + fName + getFileNameExtension();
   file = fopen(loc.c_str(), "w");
   process_started = true;
}

void xExportGmshBinary::startView(const string &commentar)
{
   clear();
   view_name = commentar;
   fprintf(file, "$PostFormat\n");
   fprintf(file, "%g %d %d\n", 1.2, 1, (int)sizeof(double));
   fprintf(file, "$EndPostFormat\n");
   fprintf(file, "$View\n");
}

void xExportGmshBinary::addDrawable(xGmshDrawable::gType t, int n, double *x, double *y, double *z, double *v)
{
   std::vector<double> *vector_ptr;
   int nbNod, sizeInfo;
   int nbSamples = n;
   switch (t)
   {
      case xGmshDrawable::NONE:
         printf("error drawable\n");
         throw 1;
      case xGmshDrawable::SP:
         nb_scalar_points++;
         vector_ptr = &scalar_points;
         nbNod = 1;
         sizeInfo = 1;
         break;
      case xGmshDrawable::SL:
         nb_scalar_lines++;
         vector_ptr = &scalar_lines;
         nbNod = 2;
         sizeInfo = 1;
         break;
      case xGmshDrawable::ST:
         nb_scalar_triangles++;
         vector_ptr = &scalar_triangles;
         nbNod = 3;
         sizeInfo = 1;
         break;
      case xGmshDrawable::SQ:
         nb_scalar_quadrangles++;
         vector_ptr = &scalar_quadrangles;
         nbNod = 4;
         sizeInfo = 1;
         break;
      case xGmshDrawable::SS:
         nb_scalar_tetrahedra++;
         vector_ptr = &scalar_tetrahedra;
         nbNod = 4;
         sizeInfo = 1;
         break;
      case xGmshDrawable::SH:
         nb_scalar_hexahedra++;
         vector_ptr = &scalar_hexahedra;
         nbNod = 8;
         sizeInfo = 1;
         break;
      case xGmshDrawable::VP:
         nb_vector_points++;
         vector_ptr = &vector_points;
         nbNod = 1;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VL:
         nb_vector_lines++;
         vector_ptr = &vector_lines;
         nbNod = 2;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VT:
         nb_vector_triangles++;
         vector_ptr = &vector_triangles;
         nbNod = 3;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VQ:
         nb_vector_quadrangles++;
         vector_ptr = &vector_quadrangles;
         nbNod = 4;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VS:
         nb_vector_tetrahedra++;
         vector_ptr = &vector_tetrahedra;
         nbNod = 4;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VH:
         nb_vector_hexahedra++;
         vector_ptr = &vector_hexahedra;
         nbNod = 8;
         sizeInfo = 3;
         break;
      case xGmshDrawable::TP:
         nb_tensor_points++;
         vector_ptr = &tensor_points;
         nbNod = 1;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TL:
         nb_tensor_lines++;
         vector_ptr = &tensor_lines;
         nbNod = 2;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TT:
         nb_tensor_triangles++;
         vector_ptr = &tensor_triangles;
         nbNod = 3;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TQ:
         nb_tensor_quadrangles++;
         vector_ptr = &tensor_quadrangles;
         nbNod = 4;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TS:
         nb_tensor_tetrahedra++;
         vector_ptr = &tensor_tetrahedra;
         nbNod = 4;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TH:
         nb_tensor_hexahedra++;
         vector_ptr = &tensor_hexahedra;
         nbNod = 8;
         sizeInfo = 9;
         break;
      default:
         throw;
   }
   int nbval = nbNod * nbSamples * sizeInfo;
   vector_ptr->insert(vector_ptr->end(), x, x + nbNod);
   vector_ptr->insert(vector_ptr->end(), y, y + nbNod);
   vector_ptr->insert(vector_ptr->end(), z, z + nbNod);
   vector_ptr->insert(vector_ptr->end(), v, v + nbval);
}

void xExportGmshBinary::endView()
{
   int one = 1;
   fprintf(file,
           "%s %d "
           "%d %d %d "
           "%d %d %d "
           "%d %d %d "
           "%d %d %d "
           "%d %d %d "
           "%d %d %d "
           "%d %d %d "
           "%d %d %d "
           "%d %d %d %d\n",
           view_name.c_str(), nb_time_steps, nb_scalar_points, nb_vector_points, nb_tensor_points, nb_scalar_lines,
           nb_vector_lines, nb_tensor_lines, nb_scalar_triangles, nb_vector_triangles, nb_tensor_triangles, nb_scalar_quadrangles,
           nb_vector_quadrangles, nb_tensor_quadrangles, nb_scalar_tetrahedra, nb_vector_tetrahedra, nb_tensor_tetrahedra,
           nb_scalar_hexahedra, nb_vector_hexahedra, nb_tensor_hexahedra, nb_scalar_prisms, nb_vector_prisms, nb_tensor_prisms,
           nb_scalar_pyramids, nb_vector_pyramids, nb_tensor_pyramids, nb_text2d, nb_text2d_chars, nb_text3d, nb_text3d_chars);
   // ascii long
   fwrite(&one, sizeof(int), 1, file);
   fwrite(&(*time_step_values.begin()), sizeof(double), nb_time_steps, file);
   fwrite(&(*scalar_points.begin()), sizeof(double), scalar_points.size(), file);
   fwrite(&(*vector_points.begin()), sizeof(double), vector_points.size(), file);
   fwrite(&(*tensor_points.begin()), sizeof(double), tensor_points.size(), file);
   fwrite(&(*scalar_lines.begin()), sizeof(double), scalar_lines.size(), file);
   fwrite(&(*vector_lines.begin()), sizeof(double), vector_lines.size(), file);
   fwrite(&(*tensor_lines.begin()), sizeof(double), tensor_lines.size(), file);
   fwrite(&(*scalar_triangles.begin()), sizeof(double), scalar_triangles.size(), file);
   fwrite(&(*vector_triangles.begin()), sizeof(double), vector_triangles.size(), file);
   fwrite(&(*tensor_triangles.begin()), sizeof(double), tensor_triangles.size(), file);
   fwrite(&(*scalar_quadrangles.begin()), sizeof(double), scalar_quadrangles.size(), file);
   fwrite(&(*vector_quadrangles.begin()), sizeof(double), vector_quadrangles.size(), file);
   fwrite(&(*tensor_quadrangles.begin()), sizeof(double), tensor_quadrangles.size(), file);
   fwrite(&(*scalar_tetrahedra.begin()), sizeof(double), scalar_tetrahedra.size(), file);
   fwrite(&(*vector_tetrahedra.begin()), sizeof(double), vector_tetrahedra.size(), file);
   fwrite(&(*tensor_tetrahedra.begin()), sizeof(double), tensor_tetrahedra.size(), file);
   fwrite(&(*scalar_hexahedra.begin()), sizeof(double), scalar_hexahedra.size(), file);
   fwrite(&(*vector_hexahedra.begin()), sizeof(double), vector_hexahedra.size(), file);
   fwrite(&(*tensor_hexahedra.begin()), sizeof(double), tensor_hexahedra.size(), file);
   fwrite(&(*scalar_prisms.begin()), sizeof(double), scalar_prisms.size(), file);
   fwrite(&(*vector_prisms.begin()), sizeof(double), vector_prisms.size(), file);
   fwrite(&(*tensor_prisms.begin()), sizeof(double), tensor_prisms.size(), file);
   fwrite(&(*scalar_pyramids.begin()), sizeof(double), scalar_pyramids.size(), file);
   fwrite(&(*vector_pyramids.begin()), sizeof(double), vector_pyramids.size(), file);
   fwrite(&(*tensor_pyramids.begin()), sizeof(double), tensor_pyramids.size(), file);
   fprintf(file, "\n$EndView\n");
}

void xExportGmshBinary::closeFile()
{
   if (!process_started) return;
   fclose(file);
   process_started = false;
}

xExportGmshBinary::~xExportGmshBinary() { closeFile(); }

void xExportGmshMPIIO::clear()
{
   view_name = "";
   std::fill(&nb_info[1], &nb_info[29], 0);
   std::fill(&sz_info[0], &sz_info[29], 0);
   time_step_values.clear();
   if (proc_id)
   {
      nb_info[info_idx::nb_time_steps] = 0;
   }
   else
   {
      nb_info[info_idx::nb_time_steps] = 1;
      time_step_values.push_back(1.0);
      sz_info[0] = 1;
   }

   scalar_points.clear();
   vector_points.clear();
   tensor_points.clear();
   scalar_lines.clear();
   vector_lines.clear();
   tensor_lines.clear();
   scalar_triangles.clear();
   vector_triangles.clear();
   tensor_triangles.clear();
   scalar_quadrangles.clear();
   vector_quadrangles.clear();
   tensor_quadrangles.clear();
   scalar_tetrahedra.clear();
   vector_tetrahedra.clear();
   tensor_tetrahedra.clear();
   scalar_hexahedra.clear();
   vector_hexahedra.clear();
   tensor_hexahedra.clear();
   scalar_prisms.clear();
   vector_prisms.clear();
   tensor_prisms.clear();
   scalar_pyramids.clear();
   vector_pyramids.clear();
   tensor_pyramids.clear();
}

xExportGmshMPIIO::xExportGmshMPIIO(MPI_Comm world, int nb_files) : xExportGmsh(world)
{
   MPI_Comm_rank(world, &proc_id);
   MPI_Comm_size(world, &nb_proc);
   int file_nb_proc = nb_proc;
   int file_proc_id = 0;
   // check nb_files
   // if less then 2 it is set to 0 which mean one unique file (comunicator is world)
   if (nb_files < 2) nb_files = 0;
   // if user ask for more then one file and there is at most 7 proc it's request is refused and
   // only one file is used
   if (nb_files && nb_proc < 8) nb_files = 0;

   // set communicator
   if (nb_files)
   {
      // To force at least 4 proc per file. This avoid user to use MPIIO like Ascii with one file
      // per proc. If this one file per proc is wanted better use Ascii/binary version.
      int nb_files_limit = nb_proc / 4;
      if (nb_files > nb_files_limit) nb_files = nb_files_limit;

      // color algo:
      // modulo: basic but not clear if it offers best performance
      // depend on computer architecture and job managing tool it may be
      // the worst case (say inter lam communication instead of intra lam communication)
      // auto balanced
      // To be tested
      // file_proc_id = proc_id % nb_files;
      // range: a little more complex but insure, a priori,  intra lam communication
      // (process are packed on one lam then on the other, ....)
      // balanced by use of larger chunk for first files up to shift * (chunck + 1) proc
      {
         const int chunk = nb_proc / nb_files;
         int shift = nb_proc % nb_files;
         int over = 1;
         if (shift && proc_id < shift * (chunk + 1))
         {
            shift = 0;
         }
         else
            over = 0;
         file_proc_id = (proc_id - shift) / (chunk + over);
      }
      // split comm with colors: 1 color per file
      MPI_Comm_split(world, file_proc_id, 0, &univ);

      // reset comunicator according to new comunicator
      MPI_Comm_rank(univ, &proc_id);
      MPI_Comm_size(univ, &nb_proc);
   }
   else
   {
      univ = world;
   }
   assert(file_nb_proc < 10000);
   sprintf(&filename_prefix.at(0), "S_%04d_F_%04d_", file_nb_proc, file_proc_id);
   MPI_Type_size(MPI_DOUBLE, &szd);
   MPI_Type_size(MPI_CHAR, &szc);
   MPI_Type_size(MPI_INT, &szi);
   head = "$PostFormat\n1.2 1 " + std::to_string(szd) + "\n$EndPostFormat\n";
   clear();
}

void xExportGmshMPIIO::openFile(const string &fName)
{
   clear();
   MPI_Status status;
   string loc = getFileNamePrefix() + fName + getFileNameExtension();
   MPI_File_open(univ, loc.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
   MPI_File_set_size(file, 0);
   if (!proc_id) MPI_File_write_at(file, 0, const_cast<char *>(head.c_str()), head.size(), MPI_CHAR, &status);
   process_started = true;
}

void xExportGmshMPIIO::startView(const string &commentar)
{
   clear();
   view_name = commentar;
}

void xExportGmshMPIIO::addDrawable(xGmshDrawable::gType t, int n, double *x, double *y, double *z, double *v)
{
   std::vector<double> *vector_ptr;
   int nbNod, sizeInfo;
   int nbSamples = n;
   info_idx idt;
   switch (t)
   {
      case xGmshDrawable::NONE:
         printf("error drawable\n");
         throw 1;
      case xGmshDrawable::SP:
         idt = info_idx::nb_scalar_points;
         vector_ptr = &scalar_points;
         nbNod = 1;
         sizeInfo = 1;
         break;
      case xGmshDrawable::SL:
         idt = info_idx::nb_scalar_lines;
         vector_ptr = &scalar_lines;
         nbNod = 2;
         sizeInfo = 1;
         break;
      case xGmshDrawable::ST:
         idt = info_idx::nb_scalar_triangles;
         vector_ptr = &scalar_triangles;
         nbNod = 3;
         sizeInfo = 1;
         break;
      case xGmshDrawable::SQ:
         idt = info_idx::nb_scalar_quadrangles;
         vector_ptr = &scalar_quadrangles;
         nbNod = 4;
         sizeInfo = 1;
         break;
      case xGmshDrawable::SS:
         idt = info_idx::nb_scalar_tetrahedra;
         vector_ptr = &scalar_tetrahedra;
         nbNod = 4;
         sizeInfo = 1;
         break;
      case xGmshDrawable::SH:
         idt = info_idx::nb_scalar_hexahedra;
         vector_ptr = &scalar_hexahedra;
         nbNod = 8;
         sizeInfo = 1;
         break;
      case xGmshDrawable::VP:
         idt = info_idx::nb_vector_points;
         vector_ptr = &vector_points;
         nbNod = 1;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VL:
         idt = info_idx::nb_vector_lines;
         vector_ptr = &vector_lines;
         nbNod = 2;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VT:
         idt = info_idx::nb_vector_triangles;
         vector_ptr = &vector_triangles;
         nbNod = 3;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VQ:
         idt = info_idx::nb_vector_quadrangles;
         vector_ptr = &vector_quadrangles;
         nbNod = 4;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VS:
         idt = info_idx::nb_vector_tetrahedra;
         vector_ptr = &vector_tetrahedra;
         nbNod = 4;
         sizeInfo = 3;
         break;
      case xGmshDrawable::VH:
         idt = info_idx::nb_vector_hexahedra;
         vector_ptr = &vector_hexahedra;
         nbNod = 8;
         sizeInfo = 3;
         break;
      case xGmshDrawable::TP:
         idt = info_idx::nb_tensor_points;
         vector_ptr = &tensor_points;
         nbNod = 1;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TL:
         idt = info_idx::nb_tensor_lines;
         vector_ptr = &tensor_lines;
         nbNod = 2;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TT:
         idt = info_idx::nb_tensor_triangles;
         vector_ptr = &tensor_triangles;
         nbNod = 3;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TQ:
         idt = info_idx::nb_tensor_quadrangles;
         vector_ptr = &tensor_quadrangles;
         nbNod = 4;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TS:
         idt = info_idx::nb_tensor_tetrahedra;
         vector_ptr = &tensor_tetrahedra;
         nbNod = 4;
         sizeInfo = 9;
         break;
      case xGmshDrawable::TH:
         idt = info_idx::nb_tensor_hexahedra;
         vector_ptr = &tensor_hexahedra;
         nbNod = 8;
         sizeInfo = 9;
         break;
      default:
         throw;
   }
   int nbval = nbSamples * sizeInfo;
   nb_info[idt]++;
   sz_info[idt] += nbNod * (nbval + 3);
   nbval *= nbNod;
   vector_ptr->insert(vector_ptr->end(), x, x + nbNod);
   vector_ptr->insert(vector_ptr->end(), y, y + nbNod);
   vector_ptr->insert(vector_ptr->end(), z, z + nbNod);
   vector_ptr->insert(vector_ptr->end(), v, v + nbval);
}

void xExportGmshMPIIO::endView()
{
   MPI_Status status;
   int last = nb_proc - 1;
   MPI_Offset offsets[29];
   MPI_Offset distg[30];
   // compute sizes accross all proc
   int sizes[29];
   MPI_Reduce(nb_info, sizes, 29, MPI_INT, MPI_SUM, last, univ);

   // compute offset from local size
   MPI_Scan(sz_info, offsets, 29, MPI_OFFSET, MPI_SUM, univ);

   // last proc of univ got sum across all proc of number of info to store per information category
   if (proc_id == last)
   {
      string nbs = "$View\n" + view_name;
      for (int i = 0; i < 29; ++i)
      {
         distg[i + 1] = offsets[i] * szd;
         nbs += " " + std::to_string(sizes[i]);
      }
      nbs += "\n";
      MPI_Offset current_offset;
      MPI_File_get_size(file, &current_offset);
      // write header + nb line information
      MPI_File_write_at(file, current_offset, const_cast<char *>(nbs.c_str()), nbs.size(), MPI_CHAR, &status);
      current_offset += nbs.size() * szc;
      // write 1
      int one = 1;
      MPI_File_write_at(file, current_offset, &one, 1, MPI_INT, &status);
      distg[0] = current_offset + szi;

      // compute once distance per information category
      for (int i = 1; i < 30; ++i) distg[i] += distg[i - 1];
   }

   // brodcast distance per information category
   MPI_Bcast(distg, 30, MPI_OFFSET, last, univ);

   // compute distance per proc per information category
   MPI_Offset dist[29];
   for (int i = 0; i < 29; ++i)
   {
      dist[i] = (offsets[i] - sz_info[i]) * szd + distg[i];
   }

   MPI_File_write_at_all(file, dist[info_idx::nb_time_steps], time_step_values.data(), nb_info[info_idx::nb_time_steps],
                         MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_scalar_points], scalar_points.data(), scalar_points.size(), MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_vector_points], vector_points.data(), vector_points.size(), MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_tensor_points], tensor_points.data(), tensor_points.size(), MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_scalar_lines], scalar_lines.data(), scalar_lines.size(), MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_vector_lines], vector_lines.data(), vector_lines.size(), MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_tensor_lines], tensor_lines.data(), tensor_lines.size(), MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_scalar_triangles], scalar_triangles.data(), scalar_triangles.size(), MPI_DOUBLE,
                         &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_vector_triangles], vector_triangles.data(), vector_triangles.size(), MPI_DOUBLE,
                         &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_tensor_triangles], tensor_triangles.data(), tensor_triangles.size(), MPI_DOUBLE,
                         &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_scalar_quadrangles], scalar_quadrangles.data(), scalar_quadrangles.size(),
                         MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_vector_quadrangles], vector_quadrangles.data(), vector_quadrangles.size(),
                         MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_tensor_quadrangles], tensor_quadrangles.data(), tensor_quadrangles.size(),
                         MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_scalar_tetrahedra], scalar_tetrahedra.data(), scalar_tetrahedra.size(),
                         MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_vector_tetrahedra], vector_tetrahedra.data(), vector_tetrahedra.size(),
                         MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_tensor_tetrahedra], tensor_tetrahedra.data(), tensor_tetrahedra.size(),
                         MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_scalar_hexahedra], scalar_hexahedra.data(), scalar_hexahedra.size(), MPI_DOUBLE,
                         &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_vector_hexahedra], vector_hexahedra.data(), vector_hexahedra.size(), MPI_DOUBLE,
                         &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_tensor_hexahedra], tensor_hexahedra.data(), tensor_hexahedra.size(), MPI_DOUBLE,
                         &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_scalar_prisms], scalar_prisms.data(), scalar_prisms.size(), MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_vector_prisms], vector_prisms.data(), vector_prisms.size(), MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_tensor_prisms], tensor_prisms.data(), tensor_prisms.size(), MPI_DOUBLE, &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_scalar_pyramids], scalar_pyramids.data(), scalar_pyramids.size(), MPI_DOUBLE,
                         &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_vector_pyramids], vector_pyramids.data(), vector_pyramids.size(), MPI_DOUBLE,
                         &status);
   MPI_File_write_at_all(file, dist[info_idx::nb_tensor_pyramids], tensor_pyramids.data(), tensor_pyramids.size(), MPI_DOUBLE,
                         &status);

   if (!proc_id) MPI_File_write_at(file, distg[29], const_cast<char *>("\n$EndView\n"), 10, MPI_CHAR, &status);
}

void xExportGmshMPIIO::closeFile()
{
   if (!process_started) return;
   MPI_File_close(&file);
   process_started = false;
}

xExportGmshMPIIO::~xExportGmshMPIIO() { closeFile(); }

xGmshDrawable::xGmshDrawable() : nbNod(0), sizeInfo(0), X(nullptr), Y(nullptr), Z(nullptr), V(nullptr), nbSamples(0), type(NONE)
{
}
xGmshDrawable::xGmshDrawable(const xGmshDrawable &in)
    : nbNod(in.nbNod), sizeInfo(in.sizeInfo), nbSamples(in.nbSamples), type(in.type)
{
   const int s_v = nbNod * nbSamples * sizeInfo;
   X = new double[nbNod];
   Y = new double[nbNod];
   Z = new double[nbNod];
   V = new double[s_v];
   memcpy(X, in.X, nbNod * sizeof(double));
   memcpy(Y, in.Y, nbNod * sizeof(double));
   memcpy(Z, in.Z, nbNod * sizeof(double));
   memcpy(V, in.V, s_v * sizeof(double));
}

xGmshDrawable::xGmshDrawable(gType t, int n, double *x, double *y, double *z, double *v) : nbSamples(n), type(t)
{
   switch (type)
   {
      case NONE:
         printf("error drawable\n");
         throw 1;
      case SP:
         nbNod = 1;
         sizeInfo = 1;
         break;
      case SL:
         nbNod = 2;
         sizeInfo = 1;
         break;
      case ST:
         nbNod = 3;
         sizeInfo = 1;
         break;
      case SS:
         nbNod = 4;
         sizeInfo = 1;
         break;
      case SQ:
         nbNod = 4;
         sizeInfo = 1;
         break;
      case SH:
         nbNod = 8;
         sizeInfo = 1;
         break;
      case VP:
         nbNod = 1;
         sizeInfo = 3;
         break;
      case VL:
         nbNod = 2;
         sizeInfo = 3;
         break;
      case VT:
         nbNod = 3;
         sizeInfo = 3;
         break;
      case VS:
         nbNod = 4;
         sizeInfo = 3;
         break;
      case VQ:
         nbNod = 4;
         sizeInfo = 3;
         break;
      case VH:
         nbNod = 8;
         sizeInfo = 3;
         break;
      case TP:
         nbNod = 1;
         sizeInfo = 9;
         break;
      case TL:
         nbNod = 2;
         sizeInfo = 9;
         break;
      case TT:
         nbNod = 3;
         sizeInfo = 9;
         break;
      case TS:
         nbNod = 4;
         sizeInfo = 9;
         break;
      case TQ:
         nbNod = 4;
         sizeInfo = 9;
         break;
      case TH:
         nbNod = 8;
         sizeInfo = 9;
         break;
      default:
         assert(0);
         break;
   }
   X = new double[nbNod];
   Y = new double[nbNod];
   Z = new double[nbNod];
   V = new double[nbNod * nbSamples * sizeInfo];
   for (int i = 0; i < nbNod; i++)
   {
      X[i] = x[i];
      Y[i] = y[i];
      Z[i] = z[i];
   }
   memcpy(V, v, nbNod * nbSamples * sizeInfo * sizeof(double));
}

xGmshDrawable::~xGmshDrawable()
{
   if (X) delete[] X;
   if (Y) delete[] Y;
   if (Z) delete[] Z;
   if (V) delete[] V;
}

/*
  std::ostream & operator << (std::ostream & ofs, const xGmshDrawable &p)
  {
  switch(p.type)
  {
  case xGmshDrawable::NONE:return ofs;
  case xGmshDrawable::SP  :ofs << "SP (";break;
  case xGmshDrawable::SL  :ofs << "SL (";break;
  case xGmshDrawable::ST  :ofs << "ST (";break;
  case xGmshDrawable::SS  :ofs << "SS (";break;
  case xGmshDrawable::SQ  :ofs << "SQ (";break;
  case xGmshDrawable::SH  :ofs << "SH (";break;
  case xGmshDrawable::VP  :ofs << "VP (";break;
  case xGmshDrawable::VL  :ofs << "VL (";break;
  case xGmshDrawable::VT  :ofs << "VT (";break;
  case xGmshDrawable::VS  :ofs << "VS (";break;
  case xGmshDrawable::VQ  :ofs << "VQ (";break;
  case xGmshDrawable::VH  :ofs << "VH (";break;
  case xGmshDrawable::TP  :ofs << "TL (";break;
  case xGmshDrawable::TL  :ofs << "TL (";break;
  case xGmshDrawable::TT  :ofs << "TT (";break;
  case xGmshDrawable::TS  :ofs << "TS (";break;
  case xGmshDrawable::TQ  :ofs << "TQ (";break;
  case xGmshDrawable::TH  :ofs << "TH (";break;
  default: assert(0); break;
  }

  for(int i=0;i<p.nbNod;i++)
  {
  if(i)ofs << "," << p.X[i] << "," << p.Y[i] << "," << p.Z[i];
  else ofs << p.X[i] << "," << p.Y[i] << "," << p.Z[i];
  }
  ofs << "){";
  for(int i=0;i<p.nbNod*p.nbSamples*p.sizeInfo;i++)
  {
  if(i) ofs << "," << p.V[i];
  else ofs << p.V[i];
  }
  ofs << "};\n";

  return ofs;
  }
*/

bool xGmshDrawable::lessthan(const xGmshDrawable &other) const
{
   // differente type give the first level of sorting
   if (type < other.type) return true;
   if (type > other.type) return false;

   // local
   int i, n = nbNod;

   // if of same type second level of sorting on nodes
   for (i = 0; i < n; ++i)
   {
      if (X[i] < other.X[i]) return true;
      if (X[i] > other.X[i]) return false;
      if (Y[i] < other.Y[i]) return true;
      if (Y[i] > other.Y[i]) return false;
      if (Z[i] < other.Z[i]) return true;
      if (Z[i] > other.Z[i]) return false;
   }

   // if of same type with same nodes third level of sorting on values
   n *= nbSamples * sizeInfo;
   for (i = 0; i < n; ++i)
   {
      if (V[i] < other.V[i]) return true;
      if (V[i] > other.V[i]) return false;
   }

   // exactely the same
   return false;
}

bool xGmshSortDrawable::operator()(const xGmshDrawable &d1, const xGmshDrawable &d2) const { return d1.lessthan(d2); }

}  // namespace xexport
