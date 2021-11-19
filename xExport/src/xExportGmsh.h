/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef _export_gmsh_H
#define _export_gmsh_H

#include <iostream>
#ifdef WITH_GZIPOUTPUT
#include <boost/iostreams/filtering_stream.hpp>
#endif
#include "xExport.h"
#include "xPoint.h"

namespace xexport
{
using xtensor::xPoint;

class xGmshDrawable
{
  public:
   enum gType
   {
      SP,
      SL,
      ST,
      SS,
      SQ,
      SH,
      VP,
      VL,
      VT,
      VS,
      VQ,
      VH,
      TP,
      TL,
      TT,
      TS,
      TQ,
      TH,
      NONE
   };

  private:
   int nbNod;
   int sizeInfo;
   double *X, *Y, *Z;
   double *V;
   int nbSamples;
   gType type;

  public:
   xGmshDrawable();
   xGmshDrawable(gType, int, double *, double *, double *, double *);
   xGmshDrawable(const xGmshDrawable &);
   virtual ~xGmshDrawable();
   template <typename O>
   friend O &operator<<(O &ofs, const xGmshDrawable &p);
   inline int getSizeInfo() const { return sizeInfo; }
   inline int getNbNod() const { return nbNod; }
   inline int getNbSamples() const { return nbSamples; }
   inline double getV(int i) const { return V[i]; }
   inline double getX(int i) const { return X[i]; }
   inline double getY(int i) const { return Y[i]; }
   inline double getZ(int i) const { return Z[i]; }
   bool lessthan(const xGmshDrawable &) const;
};

class xGmshSortDrawable
{
  public:
   bool operator()(const xGmshDrawable &d1, const xGmshDrawable &d2) const;
};

class xExportGmsh : public xExport
{
  public:
   xExportGmsh(MPI_Comm world_);
   void exportPoint(const xPoint &P1, const double &val1) override;
   void exportPoint(const xPoint &P1, const xtensor::xVector<> &val1) override;
   void exportPoint(const xPoint &P1, const xtensor::xTensor2<> &val1) override;

   void exportLine(const xPoint &P1, const xPoint &P2, const double &val1, const double &val2) override;
   void exportLine(const xPoint &P1, const xPoint &P2, const xtensor::xVector<> &val1, const xtensor::xVector<> &val2) override;
   void exportLine(const xPoint &P1, const xPoint &P2, const xtensor::xTensor2<> &val1, const xtensor::xTensor2<> &val2) override;

   void exportTriangle(const xPoint &P1, const xPoint &P2, const xPoint &P3, const double &val1, const double &val2,
                       const double &val3) override;
   void exportTriangle(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xtensor::xVector<> &val1,
                       const xtensor::xVector<> &val2, const xtensor::xVector<> &val3) override;
   void exportTriangle(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xtensor::xTensor2<> &val1,
                       const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3) override;

   void exportQuad(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const double &val1, const double &val2,
                   const double &val3, const double &val4) override;
   void exportQuad(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xtensor::xVector<> &val1,
                   const xtensor::xVector<> &val2, const xtensor::xVector<> &val3, const xtensor::xVector<> &val4) override;
   void exportQuad(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xtensor::xTensor2<> &val1,
                   const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3, const xtensor::xTensor2<> &val4) override;

   void exportTetra(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const double &val1,
                    const double &val2, const double &val3, const double &val4) override;
   void exportTetra(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xtensor::xVector<> &val1,
                    const xtensor::xVector<> &val2, const xtensor::xVector<> &val3, const xtensor::xVector<> &val4) override;
   void exportTetra(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xtensor::xTensor2<> &val1,
                    const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3, const xtensor::xTensor2<> &val4) override;

   void exportHex(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xPoint &P5, const xPoint &P6,
                  const xPoint &P7, const xPoint &P8, const double &val1, const double &val2, const double &val3,
                  const double &val4, const double &val5, const double &val6, const double &val7, const double &val8) override;
   void exportHex(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xPoint &P5, const xPoint &P6,
                  const xPoint &P7, const xPoint &P8, const xtensor::xVector<> &val1, const xtensor::xVector<> &val2,
                  const xtensor::xVector<> &val3, const xtensor::xVector<> &val4, const xtensor::xVector<> &val5,
                  const xtensor::xVector<> &val6, const xtensor::xVector<> &val7, const xtensor::xVector<> &val8) override;
   void exportHex(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xPoint &P5, const xPoint &P6,
                  const xPoint &P7, const xPoint &P8, const xtensor::xTensor2<> &val1, const xtensor::xTensor2<> &val2,
                  const xtensor::xTensor2<> &val3, const xtensor::xTensor2<> &val4, const xtensor::xTensor2<> &val5,
                  const xtensor::xTensor2<> &val6, const xtensor::xTensor2<> &val7, const xtensor::xTensor2<> &val8) override;

  protected:
   virtual void addDrawable(xGmshDrawable::gType, int, double *, double *, double *, double *) = 0;
};

class xExportGmshAscii : public xExportGmsh
{
  public:
   xExportGmshAscii(MPI_Comm world = MPI_COMM_WORLD);
   ~xExportGmshAscii() override;
   void startView(const std::string &comment) override;
   void addDrawable(xGmshDrawable::gType, int, double *, double *, double *, double *) override;  //
   void endView() override;
   void openFile(const std::string &fName) override;
   void closeFile() override;

  protected:
   std::ofstream *ofs;
};

#ifdef WITH_GZIPOUTPUT
class xExportGmshAsciigz : public xExportGmsh
{
  public:
   xExportGmshAsciigz(MPI_Comm world = MPI_COMM_WORLD);
   virtual ~xExportGmshAsciigz();
   void startView(const std::string &comment);
   virtual void addDrawable(xGmshDrawable::gType, int, double *, double *, double *, double *);  //
   void endView();
   void openFile(const std::string &fName);
   void closeFile();

  protected:
   boost::iostreams::filtering_ostream ofsf;
   std::ofstream *ofs;
};
#endif

class xExportGmshBinary : public xExportGmsh
{
  public:
   xExportGmshBinary(MPI_Comm world = MPI_COMM_WORLD);
   ~xExportGmshBinary() override;
   void startView(const std::string &comment) override;
   void addDrawable(xGmshDrawable::gType, int, double *, double *, double *, double *) override;
   void endView() override;
   void openFile(const std::string &fName) override;
   void closeFile() override;

  private:
   void clear();
   FILE *file;
   std::string view_name;
   int nb_time_steps, nb_scalar_points, nb_vector_points, nb_tensor_points, nb_scalar_lines, nb_vector_lines, nb_tensor_lines,
       nb_scalar_triangles, nb_vector_triangles, nb_tensor_triangles, nb_scalar_quadrangles, nb_vector_quadrangles,
       nb_tensor_quadrangles, nb_scalar_tetrahedra, nb_vector_tetrahedra, nb_tensor_tetrahedra, nb_scalar_hexahedra,
       nb_vector_hexahedra, nb_tensor_hexahedra, nb_scalar_prisms, nb_vector_prisms, nb_tensor_prisms, nb_scalar_pyramids,
       nb_vector_pyramids, nb_tensor_pyramids, nb_text2d, nb_text2d_chars, nb_text3d, nb_text3d_chars;

   std::vector<double> scalar_points, vector_points, tensor_points, scalar_lines, vector_lines, tensor_lines, scalar_triangles,
       vector_triangles, tensor_triangles, scalar_quadrangles, vector_quadrangles, tensor_quadrangles, scalar_tetrahedra,
       vector_tetrahedra, tensor_tetrahedra, scalar_hexahedra, vector_hexahedra, tensor_hexahedra, scalar_prisms, vector_prisms,
       tensor_prisms, scalar_pyramids, vector_pyramids, tensor_pyramids, time_step_values;
};

/*! GMSH export class that use MPIIO functionality: parallel i/o
 *
 * The constructor have 2 defaulted arguments:
 *   - the communicator witch is involved in the export. It could be anything you want. By default it
 *     is MPI_COMM_WORLD.
 *   - the number of file to generate. By default it is zero witch means that only one file will be created to
 *     collect all information generated by all world communicator process. A value of one is treated the same way
 *     as zero. Any other value greater or equal to 2 will be used such as a file is associated to 4 processes
 *     at least. Thus asking for a large number of files will always be limited by the number of process in the given
 *     communicator divide by 4:
 *        nb proc<8 => max number of file is 1
 *        nb proc=110 => max number of file is 27:
 *                    + If user have set nb_files to 20 then the first 10 files will hold export of 6 processes and
 *                      the 10 last files will hold export of 5 processes. All files use process with ascending id (i.e.
 *                      process from 0 to 5 will be hold by first file, 6 to 11 by second file, ....)
 *                    + If user have set nb_files to 40 then it is reset to 27. The first 2 files will hold export of
 *                      5 processes and the 25 last files will hold export of 4 processes.
 *
 * The constraint imposed by the class that a file must collect information at least from 4 processes avoid using
 * MPIIO for to few communication and avoid the case of 1 process witch can be more efficiently treated using
 * xExportGmshBinary or xExportGmshAscii.
 *
 * User can also fully control the number of files he generates with his export by providing a appropriate communicator,
 * setting nb_files and set appropriate file name. For example if using MPI_COMM_SELF then the export will create
 * one file per process .... which is not recommended.
 *
 * All methods are collective (except addDrawable). Thus if instance of this class are used with Export function, this
 * function must be called by all processes of the instance communicator.
 *
 */
class xExportGmshMPIIO : public xExportGmsh
{
  public:
   xExportGmshMPIIO(MPI_Comm world = MPI_COMM_WORLD, int nb_files = 0);
   ~xExportGmshMPIIO() override;
   void startView(const std::string &comment) override;
   void addDrawable(xGmshDrawable::gType, int, double *, double *, double *, double *) override;
   void endView() override;
   void openFile(const std::string &fName) override;
   void closeFile() override;

  private:
   void clear();
   MPI_Comm univ;
   MPI_File file;
   std::string view_name, head;
   int proc_id, nb_proc, szd, szc, szi;
   enum info_idx
   {
      nb_time_steps = 0,
      nb_scalar_points,
      nb_vector_points,
      nb_tensor_points,
      nb_scalar_lines,
      nb_vector_lines,
      nb_tensor_lines,
      nb_scalar_triangles,
      nb_vector_triangles,
      nb_tensor_triangles,
      nb_scalar_quadrangles,
      nb_vector_quadrangles,
      nb_tensor_quadrangles,
      nb_scalar_tetrahedra,
      nb_vector_tetrahedra,
      nb_tensor_tetrahedra,
      nb_scalar_hexahedra,
      nb_vector_hexahedra,
      nb_tensor_hexahedra,
      nb_scalar_prisms,
      nb_vector_prisms,
      nb_tensor_prisms,
      nb_scalar_pyramids,
      nb_vector_pyramids,
      nb_tensor_pyramids,
      nb_text2d,
      nb_text2d_chars,
      nb_text3d,
      nb_text3d_chars
   };
   int nb_info[29];
   MPI_Offset sz_info[29];

   std::vector<double> scalar_points, vector_points, tensor_points, scalar_lines, vector_lines, tensor_lines, scalar_triangles,
       vector_triangles, tensor_triangles, scalar_quadrangles, vector_quadrangles, tensor_quadrangles, scalar_tetrahedra,
       vector_tetrahedra, tensor_tetrahedra, scalar_hexahedra, vector_hexahedra, tensor_hexahedra, scalar_prisms, vector_prisms,
       tensor_prisms, scalar_pyramids, vector_pyramids, tensor_pyramids, time_step_values;
};

class xExportGmshAsciiSort : public xExportGmsh
{
  public:
   xExportGmshAsciiSort(MPI_Comm world = MPI_COMM_WORLD);
   ~xExportGmshAsciiSort() override;
   void startView(const std::string &comment) override;
   void addDrawable(xGmshDrawable::gType, int, double *, double *, double *, double *) override;  //
   void endView() override;
   void openFile(const std::string &fName) override;
   void closeFile() override;

  protected:
   std::ofstream *ofs;
   typedef std::multiset<xGmshDrawable, xGmshSortDrawable> xGmshDrawableSortContainer;
   xGmshDrawableSortContainer d_container;
};

#include "xExportGmsh_imp.h"

}  // namespace xexport

#endif
