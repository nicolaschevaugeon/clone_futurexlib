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
#include "mPoint.h"

namespace xfem
{

  using Trellis_Util::mPoint;




class xGmshDrawable
{
public :
  enum gType {SP,SL,ST,SS,SQ,SH,VP,VL,VT,VS,VQ,VH,TP,TL,TT,TS,TQ,TH,NONE};
private:
  int nbNod;
  int sizeInfo;
  double *X,*Y,*Z;
  double *V;
  int nbSamples;
  gType type;
public :
  xGmshDrawable();
  xGmshDrawable(gType, int, double *, double *, double *, double *);
  xGmshDrawable(const xGmshDrawable&);
  virtual ~xGmshDrawable(); 
  template <typename O> friend O & operator << (O & ofs, const xGmshDrawable &p);
  inline int getSizeInfo() const {return sizeInfo;}
  inline int getNbNod() const {return nbNod;}
  inline int getNbSamples() const {return nbSamples;}
  inline double getV(int i) const {return V[i];}
  inline double getX(int i) const {return X[i];}
  inline double getY(int i) const {return Y[i];}
  inline double getZ(int i) const {return Z[i];}
  bool lessthan (const xGmshDrawable &) const;

};

class  xGmshSortDrawable
{
public:
     bool operator()(const xGmshDrawable& d1, const xGmshDrawable& d2) const;
};

class xexport::xExportGmsh : public xexport::xExport
{

public:
  xexport::xExportGmsh();
  void exportPoint (const Trellis_Util::mPoint & P1, const double   & val1);
  void exportPoint (const Trellis_Util::mPoint & P1, const xtensor::xVector  & val1);
  void exportPoint (const Trellis_Util::mPoint & P1, const xtensor::xTensor2 & val1);

  void exportLine (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const double & val1, const double & val2);
  void exportLine (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const xtensor::xVector & val1, const xtensor::xVector & val2);
  void exportLine (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2);

  void exportTriangle (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, 
		       const double & val1, const double & val2, const double & val3 ); 
  void exportTriangle (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, 
		       const xtensor::xVector & val1, const xtensor::xVector & val2, const xtensor::xVector & val3 );
  void exportTriangle (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, 
		       const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2, const xtensor::xTensor2 & val3 );

  void exportQuad (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
		   const double & val1, const double & val2, const double & val3,  const double & val4); 
  void exportQuad (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
		   const xtensor::xVector & val1, const xtensor::xVector & val2, const xtensor::xVector & val3, const xtensor::xVector & val4 );
  void exportQuad (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
		   const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2, const xtensor::xTensor2 & val3,  const xtensor::xTensor2 & val4);

  void exportTetra (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
		    const double & val1, const double & val2, const double & val3, const double & val4 );
  void exportTetra (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
		    const xtensor::xVector & val1, const xtensor::xVector & val2, const xtensor::xVector & val3, const xtensor::xVector & val4);
  void exportTetra (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
		    const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2, const xtensor::xTensor2 & val3, const xtensor::xTensor2 & val4 );

  void exportHex (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4,
			  const Trellis_Util::mPoint & P5, const Trellis_Util::mPoint & P6, const Trellis_Util::mPoint & P7, const Trellis_Util::mPoint & P8, 
			    const double & val1, const double & val2, const double & val3, const double & val4, const double & val5, const double & val6, const double & val7, const double & val8 ) ;
  void exportHex (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4,
			  const Trellis_Util::mPoint & P5, const Trellis_Util::mPoint & P6, const Trellis_Util::mPoint & P7, const Trellis_Util::mPoint & P8, 
			    const xtensor::xVector & val1, const xtensor::xVector & val2, const xtensor::xVector & val3, const xtensor::xVector & val4, const xtensor::xVector & val5, const xtensor::xVector & val6, const xtensor::xVector & val7, const xtensor::xVector & val8 ) ;
  void exportHex (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4,
			  const Trellis_Util::mPoint & P5, const Trellis_Util::mPoint & P6, const Trellis_Util::mPoint & P7, const Trellis_Util::mPoint & P8, 
			    const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2, const xtensor::xTensor2 & val3, const xtensor::xTensor2 & val4, const xtensor::xTensor2 & val5, const xtensor::xTensor2 & val6, const xtensor::xTensor2 & val7, const xtensor::xTensor2 & val8 ) ;

protected:
  virtual void addDrawable  (xGmshDrawable::gType, int, double *, double *, double *, double *) = 0;       

};


class xexport::xExportGmshAscii : public xexport::xExportGmsh
{
public:
  
  virtual ~xexport::xExportGmshAscii();
  void startView     (const string& comment);
  virtual void addDrawable  (xGmshDrawable::gType, int, double *, double *, double *, double *);	// 
  void endView ();
  void openFile (const string& fName);
  void closeFile ();
protected:
  std::ofstream *ofs;
  
};

#ifdef WITH_GZIPOUTPUT 
class xexport::xExportGmshAsciigz : public xexport::xExportGmsh
{
public:
  
  virtual ~xexport::xExportGmshAsciigz();
  void startView     (const string& comment);
  virtual void addDrawable  (xGmshDrawable::gType, int, double *, double *, double *, double *);	// 
  void endView ();
  void openFile (const string& fName);
  void closeFile ();
protected:
  boost::iostreams::filtering_ostream ofsf;
  std::ofstream *ofs;
  
};
#endif


class xexport::xExportGmshBinary : public xexport::xExportGmsh
{

public:
  xexport::xExportGmshBinary();
  virtual ~xexport::xExportGmshBinary();
  void startView    (const string& comment);
  virtual void addDrawable  (xGmshDrawable::gType, int, double *, double *, double *, double *);
  void endView ();
  void openFile (const string& fName);
  void closeFile ();

private:
  void clear();
  FILE * file;
  std::string view_name; 
  int nb_time_steps,
  nb_scalar_points, nb_vector_points, nb_tensor_points,
  nb_scalar_lines, nb_vector_lines, nb_tensor_lines,
  nb_scalar_triangles, nb_vector_triangles, nb_tensor_triangles,
  nb_scalar_quadrangles, nb_vector_quadrangles, nb_tensor_quadrangles,
  nb_scalar_tetrahedra, nb_vector_tetrahedra, nb_tensor_tetrahedra,
  nb_scalar_hexahedra, nb_vector_hexahedra, nb_tensor_hexahedra,
  nb_scalar_prisms, nb_vector_prisms, nb_tensor_prisms,
  nb_scalar_pyramids, nb_vector_pyramids, nb_tensor_pyramids,
  nb_text2d, nb_text2d_chars, nb_text3d, nb_text3d_chars;

  std::vector<double> scalar_points, vector_points, tensor_points,
  scalar_lines, vector_lines, tensor_lines,
  scalar_triangles, vector_triangles, tensor_triangles,
  scalar_quadrangles, vector_quadrangles, tensor_quadrangles,
  scalar_tetrahedra, vector_tetrahedra, tensor_tetrahedra,
  scalar_hexahedra, vector_hexahedra, tensor_hexahedra,
  scalar_prisms, vector_prisms, tensor_prisms,
  scalar_pyramids, vector_pyramids, tensor_pyramids,
  time_step_values;

};


class xexport::xExportGmshAsciiSort : public xexport::xExportGmsh
{
public:
  
  virtual ~xexport::xExportGmshAsciiSort();
  void startView     (const string& comment);
  virtual void addDrawable  (xGmshDrawable::gType, int, double *, double *, double *, double *);	// 
  void endView ();
  void openFile (const string& fName);
  void closeFile ();
protected:
  std::ofstream *ofs;
  typedef std::multiset<xGmshDrawable,xGmshSortDrawable> xGmshDrawableSortContainer;
  xGmshDrawableSortContainer d_container;
};

#include "xExportGmsh_imp.h"

} // end of namespace

#endif
