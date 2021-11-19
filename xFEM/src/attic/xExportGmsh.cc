/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xExportGmsh.h"
#include "mPoint.h"
#include <fstream>
#include <sstream>
#include "workInProgress.h"
#ifdef WITH_GZIPOUTPUT 
#include <boost/iostreams/filter/gzip.hpp>
#endif
#ifdef PARALLEL
#include "mpi.h"
#endif


namespace xfem
{

  using Trellis_Util::mPoint;
  using namespace std;


  xexport::xExportGmsh::xexport::xExportGmsh() : xexport::xExport() { 
    filename_extension = ".pos";  
    // Do the same with mpi_com as argument
    // for now assert is commented out as it will block first para test that use there own mecanism to set names
    //assert(xtool::workInProgress());
#ifdef PARALLEL
    int rank;
    int mpisize;
    MPI_Comm_size(MPI_COMM_WORLD,&mpisize);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    std::stringstream filename_ex;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx" << std::endl;
    filename_ex <<"_PROC_"<< rank+1 << "_" << mpisize << ".pos";
    filename_extension = filename_ex.str();
#endif    
  }

void xexport::xExportGmsh::exportPoint ( const Trellis_Util::mPoint & P1, 
				const double & val1) 
{  
  double x[1] = {P1(0)};
  double y[1] = {P1(1)};
  double z[1] = {P1(2)};
  double v[1] = {val1};
  addDrawable ( xGmshDrawable :: SP , 1 , x,y,z,v ); 
}

void xexport::xExportGmsh::exportPoint ( const Trellis_Util::mPoint & P1, 
				const xtensor::xVector & val1) 
{  
  double x[1] = {P1(0)};
  double y[1] = {P1(1)};
  double z[1] = {P1(2)};
  double v[3] = {val1(0), val1(1), val1(2)};
  addDrawable ( xGmshDrawable :: VP , 1 , x,y,z,v ); 
}

void xexport::xExportGmsh::exportPoint ( const Trellis_Util::mPoint & P1, 
				const xtensor::xTensor2 & val1) 
{  
  double x[1] = {P1(0)};
  double y[1] = {P1(1)};
  double z[1] = {P1(2)};
  double v[9] = {val1(0,0),val1(0,1),val1(0,2),
		 val1(1,0),val1(1,1),val1(1,2),
		 val1(2,0),val1(2,1),val1(2,2)};
  addDrawable ( xGmshDrawable :: TP , 1 , x,y,z,v ); 
}


void xexport::xExportGmsh::exportLine ( const Trellis_Util::mPoint & P1, 
			       const Trellis_Util::mPoint & P2, 
			       const double & val1, 
			       const double & val2)
{  
  double x[2] = {P1(0),P2(0)};
  double y[2] = {P1(1),P2(1)};
  double z[2] = {P1(2),P2(2)};
  double v[2] = {val1,val2};
  addDrawable ( xGmshDrawable :: SL , 1 , x,y,z,v ); 
}

void xexport::xExportGmsh::exportLine ( const Trellis_Util::mPoint & P1, 
			      const Trellis_Util::mPoint & P2, 
			      const xtensor::xVector & val1, 
			      const xtensor::xVector & val2 )
{       
  double x[2] = {P1(0),P2(0)};
  double y[2] = {P1(1),P2(1)};
  double z[2] = {P1(2),P2(2)};
  double v[6] = {val1(0),val1(1),val1(2),val2(0),val2(1),val2(2)};
  addDrawable ( xGmshDrawable :: VL , 1 , x,y,z,v ); 
}

void xexport::xExportGmsh::exportLine (const Trellis_Util::mPoint & P1, 
			      const Trellis_Util::mPoint & P2, 
			      const xtensor::xTensor2 & val1, 
			      const xtensor::xTensor2 & val2)
{  
  double x[2] = {P1(0),P2(0)};
  double y[2] = {P1(1),P2(1)};
  double z[2] = {P1(2),P2(2)};
  double v[18] = {val1(0,0),val1(0,1),val1(0,2),
		  val1(1,0),val1(1,1),val1(1,2),
		  val1(2,0),val1(2,1),val1(2,2),
		  val2(0,0),val2(0,1),val2(0,2),
		  val2(1,0),val2(1,1),val2(1,2),
		  val2(2,0),val2(2,1),val2(2,2)};
  addDrawable ( xGmshDrawable :: TL , 1 , x,y,z,v ); 
}



void xexport::xExportGmsh::exportTriangle ( const Trellis_Util::mPoint & P1, 
				   const Trellis_Util::mPoint & P2, 
				   const Trellis_Util::mPoint & P3, 
				   const double & val1, 
				   const double & val2, 
				   const double & val3 )
{  
  double x[3] = {P1(0),P2(0),P3(0)};
  double y[3] = {P1(1),P2(1),P3(1)};
  double z[3] = {P1(2),P2(2),P3(2)};
  double v[3] = {val1,val2,val3};
  addDrawable ( xGmshDrawable :: ST , 1 , x,y,z,v ); 
}



void xexport::xExportGmsh::exportTriangle ( const Trellis_Util::mPoint & P1, 
				   const Trellis_Util::mPoint & P2, 
				   const Trellis_Util::mPoint & P3, 
				   const xtensor::xVector & val1, 
				   const xtensor::xVector & val2, 
				   const xtensor::xVector & val3 )
{  
  double x[3] = {P1(0),P2(0),P3(0)};
  double y[3] = {P1(1),P2(1),P3(1)};
  double z[3] = {P1(2),P2(2),P3(2)};
  double v[9] = {val1(0),val1(1),val1(2),
		 val2(0),val2(1),val2(2),
		 val3(0),val3(1),val3(2)};
  addDrawable ( xGmshDrawable :: VT , 1 , x,y,z,v ); 
}

void xexport::xExportGmsh::exportTriangle (const Trellis_Util::mPoint & P1, 
				  const Trellis_Util::mPoint & P2, 
				  const Trellis_Util::mPoint & P3, 
				  const xtensor::xTensor2 & val1, 
				  const xtensor::xTensor2 & val2, 
				  const xtensor::xTensor2 & val3 )
{  
  double x[3] = {P1(0),P2(0),P3(0)};
  double y[3] = {P1(1),P2(1),P3(1)};
  double z[3] = {P1(2),P2(2),P3(2)};
  double v[27] = {val1(0,0),val1(0,1),val1(0,2),
		  val1(1,0),val1(1,1),val1(1,2),
		  val1(2,0),val1(2,1),val1(2,2),
		  val2(0,0),val2(0,1),val2(0,2),
		  val2(1,0),val2(1,1),val2(1,2),
		  val2(2,0),val2(2,1),val2(2,2),
		  val3(0,0),val3(0,1),val3(0,2),
		  val3(1,0),val3(1,1),val3(1,2),
		  val3(2,0),val3(2,1),val3(2,2)};
  addDrawable ( xGmshDrawable :: TT , 1 , x,y,z,v ); 
}






void xexport::xExportGmsh::exportQuad ( const Trellis_Util::mPoint & P1, 
			       const Trellis_Util::mPoint & P2, 
			       const Trellis_Util::mPoint & P3, 
			       const Trellis_Util::mPoint & P4, 
			       const double & val1, 
			       const double & val2, 
			       const double & val3,
			       const double & val4 )
{  
  double x[4] = {P1(0),P2(0),P3(0), P4(0)};
  double y[4] = {P1(1),P2(1),P3(1), P4(1)};
  double z[4] = {P1(2),P2(2),P3(2), P4(2)};
  double v[4] = {val1,val2,val3,val4};
  addDrawable ( xGmshDrawable :: SQ , 1 , x,y,z,v ); 
}



void xexport::xExportGmsh::exportQuad ( const Trellis_Util::mPoint & P1, 
			       const Trellis_Util::mPoint & P2, 
			       const Trellis_Util::mPoint & P3,
			       const Trellis_Util::mPoint & P4,  
			       const xtensor::xVector & val1, 
			       const xtensor::xVector & val2, 
			       const xtensor::xVector & val3,
			       const xtensor::xVector & val4 )
{  
  double x[4] = {P1(0),P2(0),P3(0), P4(0)};
  double y[4] = {P1(1),P2(1),P3(1), P4(1)};
  double z[4] = {P1(2),P2(2),P3(2), P4(2)};
  double v[12] = {val1(0),val1(1),val1(2),
		  val2(0),val2(1),val2(2),
		  val3(0),val3(1),val3(2),
		  val4(0),val4(1),val4(2)};
  addDrawable ( xGmshDrawable :: VQ , 1 , x,y,z,v ); 
}

void xexport::xExportGmsh::exportQuad (const Trellis_Util::mPoint & P1, 
			      const Trellis_Util::mPoint & P2, 
			      const Trellis_Util::mPoint & P3, 
			      const Trellis_Util::mPoint & P4, 
			      const xtensor::xTensor2 & val1, 
			      const xtensor::xTensor2 & val2, 
			      const xtensor::xTensor2 & val3,
			      const xtensor::xTensor2 & val4 )
{  
  double x[4] = {P1(0),P2(0),P3(0), P4(0)};
  double y[4] = {P1(1),P2(1),P3(1), P4(1)};
  double z[4] = {P1(2),P2(2),P3(2), P4(2)};
  double v[36] = {val1(0,0),val1(0,1),val1(0,2),
		  val1(1,0),val1(1,1),val1(1,2),
		  val1(2,0),val1(2,1),val1(2,2),
		  val2(0,0),val2(0,1),val2(0,2),
		  val2(1,0),val2(1,1),val2(1,2),
		  val2(2,0),val2(2,1),val2(2,2),
		  val3(0,0),val3(0,1),val3(0,2),
		  val3(1,0),val3(1,1),val3(1,2),
		  val3(2,0),val3(2,1),val3(2,2),
		  val4(0,0),val4(0,1),val4(0,2),
		  val4(1,0),val4(1,1),val4(1,2),
		  val4(2,0),val4(2,1),val4(2,2)};
  addDrawable ( xGmshDrawable :: TQ , 1 , x,y,z,v ); 
}





void xexport::xExportGmsh::exportTetra (const Trellis_Util::mPoint & P1, 
			       const Trellis_Util::mPoint & P2, 
			       const Trellis_Util::mPoint & P3, 
			       const Trellis_Util::mPoint & P4, 
			       const double & val1, 
			       const double & val2, 
			       const double & val3, 
			       const double & val4 )
{  
  double x[4] = {P1(0),P2(0),P3(0), P4(0)};
  double y[4] = {P1(1),P2(1),P3(1), P4(1)};
  double z[4] = {P1(2),P2(2),P3(2), P4(2)};
  double v[4] = {val1,val2,val3,val4};
  addDrawable ( xGmshDrawable :: SS , 1 , x,y,z,v ); 
}


void xexport::xExportGmsh::exportTetra ( const Trellis_Util::mPoint & P1, 
				const Trellis_Util::mPoint & P2, 
				const Trellis_Util::mPoint & P3, 
				const Trellis_Util::mPoint & P4, 
				const xtensor::xVector & val1, 
				const xtensor::xVector & val2, 
				const xtensor::xVector & val3, 
				const xtensor::xVector & val4 )
{  
  double x[4] = {P1(0),P2(0),P3(0),P4(0)};
  double y[4] = {P1(1),P2(1),P3(1),P4(1)};
  double z[4] = {P1(2),P2(2),P3(2),P4(2)};
  double v[12] = {val1(0),val1(1),val1(2),
		  val2(0),val2(1),val2(2),
		  val3(0),val3(1),val3(2),
		  val4(0),val4(1),val4(2)};
  addDrawable ( xGmshDrawable :: VS , 1 , x,y,z,v ); 
}




void xexport::xExportGmsh::exportTetra (const Trellis_Util::mPoint & P1, 
			       const Trellis_Util::mPoint & P2, 
			       const Trellis_Util::mPoint & P3, 
			       const Trellis_Util::mPoint & P4, 
			       const xtensor::xTensor2 & val1, 
			       const xtensor::xTensor2 & val2, 
			       const xtensor::xTensor2 & val3, 
			       const xtensor::xTensor2 & val4 )
{  
  double x[4] = {P1(0),P2(0),P3(0),P4(0)};
  double y[4] = {P1(1),P2(1),P3(1),P4(1)};
  double z[4] = {P1(2),P2(2),P3(2),P4(2)};
  double v[36] = {val1(0,0),val1(0,1),val1(0,2),
		  val1(1,0),val1(1,1),val1(1,2),
		  val1(2,0),val1(2,1),val1(2,2),
		  val2(0,0),val2(0,1),val2(0,2),
		  val2(1,0),val2(1,1),val2(1,2),
		  val2(2,0),val2(2,1),val2(2,2),
		  val3(0,0),val3(0,1),val3(0,2),
		  val3(1,0),val3(1,1),val3(1,2),
		  val3(2,0),val3(2,1),val3(2,2),
		  val4(0,0),val4(0,1),val4(0,2),
		  val4(1,0),val4(1,1),val4(1,2),
		  val4(2,0),val4(2,1),val4(2,2)};
  addDrawable ( xGmshDrawable :: TS , 1 , x,y,z,v ); 
}


void xexport::xExportGmsh::exportHex (const Trellis_Util::mPoint & P1, 
			       const Trellis_Util::mPoint & P2, 
			       const Trellis_Util::mPoint & P3, 
			       const Trellis_Util::mPoint & P4, 
			       const Trellis_Util::mPoint & P5, 
			       const Trellis_Util::mPoint & P6, 
			       const Trellis_Util::mPoint & P7, 
			       const Trellis_Util::mPoint & P8, 
			       const double & val1, 
			       const double & val2, 
			       const double & val3, 
			       const double & val4,
			       const double & val5, 
			       const double & val6, 
			       const double & val7, 
			       const double & val8 )
{  
  double x[8] = {P1(0),P2(0),P3(0), P4(0), P5(0),P6(0),P7(0), P8(0)};
  double y[8] = {P1(1),P2(1),P3(1), P4(1), P5(1),P6(1),P7(1), P8(1)};
  double z[8] = {P1(2),P2(2),P3(2), P4(2), P5(2),P6(2),P7(2), P8(2)};
  double v[8] = {val1,val2,val3,val4,val5,val6,val7,val8};
  addDrawable ( xGmshDrawable :: SH , 1 , x,y,z,v ); 
}


void xexport::xExportGmsh::exportHex ( const Trellis_Util::mPoint & P1, 
				const Trellis_Util::mPoint & P2, 
				const Trellis_Util::mPoint & P3, 
				const Trellis_Util::mPoint & P4, 
				const Trellis_Util::mPoint & P5, 
			       const Trellis_Util::mPoint & P6, 
			       const Trellis_Util::mPoint & P7, 
			       const Trellis_Util::mPoint & P8, 
			       const xtensor::xVector & val1, 
			       const xtensor::xVector & val2, 
			       const xtensor::xVector & val3, 
			       const xtensor::xVector & val4,
			       const xtensor::xVector & val5, 
			       const xtensor::xVector & val6, 
			       const xtensor::xVector & val7, 
			       const xtensor::xVector & val8 )
{  
  double x[8] = {P1(0),P2(0),P3(0), P4(0), P5(0),P6(0),P7(0), P8(0)};
  double y[8] = {P1(1),P2(1),P3(1), P4(1), P5(1),P6(1),P7(1), P8(1)};
  double z[8] = {P1(2),P2(2),P3(2), P4(2), P5(2),P6(2),P7(2), P8(2)};
  double v[24] = {val1(0),val1(1),val1(2),
		  val2(0),val2(1),val2(2),
		  val3(0),val3(1),val3(2),
		  val4(0),val4(1),val4(2),
		  val5(0),val5(1),val5(2),
		  val6(0),val6(1),val6(2),
		  val7(0),val7(1),val7(2),
		  val8(0),val8(1),val8(2)};
  addDrawable ( xGmshDrawable :: VH , 1 , x,y,z,v ); 
}




void xexport::xExportGmsh::exportHex (const Trellis_Util::mPoint & P1, 
			       const Trellis_Util::mPoint & P2, 
			       const Trellis_Util::mPoint & P3, 
			       const Trellis_Util::mPoint & P4, 
				const Trellis_Util::mPoint & P5, 
			       const Trellis_Util::mPoint & P6, 
			       const Trellis_Util::mPoint & P7, 
			       const Trellis_Util::mPoint & P8, 
			       const xtensor::xTensor2 & val1, 
			       const xtensor::xTensor2 & val2, 
			       const xtensor::xTensor2 & val3, 
			       const xtensor::xTensor2 & val4,
			       const xtensor::xTensor2 & val5, 
			       const xtensor::xTensor2 & val6, 
			       const xtensor::xTensor2 & val7, 
			       const xtensor::xTensor2 & val8 )
{  
  double x[8] = {P1(0),P2(0),P3(0), P4(0), P5(0),P6(0),P7(0), P8(0)};
  double y[8] = {P1(1),P2(1),P3(1), P4(1), P5(1),P6(1),P7(1), P8(1)};
  double z[8] = {P1(2),P2(2),P3(2), P4(2), P5(2),P6(2),P7(2), P8(2)};
  double v[72] = {val1(0,0),val1(0,1),val1(0,2),
		  val1(1,0),val1(1,1),val1(1,2),
		  val1(2,0),val1(2,1),val1(2,2),
		  val2(0,0),val2(0,1),val2(0,2),
		  val2(1,0),val2(1,1),val2(1,2),
		  val2(2,0),val2(2,1),val2(2,2),
		  val3(0,0),val3(0,1),val3(0,2),
		  val3(1,0),val3(1,1),val3(1,2),
		  val3(2,0),val3(2,1),val3(2,2),
		  val4(0,0),val4(0,1),val4(0,2),
		  val4(1,0),val4(1,1),val4(1,2),
		  val4(2,0),val4(2,1),val4(2,2),
		  val5(0,0),val5(0,1),val5(0,2),
		  val5(1,0),val5(1,1),val5(1,2),
		  val5(2,0),val5(2,1),val5(2,2),
		  val6(0,0),val6(0,1),val6(0,2),
		  val6(1,0),val6(1,1),val6(1,2),
		  val6(2,0),val6(2,1),val6(2,2),
		  val7(0,0),val7(0,1),val7(0,2),
		  val7(1,0),val7(1,1),val7(1,2),
		  val7(2,0),val7(2,1),val7(2,2),
		  val8(0,0),val8(0,1),val8(0,2),
		  val8(1,0),val8(1,1),val8(1,2),
		  val8(2,0),val8(2,1),val8(2,2)};
  addDrawable ( xGmshDrawable :: TH , 1 , x,y,z,v ); 
}


void xexport::xExportGmshAscii::openFile (const string& fName)
{
  //  std::cout << getFileNameExtension() << std::endl;
  string loc = fName + getFileNameExtension();
  ofs = new std::ofstream(loc.c_str());
  process_started = true;
}

void xexport::xExportGmshAscii::startView (const string& commentar)
{
  (*ofs) << "View \"" << commentar << "\" {" << std::endl; 
}

void xexport::xExportGmshAscii:: addDrawable  (xGmshDrawable::gType t, int n, double *x, double *y, double *z, double *v)
{
  xGmshDrawable d (t,n,x,y,z,v);
  (*ofs) << d;
}

void xexport::xExportGmshAscii::endView ()
{
  (*ofs) << "};" << std::endl;
}

void xexport::xExportGmshAscii::closeFile ()
{
  if(!process_started)return; 
  ofs->close();
  delete ofs;
  process_started = false;
}

xexport::xExportGmshAscii::~xexport::xExportGmshAscii()
{
  closeFile();
}

#ifdef WITH_GZIPOUTPUT 
void xexport::xExportGmshAsciigz::openFile (const string& fName)
{
    // add to the chain gzip compresor filter 
    ofsf.push(boost::iostreams::gzip_compressor());

    // open a standard stream 
    string loc = fName + getFileNameExtension() + ".gz";
    ofs = new std::ofstream(loc.c_str(),ios::out|ios::binary);

    if (!ofs) throw;

    // add the standard stream to the chain
    ofsf.push(*ofs);

    process_started = true;
}

void xexport::xExportGmshAsciigz::startView (const string& commentar)
{
  ofsf << "View \"" << commentar << "\" {" << std::endl; 
}

void xexport::xExportGmshAsciigz:: addDrawable  (xGmshDrawable::gType t, int n, double *x, double *y, double *z, double *v)
{
  xGmshDrawable d (t,n,x,y,z,v);
  ofsf << d;
}

void xexport::xExportGmshAsciigz::endView ()
{
  ofsf << "};" << std::endl;
}

void xexport::xExportGmshAsciigz::closeFile ()
{
  if(!process_started)return; 
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


xexport::xExportGmshAsciigz::~xexport::xExportGmshAsciigz()
{
  closeFile();
}
#endif


void xexport::xExportGmshAsciiSort::openFile (const string& fName)
{
  string loc = fName + getFileNameExtension();
  ofs = new std::ofstream(loc.c_str());
  process_started = true;
  d_container.clear();
}

void xexport::xExportGmshAsciiSort::startView (const string& commentar)
{
  (*ofs) << "View \"" << commentar << "\" {" << std::endl; 
}

void xexport::xExportGmshAsciiSort:: addDrawable  (xGmshDrawable::gType t, int n, double *x, double *y, double *z, double *v)
{
  xGmshDrawable d (t,n,x,y,z,v);
  d_container.insert(d);
}

void xexport::xExportGmshAsciiSort::endView ()
{
  xGmshDrawableSortContainer::iterator it = d_container.begin();
  xGmshDrawableSortContainer::iterator itend = d_container.end();
  for (;it!=itend;++it)
  {
      (*ofs) << (*it);
  }
  (*ofs) << "};" << std::endl;
  d_container.clear();
}

void xexport::xExportGmshAsciiSort::closeFile ()
{
  if(!process_started)return; 
  ofs->close();
  delete ofs;
  process_started = false;
}

xexport::xExportGmshAsciiSort::~xexport::xExportGmshAsciiSort()
{
  closeFile();
}







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




void xexport::xExportGmshBinary::clear()
{
  view_name = "";
  nb_time_steps = 1;
  time_step_values.clear();
  time_step_values.push_back(1.0);
  
  
  nb_scalar_points= 0; nb_vector_points= 0; nb_tensor_points= 0;
  nb_scalar_lines= 0; nb_vector_lines= 0; nb_tensor_lines= 0;
  nb_scalar_triangles= 0; nb_vector_triangles= 0; nb_tensor_triangles= 0;
  nb_scalar_quadrangles= 0; nb_vector_quadrangles= 0; nb_tensor_quadrangles= 0;
  nb_scalar_tetrahedra= 0; nb_vector_tetrahedra= 0; nb_tensor_tetrahedra= 0;
  nb_scalar_hexahedra= 0; nb_vector_hexahedra= 0; nb_tensor_hexahedra= 0;
  nb_scalar_prisms= 0; nb_vector_prisms= 0; nb_tensor_prisms= 0;
  nb_scalar_pyramids= 0; nb_vector_pyramids= 0; nb_tensor_pyramids= 0;
  nb_text2d= 0; nb_text2d_chars= 0; nb_text3d= 0; nb_text3d_chars= 0;
  
  scalar_points.clear(); vector_points.clear(); tensor_points.clear();
  scalar_lines.clear(); vector_lines.clear(); tensor_lines.clear();
  scalar_triangles.clear(); vector_triangles.clear(); tensor_triangles.clear();
  scalar_quadrangles.clear(); vector_quadrangles.clear(); tensor_quadrangles.clear();
  scalar_tetrahedra.clear(); vector_tetrahedra.clear(); tensor_tetrahedra.clear();
  scalar_hexahedra.clear(); vector_hexahedra.clear(); tensor_hexahedra.clear();
  scalar_prisms.clear(); vector_prisms.clear(); tensor_prisms.clear();
  scalar_pyramids.clear(); vector_pyramids.clear(); tensor_pyramids.clear();

}



xexport::xExportGmshBinary::xexport::xExportGmshBinary() : 
    xexport::xExportGmsh()
{
  clear();
}

void xexport::xExportGmshBinary::openFile (const string& fName)
{
  clear();
  string loc = fName + getFileNameExtension(); 
  file = fopen(loc.c_str(), "w"); 
  process_started = true;
}

void xexport::xExportGmshBinary::startView (const string& commentar)
{
  clear();
  view_name = commentar;
  fprintf(file, "$PostFormat\n");
  fprintf(file, "%g %d %d\n", 1.2, 1, (int)sizeof(double));
  fprintf(file, "$EndPostFormat\n");
  fprintf(file, "$View\n");
}

void xexport::xExportGmshBinary:: addDrawable  (xGmshDrawable::gType t, int n, double *x, double *y, double *z, double *v)
{
 std::vector<double>* vector_ptr;
 int nbNod, sizeInfo;
 int nbSamples = n;
  switch(t)
    {
    case xGmshDrawable::NONE:printf("error drawable\n");throw 1;
    case xGmshDrawable::SL  : nb_scalar_points++;      vector_ptr = &scalar_points; nbNod = 2;sizeInfo=1;break;     
    case xGmshDrawable::ST  : nb_scalar_triangles++;   vector_ptr = &scalar_triangles; nbNod = 3;sizeInfo=1;break;
    case xGmshDrawable::SQ  : nb_scalar_quadrangles++; vector_ptr = &scalar_quadrangles; nbNod = 4;sizeInfo=1;break;
    case xGmshDrawable::SS  : nb_scalar_tetrahedra++;  vector_ptr = &scalar_tetrahedra; nbNod = 4;sizeInfo=1;break;
    case xGmshDrawable::SH  : nb_scalar_hexahedra++;  vector_ptr = &scalar_hexahedra; nbNod = 8;sizeInfo=1;break;
    case xGmshDrawable::VL  : nb_vector_lines++;       vector_ptr = &vector_lines; nbNod = 2;sizeInfo=3;break;
    case xGmshDrawable::VT  : nb_vector_triangles++;   vector_ptr = &vector_triangles; nbNod = 3;sizeInfo=3;break;
    case xGmshDrawable::VQ  : nb_vector_quadrangles++; vector_ptr = &vector_quadrangles; nbNod = 4;sizeInfo=3;break;
    case xGmshDrawable::VS  : nb_vector_tetrahedra++;  vector_ptr = &vector_tetrahedra; nbNod = 4;sizeInfo=3;break;
    case xGmshDrawable::VH  : nb_vector_hexahedra++;  vector_ptr = &vector_hexahedra; nbNod = 8;sizeInfo=3;break;
    case xGmshDrawable::TL  : nb_tensor_lines++;       vector_ptr = &tensor_lines; nbNod = 2;sizeInfo=9;break;
    case xGmshDrawable::TT  : nb_tensor_triangles++;   vector_ptr = &tensor_triangles; nbNod = 3;sizeInfo=9;break;
    case xGmshDrawable::TQ  : nb_tensor_quadrangles++; vector_ptr = &tensor_quadrangles; nbNod = 4;sizeInfo=9;break;
    case xGmshDrawable::TS  : nb_tensor_tetrahedra++;  vector_ptr = &tensor_tetrahedra; nbNod = 4;sizeInfo=9;break;
    case xGmshDrawable::TH  : nb_tensor_hexahedra++;  vector_ptr = &tensor_hexahedra; nbNod = 8;sizeInfo=9;break;
    default: throw;
    }
    int nbval = nbNod*nbSamples*sizeInfo;
    vector_ptr->insert(vector_ptr->end(), x, x+nbNod);
    vector_ptr->insert(vector_ptr->end(), y, y+nbNod);
    vector_ptr->insert(vector_ptr->end(), z, z+nbNod);
    vector_ptr->insert(vector_ptr->end(), v, v+nbval);
}

void xexport::xExportGmshBinary::endView ()
{
  int one = 1;
  fprintf(file, "%s %d "
	  "%d %d %d "
	  "%d %d %d "
	  "%d %d %d "
	  "%d %d %d "
	  "%d %d %d "
	  "%d %d %d "
	  "%d %d %d "
	  "%d %d %d "
	  "%d %d %d %d\n", 
	  view_name.c_str(), nb_time_steps,
	  nb_scalar_points, nb_vector_points, nb_tensor_points,
	  nb_scalar_lines, nb_vector_lines, nb_tensor_lines,
	  nb_scalar_triangles, nb_vector_triangles, nb_tensor_triangles,
	  nb_scalar_quadrangles, nb_vector_quadrangles, nb_tensor_quadrangles,
	  nb_scalar_tetrahedra, nb_vector_tetrahedra, nb_tensor_tetrahedra,
	  nb_scalar_hexahedra, nb_vector_hexahedra, nb_tensor_hexahedra,
	  nb_scalar_prisms, nb_vector_prisms, nb_tensor_prisms,
	  nb_scalar_pyramids, nb_vector_pyramids, nb_tensor_pyramids,
	  nb_text2d, nb_text2d_chars, nb_text3d, nb_text3d_chars);
//ascii long
  fwrite(&one, sizeof(int), 1, file);
  fwrite(&(*time_step_values.begin()), sizeof(double), nb_time_steps, file);
  fwrite(&(*scalar_points.begin()), sizeof(double), scalar_points.size(),file);
  fwrite(&(*vector_points.begin()), sizeof(double), vector_points.size(),file);
  fwrite(&(*tensor_points.begin()), sizeof(double), tensor_points.size(),file);
  fwrite(&(*scalar_triangles.begin()), sizeof(double), scalar_triangles.size(),file);
  fwrite(&(*vector_triangles.begin()), sizeof(double), vector_triangles.size(),file);
  fwrite(&(*tensor_triangles.begin()), sizeof(double), tensor_triangles.size(),file);
  fwrite(&(*scalar_quadrangles.begin()), sizeof(double), scalar_quadrangles.size(),file);
  fwrite(&(*vector_quadrangles.begin()), sizeof(double), vector_quadrangles.size(),file);
  fwrite(&(*tensor_quadrangles.begin()), sizeof(double), tensor_quadrangles.size(),file);
  fwrite(&(*scalar_tetrahedra.begin()), sizeof(double), scalar_tetrahedra.size(),file);
  fwrite(&(*vector_tetrahedra.begin()), sizeof(double), vector_tetrahedra.size(),file);
  fwrite(&(*tensor_tetrahedra.begin()), sizeof(double), tensor_tetrahedra.size(),file);
  fwrite(&(*scalar_hexahedra.begin()), sizeof(double), scalar_hexahedra.size(),file);
  fwrite(&(*vector_hexahedra.begin()), sizeof(double), vector_hexahedra.size(),file);
  fwrite(&(*tensor_hexahedra.begin()), sizeof(double), tensor_hexahedra.size(),file);
  fwrite(&(*scalar_prisms.begin()), sizeof(double), scalar_prisms.size(),file);
  fwrite(&(*vector_prisms.begin()), sizeof(double), vector_prisms.size(),file);
  fwrite(&(*tensor_prisms.begin()), sizeof(double), tensor_prisms.size(),file);
  fwrite(&(*scalar_pyramids.begin()), sizeof(double), scalar_pyramids.size(),file);
  fwrite(&(*vector_pyramids.begin()), sizeof(double), vector_pyramids.size(),file);
  fwrite(&(*tensor_pyramids.begin()), sizeof(double), tensor_pyramids.size(),file);
  fprintf(file, "\n$EndView\n");
}

void xexport::xExportGmshBinary::closeFile ()
{
  if(!process_started)return; 
  fclose(file);
  process_started = false;
}

xexport::xExportGmshBinary::~xexport::xExportGmshBinary()
{
 closeFile();
}






xGmshDrawable::xGmshDrawable()
  : nbNod(0),sizeInfo(0), X(0),Y(0),Z(0),V(0),nbSamples(0),type(NONE)
{
}
xGmshDrawable::xGmshDrawable( const  xGmshDrawable& in)
  : nbNod(in.nbNod),sizeInfo(in.sizeInfo), nbSamples(in.nbSamples),type(in.type)
{
  const int s_v =nbNod*nbSamples*sizeInfo;
  X = new double[nbNod];
  Y = new double[nbNod];
  Z = new double[nbNod];
  V = new double[s_v];
  memcpy(X,in.X,nbNod*sizeof(double));
  memcpy(Y,in.Y,nbNod*sizeof(double));
  memcpy(Z,in.Z,nbNod*sizeof(double));
  memcpy(V,in.V,s_v*sizeof(double));
}

xGmshDrawable::xGmshDrawable(gType t, int n, double *x, double *y, double *z, double *v)
  : nbSamples(n), type (t)
{
  switch(type)
    {
    case NONE:printf("error drawable\n");throw 1;
    case SP  :nbNod = 1;sizeInfo=1;break;
    case SL  :nbNod = 2;sizeInfo=1;break;
    case ST  :nbNod = 3;sizeInfo=1;break;
    case SS  :nbNod = 4;sizeInfo=1;break;
    case SQ  :nbNod = 4;sizeInfo=1;break;
    case SH  :nbNod = 8;sizeInfo=1;break;
    case VP  :nbNod = 1;sizeInfo=3;break;
    case VL  :nbNod = 2;sizeInfo=3;break;
    case VT  :nbNod = 3;sizeInfo=3;break;
    case VS  :nbNod = 4;sizeInfo=3;break;
    case VQ  :nbNod = 4;sizeInfo=3;break;
    case VH  :nbNod = 8;sizeInfo=3;break;
    case TP  :nbNod = 1;sizeInfo=9;break;  
    case TL  :nbNod = 2;sizeInfo=9;break;
    case TT  :nbNod = 3;sizeInfo=9;break;
    case TS  :nbNod = 4;sizeInfo=9;break;
    case TQ  :nbNod = 4;sizeInfo=9;break;
    case TH  :nbNod = 8;sizeInfo=9;break;
    default: assert(0); break;
    }
  X = new double[nbNod];
  Y = new double[nbNod];
  Z = new double[nbNod];
  V = new double[nbNod*nbSamples*sizeInfo];
  for(int i=0;i<nbNod;i++)
    {
      X[i] = x[i];
      Y[i] = y[i];
      Z[i] = z[i];
    }
  memcpy(V,v,nbNod*nbSamples*sizeInfo*sizeof(double));
}

xGmshDrawable::~xGmshDrawable()
{
  if(X)delete [] X;
  if(Y)delete [] Y;
  if(Z)delete [] Z;
  if(V)delete [] V;
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
  if (type<other.type) return true;
  if (type>other.type) return false;
  
  // local
  int i,n=nbNod;

  // if of same type second level of sorting on nodes
  for(i=0;i<n;++i)
  {
      if (X[i]<other.X[i]) return true;
      if (X[i]>other.X[i]) return false;
      if (Y[i]<other.Y[i]) return true;
      if (Y[i]>other.Y[i]) return false;
      if (Z[i]<other.Z[i]) return true;
      if (Z[i]>other.Z[i]) return false;
  }

  // if of same type with same nodes third level of sorting on values
  n*=nbSamples*sizeInfo;
  for(i=0;i<n;++i)
  {
      if (V[i]<other.V[i]) return true;
      if (V[i]>other.V[i]) return false;
  } 

  // exactely the same
  return false;
}


bool xGmshSortDrawable::operator()(const xGmshDrawable& d1, const xGmshDrawable& d2) const 
{
    return d1.lessthan(d2);
}


} // end of namespace



