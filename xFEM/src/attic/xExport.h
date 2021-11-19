/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _x_EXPORT_
#define _x_EXPORT_

#include<string>
#include<vector>

#include "xEval.h"
#include "xIntegrationRule.h"
#include "xEntityToEntity.h"
#include "xCommandOnGeomElem.h"
#include "mEntity.h"
#include "mVertex.h"
#include "mPoint.h"
#include "xMappingBuilderHolder.h"




namespace xfem
{
  class xtensor::xVector;
  class xtensor::xTensor2;
using std::string;

class xexport::xExport 
{
public:
  xexport::xExport() : process_started(false), nbsplit(1), nbsplitdefault(1), filterSplit(xAcceptAll())  {}
  bool processStarted() const {return process_started; }
  void setNbSplit(int n) {nbsplit=n;}
  void setNbSplitDefault(int n) {nbsplitdefault=n;}

  int  getNbSplit() const {return nbsplit;}
  int getNbSplit( AOMD::mEntity *e) const {return filterSplit(e)?nbsplit:nbsplitdefault;}
  string getFileNameExtension() const {return filename_extension; }
  void setFilterSplit(const xEntityFilter &fin){filterSplit=fin;}
  virtual ~xexport::xExport() {}
  virtual void openFile  (const string& fName) = 0;
  virtual void startView (const string& comment) = 0;

  virtual void exportPoint (const Trellis_Util::mPoint & P1, const double   & val1) = 0;
  virtual void exportPoint (const Trellis_Util::mPoint & P1, const xtensor::xVector  & val1) = 0;
  virtual void exportPoint (const Trellis_Util::mPoint & P1, const xtensor::xTensor2 & val1) = 0;
  
  virtual void exportLine (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const double & val1, const double & val2) = 0;
  virtual void exportLine (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const xtensor::xVector & val1, const xtensor::xVector & val2) = 0;
  virtual void exportLine (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2) = 0;

  virtual void exportTriangle (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, 
			       const double & val1, const double & val2, const double & val3 ) = 0; 
  virtual void exportTriangle (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, 
			       const xtensor::xVector & val1, const xtensor::xVector & val2, const xtensor::xVector & val3 ) = 0;
  virtual void exportTriangle (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, 
			       const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2, const xtensor::xTensor2 & val3 ) = 0;

  
  virtual void exportQuad (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
			    const double & val1, const double & val2, const double & val3, const double & val4 ) = 0;
  virtual void exportQuad (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
			    const xtensor::xVector & val1, const xtensor::xVector & val2, const xtensor::xVector & val3, const xtensor::xVector & val4) = 0;
  virtual void exportQuad (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
			    const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2, const xtensor::xTensor2 & val3, const xtensor::xTensor2 & val4 ) = 0;

  virtual void exportTetra (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
			    const double & val1, const double & val2, const double & val3, const double & val4 ) = 0;
  virtual void exportTetra (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
			    const xtensor::xVector & val1, const xtensor::xVector & val2, const xtensor::xVector & val3, const xtensor::xVector & val4) = 0;
  virtual void exportTetra (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4, 
			    const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2, const xtensor::xTensor2 & val3, const xtensor::xTensor2 & val4 ) = 0;

  virtual void exportHex (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4,
			  const Trellis_Util::mPoint & P5, const Trellis_Util::mPoint & P6, const Trellis_Util::mPoint & P7, const Trellis_Util::mPoint & P8, 
			    const double & val1, const double & val2, const double & val3, const double & val4, const double & val5, const double & val6, const double & val7, const double & val8 ) = 0;
  virtual void exportHex (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4,
			  const Trellis_Util::mPoint & P5, const Trellis_Util::mPoint & P6, const Trellis_Util::mPoint & P7, const Trellis_Util::mPoint & P8, 
			    const xtensor::xVector & val1, const xtensor::xVector & val2, const xtensor::xVector & val3, const xtensor::xVector & val4, const xtensor::xVector & val5, const xtensor::xVector & val6, const xtensor::xVector & val7, const xtensor::xVector & val8 ) = 0;
  virtual void exportHex (const Trellis_Util::mPoint & P1, const Trellis_Util::mPoint & P2, const Trellis_Util::mPoint & P3, const Trellis_Util::mPoint & P4,
			  const Trellis_Util::mPoint & P5, const Trellis_Util::mPoint & P6, const Trellis_Util::mPoint & P7, const Trellis_Util::mPoint & P8, 
			    const xtensor::xTensor2 & val1, const xtensor::xTensor2 & val2, const xtensor::xTensor2 & val3, const xtensor::xTensor2 & val4, const xtensor::xTensor2 & val5, const xtensor::xTensor2 & val6, const xtensor::xTensor2 & val7, const xtensor::xTensor2 & val8 ) = 0;



  //virtual void exportNode(ostream &fs, mVertex *f, int &k);
  //virtual void exportEdge(ostream &fs, mEdge *f, int &k);
  //virtual void exportFace(ostream &fs, mFace *f, int &k);
  //virtual void exportVolume(ostream &fs, mEdge *f, int &k);



  virtual void endView ()  = 0;
  virtual void closeFile ()  = 0;


protected:
  bool process_started;
  int nbsplit;
  int nbsplitdefault;
  string filename_extension;
  xEntityFilter filterSplit;

private:


};



class xSplitEdge
{
  int level;
  Trellis_Util::mPoint t[2];
  int thisEdge,nF;
public :
  xSplitEdge(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2);
  int nbEdge ();
  bool nextEdge(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2);
};

class xSplitTri
{
  int level;
  Trellis_Util::mPoint t[3];
  int thisTri,iF,nF;
public :
  xSplitTri(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2, const Trellis_Util::mPoint &p3);
  int nbTri ();
  bool nextTri(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2, Trellis_Util::mPoint &pt3);
};

class xSplitQuad
{
  int level;
  Trellis_Util::mPoint t[4];
  int thisQuad,iF,nF;
  std::vector<Trellis_Util::mPoint> v[4];

public :
  xSplitQuad(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2, const Trellis_Util::mPoint &p3, const Trellis_Util::mPoint &p4);
  int nbQuad ();
  bool nextQuad(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2, Trellis_Util::mPoint &pt3, Trellis_Util::mPoint &pt4);
};

class xSplitTet
{
  std::vector<Trellis_Util::mPoint> v[4];
  int thisTet;
public :
  xSplitTet(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2, 
	   const Trellis_Util::mPoint &p3,const Trellis_Util::mPoint &p4);
  int nbTet ();
  bool nextTet(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2, Trellis_Util::mPoint &pt3, Trellis_Util::mPoint &pt4);
};

class xSplitHex
{
  int level;
  Trellis_Util::mPoint t[8];
  int thisHex,iF,nF;
  std::vector<Trellis_Util::mPoint> v[8];
public :
  xSplitHex(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2, const Trellis_Util::mPoint &p3, const Trellis_Util::mPoint &p4, const Trellis_Util::mPoint &pt5, const Trellis_Util::mPoint &pt6, const Trellis_Util::mPoint &pt7, const Trellis_Util::mPoint &pt8);
  int nbHex ();
  bool nextHex(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2, Trellis_Util::mPoint &pt3, Trellis_Util::mPoint &pt4, Trellis_Util::mPoint &pt5, Trellis_Util::mPoint &pt6, Trellis_Util::mPoint &pt7, Trellis_Util::mPoint &pt8);
};





template <typename T>
class xexport::xPlotCommand :  public xCommandOnGeomElem
{

private:
  const xEval<T>& eval;
  T val1, val2, val3, val4, val5, val6, val7, val8;
  xexport::xExport& pexport;
  int recurSplit;
  const bool simplex;

public:
  xexport::xPlotCommand (const xEval<T>& eval_, xexport::xExport& e, const bool simplex_=true)
    : eval(eval_), pexport(e), recurSplit( e.getNbSplit()), simplex(simplex_)  {}
  
  void execute(xGeomElem* geom_integ) 
  {
    AOMD::mEntity* e = geom_integ->getEntity();
    Trellis_Util::mPoint puvw1,puvw2,puvw3,puvw4,puvw5,puvw6,puvw7,puvw8,pxyz1,pxyz2, pxyz3, pxyz4,pxyz5,pxyz6, pxyz7, pxyz8; 
    std::vector<AOMD::mEntity *> Vertex ;
    std::vector<xSplitEdge> EdgesInRefSpace ;
    std::vector<xSplitTri> trianglesInRefSpace ;
    std::vector<xSplitQuad> quadsInRefSpace ;
    std::vector<xSplitTet> tetsInRefSpace ;
    std::vector<xSplitHex> hexsInRefSpace ;
    recurSplit= pexport.getNbSplit(geom_appro->getEntity());
    switch(e->getType()) 
      {
      case AOMD::mEntity::VERTEX :
	Vertex.push_back((AOMD::mEntity*) e);
	break;
      case AOMD::mEntity::EDGE :
	EdgesInRefSpace.push_back ( xSplitEdge ( recurSplit , 
						 Trellis_Util::mPoint (-1,0,0) ,
						 Trellis_Util::mPoint (1,0,0 ) ) );
	
	
	break;
      case AOMD::mEntity::TRI :
	trianglesInRefSpace.push_back ( xSplitTri ( recurSplit , 
						    Trellis_Util::mPoint (0,0,0) ,
						    Trellis_Util::mPoint (1,0,0) ,
						    Trellis_Util::mPoint (0,1,0) ) );
	break;
      case AOMD::mEntity::TET :
	tetsInRefSpace.push_back ( xSplitTet ( recurSplit , 
					       Trellis_Util::mPoint (0,0,0) ,
					       Trellis_Util::mPoint (1,0,0) ,
					       Trellis_Util::mPoint (0,1,0) ,
					       Trellis_Util::mPoint (0,0,1) ) );
	break;
      case AOMD::mEntity::QUAD :
	quadsInRefSpace.push_back ( xSplitQuad ( recurSplit , 
						 Trellis_Util::mPoint (-1.,-1.,0.) ,
						 Trellis_Util::mPoint (1.,-1.,0.) ,
						 Trellis_Util::mPoint (1.,1.,0.) ,
						 Trellis_Util::mPoint (-1.,1.,0.) ) );

// 	trianglesInRefSpace.push_back ( xSplitTri ( recurSplit , 
// 						    Trellis_Util::mPoint (-1,-1,0) ,
// 						    Trellis_Util::mPoint (-1,1,0) ,
// 						    Trellis_Util::mPoint (1,1,0) ) );
// 	trianglesInRefSpace.push_back ( xSplitTri ( recurSplit , 
// 						    Trellis_Util::mPoint (1,1,0) ,
// 						    Trellis_Util::mPoint (1,-1,0) ,
// 						    Trellis_Util::mPoint (-1,-1,0) ) );
	break;	  
      case AOMD::mEntity::HEX :
      if(simplex){
	tetsInRefSpace.push_back ( xSplitTet ( recurSplit , 
					       Trellis_Util::mPoint (-1,1,-1),
					       Trellis_Util::mPoint (-1,1,1),
					       Trellis_Util::mPoint ( 1,-1,1),
					       Trellis_Util::mPoint (1,1,1) ) );
	tetsInRefSpace.push_back ( xSplitTet ( recurSplit , 
					       Trellis_Util::mPoint (-1,1,-1), 
					       Trellis_Util::mPoint (1,1,1),
					       Trellis_Util::mPoint ( 1,-1,1),
					       Trellis_Util::mPoint ( 1,1,-1) ));
	tetsInRefSpace.push_back ( xSplitTet ( recurSplit , 
					       Trellis_Util::mPoint (-1,1,-1), 
					       Trellis_Util::mPoint (-1,-1,1),
					       Trellis_Util::mPoint ( 1,-1,1),
					       Trellis_Util::mPoint ( -1,1,1) ));
	tetsInRefSpace.push_back ( xSplitTet ( recurSplit , 
					       Trellis_Util::mPoint (-1,1,-1), 
					       Trellis_Util::mPoint (1,1,-1),
					       Trellis_Util::mPoint (1,-1,1),
					       Trellis_Util::mPoint ( 1,-1,-1) ));
	tetsInRefSpace.push_back ( xSplitTet ( recurSplit , 
					       Trellis_Util::mPoint (-1,1,-1), 
					       Trellis_Util::mPoint (1,-1,-1),
					       Trellis_Util::mPoint (1,-1,1),
					       Trellis_Util::mPoint (-1,-1,-1) ));
	tetsInRefSpace.push_back ( xSplitTet ( recurSplit , 
					       Trellis_Util::mPoint (-1,1,-1), 
					       Trellis_Util::mPoint (-1,-1,-1),
					       Trellis_Util::mPoint ( 1,-1,1),
					       Trellis_Util::mPoint (-1,-1,1) ));
	}else{
	hexsInRefSpace.push_back ( xSplitHex ( recurSplit , 
					       Trellis_Util::mPoint (-1,-1,-1), 
					       Trellis_Util::mPoint (1,-1,-1),
					       Trellis_Util::mPoint ( 1,1,-1),
					       Trellis_Util::mPoint (-1,1,-1), 
					       Trellis_Util::mPoint (-1,-1,1), 
					       Trellis_Util::mPoint (1,-1,1),
					       Trellis_Util::mPoint ( 1,1,1),
					       Trellis_Util::mPoint (-1,1,1) ));
	}
	break;
      default : ;
      }
    for (int i=0;(unsigned)i<Vertex.size();i++)
      {
	pxyz1 = geom_integ->getXYZ();
	if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	else geom_appro->setUVW(geom_integ->getUVW());
	eval(geom_appro, geom_integ, val1);
	pexport.exportPoint(pxyz1,val1);
      }
    for (int i=0;(unsigned)i<EdgesInRefSpace.size();i++)
      {
	for (int j=0;j<EdgesInRefSpace[i].nbEdge();j++)
	  {	      
	    EdgesInRefSpace[i].nextEdge ( puvw1 , puvw2 );
	    /// 1st value
	    geom_integ->setUVW (puvw1);
	    pxyz1 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());	    
	    eval(geom_appro, geom_integ, val1);

	    /// 2nd value
	    geom_integ->setUVW (puvw2);
	    pxyz2 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val2);
	    pexport.exportLine(pxyz1,pxyz2,val1,val2);
	  }
      }
    for (int i=0;(unsigned)i<trianglesInRefSpace.size();i++)
      {
	for (int j=0;j<trianglesInRefSpace[i].nbTri();j++)
	  {	      
	    trianglesInRefSpace[i].nextTri ( puvw1 , puvw2 , puvw3);
	    /// 1st value
	    geom_integ->setUVW (puvw1);
	    pxyz1 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());   
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val1);
	    /// 2nd value
	    geom_integ->setUVW (puvw2);
	    pxyz2 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val2);
	    /// 3rd value
	    geom_integ->setUVW (puvw3);
	    pxyz3 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val3);	      
	    pexport.exportTriangle(pxyz1,pxyz2,pxyz3,val1,val2,val3);
	  }
      }
    for (int i=0;(unsigned)i<quadsInRefSpace.size();i++)
      {
	for (int j=0;j<quadsInRefSpace[i].nbQuad();j++)
	  {	      
	    quadsInRefSpace[i].nextQuad ( puvw1 , puvw2 , puvw3, puvw4);
	    /// 1st value
	    geom_integ->setUVW (puvw1);
	    pxyz1 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());   
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val1);
	    /// 2nd value
	    geom_integ->setUVW (puvw2);
	    pxyz2 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val2);
	    /// 3rd value
	    geom_integ->setUVW (puvw3);
	    pxyz3 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val3);	      
	    /// 4th value
	    geom_integ->setUVW (puvw4);
	    pxyz4 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val4);	      
	    pexport.exportQuad(pxyz1,pxyz2,pxyz3,pxyz4,val1,val2,val3,val4);
	  }
      }
    for (int i=0;(unsigned)i<tetsInRefSpace.size();i++)
      {
	for (int j=0;j<tetsInRefSpace[i].nbTet();j++)
	  {	      
	    tetsInRefSpace[i].nextTet ( puvw1 , puvw2 , puvw3, puvw4);
	    /// 1st value
	    geom_integ->setUVW (puvw1);
	    pxyz1 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val1);
	    
	    /// 2nd value
	    geom_integ->setUVW (puvw2);
	    pxyz2 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val2);
	    
	    /// 3rd value
	    geom_integ->setUVW (puvw3);
	    pxyz3 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val3);	      
	    
	    /// 4th value
	    geom_integ->setUVW (puvw4);
	    pxyz4 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val4);	    
	    pexport.exportTetra(pxyz1,pxyz2,pxyz3,pxyz4,val1,val2,val3,val4);
	    
	  }
      }
    for (int i=0;(unsigned)i<hexsInRefSpace.size();i++)
      {
	for (int j=0;j<hexsInRefSpace[i].nbHex();j++)
	  {	      
	    hexsInRefSpace[i].nextHex ( puvw1 , puvw2 , puvw3, puvw4, puvw5 , puvw6 , puvw7, puvw8);
	    /// 1st value
	    geom_integ->setUVW (puvw1);
	    pxyz1 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val1);
	    
	    /// 2nd value
	    geom_integ->setUVW (puvw2);
	    pxyz2 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val2);
	    
	    /// 3rd value
	    geom_integ->setUVW (puvw3);
	    pxyz3 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val3);	      
	    
	    /// 4th value
	    geom_integ->setUVW (puvw4);
	    pxyz4 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val4);	      

	    /// 5st value
	    geom_integ->setUVW (puvw5);
	    pxyz5 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val5);
	    
	    /// 6nd value
	    geom_integ->setUVW (puvw6);
	    pxyz6 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val6);
	    
	    /// 7rd value
	    geom_integ->setUVW (puvw7);
	    pxyz7 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val7);	      
	    
	    /// 8th value
	    geom_integ->setUVW (puvw8);
	    pxyz8 = geom_integ->getXYZ();
	    if (geom_appro->getEntity() != geom_integ->getEntity()) geom_appro->setUVWForXYZ(geom_integ->getXYZ());
	    else geom_appro->setUVW(geom_integ->getUVW());
	    eval(geom_appro, geom_integ, val8);	 
  
	    pexport.exportHex(pxyz1,pxyz2,pxyz3,pxyz4,pxyz5,pxyz6,pxyz7,pxyz8,val1,val2,val3,val4,val5,val6,val7,val8);

	  }
      }
  }
};



/// To export fields defined at gauss points
/// Still quite a bit experimental !
template <class T>
class xexport::xPlotPWCommand : public xCommandOnGeomElem
{
  
private:
  xexport::xExport& pexport;
  const xIntegrationRule& integ_rule;
  //const xEvalField<T,xFieldPointwise>& eval;
  const xEval<T>& eval;
  typedef typename T::result_type result_type;
  
public:
  //xexport::xPlotPWCommand (const xEvalField<T,xFieldPointwise>& eval_, xexport::xExport& e,const xIntegrationRule& integ_rule_)
  xexport::xPlotPWCommand (const xEval<T>& eval_, xexport::xExport& e,const xIntegrationRule& integ_rule_)
  : eval(eval_), pexport(e), integ_rule(integ_rule_)  {}
  
  void execute(xGeomElem* geom_integ)
  {
    //Il faut reconstruire un xGeomElem avec les bons PW
    AOMD::mEntity* e= geom_integ->getEntity();//Entite...
    xmapping::xMappingBuilderHolder& holder = xmapping::xMappingBuilderHolderSingleton::instance();
    Trellis_Util::Mapping* mapping = holder.getMappingBuilder()->BuildMapping(e);//Mapping...
    Trellis_Util::Integrator* integrator=integ_rule.getIntegrator(e);
    xGeomElem Ginteg(e,mapping,integrator);//...et ca nous donne le BON geomElem
    Ginteg=*geom_integ;//NORMALEMENT ON DEVRAIT POUVOIR EFFACER mais ca ne marche pas
    Ginteg.SetIntegrationPointNumberForDegree(geom_integ->getOrder());//On copie le nombre de pts d'integration en fonction du degres
    //cout<<Ginteg.GetNbIntegrationPoints()<<endl;
    
    //e->print();
    for(int k=0;k<Ginteg.GetNbIntegrationPoints();k++)//Pour chaque point d'integration
    {
      //cout<<k<<" / ";
      Ginteg.setUVW(k);//On dit quel point on veut
      Trellis_Util::mPoint p=Ginteg.getXYZ();//On en fait un point
      //typename T::result_type result_type;
      result_type val;
      //field_PW.getVal(geom_integ,geom_integ,val);
      if (geom_appro->getEntity() != Ginteg.getEntity()) geom_appro->setUVWForXYZ(Ginteg.getXYZ());
      else geom_appro->setUVW(Ginteg.getUVW());
                     
      //cout<<p(0)<<" - "<<p(1)<<" - "<<p(2)<<endl;
      eval(geom_appro, &Ginteg, val);
      pexport.exportPoint(p,val);
    }
    delete mapping;
  }
  
};


} // end of namespace



#endif
