/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _export_ensight_H
#define _export_ensight_H

#include <iostream>
#include <fstream>
#include <list>
#include "xExport.h"
#include "mPoint.h"



namespace xfem
{

//voir si il est dans xexport::xExport.

using std::list;
  using std::ofstream;
using Trellis_Util::mPoint;

 

 class xexport::xExportEnsight : public xexport::xExport
{
 public:


  //ClasseInterne pour l'enregistrement :


  class Figure{

  public:
    class Tensor{
    private:
      float tab_[9];
    public:
      Tensor();
      Tensor(float, float, float, float, float, float, float, float, float);
      Tensor(const Tensor&);
      float getTensor(int)const;
      void afficheTensor() const;
      
    };
    
    //classe interne coordonne
    class Coordon{
      //mettre la classe qui touche aux coordonnées en friends
    private:
      float x_;
      float y_;
      float z_;
      
      void setX(float a);
      void setY(float a);
      void setZ(float a);
      
    public:
      Coordon();
      Coordon(float, float, float);
      Coordon(const Coordon&);
    
      void afficheCoordon() const;
      //faire surcharge de <<
      float getX() const;
      float getY() const;
      float getZ() const;
  
      //  ~Coordon();
    };
    
  private:
    string project_name_; //nom du projet
    //int cptFigure_;
    string current_variable_name_;
    short typeElement_;
    short typeVariable_;
    
  public :
    list<Coordon> liste;//Pour les coordonnées.
    /*pour la dependance du temps on vérifiera la longeur de la liste
      de Coordon par rapport à la taille des variables, et si on est pas égal, alors c'est en fonction du temps
      et on pourra trouver la dépendace ( nombre de pas de temps ) pa rapport a un calcul */
    list<float> scalaire;
    list<Coordon> vecteur;
    list<Tensor> tensor;
    
    Figure();
    Figure(string);// on donne juste le nom du project
    //Figure(const Figure&);
    
    string getProjectName() const;
    short getTypeElem() const;
    //int getCptFigure() const;
    string getVariableName() const;
    void setVariableName(const string& s);
    short getTypeVariable() const;
    
    
    void addElem(float, float, float);
    void addScalar(float);
    void addVect(float, float, float);
    void addTensor(float, float, float, float, float, float, float, float, float);
    
    void addtypeVariable(short);
    void incrementCptFigure();

    void addTypeElem(short);
    void addNomFigure(const string&);
    
    
    void afficheToutCoordon() const;
    //classe interne tensor.


  public:
    
  };
  //FIn classe interne
  

  xexport::xExportEnsight():xexport::xExport(), writeGeo(false) {};
  
  void exportPoint(const Trellis_Util::mPoint&, const double&){
    std::cout << "xexport::xExportEnsight::exportPoint(const Trellis_Util::mPoint&, const double&) not implemented yet" << std::endl;
    throw;
  }
  
  void exportPoint(const Trellis_Util::mPoint&, const xtensor::xVector&){
    std::cout << "xexport::xExportEnsight::exportPoint(const Trellis_Util::mPoint&, const xtensor::xVector&) not implemented yet" << std::endl;
    throw;
  }

  void exportPoint(const Trellis_Util::mPoint&, const xtensor::xTensor2&){
    std::cout << "xexport::xExportEnsight::exportPoint(const Trellis_Util::mPoint&, const xtensor::xTensor2&) not implemented yet" << std::endl;
    throw;
  }

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
  Figure fig;
  ofstream fcase;
  ofstream fgeo;
  ofstream fscl;
  ofstream fvct;
  ofstream ftens;
  bool writeGeo;
  virtual void exportGeometry()=0;
  virtual void exportVariableScl()=0;
  virtual void exportVariableVect()=0;
  virtual void exportVariableTensor()=0;
  void exportVariableToEnsight();
  virtual void openLoadFic(std::ofstream&)=0;

  void writeHeaderCase();
  string giveElemEnsight(short type);
  
  string giveExtension(short a);
  bool isLoad() const;
  void firstOn(short, short);
  

};


  class xexport::xExportEnsightAscii : public xexport::xExportEnsight
  {
  public:
    
    virtual ~xexport::xExportEnsightAscii();
    void startView(const string& comment);
    
    void endView ();
    void openFile (const string& fName);
    void closeFile ();
  protected:
    //mettre en virtual dans xexport::xExportEnsight.
    void exportGeometry();
    void exportVariableScl();
    void exportVariableVect();
    void exportVariableTensor();
    void openLoadFic(std::ofstream&);    
  };

 class xexport::xExportEnsightBinary : public xexport::xExportEnsight
{
  public:
    
    virtual ~xexport::xExportEnsightBinary();
    void startView(const string& comment);
    
    void endView ();
    void openFile (const string& fName);
    void closeFile ();
  protected:
    //mettre en virtual dans xexport::xExportEnsight.
    void exportGeometry();
    void exportVariableScl();
    void exportVariableVect();
    void exportVariableTensor();
    void openLoadFic(std::ofstream&);
};


}
#endif
