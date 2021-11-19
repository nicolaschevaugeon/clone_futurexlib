/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include "xExportEnsight.h"

#include <cassert>
#include <fstream>

#include "workInProgress.h"
#include "xPoint.h"

#define CST_CHAR sizeof(char) * 80
#define EXT_SCL ".scl"
#define EXT_VCT ".vct"
#define EXT_TENS ".tens"
#define EXT_CASE ".case"
#define EXT_GEO ".geo"
//#define FORMAT_ASCII_GEO fgeo.fill(' ');fgeo.width(10);

namespace xexport
{
using xtensor::xPoint;
using namespace std;

xExportEnsight::xExportEnsight(MPI_Comm world) : xExport(world), writeGeo(false)
{
   // every thing must be checked
   assert(xtool::workInProgress());
}

bool xExportEnsight::isLoad() const { return fig.getTypeVariable() == 0; }

void xExportEnsight::firstOn(short var, short elem)
{
   fig.addtypeVariable(var);
   fig.addTypeElem(elem);
}

string xExportEnsight::giveElemEnsight(short type)
{
   switch (type)
   {
      case 4:
         return "tetra4";
      case 3:
         return "tria3";
      case 1:
         return "point";
      case 2:
         return "bar2";
      case 5:
         return "quad4";
      case 8:
         return "hexa8";
      default:
         throw;
   }
}

string xExportEnsight::giveExtension(short a)
{
   switch (a)
   {
      case 1:
         return EXT_SCL;
      case 2:
         return EXT_VCT;
      case 3:
         return EXT_TENS;
      default:
         throw;
   }
}

void xExportEnsight::exportLine(const xPoint &P1, const xPoint &P2, const double &val1, const double &val2)
{
   if (isLoad()) firstOn(1, 2);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
   }
   fig.addScalar(val1);
   fig.addScalar(val2);
}

void xExportEnsight::exportLine(const xPoint &P1, const xPoint &P2, const xtensor::xVector<> &val1,
                                const xtensor::xVector<> &val2)
{
   if (isLoad()) firstOn(2, 2);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
   }
   fig.addVect(val1(0), val1(1), val1(2));
   fig.addVect(val2(0), val2(1), val2(2));
}

void xExportEnsight::exportLine(const xPoint &P1, const xPoint &P2, const xtensor::xTensor2<> &val1,
                                const xtensor::xTensor2<> &val2)
{
   if (isLoad()) firstOn(3, 2);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
   }
   fig.addTensor(val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2));
   fig.addTensor(val2(0, 0), val2(0, 1), val2(0, 2), val2(1, 0), val2(1, 1), val2(1, 2), val2(2, 0), val2(2, 1), val2(2, 2));
}

void xExportEnsight::exportTriangle(const xPoint &P1, const xPoint &P2, const xPoint &P3, const double &val1, const double &val2,
                                    const double &val3)
{
   if (isLoad()) firstOn(1, 3);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
      fig.addElem(P3(0), P3(1), P3(2));
   }
   fig.addScalar(val1);
   fig.addScalar(val2);
   fig.addScalar(val3);
}
void xExportEnsight::exportTriangle(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xtensor::xVector<> &val1,
                                    const xtensor::xVector<> &val2, const xtensor::xVector<> &val3)
{
   if (isLoad()) firstOn(2, 3);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
      fig.addElem(P3(0), P3(1), P3(2));
   }
   fig.addVect(val1(0), val1(1), val1(2));
   fig.addVect(val2(0), val2(1), val2(2));
   fig.addVect(val3(0), val3(1), val3(2));
}
void xExportEnsight::exportTriangle(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xtensor::xTensor2<> &val1,
                                    const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3)
{
   if (isLoad()) firstOn(3, 3);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
      fig.addElem(P3(0), P3(1), P3(2));
   }
   fig.addTensor(val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2));
   fig.addTensor(val2(0, 0), val2(0, 1), val2(0, 2), val2(1, 0), val2(1, 1), val2(1, 2), val2(2, 0), val2(2, 1), val2(2, 2));
   fig.addTensor(val3(0, 0), val3(0, 1), val3(0, 2), val3(1, 0), val3(1, 1), val3(1, 2), val3(2, 0), val3(2, 1), val3(2, 2));
}

void xExportEnsight::exportQuad(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const double &val1,
                                const double &val2, const double &val3, const double &val4)
{
   if (isLoad()) firstOn(1, 5);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
      fig.addElem(P3(0), P3(1), P3(2));
      fig.addElem(P4(0), P4(1), P4(2));
   }
   fig.addScalar(val1);
   fig.addScalar(val2);
   fig.addScalar(val3);
   fig.addScalar(val4);
}
void xExportEnsight::exportQuad(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4,
                                const xtensor::xVector<> &val1, const xtensor::xVector<> &val2, const xtensor::xVector<> &val3,
                                const xtensor::xVector<> &val4)
{
   if (isLoad()) firstOn(2, 5);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
      fig.addElem(P3(0), P3(1), P3(2));
      fig.addElem(P4(0), P4(1), P4(2));
   }
   fig.addVect(val1(0), val1(1), val1(2));
   fig.addVect(val2(0), val2(1), val2(2));
   fig.addVect(val3(0), val3(1), val3(2));
   fig.addVect(val4(0), val4(1), val4(2));
}
void xExportEnsight::exportQuad(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4,
                                const xtensor::xTensor2<> &val1, const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3,
                                const xtensor::xTensor2<> &val4)
{
   if (isLoad()) firstOn(3, 5);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
      fig.addElem(P3(0), P3(1), P3(2));
      fig.addElem(P4(0), P4(1), P4(2));
   }
   fig.addTensor(val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2));
   fig.addTensor(val2(0, 0), val2(0, 1), val2(0, 2), val2(1, 0), val2(1, 1), val2(1, 2), val2(2, 0), val2(2, 1), val2(2, 2));
   fig.addTensor(val3(0, 0), val3(0, 1), val3(0, 2), val3(1, 0), val3(1, 1), val3(1, 2), val3(2, 0), val3(2, 1), val3(2, 2));
   fig.addTensor(val4(0, 0), val4(0, 1), val4(0, 2), val4(1, 0), val4(1, 1), val4(1, 2), val4(2, 0), val4(2, 1), val4(2, 2));
}

void xExportEnsight::exportTetra(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const double &val1,
                                 const double &val2, const double &val3, const double &val4)
{
   if (isLoad()) firstOn(1, 4);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
      fig.addElem(P3(0), P3(1), P3(2));
      fig.addElem(P4(0), P4(1), P4(2));
   }
   fig.addScalar(val1);
   fig.addScalar(val2);
   fig.addScalar(val3);
   fig.addScalar(val4);
}
void xExportEnsight::exportTetra(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4,
                                 const xtensor::xVector<> &val1, const xtensor::xVector<> &val2, const xtensor::xVector<> &val3,
                                 const xtensor::xVector<> &val4)
{
   if (isLoad()) firstOn(2, 4);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
      fig.addElem(P3(0), P3(1), P3(2));
      fig.addElem(P4(0), P4(1), P4(2));
   }
   fig.addVect(val1(0), val1(1), val1(2));
   fig.addVect(val2(0), val2(1), val2(2));
   fig.addVect(val3(0), val3(1), val3(2));
   fig.addVect(val4(0), val4(1), val4(2));
}
void xExportEnsight::exportTetra(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4,
                                 const xtensor::xTensor2<> &val1, const xtensor::xTensor2<> &val2,
                                 const xtensor::xTensor2<> &val3, const xtensor::xTensor2<> &val4)
{
   if (isLoad()) firstOn(3, 4);
   if (writeGeo)
   {
      fig.addElem(P1(0), P1(1), P1(2));
      fig.addElem(P2(0), P2(1), P2(2));
      fig.addElem(P3(0), P3(1), P3(2));
      fig.addElem(P4(0), P4(1), P4(2));
   }
   fig.addTensor(val1(0, 0), val1(0, 1), val1(0, 2), val1(1, 0), val1(1, 1), val1(1, 2), val1(2, 0), val1(2, 1), val1(2, 2));
   fig.addTensor(val2(0, 0), val2(0, 1), val2(0, 2), val2(1, 0), val2(1, 1), val2(1, 2), val2(2, 0), val2(2, 1), val2(2, 2));
   fig.addTensor(val3(0, 0), val3(0, 1), val3(0, 2), val3(1, 0), val3(1, 1), val3(1, 2), val3(2, 0), val3(2, 1), val3(2, 2));
   fig.addTensor(val4(0, 0), val4(0, 1), val4(0, 2), val4(1, 0), val4(1, 1), val4(1, 2), val4(2, 0), val4(2, 1), val4(2, 2));
}
void xExportEnsight::exportHex(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xPoint &P5,
                               const xPoint &P6, const xPoint &P7, const xPoint &P8, const double &val1, const double &val2,
                               const double &val3, const double &val4, const double &val5, const double &val6, const double &val7,
                               const double &val8)
{
   cout << "Not Coded\n";
   assert(0);
}
void xExportEnsight::exportHex(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xPoint &P5,
                               const xPoint &P6, const xPoint &P7, const xPoint &P8, const xtensor::xVector<> &val1,
                               const xtensor::xVector<> &val2, const xtensor::xVector<> &val3, const xtensor::xVector<> &val4,
                               const xtensor::xVector<> &val5, const xtensor::xVector<> &val6, const xtensor::xVector<> &val7,
                               const xtensor::xVector<> &val8)
{
   cout << "Not Coded\n";
   assert(0);
}
void xExportEnsight::exportHex(const xPoint &P1, const xPoint &P2, const xPoint &P3, const xPoint &P4, const xPoint &P5,
                               const xPoint &P6, const xPoint &P7, const xPoint &P8, const xtensor::xTensor2<> &val1,
                               const xtensor::xTensor2<> &val2, const xtensor::xTensor2<> &val3, const xtensor::xTensor2<> &val4,
                               const xtensor::xTensor2<> &val5, const xtensor::xTensor2<> &val6, const xtensor::xTensor2<> &val7,
                               const xtensor::xTensor2<> &val8)
{
   cout << "Not Coded\n";
   assert(0);
}

void xExportEnsight::writeHeaderCase()
{
   string project_name(fig.getProjectName());
   fcase << "#" << endl << "# " << project_name << endl << "# " << endl;
   fcase << "FORMAT" << endl;
   fcase << "type:\t\t\tensight gold" << endl;
   fcase << endl;
   fcase << "GEOMETRY" << endl;
   fcase << "model:\t\t\t" << project_name << ".geo" << endl;

   short typeV = fig.getTypeVariable();
   if (typeV != 0)
   {
      fcase << endl << "VARIABLE" << endl;
      if (typeV == 1) fcase << "scalar per node";
      if (typeV == 2) fcase << "vector per node";
      if (typeV == 3) fcase << "tensor asym per node";

      fcase << ":\t" << fig.getVariableName() << "\t" << project_name;
      if (project_name != fig.getVariableName()) fcase << "_" << fig.getVariableName();
      fcase << giveExtension(fig.getTypeVariable()) << endl;
   }
}

xExportEnsightAscii::xExportEnsightAscii(MPI_Comm world) : xExportEnsight(world) {}
// Fonction qui ouvre un fichier en fonction de sa variable
void xExportEnsightAscii::openLoadFic(ofstream &fic)
{
   string nom(fig.getProjectName());

   if (fig.getVariableName() != "" && fig.getVariableName() != fig.getProjectName())
   {
      nom += "_";
      nom += fig.getVariableName();
   }
   nom += giveExtension(fig.getTypeVariable());

   fic.open(nom.c_str(), fstream::out);
}

void xExportEnsightAscii::exportGeometry()
{
   fgeo << "fichier : " << fig.getProjectName() << endl;
   fgeo << endl;
   fgeo << "node id assign" << endl;
   fgeo << "element id assign" << endl;
   fgeo << "part" << endl;
   fgeo.fill(' ');
   fgeo.width(10);
   fgeo << 1 << endl;
   fgeo << "partie 1" << endl;
   fgeo << "coordinates" << endl;
   fgeo.fill(' ');
   fgeo.width(10);
   int tei = fig.liste.size();
   fgeo << tei << endl;
   // on rentre les coordonees x puis y puis z;
   int cpt;
   list<Figure::Coordon>::const_iterator i;
   float var;
   for (cpt = 1; cpt < 4; cpt++)
   {
      i = fig.liste.begin();
      for (; i != fig.liste.end(); ++i)
      {
         if (cpt == 1) var = i->getX();
         if (cpt == 2) var = i->getY();
         if (cpt == 3) var = i->getZ();
         fgeo.setf(ios::scientific, ios::floatfield);
         fgeo.precision(5);
         fgeo.fill(' ');
         fgeo.width(12);
         fgeo << var << endl;
      }
   }
   short typeEnsight = fig.getTypeElem();
   string nomElem(giveElemEnsight(typeEnsight));
   fgeo << nomElem << endl;
   if (typeEnsight == 5) typeEnsight--;
   fgeo.fill(' ');
   fgeo.width(10);
   int taille = fig.liste.size();
   fgeo << (taille / typeEnsight) << endl;
   int connect = 0;
   int compt = 1;

   while (connect < (taille))
   {
      for (compt = 0; compt < typeEnsight; compt++)
      {
         fgeo.fill(' ');
         fgeo.width(10);
         fgeo << ++connect;
      }
      fgeo << endl;
   }
   fgeo.close();
   writeGeo = false;
}

void xExportEnsightAscii::openFile(const string &fName)
{
   // on ouvre tous les fichier
   string nfcase = fName + EXT_CASE;
   string nfgeo = fName + EXT_GEO;
   fcase.open(nfcase.string::c_str());
   fgeo.open(nfgeo.string::c_str());
   writeGeo = true;
   fig.addNomFigure(fName);
   process_started = true;
}

void xExportEnsightAscii::startView(const string &variable_name)
{
   // ici la figure n'est pas rempi donc j'ecris juste un bout du fichier case
   fig.setVariableName(variable_name);
   // cout<<"Chargement de la figure en memoire"<<endl;
   // ou
   //  cout<<"cahrgement d'une nouvelle variable de la figure"<<endl;
}

void xExportEnsightAscii::endView()
{
   // cout<<"debut de l'ecriture"<<endl;
   // writeHeaderCase();
   //    exportGeometry();
   if (writeGeo)
   {
      // Ici on fait les autres traitement.
      // Export de la geometrie
      writeHeaderCase();
      exportGeometry();
      // on commence l'export des variables
      // cout<<"******************************************exportvariable"<<endl;
      exportVariableToEnsight();
      fig.liste.clear();
   }
   else
   {
      // cout<<"je passe la par contre"<<endl;
      // on exporte la variable avec uniquement une nouvelle ligne au .case.
      short typeV = fig.getTypeVariable();
      if (typeV == 1) fcase << "scalar per node";
      if (typeV == 2) fcase << "vector per node";
      if (typeV == 3) fcase << "tensor asym per node";

      fcase << ":\t" << fig.getVariableName() << "\t" << fig.getProjectName() << "_" << fig.getVariableName()
            << giveExtension(fig.getTypeVariable()) << endl;

      exportVariableToEnsight();
   }
   // ici on met a zero les variables et la geo
   // vidage de la figure

   fig.scalaire.clear();
   fig.vecteur.clear();
   fig.tensor.clear();
   // mise a zero de la variable pour une seconde si necessaire.
   fig.addtypeVariable(0);
   // fig.addTypeElem(0);
}

void xExportEnsightAscii::closeFile()
{
   fcase.close();
   process_started = false;
}

xExportEnsightAscii::~xExportEnsightAscii()
{
   // on detruit le reste des liste
   // mais deja fait...
}

// varIables
void xExportEnsight::exportVariableToEnsight()
{
   switch (fig.getTypeVariable())
   {
      case 1:
         exportVariableScl();
         break;
      case 2:
         exportVariableVect();
         break;
      case 3:
         exportVariableTensor();
         break;
   }
}

void xExportEnsightAscii::exportVariableTensor()
{
   openLoadFic(ftens);
   ftens << "tensor" << endl;
   ftens << "part" << endl;
   ftens.fill(' ');
   ftens.width(10);

   ftens << 1 << endl;
   ftens << "coordinates" << endl;
   int cpt;
   list<Figure::Tensor>::const_iterator i;

   for (cpt = 0; cpt < 9; cpt++)
   {
      // cout<<cpt<<"  ";
      i = fig.tensor.begin();
      for (; i != fig.tensor.end(); ++i)
      {
         ftens.setf(ios::scientific, ios::floatfield);
         ftens.precision(5);
         ftens.fill(' ');
         ftens.width(12);

         ftens << i->Figure::Tensor::getTensor(cpt) << endl;
      }
   }
   ftens.close();
}

void xExportEnsightAscii::exportVariableVect()
{
   openLoadFic(fvct);
   fvct << "vecteur" << endl;
   fvct << "part" << endl;
   fvct.fill(' ');
   fvct.width(10);

   fvct << 1 << endl;
   fvct << "coordinates" << endl;
   int cpt;
   list<Figure::Coordon>::const_iterator i;
   float var;
   for (cpt = 1; cpt < 4; cpt++)
   {
      i = fig.vecteur.begin();
      for (; i != fig.vecteur.end(); ++i)
      {
         if (cpt == 1) var = i->Figure::Coordon::getX();
         if (cpt == 2) var = i->Figure::Coordon::getY();
         if (cpt == 3) var = i->Figure::Coordon::getZ();

         fvct.setf(ios::scientific, ios::floatfield);
         fvct.precision(5);
         fvct.fill(' ');
         fvct.width(12);

         fvct << var << endl;
      }
   }
   fvct.close();
}

void xExportEnsightAscii::exportVariableScl()
{
   openLoadFic(fscl);
   fscl << "scalaire" << endl;
   fscl << "part" << endl;
   fscl.fill(' ');
   fscl.width(10);
   fscl << 1 << endl;
   fscl << "coordinates" << endl;
   list<float>::const_iterator i = fig.scalaire.begin();
   for (; i != fig.scalaire.end(); ++i)
   {
      fscl.setf(ios::scientific, ios::floatfield);
      fscl.precision(5);
      fscl.fill(' ');
      fscl.width(12);
      // fscl<<i->getSc()<<endl;
      fscl << *i << endl;
   }
   fscl.close();
}

// Debut implementation classe interne
xExportEnsight::Figure::Figure() : project_name_(""), current_variable_name_(""), typeElement_(0), typeVariable_(0) {}

xExportEnsight::Figure::Figure(string nom) : project_name_(nom), current_variable_name_(""), typeElement_(0), typeVariable_(0) {}

string xExportEnsight::Figure::getProjectName() const { return project_name_; }

short xExportEnsight::Figure::getTypeElem() const { return typeElement_; }

string xExportEnsight::Figure::getVariableName() const { return current_variable_name_; }

void xExportEnsight::Figure::setVariableName(const string &s) { current_variable_name_ = s; }

short xExportEnsight::Figure::getTypeVariable() const { return typeVariable_; }

void xExportEnsight::Figure::addElem(float x, float y, float z)
{
   list<Figure::Coordon>::iterator iter = liste.begin();
   liste.insert(iter, Coordon(x, y, z));
}

void xExportEnsight::Figure::addScalar(float s)
{
   list<float>::iterator iter = scalaire.begin();
   scalaire.insert(iter, s);
}

void xExportEnsight::Figure::addVect(float x, float y, float z)
{
   list<Figure::Coordon>::iterator iter = vecteur.begin();
   vecteur.insert(iter, Coordon(x, y, z));
}

void xExportEnsight::Figure::addTensor(float a, float z, float e, float r, float t, float y, float u, float i, float o)
{
   list<Figure::Tensor>::iterator iter = tensor.begin();
   tensor.insert(iter, Tensor(a, z, e, r, t, y, u, i, o));
}

void xExportEnsight::Figure::addTypeElem(short bozo) { typeElement_ = bozo; }
void xExportEnsight::Figure::addtypeVariable(short bozo) { typeVariable_ = bozo; }

void xExportEnsight::Figure::addNomFigure(const string &titi) { project_name_ = titi; }

void xExportEnsight::Figure::afficheToutCoordon() const
{
   list<Figure::Coordon>::const_iterator i = liste.begin();
   // oli::iterator il = l.end();
   for (; i != liste.end(); ++i)
   {
      i->afficheCoordon();
      // cout<<*i<<" ";
   }
   cout << endl;
}

// Debut classe interne Tensor
xExportEnsight::Figure::Tensor::Tensor()
{
   int i;
   for (i = 0; i < 9; i++)
   {
      tab_[i] = 0;
   }
}
xExportEnsight::Figure::Tensor::Tensor(float a, float z, float e, float r, float t, float y, float u, float i, float o)
{
   tab_[0] = a;
   tab_[1] = z;
   tab_[2] = e;
   tab_[3] = r;
   tab_[4] = t;
   tab_[5] = y;
   tab_[6] = u;
   tab_[7] = i;
   tab_[8] = o;
   // tab_[9]={a,z,e,r,t,y,u,i,o};
}

xExportEnsight::Figure::Tensor::Tensor(const Tensor &ct)
{
   int i;
   for (i = 0; i < 9; i++)
   {
      tab_[i] = ct.tab_[i];
   }
}
float xExportEnsight::Figure::Tensor::getTensor(int i) const { return tab_[i]; }
void xExportEnsight::Figure::Tensor::afficheTensor() const
{
   int i;
   for (i = 0; i < 9; i++)
   {
      cout << i << " : " << tab_[i];
   }
   cout << endl;
}

// classe interne Coordon
xExportEnsight::Figure::Coordon::Coordon() : x_(0), y_(0), z_(0) {}

xExportEnsight::Figure::Coordon::Coordon(float cx, float cy, float cz) : x_(cx), y_(cy), z_(cz) {}
xExportEnsight::Figure::Coordon::Coordon(const Coordon &ct) : x_(ct.x_), y_(ct.y_), z_(ct.z_)
{
   // x_(ct.x_);
   //  y_(ct.y_);
   //  z_(ct.z_);
}

float xExportEnsight::Figure::Coordon::getX() const { return x_; }
float xExportEnsight::Figure::Coordon::getY() const { return y_; }
float xExportEnsight::Figure::Coordon::getZ() const { return z_; }

// Coordon::~Coordon(){
// delete this;

//}
void xExportEnsight::Figure::Coordon::setX(float a) { x_ = a; }
void xExportEnsight::Figure::Coordon::setY(float a) { y_ = a; }
void xExportEnsight::Figure::Coordon::setZ(float a) { z_ = a; }

void xExportEnsight::Figure::Coordon::afficheCoordon() const
{
   std::cout << "x : " << this->getX() << " y : " << this->getY() << " z : " << this->getZ() << std::endl;
}

//************** Begin def class xExportEnsightBinary ******************

xExportEnsightBinary::xExportEnsightBinary(MPI_Comm world) : xExportEnsight(world) {}
void xExportEnsightBinary::openLoadFic(ofstream &fic)
{
   string nom(fig.getProjectName());

   if (fig.getVariableName() != "" && fig.getVariableName() != fig.getProjectName())
   {
      nom += "_";
      nom += fig.getVariableName();
   }
   nom += giveExtension(fig.getTypeVariable());

   fic.open(nom.c_str(), fstream::out | fstream::binary);
}

void xExportEnsightBinary::openFile(const string &fName)
{
   // on ouvre tous les fichier
   string nfcase = fName + EXT_CASE;
   string nfgeo = fName + EXT_GEO;
   fcase.open(nfcase.string::c_str());
   fgeo.open(nfgeo.string::c_str(), fstream::out | fstream::binary);
   writeGeo = true;
   fig.addNomFigure(fName);
   process_started = true;
}

void xExportEnsightBinary::startView(const string &variable_name)
{
   // ici la figure n'est pas rempi donc j'ecris juste un bout du fichier case
   fig.setVariableName(variable_name);

   // cout<<"Chargement de la figure en memoire"<<endl;
   // ou
   //  cout<<"cahrgement d'une nouvelle variable de la figure"<<endl;
}
void xExportEnsightBinary::endView()
{
   // cout<<"debut de l'ecriture"<<endl;
   // writeHeaderCase();
   //    exportGeometry();
   if (writeGeo)
   {
      // Ici on fait les autres traitement.
      // Export de la geometrie
      writeHeaderCase();
      exportGeometry();
      // on commence l'export des variables
      // cout<<"******************************************exportvariable"<<endl;
      exportVariableToEnsight();
      fig.liste.clear();
   }
   else
   {
      // cout<<"je passe la par contre"<<endl;
      // on exporte la variable avec uniquement une nouvelle ligne au .case.
      short typeV = fig.getTypeVariable();
      if (typeV == 1) fcase << "scalar per node";
      if (typeV == 2) fcase << "vector per node";
      if (typeV == 3) fcase << "tensor asym per node";

      fcase << ":\t" << fig.getVariableName() << "\t" << fig.getProjectName() << "_" << fig.getVariableName()
            << giveExtension(fig.getTypeVariable()) << endl;

      exportVariableToEnsight();
   }
   // ici on met a zero les variables et la geo
   // vidage de la figure

   fig.scalaire.clear();
   fig.vecteur.clear();
   fig.tensor.clear();
   // mise a zero de la variable pour une seconde si necessaire.
   fig.addtypeVariable(0);
   // fig.addTypeElem(0);
}

void xExportEnsightBinary::closeFile()
{
   fcase.close();
   process_started = false;
}

xExportEnsightBinary::~xExportEnsightBinary()
{
   // on detruit le reste des liste
   // mais deja fait...
}

void xExportEnsightBinary::exportGeometry()
{
   // On va creer des buffer
   char bufChar[80];
   memset(&bufChar, '\0', 80);
   int bufInt;
   float bufFloat;
   memcpy(&bufChar, "C Binary", 8);
   fgeo.write((char *)&bufChar, CST_CHAR);
   memcpy(&bufChar, "fichier : ", 10);
   fgeo.write((char *)&bufChar, 10);
   fgeo.write((fig.getProjectName()).c_str(), sizeof(char) * 70);
   fgeo.write("", CST_CHAR);
   fgeo.write("node id assign\n", CST_CHAR);
   fgeo.write("element id assign", CST_CHAR);

   fgeo.write("part", CST_CHAR);
   bufInt = 1;
   fgeo.write((char *)&bufInt, sizeof(int));

   fgeo.write("partie 1", CST_CHAR);
   fgeo.write("coordinates", CST_CHAR);
   bufInt = fig.liste.size();
   fgeo.write((char *)&bufInt, sizeof(int));

   int cpt;
   list<Figure::Coordon>::const_iterator i;
   for (cpt = 1; cpt < 4; cpt++)
   {
      i = fig.liste.begin();
      for (; i != fig.liste.end(); ++i)
      {
         if (cpt == 1) bufFloat = i->getX();
         if (cpt == 2) bufFloat = i->getY();
         if (cpt == 3) bufFloat = i->getZ();
         fgeo.write((char *)&bufFloat, sizeof(float));
      }
   }
   short typeEnsight = fig.getTypeElem();
   memset(&bufChar, '\0', 80);
   memcpy(&bufChar, giveElemEnsight(typeEnsight).string::c_str(), giveElemEnsight(typeEnsight).length());
   fgeo.write((char *)&bufChar, CST_CHAR);
   if (typeEnsight == 5) typeEnsight--;
   int taille = fig.liste.size();
   bufInt = taille / typeEnsight;
   fgeo.write((char *)&bufInt, sizeof(int));
   int connect = 0;
   int compt = 1;

   while (connect < (taille))
   {
      for (compt = 0; compt < typeEnsight; compt++)
      {
         ++connect;
         fgeo.write((char *)&connect, sizeof(int));
      }
   }

   fgeo.close();
   writeGeo = false;
}

void xExportEnsightBinary::exportVariableTensor()
{
   openLoadFic(ftens);
   float ff;
   int bufInt;
   ftens.write("tensor", CST_CHAR);
   ftens.write("part", CST_CHAR);
   bufInt = 1;
   ftens.write((char *)&bufInt, sizeof(int));
   ftens.write("coordinates", CST_CHAR);
   int cpt;
   list<Figure::Tensor>::const_iterator i;
   for (cpt = 0; cpt < 9; cpt++)
   {
      i = fig.tensor.begin();
      for (; i != fig.tensor.end(); ++i)
      {
         ff = i->Figure::Tensor::getTensor(cpt);
         ftens.write((char *)&ff, sizeof(float));
      }
   }
   ftens.close();
}

void xExportEnsightBinary::exportVariableVect()
{
   openLoadFic(fvct);
   float ff;
   fvct.write("vecteur", CST_CHAR);
   fvct.write("part", CST_CHAR);
   int buff;
   buff = 1;
   fvct.write((char *)&buff, sizeof(int));
   fvct.write("coordinates", CST_CHAR);
   int cpt;
   list<Figure::Coordon>::const_iterator i;
   for (cpt = 1; cpt < 4; cpt++)
   {
      i = fig.vecteur.begin();
      for (; i != fig.vecteur.end(); ++i)
      {
         if (cpt == 1) ff = i->Figure::Coordon::getX();
         if (cpt == 2) ff = i->Figure::Coordon::getY();
         if (cpt == 3) ff = i->Figure::Coordon::getZ();

         fvct.write((char *)&ff, sizeof(float));
      }
   }
   fvct.close();
}

void xExportEnsightBinary::exportVariableScl()
{
   openLoadFic(fscl);
   fscl.write("scalaire", CST_CHAR);
   fscl.write("part", CST_CHAR);
   int buff = 1;
   float bufFloat;
   fscl.write((char *)&buff, sizeof(int));
   fscl.write("coordinates", CST_CHAR);
   list<float>::const_iterator i = fig.scalaire.begin();
   for (; i != fig.scalaire.end(); ++i)
   {
      bufFloat = *i;
      fscl.write((char *)&bufFloat, sizeof(float));
   }
   fscl.close();
}

}  // namespace xexport
