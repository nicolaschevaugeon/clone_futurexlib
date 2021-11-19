/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef _export_ensight_H
#define _export_ensight_H

#include <fstream>
#include <iostream>
#include <list>

#include "xExport.h"

namespace xexport
{
// voir si il est dans xExport.

class xExportEnsight : public xExport
{
  public:
   // ClasseInterne pour l'enregistrement :

   class Figure
   {
     public:
      class Tensor
      {
        private:
         float tab_[9];

        public:
         Tensor();
         Tensor(float, float, float, float, float, float, float, float, float);
         Tensor(const Tensor&);
         float getTensor(int) const;
         void afficheTensor() const;
      };

      // classe interne coordonne
      class Coordon
      {
         // mettre la classe qui touche aux coordonnées en friends
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
         // faire surcharge de <<
         float getX() const;
         float getY() const;
         float getZ() const;

         //  ~Coordon();
      };

     private:
      std::string project_name_;  // nom du projet
      // int cptFigure_;
      std::string current_variable_name_;
      short typeElement_;
      short typeVariable_;

     public:
      std::list<Coordon> liste;  // Pour les coordonnées.
      /*pour la dependance du temps on vérifiera la longeur de la std::liste
        de Coordon par rapport à la taille des variables, et si on est pas égal, alors c'est en fonction du temps
        et on pourra trouver la dépendace ( nombre de pas de temps ) pa rapport a un calcul */
      std::list<float> scalaire;
      std::list<Coordon> vecteur;
      std::list<Tensor> tensor;

      Figure();
      Figure(std::string);  // on donne juste le nom du project
      // Figure(const Figure&);

      std::string getProjectName() const;
      short getTypeElem() const;
      // int getCptFigure() const;
      std::string getVariableName() const;
      void setVariableName(const std::string& s);
      short getTypeVariable() const;

      void addElem(float, float, float);
      void addScalar(float);
      void addVect(float, float, float);
      void addTensor(float, float, float, float, float, float, float, float, float);

      void addtypeVariable(short);
      void incrementCptFigure();

      void addTypeElem(short);
      void addNomFigure(const std::string&);

      void afficheToutCoordon() const;
      // classe interne tensor.

     public:
   };
   // FIn classe interne

   xExportEnsight(MPI_Comm world);

   void exportPoint(const xtensor::xPoint&, const double&) override
   {
      std::cout << "xExportEnsight::exportPoint(const xtensor::xPoint&, const double&) not implemented yet" << std::endl;
      throw;
   }

   void exportPoint(const xtensor::xPoint&, const xtensor::xVector<>&) override
   {
      std::cout << "xExportEnsight::exportPoint(const xtensor::xPoint&, const xtensor::xVector&) not implemented yet"
                << std::endl;
      throw;
   }

   void exportPoint(const xtensor::xPoint&, const xtensor::xTensor2<>&) override
   {
      std::cout << "xExportEnsight::exportPoint(const xtensor::xPoint&, const xtensor::xTensor2&) not implemented yet"
                << std::endl;
      throw;
   }

   void exportLine(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const double& val1, const double& val2) override;
   void exportLine(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xVector<>& val1,
                   const xtensor::xVector<>& val2) override;
   void exportLine(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xTensor2<>& val1,
                   const xtensor::xTensor2<>& val2) override;

   void exportTriangle(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const double& val1,
                       const double& val2, const double& val3) override;
   void exportTriangle(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3,
                       const xtensor::xVector<>& val1, const xtensor::xVector<>& val2, const xtensor::xVector<>& val3) override;
   void exportTriangle(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3,
                       const xtensor::xTensor2<>& val1, const xtensor::xTensor2<>& val2,
                       const xtensor::xTensor2<>& val3) override;

   void exportQuad(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const xtensor::xPoint& P4,
                   const double& val1, const double& val2, const double& val3, const double& val4) override;
   void exportQuad(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const xtensor::xPoint& P4,
                   const xtensor::xVector<>& val1, const xtensor::xVector<>& val2, const xtensor::xVector<>& val3,
                   const xtensor::xVector<>& val4) override;
   void exportQuad(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const xtensor::xPoint& P4,
                   const xtensor::xTensor2<>& val1, const xtensor::xTensor2<>& val2, const xtensor::xTensor2<>& val3,
                   const xtensor::xTensor2<>& val4) override;

   void exportTetra(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const xtensor::xPoint& P4,
                    const double& val1, const double& val2, const double& val3, const double& val4) override;
   void exportTetra(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const xtensor::xPoint& P4,
                    const xtensor::xVector<>& val1, const xtensor::xVector<>& val2, const xtensor::xVector<>& val3,
                    const xtensor::xVector<>& val4) override;
   void exportTetra(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const xtensor::xPoint& P4,
                    const xtensor::xTensor2<>& val1, const xtensor::xTensor2<>& val2, const xtensor::xTensor2<>& val3,
                    const xtensor::xTensor2<>& val4) override;
   void exportHex(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const xtensor::xPoint& P4,
                  const xtensor::xPoint& P5, const xtensor::xPoint& P6, const xtensor::xPoint& P7, const xtensor::xPoint& P8,
                  const double& val1, const double& val2, const double& val3, const double& val4, const double& val5,
                  const double& val6, const double& val7, const double& val8) override;
   void exportHex(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const xtensor::xPoint& P4,
                  const xtensor::xPoint& P5, const xtensor::xPoint& P6, const xtensor::xPoint& P7, const xtensor::xPoint& P8,
                  const xtensor::xVector<>& val1, const xtensor::xVector<>& val2, const xtensor::xVector<>& val3,
                  const xtensor::xVector<>& val4, const xtensor::xVector<>& val5, const xtensor::xVector<>& val6,
                  const xtensor::xVector<>& val7, const xtensor::xVector<>& val8) override;
   void exportHex(const xtensor::xPoint& P1, const xtensor::xPoint& P2, const xtensor::xPoint& P3, const xtensor::xPoint& P4,
                  const xtensor::xPoint& P5, const xtensor::xPoint& P6, const xtensor::xPoint& P7, const xtensor::xPoint& P8,
                  const xtensor::xTensor2<>& val1, const xtensor::xTensor2<>& val2, const xtensor::xTensor2<>& val3,
                  const xtensor::xTensor2<>& val4, const xtensor::xTensor2<>& val5, const xtensor::xTensor2<>& val6,
                  const xtensor::xTensor2<>& val7, const xtensor::xTensor2<>& val8) override;

  protected:
   Figure fig;
   std::ofstream fcase;
   std::ofstream fgeo;
   std::ofstream fscl;
   std::ofstream fvct;
   std::ofstream ftens;
   bool writeGeo;
   virtual void exportGeometry() = 0;
   virtual void exportVariableScl() = 0;
   virtual void exportVariableVect() = 0;
   virtual void exportVariableTensor() = 0;
   void exportVariableToEnsight();
   virtual void openLoadFic(std::ofstream&) = 0;

   void writeHeaderCase();
   std::string giveElemEnsight(short type);

   std::string giveExtension(short a);
   bool isLoad() const;
   void firstOn(short, short);
};

class xExportEnsightAscii : public xExportEnsight
{
  public:
   xExportEnsightAscii(MPI_Comm world = MPI_COMM_WORLD);
   ~xExportEnsightAscii() override;
   void startView(const std::string& comment) override;

   void endView() override;
   void openFile(const std::string& fName) override;
   void closeFile() override;

  protected:
   // mettre en virtual dans xExportEnsight.
   void exportGeometry() override;
   void exportVariableScl() override;
   void exportVariableVect() override;
   void exportVariableTensor() override;
   void openLoadFic(std::ofstream&) override;
};

class xExportEnsightBinary : public xExportEnsight
{
  public:
   xExportEnsightBinary(MPI_Comm world = MPI_COMM_WORLD);
   ~xExportEnsightBinary() override;
   void startView(const std::string& comment) override;

   void endView() override;
   void openFile(const std::string& fName) override;
   void closeFile() override;

  protected:
   // mettre en virtual dans xExportEnsight.
   void exportGeometry() override;
   void exportVariableScl() override;
   void exportVariableVect() override;
   void exportVariableTensor() override;
   void openLoadFic(std::ofstream&) override;
};

}  // namespace xexport
#endif
