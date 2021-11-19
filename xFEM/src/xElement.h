/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XELEMENT_H
#define _XELEMENT_H

#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

// aomd
#include "AOMD_Internals.h"
#include "mEntity.h"
// xtool
#include "xSingleton.h"
// xtensor
#include "xTensor2.h"
#include "xVector.h"

namespace xfem
{
class xElementInfoManager;
typedef xtool::xSingleton<xElementInfoManager> xElementInfoManagerSingleton;

// this a simplex element
// We assume jac and detjac constant over the element
class xElement
{
  public:
   //  ... Constructor
   xElement(const AOMD::mEntity* e);
   // const mis
   //  ... calcul d un champ scalaire interpole avec les fonctions de formes
   double getInterpoSca(const std::vector<double>& scas);
   //  ... calcul d un champ vectoriel interpole avec les fonctions des formes
   xtensor::xVector<> getInterpoVec(const std::vector<xtensor::xVector<>>& vecs);
   //  ... calcul du gradient d un champ scalaire interpole avec les fonctions de formes
   xtensor::xVector<> getGradInterpoSca(const std::vector<double>& scas);

   //  ... calcul du gradient d un champ vectoriel interpole avec les fonctions de formes

   void getGradInterpoVec(const std::vector<xtensor::xVector<>>& vecs, xtensor::xTensor2<>& ten);

   // fin const mis

   //  ... determination des coordonnees locales, de la matrice jacobienne, de son inverse
   //      et de la transposee de son inverse
   void xyz2uvw(const xtensor::xPoint& xyz);
   // JE NE SAIS PLUS A QUOI SERT CETTE FONCTION !!!!!!!!!!
   double getVar(std::vector<double>& vals);
   //  ... determination des fonctions de formes en coordonnees locales
   void getFF(std::vector<double>& ff);
   //  ... determination du gradient des fonctions de formes (exprime dans la base globale)
   void getGradFF(std::vector<xtensor::xVector<>>& gradff);

   //  ... volume de l element
   double getVolume();

   void setUvw(const xtensor::xPoint&);
   xtensor::xPoint getUvw() const;

  private:
   std::vector<xtensor::xPoint> p;
   const int dim;
   xtensor::xPoint uvw;
   xtensor::xTensor2<> jac;
   xtensor::xTensor2<> invjac;
   xtensor::xTensor2<> tinvjac;
   double detjac;
   static void getFFEdge(const double& u, std::vector<double>& ff);
   static void getFFTri(const double& u, const double& v, std::vector<double>& ff);
   static void getFFQuadSimplex(const double& u, const double& v, std::vector<double>& ff);
   static void getFFTet(const double& u, const double& v, const double& w, std::vector<double>& ff);
   static void getFFHexSimplex(const double& u, const double& v, const double& w, std::vector<double>& ff);
   AOMD::mEntity::mType type;
};

class xElementInfo
{
   typedef std::map<int, xtensor::xVector<>>::value_type value;

  public:
   //  ... definition des informations associees a un element
   xElementInfo(AOMD::mEntity* e);
   double getVolume() const { return volume; }
   xtensor::xVector<> getGradFFNode(int i) const { return grad_ff_node[i]; }

  private:
   xElementInfo();
   double volume;
   std::vector<xtensor::xVector<>> grad_ff_node;
   AOMD::mEntity* elem;
   mType type;
};

class xElementInfoManager
{
  public:
   typedef std::shared_ptr<xElementInfo> xElementInfoPtr;

  private:
   typedef std::unordered_map<AOMD::mEntity*, xElementInfoPtr, AOMD::EntityHashKey, AOMD::EntityEqualKey> info_t;
   typedef info_t::const_iterator iterator;
   typedef info_t::value_type value;

  public:
   //~xElementInfoManager();
   xElementInfoPtr getInfo(AOMD::mEntity* e);
   void clear();

  private:
   info_t infos;
};

}  // namespace xfem

#endif
