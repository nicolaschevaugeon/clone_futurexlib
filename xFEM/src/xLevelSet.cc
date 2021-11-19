/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xLevelSet.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "mAttachableDataContainer.h"
#include "mEdge.h"
#include "mEntity.h"
#include "mFace.h"
#include "mIterator.h"
#include "mTet.h"
#include "mVertex.h"
#include "xDebug.h"
#include "xElement.h"
#include "xIntegrationRule.h"
#include "xLevelSetOperators.h"
#include "xMesh.h"
#include "xRegion.h"

// eXlibris_tools
#include "xDataExchanger.h"

// to be removed when fully parallel
#include "workInProgress.h"

using AOMD::mEntity;
using AOMD::mVertex;

namespace xfem
{
xLevelSet::~xLevelSet()
{
   if (keys_container) delete keys_container;
}

xLevelSet::xLevelSet()
    : key_up_to_date(false), grad_up_to_date(false), curv_up_to_date(false), true_curv_up_to_date(false), keys_container(nullptr)
{
}

xLevelSet::xLevelSet(const xRegion& s, const double& val)
    : key_up_to_date(false),
      grad_up_to_date(false),
      curv_up_to_date(false),
      true_curv_up_to_date(false),
      support(s),
      keys_container(nullptr)
{
   for (xIter it = support.begin(0); it != support.end(0); ++it)
   {
      mVertex* v = (mVertex*)*it;
      ls[v] = val;
   }
   keys_container =
       new xtool::xKeyContainerSendAndRecv<xLevelSetKeyManager::information_key_t>(support.getPartitionManager().getComm());
}

xLevelSet::xLevelSet(const xRegion& s, const xPointToDouble& func)
    : key_up_to_date(false),
      grad_up_to_date(false),
      curv_up_to_date(false),
      true_curv_up_to_date(false),
      support(s),
      keys_container(nullptr)
{
   load(func);
   keys_container =
       new xtool::xKeyContainerSendAndRecv<xLevelSetKeyManager::information_key_t>(support.getPartitionManager().getComm());
}

xLevelSet::xLevelSet(const xLevelSet& in) : xLevelSet()
{
   this->operator=(in);
   return;
}

void xLevelSet::load(const xPointToDouble& func)
{
   for (xIter it = support.begin(0); it != support.end(0); ++it)
   {
      AOMD::mVertex* v = (AOMD::mVertex*)*it;
      ls[v] = func(v->point());
   }
}

void xLevelSet::loadOnElement(const xPointToDouble& func, AOMD::mEntity* e)
{
   for (int i = 0; i < e->size(0); ++i)
   {
      AOMD::mVertex* v = (AOMD::mVertex*)e->get(0, i);
      ls[v] = func(v->point());
   }
}

xLevelSet& xLevelSet::operator=(const xLevelSet& other)
{
   if (keys_container)
   {
      delete keys_container;
      keys_container = nullptr;
   }
   support = other.support;

   for (xIter it = support.begin(0); it != support.end(0); ++it)
   {
      mVertex* v = (mVertex*)*it;
      ls[v] = other(v);
   }
   keys_container =
       new xtool::xKeyContainerSendAndRecv<xLevelSetKeyManager::information_key_t>(support.getPartitionManager().getComm());

   key_up_to_date = false;
   grad_up_to_date = false;
   curv_up_to_date = false;
   true_curv_up_to_date = false;
   return *this;
}

void xLevelSet::complement()
{
   for (auto& keyval : ls) keyval.second = -keyval.second;
}

double& xLevelSet::operator()(mEntity* v)
{
   // cout << "arg" ;
   key_up_to_date = false;
   grad_up_to_date = false;
   curv_up_to_date = false;
   true_curv_up_to_date = false;
   return ls[v];
}

bool xLevelSet::isDefinedAt(mEntity* v) const { return ls.count(v); }

bool xLevelSet::isDefinedOnElement(mEntity* e) const
{
   for (int i = 0; i < e->size(0); ++i)
   {
      mEntity* v = e->get(0, i);
      if (!isDefinedAt(v)) return false;
   }
   return true;
}

int xLevelSet::side_of(const xGeomElem* geo_appro, const xGeomElem* geo_integ) const
{
   const bool debug = false;
   unsigned int side_tag;
   side_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
   int orientation = (geo_integ->getEntity())->getAttachedInt(side_tag);

   if (orientation == 0)
   {
      if (debug) cout << "in sideof : in domain\n";
      mEntity* e = geo_appro->getEntity();
      // mEntity* e_sub;
      xElement elem(e);
      elem.setUvw(geo_appro->getUVW());
      std::vector<double> vals = getVals(e);
      double f = elem.getInterpoSca(vals);
      if (fabs(f) < 1e-3)
      {
         if (geo_integ != nullptr)
         {
            // e_sub = geo_integ->getEntity();
            elem.xyz2uvw(geo_integ->getCDGxyz());
         }
         else
         {
            // e_sub = geo_appro->getEntity();
            elem.setUvw(geo_appro->getCDGuvw());
         }

         std::vector<double> vals = getVals(e);
         f = elem.getInterpoSca(vals);
      }

      if (f >= 0)
         return (1);
      else
         return (-1);
   }
   else
   {
      if (debug) cout << "in sideof : on slice\n";
      if (debug)
      {
         cout << "Entity in sideof :";
         geo_integ->getEntity()->print();
         cout << endl;
      }
      if (debug) cout << "Orientation in sideof\n" << orientation << endl;
      if (orientation >= 0)
         return (1);
      else
         return (-1);
   }
}

double& xLevelSet::at(const AOMD::mVertex& v)
{
   key_up_to_date = false;
   grad_up_to_date = false;
   curv_up_to_date = false;
   true_curv_up_to_date = false;
   return ls.at(const_cast<mVertex*>(&v));
}

const double& xLevelSet::at(const AOMD::mVertex& v) const { return ls.at(const_cast<mVertex*>(&v)); }

double xLevelSet::operator()(const mEntity* pe) const
{
   const mVertex* pv = dynamic_cast<const mVertex*>(pe);
   if (pv) return at(*pv);
   throw;
}

double xLevelSet::operator()(const mVertex& v) const { return at(v); }

const xRegion& xLevelSet::getSupport() const { return support; }

void xLevelSet::setSupport(const xRegion& m, const double& val)
{
   support = m;
   clear();
   for (mEntity* pe : support.range(0))
   {
      mVertex* pv = static_cast<mVertex*>(pe);
      ls[pv] = val;
   }
}

void xLevelSet::reduceSupport(const xRegion& m)
{
   support = m;
   key_up_to_date = false;
   if (keys_container) delete keys_container;
   keys_container =
       new xtool::xKeyContainerSendAndRecv<xLevelSetKeyManager::information_key_t>(support.getPartitionManager().getComm());
}

void xLevelSet::clear()
{
   ls.clear();
   grad.clear();
   grad_up_to_date = false;
   curv.clear();
   curv_up_to_date = false;
   true_curv.clear();
   true_curv_up_to_date = false;
   key_up_to_date = false;
   if (keys_container) delete keys_container;
   keys_container =
       new xtool::xKeyContainerSendAndRecv<xLevelSetKeyManager::information_key_t>(support.getPartitionManager().getComm());
}

void xLevelSet::accept(xLevelSetModifier& visitor) { visitor.visit(*this, getSupport()); }
void xLevelSet::accept(xLevelSetModifier& visitor, const xRegion& target) { visitor.visit(*this, target); }

void xLevelSet::accept(xLevelSetInspector& visitor) const { visitor.visit(*this, getSupport()); }
void xLevelSet::accept(xLevelSetInspector& visitor, const xRegion& target) const { visitor.visit(*this, target); }

void xLevelSet::accept(xLevelSetCreator& visitor, xLevelSet& f) const { visitor.visit(*this, f); }

std::vector<double> xLevelSet::getVals(const mEntity* e) const { return getVals(*e); }

std::vector<double> xLevelSet::getVals(const mEntity& e) const
{
   const bool debug = xdebug_flag;
   std::vector<double> vals;
   vals.reserve(e.size(0));
   if (debug)
   {
      cout << " e in getvals " << endl;
      e.print();
   }
   if (e.getLevel() > 0)
   {
      for (int i = 0; i < e.size(0); i++)
      {
         mVertex& v = static_cast<mVertex&>(*e.get(0, i));
         lstype::const_iterator itend = ls.end();
         lstype::const_iterator it = ls.find(&v);
         if (it != itend)
            vals.push_back(it->second);
         else
         {
            std::cerr << "can't find vertex " << &v << "in level set " << this;
            throw 1;
         }
      }
   }
   else
   {
      lstype::const_iterator itend = ls.end();
      lstype::const_iterator it = ls.find(const_cast<AOMD::mEntity*>(&e));
      if (it != itend)
         vals.push_back(it->second);
      else
      {
         std::cerr << "can't find vertex " << &e << "in level set " << this;
         throw 1;
      }
   }
   if (debug)
   {
      cout << " ls values at the nodes " << endl;
      cout << " number of nodes " << vals.size() << endl;
      std::copy(vals.begin(), vals.end(), std::ostream_iterator<double>(cout, " "));
      cout << endl;
   }
   return vals;
}

// get values of the level set for the nodes given in argument
// v is a vector of pointer to mEntity entity
std::vector<double> xLevelSet::getVals(std::vector<mEntity*>& v) const
{
   std::vector<double> vals;
   std::vector<mEntity*>::iterator itv = v.begin();
   std::vector<mEntity*>::const_iterator itve = v.end();
   lstype::const_iterator itl;
   lstype::const_iterator itle = ls.end();

   for (; itv != itve; ++itv)
   {
      itl = ls.find(*itv);
      if (itl != itle)
         vals.push_back(itl->second);
      else
      {
         mVertex* itvv = (mVertex*)*itv;
         std::cerr << "can't find vertex " << itvv << "in level set " << this;
         throw 1;
      }
   }
   return vals;
}

xtensor::xVector<> xLevelSet::getGrad(mEntity* e) const
{
   //  ... cette fonction calcule le gradient du xLevelSet sur tous les noeuds du maillage
   const bool debug = xdebug_flag;
   if (debug) cout << " e in getgrad is " << endl;

   if (debug) e->print();
   if (debug) cout << " bool grad_up_to_date is " << grad_up_to_date << endl;
   if (!grad_up_to_date) compute_grad();
   auto it = grad.find(e);
   if (it == grad.end()) throw 1;
   xtensor::xVector<> out(grad.find(e)->second);
   return out;
}

double xLevelSet::getVal(mEntity* e, const xtensor::xPoint& uvw) const
{
   //  ... cette fonction calcule EN UN POINT D UN ELEMENT la valeur du champ du xLevelSet
   //      a partir des fonctions de forme et de la valeur aux noeuds de l element
   //      de ce meme xLevelSet
   // const bool debug = xdebug_flag;  const int nnodes = e->size(0);
   std::vector<double> vals = getVals(e);

   switch (e->getType())
   {
      case mEntity::VERTEX:
         return vals[0];
      case mEntity::EDGE:
      {
         const double u = uvw(0);
         return vals[0] * (0.5 * (1. - u)) + vals[1] * (0.5 * (1. + u));
      }
      case mEntity::TRI:
      {
         const double u = uvw(0);
         const double v = uvw(1);
         return vals[0] * (1. - u - v) + vals[1] * u + vals[2] * v;
      }
      case mEntity::TET:
      {
         const double u = uvw(0);
         const double v = uvw(1);
         const double w = uvw(2);
         return vals[0] * (1. - u - v - w) + vals[1] * u + vals[2] * v + vals[3] * w;
      }
      case mEntity::QUAD:
      {
         xElement elem(e);
         elem.setUvw(uvw);
         double out = elem.getInterpoSca(vals);
         return out;
      }
      case mEntity::HEX:
      {
         xElement elem(e);
         elem.setUvw(uvw);
         double out = elem.getInterpoSca(vals);
         return out;
      }
      default:
         throw;
   }
}

bool xLevelSet::getVal(mEntity* e, const xtensor::xPoint& uvw, double& val) const
{
   try
   {
      val = getVal(e, uvw);
      return true;
   }
   catch (...)
   {
      val = 0.;
      return false;
   }
}

xtensor::xVector<> xLevelSet::getGrad(mEntity* e, const xtensor::xPoint& uvw) const
{
   //  ... cette fonction calcule EN UN POINT D UN ELEMENT la valeur du gradient du xLevelSet
   //      a partir des fonctions de forme et de la valeur aux noeuds de l element du gradient
   //      de ce meme xLevelSet
   if (!grad_up_to_date) compute_grad();
   // xElement elem(e);
   //  ... determination des coordonnees locales du noeud considere
   // elem.xyz2uvw(p);
   // elem.setUvw(uvw);
   //  get des gradients aux noeuds
   const int nnodes = e->size(0);
   std::vector<xtensor::xVector<>> grad_v(nnodes);
   for (int j = 0; j < nnodes; j++)
   {
      grad_v[j] = getGrad(e->get(0, j));
   }
   switch (e->getType())
   {
      case mEntity::VERTEX:
         return grad_v[0];
      case mEntity::EDGE:
      {
         const double u = uvw(0);
         return grad_v[0] * (0.5 * (1. - u)) + grad_v[1] * (0.5 * (1. + u));
      }
      case mEntity::TRI:
      {
         const double u = uvw(0);
         const double v = uvw(1);
         return grad_v[0] * (1. - u - v) + grad_v[1] * u + grad_v[2] * v;
      }
      case mEntity::TET:
      {
         const double u = uvw(0);
         const double v = uvw(1);
         const double w = uvw(2);
         return grad_v[0] * (1. - u - v - w) + grad_v[1] * u + grad_v[2] * v + grad_v[3] * w;
      }
      case mEntity::QUAD:
      {
         xElement elem(e);
         //  ... determination des coordonnees locales du noeud considere
         elem.setUvw(uvw);
         xtensor::xVector<> out(elem.getInterpoVec(grad_v));
         return out;
      }
      case mEntity::HEX:
      {
         xElement elem(e);
         //  ... determination des coordonnees locales du noeud considere
         elem.setUvw(uvw);
         xtensor::xVector<> out(elem.getInterpoVec(grad_v));
         return out;
      }
      default:
         throw;
   }
}

bool xLevelSet::getGrad(mEntity* e, const xtensor::xPoint& uvw, xtensor::xVector<>& vals) const
{
   try
   {
      vals = getGrad(e, uvw);
      return true;
   }
   catch (...)
   {
      return false;
   }
}

xtensor::xTensor2<> xLevelSet::getCurv(mEntity* e) const
{
   if (!curv_up_to_date) compute_curv();
   return (curv.find(e)->second);
}

xtensor::xTensor2<> xLevelSet::getTrueCurv(mEntity* e) const
{
   if (!true_curv_up_to_date) compute_true_curv();
   return (true_curv.find(e)->second);
}

void xLevelSet::create_key() const
{
   xLevelSetKeyManager km(support.getPartitionManager());
   keys_container->accumulateKeysAllGather(ls.begin(), ls.end(), km);
   key_up_to_date = true;
}
//  rentre des valeurs dans grad_region, grad_vertex
void xLevelSet::compute_grad() const
{
   const bool debug = xdebug_flag;
   if (debug) cout << " support in compute_grad is of size dim=3 " << support.size(3) << endl;
   for (xIter it = support.begin(); it != support.end(); ++it)
   {
      mEntity* e = *it;
      if (debug) cout << " compute_grad for element " << endl;
      if (debug) e->print();
      std::vector<double> ls_v = getVals(e);
      if (debug) cout << " nodal level set values " << endl;
      if (debug) std::copy(ls_v.begin(), ls_v.end(), std::ostream_iterator<double>(cout, " "));
      xElement elem(e);
      grad[e] = elem.getGradInterpoSca(ls_v);
   }

   // tempory container
   containers_t<int> nbconnect;

   for (xIter it = support.begin(0); it != support.end(0); ++it)
   {
      mEntity* v = *it;
      grad[v] = xtensor::xVector<>(0., 0., 0.);
      nbconnect[v] = 0;
   }

   for (xIter it = support.begin(); it != support.end(); ++it)
   {
      // for(mIteratorPtr it(support.getTopIter()); !it->done(); it->next() ) {
      mEntity* e = *it;
      for (int j = 0; j < e->size(0); ++j)
      {
         mEntity* v = e->get(0, j);
         grad[v] += grad[e];
         nbconnect[v] += 1;
      }
   }

   // create key if not already created
   if (!key_up_to_date) create_key();

   // sum contribution from other proc
   xLevelSetGradInfoManager im(grad, nbconnect);
   exchangeInformation(*keys_container, im);

   for (xIter it = support.begin(0); it != support.end(0); ++it)
   {
      // for(mIteratorPtr it(support.getIter(0)); !it->done(); it->next() ) {
      mEntity* v = *it;
      grad[v] /= (double)nbconnect[v];
   }
   grad_up_to_date = true;
   return;
}

//  rentre des valeurs dans grad
void xLevelSet::compute_curv() const
{
   // let grad be the field gradient
   // curv is the matrix [grad,1  grad,2 grad,3]
   if (!grad_up_to_date) compute_grad();
   xtensor::xTensor2<> ten;
   for (xIter it = support.begin(); it != support.end(); ++it)
   {
      // for(mIteratorPtr it(support.getTopIter()); !it->done(); it->next() ) {
      //  ... determination de l element courant
      mEntity* e = *it;
      std::vector<xtensor::xVector<>> grad_v;
      for (int j = 0; j < e->size(0); j++)
      {
         mVertex* v = (mVertex*)e->get(0, j);
         grad_v.push_back(getGrad(v));
      }
      xElement elem(e);
      elem.getGradInterpoVec(grad_v, ten);
      // curv must be symmetric it would be good to check
      curv[e] = ten;
   }
   curv_up_to_date = true;
   return;
}

void xLevelSet::compute_true_curv() const
{
   if (!grad_up_to_date) compute_grad();
   xtensor::xTensor2<> ten;
   for (xIter it = support.begin(); it != support.end(); ++it)
   {
      mEntity* e = *it;
      std::vector<xtensor::xVector<>> grad_v;
      for (int j = 0; j < e->size(0); j++)
      {
         mVertex* v = (mVertex*)e->get(0, j);
         grad_v.push_back(getGrad(v).norm());
      }
      xElement elem(e);
      elem.getGradInterpoVec(grad_v, ten);
      true_curv[e] = ten;
   }
   true_curv_up_to_date = true;
   return;
}

void xLevelSet::restrictTo(xRegion new_support, double init_val)
{
   for (xIter it = new_support.begin(0); it != new_support.end(0); ++it)
   {
      // for(mIteratorPtr it(new_support.getIter(0)); !it->done(); it->next() ) {
      mEntity* v = *it;
      if (!isDefinedAt(v)) ls[v] = init_val;
   }
   for (xIter it = support.begin(0); it != support.end(0); ++it)
   {
      // for(mIteratorPtr it(support.getIter(0)); !it->done(); it->next() ) {
      mEntity* v = *it;
      if (!new_support.find(v))
      {
         ls.erase(v);
      }
   }
   support = new_support;
}

void xLevelSet::printDebug() const
{
   for (AOMD::mEntity* pe : support.range(0))
   {
      AOMD::mVertex* pv = static_cast<AOMD::mVertex*>(pe);
      cout << "At Point :" << pv->point() << " Value ls is :" << ls.find(pv)->second << endl;
   }
}

void xLevelSet::exportMatlab(ostream& fout, const std::string& field_name, int level) const
{
   fout << field_name << "= [ ";
   for (mEntity* e : support.range(level)) fout << ls.find(e)->second << " ";
   fout << "];\n";
   return;
}

void xLevelSet::takeTraceOn(const xMesh& target_mesh, const xMesh::datamanager_t<AOMD::mEntity*>& was_created_by,
                            const xMesh::datamanager_t<double>& r_on_edge, xLevelSet& trace) const
{
   const bool debug = xdebug_flag;
   if (debug) cout << "takeTraceOn : entree\n\n";
   trace.setSupport(&target_mesh);
   for (mEntity* pv : target_mesh.range(0))
   {
      mVertex& v = static_cast<mVertex&>(*pv);
      mEntity& e = *was_created_by.at(v);
      switch (e.getLevel())
      {
         case 0:
         {
            trace.at(v) = (*this)(static_cast<mVertex&>(e));
            break;
         }
         case 1:
         {
            // we need to interpolate
            const double r = r_on_edge.at(v);
            const double lsv0 = (*this)(static_cast<mVertex*>(e.get(0, 0)));
            const double lsv1 = (*this)(static_cast<mVertex*>(e.get(0, 1)));
            trace.at(v) = lsv0 + r * (lsv1 - lsv0);
            break;
         }
         default:
            throw;
      }
   }
}

/// enables to evaluate x_field on the nodes of mesh_new
/// where mesh_new is necessarily a submesh
///
/// was_duplicated_from() was_created_by r_on_edge
void xLevelSet::interpolateTo(const xMesh& mesh_new, const xMesh::datamanager_t<AOMD::mEntity*>& was_duplicated_from,
                              const xMesh::datamanager_t<AOMD::mEntity*>& was_created_by, xLevelSet& lsnew) const
{
   const bool debug = xdebug_flag;
   if (debug) cout << "interpolateTo : ENTREE\n\n";
   lsnew.setSupport(&mesh_new);
   for (mEntity* pe : mesh_new.range(0))
   {
      mVertex& v = static_cast<mVertex&>(*pe);
      mVertex& vd = static_cast<mVertex&>(*was_duplicated_from.at(v));
      if (debug)
      {
         cout << "noeud courant du submesh " << v.point() << "****\n";
         cout << " id is " << v.getId() << endl;
         cout << "noeud duplicateur " << vd.point() << "****\n";
      }
      mEntity* const* ppec = was_created_by.getData(vd);
      if (!ppec)
      {
         lsnew.at(v) = (*this)(static_cast<mVertex&>(vd));
         if (debug) cout << "valeur de lsnew en " << v.point() << "est " << lsnew(v) << "****\n";
      }
      else
      {
         mEntity& ec = **ppec;
         if (debug)
         {
            cout << "level de e = " << ec.getLevel() << "*****\n";
            ec.print();
         }
         if (ec.getLevel() == 0)
            lsnew.at(v) = (*this)(static_cast<mVertex&>(ec));
         else  // e is then an edge
         {
            xElement el(&ec);
            el.xyz2uvw(vd.point());
            const double lsv0 = (*this)(static_cast<mVertex&>(*ec.get(0, 0)));
            const double lsv1 = (*this)(static_cast<mVertex&>(*ec.get(0, 1)));
            lsnew.at(v) = el.getInterpoSca(std::vector<double>{lsv0, lsv1});
            if (debug) cout << "valeur de lsv0 en  " << v.point() << " est " << lsv0 << "****\n";
            if (debug) cout << "valeur de lsv1 en  " << v.point() << " est " << lsv1 << "****\n";
            if (debug) cout << "valeur de lsnew en " << v.point() << " est " << lsnew(v) << "****\n";
         }
      }
   }
   if (debug) cout << "interpolateTo : SORTIE\n\n";
   return;
}

void xLevelSet::forceGradComputation()
{
   if (!grad_up_to_date) compute_grad();
}
//---------------------------------------
xLevelSet::xLevelSetKeyManager::xLevelSetKeyManager(const partmanAOMD_t& _partman) : partman(_partman) {}

//---------------------------------------
xtool::xConstPartitionObject<AOMD::mEntity> xLevelSet::xLevelSetKeyManager::getConstPartitionObject(const data_t& o) const
{
   AOMD::mEntity& e = *o.first;
   return partman.getConstPartitionObject(e);
}

//---------------------------------------
auto xLevelSet::xLevelSetKeyManager::localObjectKey(const data_t& o) const -> information_key_t { return o.first; }

//---------------------------------------
auto xLevelSet::xLevelSetKeyManager::remoteObjectKey(const xtool::xRemoteObject<AOMD::mEntity>& ro, const data_t& lo) const
    -> information_key_t
{
   return (AOMD::mEntity*)(ro.getObjectAddress());
}
//---------------------------------------
xLevelSet::xLevelSetGradInfoManager::xLevelSetGradInfoManager(xLevelSet::containers_t<xtensor::xVector<>>& grad_,
                                                              xLevelSet::containers_t<int>& nbconnect_)
    : s(3 * sizeof(MPI_DOUBLE) + sizeof(MPI_INT)), grad(grad_), nbconnect(nbconnect_)
{
}
//---------------------------------------
void xLevelSet::xLevelSetGradInfoManager::getInfo(information_key_t key, xtool::xMpiInputBuffer& buff, int sendto)
{
   auto k = const_cast<AOMD::mEntity*>(key);
   buff.pack(grad[k].data(), 3, MPI_DOUBLE);
   buff.pack(&nbconnect[k], 1, MPI_INT);
}
//---------------------------------------
void xLevelSet::xLevelSetGradInfoManager::setInfo(information_key_t key, const xtool::xMpiOutputBuffer& buff, int receivedfrom)
{
   xtensor::xVector<> g;
   int nbc;
   auto k = const_cast<AOMD::mEntity*>(key);
   buff.unPack(g.data(), 3, MPI_DOUBLE);
   grad[k] += g;
   buff.unPack(&nbc, 1, MPI_INT);
   nbconnect[k] += nbc;
}
//---------------------------------------
size_t xLevelSet::xLevelSetGradInfoManager::getApproxDataSize() { return s; }

//---------------------------------------
}  // namespace xfem
