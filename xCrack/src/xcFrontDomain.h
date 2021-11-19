/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xcFrontPartDomain_
#define _xcFrontPartDomain_

#include <functional>
#include <string>

// AOMD
#include "ParUtil.h"
#include "mExchangeData.h"
// xfem
#include "xLevelSetOperators.h"
#include "xMesh.h"
#include "xSubMesh.h"
// xcrack
#include "xcFrontSpaces.h"
// xExport
#include "xAlgorithm.h"
#include "xExportGmsh.h"

class xcFrontPartBase;
class xParseData;

namespace xfem
{
class xIntegrationRule;
class xSpaceRegular;
}  // namespace xfem

enum processType
{
   CRACKFRONT,
   DAMAGEFRONT
};

/// This class define a domain, i.e. a xsubmesh and a levelset for parametrisation, for a given frontpart
/*!
    notes : parameters for the domain creation and the verobosity of the function are taken from the one used to construct    the
 object of the class xcFrontPart. (true on objet of the clas xParseData.) The parameters should include int : "front_verbosity". 0
 no message. 1 more message to follow what's happening. and, depending on the value of processType:

    Damage front : The front is 1d, the domain is 2d.
        double :  "front_domain_rho_geo"
              "front_domain_fr_distance"



    front_domain_xx are used to describe a domain around a front. The details are a bit subtle,
    but basically, _rho_geo refer to a distance in the damaged part,
    and _fr_distance refer to a distance in the non damaged part.

    CrackFront : The front represent the crack tip two dimension (a point)
                 or a curve (the crack front) in tree dimension.
         The domain is disk around the point in 2d, a cylinder around the point in 3d.
         The cylinder might have a cylindrical hole (the core").
         double :  "front_domain_rho_geo"
         int :     "front_domain_nb_layers_core"
                 int :     "front_domain_nb_layers_cylinder"

         The cylinder is constructed by first adding _nb_layers_core element layers around the elements crossed by the front. This
 give a "minimum" size to the cylinder. Then more layer are added, until no element in the added layer have at least one node at a
 distance less                  than _rho_geo. Then some _nb_layers are removed from the cylinder, starting at the core

 !*/
class xcFrontDomain
{
  public:
   friend class xcFrontDomainManager;
   xcFrontDomain(xcFrontPartBase *_part, int _physical_process);
   virtual ~xcFrontDomain();
   const std::string &getFrontName() const { return front_part_name; }
   const std::string subset_for_q_label;
   const std::string subset_for_int_label;
   const xcFrontPartBase *getFrontPart() const { return part; }
   const xParseData &getParameters() const { return parameters; }
   const xfem::xLevelSet &getFr() const { return fr; }
   virtual int get_nb_modes() { return nb_modes; }
   virtual void set_nb_modes(int n) { nb_modes = n; }
   virtual xcApproxFunctionModalPolynomeHermiteUniform *getNewApproxFunction(int imod, int dim) = 0;

  private:
   void createDomainAndRadialFunction();
   virtual void extendParametrization() = 0;

  protected:
   int nb_modes;
   std::function<double(AOMD::mVertex *)> &front_distance;
   const std::string front_part_name;
   const xcFrontPartBase *part;
   int physical_process;
   xfem::xLevelSet fr;  // alpha field
   const xParseData &parameters;
   const xfem::xMesh &mesh, &front_mesh;
};

///---------------------------------------------------------- Specialization for point
class xcFrontDomainPoint : public xcFrontDomain
{
  public:
   friend class xcFrontDomainManager;
   xcFrontDomainPoint(xcFrontPartBase *_part, int _physical_process);
   xcApproxFunctionModalPolynomeHermiteUniform *getNewApproxFunction(int imod, int dim) override { return nullptr; }

  private:
   void extendParametrization() override;
};

///---------------------------------------------------------- Specialization for line open
class xcFrontDomainLineOpen : public xcFrontDomain
{
  public:
   friend class xcFrontDomainManager;
   xcFrontDomainLineOpen(xcFrontPartBase *_part, int _physical_process);
   const xfem::xLevelSet &getLss3d() const { return lss3d; }
   xcApproxFunctionModalPolynomeHermiteUniform *getNewApproxFunction(int imod, int dim) override;

  private:
   void extendParametrization() override;
   xfem::xLevelSet lss3d;
};

///---------------------------------------------------------- Specialization for line closed
class xcFrontDomainLineClosed : public xcFrontDomain
{
  public:
   friend class xcFrontDomainManager;
   xcFrontDomainLineClosed(xcFrontPartBase *_part, int _physical_process);
   const xfem::xLevelSet &getLss3dCos() const { return lss3dCos; }
   const xfem::xLevelSet &getLss3dSin() const { return lss3dSin; }
   xcApproxFunctionModalPolynomeHermiteUniform *getNewApproxFunction(int imod, int dim) override;

  private:
   void extendParametrization() override;
   xfem::xLevelSet lss3dCos, lss3dSin;
};

class xcSetRadialFunction : public xfem::xLevelSetModifier
{
  public:
   xcSetRadialFunction();
   void visit(xfem::xLevelSet &f, xfem::xRegion target) override;
};

///----------------------------------------------------------
class xcElementsAlongFrontCreator : public xfem::xSubMeshCreator, public AOMD::AOMD_DataExchanger
{
  public:
   xcElementsAlongFrontCreator(const xfem::xMesh &_front, const std::string &_name, std::function<double(AOMD::mVertex *)> d)
       : xfem::xSubMeshCreator(), front_mesh(_front), front_part_name(_name), front_distance(d)
   {
   }
   void create(const xfem::xMesh &, const std::string &name) override;
   // implementation of base class AOMD_DataExchanger member
   int tag() const override;
   void *AP_alloc_and_fill_buffer(AOMD::mEntity *e, AOMD::AOMD_SharedInfo &si, int tag) override;
   void receiveData(int pid, void *buf) override;

  private:
   const xfem::xMesh &front_mesh;
   xfem::xMesh *front_mesh_global;
   const std::string front_part_name;
   std::function<double(AOMD::mVertex *)> front_distance;
   std::string name;
   std::list<AOMD::mEntity *> toclean;
   mutable xinterface::aomd::xAttachedDataManagerAOMD<int> added;
};

///----------------------------------------------------------

class xcElementsAtLcDistanceFrom1DFront : public xfem::xSubMeshCreator
{
  public:
   xcElementsAtLcDistanceFrom1DFront(const xfem::xMesh &_front, const std::string &_name, double _ann_distance,
                                     double _fr_distance, std::function<double(AOMD::mVertex *)> dist_fct)
       : xSubMeshCreator(),
         front_mesh(_front),
         front_part_name(_name),
         ann_distance(_ann_distance),
         fr_distance(_fr_distance),
         front_distance(dist_fct)
   {
      std::cout << ann_distance << " " << fr_distance << std::endl;
   }
   void create(const xfem::xMesh &, const std::string &name) override;
   // implementation of base class AOMD_DataExchanger member
   /*		int tag() const;
       void * AP_alloc_and_fill_buffer (AOMD::mEntity *e, AOMD::AOMD_SharedInfo &si, int tag);
       void receiveData (int pid, void *buf);*/
  private:
   const xfem::xMesh &front_mesh;
   xfem::xMesh *front_mesh_global;
   const std::string front_part_name;
   double ann_distance, fr_distance;
   std::string name;
   std::function<double(AOMD::mVertex *)> front_distance;
   //		std::list<AOMD::mEntity *> toclean;
};

#endif
