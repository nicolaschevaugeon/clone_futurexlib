/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef _xcFront_
#define _xcFront_

#include <functional>
#include <string>
// AOMD
#include "mExchangeData.h"
// xfem
#include "xField.h"
#include "xLevelSet.h"
#include "xLevelSetOperators.h"
#include "xRegion.h"
#include "xSpace.h"
#include "xSubMesh.h"

struct xParseData;

class xcSpaceModalPolynomeHermiteUniform;
namespace xfem
{
class xIntegrationRule;
class xSpaceRegular;
class xMesh;
}  // namespace xfem
namespace AOMD
{
class mVertex;
}

namespace xcrack
{
class xcInteractionIntegralsOnCrack;
class xcFrontPart;
class xcFront
{
  public:
   xcFront(const xfem::xMesh& front_mesh, const std::string& _frontname, const xfem::xMesh& mesh,
           std::function<double(AOMD::mVertex*)> front_distance, const xParseData& _parameters, bool _damage_type = false);
   ~xcFront();
   /// separate the crack fronts in simply connected parts and create the parts objects
   void createFronts();
   /// create domains to compute the configurational forces
   //  void createDomains();
   /// create space
   void createSpace();
   void reLoadSpace();
   xfem::xSpaceComposite& getSpace();
   xfem::xSpaceComposite& getSpace1D();
   /// export domains
   void exportDomainForIntegral() const;
   typedef std::list<xcFrontPart*>::iterator iterator;
   typedef std::list<xcFrontPart*>::const_iterator const_iterator;
   size_t size_front_parts() const;
   iterator begin_front_parts();
   iterator end_front_parts();
   const_iterator begin_front_parts() const;
   const_iterator end_front_parts() const;
   iterator begin() { return begin_front_parts(); }
   iterator end() { return end_front_parts(); }
   const_iterator begin() const { return begin_front_parts(); }
   const_iterator end() const { return end_front_parts(); }

   void clear_front_parts();
   const xfem::xMesh& getFrontMesh() const { return *front_mesh_global; }
   xfem::xSubMesh& createGlobalSubMesh();
   void extendParametrization();

  private:
   const xfem::xMesh& front_mesh;
   xfem::xMesh* front_mesh_global;
   const std::string front_name;
   const xfem::xMesh& mesh;
   std::function<double(AOMD::mVertex*)> front_distance;
   const xParseData& parameters;
   std::list<xcFrontPart*> fronts;
   xfem::xSpaceComposite space_modal_times_fr;
   xfem::xSpaceComposite front_modal_space_1d;
   bool damage_type;
};

/// This class define a part of the  front, simply connected or a loop.
/// it stors the space for the virtual speed as well as the domain for integration
class xcFrontPart
{
  public:
   friend class xcInteractionIntegralsOnCrack;
   enum frontType
   {
      Point,
      LineOpen,
      LineClose,
      None
   };
   xcFrontPart(const xfem::xMesh& front_mesh, const std::string& _frontpartname, const xfem::xMesh& mesh,
               std::function<double(AOMD::mVertex*)> front_distance, const xParseData& parameters, bool _damage_type = false);
   virtual ~xcFrontPart();
   const std::string& getFrontName() const { return front_part_name; }
   virtual void resetModes(double minwidth);
   virtual xfem::xApproxFunction* getModalApproxFunction(int imod, int dim);
   void exportDomainForIntegral() const;
   xfem::xSpaceXFEM* getSpace();
   xcSpaceModalPolynomeHermiteUniform* getSpace1D();

   const std::string subset_for_q_label;
   const std::string subset_for_int_label;
   const xfem::xRegion& getFrontRegion() const { return front_region; }
   double getFrontLenght() const { return totallenght; }
   virtual int get_nb_modes() { return nb_modes; }
   virtual double get_mode_width() { return mode_width; }
   virtual void extendParametrization();
   std::list<AOMD::mVertex*> local_vertices;

  private:
   void createDomainsAndRadialFunction();
   void deleteDomains();
   virtual void setParametrisation() = 0;
   virtual void deleteSpace();

  protected:
   const xfem::xMesh& front_mesh;
   const std::string front_part_name;
   const xfem::xRegion front_region;
   const xfem::xMesh& mesh;
   std::function<double(AOMD::mVertex*)> front_distance;
   frontType fronttype;
   const xParseData& parameters;
   bool damage_type;
   int nb_modes;
   xfem::xSpaceXFEM* modal_space;
   xcSpaceModalPolynomeHermiteUniform* modal_space_1d;
   xfem::xLevelSet fr;  // alpha field
   double totallenght, mode_width;
};

/// Specialisation for front that are just point (2D case)
class xcFrontPartPoint : public xcFrontPart
{
  public:
   xcFrontPartPoint(const xfem::xMesh& front_mesh, const std::string& _frontname, const xfem::xMesh& mesh,
                    std::function<double(AOMD::mVertex*)> _front_distance, const xParseData& parameters,
                    bool _damage_type = false);
   void createSpace();
   void setParametrisation() override;

  private:
};

/// specialisation for front part that are open (for example a line segment)
class xcFrontPartLineOpen : public xcFrontPart
{
  public:
   xcFrontPartLineOpen(const xfem::xMesh& front, const std::string& _frontname, const xfem::xMesh& mesh,
                       std::function<double(AOMD::mVertex*)> _front_distance, const xParseData& _parameters,
                       AOMD::mVertex* _start, AOMD::mVertex* _end, bool _damage_type = false);
   void setFrontParametrisation();
   void setParametrisation() override;
   void createSpace();
   void resetModes(double minwidth) override;
   xfem::xApproxFunction* getModalApproxFunction(int imod, int dim) override;
   const xfem::xLevelSet& getLss1d() const;
   void extendParametrization() override;
   xfem::xLevelSet lss3d;
   xfem::xLevelSet lss1d;

  private:
   AOMD::mVertex* start;
   AOMD::mVertex* end;
};

/// Specialisation for front that are closed (periodical, for example a circle)
class xcFrontPartLineClose : public xcFrontPart
{
  public:
   xcFrontPartLineClose(const xfem::xMesh& front, const std::string& _frontname, const xfem::xMesh& mesh,
                        std::function<double(AOMD::mVertex*)> _front_distance, const xParseData& _parameters,
                        AOMD::mVertex* _start, bool _damage_type = false);
   void setFrontParametrisation();
   void setParametrisation() override;
   void createSpace();
   void resetModes(double minwidth) override;
   xfem::xApproxFunction* getModalApproxFunction(int imod, int dim) override;
   const xfem::xLevelSet& getLss1dCos() const;
   const xfem::xLevelSet& getLss1dSin() const;
   void extendParametrization() override;
   xfem::xLevelSet lss3dCos;
   xfem::xLevelSet lss3dSin;
   xfem::xLevelSet lss1dCos;
   xfem::xLevelSet lss1dSin;

  private:
   AOMD::mVertex* start;
};

std::map<std::string, std::pair<AOMD::mVertex*, AOMD::mVertex*>> Separate1dMeshBranch(xfem::xMesh& mesh_crack_front,
                                                                                      const std::string& front_base_name);

void CreateSubsetforConfigForce(const xfem::xMesh& front_mesh, const std::string& front_part_name, const xfem::xMesh& mesh,
                                double rho_cylinder, int nb_layers_cylinder, int nb_layers_core, const std::string& subsetfor_q,
                                const std::string& subsetfor_domain_int);

/// fill lss_1D with a curvilinar parameter along the mesh crack_front.
/*! the Mesh is supposed to be a 1 d mesh. it can also be a 0d mesh (a node) then the parameter is set to 0 for this node
    each part of the front is parametrized separately. the function return a map of bool, indexed by the part number that store if
   the part is open (false ) or closed (true)
  */

std::list<AOMD::mVertex*> Parametrize1dMeshBranch(AOMD::mVertex* vstart, AOMD::mVertex* vend, xfem::xLevelSet& ls,
                                                  double& totallenght);
std::list<AOMD::mVertex*> Parametrize1dMeshLoop(AOMD::mVertex* vstart, xfem::xLevelSet& lscos, xfem::xLevelSet& lssin,
                                                double& totallenght);
void ExportSifs(const std::string& filename, const xfem::xField<>& j_modal, xfem::xRegion&);

class xcElementsAlongFrontCreator : public xfem::xSubMeshCreator, public AOMD::AOMD_DataExchanger
{
  public:
   xcElementsAlongFrontCreator(const xfem::xMesh& _front, const std::string& _name, std::function<double(AOMD::mVertex*)> d)
       : xfem::xSubMeshCreator(), front_mesh(_front), front_part_name(_name), front_distance(d)
   {
   }
   void create(const xfem::xMesh&, const std::string& name) override;
   // implementation of base class AOMD_DataExchanger member
   int tag() const override;
   void* AP_alloc_and_fill_buffer(AOMD::mEntity* e, AOMD::AOMD_SharedInfo& si, int tag) override;
   void receiveData(int pid, void* buf) override;

  private:
   const xfem::xMesh& front_mesh;
   xfem::xMesh* front_mesh_global;
   const std::string front_part_name;
   std::function<double(AOMD::mVertex*)> front_distance;
   std::string name;
   std::list<AOMD::mEntity*> toclean;
   mutable xinterface::aomd::xAttachedDataManagerAOMD<int> added;
};

class xcSetRadialFunction : public xfem::xLevelSetModifier
{
  public:
   xcSetRadialFunction();
   void visit(xfem::xLevelSet& f, xfem::xRegion target) override;
};

class xcApproxFunctionLevelSet : public xfem::xApproxFunction
{
  public:
   xcApproxFunctionLevelSet(const xfem::xLevelSet& l) : ls(l) {}
   void getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double&) const override;
   void getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   std::string name() override { return "xcApproxFunctionLevelSet"; }

  private:
   const xfem::xLevelSet& ls;
};

class xcApproxFunctionGradLevelSet : public xfem::xApproxFunction
{
  public:
   xcApproxFunctionGradLevelSet(const xfem::xLevelSet& l) : ls(l) {}
   std::string name() override { return "xcApproxFunctionGradLevelSet"; }
   void getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;

  private:
   const xfem::xLevelSet& ls;
};

class xcApproxFunctionNormedGradLevelSet : public xfem::xApproxFunction
{
  public:
   xcApproxFunctionNormedGradLevelSet(const xfem::xLevelSet& l) : ls(l) {}
   std::string name() override { return "xcApproxFunctionNormedGradLevelSet"; }
   void getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   void getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xTensor2<>&) const override;

  private:
   const xfem::xLevelSet& ls;
};
}  // namespace xcrack

#endif
