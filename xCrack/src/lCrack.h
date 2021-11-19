/*
     This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef ___LCRACK__H
#define ___LCRACK__H
#include <functional>

// aomd
#include "mEdge.h"
#include "mEntity.h"
// xtensor
#include "xTensor2.h"
#include "xVector.h"
// xfem
#include "xElement.h"
#include "xEntityToEntity.h"
#include "xExportGmsh.h"
#include "xField.h"
#include "xGeomElem.h"
#include "xLevelSet.h"
#include "xLevelSetOperators.h"
#include "xPointToDouble.h"
#include "xSubMesh.h"
#include "xVectorField.h"
#include "xcCrackBase.h"

// parallel
// CLEAN 16/11/2021 #include "xParallel.h"

// A crack is defined by the 0 isovalue of the 1st level set; and the negative value of the second level set
// The crack owns both level set representing it.
// To create a crack, a mesh is given as well as a two functors yielding required
// level set values at each node.

namespace xcrack
{
class xcEvalJintTerms;
/// struct to define parameter for xcrack instance
struct cParam
{
   cParam()
       : max_iter_reinit(100),
         max_iter_ortho(100),
         max_iter_extend(200),
         coeff_dt_virt(1),
         coeff_dt_phys(0.5),
         coeff_real_dt(10.0),
         tolerance(5e-5),
         tolfit(1e-3)
   {
   }
   int max_iter_reinit;
   int max_iter_ortho;
   int max_iter_extend;
   double coeff_dt_virt;
   double coeff_dt_phys;
   double coeff_real_dt;
   double tolerance;
   double tolfit;
};

//! class to describ a crack.
/*! lst refer to the level set describing the crack surface.
  lsn refer to the level set describing the surface where the crack front lies.
  lsn is supposed to be orthogonal to lst.
  lss refer to an curvilinear absice along the crack tipe.*/

class lCrack : public xcCrackBase
{
  public:
   /// lCrack constructor from the mesh. All the levelset function must be added later.
   /*! additional bool parameter fite :
     if set to true, will move the iso zero to vertex if close to it when calling post_initfunc.
     set to true by default.
   */
   lCrack(xfem::xMesh* m, bool fite = true, std::string _label = "");

   /// lCrack constructor from the mesh and 2 levelst function
   /*! additional bool parameter fite :
     if set to true, will move the iso zero to vertex if close to it when calling post_initfunc.
     set to true by default.
   */
   lCrack(xfem::xMesh* m, const xfem::xPointToDouble& lsfn, const xfem::xPointToDouble& lsft, bool fite = true,
          std::string _label = "");

   /// lCrack destructor
   ~lCrack() override;

   /// A filter that return true if the entities in entry has its support completly cut by the crack.
   /*!
     Primarly used for Heaviside enrichment
     Inherited from xcCrackBase
   !*/
   bool supportsCutthruByCrack(AOMD::mEntity* e) const override;
   xfem::xEntityFilter supportsCutthruByCrackFilter() const override;
   bool supportInCylinderAroundFront(double r, AOMD::mEntity* e) const;
   xfem::xEntityFilter supportInCylinderAroundFrontFilter(double r) const;

   /// This evaluator is meant to decide if an entity is above or below the crack surface.
   /*!
     Return one on top of the crack surface (lsn > 0), zero or -1 otherwise.
     It uses the usual xEval convention, geo_appro store the entity on which the approximation is defined, and geo_integ store the
   entity on which the information is evaluated. Inherited from xcCrackBase
   !*/
   int sideOf(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ) const override;

   /// Get distance from  the lip and angle;
   /*!
     From  emesh entity e and coordinate in the reference frame on this element, return r = distance to the crack lip and
      theta direction angle, 0 behing the crack plan
      Inherited from xcCrackBase
   !*/
   void getLocalCoords(AOMD::mEntity* e, xtensor::xPoint& uvw, double& r, double& theta) const override;
   std::vector<double> getLstValues(AOMD::mEntity* e) const override;
   std::vector<double> getLsnValues(AOMD::mEntity* e) const override;
   /// Transform the coordinate of a vector expressed in the crack coordinate system (t, n, t^n) to the global coordinate system
   /*!
     Entity e and uvw (coordinate of the point on which the vector is computed in the reference element of e)
     are used to compute the frame change. The "smoothed gradient" is use here
      Inherited from xcCrackBase
   !*/
   void localToGlobal(AOMD::mEntity* e, const xtensor::xPoint& uvw, xtensor::xVector<>& local) const override;

   /// Transform the coordinate of a tensor2 expressed in the crack coordinate system (t, n, t^n) to the global coordinate system
   /*!
     Entity e and uvw (coordinate of the point on which the Tensor  is computed in the reference element of e)
     are used to compute the frame change.
      Inherited from xcCrackBase
   !*/
   void localToGlobal(AOMD::mEntity* e, const xtensor::xPoint& uvw, xtensor::xTensor2<>& local) const override;

   void post_initfunc();
   /// set the two  xLevelSet describing the crack. first the ls normal then the levelset tangent. The gradient is computed for
   /// both (this involve communication)
   /*! additional parameter fite set to true by default. See constructor for definition.
    */
   void setLevelSets(xfem::xLevelSet& lsfn, xfem::xLevelSet& lsft, bool fite = true);
   double getFrontDistance(AOMD::mVertex* v) const override { return sqrt(pow(lst(v), 2) + pow(lsn(v), 2)); }
   // Name should be changed here ..
   /// return a pointer to xLevelSet tangent to the crack surface.
   xfem::xLevelSet const* getFieldt() const override { return &lst; }
   /// return a pointer to the xLevelSet where the crack lipe is.
   xfem::xLevelSet const* getFieldn() const override { return &lsn; }
   /// return a pointer to the xLevelSet giving a curvilinear abscisse along the cracklip of the closest point on the crack lip.
   xfem::xLevelSet const* getFields() const { return &lss; }
   xfem::xLevelSet const* getField_cos_s() const { return &lss_cos; }
   xfem::xLevelSet const* getField_sin_s() const { return &lss_sin; }
   xfem::xVectorField const* getV3D() const { return &v3D; }
   xfem::xLevelSet const* getVn() const { return &v_oplane3D; }
   xfem::xLevelSet const* getVt() const { return &v_iplane3D; }

   /// Return true if e is touched by the crack
   /*! renvoie true si le support de l'entite e est touche par la fissure.
     elements : la fissure touche l'element (par un noeud ou une face ou encore une arete), ou le traverse totalement ou non
     noeuds : la fissure touche le support du noeud ( idem def precedente)
     arretes ou face : idem
   !*/
   bool supports_touched_by_crack(AOMD::mEntity* e) const;

   /// renvoie les entites dont le support est touch�par la (les) pointe(s)s de fissure
   /*!
     meme definitions que pour touched_by_crack
   !*/
   bool supports_touched_by_front(AOMD::mEntity* e) const;

   /// renvoie les entites dont le support est effectivement coupe par la fissure (totalement ou pas)
   /*! ne renvoie par exemple pas tous les elements qui sont touched_by_crack.
     elements : revoie ceux qui sont traverses par la fissure
     noeuds : renvoie ceux dont le support est traverse par la fissure
   !*/
   bool supports_cut_by_crack(AOMD::mEntity* e) const;

   bool side_of_crack_surface(AOMD::mEntity* e) const;

   /// renvoie le support concerne par l'integration d'ordre plus eleve pour les fonctions de forme asymptotiques
   bool supports_integrate_front(AOMD::mEntity* e) const;

   /// return the dimension of the mesh support of the crack
   int dim() const;

   //  xMesh const *getMesh(void) const {return mesh;}
   //  xMesh const *getMeshCrackFront(void) const {return mesh_crack_front;}
   xfem::xMesh* getMesh() override { return mesh; }
   const xfem::xMesh& getMesh() const override { return *mesh; }
   xfem::xMesh* getMeshCrackFront() override { return mesh_crack_front; }
   const xfem::xMesh& getMeshCrackFront() const override { return *mesh_crack_front; };

   xfem::xMesh* getMeshCrackSurface() { return mesh_crack_surface; }

   /// From one element on which lsn and lst must be defined, return the crack tip bases and their gradient.
   void getLocalCurv(AOMD::mEntity* e, const xtensor::xPoint& uvw, xtensor::xVector<>& e1, xtensor::xVector<>& e2,
                     xtensor::xVector<>& e3, xtensor::xTensor2<>& curv1, xtensor::xTensor2<>& curv2,
                     xtensor::xTensor2<>& curv3) const override;

   void getLocalSmoothOrthoAxis(AOMD::mEntity* e, const xtensor::xPoint& uvw, xtensor::xPoint& local_coords,
                                xtensor::xVector<>& eo1, xtensor::xVector<>& eo2, xtensor::xVector<>& eo3) const;

   void getCrackAxis(AOMD::mEntity* e, xtensor::xTensor2<>& base) const override;
   void getCrackOrthoAxis(AOMD::mEntity* e, xtensor::xTensor2<>& base) const override;

   void getLocalSmoothOrthoAxis(AOMD::mEntity* e, xtensor::xVector<>& eo1, xtensor::xVector<>& eo2,
                                xtensor::xVector<>& eo3) const;

   void setLocalInfos(AOMD::mEntity* e, const xtensor::xPoint& uvw) const;

   void extractSIFs(const xfem::xField<>& disp_l, std::string filename, bool autop = true, double ratiop = 2.0, double hp = 1.0,
                    double r1p = 1.0, double r2p = 1.0);  // output file with the sifs

   void getJint3D(const xfem::xField<>& disp_l, xfem::xValueManagerDist<double>* DoubleManager,
                  xfem::xIntegrationRule& integrator, double width, double radius, std::string filename, int layers = 1) const;

   void getMatF3D(xfem::xField<>& disp_l, xfem::xValueManagerDist<double>* DoubleManager, xfem::xIntegrationRule& integrator,
                  double width, double radius, std::string filename, int layers = 1);

   void getJint3DParks(xfem::xField<>& disp_l, xfem::xIntegrationRule& integrator_vol, xfem::xIntegrationRule& integrator_bnd,
                       const double& rho, const int& nb_internodes, std::string filename,
                       const std::function<double(const xtensor::xPoint&)>& p_to_s);

   /*void getJint3DModalParks(xField<>& disp_l,
                            xIntegrationRule& integrator_vol,
                            xIntegrationRule& integrator_bnd,
                            const double& rho_cylinder,
                             int nb_layers_cylinder,
                            int nb_layers_core,
                            int nb_modes,
                            std::string filename);

   void getJint3DModalParks(xField<>& disp_l,
                            xIntegrationRule& integrator_vol,
                            xIntegrationRule& integrator_bnd,
                            const double& rho_cylinder,
                            int nb_layers_cylinder,
                            int nb_layers_core,
                            int nb_modes,
                            std::string filename,
                            const xcEvalJintTerms& terms_evaluator
                            );

   void getJint3DModalParksV2(xField<>& disp_l,
                              xIntegrationRule& integrator_vol,
                              xIntegrationRule& integrator_bnd,
                              const double& rho_cylinder,
                              int nb_layers_cylinder,
                              int nb_layers_core,
                              int nb_modes,
                              std::string filenamebase
                              );

   void getJint3DModalParksV2(
                              xIntegrationRule& integrator_vol,
                              xIntegrationRule& integrator_bnd,
                              const double& rho_cylinder,
                              int nb_layers_cylinder,
                              int nb_layers_core,
                              int nb_modes,
                              std::string filenamebase,
                              const xcEvalJintTerms& terms_evaluator
                              );
   */

   void printSifs(std::string filename, const xfem::xValueManagerDist<double>& double_manager,
                  const std::function<double(const xtensor::xPoint&)>& p_to_s) const;
   void printSifsModal(const std::string& filename, const xfem::xField<>& J_modal, AOMD::mEntity*);

   void extractMeshCrackFrontDofsForJint3DParks(
       xfem::xMesh* mesh_crack_front_dofs,
       std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d_s_loc,
       std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e1d_dof_length,
       std::unordered_map<AOMD::mEntity*, AOMD::mEntity*, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d_dof_e0d,
       std::unordered_map<AOMD::mEntity*, AOMD::mEntity*, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d_e1d_dof,
       int nb_internodes) const;

   void set3DElementsToFrontRelation(
       const xfem::xRegion& domain_for_integral,
       const std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d_s_loc,
       const std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e1d_dof_length,

       const std::unordered_map<AOMD::mEntity*, AOMD::mEntity*, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d_e1d_dof,
       std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d3d_s_loc,
       std::unordered_map<AOMD::mEntity*, AOMD::mEntity*, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d3d_e1d_dof) const;

   double getClosestSegmentFromTo(AOMD::mVertex* v, xfem::xMesh* mesh1d, AOMD::mEntity*& edge_closest,
                                  double& percent_position) const;

   void setSlocLevelSet(
       const xfem::xRegion& domain_for_integral, xfem::xLevelSet& sloc_level_set,
       const std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d3d_s_loc,
       const std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e1d_dof_length,

       const std::unordered_map<AOMD::mEntity*, AOMD::mEntity*, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d3d_e1d_dof,

       std::unordered_map<AOMD::mEntity*, AOMD::mEntity*, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d3d_e1d_dof_other,
       bool& sloc_ok) const;

   void setAllSlocLevelSets(
       const xfem::xRegion& domain_for_integral,
       std::unordered_map<AOMD::mEntity*, xfem::xLevelSet*, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d_dof_sloc_level_set,
       const std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d3d_s_loc,
       const std::unordered_map<AOMD::mEntity*, double, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e1d_dof_length,
       const std::unordered_map<AOMD::mEntity*, AOMD::mEntity*, AOMD::EntityHashKey, AOMD::EntityEqualKey>& e0d3d_e1d_dof) const;

   void cutElementsOnTransitions(const xfem::xLevelSet& sloc_level_set) const;
   void setQDirectionParks(xfem::xVectorField& q_dir) const;
   static void orderEdges(AOMD::mEdge*& ed1, AOMD::mEdge*& ed2);
   void setLssOnFront(double mins = -1., double maxs = 1.);

  public:
   mutable double Jscalar;
   /// return a reference to the current time step
   int& tstep();
   /// modify lst et lsn so that they fit to vertices
   void fittovertices();
   double propagate();
   double lspropagate();
   double lsnextend();
   double reinit();
   double reinit2D();
   double reortho();
   double modifyVoplane();
   void createNarrowBand();
   void createSupportInfo();
   void cutMeshAndSlice();
   void fromV1DToVioplane();
   void fromV1DToVioplaneOnly();

   void updateForNarrowBand();
   void exportV3D(string s);
   void exportVoplane3D(string s);
   void exportViplane3D(string s);
   void exportVoplane1D(string s);
   void exportViplane1D(string s);
   void exportlsn(string s);
   void exportlst(string s);
   void exportmeshcracksurface(string s);
   void exportmeshcracktip(string s);
   void exportlstcracksurface(string s);

   void setV1Ddebug(std::function<xtensor::xVector<>(xtensor::xPoint)> func);
   void extendVelocity();   // goes from v1D to v_iplane3D and v_oplane3D and V3D
   void extendVelocity2();  // goes from v1D to v_iplane3D and v_oplane3D

   void setV1Ddebug(xtensor::xVector<> Vec);
   void setVioplane3Ddebug(double iplane, double oplane);
   void getAsymptoticFields(const xfem::xGeomElem* geo_appro, int mode, xtensor::xTensor2<>& stress_aux,
                            xtensor::xTensor2<>& grad_disp_aux) const override;
   bool checkIfClosed()  // return true if the crack front is closed, false other wise
   {
      if (!lss_created) setLssOnFront();
      return is_closed;
   }

  private:
   xfem::xPointToDouble lsft;
   xfem::xPointToDouble lsfn;

  protected:
   bool fit;
   xfem::xMesh* mesh;    // all the domain
   xfem::xLevelSet lsn;  // defining the support for the crack surface )
   xfem::xLevelSet lst;  // defining the front of the crack
   bool lss_created;
   xfem::xLevelSet lss;      // curvilinear abs
   xfem::xLevelSet lss_cos;  // cosine of the curvilinear abs (for closed crack)
   xfem::xLevelSet lss_sin;  // sine of the curvilinear abs (for closed crack)

   xfem::xLevelSet lss_1D;      // curvilinear abs defined on crack front
   xfem::xLevelSet lss_1D_cos;  // cosine of the curvilinear abs (for closed crack)
   xfem::xLevelSet lss_1D_sin;  // sine of the curvilinear abs (for closed crack)
   bool is_closed;              // flag set to 1 by setlss  if the crack front is closed

   xfem::xMesh* mesh_crack_surface;    // mesh of the crack surface (surface mesh (3d) or linear mesh (in 2d) )
   xfem::xMesh* mesh_crack_front;      // mesh of the crack tip ( linear (3d) or just nodes (in 2D) )
   xfem::xLevelSet lst_crack_surface;  // interpolated normal level set on the mesh of the surface of the crack
   xfem::xRegion narrow_band;

   mutable xfem::xVectorField v1D;
   xfem::xVectorField v3D;

   xfem::xLevelSet v_iplane1D;
   xfem::xLevelSet v_oplane1D;
   xfem::xLevelSet v_x1D;
   xfem::xLevelSet v_y1D;
   xfem::xLevelSet v_z1D;
   xfem::xLevelSet v_iplane3D;
   xfem::xLevelSet v_oplane3D;
   xfem::xLevelSet v_x3D;
   xfem::xLevelSet v_y3D;
   xfem::xLevelSet v_z3D;

   // ne sert que pour les export
   int timestep;

   AOMD::mEntity* otherEdgeOnVertex(AOMD::mVertex* v, AOMD::mEntity* e) const;
   AOMD::mVertex* otherVertexOnEdge(AOMD::mEntity* e, AOMD::mVertex* v) const;
   double lengthOfEdge(AOMD::mEntity* e) const;
   void set_Subset_Labels();

  public:
   // local infos stored at the current point

   // local coordinates of the point with global coordinates p
   // seul the two first coordinates make sense
   mutable xtensor::xPoint local;

   // these are the normal to the iso-level sets
   // they are one form, i.e dx^1, dx^2, dx^3
   // for instance n1 is the derivative of the first level-set (the front one)
   // with respect to the global coordinates
   // THESE AXIS ARE OBTAINED FROM the smoothed gradient averaged at the nodes
   mutable xtensor::xVector<> n1, n2, n3;
   mutable double n1_norm, n2_norm, n3_norm;
   // normalized vectors
   mutable xtensor::xVector<> no1, no2, no3;

   // the matrix below are composed by placing the basis vectors
   // n1, n2, n3 or no1, no2, no3
   // either as row or column
   mutable xtensor::xTensor2<> n123_row, n123_col;
   mutable xtensor::xTensor2<> no123_row, no123_col;

   // the metric obtained by ni dot nj is called metricn
   // the metric obtained by ti dot tj is called metricn
   // these two matrices are inverse of each other
   mutable xtensor::xTensor2<> metricn;

   // curv1 is the second derivatives of the first (front) level-sets
   // with respect to the physical coordinates
   // curv2 is for the second derivatives and
   mutable xtensor::xTensor2<> curv1, curv2, curv3;
   // derivatives of the normalized vectors
   mutable xtensor::xTensor2<> curvo1, curvo2, curvo3;

   cParam Param;

  public:
   //  std::string label;
   std::string touched_by_crack_label;
   std::string touched_by_front_label;
   std::string cut_by_crack_label;
   std::string cut_by_crack_to_delete_label;
   std::string cutthru_by_crack_label;
   std::string integrate_front_label;
   std::string narrow_band_label;
   std::string all_minus_touched_by_front_label;

  private:
   xexport::xExportGmshAscii pexport;
};

class xcSupportsTouchedByCrackCreator : public xfem::xSubMeshCreator
{
  public:
   xcSupportsTouchedByCrackCreator(const xfem::xLevelSet& c);
   void create(const xfem::xMesh&, const string& name) override;

  private:
   const xfem::xLevelSet& crack2D;
};

class xcSupportsTouchedByFrontCreator : public xfem::xSubMeshCreator
{
  public:
   xcSupportsTouchedByFrontCreator(xfem::xMesh* m1D);
   void create(const xfem::xMesh&, const string& name) override;

  private:
   xfem::xMesh* mesh_front1D;
};

class xcSupportsIntegrateCreator : public xfem::xSubMeshCreator
{
  public:
   xcSupportsIntegrateCreator(const string& name_tbf_);
   void create(const xfem::xMesh&, const string& name) override;

  private:
   string name_tbf;
};

class xcSupportsCutByCrackCreator : public xfem::xSubMeshCreator
{
  public:
   xcSupportsCutByCrackCreator(const xfem::xLevelSet& c);
   void create(const xfem::xMesh&, const string& name) override;

  private:
   const xfem::xLevelSet& crack2D;
};

class xcSupportsCutByCrackToDeleteCreator : public xfem::xSubMeshCreator
{
  public:
   xcSupportsCutByCrackToDeleteCreator(const xfem::xLevelSet& c);
   void create(const xfem::xMesh&, const string& name) override;

  private:
   const xfem::xLevelSet& crack2D;
};

class xcSupportsCutthruByCrackCreator : public xfem::xSubMeshCreator
{
  public:
   xcSupportsCutthruByCrackCreator(const xfem::xLevelSet& c, const string& ncut, const string& ndel);
   void create(const xfem::xMesh&, const string& name) override;

  private:
   const xfem::xLevelSet& crack2D;
   string name_cut;
   string name_del;
};

struct tipfilterlcrack
{
   tipfilterlcrack(const lCrack& _crack, double _r) : crack(_crack), r(_r){};
   bool operator()(AOMD::mEntity* e) const { return crack.supportInCylinderAroundFront(r, e); }

  private:
   const lCrack& crack;
   double r;
};

// Union

class xcStrictNegativeSideCreator : public xfem::xSubMeshCreator
{
  public:
   xcStrictNegativeSideCreator(const xfem::xLevelSet& c);
   void create(const xfem::xMesh&, const string& name) override;

  private:
   const xfem::xLevelSet& field;
};

class xcNegativeSideCreator : public xfem::xSubMeshCreator
{
  public:
   xcNegativeSideCreator(const xfem::xLevelSet& c);
   void create(const xfem::xMesh&, const string& name) override;

  private:
   const xfem::xLevelSet& field;
};

/// sub_mesh creator : using element in the sub "supports_touched_by_front" create a sub_mesh "name".
/*!
  For each element in supports_touched_by_front, find the closest node and add all element sharing this
  node in sub_mesh "name".
  !*/
class xcElementsAlongCrackFrontCreator : public xfem::xSubMeshCreator
{
  public:
   xcElementsAlongCrackFrontCreator(const lCrack& _c);
   void create(const xfem::xMesh&, const string& name) override;

  private:
   const lCrack& c;
};

class xcModifyOutOfPlaneCrackVelocity : public xfem::xLevelSetModifier
{
  public:
   xcModifyOutOfPlaneCrackVelocity(const xfem::xLevelSet& front3D, const xfem::xLevelSet& plane3D,
                                   const xfem::xLevelSet& v_iplane3D, string s, double d)
       : front(front3D), plane(plane3D), v_i(v_iplane3D), dt(d), supports_cuthru_by_crack(s)
   {
   }
   void visit(xfem::xLevelSet& v_oplane3D, xfem::xRegion target) override;

  private:
   const xfem::xLevelSet &front, plane, v_i;
   double dt;
   string supports_cuthru_by_crack;
};

class xcValKeyExtendAndSetRefe
{
  public:
   xcValKeyExtendAndSetRefe(const std::string& e_, int ref_) : extension(e_), ref(ref_) {}
   void operator()(xfem::xValKey& key)
   {
      std::string geom_new = xfem::xKeyInfo::getGeomName(key.getGeom()) + extension;
      key.setGeom(xfem::xKeyInfo::getGeomId(geom_new));
      key.setRefe(ref);
   }

  private:
   std::string extension;
   int ref;
};

}  // namespace xcrack

#endif
