/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xcFrontSpaces_
#define _xcFrontSpaces_

#include "xApproxFunction.h"
#include "xSpace.h"

namespace xfem
{
class xIntegrationRule;
class xLevelSet;
class xSpaceRegular;
}  // namespace xfem

class xcApproxFunctionModalConstant : public xfem::xApproxFunction
{
  public:
   xcApproxFunctionModalConstant(const xfem::xLevelSet* _ls) : ls(_ls) {}
   std::string name() override { return "xcApproxFunctionModalConstant"; }
   void getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double&) const override;
   void getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;

  private:
   const xfem::xLevelSet* ls;
};

class ModalFourier
{
  public:
   ModalFourier(const xfem::xLevelSet& cos_s_modal_, const xfem::xLevelSet& sin_s_modal_, int _nbmodes)
       : ls_cos_s(cos_s_modal_),
         ls_sin_s(sin_s_modal_),
         nb_modes(_nbmodes),
         eintheta(nb_modes, std::complex<double>(0., 0.)),
         gx(nb_modes, std::complex<double>(0., 0.)),
         gy(nb_modes, std::complex<double>(0., 0.)),
         gz(nb_modes, std::complex<double>(0., 0.)),
         e_current_v(nullptr),
         e_current_g(nullptr)
   {
   }

   void getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ) const;
   void getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ) const;
   friend class xcApproxFunctionModalFourier;

  private:
   const xfem::xLevelSet& ls_cos_s;
   const xfem::xLevelSet& ls_sin_s;
   const int nb_modes;
   mutable std::vector<std::complex<double>> eintheta;
   mutable std::vector<std::complex<double>> gx;
   mutable std::vector<std::complex<double>> gy;
   mutable std::vector<std::complex<double>> gz;
   mutable const AOMD::mEntity* e_current_v;
   mutable xtensor::xPoint uvw_current_v;
   mutable const AOMD::mEntity* e_current_g;
   mutable xtensor::xPoint uvw_current_g;
   bool equal_uvw(const xtensor::xPoint& a, const xtensor::xPoint& b) const
   {
      return (a(0) == b(0) && a(1) == b(1) && a(2) == b(2));
   }
};

/*class xcApproxFunctionModalFourier : public xApproxFunction {
    public:
        enum cos_sin_t { COS, SIN };
        xcApproxFunctionModalFourier(const xLevelSet& cos_s_modal_,const xLevelSet& sin_s_modal_,   int i, cos_sin_t t);
        std::string name(){
            return "xcApproxFunctionModalFourier";
        }
        void  getVal(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, double&) const;
        void getGrad(const xfem::xGeomElem* geo_appro,const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const;

    private:
        const xLevelSet& cos_s_modal;
        const xLevelSet& sin_s_modal;

        int ith;
        cos_sin_t cos_or_sin;
        };*/

class xcApproxFunctionModalFourier : public xfem::xApproxFunction
{
  public:
   enum cos_sin_t
   {
      COS,
      SIN
   };
   xcApproxFunctionModalFourier(const ModalFourier& modalfourrier, int i, cos_sin_t t);
   std::string name() override { return "xcApproxFunctionModalFourier"; }
   void getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double&) const override;
   void getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;

  private:
   const ModalFourier& modalfourier;
   int ith;
   cos_sin_t cos_or_sin;
};

/// ApproxFunction Legendre(1d), defined between -1 and 1.
class xcApproxFunctionModalLegendre : public xfem::xApproxFunction
{
  public:
   xcApproxFunctionModalLegendre(const xfem::xLevelSet& s_modal_, int i);
   std::string name() override { return "xcApproxFunctionModalLegendre"; }
   void getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double&) const override;
   void getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;

  private:
   const xfem::xLevelSet& s_modal;
   int ith;
};

/// ApproxFunction 1d, set of hermite polynomes, order 3 in two parts
// (to get C^1 functions on the whole front AND partition of unity), defined between [-1,0] and [0,1].
// with uniform repartition
class xcApproxFunctionModalPolynomeHermiteUniform : public xfem::xApproxFunction
{
  public:
   xcApproxFunctionModalPolynomeHermiteUniform(const xfem::xLevelSet& s_modal_, int i, int nb_mod, double mode_width);
   xcApproxFunctionModalPolynomeHermiteUniform(const xfem::xLevelSet& _s_modal_cos, const xfem::xLevelSet& _s_modal_sin, int i,
                                               int nb_mod, double mode_width);
   std::string name() override { return "xcApproxFunctionModalPolynomeHermiteUniform"; }
   void getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double&) const override;
   void getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;

  private:
   const xfem::xLevelSet& s_modal;
   const xfem::xLevelSet& s_modal_sin;
   int ith;
   int total_number_of_modes;
   bool line_open;
   const double mode_width;
};

class xcSpaceModalConstant : public xfem::xSpaceRegular
{
  public:
   xcSpaceModalConstant(const std::string& a, TensorialType_t b, AOMD::mEntity* _keyedentity,
                        const xfem::xLevelSet* _ls = nullptr);
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  private:
   AOMD::mEntity* keyedentity;
   const xfem::xLevelSet* ls;
};

class xcSpaceModalFourier : public xfem::xSpaceRegular
{
  public:
   xcSpaceModalFourier(const std::string& a, TensorialType_t b, int nb_mod, const xfem::xLevelSet& _ls_cos_s,
                       const xfem::xLevelSet& _ls_sin_s);
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;
   int get_nb_modes() { return nb_modes; }

  private:
   const int nb_modes;

   const xfem::xLevelSet& ls_cos_s;
   const xfem::xLevelSet& ls_sin_s;
   ModalFourier modalfourier;
   std::vector<xfem::xValKey::ids_size_t> MODE_SIN_KEY;
   std::vector<xfem::xValKey::ids_size_t> MODE_COS_KEY;
};

class xcSpaceModalLegendre : public xfem::xSpaceRegular
{
  public:
   xcSpaceModalLegendre(const std::string& a, TensorialType_t b, int nb_mod, const xfem::xLevelSet& _curvi_coord_on_front_part);
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;
   int get_nb_modes() { return nb_modes; }

  private:
   const int nb_modes;
   const xfem::xLevelSet& curvi_coord_on_front_part;
   std::vector<xfem::xValKey::ids_size_t> MODE_LEGENDRE_KEY;
};

class xcSpaceModalPolynomeHermiteUniform : public xfem::xSpaceRegular
{
  public:
   xcSpaceModalPolynomeHermiteUniform(const std::string& a, TensorialType_t b, int nb_mod,
                                      const xfem::xLevelSet& _curvi_coord_on_front_part, double _mode_width);
   xcSpaceModalPolynomeHermiteUniform(const std::string& a, TensorialType_t b, int nb_mod,
                                      const xfem::xLevelSet& _curvi_coord_on_front_part,
                                      const xfem::xLevelSet& _curvi_coord_on_front_part_sin, double _mode_width);
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;
   int get_nb_modes() { return nb_modes; }

  private:
   const int nb_modes;
   const xfem::xLevelSet& curvi_coord_on_front_part;
   const xfem::xLevelSet& curvi_coord_on_front_part_sin;
   std::vector<xfem::xValKey::ids_size_t> MODE_POLYNOME_HERMITE_UNIFORM_KEY;
   bool line_open;
   const double mode_width;
};

#endif
