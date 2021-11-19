/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#include <complex>
#include <iostream>

// xcrack
#include "xcFront.h"
#include "xcFrontSpaces.h"
// xfem
#include "xGatherMesh.h"
#include "xMesh.h"
#include "xParseData.h"
// xExport
#include "xExportAlgorithm.h"
#include "xExportGmsh.h"

// AOMD
#include "ParUtil.h"
#include "mEdge.h"

#ifdef PARALLEL
#include "autopack.h"
#endif

#include <cstdlib>

using std::cout;
using std::endl;
using namespace AOMD;
using namespace xfem;

namespace xcrack
{
xcFront::xcFront(const xMesh &_front_mesh, const string &_frontname, const xMesh &_mesh, std::function<double(mVertex *)> f,
                 const xParseData &_parameters, bool _damage_type)
    : front_mesh(_front_mesh),
      front_name(_frontname),
      mesh(_mesh),
      front_distance(f),
      parameters(_parameters),
      damage_type(_damage_type)
{
   createFronts();
   createSpace();
}

xcFront::~xcFront()
{
   for (xcFrontPart *fp : fronts)
   {
      ((*front_mesh_global)).deleteSubMesh(fp->getFrontName());
      delete (fp);
   }
   fronts.clear();
   delete front_mesh_global;
   front_mesh_global = nullptr;
}

// a mesh of the whole front integration region is required
// to solve multiple fronts at once
xSubMesh &xcFront::createGlobalSubMesh()
{
   xSubMesh &globalSubMesh = mesh.createSubMesh("globalsubmesh");
   // add other frontparts submeshes
   for (xcFrontPart *currentPart : fronts)
   {
      const string subset_name = currentPart->subset_for_int_label;
      xSubMesh *mysubmesh = (const_cast<xSubMesh *>(&mesh.getSubMesh(subset_name)));
      for (mEntity *e : mysubmesh->range(2))
      {
         globalSubMesh.add(e);
         // seules les faces sont ajoutées par defaut... manque les noeuds...
         int nb_vertices_on_face = e->size(0);
         for (int inode = 0; inode < nb_vertices_on_face; inode++)
         {
            globalSubMesh.add(e->get(0, inode));
         }
      }
   }
   return globalSubMesh;
}

// extends the parametrisation of one front part to the other frontparts domains for integration
// i.e. extends the levelsets supports
// NB: the values of the modes in extended domains doesn't matter since multiplied by the value of fr = zero.
void xcFront::extendParametrization()
{
   for (xcFrontPart *currentPart : fronts) currentPart->extendParametrization();
}

xcFront::iterator xcFront::begin_front_parts() { return fronts.begin(); }
xcFront::iterator xcFront::end_front_parts() { return fronts.end(); }
xcFront::const_iterator xcFront::begin_front_parts() const { return fronts.begin(); }
xcFront::const_iterator xcFront::end_front_parts() const { return fronts.end(); }

size_t xcFront::size_front_parts() const { return fronts.size(); }

void xcFront::createFronts()
{
   const bool debug = false;
   int verbose = parameters.getInt("sifs_verbosity");
   xMesh &front_mesh_mod = (const_cast<xMesh &>(front_mesh));
   AOMD::mMesh &mmesh = front_mesh_mod.getMesh();
   // front_mesh_mod.modifyAllState();
   mmesh.modifyState(3, 2, true);
   mmesh.modifyState(3, 1, true);
   mmesh.modifyState(3, 0, true);
   mmesh.modifyState(2, 1, true);
   mmesh.modifyState(2, 0, true);
   mmesh.modifyState(1, 0, true);
   mmesh.modifyState(0, 1, true);
   mmesh.modifyState(0, 2, true);
   mmesh.modifyState(0, 3, true);
   mmesh.modifyState(1, 2, true);
   mmesh.modifyState(1, 3, true);
   mmesh.modifyState(2, 3, true);

   front_mesh_global = new xMesh;
   if (verbose) cout << "broadcasting_front_mesh" << endl;
   xEntityCopyCallBackAttachEntity<xMesh::datamanager_t> callback(xMesh::get_const_was_created_by(), xMesh::get_was_created_by());
   BroadCastMesh(&front_mesh_mod, front_mesh_global, &callback);

   if (verbose) cout << "separting front_mesh in connected parts" << endl;
   std::map<std::string, std::pair<mVertex *, mVertex *>> frontParts = Separate1dMeshBranch(*front_mesh_global, front_name);
   if (debug) cout << " xcInteractionIntegralsOnCrack::createFront : nb front " << frontParts.size() << endl;
   if (verbose) cout << "creating xcFronts" << endl;
   for (const auto &s_vv : frontParts)
   {
      xcFrontPart::frontType fronttype;
      int dim = front_mesh_global->dim();
      if (dim == 0)
         fronttype = xcFrontPart::Point;
      else if (dim == 1)
      {
         if (s_vv.second.first == s_vv.second.second)
            fronttype = xcFrontPart::LineClose;
         else
            fronttype = xcFrontPart::LineOpen;
      }
      else
         abort();
      std::string frontsubname = s_vv.first;
      mVertex *start = s_vv.second.first;
      mVertex *end = s_vv.second.second;
      switch (fronttype)
      {
         case xcFrontPart::Point:
         {
            fronts.push_back(
                new xcFrontPartPoint(*front_mesh_global, frontsubname, mesh, front_distance, parameters, damage_type));
            break;
         }
         case xcFrontPart::LineOpen:
         {
            fronts.push_back(new xcFrontPartLineOpen(*front_mesh_global, frontsubname, mesh, front_distance, parameters, start,
                                                     end, damage_type));
            break;
         }
         case xcFrontPart::LineClose:
         {
            fronts.push_back(
                new xcFrontPartLineClose(*front_mesh_global, frontsubname, mesh, front_distance, parameters, start, damage_type));
            break;
         }
         case xcFrontPart::None:
         {
            cout << "warning can't create  front part with type xcFrontPart::None" << endl;
            break;
         }
      }
   }
   if (verbose) cout << "done xcFronts" << endl;
};

void xcFront::createSpace()
{
   for (xcFrontPart *pf : fronts)
   {
      xfem::xSpaceXFEM *sp = (pf->getSpace());
      space_modal_times_fr.insert(*sp);
   }
}

// PEB: modes (and spaces) should not be always created in the constructor
// since several parameters might still be needed: variable number of modes,
// modes depending on a specific coordinate, etc .
// resetModes and reLoadSpace allows one to redefine the spaces

// DEPRECATED !!!!!!!!!!!!!!!!!!!!!!!!!!

void xcFront::reLoadSpace()
{
   // modal spaces have changed, need to update space_modal_times_fr
   space_modal_times_fr.clear();
   front_modal_space_1d.clear();
   for (xcFrontPart *pf : fronts)
   {
      space_modal_times_fr.insert(*(pf->getSpace()));
      front_modal_space_1d.insert(*(pf->getSpace1D()));
   }
}

xSpaceComposite &xcFront::getSpace() { return space_modal_times_fr; }

xSpaceComposite &xcFront::getSpace1D() { return front_modal_space_1d; }

xcFrontPart::xcFrontPart(const xMesh &_front_mesh, const std::string &_frontpartname, const xMesh &_mesh,
                         std::function<double(mVertex *)> _front_distance, const xParseData &_parameters, bool _damage_type)
    : subset_for_q_label(_frontpartname + "_subset_for_q"),
      subset_for_int_label(_frontpartname + "_subset_for_int"),
      front_mesh(_front_mesh),
      front_part_name(_frontpartname),
      front_region(&_front_mesh, _frontpartname),
      mesh(_mesh),
      front_distance(_front_distance),
      fronttype(None),
      parameters(_parameters),
      damage_type(_damage_type),
      modal_space(nullptr),
      modal_space_1d(nullptr)
{
   int verbose = parameters.getInt("sifs_verbosity");
   if (verbose) cout << "begin_create_domain" << endl;
   createDomainsAndRadialFunction();
   if (verbose) cout << "end_create_domain" << endl;
};

xSpaceXFEM *xcFrontPart::getSpace() { return modal_space; }

xcSpaceModalPolynomeHermiteUniform *xcFrontPart::getSpace1D() { return modal_space_1d; }

void xcFrontPart::deleteDomains()
{
   const_cast<xMesh &>(mesh).deleteSubMesh(subset_for_q_label);
   const_cast<xMesh &>(mesh).deleteSubMesh(subset_for_int_label);
}

void xcFrontPart::resetModes(double minwidth)
{
   std::cerr << " resetModes is not implemented yet for this front type ! " << endl;
}
xApproxFunction *xcFrontPart::getModalApproxFunction(int imod, int dim)
{
   std::cerr << " getModalApproxFunction is not implemented yet for this front type ! " << endl;
   return nullptr;
}

void xcFrontPart::exportDomainForIntegral() const
{
   int verbose = parameters.getInt("sifs_verbosity");
   if (verbose) cout << "exportDomainForIntegral" << endl;
   mesh.getSubMesh(subset_for_int_label).exportGmsh(front_part_name + "_domain_for_integral.msh", 1);
}

xcFrontPart::~xcFrontPart() { deleteDomains(); }

void xcFrontPart::deleteSpace()
{
   if (modal_space)
   {
      delete modal_space;
      modal_space = nullptr;
   }
   if (modal_space_1d)
   {
      delete modal_space_1d;
      modal_space_1d = nullptr;
   }
};

void xcFrontPart::extendParametrization()
{
   std::cerr << " reset Parametrization not implemented yet for this front type ! " << endl;
}

xcFrontPartPoint::xcFrontPartPoint(const xMesh &_front_mesh, const std::string &_frontname, const xMesh &_mesh,
                                   std::function<double(mVertex *)> _front_distance, const xParseData &_parameters,
                                   bool _damage_type)
    : xcFrontPart(_front_mesh, _frontname, _mesh, _front_distance, _parameters, _damage_type)
{
   fronttype = Point;
   xValKeyExtend key_modifier("");
   xcApproxFunctionLevelSet bla(fr);
   // xcApproxFunction *pbla = &bla;
   modal_space = new xSpaceXFEM(xcSpaceModalConstant("J_modal", xSpace::SCALAR, (*front_region.begin(0))), bla, key_modifier);
}

void xcFrontPartPoint::setParametrisation(){};

xcFrontPartLineOpen::xcFrontPartLineOpen(const xMesh &_front_mesh, const std::string &_frontname, const xMesh &_mesh,
                                         std::function<double(mVertex *)> _front_distance, const xParseData &_parameters,
                                         mVertex *_start, mVertex *_end, bool _damage_type)
    : xcFrontPart(_front_mesh, _frontname, _mesh, _front_distance, _parameters, _damage_type), start(_start), end(_end)
{
   fronttype = LineOpen;
   int verbose = parameters.getInt("sifs_verbosity");
   if (verbose) cout << "start parametrisation" << endl;
   setParametrisation();
   if (verbose) cout << "parametrisation done" << endl;
   int nb_modes = parameters.getInt("sifs_nb_modes");
   xValKeyExtend key_modifier("_mult_fr_for_" + getFrontName());
   modal_space = new xSpaceXFEM(xcSpaceModalLegendre("J_modal", xSpace::SCALAR, nb_modes, lss3d), xcApproxFunctionLevelSet(fr),
                                key_modifier);
};

void xcFrontPartLineOpen::setFrontParametrisation()
{
   lss1d.setSupport(front_region);
   local_vertices = Parametrize1dMeshBranch(start, end, lss1d, totallenght);
}

void xcFrontPartLineOpen::extendParametrization()
{
   lss1d.setSupport(xRegion(&front_mesh));
   local_vertices = Parametrize1dMeshBranch(start, end, lss1d, totallenght);
   lss3d.setSupport(xRegion(&mesh, "globalsubmesh"));
   xVelocityInitialization lss_i(lss1d, xRegion(&mesh, subset_for_q_label));
   lss3d.accept(lss_i, xRegion(&mesh, subset_for_q_label));
   fr.setSupport(xRegion(&mesh, "globalsubmesh"));
   xcSetRadialFunction set_fr;
   fr.accept(set_fr, xRegion(&mesh, subset_for_q_label));

   int verbose = parameters.getInt("sifs_verbosity");
   if (verbose)
   {
      xexport::xExportGmshAscii pexport;
      xIntegrationRuleBasic intr;
      std::stringstream outfilename;
      if (ParUtil::Instance()->size() > 1)
         outfilename << front_part_name << "_fr_extend_" << ParUtil::Instance()->size() << "_" << ParUtil::Instance()->rank() + 1;
      else
         outfilename << front_part_name << "_fr_extend";
      Export(fr, pexport, outfilename.str());
   }
}

void xcFrontPartLineOpen::setParametrisation()
{
   int verbose = parameters.getInt("sifs_verbosity");
   if (verbose) cout << "set 1d Parametrisation " << endl;
   setFrontParametrisation();
   if (verbose) cout << "done 1d Parametrisation " << endl;

   lss3d.setSupport(xRegion(&mesh, subset_for_q_label));
   xVelocityInitialization lss_i(lss1d, xRegion(&mesh, subset_for_q_label));
   lss3d.accept(lss_i, xRegion(&mesh, subset_for_q_label));
   xexport::xExportGmshAscii pexport;
   xIntegrationRuleBasic intr;
   std::stringstream outfilename("lss3d");
   if (ParUtil::Instance()->size() > 1)
      outfilename << "_" << ParUtil::Instance()->size() << "_" << ParUtil::Instance()->rank() + 1 << ".pos";

   Export(lss3d, pexport, outfilename.str());
   if (verbose) cout << "export lss3d" << endl;
}

void xcFrontPartLineOpen::resetModes(double minwidth)
{
   // computes the number of modes to ensure a minimum mode width of minwidth
   // for the "spline" mode type
   nb_modes = (int)(ceil(2 * getFrontLenght() / minwidth + 1));
   mode_width = 2 * getFrontLenght() / (nb_modes - 1);

   // delete possible previous modal spaces
   if (modal_space) delete modal_space;
   if (modal_space_1d) delete modal_space_1d;

   // create new modal space
   xValKeyExtend key_modifier("_mult_fr_for_" + getFrontName());
   modal_space = new xSpaceXFEM(
       xcSpaceModalPolynomeHermiteUniform("Modal_Damage_Int", xSpace::SCALAR, nb_modes, lss3d, mode_width / totallenght * 2),
       xcApproxFunctionLevelSet(fr), key_modifier);
   modal_space_1d = new xcSpaceModalPolynomeHermiteUniform("Modal_Damage_One_Dimensional", xSpace::SCALAR, nb_modes, lss1d,
                                                           mode_width / totallenght * 2);
}

xApproxFunction *xcFrontPartLineOpen::getModalApproxFunction(int imod, int dim)
{
   xApproxFunction *modalpoly;
   switch (dim)
   {
      case 1:
         modalpoly = new xcApproxFunctionModalPolynomeHermiteUniform(lss1d, imod, nb_modes, mode_width / getFrontLenght() * 2.);
         break;
      case 2:
         modalpoly = new xcApproxFunctionModalPolynomeHermiteUniform(lss3d, imod, nb_modes, mode_width / getFrontLenght() * 2.);
         break;
      default:
         cerr << "wrong dimension in xcFrontPartLineOpen::getModalApproxFunction" << endl;
         throw;
   }
   return modalpoly;
}

const xLevelSet &xcFrontPartLineOpen::getLss1d() const { return lss1d; }

xcFrontPartLineClose::xcFrontPartLineClose(const xMesh &_front_mesh, const std::string &_frontname, const xMesh &_mesh,
                                           std::function<double(mVertex *)> _front_distance, const xParseData &_parameters,
                                           mVertex *_start, bool _damage_type)
    : xcFrontPart(_front_mesh, _frontname, _mesh, _front_distance, _parameters, _damage_type), start(_start)
{
   fronttype = LineClose;
   setParametrisation();
   int nb_modes = parameters.getInt("sifs_nb_modes");
   xValKeyExtend key_modifier("");
   modal_space = new xSpaceXFEM(xcSpaceModalFourier("J_modal", xSpace::SCALAR, nb_modes, lss3dCos, lss3dSin),
                                xcApproxFunctionLevelSet(fr), key_modifier);
};

void xcFrontPartLineClose::setFrontParametrisation()
{
   lss1dCos.setSupport(front_region);
   lss1dSin.setSupport(front_region);
   local_vertices = Parametrize1dMeshLoop(start, lss1dCos, lss1dSin, totallenght);
}

void xcFrontPartLineClose::setParametrisation()
{
   setFrontParametrisation();
   lss3dCos.setSupport(xRegion(&mesh, subset_for_q_label));
   xVelocityInitialization lsscos_i(lss1dCos, xRegion(&mesh, subset_for_q_label));
   lss3dCos.accept(lsscos_i, xRegion(&mesh, subset_for_q_label));

   lss3dSin.setSupport(xRegion(&mesh, subset_for_q_label));
   xVelocityInitialization lsssin_i(lss1dSin, xRegion(&mesh, subset_for_q_label));
   lss3dSin.accept(lsssin_i, xRegion(&mesh, subset_for_q_label));
}

void xcFrontPartLineClose::extendParametrization()
{
   lss1dCos.setSupport(xRegion(&front_mesh));
   lss1dSin.setSupport(xRegion(&front_mesh));
   local_vertices = Parametrize1dMeshLoop(start, lss1dCos, lss1dSin, totallenght);
   lss3dCos.setSupport(xRegion(&mesh, "globalsubmesh"));
   xVelocityInitialization lsscos_i(lss1dCos, xRegion(&mesh, subset_for_q_label));
   lss3dCos.accept(lsscos_i, xRegion(&mesh, subset_for_q_label));
   lss3dSin.setSupport(xRegion(&mesh, "globalsubmesh"));
   xVelocityInitialization lsssin_i(lss1dSin, xRegion(&mesh, subset_for_q_label));
   lss3dSin.accept(lsssin_i, xRegion(&mesh, subset_for_q_label));
   fr.setSupport(xRegion(&mesh, "globalsubmesh"));
   xcSetRadialFunction set_fr;
   fr.accept(set_fr, xRegion(&mesh, subset_for_q_label));

   int verbose = parameters.getInt("sifs_verbosity");
   if (verbose)
   {
      xexport::xExportGmshAscii pexport;
      xIntegrationRuleBasic intr;
      std::stringstream outfilename;
      if (ParUtil::Instance()->size() > 1)
         outfilename << front_part_name << "_fr_extend_" << ParUtil::Instance()->size() << "_" << ParUtil::Instance()->rank() + 1;
      else
         outfilename << front_part_name << "_fr_extend";
      Export(fr, pexport, outfilename.str());
   }
}

const xLevelSet &xcFrontPartLineClose::getLss1dCos() const { return lss1dCos; }
const xLevelSet &xcFrontPartLineClose::getLss1dSin() const { return lss1dSin; }

void xcFrontPartLineClose::resetModes(double minwidth)
{
   // computes the number of modes to ensure a minimum mode width of minwidth
   // for the "spline" mode type
   nb_modes = (int)(ceil(2 * getFrontLenght() / minwidth));
   int moduloo = nb_modes % 2;
   if (moduloo != 0)
   {
      // cout << "Warning: number of modes modified to be an even number, using (number_of_modes + 1)" << endl;
      nb_modes += 1;
   }
   mode_width = 2 * getFrontLenght() / nb_modes;

   // delete possible previous modal spaces
   if (modal_space) delete modal_space;
   if (modal_space_1d) delete modal_space_1d;

   // create new modal space
   xValKeyExtend key_modifier("_mult_fr_for_" + getFrontName());
   modal_space = new xSpaceXFEM(xcSpaceModalPolynomeHermiteUniform("Modal_Damage_Int", xSpace::SCALAR, nb_modes, lss3dCos,
                                                                   lss3dSin, mode_width / totallenght * 2),
                                xcApproxFunctionLevelSet(fr), key_modifier);
   modal_space_1d = new xcSpaceModalPolynomeHermiteUniform("Modal_Damage_One_Dimensional", xSpace::SCALAR, nb_modes, lss1dCos,
                                                           lss1dSin, mode_width / totallenght * 2);
}

xApproxFunction *xcFrontPartLineClose::getModalApproxFunction(int imod, int dim)
{
   xApproxFunction *modalpoly;
   switch (dim)
   {
      case 1:
         modalpoly = new xcApproxFunctionModalPolynomeHermiteUniform(lss1dCos, lss1dSin, imod, nb_modes,
                                                                     mode_width / getFrontLenght() * 2.);
         break;
      case 2:
         modalpoly = new xcApproxFunctionModalPolynomeHermiteUniform(lss3dCos, lss3dSin, imod, nb_modes,
                                                                     mode_width / getFrontLenght() * 2.);
         break;
      default:
         cerr << "wrong dimension in xcFrontPartLineClose::getModalApproxFunction" << endl;
         throw;
   }
   return modalpoly;
}

void xcFrontPart::createDomainsAndRadialFunction()
{
   // subset_for_q is a region in the mesh created to define the virtual velocity q
   // it is created by first taking the xUpperAdjacencyRecursive() operator on the
   // element of the front
   // then adding nb_layers_core     --> the set of elements obtained at this point is called foo
   // then adding nb_layers_cylinder
   // then adding layers provided the element has at least one node whose distance to the front is less
   // than rho_cylinder.
   // subset_for_int is defined as subset_for_q minus foo
   // remark : the region "cilynder_core" is a temporary region used to produced the above.

   double rho_cylinder = parameters.getDouble("sifs_rho_geo");
   int nb_layers_cylinder = parameters.getInt("sifs_nb_layers_cylinder");
   int nb_layers_core = parameters.getInt("sifs_nb_layers_core");
   int verbose = parameters.getInt("sifs_verbosity");

   if (verbose)
      cout << "parameters forConfigForce    Subset :" << rho_cylinder << " " << nb_layers_cylinder << " " << nb_layers_core
           << endl;

   // int dim_front = front_mesh.dim();
   int dim_mesh = mesh.dim();

   // (const_cast<xMesh &>(mesh)).modifyState(0, dim_mesh, true);
   (const_cast<xMesh &>(mesh)).getMesh().modifyState(0, dim_mesh, true);

   mesh.createSubMesh(subset_for_int_label);
   xcElementsAlongFrontCreator seed_creator(front_mesh, front_part_name, front_distance);
   xSubMesh &cylinder_core = mesh.createSubMesh("cylinder_core", seed_creator);

   stringstream fname;
   if (verbose) cylinder_core.exportGmsh(front_part_name + "_cylinder_core.msh");

   xAddLayerModifier add_core_layers(nb_layers_core);
   add_core_layers.modify(mesh, "cylinder_core");
   if (verbose) cylinder_core.exportGmsh(front_part_name + "_cylinder_core_plus_layer.msh");

   xAddLayerCreator add_cylinder_topo(nb_layers_cylinder, "cylinder_core");
   xSubMesh &subset_for_q = mesh.createSubMesh(subset_for_q_label, add_cylinder_topo);
   if (verbose) subset_for_q.exportGmsh(front_part_name + "_cylinder_topo.msh");

   int sizesub_loc = subset_for_q.size(dim_mesh);
   int sizesub = sizesub_loc;
#ifdef PARALLEL
   MPI_Allreduce(&sizesub_loc, &sizesub, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD /* ierr*/);
#endif
   int sizesubnew = sizesub;

   xAcceptLessThanForOneNode filter(front_distance, rho_cylinder);
   xAddLayerModifier add_cylinder_geo(1, 0, filter);
   do
   {
      sizesub = sizesubnew;
      // cout << "size "<< sizesub << endl;
      add_cylinder_geo.modify(mesh, subset_for_q_label);
      int sizesubnew_loc = subset_for_q.size(dim_mesh);
      sizesubnew = sizesubnew_loc;
#ifdef PARALLEL
      MPI_Allreduce(&sizesubnew_loc, &sizesubnew, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD /* ierr*/);
#endif
   } while (sizesubnew != sizesub);
   if (verbose) subset_for_q.exportGmsh(front_part_name + "_cylinder_geo.msh");

   // subset_for_q minus cylinder_core : subset_for_int

   if (damage_type)
   {  // damage case: nb_layers_core considered zero, the domain for integration is subset_for_q
      xCopyCreator copycreator(subset_for_q_label);
      xSubMesh &subset_for_int = mesh.createSubMesh(subset_for_int_label, copycreator);
      if (verbose) subset_for_int.exportGmsh(front_part_name + "_subset_for_int_damage.msh");
   }
   else
   {
      xFirstMinusSecondCreator minus_creator(subset_for_q_label, "cylinder_core");
      xSubMesh &subset_for_int = mesh.createSubMesh(subset_for_int_label, minus_creator);
      if (verbose) subset_for_int.exportGmsh(front_part_name + "_subset_for_int.msh");
   }

   mesh.deleteSubMesh("cylinder_core");
   fr.setSupport(xRegion(&mesh, subset_for_q_label));
   xcSetRadialFunction set_fr;
   fr.accept(set_fr);

   if (verbose)
   {
      xexport::xExportGmshAscii pexport;
      xIntegrationRuleBasic intr;
      std::stringstream outfilename;
      if (ParUtil::Instance()->size() > 1)
         outfilename << front_part_name << "_fr_" << ParUtil::Instance()->size() << "_" << ParUtil::Instance()->rank() + 1;
      else
         outfilename << front_part_name << "_fr";
      cout << "##### export fr " << outfilename.str() << endl;
      Export(fr, pexport, outfilename.str());
   }
}

std::map<std::string, std::pair<mVertex *, mVertex *>> Separate1dMeshBranch(xfem::xMesh &front_mesh,
                                                                            const std::string &front_base_name)
{
   // NB: when considering multiple frontparts merging, a node might be connected to more than 2 edges
   // infinte loops are then observed with the algorithm below.
   // Quick Fix used: slightly shift the level set before calling xcFront...
   const bool debug = true;
   xinterface::aomd::xAttachedDataManagerAOMD<int> visited;

   std::string subsetnamebase(front_base_name + "_part_");
   std::map<std::string, std::pair<mVertex *, mVertex *>> parts;
   if (debug) cout << front_mesh.dim() << endl;
   switch (front_mesh.dim())
   {
      case 0:
      {
         int partnumber = 1;
         for (AOMD::mEntity *pv : front_mesh.range(0))
         {
            std::stringstream subsetname;
            subsetname << subsetnamebase << partnumber;
            xSubMesh &subset = front_mesh.createSubMesh(subsetname.str());
            subset.add(pv);
            parts.insert(std::make_pair(subsetname.str(), std::make_pair((AOMD::mVertex *)(pv), (AOMD::mVertex *)nullptr)));
            ++partnumber;
         }
         break;
      }
      case 1:
      {
         std::list<mVertex *> totreatv;
         front_mesh.getMesh().modifyState(0, 1, true);
         for (AOMD::mEntity *pv : front_mesh.range(0)) totreatv.push_back(static_cast<mVertex *>(pv));
         int partnumber = 0;
         if (debug) cout << "Separate 1d mesh Branch start " << endl;
         while (totreatv.size() != 0)
         {  // continue untill every node has been stored in a front part
            ++partnumber;
            std::stringstream subsetname;
            subsetname << subsetnamebase << partnumber;
            xSubMesh &subset = front_mesh.createSubMesh(subsetname.str());
            // cout << "totreat " <<totreatv.size() << endl;
            std::list<mVertex *>::iterator itv = totreatv.begin();
            mVertex *starttmp = (*itv);
            mVertex *start = nullptr;
            mVertex *theend = nullptr;
            mVertex *current = starttmp;
            mEdge *previousedge = nullptr;
            bool isloop = false;
            subset.add(starttmp);
            visited.setData(*starttmp) = 1;
            if (starttmp->size(1) == 1)
            {
               if (debug) cout << "Separate 1d mesh Branch start on boundary" << endl;
               start = starttmp;
               mEdge *currentedge = static_cast<mEdge *>(current->get(1, 0));
               if (currentedge == previousedge)
               {
                  currentedge = static_cast<mEdge *>(current->get(1, 1));
               }
               current = static_cast<mVertex *>(E_otherVertex(currentedge, current));
               subset.add(current);
               subset.add(currentedge);
               visited.setData(*current) = 1;
               previousedge = currentedge;
               while (current->size(1) > 1)
               {
                  mEdge *currentedge = static_cast<mEdge *>(current->get(1, 0));
                  if (currentedge == previousedge)
                  {
                     currentedge = static_cast<mEdge *>(current->get(1, 1));
                  }
                  current = static_cast<mVertex *>(E_otherVertex(currentedge, current));

                  subset.add(current);
                  subset.add(currentedge);
                  visited.setData(*current) = 1;
                  previousedge = currentedge;
               }
               theend = current;
            }
            else
            {
               if (debug) cout << "Separate 1d mesh Branch start on a middle" << endl;
               while (current->size(1) > 1)
               {
                  mEdge *currentedge = static_cast<mEdge *>(current->get(1, 0));
                  if (currentedge == previousedge)
                  {
                     currentedge = static_cast<mEdge *>(current->get(1, 1));
                  }
                  current = static_cast<mVertex *>(E_otherVertex(currentedge, current));
                  subset.add(currentedge);
                  if (current == starttmp) break;
                  subset.add(current);
                  visited.setData(*current) = 1;
                  previousedge = currentedge;
               }
               if (current == starttmp)
               {
                  start = starttmp;
                  theend = starttmp;
               }
               else
               {
                  if (debug) cout << "Separate 1d mesh Branch going the other way" << endl;
                  theend = current;
                  previousedge = nullptr;
                  current = starttmp;
                  while (current->size(1) > 1)
                  {
                     mEdge *currentedge = static_cast<mEdge *>(current->get(1, 1));
                     if (currentedge == previousedge)
                     {
                        currentedge = static_cast<mEdge *>(current->get(1, 0));
                     }
                     current = static_cast<mVertex *>(E_otherVertex(currentedge, current));
                     visited.setData(*current) = 1;
                     subset.add(current);
                     subset.add(currentedge);
                     previousedge = currentedge;
                  }
                  start = current;
               }
            }
            isloop = (theend == start);
            if (debug) cout << "Separate 1d mesh Branch find branch " << start << " " << theend << " " << isloop << endl;
            parts.insert(std::make_pair(subsetname.str(), std::make_pair(start, theend)));
            totreatv.erase(std::remove_if(totreatv.begin(), totreatv.end(),
                                          [&visited](AOMD::mVertex *pv) { return bool(visited.getData(*pv)); }),
                           totreatv.end());

            // front_mesh.modifyState_sub(1,0,true, subsetname.str());
            subset.modifyState(0, 1, true);
         }
         break;
      }
      default:
      {
         cout << "Error : front_mesh should be of dimension 0 or 1 " << endl;
      }
   }
   return parts;
}

std::list<mVertex *> Parametrize1dMeshBranch(mVertex *vstart, mVertex *vend, xLevelSet &ls, double &totallenght)
{
   xinterface::aomd::xAttachedDataManagerAOMD<double> s;
   mVertex *vcurrent = vstart;
   s.setData(*vcurrent) = 0.;
   mEdge *ecurrent = static_cast<mEdge *>(vcurrent->get(1, 0));
   totallenght = 0.;
   std::list<mVertex *> visitednode;
   visitednode.push_back(vcurrent);
   while (vcurrent != vend)
   {
      mVertex *vnext = static_cast<mVertex *>(E_otherVertex(ecurrent, vcurrent));
      xtensor::xVector<> v01(vnext->point(), vcurrent->point());
      totallenght += v01.mag();
      s.setData(*vnext) = totallenght;
      mEdge *enext = static_cast<mEdge *>(vnext->get(1, 0));
      if (enext == ecurrent)
      {
         if (vnext->size(1) > 1) enext = static_cast<mEdge *>(vnext->get(1, 1));
      }
      ecurrent = enext;
      vcurrent = vnext;
      visitednode.push_back(vcurrent);
   }

   for (AOMD::mVertex *pv : visitednode) ls(pv) = s.at(*pv) * 2. / totallenght - 1.;
   return visitednode;
}

std::list<mVertex *> Parametrize1dMeshLoop(mVertex *vstart, xLevelSet &lscos, xLevelSet &lssin, double &totallenght)
{
   xinterface::aomd::xAttachedDataManagerAOMD<double> s;
   mVertex *vcurrent = vstart;
   s.setData(*vcurrent) = 0.;
   mEdge *ecurrent = static_cast<mEdge *>(vcurrent->get(1, 0));
   totallenght = 0.;
   std::list<mVertex *> visitednode;
   visitednode.push_back(vcurrent);
   while (1)
   {
      mVertex *vnext = static_cast<mVertex *>(E_otherVertex(ecurrent, vcurrent));
      xtensor::xVector<> v01(vnext->point(), vcurrent->point());
      totallenght += v01.mag();
      if (vnext == vstart) break;
      mEdge *enext = static_cast<mEdge *>(vnext->get(1, 0));
      if (enext == ecurrent)
      {
         enext = static_cast<mEdge *>(vnext->get(1, 1));
      }
      ecurrent = enext;
      vcurrent = vnext;
      s.setData(*vcurrent) = totallenght;
      visitednode.push_back(vcurrent);
   }

   for (AOMD::mEntity *pv : visitednode)
   {
      const double theta = s.at(*pv) * 2 * M_PI / totallenght;
      lscos(pv) = cos(theta);
      lssin(pv) = sin(theta);
   }
   return visitednode;
}

void xcElementsAlongFrontCreator::create(const xMesh &_m, const string &_name)
{
   name = _name;
   xMesh &m = const_cast<xMesh &>(_m);
   xSubMesh &sub = m.getSubMesh(name);
   xSubMesh &front_part = front_mesh.getSubMesh(front_part_name);

   int dim_front = front_mesh.dim();
   int dim_mesh = m.dim();
   assert((dim_mesh - dim_front) == 1 || (dim_mesh - dim_front) == 2);
   for (mEntity *e : front_part.range(dim_front))
   {
      mVertex *vmin = nullptr;
      double rmin = 1.;  // fake value to avoid -Wmaybe-uninitialized below. Anyway when j= 0 rmin is set with r ...
      e = xCreatorRecursive()(e);
      // cout<< e << " in xcElementsAlongFrontCreator  "  << endl;
      if (e != nullptr)
      {
         // cout<< e << " " << e->size(0) << endl;
         //	  if ((dim_mesh - dim_front) == 2)  e = xUpperCreator()(e);
         for (int j = 0; j < e->size(0); j++)
         {
            mVertex *v = static_cast<mVertex *>(e->get(0, j));
            double r = front_distance(v);
            if (j == 0 || r < rmin)
            {
               vmin = v;
               rmin = r;
            }
         }
         if (vmin == nullptr) vmin = dynamic_cast<mVertex *>(e);  // if e->size(0)==0, e is a vertex
         added.setData(*vmin) = 1;
         toclean.push_back(vmin);
      }
   }

   m.getMesh().exchangeDataOnPartBdrys(const_cast<xcElementsAlongFrontCreator &>(*this), true);
   for (mEntity *v : toclean)
   {
      for (int j = 0; j < v->size(dim_mesh); j++)
      {
         mEntity *ee = v->get(dim_mesh, j);
         sub.add(ee);
      }
      added.deleteData(*v);
   }
   toclean.clear();
   ParUtil::Instance()->Barrier(33, "b2");
   sub.modifyAllState();
   added.clear();
}

int xcElementsAlongFrontCreator::tag() const { return 99996; };

void *xcElementsAlongFrontCreator::AP_alloc_and_fill_buffer(mEntity *e, AOMD::AOMD_SharedInfo &si, int tag)
{
#ifdef PARALLEL
   int buffer_size = sizeof(mEntity *) + sizeof(int);
   void *buf = AP_alloc(si.pid(), this->tag(), buffer_size);
   mEntity **ebuf = reinterpret_cast<mEntity **>(buf);
   *(ebuf++) = si.getRemotePointer();
   int *ibuf = reinterpret_cast<int *>(ebuf);
   *ibuf = added.at(*e);
   return buf;
#else
   throw;
   return nullptr;
#endif
};
void xcElementsAlongFrontCreator::receiveData(int pid, void *buf)
{
   //#ifdef PARALLEL
   mEntity **ebuf = reinterpret_cast<mEntity **>(buf);
   mEntity *vmin = *(ebuf++);
   int *tag = reinterpret_cast<int *>(ebuf);

   if (*tag)
   {
      toclean.push_back(vmin);
      added.setData(*vmin) = 1;
   }
   //#endif
};

xcSetRadialFunction::xcSetRadialFunction(){};
void xcSetRadialFunction::visit(xLevelSet &f, xRegion target)
{
   for (mEntity *pe : target.range(0)) f(pe) = 1;

   xAcceptOnBoundaryOfSubMesh filter_bnd_region(*(target.getSubMesh()));
   //    xAcceptOnBoundary filter_bnd;
   int dim = target.dim();
   for (mEntity *e : target.range(dim - 1))
   {
      if (filter_bnd_region(e) && (GEN_type(e->getClassification()) != dim - 1))
      {
         for (int i = 0; i < e->size(0); ++i)
         {
            mEntity *v = e->get(0, i);
            f(v) = 0.;
         }
      }
   }

   int mpi_s;
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_s);
   if (mpi_s > 1) {
    std::cout << "ERROR in :" <<  __FILE__ << __LINE__ << "if parallel, ls need to be unifyed "<< std::endl;
    throw;
     //xParallelLevelSetUnifier up;
     //up.visit(f, target);
   }
   return;
}

void xcApproxFunctionLevelSet::getVal(const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ, double &res) const
{
   mEntity *e = geo_appro->getEntity();
   if (ls.getSupport().IsInRegion(e))
      res = ls.getVal(geo_appro->getEntity(), geo_appro->getUVW());
   else
      res = 0.;
}

void xcApproxFunctionLevelSet::getGrad(const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ,
                                       xtensor::xVector<> &res) const
{
   mEntity *e = geo_appro->getEntity();
   if (ls.getSupport().IsInRegion(e))
      res = ls.getGrad(geo_appro->getEntity());
   else
      res = xtensor::xVector<>(0., 0., 0.);
}

void xcApproxFunctionGradLevelSet::getVal(const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ,
                                          xtensor::xVector<> &res) const
{
   mEntity *e = geo_appro->getEntity();
   if (ls.getSupport().IsInRegion(e))
      res = ls.getGrad(geo_appro->getEntity(), geo_appro->getUVW());
   else
      res = xtensor::xVector<>(0., 0., 0.);
}

void xcApproxFunctionGradLevelSet::getGrad(const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ,
                                           xtensor::xTensor2<> &res) const
{
   mEntity *e = geo_appro->getEntity();
   if (ls.getSupport().IsInRegion(e))
      res = ls.getCurv(geo_appro->getEntity());
   else
      res = xtensor::xTensor2<>(0.);
}

void xcApproxFunctionNormedGradLevelSet::getVal(const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ,
                                                xtensor::xVector<> &res) const
{
   mEntity *e = geo_appro->getEntity();
   if (ls.getSupport().IsInRegion(e))
      res = (ls.getGrad(geo_appro->getEntity(), geo_appro->getUVW())).norm();
   else
      res = xtensor::xVector<>(0., 0., 0.);
}

void xcApproxFunctionNormedGradLevelSet::getGrad(const xfem::xGeomElem *geo_appro, const xfem::xGeomElem *geo_integ,
                                                 xtensor::xTensor2<> &res) const
{
   mEntity *e = geo_appro->getEntity();
   if (ls.getSupport().IsInRegion(e))
      res = ls.getCurv(geo_appro->getEntity());
   else
      res = xtensor::xTensor2<>(0.);
}

}  // end of namespace xcrack
