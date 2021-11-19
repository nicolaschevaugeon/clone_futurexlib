/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

// xfem
#ifdef HAVE_ANN
#include "xNearestNeighborInterface.h"
#endif
#include "xMesh.h"
#include "xParseData.h"
// xexport
#include "xExportAlgorithm.h"
// xcrack
#include "xcFrontDomain.h"
#include "xcFrontPart.h"
#include "xcFrontSpaces.h"

using namespace AOMD;
using namespace xfem;
using namespace std;

xcFrontDomain::xcFrontDomain(xcFrontPartBase* _part, int _physical_process)
    : subset_for_q_label(_part->getFrontName() + "_subset_for_q"),
      subset_for_int_label(_part->getFrontName() + "_subset_for_int"),
      front_distance(_part->getDistanceFunction()),
      front_part_name(_part->getFrontName()),
      part(_part),
      physical_process(_physical_process),
      parameters(_part->getParameters()),
      mesh(_part->getMesh()),
      front_mesh(_part->getFrontMesh())
{
   createDomainAndRadialFunction();
};

xcFrontDomain::~xcFrontDomain() = default;
;

void xcFrontDomain::createDomainAndRadialFunction()
{
   if (!(physical_process == DAMAGEFRONT))
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

      double rho_cylinder = parameters.getDouble("front_domain_rho_geo");
      int nb_layers_cylinder = parameters.getInt("front_domain_nb_layers_cylinder");
      int nb_layers_core = parameters.getInt("front_domain_nb_layers_core");
      int verbose = parameters.getInt("front_verbosity");

      if (verbose)
         cout << "parameters forConfigForce    Subset :" << rho_cylinder << " " << nb_layers_cylinder << " " << nb_layers_core
              << endl;

      // int dim_front = front_mesh.dim();
      int dim_mesh = mesh.dim();

      (const_cast<xMesh&>(mesh)).getMesh().modifyState(0, dim_mesh, true);

      mesh.createSubMesh(subset_for_int_label);
      xcElementsAlongFrontCreator seed_creator(front_mesh, front_part_name, front_distance);
      xSubMesh& cylinder_core = mesh.createSubMesh("cylinder_core", seed_creator);

      stringstream fname;
      if (verbose) cylinder_core.exportGmsh(front_part_name + "_cylinder_core.msh");

      xAddLayerModifier add_core_layers(nb_layers_core);
      add_core_layers.modify(mesh, "cylinder_core");
      if (verbose) cylinder_core.exportGmsh(front_part_name + "_cylinder_core_plus_layer.msh");

      xAddLayerCreator add_cylinder_topo(nb_layers_cylinder, "cylinder_core");
      xSubMesh& subset_for_q = mesh.createSubMesh(subset_for_q_label, add_cylinder_topo);
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

      if (physical_process == DAMAGEFRONT)
      {  // damage case: nb_layers_core considered zero, the domain for integration is subset_for_q
         xCopyCreator copycreator(subset_for_q_label);
         xSubMesh& subset_for_int = mesh.createSubMesh(subset_for_int_label, copycreator);
         if (verbose) subset_for_int.exportGmsh(front_part_name + "_subset_for_int_damage.msh");
      }
      else
      {  // physical_process==CRACKFRONT
         xFirstMinusSecondCreator minus_creator(subset_for_q_label, "cylinder_core");
         xSubMesh& subset_for_int = mesh.createSubMesh(subset_for_int_label, minus_creator);
         if (verbose) subset_for_int.exportGmsh(front_part_name + "_subset_for_int.msh");
      }

      mesh.deleteSubMesh("cylinder_core");
      fr.setSupport(xRegion(&mesh, subset_for_q_label));
      // cout << "###########################################################" << endl;
      xcSetRadialFunction set_fr;
      fr.accept(set_fr);

      if (verbose)
      {
         xexport::xExportGmshAscii pexport;
         xIntegrationRuleBasic intr;
         std::stringstream outfilename;
         if (AOMD::ParUtil::Instance()->size() > 1)
            outfilename << front_part_name << "_fr_" << AOMD::ParUtil::Instance()->size() << "_"
                        << AOMD::ParUtil::Instance()->rank() + 1;
         else
            outfilename << front_part_name << "_fr";
         cout << "##### export fr " << outfilename.str() << endl;
         xexport::Export(fr, pexport, outfilename.str());
      }
   }
   else
   {  // ----------------------------------- DAMAGE FRONT USING ANN
      double fr_distance = parameters.getDouble("front_domain_fr_distance");
      double rho_cylinder = parameters.getDouble("front_domain_rho_geo");
      int verbose = parameters.getInt("front_verbosity");

      if (verbose) cout << "parameters forConfigForce    Subset :" << rho_cylinder << endl;

      int dim_mesh = mesh.dim();
      (const_cast<xMesh&>(mesh)).getMesh().modifyState(0, dim_mesh, true);

      mesh.createSubMesh(subset_for_q_label);
      xcElementsAtLcDistanceFrom1DFront ann_creator(front_mesh, front_part_name, rho_cylinder, fr_distance, front_distance);
      xSubMesh& subset_for_q = mesh.createSubMesh(subset_for_q_label, ann_creator);
      if (verbose) subset_for_q.exportGmsh(front_part_name + "_cylinder_topo.msh");

      xCopyCreator copycreator(subset_for_q_label);
      xSubMesh& subset_for_int = mesh.createSubMesh(subset_for_int_label, copycreator);
      if (verbose) subset_for_int.exportGmsh(front_part_name + "_subset_for_int_damage.msh");

      fr.setSupport(xRegion(&mesh, subset_for_q_label));
      subset_for_q.updateBoundary();
      xcSetRadialFunction set_fr;
      // cout <<"***********************************" << endl;
      fr.accept(set_fr);
   }
}

///------------------------------------------------------------

xcSetRadialFunction::xcSetRadialFunction(){};
void xcSetRadialFunction::visit(xLevelSet& f, xRegion target)
{
   for (mEntity* pe : target.range(0)) f(pe) = 1.;

   xAcceptOnBoundaryOfSubMesh filter_bnd_region(*(target.getSubMesh()));
   //    xAcceptOnBoundary filter_bnd;
   int dim = target.dim();
   for (mEntity* e : target.range(dim - 1))
   {
      if (filter_bnd_region(e) && (GEN_type(e->getClassification()) != dim - 1))
      {
         for (int i = 0; i < e->size(0); ++i)
         {
            mEntity* v = e->get(0, i);
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

///------------------------------------------------------------

xcFrontDomainPoint::xcFrontDomainPoint(xcFrontPartBase* _part, int _physical_process) : xcFrontDomain(_part, _physical_process) {}

void xcFrontDomainPoint::extendParametrization() {}

///------------------------------------------------------------

xcFrontDomainLineOpen::xcFrontDomainLineOpen(xcFrontPartBase* _part, int _physical_process)
    : xcFrontDomain(_part, _physical_process)
{
   lss3d.setSupport(xRegion(&mesh, subset_for_q_label));
   xVelocityInitialization lss_i(((xcFrontPartLineOpen*)part)->getLss1d(), xRegion(&mesh, subset_for_q_label));
   lss3d.accept(lss_i, xRegion(&mesh, subset_for_q_label));

   int verbose = parameters.getInt("front_verbosity");
   if (verbose)
   {
      xexport::xExportGmshAscii pexport;
      xIntegrationRuleBasic intr;
      std::stringstream outfilename("lss3d");
      if (AOMD::ParUtil::Instance()->size() > 1)
         outfilename << "_" << AOMD::ParUtil::Instance()->size() << "_" << AOMD::ParUtil::Instance()->rank() + 1 << ".pos";
      xexport::Export(lss3d, pexport, outfilename.str());
      cout << "export lss3d" << endl;
   }
}

void xcFrontDomainLineOpen::extendParametrization()
{
   ((xcFrontPartLineOpen*)part)->restrict_ls1d_support(xRegion(&front_mesh));
   lss3d.restrictTo(xRegion(&mesh, "globalsubmesh"), 0.);
   fr.restrictTo(xRegion(&mesh, subset_for_q_label), 0.);

   /*	lss1d.setSupport(xRegion(&front_mesh));
           local_vertices = Parametrize1dMeshBranch(start, end, lss1d, totallenght);
           lss3d.setSupport(xRegion(&mesh, "globalsubmesh"));
           xVelocityInitialization lss_i(lss1d, xRegion(&mesh, subset_for_q_label));
           lss3d.accept( lss_i, xRegion(&mesh, subset_for_q_label));
           fr.setSupport(xRegion(&mesh, "globalsubmesh"));
           xcSetRadialFunction set_fr;
           fr.accept(set_fr, xRegion(&mesh, subset_for_q_label));  */

   int verbose = parameters.getInt("front_verbosity");
   if (verbose)
   {
      xexport::xExportGmshAscii pexport;
      xIntegrationRuleBasic intr;
      std::stringstream outfilename;
      std::stringstream outfilename_ls3d;
      if (AOMD::ParUtil::Instance()->size() > 1)
      {
         outfilename << front_part_name << "_fr_extend_" << AOMD::ParUtil::Instance()->size() << "_"
                     << AOMD::ParUtil::Instance()->rank() + 1;
         outfilename_ls3d << front_part_name << "_ls3d_extended_" << AOMD::ParUtil::Instance()->size() << "_"
                          << AOMD::ParUtil::Instance()->rank() + 1;
      }
      else
      {
         outfilename << front_part_name << "_fr_extend";
         outfilename_ls3d << front_part_name << "_ls3d_extended";
      }
      xexport::Export(fr, pexport, outfilename.str());
      xexport::Export(lss3d, pexport, outfilename_ls3d.str());
      cout << "export extended lss3d and fr" << endl;
   }
};

xcApproxFunctionModalPolynomeHermiteUniform* xcFrontDomainLineOpen::getNewApproxFunction(int imod, int dim)
{
   double mode_width = 4. / (nb_modes - 1);
   switch (dim)
   {
      case 1:
         return (new xcApproxFunctionModalPolynomeHermiteUniform(((xcFrontPartLineOpen*)part)->getLss1d(), imod, nb_modes,
                                                                 mode_width));
         break;
      case 2:
         return (new xcApproxFunctionModalPolynomeHermiteUniform(lss3d, imod, nb_modes, mode_width));
         break;
      default:
         throw;
   }
}

///------------------------------------------------------------

xcFrontDomainLineClosed::xcFrontDomainLineClosed(xcFrontPartBase* _part, int _physical_process)
    : xcFrontDomain(_part, _physical_process)
{
   lss3dCos.setSupport(xRegion(&mesh, subset_for_q_label));
   xVelocityInitialization lsscos_i(((xcFrontPartLineClose*)part)->getLss1dCos(), xRegion(&mesh, subset_for_q_label));
   lss3dCos.accept(lsscos_i, xRegion(&mesh, subset_for_q_label));

   lss3dSin.setSupport(xRegion(&mesh, subset_for_q_label));
   xVelocityInitialization lsssin_i(((xcFrontPartLineClose*)part)->getLss1dSin(), xRegion(&mesh, subset_for_q_label));
   lss3dSin.accept(lsssin_i, xRegion(&mesh, subset_for_q_label));

   int verbose = parameters.getInt("front_verbosity");
   if (verbose)
   {
      xexport::xExportGmshAscii pexport;
      xIntegrationRuleBasic intr;
      std::stringstream outfilenameCos("lss3dCos");
      std::stringstream outfilenameSin("lss3dSin");
      if (AOMD::ParUtil::Instance()->size() > 1)
      {
         outfilenameCos << "_" << AOMD::ParUtil::Instance()->size() << "_" << AOMD::ParUtil::Instance()->rank() + 1 << ".pos";
         outfilenameSin << "_" << AOMD::ParUtil::Instance()->size() << "_" << AOMD::ParUtil::Instance()->rank() + 1 << ".pos";
      }
      xexport::Export(lss3dCos, pexport, outfilenameCos.str());
      xexport::Export(lss3dSin, pexport, outfilenameSin.str());
      cout << "export lss3d Cos and Sin" << endl;
   }
}

void xcFrontDomainLineClosed::extendParametrization()
{
   ((xcFrontPartLineClose*)part)->restrict_ls1d_support(xRegion(&front_mesh));
   lss3dCos.restrictTo(xRegion(&mesh, "globalsubmesh"), 0.);
   lss3dSin.restrictTo(xRegion(&mesh, "globalsubmesh"), 0.);
   fr.restrictTo(xRegion(&mesh, subset_for_q_label), 0.);

   /*	lss1d.setSupport(xRegion(&front_mesh));
           local_vertices = Parametrize1dMeshBranch(start, end, lss1d, totallenght);
           lss3d.setSupport(xRegion(&mesh, "globalsubmesh"));
           xVelocityInitialization lss_i(lss1d, xRegion(&mesh, subset_for_q_label));
           lss3d.accept( lss_i, xRegion(&mesh, subset_for_q_label));
           fr.setSupport(xRegion(&mesh, "globalsubmesh"));
           xcSetRadialFunction set_fr;
           fr.accept(set_fr, xRegion(&mesh, subset_for_q_label));  */

   int verbose = parameters.getInt("front_verbosity");
   if (verbose)
   {
      xexport::xExportGmshAscii pexport;
      xIntegrationRuleBasic intr;
      std::stringstream outfilename;
      std::stringstream outfilename_ls3dcos;
      std::stringstream outfilename_ls3dsin;
      if (AOMD::ParUtil::Instance()->size() > 1)
      {
         outfilename << front_part_name << "_fr_extend_" << AOMD::ParUtil::Instance()->size() << "_"
                     << AOMD::ParUtil::Instance()->rank() + 1;
         outfilename_ls3dcos << front_part_name << "_ls3dcos_extended_" << AOMD::ParUtil::Instance()->size() << "_"
                             << AOMD::ParUtil::Instance()->rank() + 1;
         outfilename_ls3dsin << front_part_name << "_ls3dsin_extended_" << AOMD::ParUtil::Instance()->size() << "_"
                             << AOMD::ParUtil::Instance()->rank() + 1;
      }
      else
      {
         outfilename << front_part_name << "_fr_extend";
         outfilename_ls3dcos << front_part_name << "_ls3dcos_extended";
         outfilename_ls3dsin << front_part_name << "_ls3dsin_extended";
      }
      xexport::Export(fr, pexport, outfilename.str());
      xexport::Export(lss3dCos, pexport, outfilename_ls3dcos.str());
      xexport::Export(lss3dSin, pexport, outfilename_ls3dsin.str());
      cout << "export extended lss3d's and fr" << endl;
   }
}

xcApproxFunctionModalPolynomeHermiteUniform* xcFrontDomainLineClosed::getNewApproxFunction(int imod, int dim)
{
   double mode_width = 4. / nb_modes;
   switch (dim)
   {
      case 1:
         return (new xcApproxFunctionModalPolynomeHermiteUniform(((xcFrontPartLineClose*)part)->getLss1dCos(),
                                                                 ((xcFrontPartLineClose*)part)->getLss1dSin(), imod, nb_modes,
                                                                 mode_width));
         break;
      case 2:
         return (new xcApproxFunctionModalPolynomeHermiteUniform(lss3dCos, lss3dSin, imod, nb_modes, mode_width));
         break;
      default:
         throw;
   }
}

///------------------------------------------------------------

void xcElementsAlongFrontCreator::create(const xMesh& _m, const string& _name)
{
   name = _name;
   xMesh& m = const_cast<xMesh&>(_m);
   xSubMesh& sub = m.getSubMesh(name);
   xSubMesh& front_part = front_mesh.getSubMesh(front_part_name);

   int dim_front = front_mesh.dim();
   int dim_mesh = m.dim();
   assert((dim_mesh - dim_front) == 1 || (dim_mesh - dim_front) == 2);
   for (mEntity* e : front_part.range(dim_front))
   {
      mVertex* vmin = nullptr;
      double rmin = 1.;  // fake value to avoid -Wmaybe-uninitialized below. Anyway when j= 0 rmin is set with r ...
      e = xCreatorRecursive()(e);
      if (e)
      {
         for (int j = 0; j < e->size(0); j++)
         {
            mVertex* v = static_cast<mVertex*>(e->get(0, j));
            const double r = front_distance(v);
            if (j == 0 || r < rmin)
            {
               vmin = v;
               rmin = r;
            }
         }
         if (vmin == nullptr) vmin = dynamic_cast<mVertex*>(e);  // if e->size(0)==0, e is a vertex
         added.setData(*vmin) = 1;
         toclean.push_back(vmin);
      }
   }

   m.getMesh().exchangeDataOnPartBdrys(const_cast<xcElementsAlongFrontCreator&>(*this), true);
   for (mEntity* v : toclean)
   {
      for (int j = 0; j < v->size(dim_mesh); j++)
      {
         mEntity* ee = v->get(dim_mesh, j);
         sub.add(ee);
      }
      added.deleteData(*v);
   }
   toclean.clear();
   AOMD::ParUtil::Instance()->Barrier(33, "b2");
   sub.modifyAllState();
   added.clear();
}

int xcElementsAlongFrontCreator::tag() const { return 99996; };

void* xcElementsAlongFrontCreator::AP_alloc_and_fill_buffer(mEntity* e, AOMD::AOMD_SharedInfo& si, int tag)
{
#ifdef PARALLEL
   int buffer_size = sizeof(mEntity*) + sizeof(int);
   void* buf = AP_alloc(si.pid(), this->tag(), buffer_size);
   mEntity** ebuf = reinterpret_cast<mEntity**>(buf);
   *(ebuf++) = si.getRemotePointer();
   int* ibuf = reinterpret_cast<int*>(ebuf);
   *ibuf = added.setData(*e);
   return buf;
#else
   throw;
   return nullptr;
#endif
};

void xcElementsAlongFrontCreator::receiveData(int pid, void* buf)
{
   //#ifdef PARALLEL
   mEntity** ebuf = reinterpret_cast<mEntity**>(buf);
   mEntity* vmin = *(ebuf++);
   int* tag = reinterpret_cast<int*>(ebuf);

   if (*tag)
   {
      toclean.push_back(vmin);
      added.setData(*vmin) = 1;
   }
   //#endif
};

///------------------------------------------------------------

void xcElementsAtLcDistanceFrom1DFront::create(const xMesh& _m, const string& _name)
{
#ifdef HAVE_ANN
   // First, consider only the nodes at a distance ann_distance of the frontmesh for this front part
   // then, use the front_distance function for a second filter
   // The front_distance only is not sufficient to separate the integration zones in case of multiple damage zones...
   name = _name;
   xMesh& m = const_cast<xMesh&>(_m);
   xSubMesh& sub = m.getSubMesh(name);
   xSubMesh& front_part = front_mesh.getSubMesh(front_part_name);
   xNearestNeighborInterface<xIter> annint(front_part.begin(0), front_part.end(0));
   for (mEntity* e : m.range(2))
   {  // for every face in the mesh
      int nb_vertices_on_face = e->size(0);
      double dmin = ann_distance + 1., d;
      for (int inode = 0; inode < nb_vertices_on_face; inode++)
      {
         mVertex* v = (mVertex*)e->get(0, inode);
         // mVertex *nearest_on_front = annint.nearestvertexto(*v, d );
         annint.nearestvertexto(*v, d);
         if (inode == 0)
            dmin = d;
         else
            dmin = min(dmin, d);
         // cout << "inode " << inode << " d " << d << " dmin " << dmin << " ann_distance " << ann_distance << endl;
      }
      double dmin2 = fr_distance + 1.;
      if (dmin <= ann_distance)
      {
         for (int inode = 0; inode < nb_vertices_on_face; inode++)
         {
            mVertex* v = (mVertex*)e->get(0, inode);
            double d = front_distance(v);
            if (inode == 0)
               dmin2 = d;
            else
               dmin2 = min(dmin2, d);
         }
         if (dmin2 <= fr_distance)
         {
            sub.add(e);  // ajout de la face
            for (int inode = 0; inode < nb_vertices_on_face; inode++)
            {
               sub.add(e->get(0, inode));  // ajout des noeuds
            }
         }
      }
   }
   sub.modifyAllState();
#else
   std::cerr << "Approximate Nearest Neighbor Library not defined. Please compile both the Xfem and Xcrack codes with the "
                "-HAVE_ANN=1 flag."
             << std::endl;
   throw;
#endif
   // an old parallele version below. need to be updated
   /*	const unsigned int  added_tag  = AOMD::AOMD_Util::Instance()->loOOokupMeshDataId("added_entity_tag");
   int dim_front = front_mesh.dim();
   int dim_mesh  = m.dim();
   int dim_diff = dim_mesh - dim_front;
   assert(dim_diff == 1 || dim_diff == 2);
   for(xIter it=front_part.begin(dim_front); it!= front_part.end(dim_front); ++it)
   {
       mEntity* e = *it;
       mVertex* vmin=0;
       double rmin;
       e = xCreatorRecursive()(e);
       if (e!=0){
           //cout<< e << " " << e->size(0) << endl;
           //	  if (dim_diff == 2)  e = xUpperCreator()(e);
           for (int j = 0; j < e->size(0); j++)
           {
               mVertex* v = (mVertex*) e->get(0,j);
               double r = front_distance(v);
               if (j == 0 || r < rmin) {vmin = v; rmin = r;}
           }
           if (vmin==0) vmin=dynamic_cast<mVertex * > (e); //if e->size(0)==0, e is a vertex
           vmin->attachInt(added_tag, 1);
           toclean.push_back(vmin);
       }
   }

   m.exchangeDataOnPartBdrys(const_cast< xcElementsAtLcDistanceFrom1DFront& >( *this), true);
   for (std::list<mEntity *>::iterator it = toclean.begin(); it!=toclean.end(); ++it){
       mEntity * v = *it;
       for (int j = 0; j < v->size(dim_mesh); j++) {
           mEntity* ee = v->get(dim_mesh,j);
           sub.add(ee);
       }
       (*it)->deleteData(added_tag);
   }
   toclean.clear();

   //cout << "almost" << endl;
   AOMD::ParUtil::Instance()->Barrier(33, "b2");
   //cout << "done" << endl;
   sub.modifyAllState();*/
}
