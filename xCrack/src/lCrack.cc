/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/
// std
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
// boost
#include <boost/timer.hpp>
// Trellis
#include "ParUtil.h"
#include "mAOMD.h"
// xtensor
#include "xTensor2.h"
#include "xVector.h"
// xfem
#include "xAlgorithm.h"
#include "xEntityFilter.h"
#include "xField.h"
#include "xMaterial.h"
#include "xMaterialManager.h"
#include "xMesh.h"
#include "xMeshCut.h"
#include "xPointToDouble.h"
#include "xSimpleGeometry.h"
// xexport
#include "xExportAlgorithm.h"
#include "xExportGmsh.h"
// xcrack
#include "CrackPostpro.h"
#include "lCrack.h"

using namespace std;
using namespace xfem;
using namespace AOMD;
using AOMD::ParUtil;

namespace xcrack
{
const bool do_reortho = true;
const bool do_reinit = true;
// anthony
const bool do_narrow = false;

const bool iter_lsn = false;

int& lCrack::tstep() { return timestep; }

lCrack::lCrack(xMesh* m, bool fite, std::string _label)
    : fit(fite),
      mesh(m),
      lsn(m),
      lst(m),
      lss_created(false),
      is_closed(false),
      mesh_crack_surface(nullptr),
      mesh_crack_front(nullptr),
      timestep(0)
{
   label = _label;
   set_Subset_Labels();
}

lCrack::lCrack(xMesh* m, const xPointToDouble& lsfn, const xPointToDouble& lsft, bool fite, std::string _label)
    : fit(fite),
      mesh(m),
      lsn(m),
      lst(m),
      lss_created(false),
      is_closed(false),
      mesh_crack_surface(nullptr),
      mesh_crack_front(nullptr),
      timestep(0)
{
   label = _label;
   set_Subset_Labels();
   // std::cout <<" before and after load" << std::endl;
   // system("ps v  -C CrackComput-O ");
   lsn.load(lsfn);
   lst.load(lsft);
   // system("ps v  -C CrackComput-O ");
   // mesh_crack_surface=NULL;
   // mesh_crack_front=NULL;
   post_initfunc();
   // std::cout <<" postinit" << std::endl;
   // system("ps v  -C CrackComput-O ");
}

void lCrack::set_Subset_Labels()
{
   touched_by_crack_label = "supports_touched_by_crack_" + label;
   touched_by_front_label = "supports_touched_by_front_" + label;
   cut_by_crack_label = "supports_cut_by_crack_" + label;
   cut_by_crack_to_delete_label = "supports_cut_by_crack_to_delete_" + label;
   cutthru_by_crack_label = "supports_cutthru_by_crack_" + label;
   integrate_front_label = "supports_integrate_front_" + label;
   narrow_band_label = "narrow_band_" + label;
   all_minus_touched_by_front_label = "all_minus_supports_touched_by_front_" + label;
}

void lCrack::post_initfunc()
{
   if (fit) fittovertices();
   bool debug = false;
   lsn.forceGradComputation();
   lst.forceGradComputation();

   if (debug) std::cout << " postinit 0" << std::endl;
   // system("ps v  -C CrackComput-O ");
   cutMeshAndSlice();
   if (debug) std::cout << " postinit 1" << std::endl;
   // system("ps v  -C CrackComput-O ");
   createSupportInfo();
   if (debug) std::cout << " postinit 2" << std::endl;
   // system("ps v  -C CrackComput-O ");
   createNarrowBand();
   if (debug) std::cout << " postinit 3" << std::endl;
   // system("ps v  -C CrackComput-O ");
}

void lCrack::setLevelSets(xLevelSet& lsfn, xLevelSet& lsft, bool fite)
{
   fit = fite;
   lsn = lsfn;
   lst = lsft;
   if (fit) fittovertices();
   cutMeshAndSlice();
   createSupportInfo();
   createNarrowBand();
   lss_created = false;
   is_closed = false;
   lss.clear();
   lss_cos.clear();
   lss_sin.clear();
   lss_1D.clear();
   lss_1D_cos.clear();
   lss_1D_sin.clear();
   v_iplane1D.clear();
   v_oplane1D.clear();
   v_x1D.clear();
   v_y1D.clear();
   v_z1D.clear();
   v_iplane3D.clear();
   v_oplane3D.clear();
   v_x3D.clear();
   v_y3D.clear();
   v_z3D.clear();
   v1D.clear();
   v3D.clear();
}

double lCrack::propagate()
{
   boost::timer clock;
   double timefromV1DToVioplane = 0.;
   double timeextendVelocity = 0.;
   double timelspropagate = 0.;
   double timereinit = 0.;
   double timereortho = 0.;
   double timemeshstuff = 0.;

   tstep()++;

   fromV1DToVioplane();
   timefromV1DToVioplane += clock.elapsed();
   clock.restart();
   // exportViplane1D("v_iplane1D");
   // exportVoplane1D("v_oplane1D");
   std::cout << "extendVelocity()" << std::endl;
   extendVelocity();
   timeextendVelocity += clock.elapsed();
   clock.restart();

   if (iter_lsn)
   {
      modifyVoplane();
      exportVoplane3D("v_oplane3D_aftermodify");
   }
   std::cout << "lspropagate()" << std::endl;
   double dt = lspropagate();
   timelspropagate += clock.elapsed();
   clock.restart();
   // exportlst("lst_prop");
   // exportlsn("lsn_prop");

#ifdef PARALLEL
   ParUtil::Instance()->Barrier(__LINE__, __FILE__);
#endif
   std::cout << "reinit     ()" << std::endl;
   // reinit();
   timereinit += clock.elapsed();
   clock.restart();

#ifdef PARALLEL
   ParUtil::Instance()->Barrier(__LINE__, __FILE__);
#endif
   std::cout << "reortho    ()" << std::endl;
   // reortho();
   timereortho += clock.elapsed();
   clock.restart();

#ifdef PARALLEL
   ParUtil::Instance()->Barrier(__LINE__, __FILE__);
#endif

   //  exportlst("lst_reortho");
   //  exportlsn("lsn_reinit");

   if (fit) fittovertices();
   cutMeshAndSlice();
   createSupportInfo();
   createNarrowBand();
   timemeshstuff += clock.elapsed();
   clock.restart();

   double timetotal = timefromV1DToVioplane + timeextendVelocity + timelspropagate + timereinit + timereortho + timemeshstuff;
   std::cout << "time fromV1DToVioplane() : " << timefromV1DToVioplane << " " << timefromV1DToVioplane / timetotal << std::endl;
   std::cout << "time extendVelocity()    : " << timeextendVelocity << " " << timeextendVelocity / timetotal << std::endl;
   std::cout << "time lspropagate()       : " << timelspropagate << " " << timelspropagate / timetotal << std::endl;
   std::cout << "time reinit()            : " << timereinit << " " << timereinit / timetotal << std::endl;
   std::cout << "time reortho()           : " << timereortho << " " << timereortho / timetotal << std::endl;
   std::cout << "time meshstuff()         : " << timemeshstuff << " " << timemeshstuff / timetotal << std::endl;
   std::cout << "time total : " << timetotal << std::endl;

   return dt;
}

void lCrack::fittovertices()
{
   xFitToVertices fit(Param.tolfit);
   lsn.accept(fit);
   lst.accept(fit);
}

double lCrack::reinit2D()
{
   xL2norm l2n;
   double dtvirt = xCFL::getDt(mesh_crack_surface);

#ifdef PARALLEL
   {
      double dtvirtMin = 0.;
      MPI_Allreduce(&dtvirt, &dtvirtMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      dtvirt = dtvirtMin;
   }
#endif

   cout << "dtvirt=" << dtvirt << endl;
   xPilotError statio(l2n, Param.tolerance, Param.max_iter_reinit * 100, dtvirt * Param.coeff_dt_virt);
   //  initialisations du pilot pour la reorthogonalisation

   xReInitOperator reinit;
   xEvolveToStationary rei(reinit, statio);
   if (do_reinit)
   {
      cout << "reinit" << endl;
      lst_crack_surface.accept(rei);
      //    crack3D.accept ( rei , lSupport(mesh, "all_minus_supports_touched_by_front"));
   }
   return dtvirt;
}

double lCrack::reinit()
{
   //  ... reinitialisation des deux levelsets
   // valeurs bloquees sur le front
   //  cFirstMinusSecondCreator cmin("all", "supports_touched_by_front");
   //  mesh->createSubsetEntities("all_minus_supports_touched_by_front", cmin);
   // valeurs bloquees ou non bloquees sur le front

   xL2norm l2n;
   double dtvirt = xCFL::getDt(narrow_band);

#ifdef PARALLEL
   {
      double dtvirtMin = 0.;
      MPI_Allreduce(&dtvirt, &dtvirtMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      dtvirt = dtvirtMin;
   }
#endif

   cout << "dtvirt=" << dtvirt << endl;
   xPilotError statio(l2n, Param.tolerance, Param.max_iter_reinit, dtvirt * Param.coeff_dt_virt);
   //  initialisations du pilot pour la reorthogonalisation

   xReInitOperator reinit;
   xEvolveToStationary rei(reinit, statio);
   if (do_reinit)
   {
      cout << "reinit" << endl;
      lsn.accept(rei, narrow_band);
      //    crack3D.accept ( rei , lSupport(mesh, "all_minus_supports_touched_by_front"));
   }
   return dtvirt;
}

double lCrack::reortho()
{
   xL2norm l2;
   double dtvirt = xCFL::getDt(narrow_band);

#ifdef PARALLEL
   {
      double dtvirtMin = 0.;
      MPI_Allreduce(&dtvirt, &dtvirtMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      dtvirt = dtvirtMin;
   }
#endif

   xPilotError statio(l2, Param.tolerance, Param.max_iter_ortho, dtvirt * Param.coeff_dt_virt * 2.);
   //  ... re-orthogonalisation (avec etape de reinitialisation)
   //  std::list<xSpaceOperator*> ortho_reinit;
   //  ortho_reinit.push_back(new xOrthoOperator(lsn));
   xOrthoOperator orthogo(lsn);
   // if (do_reinit)
   //    ortho_reinit.push_back(new xReInitOperator());
   //  xEvolveToStationary ortho(ortho_reinit, statio2);

   xEvolveToStationary ortho(orthogo, statio);
   // xEvolveToStationaryT <xGaussSeidel <xOrthoOperator >  > ortho(orthogo, statio);

   if (do_reortho)
   {
      cout << "reortho" << endl;
      lst.accept(ortho, narrow_band);
   }
   return dtvirt;
}

double lCrack::lsnextend()
{
   double dtvirt = xCFL::getDt(narrow_band);

#ifdef PARALLEL
   {
      double dtvirtMin = 0.;
      MPI_Allreduce(&dtvirt, &dtvirtMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      dtvirt = dtvirtMin;
   }
#endif

   xL2norm l2;
   // xPilotError  statio(l2, Param.tolerance, Param.max_iter_extend, dtvirt*Param.coeff_dt_virt);
   xPilotError statio(l2, Param.tolerance, 3000, dtvirt * Param.coeff_dt_virt);
   xLSExtensionOperator ext(lst, v3D);
   xEvolveToStationary extension(ext, statio);
   lsn.accept(extension, narrow_band);
   return dtvirt;
}

double lCrack::lspropagate()
{
   int nbsteps;
   double dtstep;

   if (iter_lsn)
   {
      // propagation de lsn (plan de la fissure)
      double dto = xCFL::getDt(v_oplane3D) * Param.coeff_dt_phys;

#ifdef PARALLEL
      {
         double dtoMin = 0.;
         MPI_Allreduce(&dto, &dtoMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
         dto = dtoMin;
      }
#endif

      nbsteps = (int)std::floor(1.0 / dto);
      if (nbsteps == 0) nbsteps = 1;
      //  nbsteps=1;
      dtstep = 1.0 / nbsteps;
      do
      {
         nbsteps--;
         xPilotOneStep pilot_o(dtstep);
         xPropagationOperator pro_o(v_oplane3D);
         xTimeIntegration tim_o(pro_o, pilot_o);
         cout << "propag O" << endl;
         lsn.accept(tim_o);
      } while (nbsteps > 0);
   }

   else
      lsnextend();

   // propagation de lst (position du front sur le plan de la fissure)
   nbsteps = int(Param.coeff_real_dt);
   if (nbsteps == 0) nbsteps = 1;
   double dti = xCFL::getDt(v_iplane3D) * Param.coeff_dt_phys;

#ifdef PARALLEL
   {
      double dtiMin = 0.;
      MPI_Allreduce(&dti, &dtiMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      dti = dtiMin;
   }
#endif

   dtstep = dti;
   double dt = 0;

   do
   {
      nbsteps--;
      dt += dtstep;
      xPilotOneStep pilot_i(dtstep);
      xPropagationOperator pro_i(v_iplane3D);
      xTimeIntegration tim_i(pro_i, pilot_i);
      lst.accept(tim_i, narrow_band);
      // cout <<"propag I" << endl;
   } while (nbsteps > 0);
   cout << dt << endl;
   return dt;
}

double lCrack::modifyVoplane()
{
   //  ... determination du pas de temps critique
   double dti = xCFL::getDt(v_iplane3D);
   double dto = xCFL::getDt(v_oplane3D);

#ifdef PARALLEL
   {
      double dtiMin = 0., dtoMin = 0.;
      MPI_Allreduce(&dti, &dtiMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&dto, &dtoMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      dti = dtiMin;
      dto = dtoMin;
   }
#endif

   //  pas de temps critique optimal
   double dt = min(dti, dto);

   //  cout << "dtiav=" << dti << " dtoav=" << dto << endl;
   //  pas de temps critique relatif a i
   //  double dt = dti;
   dt = 1.0;
   xcModifyOutOfPlaneCrackVelocity velc_i(lst, lsn, v_iplane3D, cutthru_by_crack_label, dt);

   v_oplane3D.accept(velc_i, narrow_band);
   return dt;
}

void lCrack::createNarrowBand()
{
   if (!do_narrow)
   {
      xAllCreator cc;
      mesh->createSubMesh(narrow_band_label, cc);
      narrow_band = xRegion(mesh, narrow_band_label);
   }
   else
   {
      // create a narrow band
      xAddLayerCreator adl(1, touched_by_crack_label);
      mesh->createSubMesh(narrow_band_label, adl);
      narrow_band = xRegion(mesh, narrow_band_label);
      // reduce front3D and crack3D to a narrow band
      lst.restrictTo(narrow_band, 0.0);
      lsn.restrictTo(narrow_band, 0.0);
      updateForNarrowBand();
   }
}

void lCrack::extendVelocity2()
{
   //  on initialise le champ v_iplane3D a partir du champ v_iplane1D using distance to front
   //   xRegion reg(mesh, touched_by_front_label);
   xRegion reg(mesh, narrow_band_label);
   v_iplane3D.setSupport(narrow_band);
   v_oplane3D.setSupport(narrow_band);
   xVelocityInitialization vel_i(v_iplane1D, reg);
   v_iplane3D.accept(vel_i, narrow_band);

   xVelocityInitialization vel_o(v_oplane1D, reg);
   v_oplane3D.accept(vel_o, narrow_band);
   // v_iplane3D.exportGmsh("vi");
   // v_oplane3D.exportGmsh("vo");
   xexport::Export(v_iplane3D, pexport, "vi");
   xexport::Export(v_oplane3D, pexport, "vo");

   /*
 // Smoothing out
 xL2norm        l2n3;
 double dtvirt = xCFL::getDt(narrow_band) ;
 // xPilotError  statio3(l2n3, Param.tolerance, Param.max_iter_extend, dtvirt*Param.coeff_dt_virt);
 xPilotError  statio3(l2n3, Param.tolerance, 10, dtvirt*Param.coeff_dt_virt);
 std::list<xSpaceOperator*> evo_evi;
 xExtensionOperator extlsn(lsn);
 xExtensionOperator extlst(lst);
 evo_evi.push_back(&extlsn);
 evo_evi.push_back(&extlst);
 xEvolveToStationaryT <xGaussSeidel <xExtensionOperator >  > evol(evo_evi, statio3);
 xFirstMinusSecondCreator cminus(narrow_band_label, touched_by_front_label);
 mesh->createSubsetEntities(all_minus_touched_by_front_label, cminus);
 //xRegion extendreg (mesh, all_minus_touched_by_front_label);
 xRegion extendreg (mesh, narrow_band_label);
 v_iplane3D.accept ( evol, extendreg);
 v_oplane3D.accept ( evol, extendreg );
 */
}

void lCrack::extendVelocity()
{
   v_iplane3D.setSupport(narrow_band);
   v_oplane3D.setSupport(narrow_band);
   v_x3D.setSupport(narrow_band);
   v_y3D.setSupport(narrow_band);
   if (dim() == 3) v_z3D.setSupport(narrow_band);
   v3D.setSupport(narrow_band);

   //  ... initialisations du pilot pour l extension
   //  double tol    = 1.e-6;
   //  int itmax3     = 20;
   xL2norm l2n3;
   double dtvirt = xCFL::getDt(narrow_band);

#ifdef PARALLEL
   {
      double dtvirtMin = 0.;
      MPI_Allreduce(&dtvirt, &dtvirtMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      dtvirt = dtvirtMin;
   }
#endif

   // out << "dtvirt in extend velocity " << dtvirt << endl;
   xPilotError statio3(l2n3, Param.tolerance, Param.max_iter_extend, dtvirt * Param.coeff_dt_virt);

   //  on initialise le champ v_iplane3D a partir du champ v_iplane1D
   // xRegion reg(mesh, touched_by_front_label);
   xRegion reg(mesh, narrow_band_label);
   xVelocityInitialization vel_i(v_iplane1D, reg);

#ifdef PARALLEL
   if (reg.size() == 0) cout << "BEGIN=END ID=" << ParUtil::Instance()->rank() << endl;
#endif

   v_iplane3D.accept(vel_i, narrow_band);

   std::list<xSpaceOperator*> evo_evi;
   xExtensionOperator extlsn(lsn);
   xExtensionOperator extlst(lst);

   evo_evi.push_back(&extlsn);
   evo_evi.push_back(&extlst);

   //  evo_evi.push_back(new xExtensionRadialOperator(lsn,lst));
   xEvolveToStationaryT<xGaussSeidel<xExtensionOperator>> evol(evo_evi, statio3);
   // xEvolveToStationary evol(evo_evi, statio3);

   xFirstMinusSecondCreator cminus(narrow_band_label, touched_by_front_label);
   mesh->createSubMesh(all_minus_touched_by_front_label, cminus);
   // xRegion extendreg (mesh, all_minus_touched_by_front_label);
   xRegion extendreg(mesh, narrow_band_label);
   v_iplane3D.accept(evol, extendreg);

   //  exportViplane3D("v_iplane3D");
   //  v_iplane3D.printDebug();

   if (iter_lsn)
   {
      //  on initialise le champ v_oplane3D a partir du champ v_oplane1D
      xVelocityInitialization vel_o(v_oplane1D, reg);
      v_oplane3D.accept(vel_o);
      // exportVoplane3D("v_oplane3D_b");

      v_oplane3D.accept(evol, extendreg);
      // exportVoplane3D("v_oplane3D");
   }
   else
   {
      xVelocityInitialization vel_x(v_x1D, reg);
      v_x3D.accept(vel_x);
      xVelocityInitialization vel_y(v_y1D, reg);
      v_y3D.accept(vel_y);

      v_x3D.accept(evol, extendreg);
      v_y3D.accept(evol, extendreg);

      if (dim() == 3)
      {
         xVelocityInitialization vel_z(v_z1D, reg);
         v_z3D.accept(vel_z);
         v_z3D.accept(evol, extendreg);
      }

      for (mEntity* pe : narrow_band.range(0))
      {
         mVertex& v = static_cast<mVertex&>(*pe);
         if (dim() == 3)
         {
            xtensor::xVector<> vec(v_x3D(v), v_y3D(v), v_z3D(v));
            v3D(&v) = vec;
         }
         else
         {
            xtensor::xVector<> vec(v_x3D(v), v_y3D(v), 0.);
            v3D(&v) = vec;
         }
      }
      // exportV3D("v_3D");
   }
}

void lCrack::fromV1DToVioplaneOnly()
{
   v_iplane1D.clear();
   v_iplane1D.setSupport(mesh_crack_front);
   v_oplane1D.clear();
   v_oplane1D.setSupport(mesh_crack_front);
   for (mEntity* pe : mesh_crack_front->range(0))
   {
      mVertex* v = static_cast<mVertex*>(pe);
      xtensor::xPoint p = v->point();
      mEntity* e_mesh2d = xMesh::get_const_was_created_by().at(*v);
      mEntity* e_mesh3d = xMesh::get_const_was_created_by().at(*e_mesh2d);
      mEntity* e;
      if (e_mesh3d->getLevel() != mesh->dim())
         e = e_mesh3d->get(mesh->dim(), 0);
      else
         e = e_mesh3d;
      // peut-etre couteux a revoir
      // new tLocalInfos(p, e);
      xfem::xGeomElem elem(e);
      elem.setUVWForXYZ(p);
      setLocalInfos(e, elem.getUVW());
      cout << "V1D " << v1D(v) << endl;
      cout << "no1 " << no1 << endl;
      cout << "no2 " << no2 << endl;
      cout << "Vi " << v1D(v) * no1 << endl;
      cout << "Vo " << v1D(v) * no2 << endl;

      v_iplane1D(v) = v1D(v) * no1;
      v_oplane1D(v) = v1D(v) * no2;
   }
}

void lCrack::fromV1DToVioplane()
{
   //  ... CONSTRUCTION DES CHAMPS SCALAIRES 1D A PARTIR D UN CHAMP VECTORIEL 1D
   //  ... on initialise les champs des vitesses scalaires 1D a parti du champ vectoriel 1D
   cout << "fromV1DToVioplane()" << endl;
   v_iplane1D.clear();
   v_iplane1D.setSupport(mesh_crack_front);
   v_oplane1D.clear();
   v_oplane1D.setSupport(mesh_crack_front);

   v_x1D.clear();
   v_x1D.setSupport(mesh_crack_front);
   v_y1D.clear();
   v_y1D.setSupport(mesh_crack_front);
   v_z1D.clear();
   v_z1D.setSupport(mesh_crack_front);

   for (mEntity* pe : mesh_crack_front->range(0))
   {
      mVertex* v = static_cast<mVertex*>(pe);
      xtensor::xPoint p = v->point();
      mEntity* e_mesh2d = xMesh::get_const_was_created_by().at(*v);
      mEntity* e_mesh3d = xMesh::get_const_was_created_by().at(*e_mesh2d);
      mEntity* e;
      if (e_mesh3d->getLevel() != mesh->dim())
         e = e_mesh3d->get(mesh->dim(), 0);
      else
         e = e_mesh3d;
      // peut-etre couteux a revoir
      // new tLocalInfos(p, e);
      xfem::xGeomElem elem(e);
      elem.setUVWForXYZ(p);
      setLocalInfos(e, elem.getUVW());
      cout << v1D(v) << endl;
      v_iplane1D(v) = v1D(v) * no1;
      v_oplane1D(v) = v1D(v) * no2;
      v_x1D(v) = v1D(v)[0];
      v_y1D(v) = v1D(v)[1];
      v_z1D(v) = v1D(v)[2];

      //    cout << v_iplane1D(v) <<" " <<  v_oplane1D(v) << endl;
   }
}

void lCrack::updateForNarrowBand()
{
   /*
 //definition des domaines
 //
 //reinitialisation du front  (to improve front1d should be forzen)
 xL2norm       l2n;
 double dtvirt = xCFL::getDt(narrow_band);
 xPilotError  statio(l2n, Param.tolerance, Param.max_iter_reinit, dtvirt*Param.coeff_dt_virt);

 xReInitOperator         reinit;
 xEvolveToStationary rei(reinit, statio);
 lst.accept ( rei );

 //
 xL2norm        l2n2;
 xPilotError  statio2(l2n2, Param.tolerance, Param.max_iter_ortho, dtvirt*Param.coeff_dt_virt);
 xOrthoOperator  ortho_wrt_front(lst);
 xEvolveToStationary orthof(ortho_wrt_front, statio2);

 //anthony
 xL2norm        l2n3;
 xPilotError  statio3(l2n3, Param.tolerance, Param.max_iter_extend, dtvirt*Param.coeff_dt_virt);
 xExtensionOperator  extend_wrt_front(lst);
 xEvolveToStationary extendf(extend_wrt_front, statio3);

 //anthony
 cNegativeSideCreator cneg(lst);
 mesh->createSubsetEntities("supports_negative_front", cneg);
 cFirstMinusSecondCreator cminus("narrow_band", "supports_negative_front");
 mesh->createSubsetEntities("all_minus_negative_front", cminus);
 //  crack3D.accept ( orthof , lSupport(mesh, "all_minus_negative_front"));
 lsn.accept ( extendf , xRegion(mesh, "all_minus_negative_front"));
 //  crack3D.accept ( orthof );

 //
 lsn.accept ( rei );

 //
 xOrthoOperator  ortho_wrt_crack(lsn);
 xEvolveToStationary orthoc(ortho_wrt_crack, statio2);
 lst.accept ( orthoc );
 */
}

void lCrack::getLocalCurv(mEntity* e, const xtensor::xPoint& uvw, xtensor::xVector<>& e1, xtensor::xVector<>& e2,
                          xtensor::xVector<>& e3, xtensor::xTensor2<>& curv1, xtensor::xTensor2<>& curv2,
                          xtensor::xTensor2<>& curv3) const
{
   const bool debug = false;
   curv1 = lst.getCurv(e);
   curv2 = lsn.getCurv(e);
   // curv3 is the gradient of  e1 x e2
   // curv3[.][0] is the gradient of e1 x e2 with respect to zero
   // curv3[.][0] =  e1,0 x e2 + e1 x e2,0 = curv1[.][0] x e2 + e1 x curv2[.][0]
   // so curv3[.][j] is
   // curv3[.][j] = curv1[.][j] x e2 + e1 x curv2[.][j]
   // curv3[.][j] = temp1 x e2 + e1 x temp2
   e1 = lst.getGrad(e, uvw);
   e2 = lsn.getGrad(e, uvw);
   e3 = e1 % e2;
   xtensor::xVector<> temp1, temp2;
   for (int j = 0; j < 3; j++)
   {
      for (int i = 0; i < 3; i++)
      {
         temp1(i) = curv1(i, j);
         temp2(i) = curv2(i, j);
      }
      xtensor::xVector<> res = temp1 % e2 + e1 % temp2;
      for (int i = 0; i < 3; i++) curv3(i, j) = res(i);
   }
   if (debug)
   {
      cout << "results in getLocalCurv" << endl;
      cout << e1 << " e1 " << endl;
      cout << " e2 " << e2 << endl;
      cout << " e3 " << e3 << endl;
      cout << " curv1 " << curv1 << endl;
      cout << " curv2 " << curv2 << endl;
      cout << " curv3 " << curv3 << endl;
   }
}

void lCrack::getLocalSmoothOrthoAxis(mEntity* e, const xtensor::xPoint& uvw, xtensor::xPoint& local_coords,
                                     xtensor::xVector<>& eo1, xtensor::xVector<>& eo2, xtensor::xVector<>& eo3) const
{
   const bool debug = false;
   xElement elem(e);
   elem.setUvw(uvw);
   std::vector<double> c = lsn.getVals(e);
   std::vector<double> f = lst.getVals(e);
   local_coords(0) = elem.getInterpoSca(f);
   local_coords(1) = elem.getInterpoSca(c);
   local_coords(2) = 0.;  // this value is in fact never used and
   xtensor::xVector<> e1, e2;
   // new p -> uvw
   eo1 = lst.getGrad(e, uvw);
   eo1.norm();
   e2 = lsn.getGrad(e, uvw);
   eo2 = e2 - eo1 * (e2 * eo1);
   eo2.norm();
   eo3 = eo1 % eo2;

   if (debug)
   {
      cout << "results in getLocalSmoothOrthoAxis" << endl;
      cout << " local_coords " << local_coords << endl;
      cout << " eo1 " << eo1 << endl;
      cout << " eo2 " << eo2 << endl;
      cout << " eo3 " << eo3 << endl;
   }

#if DEBUGNIC >= 1
   printf("axis smooth ortho at point %12.5e %12.5e %12.5e are\n", p(0), p(1), p(2));
   printf("%12.5e %12.5e %12.5e\n", eo1(0), eo1(1), eo1(2));
   printf("%12.5e %12.5e %12.5e\n", eo2(0), eo2(1), eo2(2));
   printf("%12.5e %12.5e %12.5e\n", eo3(0), eo3(1), eo3(2));
   printf("local is %12.5e %12.5e %12.5e are\n", local(0), local(1), local(2));
#endif
   return;
}

void lCrack::getLocalSmoothOrthoAxis(mEntity* e, xtensor::xVector<>& eo1, xtensor::xVector<>& eo2, xtensor::xVector<>& eo3) const
{
   xtensor::xVector<> e1, e2;
   eo1 = lst.getGrad(e);
   eo1.norm();
   e2 = lsn.getGrad(e);
   eo2 = e2 - eo1 * (e2 * eo1);
   eo2.norm();
   eo3 = eo1 % eo2;
}

void lCrack::setLocalInfos(mEntity* e, const xtensor::xPoint& uvw) const
{
   const bool debug = false;
   //
   // Computation of all the local quantities
   //
   xElement elem(e);
   // new beg elem.xyz2uvw(p);
   elem.setUvw(uvw);
   // new end
   std::vector<double> f = lst.getVals(e);
   std::vector<double> c = lsn.getVals(e);
   local(0) = elem.getInterpoSca(f);
   local(1) = elem.getInterpoSca(c);
   local(2) = 0.;
   // new p -> uvw
   n1 = lst.getGrad(e, uvw);
   n1_norm = std::sqrt(n1 * n1);
   no1 = n1 / n1_norm;
   n2 = lsn.getGrad(e, uvw);
   n2_norm = std::sqrt(n2 * n2);
   no2 = n2 / n2_norm;
   n3 = (n1 % n2);
   n3_norm = std::sqrt(n3 * n3);
   no3 = n3 / n3_norm;

   for (int i = 0; i < 3; i++)
   {
      n123_row(0, i) = n1(i);
      n123_row(1, i) = n2(i);
      n123_row(2, i) = n3(i);
      n123_col(i, 0) = n1(i);
      n123_col(i, 1) = n2(i);
      n123_col(i, 2) = n3(i);
      no123_row(0, i) = no1(i);
      no123_row(1, i) = no2(i);
      no123_row(2, i) = no3(i);
      no123_col(i, 0) = no1(i);
      no123_col(i, 1) = no2(i);
      no123_col(i, 2) = no3(i);
   }

   metricn = n123_row * n123_col;

   curv1 = lst.getCurv(e);
   curv2 = lsn.getCurv(e);
   // curv3 is the gradient of  n1 x n2;
   // first we compute the derivatives of n1 X n2
   // curv3[.][0] is the gradient of n1 x n2 with respect to zero
   // curv3[.][0] =  n1,0 x n2 + n1 x n2,0 = curv1[.][0] x n2 + n1 x curv2[.][0]
   // so curv3[.][j] is
   // curv3[.][j] = curv1[.][j] x n2 + n1 x curv2[.][j]
   // curv3[.][j] = temp1 x n2 + n1 x temp2
   xtensor::xVector<> temp1, temp2;
   for (int j = 0; j < 3; j++)
   {
      for (int i = 0; i < 3; i++)
      {
         temp1(i) = curv1(i, j);
         temp2(i) = curv2(i, j);
      }
      xtensor::xVector<> res = (temp1 % n2) + (n1 % temp2);
      for (int i = 0; i < 3; i++) curv3(i, j) = res(i);
   }

   // symmetrizing??
   //  curv1.symmetrize(); curv2.symmetrize(); curv3.symmetrize();

   // derivatives of the normalized vectors
   temp1 = no1 * curv1;
   temp2 = no2 * curv2;
   for (int j = 0; j < 3; j++)
   {
      for (int i = 0; i < 3; i++)
      {
         curvo1(i, j) = (curv1(i, j) - no1(i) * temp1(j)) / n1_norm;
         curvo2(i, j) = (curv2(i, j) - no2(i) * temp2(j)) / n2_norm;
      }
   }
   for (int j = 0; j < 3; j++)
   {
      for (int i = 0; i < 3; i++)
      {
         temp1(i) = curvo1(i, j);
         temp2(i) = curvo2(i, j);
      }
      xtensor::xVector<> res = (temp1 % no2) + (no1 % temp2);
      for (int i = 0; i < 3; i++) curvo3(i, j) = res(i);
   }

   if (debug)
   {
      cout << "results of setLocalInfos" << endl;
      cout << "local  position " << local << endl;
      cout << "n1 " << n1 << " n2 " << n2 << " n3 " << n3 << endl;
      cout << "n123_row" << endl << n123_row << endl;
      cout << "n123_col" << endl << n123_col << endl;
      cout << "no123_row" << endl << no123_row << endl;
      cout << "no123_col" << endl << no123_col << endl;

      // checks
      cout << "metric n dot n" << endl << metricn << endl;
      // end checks

      cout << "curv1 " << endl << curv1 << endl;
      cout << "curv2 " << endl << curv2 << endl;
      cout << "curv3 " << endl << curv3 << endl;
      cout << "curvo1 " << endl << curvo1 << endl;
      cout << "curvo2 " << endl << curvo2 << endl;
      cout << "curvo3 " << endl << curvo3 << endl;
   }

   return;
}

void lCrack::getLocalCoords(mEntity* e, xtensor::xPoint& uvw, double& r, double& theta) const
{
   xElement elem(e);
   elem.setUvw(uvw);
   std::vector<double> valst = lst.getVals(e);
   std::vector<double> valsn = lsn.getVals(e);
   /*  double Vn=elem.getInterpoSca(valst);
   double Vt=elem.getInterpoSca(valsn);
   //    printf(" n %f t %f ",Vn,Vt);
   double norm;
   xtensor::xVector<> Gn=lst.getGrad(e);
   xtensor::xVector<> Gt=lsn.getGrad(e);

   norm=sqrt(Gt*Gt);
   xtensor::xVector<> Gnormt=Gt/norm;
   xtensor::xVector<> Ft=(Gn%Gt);
   //    norm=sqrt(Ft*Ft);
   //    mVector Fnormt=Ft/norm;
   xtensor::xVector<> Fn=Gt%Ft;
   norm=sqrt(Fn*Fn);
   xtensor::xVector<> Fnormn=Fn/norm;

   double x=Vn/(Gn*Fnormn);
   double y=Vt/(Gt*Gnormt);
*/
   double x = elem.getInterpoSca(valst);
   double y = elem.getInterpoSca(valsn);

   //    printf(" x %f y %f ",x,y);
   theta = atan2(y, x);
   r = sqrt(x * x + y * y);
   //    printf(" r %f th %f \n",r,theta);
}

std::vector<double> lCrack::getLstValues(mEntity* e) const { return lst.getVals(e); }

std::vector<double> lCrack::getLsnValues(mEntity* e) const { return lsn.getVals(e); }

// on which side of the crack is made the interpolation
int lCrack::sideOf(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ) const
{
   if (overrideside) return overrideside;
   mEntity* e = geo_appro->getEntity();
   xElement elem(e);
   elem.setUvw(geo_appro->getUVW());
   std::vector<double> vals = getFieldn()->getVals(e);
   double f = elem.getInterpoSca(vals);
   if (fabs(f) < 1e-3)
   {
      if (geo_integ != nullptr)
      {
         elem.xyz2uvw(geo_integ->getCDGxyz());
      }
      else
      {
         elem.setUvw(geo_appro->getCDGuvw());
      }
      std::vector<double> vals = getFieldn()->getVals(e);
      f = elem.getInterpoSca(vals);
   }
   return f >= 0 ? 1 : -1;
}

void lCrack::cutMeshAndSlice()
{
   const bool debug = false;
   if (mesh_crack_surface) delete mesh_crack_surface;
   if (mesh_crack_front) delete mesh_crack_front;
   mesh_crack_surface = nullptr;
   mesh_crack_front = nullptr;

   mesh_crack_surface = new xMesh;
   mesh_crack_front = new xMesh;

   // the following create the crack surface mesh
   if (debug) cout << " creating mesh_crack_surface in CutMeshAndSlice " << endl;
   // void xMesh::cutMesh(const xLevelSet& ls, xMesh* x_interface, xEntityToEntity classify_in,
   //              EntityToEntity classify_out, bool create_partition, bool keep_old_partition, bool recursive)
   // mesh->cutMesh(lsn, mesh_crack_surface, xtool::xIdentity<mEntity*>(), xtool::xIdentity<mEntity*>(),  true, true);

   xMesh::datamanager_t<AOMD::mEntity*>& was_created_by = xfem::xMesh::get_was_created_by();
   // xMesh::datamanager_t< double > &r_on_edge           = xfem::xMesh::get_r_on_edge();
   xMesh::datamanager_t<double> r_on_edge;
   xMesh::datamanager_t<AOMD::mEntity*>& is_duplicated_in = xfem::xMesh::get_is_duplicated_in();
   xMesh::datamanager_t<AOMD::mEntity*>& was_duplicated_from = xfem::xMesh::get_was_duplicated_from();
   xMesh::datamanager_t<xMesh>& partition = xfem::xMesh::get_partition();
   xMesh::datamanager_t<AOMD::mEntity*>& is_in_partition_of = xfem::xMesh::get_is_in_partition_of();

   xcut::cutMesh(*mesh, lsn, *mesh_crack_surface, was_created_by, r_on_edge, is_duplicated_in, was_duplicated_from, partition,
                 is_in_partition_of, xtool::xIdentity<mEntity*>(), xtool::xIdentity<mEntity*>(), true, true, false);
   // on the following are variants that might be usefull
   // when the ls is not defined on the full mesh :
   //  mesh->cutMesh(lsn,mesh_crack_surface, xAcceptInsideLevelSet(lst),
   //                xtool::xIdentity<mEntity*>(), xtool::xIdentity<mEntity*>(),  true, true);
   //  mesh->cutMesh(lsn,mesh_crack_surface, xAcceptSupport(lst,mesh));
   if (debug) cout << " done creating mesh_crack_surface in CutMeshAndSlice " << endl;
   lst.takeTraceOn(*mesh_crack_surface, was_created_by, r_on_edge, lst_crack_surface);
   if (fit)
   {
      xFitToVertices fit(Param.tolfit);
      lst_crack_surface.accept(fit);
   }

   // the following create the crack front mesh
   if (debug) cout << " creating mesh_crack_front in CutMeshAndSlice " << endl;
   // mesh_crack_surface->cutMesh(lst_crack_surface,mesh_crack_front,
   //                            xtool::xIdentity<mEntity*>(),xtool::xIdentity<mEntity*>(),false);
   //
   xcut::createIsoZeroMeshFromLevelSet(lst_crack_surface, *mesh_crack_front, was_created_by, r_on_edge);

   return;
}

bool lCrack::supportInCylinderAroundFront(double r, AOMD::mEntity* e) const
{
   std::vector<AOMD::mEntity*> pvs;
   if (e->getType() == 0)
      pvs.push_back(e);
   else
   {
      int nnode = e->size(0);
      for (int i = 0; i < nnode; ++i) pvs.push_back(e->get(0, i));
   }
   for (auto pv : pvs)
   {
      double x, y;
      try
      {
         x = lst(pv);
      }
      catch (...)
      {
         x = r + 1.;
      }
      try
      {
         y = lsn(pv);
      }
      catch (...)
      {
         y = r + 1.;
      }
      if ((x * x + y * y) <= (r * r)) return true;
   }
   return false;
}

void lCrack::createSupportInfo()
{
   // build the mesh subspace for
   //  touched by front "supports_touched_by_crack"
   //  touched by crack "supports_touched_by_front"
   //  cut by crack     "supports_cut_by_crack"
   //  cuthru by crack  "supports_cuthru_by_crack"
   //
   // the results will be independent
   // of the narrow band choice provided it is big enough!!
   //
   mesh->deleteSubMesh(touched_by_crack_label);
   mesh->deleteSubMesh(touched_by_front_label);
   mesh->deleteSubMesh(cut_by_crack_label);
   mesh->deleteSubMesh(cutthru_by_crack_label);
   mesh->deleteSubMesh(integrate_front_label);
   xUnifySubMeshAccrossProcess submesunifyer;

   xcSupportsTouchedByCrackCreator c1(lst_crack_surface);
   mesh->createSubMesh(touched_by_crack_label, c1);
   submesunifyer.modify(*mesh, touched_by_crack_label);

   xcSupportsTouchedByFrontCreator c2(mesh_crack_front);
   mesh->createSubMesh(touched_by_front_label, c2);
   submesunifyer.modify(*mesh, touched_by_front_label);

   xcSupportsCutByCrackCreator c3(lst_crack_surface);
   auto& cut_by_crack = mesh->createSubMesh(cut_by_crack_label, c3);
   submesunifyer.modify(*mesh, cut_by_crack_label);

   xcSupportsCutByCrackToDeleteCreator c3bis(lst_crack_surface);
   auto& cut_by_crack_to_delete = mesh->createSubMesh(cut_by_crack_to_delete_label, c3bis);
   submesunifyer.modify(*mesh, cut_by_crack_to_delete_label);

   xSubMesh& cutthru_by_crack = mesh->createSubMesh(cutthru_by_crack_label);
   for (int d = 0; d <= 3; d++)
      for (mEntity* pe : cut_by_crack.range(d)) cutthru_by_crack.add(pe);

   // and then delete in cutthru_by_crack entities in cut_by_crack_to_delete
   for (int d = 0; d <= 3; d++)
      for (mEntity* pe : cut_by_crack_to_delete.range(d)) cutthru_by_crack.del(pe);

   xcSupportsIntegrateCreator c5(touched_by_front_label);
   mesh->createSubMesh(integrate_front_label, c5);
   return;
}

bool lCrack::supports_touched_by_crack(mEntity* e) const { return (mesh->getSubMesh(touched_by_crack_label).find(e)); }
bool lCrack::supports_touched_by_front(mEntity* e) const { return (mesh->getSubMesh(touched_by_front_label).find(e)); }
bool lCrack::supports_cut_by_crack(mEntity* e) const { return (mesh->getSubMesh(cut_by_crack_label).find(e)); }
bool lCrack::supportsCutthruByCrack(mEntity* e) const { return (mesh->getSubMesh(cutthru_by_crack_label).find(e)); }

xEntityFilter lCrack::supportsCutthruByCrackFilter() const
{
   return bind1st(mem_fun(&xcrack::lCrack::supportsCutthruByCrack), this);
}

bool lCrack::supports_integrate_front(mEntity* e) const { return (mesh->getSubMesh(integrate_front_label).find(e)); }

bool lCrack::side_of_crack_surface(mEntity* e) const { return (lsn(e) > 0); }

void lCrack::getCrackAxis(mEntity* e, xtensor::xTensor2<>& base) const
{
   const xtensor::xVector<> g1 = lst.getGrad(e);
   const xtensor::xVector<> g2 = lsn.getGrad(e);
   const xtensor::xVector<> g3 = g1 % g2;
   base = xtensor::xTensor2<>(g1, g2, g3);
}

void lCrack::getCrackOrthoAxis(mEntity* e, xtensor::xTensor2<>& base) const
{
   const xtensor::xVector<> eo1 = (lst.getGrad(e)).norm();
   xtensor::xVector<> eo2 = lsn.getGrad(e);
   eo2 = (eo2 - eo1 * (eo2 * eo1)).norm();
   const xtensor::xVector<> eo3 = eo1 % eo2;
   base = xtensor::xTensor2<>(eo1, eo2, eo3);
}

void lCrack::localToGlobal(mEntity* e, const xtensor::xPoint& uvw, xtensor::xVector<>& local) const
{
   const bool debug = false;
   xtensor::xVector<> e1 = lst.getGrad(e, uvw);
   xtensor::xVector<> e2 = lsn.getGrad(e, uvw);
   xtensor::xVector<> e3 = e1 % e2;
   if (debug)
   {
      cout << "in localtoGlobal" << endl;
      cout << " local " << local << endl;
      cout << " e1 " << e1 << endl;
      cout << " e2 " << e2 << endl;
      cout << " e3 " << e3 << endl;
   }
   xtensor::xVector<> temp = local;
   local(0) = temp(0) * e1(0) + temp(1) * e2(0) + temp(2) * e3(0);
   local(1) = temp(0) * e1(1) + temp(1) * e2(1) + temp(2) * e3(1);
   local(2) = temp(0) * e1(2) + temp(1) * e2(2) + temp(2) * e3(2);
   if (debug)
   {
      cout << " global " << local << endl;
   }
   return;
}

void lCrack::localToGlobal(mEntity* e, const xtensor::xPoint& uvw, xtensor::xTensor2<>& local) const
{
   xtensor::xVector<> e1 = lst.getGrad(e, uvw);
   xtensor::xVector<> e2 = lsn.getGrad(e, uvw);
   xtensor::xVector<> e3 = e1 % e2;
   xtensor::xTensor2<> temp(local);
   xtensor::xTensor2<> Q(e1, e2, e3);
   Q.transpose();
   temp = Q * temp;
   xtensor::xTensor2<> temp2 = Q.invert();
   local = temp * temp2;
   return;
}

xcSupportsTouchedByCrackCreator::xcSupportsTouchedByCrackCreator(const xLevelSet& c) : xSubMeshCreator(), crack2D(c) {}

void xcSupportsTouchedByCrackCreator::create(const xMesh& m, const string& name)
{
   // first find the element touching
   int dim = m.dim();
   xSubMesh& sub = m.getSubMesh(name);
   xRegion& slice = const_cast<xRegion&>(crack2D.getSupport());
   for (mEntity* pe : slice.range(0))
   {
      mVertex* v = static_cast<mVertex*>(pe);
      if (crack2D(v) <= 0.0)
      {
         mEntity* up = xMesh::get_const_was_created_by().at(*v);
         sub.add(up);
         for (int d = up->getLevel() + 1; d <= dim; d++)
         {
            for (int j = 0; j < up->size(d); j++)
            {
               mEntity* e = up->get(d, j);
               sub.add(e);
            }
         }
      }
   }

   std::set<mEntity*> to_insert;
   for (mEntity* e : sub.range(dim))
   {
      for (int d = 0; d < dim; d++)
      {
         for (int j = 0; j < e->size(d); j++)
         {
            mEntity* down = e->get(d, j);
            to_insert.insert(down);
         }
      }
   }
   for (mEntity* pe : to_insert) sub.add(pe);
}

xcSupportsTouchedByFrontCreator::xcSupportsTouchedByFrontCreator(xMesh* m1D) : xSubMeshCreator(), mesh_front1D(m1D) {}

void xcSupportsTouchedByFrontCreator::create(const xMesh& m, const string& name)
{
   // first find the element touching
   int dim = m.dim();
   xSubMesh& sub = m.getSubMesh(name);

   for (int d = 0; d <= 1; d++)
   {
      for (mEntity* e1d : mesh_front1D->range(d))
      {
         mEntity* e2d = xMesh::get_const_was_created_by().at(*e1d);
         mEntity* e3d = xMesh::get_const_was_created_by().at(*e2d);
         if (e3d->getLevel() == dim)
            sub.add(e3d);
         else
         {
            for (int j = 0; j < e3d->size(dim); j++)
            {
               mEntity* e = e3d->get(dim, j);
               sub.add(e);
            }
         }
      }
   }

   std::set<mEntity*> to_insert;
   for (mEntity* e : sub.range(dim))
   {
      for (int d = 0; d < dim; d++)
      {
         for (int j = 0; j < e->size(d); j++)
         {
            mEntity* down = e->get(d, j);
            to_insert.insert(down);
         }
      }
   }
   for (mEntity* pe : to_insert) sub.add(pe);
}

xcSupportsIntegrateCreator::xcSupportsIntegrateCreator(const string& name_tbf_) : xSubMeshCreator(), name_tbf(name_tbf_) {}

void xcSupportsIntegrateCreator::create(const xMesh& m, const string& name)
{
   int dim = m.dim();
   xSubMesh& sub = m.getSubMesh(name);
   xSubMesh& sub_tbf = m.getSubMesh(name_tbf);
   std::set<mEntity*> to_insert;

   for (int d = dim; d >= 0; d--)
   {
      for (mEntity* e : sub_tbf.range(d))
      {
         sub.add(e);
         to_insert.insert(e);
      }
   }

   for (mEntity* e : sub.range(dim))
   {
      for (int d = 0; d < dim; d++)
      {
         for (int j = 0; j < e->size(d); j++)
         {
            mEntity* down = e->get(d, j);
            if (to_insert.find(down) == to_insert.end())
            {
               to_insert.insert(down);
               sub.add(down);
            }
         }
      }
   }

   for (int d = 0; d < dim; d++)
   {
      for (mEntity* e : sub.range(d))
      {
         for (int k = dim; k > d; k--)
         {
            for (int j = 0; j < e->size(k); j++)
            {
               mEntity* up = e->get(k, j);
               if (to_insert.find(up) == to_insert.end())
               {
                  to_insert.insert(up);
                  sub.add(up);
               }
            }
         }
      }
   }
}

xcSupportsCutByCrackCreator::xcSupportsCutByCrackCreator(const xLevelSet& c) : xSubMeshCreator(), crack2D(c) {}

void xcSupportsCutByCrackCreator::create(const xMesh& m, const string& name)
{
   // first find the element touching
   xSubMesh& sub = m.getSubMesh(name);
   xRegion slice = crack2D.getSupport();
   int dim1 = slice.dim();
   for (int dim = 0; dim <= dim1; dim++)
   {
      for (mEntity* e : slice.range(dim))
      {
         double fmin;
         if (dim > 0)
         {
            std::vector<double> f = crack2D.getVals(e);
            fmin = *std::min_element(f.begin(), f.end());
         }
         else
         {
            mVertex* v = (mVertex*)e;
            fmin = crack2D(v);
         }

         if (fmin < 0.0)
         {
            mEntity* up = xMesh::get_const_was_created_by().at(*e);
            // note that up may be a face or an element
            // for an element we add it directly
            sub.add(up);
         }
      }
   }
   // tous les enfants sont a ajouter
   sub.modifyAllState();
}

// the idea is to start with the cut set
// and to discard the ones that are not cuthru

// loop over the elements of crack2D
//  this will give all the element cut but not cuthru
//                 take all the subentities of these not cuthru
//  and also give all the faces cut but not cuthru
//                 take all the subentities of these not cuthru

// remove from the cut set all the ones that are not cuthru

xcSupportsCutByCrackToDeleteCreator::xcSupportsCutByCrackToDeleteCreator(const xLevelSet& c) : crack2D(c) {}

void xcSupportsCutByCrackToDeleteCreator::create(const xMesh& m, const string& name)
{
   const bool debug = false;
   // find the element touching
   xRegion slice = crack2D.getSupport();
   xSubMesh& sub = m.getSubMesh(name);
   for (mEntity* e : slice.range())
   {
      std::vector<double> f = crack2D.getVals(e);
      if (debug)
      {
         cout << " crack2D is" << endl;
         for (xIter itd = slice.begin(); itd != slice.end(); itd++)
         {
            mEntity* ed = *itd;
            ed->print();
         }
         cout << " crack2D.getVals(e)" << endl;
         e->print();
         for (int i0 = 0; i0 < e->size(0); i0++)
         {
            mEntity* eg0 = e->get(0, i0);
            e->get(0, i0)->print();
            std::vector<double> fgO = crack2D.getVals(eg0);
            cout << "std::vector<double> fgO = crack2D.getVals(eg0); taille" << endl;
            cout << fgO.size() << endl;
            // mVertex * veg0 = (mVertex *)eg0;
            cout << "crack2D find(eg0)" << endl;
            // cout<<crack2D.find(veg0)->second<<endl;
         }
         cout << "f.size() " << f.size() << endl;
         for (size_t i = 0; i < f.size(); i++)
         {
            cout << "value of f(" << i << ") = ";
            cout << f[i] << endl;
         }
      }

      double fmax = *std::max_element(f.begin(), f.end());
      if (fmax > 0.0)
      {
         mEntity* up = xMesh::get_const_was_created_by().at(*e);
         sub.add(up);  // note that up may be a face or an element
         for (int d = up->getLevel() - 1; d >= 0; d--)
         {
            for (int j = 0; j < up->size(d); j++)
            {
               mEntity* e = up->get(d, j);
               sub.add(e);
            }
         }
      }
   }
}

// the idea is to start with the cut set
// and to discard the ones that are not cuthru
xcSupportsCutthruByCrackCreator::xcSupportsCutthruByCrackCreator(const xLevelSet& c, const string& ncut, const string& ndel)
    : xSubMeshCreator(), crack2D(c), name_cut(ncut), name_del(ndel)
{
}

void xcSupportsCutthruByCrackCreator::create(const xMesh& m, const string& name)
{
   xSubMesh& sub = m.getSubMesh(name);
   xSubMesh& sub_cut = m.getSubMesh(name_cut);
   xSubMesh& sub_del = m.getSubMesh(name_del);
   // first copy name_cut into name
   for (int d = 0; d <= 3; d++)
      for (mEntity* pe : sub_cut.range(d)) sub.add(pe);
   // and then delete name_del
   for (int d = 0; d <= 3; d++)
      for (mEntity* pe : sub_del.range(d)) sub.del(pe);
}

xcStrictNegativeSideCreator::xcStrictNegativeSideCreator(const xLevelSet& c) : xSubMeshCreator(), field(c) {}

void xcStrictNegativeSideCreator::create(const xMesh& m, const string& name)
{
   xRegion slice = field.getSupport();
   xSubMesh& sub = m.getSubMesh(name);
   for (mEntity* e : slice)
   {
      std::vector<double> f = field.getVals(e);
      double fmax = *std::max_element(f.begin(), f.end());
      if (fmax <= 0.0) sub.add(e);
   }
}

xcNegativeSideCreator::xcNegativeSideCreator(const xLevelSet& c) : xSubMeshCreator(), field(c) {}

void xcNegativeSideCreator::create(const xMesh& m, const string& name)
{
   xRegion slice = field.getSupport();
   xSubMesh& sub = m.getSubMesh(name);
   for (mEntity* e : slice.range(0))
   {
      if (field(e) <= 0.0) sub.add(e);
   }

   for (mEntity* e : slice)
   {
      std::vector<double> f = field.getVals(e);
      double fmin = *std::min_element(f.begin(), f.end());
      if (fmin <= 0.0) sub.add(e);
   }
}

xcElementsAlongCrackFrontCreator::xcElementsAlongCrackFrontCreator(const lCrack& _c) : xSubMeshCreator(), c(_c) {}

void xcElementsAlongCrackFrontCreator::create(const xMesh& m, const string& name)
{
   // assert(m->dim() == 3);
   int dim = m.dim();
   const xLevelSet& lsnm = *c.getFieldn();
   const xLevelSet& lstm = *c.getFieldt();
   xSubMesh& sub_touched_by_front = m.getSubMesh(c.touched_by_front_label);
   xSubMesh& sub = m.getSubMesh(name);
   for (mEntity* e : sub_touched_by_front.range(dim))
   {
      mVertex* vmin = nullptr;
      double rmin = 1.;  // fake value to avoid -Wmaybe-uninitialized below. Anyway when j= 0 rmin is set with r ...
      for (int j = 0; j < e->size(0); j++)
      {
         mVertex* v = (mVertex*)e->get(0, j);
         double a = lsnm(v);
         double b = lstm(v);
         double r = sqrt(a * a + b * b);
         if (j == 0 || r < rmin)
         {
            vmin = v;
            rmin = r;
         }
      }
      for (int j = 0; j < vmin->size(dim); j++)
      {
         mEntity* ee = vmin->get(dim, j);
         sub.add(ee);
      }
   }
   sub.modifyAllState();
}

void xcModifyOutOfPlaneCrackVelocity::visit(xLevelSet& v_o, xRegion target)
{
   // xField valpos;
   // xRegion r;

   //  ... initialisation de la valeur positive de psi (constante pour tout le calcul)
   //  for (xIter it=target.begin(0); it!=target.end(0); ++it)
   //  {
   //  ... on va chercher le numero du noeud courant
   //    mVertex *v = (mVertex *) *it;
   //    if (front(v)<=0.) {valpos(v) = 0.;}
   //    else {valpos(v) = 1.;}
   //     {valpos(v) = 1.;}
   //  }
   /*
 for (int d = 2; d <= 3; d++)
 {
   for (xIter it=
 v_o.getSupport().getMesh()->begin_sub(d,supports_cuthru_by_crack);it!=v_o.getSupport().getMesh()->end_sub(d,supports_cuthru_by_crack);
 ++it)
   {
     mEntity *e =  *it;
     for (int j = 0; j < e->size(0); j++)
     {
   mEntity* v = e->get(0,j);
       for (int k = 0; k < v->size(3); k++)
   {
         mEntity* e2 = v->get(3, k);
     for (int l = 0; l < e2->size(0); l++)
     {
       mEntity* v2 = e2->get(0, l);
       valpos(v2) = 0.0;
     }
   }
     }
   }
 }*/
   //  ... reste finalement a multiplier le champ des vitesses par
   //        (1./vpsi/dt*psi*valeur_positive(psi)) ................................

   for (mEntity* pe : target.range(0))
   {
      //  ... on va chercher le numero du noeud courant
      mVertex* v = static_cast<mVertex*>(pe);
      double p = plane(v);
      double f = front(v);
      double vi = v_i(v);
      double vo = v_o(v);
      double th = atan2(p, f);
      double alpha = atan2(vo, vi) / 2;
      alpha = 0;
      double valpos;
      if ((f + p * tan(alpha)) <= 0.)
      {
         valpos = 0.;
      }
      else
      {
         valpos = 1.;
      }
      v_o(v) = vo * ((f + p * tan(alpha)) * valpos / vi / dt) * cos(th);
   }
}

void lCrack::setV1Ddebug(xtensor::xVector<> Vec)
{
   v1D.setSupport(mesh_crack_front);
   v1D.clear();
   int dimf = mesh_crack_front->dim();
   for (mEntity* pe : mesh_crack_front->range(dimf)) v1D(pe) = Vec;
   for (mEntity* pe : mesh_crack_front->range(0)) v1D(static_cast<mVertex*>(pe)) = Vec;
   // fromV1DToVioplane();
   fromV1DToVioplaneOnly();
}

void lCrack::setVioplane3Ddebug(double iplane, double oplane)
{
   v_iplane3D.clear();
   v_oplane3D.clear();

   v_iplane3D.setSupport(narrow_band);
   v_oplane3D.setSupport(narrow_band);
   for (mEntity* pe : narrow_band.range(0))
   {
      mVertex* v = static_cast<mVertex*>(pe);
      v_iplane3D(v) = iplane;
      v_oplane3D(v) = oplane;
   }
}

void lCrack::exportV3D(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;
   // Export(v3D, pexport, s);
   v3D.exportGmsh(s);
}

void lCrack::exportVoplane3D(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;
   // v_oplane3D.exportGmsh(s);
   Export(v_oplane3D, pexport, s);
}

void lCrack::exportViplane3D(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;
   // v_iplane3D.exportGmsh(s);
   xexport::Export(v_iplane3D, pexport, s);
}

void lCrack::exportVoplane1D(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;
   // v_oplane1D.exportGmsh(s);
   xexport::Export(v_oplane1D, pexport, s);
}

void lCrack::exportViplane1D(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;
   // v_iplane1D.exportGmsh(s);
   xexport::Export(v_iplane1D, pexport, s);
}

void lCrack::exportlsn(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;
   // lsn.exportGmsh(s);
   xexport::Export(lsn, pexport, s);
}

void lCrack::exportlst(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;
   // lst.exportGmsh(s);
   xexport::Export(lst, pexport, s);
}

void lCrack::exportmeshcracksurface(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;

#ifdef PARALLEL
   stringstream oufilename;
   oufilename << "_" << ParUtil::Instance()->size() << "_" << ParUtil::Instance()->rank() + 1;
   s = s + oufilename.str() + ".msh";
#endif

#ifndef PARALLEL
   s = s + ".msh";
#endif

   AOMD_Util::Instance()->ex_port(s.c_str(), &mesh_crack_surface->getMesh());
}

void lCrack::exportmeshcracktip(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;

#ifdef PARALLEL
   stringstream oufilename;
   oufilename << "_" << ParUtil::Instance()->size() << "_" << ParUtil::Instance()->rank() + 1;
   s = s + oufilename.str() + ".msh";
#endif

#ifndef PARALLEL
   s = s + ".msh";
#endif

   AOMD_Util::Instance()->ex_port(s.c_str(), &mesh_crack_front->getMesh());
}

void lCrack::exportlstcracksurface(string s)
{
   char s2[10];
   sprintf(s2, "-%d", timestep);
   s = s + s2;
   // lst_crack_surface.exportGmsh(s);
   xexport::Export(lst_crack_surface, pexport, s);
}

int lCrack::dim() const { return mesh->dim(); }
void lCrack::setLssOnFront(double mins, double maxs)
{
   const bool debug = true;
   // assumptions : there is only one front.
   // it is not a closed curve.
   // int created_by_tag=mesh->get_was_created_by_tag();
   bool isloop = false;
   xinterface::aomd::xAttachedDataManagerAOMD<int> visited;
   double lenght = 0.;
   if (dim() == 3)
   {
      mesh_crack_front->getMesh().modifyState(0, 1, true);
      lss_1D.setSupport(mesh_crack_front);

      xIter it = mesh_crack_front->begin(1);
      mEntity* edge_seed = *it;
      mVertex* v0 = (mVertex*)edge_seed->get(0, 0);
      mVertex* v1 = (mVertex*)edge_seed->get(0, 1);
      xtensor::xPoint p_mid = (v0->point() + v1->point()) * 0.5;
      mEntity* e_mesh2d = xMesh::get_const_was_created_by().at(*edge_seed);
      mEntity* e_mesh3d = xMesh::get_const_was_created_by().at(*e_mesh2d);
      mEntity* e;
      if (e_mesh3d->getLevel() != mesh->dim())
         e = e_mesh3d->get(mesh->dim(), 0);
      else
         e = e_mesh3d;
      xfem::xGeomElem elem(e);
      elem.setUVWForXYZ(p_mid);
      setLocalInfos(e, elem.getUVW());
      lss_1D(v0) = 0.0;
      xtensor::xVector<> v01(v0->point(), v1->point());
      lenght = v01.mag();

      lss_1D(v1) = 0.0 + lenght * ((v01 * n3) > 0 ? 1. : -1.);
      visited.setData(*v0) = 1;
      visited.setData(*v1) = 1;
      mVertex* seed_pos;
      mVertex* seed_neg;
      if (v01 * n3 > 0.)
      {
         seed_pos = v1;
         seed_neg = v0;
      }
      else
      {
         seed_pos = v0;
         seed_neg = v1;
      }

      // positive
      mEntity* e_current = edge_seed;
      mVertex* v_current = seed_pos;
      mEntity* e_next;
      mVertex* v_next;
      while ((e_next = otherEdgeOnVertex(v_current, e_current)))
      {
         //	if (debug) cout << "new pos found " << endl;
         v_next = otherVertexOnEdge(e_next, v_current);
         lenght += lengthOfEdge(e_next);
         double abss = lss_1D(v_current) + lengthOfEdge(e_next);
         // std::cout << abss << std::endl;
         if (!visited.getData(*v_next))
         {
            visited.setData(*v_next) = 1;
            lss_1D(v_next) = abss;
            v_current = v_next;
            e_current = e_next;
         }
         else
         {
            isloop = true;
            break;
         }
      }
      if (!isloop)
      {
         // negative
         e_current = edge_seed;
         v_current = seed_neg;
         while ((e_next = otherEdgeOnVertex(v_current, e_current)))
         {
            // if (debug) cout << "new neg found " << endl;
            v_next = otherVertexOnEdge(e_next, v_current);
            lenght += lengthOfEdge(e_next);
            double abss = lss_1D(v_current) - lengthOfEdge(e_next);
            // std::cout << abss << std::endl;
            // if (!v_next->getAttachedInt(visited_tag)){
            visited.setData(*v_next) = 1;
            lss_1D(v_next) = abss;
            v_current = v_next;
            e_current = e_next;
            //}
         }
      }
      // isloop = false;
      is_closed = isloop;
      //    cleaning up the tags
      for (mEntity* pe : mesh_crack_front->range(0)) visited.deleteData(*pe);
      // scaling the parametrization between mins and maxs for open crack
      xGetMinMaxLevelSetInspector cl;
      if (!isloop)
      {
         lss_1D.accept(cl);
         if (debug) cout << " ls front  min is " << cl.getMin() << " and max is " << cl.getMax() << endl;
         double s = (maxs - mins) / (cl.getMax() - cl.getMin());
         if (debug) std::cout << " s and lenght" << s << " " << lenght << std::endl;
         xScaleLevelSetModifier scale(s);
         lss_1D.accept(scale);
         xShiftLevelSetModifier shift(mins - s * cl.getMin());
         lss_1D.accept(shift);
      }
      else
      {
         xScaleLevelSetModifier scale(1. / lenght * 2 * M_PI);
         lss_1D.accept(scale);
      }
      if (debug) std::cout << "loop " << isloop << "lenght  " << lenght << std::endl;

      if (isloop)
      {
         lss_1D_cos.setSupport(mesh_crack_front);
         lss_1D_sin.setSupport(mesh_crack_front);
         for (mEntity* pe : mesh_crack_front->range(0))
         {
            lss_1D_cos(pe) = cos(lss_1D(pe));
            lss_1D_sin(pe) = sin(lss_1D(pe));
         }
      }

      lss.setSupport(narrow_band);
      //  on initialise le champ v_iplane3D a partir du champ v_iplane1D
      // xRegion reg(mesh,"supports_touched_by_front");
      xRegion reg(mesh, narrow_band_label);
      xVelocityInitialization lss_i(lss_1D, reg);
      lss.accept(lss_i, narrow_band);
      if (isloop)
      {
         lss_cos.setSupport(narrow_band);
         lss_sin.setSupport(narrow_band);
         xVelocityInitialization lss_i_cos(lss_1D_cos, reg);
         xVelocityInitialization lss_i_sin(lss_1D_sin, reg);
         lss_cos.accept(lss_i_cos, narrow_band);
         lss_sin.accept(lss_i_sin, narrow_band);
      }

      // trop lent pour l'instant
#if 0
        //  ... initialisations du pilot pour l extension
        //  double tol    = 1.e-6;
        //  int itmax3     = 20;
        xL2norm        l2n3;
        double dtvirt = xCFL::getDt(narrow_band);
        cout << "dtvirt in extend velocity " << dtvirt << endl;
        int max_iter = 150;
        //int max_iter = Param.max_iter_extend;
        xPilotError  statio3(l2n3, Param.tolerance, max_iter, dtvirt*Param.coeff_dt_virt);
        //extension a tout le domaine
        std::list<xSpaceOperator*> evo_lss;
        evo_lss.push_back(new xExtensionOperator(lsn));
        evo_lss.push_back(new xExtensionOperator(lst));
        xEvolveToStationary evol(evo_lss, statio3);
        cFirstMinusSecondCreator cminus(narrow_band_label,  touched_by_front_label);
        mesh->createSubsetEntities(all_minus_touched_by_front_label, cminus);
        lss.accept ( evol, xRegion(mesh, all_minus_touched_by_front_label));
#endif

      lss_created = true;
   }

   else
   {  // 2D case, no paramter needed, every thing set to one
      lss_1D.setSupport(mesh_crack_front, 1.);
      // lss_1D.exportGmsh("lss_1D");
      xexport::Export(lss_1D, pexport, "lss_1D");
      // extension au noeud proche
      lss.setSupport(narrow_band);
      xRegion reg(mesh, narrow_band_label);
      xVelocityInitialization lss_i(lss_1D, reg);
      lss.accept(lss_i, narrow_band);
      // lss.exportGmsh("lss");
      Export(lss, pexport, "lss");
      lss_created = true;
   }
}

mEntity* lCrack::otherEdgeOnVertex(mVertex* v, mEntity* e) const
{
   const bool debug = false;
   if (debug) cout << "v->size(1) " << v->size(1) << endl;
   assert((v->size(1) == 1) || (v->size(1) == 2));
   if (v->size(1) == 1) return nullptr;
   if (v->get(1, 0) == e) return v->get(1, 1);
   return v->get(1, 0);
}

mVertex* lCrack::otherVertexOnEdge(mEntity* e, mVertex* v) const
{
   if ((mVertex*)e->get(0, 0) == v) return (mVertex*)e->get(0, 1);
   return (mVertex*)e->get(0, 0);
}

double lCrack::lengthOfEdge(mEntity* e) const
{
   xtensor::xVector<> xvec(((mVertex*)e->get(0, 0))->point(), ((mVertex*)e->get(0, 1))->point());
   return xvec.mag();
}

void lCrack::getAsymptoticFields(const xfem::xGeomElem* geo_appro, int mode, xtensor::xTensor2<>& stress_aux,
                                 xtensor::xTensor2<>& grad_disp_aux) const
{
   const bool debug = false;

   // const double PI = 3.141592653589793238462;
   xMaterial* mat = xMaterialManagerSingleton::instance().getMaterial(geo_appro);
   const xTensors* properties = mat->getProperties();
   std::vector<double> MatInfo(2);
   MatInfo[0] = properties->scalar("YOUNG_MODULUS");
   MatInfo[1] = properties->scalar("POISSON_RATIO");

   setLocalInfos(geo_appro->getEntity(), geo_appro->getUVW());

   xtensor::xTensor2<> Stressl, Epsl, Gradl0Stress, Gradl1Stress, GradDisp;
   xtensor::xVector<> Displ, Gradl0Disp, Gradl1Disp, Gradl00Disp;
   xtensor::xVector<> Gradl01Disp, Gradl11Disp;
   GetMech3DAuxFields(local, MatInfo, Stressl, Epsl, Displ, Gradl0Disp, Gradl1Disp, Gradl0Stress, Gradl1Stress, Gradl00Disp,
                      Gradl11Disp, Gradl01Disp, mode);

   xtensor::xTensor2<>& ri = n123_row;
   xtensor::xTensor2<>& rio = no123_row;
   xtensor::xTensor2<>& leo = no123_col;
   xtensor::xTensor2<>& curvo0 = curvo1;
   xtensor::xTensor2<>& curvo1 = curvo2;
   xtensor::xTensor2<>& curvo2 = curvo3;

   stress_aux = (leo * (Stressl * rio));

   for (int a = 0; a < 3; a++)
   {
      GradDisp(a, 0) = Gradl0Disp(a) * ri(0, 0) + Gradl1Disp(a) * ri(1, 0);
      GradDisp(a, 1) = Gradl0Disp(a) * ri(0, 1) + Gradl1Disp(a) * ri(1, 1);
      GradDisp(a, 2) = Gradl0Disp(a) * ri(0, 2) + Gradl1Disp(a) * ri(1, 2);
   }

   // we now transform the classical derivatives into
   // covariant derivatives
   // page 36 and 37 BookI
   grad_disp_aux = (leo * GradDisp);
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         grad_disp_aux(i, j) += Displ(0) * curvo0(i, j) + Displ(1) * curvo1(i, j) + Displ(2) * curvo2(i, j);
      }
   }

   if (debug) cout << "grad_disp_aux =   " << grad_disp_aux << endl;

   return;
};

xEntityFilter lCrack::supportInCylinderAroundFrontFilter(double r) const { return tipfilterlcrack(*this, r); };

void lCrack::setV1Ddebug(std::function<xtensor::xVector<>(xtensor::xPoint)> func)
{
   mesh_crack_front->getMesh().modifyState(0, 1, true);
   for (mEntity* pe : mesh_crack_front->range(0))
   {
      AOMD::mVertex* v = static_cast<AOMD::mVertex*>(pe);
      v1D(v) = func(v->point());
      std::cout << v1D(v) << std::endl;
   }

   for (mEntity* pe : mesh_crack_front->range(1))
   {
      xtensor::xVector<> vv1 = v1D(static_cast<AOMD::mVertex*>(pe->get(0, 0)));
      xtensor::xVector<> vv2 = v1D(static_cast<AOMD::mVertex*>(pe->get(0, 1)));
      v1D(pe) = (vv1 + vv2) / 2;
      std::cout << v1D(pe) << std::endl;
   }
   /*	{
   string s="speedfront";
   char s2[10];
   sprintf(s2,"-%d",timestep);
   s=s+s2;
   v1D.exportGmsh(s);
   }
 */
   fromV1DToVioplane();
}

lCrack::~lCrack()
{
   if (mesh_crack_surface) delete mesh_crack_surface;
   if (mesh_crack_front) delete mesh_crack_front;
}

}  // namespace xcrack
