/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include "xMPIEnv.h"
#include "xRegion.h"
#include "xData.h"
#include "xExportGmsh.h"
#include "xCommandOnGeomElem.h"

#include "FMSK_fastmarching.h"
#include "FMSK_fastmarching_updater.h"
#include "meshinterfacexRegion.h"
#include "FMSK_skeleton_enriched_levelset.h"

 

/// I redefined this here to avoid the dependencies implies by including the function definition file from xfem library.
template <class ITER>
void ApplyCommandOnIntegrationRule(xfem::xCommandOnGeomElem& command, 
				   const xfem::xIntegrationRule& integration_rule, 
				   ITER it, ITER end, 
				   xfem::xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
  const bool debug = xfem::xdebug_flag;
  for ( ; it != end; ++it)
    {
      AOMD::mEntity* e_integ = *it;
      AOMD::mEntity* e_appro = integ2appro(e_integ);
      if (debug) std::cout<<"ApplyCommandOnIntegrationRule case #\n";
      if (debug) e_integ->print();
      if (debug) std::cout<<"E_appro:\n";
      if (debug) e_appro->print();
      xfem::xGeomElem geom_appro(e_appro);   
      xfem::xGeomElem geom_integ(e_integ); 
      command.setIntegElem(&geom_integ); 
      command.openApproxElem(&geom_appro);
      integration_rule.accept(command, e_integ);
      command.closeApproxElem(&geom_appro);
    }
}

/// This is a copy of xexport::xExport from xAlgorithm.h .... Because I don't want to be dependant on all the dependencies in xAlgorithm at ths point ...
template <typename T, class ITER>
void  mmExport(const xfem::xEval<T>& eval,
	       xexport::xExport& pexport,
	       const std::string& fieldName, 
	       const xfem::xIntegrationRule& integration_rule, 
	       ITER it, 
	       ITER end, const bool simplex=true, 
	       xfem::xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>()) 
{      

  if (!pexport.processStarted()) 
    {
      std::string fs = fieldName;
      pexport.openFile  ( fs.c_str() );
      pexport.startView ( fieldName.c_str() );
      xexport::xPlotCommand<T> plot(eval, pexport, simplex);
      ApplyCommandOnIntegrationRule(plot, integration_rule, it, end, integ2appro);    
      pexport.endView ();
      pexport.closeFile();   
    }
  else
    {
      pexport.startView ( fieldName.c_str() );
      xexport::xPlotCommand<T> plot(eval, pexport, simplex);
      ApplyCommandOnIntegrationRule(plot, integration_rule, it, end, integ2appro);    
      pexport.endView();
    }	  
  return;
}





struct options{
  std::string inputfilename;
};
  
options get_options(int argc, char *argv[]){
  cout <<"Usage : " << argv[0] << " inputfile" << endl;
  options opt;
  if (argc > 1){
    opt.inputfilename = argv[1];
  }
  return opt;
}
  

int main(int argc, char *argv[]){
  xtool::xMPIEnv::init(argc,argv);
  const double shock_treshold = M_PI/5.;
  const double fit_tol = 1.e-2;
 
  xfem::xData data;
  string pname;
  auto opt  = get_options(argc, argv);
  data.ReadInfo(opt.inputfilename.c_str());
  data.ReadMesh();
  data.ReadZones();

  xfem::xMesh &m = *(data.mesh);
  xfem::xRegion all(&m);
  typedef xfastmarching::meshinterfacexRegion meshinterface;
  meshinterface mi(all);
    
  xfastmarching::skeleton::vertexvectorstorage gls(mi);
  xfastmarching::skeleton::lstype ls(mi);
   
  for (auto env : *(data.PhysicalEnv)){
    std::cout <<  env.Phys << " " << env.Type << " "<< env.Entity << " "<< env.getDimension() << " "   << env.getValue()  << std::endl;
    if (env.Phys == "CUSTOM")
      if (env.Type == FIX){
	xfem::xClassRegion bc(&m, env.Entity, env.getDimension());
	for (auto  e : xtool::make_range(bc.begin(1),  bc.end(1))){
	  AOMD::mVertex &v0 = static_cast<AOMD::mVertex & >(*e->get(0,0));
	  AOMD::mVertex &v1 = static_cast<AOMD::mVertex & >(*e->get(0,1));
	  ls.set(v0, 0.);
	  ls.set(v1, 0.);
	}
      }
  }
  std::list<const AOMD::mVertex * >  known_init = ls.getSupport();
  std::function<double (const  meshinterface::vertex &) >  Ffunc = [](const  meshinterface::vertex & v ){return 1.;};
  std::list< const AOMD::mVertex * > final;
  std::list< const AOMD::mVertex * > knownv;
  xfastmarching::skeleton::FMupdaterGeneric2<meshinterface, xfastmarching::vector3d<double >, double , int> updater(mi, ls, Ffunc, 0., shock_treshold, gls);
  std::cout << "Start First FM"<< std::endl;
  fmeik_intern2( ls.getSupport().begin(), ls.getSupport().end(), knownv.begin(), knownv.end(), std::back_inserter(final), updater);
  
  // fmeik( mi, ls, ls.getSupport().begin(), ls.getSupport().end(), knownv.begin(), knownv.end(), std::back_inserter(final),  Ffunc,  0., gls );

  
  xfastmarching::skeleton::exportGMSHMine exp(m, mi, "results");
  std::function<double ( const AOMD::mVertex &v ) > fmls = [&ls ](const AOMD::mVertex &v ){ double val =100.; ls.get(v, val); return val; };
  
  
  std::function< std::array<double, 3 > ( const AOMD::mVertex &v ) > fmgls = [&gls ](const AOMD::mVertex &v ){
    xfastmarching::vector3d<double > g{0.,0.,0.};  gls.get(v, g);
    std::array<double, 3 > ret{g[0], g[1], g[2]};;
    return ret;
  };
  
  exp.appendVertexScalarView("ls", fmls);
  exp.appendVertexVectorView("gls", fmgls);
  
    
  
  {
    std::cout << "Starting Enriched ls " << std::endl;
    //auto  known_vrange = xtool::make_range( ls.getSupport().begin(), ls.getSupport().end() );
    std::list <const AOMD::mVertex * > trial_vrange;
    std::function <double (const AOMD::mVertex &) > input_ls    = [](const AOMD::mVertex & v){return 0.;};
    std::function <double (const AOMD::mVertex &) > input_trans = [](const AOMD::mVertex & v){return 0.;};
    xfastmarching::skeleton::skeleton_enriched_ls <meshinterface,  decltype(input_ls), decltype(input_trans)> enriched_ls(mi, known_init, trial_vrange, input_ls, input_trans, shock_treshold, fit_tol );
    std::cout << " enriched Mesh size " << enriched_ls.get_enriched_mesh().size(0) << " " << enriched_ls.get_enriched_mesh().size(1) << std::endl;   
  
    std::cout << "Done Enriched LS" << std::endl;


    //Exporting the enriched level set
    xexport::xExportGmshAscii pexport;
    xfem::xIntegrationRulePartition integ(0, enriched_ls.enriched_mesh.getGetPartition() );
    mmExport(enriched_ls.getEvalEnrLs(), pexport, "Enriched_Level_Set", integ, all.begin(), all.end() );
    
    // Exporting the skeleton
    xfem::xIntegrationRuleBasic integ0;
    xfem::xEvalConstant<double > one(1.);
    class eval_id : public xfem::xEval<double >{
    public:
      eval_id( ){}
      void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, result_type& res) const{
	res =  geo_appro->getEntity()->getId();
      }
    };
    
    class eval_parent_id : public xfem::xEval<double >{
    public:
      eval_parent_id(const xfastmarching::skeleton::embeded_mesh_fit_to_vertex & _em ):em(_em){}
      void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, result_type& res) const{
	auto pe = geo_appro->getEntity();
	const AOMD::mEntity &parent_e =  em.get_parent(*pe);
	res = parent_e.getId();
      }
    private:
      const xfastmarching::skeleton::embeded_mesh_fit_to_vertex & em;
    };
    
    eval_parent_id evalpid(enriched_ls.get_enriched_embeded_mesh() );
    eval_id        evalid;
    mmExport( evalpid, pexport, "Skeleton_pid", integ0, enriched_ls.getSkeletonRegion().begin(), enriched_ls.getSkeletonRegion().end() );
    mmExport( evalid, pexport, "Mesh_id_v", integ0, m.begin(0), m.end(0) );
    mmExport( evalid, pexport, "Mesh_id_e", integ0, m.begin(1), m.end(1) ); 
    mmExport( evalid, pexport, "Mesh_id_f", integ0, m.begin(2), m.end(2) ); 
    
    

  }



  
  return xtool::xMPIEnv::finalize();
}
