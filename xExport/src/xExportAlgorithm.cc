#include "xExportAlgorithm.h"
#include "xExportGmsh.h"
#include "mpi.h"

namespace xexport{
void  Export(const xfem::xLevelSet& ls,
             xexport::xExport& pexport,
             const std::string& fieldName, const bool simplex)
{
  xfem::xEvalLevelSet<xtool::xIdentity<double> > eval_ls(ls);
  const xfem::xRegion& support = ls.getSupport();
  Export(eval_ls, pexport, fieldName, xfem::xIntegrationRuleBasic(0) , support.begin(), support.end(), simplex);
}

void Export(const xfem::xMesh * mesh,  xexport::xExport& pexport, const std::string& name)
{
  std::stringstream expname;
  expname << name;

  int proc_id, nb_proc;
  auto world = mesh->getPartitionManager().getComm();
  MPI_Comm_rank(world, &proc_id);
  MPI_Comm_size(world, &nb_proc);

  if (nb_proc> 1){
      expname  << "_"<< nb_proc << "_" <<  proc_id+1;
  }
  expname  << ".msh";
	AOMD::AOMD_Util::Instance()->ex_port(expname.str().c_str() ,  &mesh->getMesh() );
}

void Export(const xfem::xMesh & mesh,  xexport::xExport& pexport, const std::string& name)
{
  Export(&mesh, pexport, name);
}
}// end namespace xexport

