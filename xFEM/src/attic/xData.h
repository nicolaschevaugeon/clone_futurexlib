/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _DATA_H
#define _DATA_H

#include <map>
#include <string>

#include "xBoundary.h"
#include "xEnv.h"
#include "xMaterialManager.h"
#include "xMesh.h"
#include "xPointToDouble.h"

namespace xfem
{
class xZone;
class xZoneContainer;
class xBoundaryContainer;
class xPhysicalEnv;
class xMesh;

static const int NB_CHAR_MAX = 256;  // 48 was a bit limiting because it is used also for filenames...
static const int NB_MAX_ZONE_MAT = 1000;
static const int NB_MAX_ZONE_ENRICH = 1000;
static const int NB_MAX_ZONE_POSTPRO = 1000;
static const int NB_MAX_ENTITY = 1000;
// New for multiple crack tip postpro
static const int NB_MAX_FRONT_POSTPRO = 400;
// Size of the J-Domain Box
static const int CELLX = 2;
static const int CELLY = 2;
static const int CELLZ = 2;

class xData
{
  public:
   // BASIC CASE
   // WHAT IS ALWAYS NEEDED
   xData(void);
   ~xData(void);
   char* Getanalysis(void);
   char* GetProcedure(void);
   char* GetProcedureParametersFile(void);

   char* GetCrackDefinitionFile(void);
   void ReadInfo(const char* filename);
   /// Load the Mesh and store it in public member mesh.
   /**
      if a mesh was previously loaded and not cleaned up, the pointer mesh is not null.
      A warning will occur and the previous mesh is deleted.
      If running in Parallel with more than one processor,
      the file name set by setMeshfile is considered as a base name, the file that will be
      loaded by the current processor is a file with the current process number attached to it.
      For example if mesh_file is "mesh.msh" and the processor is proc 3,
      the file that will be read is mesh_3.msh"
      At the end of reading, there is a parallel barrier and then a boundary_link setup call to set up the communications paterns.
    **/
   void ReadMesh(void);
   void ReadZones(void);
   char* GetMeshFile(void);
   void SetMeshFile(const char* _mesh_file_name);
   char* GetSolverName(void);
   char* GetSolverParamFile(void);

   //
   // BASIC Informations ALWAYS NEEDED
   //
   int Dimension;  // Dimension of the problem: Derived from the mesh
   xMesh* mesh;
   xPhysicalEnv* PhysicalEnv;
   double LoadFactor;
   xZoneContainer* zones;  // material in each zone
   double time, dt;
   xMaterialManager& MaterialManager;

   xBoundaryContainer CurvesGroups;  // for the J integral in the case of
   // truly modeled cracks

   //  std::map<string, boost::function<double (const double&, const double&, const double&)*> level_sets;

   int do_postpro;
   int do_postpro_manual;
   int do_postpro_automatic;

   xBoundary allGroups;  // All piece of boundary for a given subdomain
   bool ComputeExactError;

   // to compare two solutions
   char formulation_one[NB_CHAR_MAX];
   char formulation_two[NB_CHAR_MAX];
   char solution_one_file[NB_CHAR_MAX];
   char solution_two_file[NB_CHAR_MAX];
   int solution_one_type;
   int solution_two_type;

   string export_format;
   xPointToDouble getLevelSet(string string_name) { return levelset[string_name]; }

  private:
   char analysis[NB_CHAR_MAX];
   char procedure[NB_CHAR_MAX];
   char regime[NB_CHAR_MAX];
   char mesh_file_type[NB_CHAR_MAX];
   char mesh_file[NB_CHAR_MAX];
   char geom_file[NB_CHAR_MAX];
   char geom_type[NB_CHAR_MAX];
   char integration_type[NB_CHAR_MAX];
   // char export_format[NB_CHAR_MAX];
   char formulation_param_file[NB_CHAR_MAX];
   char procedure_param_file[NB_CHAR_MAX];
   char crack_definition_file[NB_CHAR_MAX];
   char solver_name[NB_CHAR_MAX];
   char solver_param_file[NB_CHAR_MAX];
   char partition_file[NB_CHAR_MAX];
   char result_file[NB_CHAR_MAX];

   char reference_solution_file[NB_CHAR_MAX];
   char store_solution_file[NB_CHAR_MAX];

   int nb_zone_enrich;
   int zone_enrich_num[NB_MAX_ZONE_ENRICH];

   int nb_zone_postpro;
   int zone_postpro_num[NB_MAX_ZONE_POSTPRO];

   static char keyword_list[][NB_CHAR_MAX];
   static int keyword_sorted;
   static char analysis_list[][NB_CHAR_MAX];
   static int analysis_sorted;
   static char procedure_list[][NB_CHAR_MAX];
   static int procedure_sorted;
   static char regime_list[][NB_CHAR_MAX];
   static int regime_sorted;
   static char mesh_file_type_list[][NB_CHAR_MAX];
   static int mesh_file_type_sorted;
   static char geom_file_type_list[][NB_CHAR_MAX];
   static int geom_file_type_sorted;
   static char exact_solution_list[][NB_CHAR_MAX];
   static int exact_solution_sorted;

   typedef std::map<int, std::pair<std::string, std::string>> int_zones_t;
   int_zones_t int_zones;
   typedef std::map<std::string, std::pair<std::string, std::string>> string_zones_t;
   string_zones_t string_zones;
   // celine begin

   //  typedef boost::function<double (double&, double&, double&)> type_ls;
   //  std::map<std::string, type_ls > map_string_ls;

   std::map<std::string, xPointToDouble> levelset;

   // celine end
};

}  // namespace xfem

#endif
