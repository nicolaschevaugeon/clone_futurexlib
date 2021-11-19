/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xData.h"

#include <cassert>
#include <cstdio>
#include <sstream>
#include <string>
#include <unordered_map>

#include "mAOMD.h"
#include "xEnv.h"
#include "xFiniteElement.h"
#include "xMesh.h"


namespace xfem
{
char xData::keyword_list[][NB_CHAR_MAX] = {"PROCEDURE",
                                           "ASS",
                                           "COUPLE_X",
                                           "COUPLE_Y",
                                           "CRACK",
                                           "CRACK_DIS",
                                           "DISPLACEMENT",
                                           "DISPLACEMENT_X",
                                           "DISPLACEMENT_Y",
                                           "DISPLACEMENT_Z",
                                           "DISPLACEMENT_R",
                                           "DISPLACEMENT_TETA",
                                           "DISPLACEMENT_ZC",
                                           "KINEMATIC_X",
                                           "KINEMATIC_Y",
                                           "KINEMATIC_Z",
                                           "CONTACT"
                                           "IMPROVEMENT",
                                           "IMPROVEMENT_X",
                                           "IMPROVEMENT_Y",
                                           "IMPROVEMENT_Z",
                                           "MULTIPLIER",
                                           "MULTIPLIER_X",
                                           "MULTIPLIER_Y",
                                           "MULTIPLIER_Z",
                                           "ANY",
                                           "ANY_X",
                                           "ANY_Y",
                                           "ANY_Z",
                                           "VECTOR_POTENTIAL_X",
                                           "VECTOR_POTENTIAL_Y",
                                           "VECTOR_POTENTIAL_Z",
                                           "VELOCITY_X",
                                           "VELOCITY_Y",
                                           "VELOCITY_Z",
                                           "VELOCITY",
                                           "ACCELERATION_X",
                                           "ACCELERATION_Y",
                                           "ACCELERATION_Z",
                                           "DISPLACEMENT_T",
                                           "DO_EIGEN_ANALYSIS",
                                           "DO_ENRICH_CHANNELS",
                                           "DO_ENRICH_CHANNEL_PRESSURE",
                                           "DO_ENRICH_MANUAL",
                                           "DO_ENRICH_AUTOMATIC",
                                           "DO_POSTPRO_MANUAL",
                                           "DO_POSTPRO_AUTOMATIC",
                                           "DOF_GROUP_SCA_1",
                                           "ELASTIC_SUPPORT",
                                           "ENDIF",
                                           "FIX",
                                           "FIX_AND_MEASURE",
                                           "FORMULATION_ONE",
                                           "FORMULATION_TWO",
                                           "GEOM_FILE",
                                           "GEOM_TYPE",
                                           "INTEGRATION_TYPE",
                                           "BC_VOLUME",
                                           "BC_SURFACE",
                                           "BC_LINE",
                                           "BC_POINT",
                                           "IC_VOLUME",
                                           "IC_SURFACE",
                                           "IC_LINE",
                                           "IC_POINT",
                                           "RIGID_LINE",
                                           "RIGID_SURFACE",
                                           "CONTACT_LINE",
                                           "GROUP_ON_CRACK",
                                           "HOLE",
                                           "PARTICLE",
                                           "PARTICLE_VEL",
                                           "PARTICLE_PRESS",
                                           "PARTICLE_WALL_GAP_VEL",
                                           "PARTICLE_WALL_GAP_PRESS",
                                           "IF",
                                           "LINE",
                                           "MATERIAL_PAIR",
                                           "MAT_CLASS",
                                           "MAT_PARAM",
                                           "MESH_FILE_TYPE",
                                           "MESH_FILE",
                                           "NEAR_TIP",
                                           "NEAR_REG_TIP",
                                           "NEAR_JUNCTION",
                                           "FORMULATION_PARAM_FILE",
                                           "PROCEDURE_PARAM_FILE",
                                           "PRESSURE",
                                           "POTENTIAL",
                                           "POTENTIAL_WEAK",
                                           "NORMAL_VELOCITY",
                                           "PARTITION_FILE",
                                           "POINT",
                                           "RESULT_FILE",
                                           "ROTATION_X",
                                           "ROTATION_Y",
                                           "SOLUTION_ONE_FILE",
                                           "SOLUTION_TWO_FILE",
                                           "SOLVER_NAME",
                                           "SOLVER_PARAM_FILE",
                                           "SURFACIC_HEAT_FLUX",
                                           "TEMPERATURE",
                                           "TEMPERATURE_MAC",
                                           "TRACTION_X",
                                           "TRACTION_Y",
                                           "TRACTION_Z",
                                           "ELAST_TRACTION_X",
                                           "ELAST_TRACTION_Y",
                                           "ELAST_TRACTION_Z",
                                           "ZONE",
                                           "ZONE_ENRICH",
                                           "ZONE_ENV",
                                           "ZONE_POSTPRO",
                                           "ZONE_NAME",
                                           "REFERENCE_SOLUTION_FILE",
                                           "STORE_SOLUTION_FILE",
                                           "STRESS",
                                           "SING_STRESS",
                                           "STRESS_XX",
                                           "STRESS_XY",
                                           "STRESS_XZ",
                                           "STRESS_YY",
                                           "STRESS_YZ",
                                           "STRESS_ZZ",
                                           "COMPUTE_EXACT_ERROR",
                                           "EXACT_SOLUTION_NAME",
                                           "NB_EVO",
                                           "EVO",
                                           "SHEARBAND",
                                           "SHEARBAND_DIS",
                                           "SHEARBAND_DIS_RIGID",
                                           "TANGENT_STIFFNESS",
                                           "NORMAL_STIFFNESS",
                                           "FRICTION",
                                           "WALL",
                                           "DENSITY",
                                           "EXPORT_FORMAT",
                                           "ELECTRIC_DISPLACEMENT",
                                           "CUSTOM",
                                           "CRACK_DEFINITION_FILE",
                                           "BC_VOLUME_NAME",
                                           "BC_SURFACE_NAME",
                                           "BC_LINE_NAME",
                                           "BC_POINT_NAME"};

int xData::keyword_sorted = 0;

char xData::procedure_list[][NB_CHAR_MAX] = {"crack_growth", "cohesive_crack_growth"};
int xData::procedure_sorted = 0;
//
char xData::exact_solution_list[][NB_CHAR_MAX] = {"patch_test",          "contact2d",       "hole",         "inclusion",
                                                  "bimaterial1d",        "stokescyl",       "crack_mode_I", "bimaterial3d",
                                                  "spherical_inclusion", "spherical_cavity"};
int xData::exact_solution_sorted = 0;

static int cmp(const void *vp, const void *vq) { return strcmp((const char *)vp, (const char *)vq); }

// DEFAULT CONSTRUCTOR
xData::xData()
    : mesh(nullptr),
      internal_read(false),
      entities_id(),
      PhysicalEnv(nullptr),
      MaterialManager(xMaterialManagerSingleton::instance())
{
   Dimension = 0;
   ComputeExactError = false;

   time = 1.0;
   dt = 1.0;

   PhysicalEnv = new xPhysicalEnv;

   strcpy(procedure, "");
   strcpy(crack_definition_file, "");
   strcpy(mesh_file_type, "");
   strcpy(mesh_file, "");
   strcpy(part_file, "");
   export_format = "ascii";
   strcpy(formulation_param_file, "");
   strcpy(solver_param_file, "");
   strcpy(solver_name, "");
   strcpy(partition_file, "");

   strcpy(reference_solution_file, "");
   strcpy(store_solution_file, "");

   if (!keyword_sorted)
   {
      keyword_sorted = 1;
      qsort(keyword_list, sizeof(keyword_list) / NB_CHAR_MAX, NB_CHAR_MAX, cmp);
   }

   if (!exact_solution_sorted)
   {
      exact_solution_sorted = 1;
      qsort(exact_solution_list, sizeof(exact_solution_list) / NB_CHAR_MAX, NB_CHAR_MAX, cmp);
   }

   if (!procedure_sorted)
   {
      procedure_sorted = 1;
      qsort(procedure_list, sizeof(procedure_list) / NB_CHAR_MAX, NB_CHAR_MAX, cmp);
   }

   xMaterialManagerSingleton::instance().registerMaterial("elastic", xMaterialCreator<xElastic>());

   return;
}

// DESTRUCTOR
xData::~xData()
{
   // must be clean before mesh
   entities_id.clear();
   // std::cout << "xDATA DESTRUCTOR" << std::endl;
   // a bug needs to be fix in aomd because
   // the destructor of mesh makes aomd crash
   // I think the bug is not in AOMD ... In fact it seems the destructor break done when deleling attachedData related to
   // xcut::xPhysSurfByTagging.
   delMesh();
   if (PhysicalEnv) delete PhysicalEnv;
}

xMesh *xData::getMesh()
{
   if (mesh)
      return mesh;
   else
   {
      throw -56485;
   }
};

void xData::setMesh(xMesh *pm)
{
   delMesh();
   mesh = pm;
};

void xData::delMesh()
{
   if (internal_read && mesh) delete mesh;
   mesh = nullptr;
   internal_read = false;
};

char *xData::GetMeshFile() { return mesh_file; }
char *xData::GetPartFile()
{
   if (part_file == std::string(""))
   {
      if (mesh_file == std::string(""))
      {
         cerr << "nor part_file or mesh_file name set in xData. Fixe mesh_file to be able to call GetPartFile" << endl;
         throw -456;
      }
      if (strncmp(&mesh_file[strlen(mesh_file) - 4], ".msh", 4))
      {
         cerr << mesh_file << " do not have the .msh extention as expected" << endl;
         throw -987;
      }
      string partition_name;
      partition_name.assign(mesh_file, strlen(mesh_file) - 4);
      partition_name += ".part";
      strcpy(part_file, partition_name.c_str());
   }
   return part_file;
}

char *xData::GetSolverName() { return solver_name; }

char *xData::GetSolverParamFile() { return solver_param_file; }
char *xData::GetProcedureParametersFile() { return procedure_param_file; }

// char * xData::Getanalysis(void)
// {
//   return analysis;
// }

char *xData::GetProcedure() { return procedure; }

char *xData::GetCrackDefinitionFile() { return crack_definition_file; }

void xData::ReadMesh(xinterface::aomd::xAttachedDataManagerAOMD<int> &target_proc, MPI_Comm world)
{
   ReadMesh(world, true);
   int proc_id;
   MPI_Comm_rank(world, &proc_id);
   string partition_name;
   partition_name.assign(mesh_file, strlen(mesh_file) - 4);
   partition_name += ".part";
   ifstream f(partition_name);
   if (f.good())
   {
      f.close();
      strcpy(part_file, partition_name.c_str());
      int dim = mesh->dim();
      std::unordered_map<int, int> id_to_proc;
      id_to_proc.reserve(mesh->size(dim));
      int e_id, p_id;
      ifstream f(partition_name);
      f >> e_id >> p_id;
      while (!f.eof())
      {
         id_to_proc.insert(std::make_pair(e_id, p_id));
         f >> e_id >> p_id;
      }
      for (auto *e : xtool::make_range(mesh->begin(dim), mesh->end(dim)))
      {
         int *pe_id = entities_id.getData(*e);
         assert(pe_id);
         target_proc.setData(*e) = id_to_proc.at(*pe_id);
      }
   }
}
void xData::ReadMesh(MPI_Comm world, bool read_entid)
{
   if (mesh_file == std::string(""))
   {
      cerr << "mesh_file name not set in xData" << endl;
      return;
   }
   if (strncmp(&mesh_file[strlen(mesh_file) - 4], ".msh", 4))
   {
      cerr << mesh_file << " do not have the .msh extention as expected" << endl;
      return;
   }
   int proc_id;
   MPI_Comm_rank(world, &proc_id);

   // if there is a file with sequential name only proc 0 reads it
   ifstream f(mesh_file);
   if (f.good())
   {
      f.close();
      if (mesh != nullptr)
      {
         std::cout << "Proc " << proc_id << " WARNING : non-null xMesh pointer reinitialized!" << endl;
         delMesh();
      }
      cout << " Read mesh : " << mesh_file << endl;
      if (proc_id)
         mesh = new xMesh(world);
      else
      {
         if (read_entid)
            mesh = new xMesh(mesh_file, entities_id, world);
         else
            mesh = new xMesh(mesh_file, world);
      }
      internal_read = true;
   }
   // if there is no file with sequential name try composed name
   else
   {
      string meshPartitionName;
      meshPartitionName.assign(mesh_file, strlen(mesh_file) - 4);
      meshPartitionName += "_proc" + std::to_string(proc_id) + ".msh";

      cout << " Read mesh : " << meshPartitionName << endl;

      if (mesh != nullptr)
      {
         std::cout << "Proc " << proc_id << " WARNING : non-null xMesh pointer reinitialized!" << endl;
         delMesh();
      }
      ifstream rf(meshPartitionName.c_str());
      if (!rf.good())
      {
         cout << "no mesh file with name " << meshPartitionName << "found " << endl;
         mesh = new xMesh(world);
      }
      else
      {
         rf.close();
         strcpy(mesh_file, meshPartitionName.c_str());
         if (read_entid)
            mesh = new xMesh(meshPartitionName, entities_id, world);
         else
            mesh = new xMesh(meshPartitionName, world);
      }
      internal_read = true;
   }
}

void xData::ReadZones()
{
   for (int_zones_t::const_iterator it = int_zones.begin(); it != int_zones.end(); ++it)
   {
      MaterialManager.createZone(it->first, it->second.first, it->second.second);
   }
   for (string_zones_t::const_iterator it = string_zones.begin(); it != string_zones.end(); ++it)
   {
      MaterialManager.createZone(it->first, it->second.first, it->second.second);
   }
}

inline void check_active(int active, const char *s)
{
   if (active != 0)
   {
      printf("Error in DAT file: you forgot to close a brace for the %s info\n", s);
      assert(1 == 0);
   }
   return;
}
void xData::SetMeshFile(const char *_mesh_file_name) { strcpy(mesh_file, _mesh_file_name); }

void xData::ReadInfo(const char *filename)
{
   FILE *fp = fopen(filename, "r");
   if (fp == nullptr)
   {
      fprintf(stderr, "The data file %s cannot be opened\n", filename);
      assert(1 == 0);
      throw;
   }
   char key[256];
   char c;
   int i;
   //  int tmpScf;
   int geom, entity, val_ass;
   // ajout celine 16 nov 2007
   char entity_name[NB_CHAR_MAX];
   //  std::string entity_name_s(entity_name);
   // end ajout celine
   string bc_name;
   string phys;
   double val_fix;
   xEnv env;
   xPhysicalEnv *InfoEnv = nullptr;

   PhysicalEnv->clear();

   char mat_class[NB_CHAR_MAX], mat_param[NB_CHAR_MAX], zone_name[NB_CHAR_MAX];
   // for the load evolution
   int nb_evo = 0;
   double factor, t;
   std::map<double, double> evo;

   const int NOTHING_ACTIVE = 0;
   const int BC_VOLUME_ACTIVE = 2;
   const int BC_SURFACE_ACTIVE = 3;
   const int BC_LINE_ACTIVE = 4;
   const int BC_POINT_ACTIVE = 5;
   const int ZONE_ENV_ACTIVE = 6;
   const int ZONE_ACTIVE = 10;
   const int ZONE_NAME_ACTIVE = 11;
   const int BC_VOLUME_NAME_ACTIVE = 12;
   const int BC_SURFACE_NAME_ACTIVE = 13;
   const int BC_LINE_NAME_ACTIVE = 14;
   const int BC_POINT_NAME_ACTIVE = 15;
   const int IC_VOLUME_ACTIVE = 16;
   const int IC_SURFACE_ACTIVE = 17;
   const int IC_LINE_ACTIVE = 18;
   const int IC_POINT_ACTIVE = 19;
   const int RIGID_SURFACE_ACTIVE = 20;
   const int RIGID_LINE_ACTIVE = 21;
   const int CONTACT_LINE_ACTIVE = 22;

   int active = NOTHING_ACTIVE;

   while ((i = fscanf(fp, "%[1234567890A-Z_]", key)) != EOF)
   {
      if (i == 0)
      {
         c = fgetc(fp);
         // printf("c: %c\n", c);

         if (c != '#' && c != '\n' && c != ' ' && c != '\t' && c != '}')
         {
            fprintf(stderr, "The following character is not known : %c\n", c);
            assert(1 == 0);
         }
         if (c == '#') fscanf(fp, "%*[^\n] \n");

         if (c == '}')
         {
            if (active == BC_LINE_ACTIVE)
            {
            }
            else if (active == BC_SURFACE_ACTIVE)
            {
            }
            else if (active == BC_VOLUME_ACTIVE)
            {
            }
            else if (active == BC_POINT_ACTIVE)
            {
            }
            // ajout celine 21 nov 2007
            else if (active == BC_LINE_NAME_ACTIVE)
            {
            }
            else if (active == BC_SURFACE_NAME_ACTIVE)
            {
            }
            else if (active == BC_VOLUME_NAME_ACTIVE)
            {
            }
            else if (active == BC_POINT_NAME_ACTIVE)
            {
            }
            // end celine 21 nov 2007
            else if (active == IC_LINE_ACTIVE)
            {
            }
            else if (active == IC_SURFACE_ACTIVE)
            {
            }
            else if (active == IC_VOLUME_ACTIVE)
            {
            }
            else if (active == IC_POINT_ACTIVE)
            {
            }
            else if (active == RIGID_LINE_ACTIVE)
            {
            }
            else if (active == RIGID_SURFACE_ACTIVE)
            {
            }
            else if (active == CONTACT_LINE_ACTIVE)
            {
            }
            else if (active == ZONE_ENV_ACTIVE)
            {
            }
            else if (active == ZONE_ACTIVE)
            {
               std::string mat_class_s(mat_class);
               std::string mat_param_s(mat_param);
               int_zones.insert(std::make_pair(entity, std::make_pair(mat_class_s, mat_param_s)));
               // MaterialManager.createZone(entity, mat_class, mat_param);
            }
            else if (active == ZONE_NAME_ACTIVE)
            {
               // MaterialManager.createZone(zone_name, mat_class, mat_param);
               std::string mat_class_s(mat_class);
               std::string mat_param_s(mat_param);
               std::string zone_name_s(zone_name);
               string_zones.insert(std::make_pair(zone_name_s, std::make_pair(mat_class_s, mat_param_s)));
            }
            else if (active == NOTHING_ACTIVE)
            {
               printf("Error in DAT file :you closed a braced without opening one\n");
               assert(0);
            }
            else
            {
               assert(0);
            }
            // we clear things before reading the new FOO = 3 { .... }
            active = NOTHING_ACTIVE;
            entity = -1;
         }
      }
      else
      {
         //
         // Check if the key is valid
         //
         if (bsearch((const void *)key, keyword_list, sizeof(xData::keyword_list) / NB_CHAR_MAX, NB_CHAR_MAX, cmp) == nullptr)
         {
            fprintf(stderr, "The keyword %s is not known\n", key);
            assert(1 == 0);
         }

         if (strcmp(key, "ENDIF") == 0)
         {
            // nothing to do
         }

         //
         // action depending on the key
         //
         if (strcmp(key, "PROCEDURE") == 0)
         {
            fscanf(fp, "%*[\n= ] %s[a-z]", procedure);
            if (bsearch((const void *)procedure, procedure_list, sizeof(xData::procedure_list) / NB_CHAR_MAX, NB_CHAR_MAX, cmp) ==
                nullptr)
            {
               printf("The type of procedure %s is not known\n", procedure);
               assert(1 == 0);
            }
         }
         if (strcmp(key, "MESH_FILE_TYPE") == 0)
         {
            fscanf(fp, "%*[\n= ] %s[a-z]", mesh_file_type);
         }
         if (strcmp(key, "MESH_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", mesh_file);
         if (strcmp(key, "EXPORT_FORMAT") == 0)
         {
            char tmp[NB_CHAR_MAX];
            fscanf(fp, "%*[\n= ] %s[a-z]", tmp);
            export_format = tmp;
            if ((export_format != "binary") && (export_format != "ascii"))
            {
               fprintf(stderr, "EXPORT_FOMAT is %s\n", export_format.c_str());
               fprintf(stderr, "EXPORT_FORMAT should be ascii or binary\n");
               assert(0);
            }
         }
         if (strcmp(key, "GEOM_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", geom_file);
         if (strcmp(key, "GEOM_TYPE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", geom_type);
         if (strcmp(key, "INTEGRATION_TYPE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", integration_type);
         if (strcmp(key, "FORMULATION_PARAM_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", formulation_param_file);
         if (strcmp(key, "PROCEDURE_PARAM_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", procedure_param_file);
         if (strcmp(key, "CRACK_DEFINITION_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", crack_definition_file);

         if (strcmp(key, "SOLVER_PARAM_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", solver_param_file);
         if (strcmp(key, "SOLVER_NAME") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", solver_name);

         if (strcmp(key, "PARTITION_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", partition_file);
         if (strcmp(key, "RESULT_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", result_file);

         if (strcmp(key, "COMPUTE_EXACT_ERROR") == 0)
         {
            ComputeExactError = true;
         }

         // STRESS needed
         if (strcmp(key, "COMPUTE_EXACT_ERROR") == 0) ComputeExactError = true;

         if (strcmp(key, "REFERENCE_SOLUTION_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", reference_solution_file);
         if (strcmp(key, "STORE_SOLUTION_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", store_solution_file);

         if (strcmp(key, "SOLUTION_ONE_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", solution_one_file);
         if (strcmp(key, "SOLUTION_TWO_FILE") == 0) fscanf(fp, "%*[\n= ] %s[a-z]", solution_two_file);

         if (strcmp(key, "ZONE") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            check_active(active, "ZONE");
            active = ZONE_ACTIVE;
         }

         if (strcmp(key, "ZONE_NAME") == 0)
         {
            fscanf(fp, "%s%*[\n= {]", zone_name);
            check_active(active, "ZONE_NAME");
            active = ZONE_NAME_ACTIVE;
         }

         if (strcmp(key, "MAT_CLASS") == 0) fscanf(fp, "%*[\n= ] %[a-zA-Z0123456789_]", mat_class);
         if (strcmp(key, "MAT_PARAM") == 0)
         {
            fscanf(fp, "%*[\n= ] %[a-zA-Z0123456789_./*]", mat_param);
         }

         if (strcmp(key, "BC_VOLUME") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "BC_VOLUME");
            geom = BC_VOLUME;
            active = BC_VOLUME_ACTIVE;
         }
         if (strcmp(key, "IC_VOLUME") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "IC_VOLUME");
            geom = IC_VOLUME;
            active = IC_VOLUME_ACTIVE;
         }
         if (strcmp(key, "BC_VOLUME_NAME") == 0)
         {
            fscanf(fp, "%s%*[\n= {]", entity_name);
            InfoEnv = PhysicalEnv;
            check_active(active, "BC_VOLUME_NAME");
            geom = BC_VOLUME_NAME;
            active = BC_VOLUME_NAME_ACTIVE;
            env.set(entity_name, geom);
         }
         if (strcmp(key, "BC_SURFACE") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "BC_SURFACE");
            geom = BC_SURFACE;
            active = BC_SURFACE_ACTIVE;
            env.set(entity, geom);
         }
         if (strcmp(key, "IC_SURFACE") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "IC_SURFACE");
            geom = IC_SURFACE;
            active = IC_SURFACE_ACTIVE;
            env.set(entity, geom);
         }
         if (strcmp(key, "RIGID_SURFACE") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "RIGID_SURFACE");
            geom = RIGID_SURFACE;
            active = RIGID_SURFACE_ACTIVE;
            env.set(entity, geom);
         }
         if (strcmp(key, "BC_SURFACE_NAME") == 0)
         {
            fscanf(fp, "%s%*[\n= {]", entity_name);
            InfoEnv = PhysicalEnv;
            check_active(active, "BC_SURFACE_NAME");
            geom = BC_SURFACE_NAME;
            active = BC_SURFACE_NAME_ACTIVE;
            cout << entity_name << endl;
            env.set(entity_name, geom);
         }

         if (strcmp(key, "BC_LINE") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "BC_LINE");
            geom = BC_LINE;
            active = BC_LINE_ACTIVE;
         }
         if (strcmp(key, "IC_LINE") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "IC_LINE");
            geom = IC_LINE;
            active = IC_LINE_ACTIVE;
         }
         if (strcmp(key, "RIGID_LINE") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "RIGID_LINE");
            geom = RIGID_LINE;
            active = RIGID_LINE_ACTIVE;
         }
         if (strcmp(key, "CONTACT_LINE") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "CONTACT_LINE");
            geom = CONTACT_LINE;
            active = CONTACT_LINE_ACTIVE;
         }
         if (strcmp(key, "BC_LINE_NAME") == 0)
         {
            fscanf(fp, "%s%*[\n= {]", entity_name);
            InfoEnv = PhysicalEnv;
            check_active(active, "BC_LINE_NAME");
            geom = BC_LINE_NAME;
            active = BC_LINE_NAME_ACTIVE;
            env.set(entity_name, geom);
         }

         if (strcmp(key, "BC_POINT") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "BC_POINT");
            geom = BC_POINT;
            active = BC_POINT_ACTIVE;
         }
         if (strcmp(key, "IC_POINT") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            check_active(active, "IC_POINT");
            geom = IC_POINT;
            active = IC_POINT_ACTIVE;
         }

         if (strcmp(key, "BC_POINT_NAME") == 0)
         {
            fscanf(fp, "%s%*[\n= {]", entity_name);
            InfoEnv = PhysicalEnv;
            check_active(active, "BC_POINT_NAME");
            geom = BC_POINT_NAME;
            active = BC_POINT_NAME_ACTIVE;
            env.set(entity_name, geom);
         }

         // JF NEW
         if (strcmp(key, "ZONE_ENV") == 0)
         {
            fscanf(fp, "%d%*[\n= {]", &entity);
            InfoEnv = PhysicalEnv;
            geom = ZONE_ENV;
            check_active(active, "ZONE_ENV");
            active = ZONE_ENV_ACTIVE;
         }

         // JF END

         if (strcmp(key, "DISPLACEMENT_X") == 0) phys = key;
         if (strcmp(key, "DISPLACEMENT_Y") == 0) phys = key;
         if (strcmp(key, "DISPLACEMENT_Z") == 0) phys = key;
         if (strcmp(key, "DISPLACEMENT_R") == 0) phys = key;
         if (strcmp(key, "DISPLACEMENT_TETA") == 0) phys = key;
         if (strcmp(key, "DISPLACEMENT_ZC") == 0) phys = key;
         if (strcmp(key, "KINEMATIC_X") == 0) phys = key;
         if (strcmp(key, "KINEMATIC_Y") == 0) phys = key;
         if (strcmp(key, "KINEMATIC_Z") == 0) phys = key;
         if (strcmp(key, "CONTACT") == 0) phys = key;
         if (strcmp(key, "MULTIPLIER_X") == 0) phys = key;
         if (strcmp(key, "MULTIPLIER_Y") == 0) phys = key;
         if (strcmp(key, "MULTIPLIER_Z") == 0) phys = key;
         if (strcmp(key, "IMPROVEMENT_X") == 0) phys = key;
         if (strcmp(key, "IMPROVEMENT_Y") == 0) phys = key;
         if (strcmp(key, "IMPROVEMENT_Z") == 0) phys = key;
         if (strcmp(key, "ANY_X") == 0) phys = key;
         if (strcmp(key, "ANY_Y") == 0) phys = key;
         if (strcmp(key, "ANY_Z") == 0) phys = key;
         if (strcmp(key, "VECTOR_POTENTIAL_X") == 0) phys = key;
         if (strcmp(key, "VECTOR_POTENTIAL_Y") == 0) phys = key;
         if (strcmp(key, "VECTOR_POTENTIAL_Z") == 0) phys = key;

         if (strcmp(key, "ROTATION_X") == 0) phys = key;
         if (strcmp(key, "ROTATION_Y") == 0) phys = key;
         if (strcmp(key, "COUPLE_X") == 0) phys = key;
         if (strcmp(key, "COUPLE_Y") == 0) phys = key;
         if (strcmp(key, "TRACTION_X") == 0) phys = key;
         if (strcmp(key, "TRACTION_Y") == 0) phys = key;
         if (strcmp(key, "TRACTION_Z") == 0) phys = key;
         if (strcmp(key, "ELAST_TRACTION_X") == 0) phys = key;
         if (strcmp(key, "ELAST_TRACTION_Y") == 0) phys = key;
         if (strcmp(key, "ELAST_TRACTION_Z") == 0) phys = key;

         if (strcmp(key, "STRESS") == 0) phys = key;
         if (strcmp(key, "STRESS_XX") == 0) phys = key;
         if (strcmp(key, "STRESS_XY") == 0) phys = key;
         if (strcmp(key, "STRESS_XZ") == 0) phys = key;
         if (strcmp(key, "STRESS_YY") == 0) phys = key;
         if (strcmp(key, "STRESS_YZ") == 0) phys = key;
         if (strcmp(key, "STRESS_ZZ") == 0) phys = key;
         if (strcmp(key, "SING_STRESS") == 0) phys = key;
         if (strcmp(key, "DISPLACEMENT") == 0) phys = key;
         if (strcmp(key, "MULTIPLIER") == 0) phys = key;
         if (strcmp(key, "IMPROVEMENT") == 0) phys = key;
         if (strcmp(key, "ANY") == 0) phys = key;
         if (strcmp(key, "PRESSURE") == 0) phys = key;
         if (strcmp(key, "POTENTIAL") == 0) phys = key;
         if (strcmp(key, "POTENTIAL_WEAK") == 0) phys = key;
         if (strcmp(key, "ELECTRIC_DISPLACEMENT") == 0) phys = key;
         if (strcmp(key, "TEMPERATURE") == 0) phys = key;
         if (strcmp(key, "TEMPERATURE_MAC") == 0) phys = key;
         if (strcmp(key, "SURFACIC_HEAT_FLUX") == 0) phys = key;
         if (strcmp(key, "NORMAL_VELOCITY") == 0) phys = key;
         if (strcmp(key, "ELASTIC_SUPPORT") == 0) phys = key;

         if (strcmp(key, "VELOCITY_X") == 0) phys = key;
         if (strcmp(key, "VELOCITY_Y") == 0) phys = key;
         if (strcmp(key, "VELOCITY_Z") == 0) phys = key;
         if (strcmp(key, "VELOCITY") == 0) phys = key;

         if (strcmp(key, "ACCELERATION_X") == 0) phys = key;
         if (strcmp(key, "ACCELERATION_Y") == 0) phys = key;
         if (strcmp(key, "ACCELERATION_Z") == 0) phys = key;
         if (strcmp(key, "CUSTOM") == 0) phys = key;

         if (strcmp(key, "FIX") == 0)
         {
            //	    std:: cout << phys<< " "  << geom <<" "  << entity << " "<< val_fix << std::endl;
            fscanf(fp, "%*[\n= ] %lf", &val_fix);
            env.defineFixed(phys, geom, entity, val_fix);
            InfoEnv->add(env);
            env.clear();
         }
         if (strcmp(key, "FIX_AND_MEASURE") == 0)
         {
            fscanf(fp, "%*[\n= ] %lf", &val_fix);
            env.defineFixedAndMeas(phys, geom, entity, val_fix);
            InfoEnv->add(env);
            env.clear();
         }
         if (strcmp(key, "ASS") == 0)
         {
            fscanf(fp, "%*[\n= ] %d", &val_ass);
            env.defineAssociate(phys, geom, entity, val_ass);
            InfoEnv->add(env);
            env.clear();
         }
         if (strcmp(key, "NB_EVO") == 0)
         {
            fscanf(fp, "%*[\n= ] %d", &nb_evo);
         }
         if (strcmp(key, "EVO") == 0)
         {
            fscanf(fp, "%*[\n= ]");
            evo.clear();
            for (i = 0; i < nb_evo; i++)
            {
               fscanf(fp, "%*[( ] %lf %*[ ,] %lf %*[ )]", &t, &factor);
               if (evo.find(t) == evo.end())
               {
                  evo[t] = factor;
               }
               else
               {
                  fprintf(stderr,
                          "You cannot have two factor for the same time\n"
                          "i.e. the loading evolution must be continuous\n"
                          "the time %12.5e is there twice\n",
                          t);
               }
            }
            env.setEvolution(evo);
         }
      }
   }

   fclose(fp);

   return;
}

}  // namespace xfem
