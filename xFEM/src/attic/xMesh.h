/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _AOMDEXT_H_
#define _AOMDEXT_H_

#include <string>
#include <map>
#include <memory>
#include <utility>
#include "mAOMD.h"
#include "GEntity.h"
#include "mMesh.h"
#include "mVertex.h"
#include "mIterator.h"
#include "mExchangeData.h"
#include "AOMD_SharedInfo.h"
#include "xEntityFilter.h"
#include "xEntityToEntity.h"

class Octree;

namespace xfem
{
  

typedef AOMD::mMesh::iter    xIter;
typedef AOMD::mMesh::iterall xIterall;
typedef std::set<AOMD::mEntity*>   xPartition;
typedef AOMD::mClassIterator xClassIter;

class xSubMeshCreator;
class xRegularGrid;
class xRegion;
class mPyramid;
class xLevelSet;
class xSubMesh;

using std::map;
using std::vector;
using std::set;
using std::string;

/// This function return the "source " of e :
/*!  
    It check for the tag was_duplicated_from. 
    If it exist it call it self on the    entity pointed by the tag.
    else it return the input entity.
    In short, the function return e, if e is not the copy of an other entity,
    other wise, it goes up the hierachies of copies to return the "source" entity.
*/
 AOMD::mEntity *getSource(AOMD::mEntity * e);

 std::pair<Trellis_Util::mPoint, Trellis_Util::mPoint> compute_bounding_box(AOMD::mEntity * e) ;

/// the xMesh adds a functionality to the mMesh class of the Aomd package.
/*! 
  It adds the functionnality to create subset of a mesh.
  a subset is a group of AOMD::mEntity *, for which a name is associated.
  xMesh offer the possibility to iterate on a subset.
!*/
class xMesh  : public AOMD::mMesh {
public: 
  xMesh(int id = 1);
  xMesh(const string& filename);
  void modifyAllState();
  void modifyAllStateFalse();
  ~xMesh();

  /// return 2 points forming the bounding box of the mesh
  void compute_bounding_box(Trellis_Util::mPoint& min,Trellis_Util::mPoint &max) const;
  

    
  /// return the dimension (0, 1, 2, 3 ) of the highest level entity in the mesh (local to the proc)
  int dim() const;

  /// return the dimension (0, 1, 2, 3 ) of the highest level entity in the mesh (global  : max for the mesh on all proc)
  int dim_global() const;
  
  /// Mirror vertex are used for periodic boundary condition. 
  bool isMirror(AOMD::mEntity*); 
  AOMD::mMirrorVertex* getMirrorVertex(AOMD::mEntity*);
  
  //change the del function of the base class
  void del (AOMD::mEntity *);
  
  // new functions to deal with subsets   
  xSubMesh & createSubMesh(const string& subMeshName) const;
  xSubMesh & createSubMesh(const string& subMeshName, xSubMeshCreator& creator) const;
  void deleteSubMesh(const string& subset) const;
  void deleteAllSubMesh() const;
  xSubMesh & getSubMesh(const string& sub) const;
  
 
  //return the support of the entity e
  //the level of all the elements in l
  //is the same as dim(), all the elements in l are leaves
  //lookupSupport add things to l
  //So make sure it is empty before calling it
  //lookupSupport takes into account the fact that the entity can be a mirror
  void lookupSupport(AOMD::mEntity* e, std::set<AOMD::mEntity*> &l) ;


  /// check mesh : calculate bound and average "volume" value of elements.
  //! if there is a problem mim value are fare less then average.
  //! fact is a thresold to print entity with a "volume" less then 1/fact times the average "volume"
  std::pair<double,double> checkMesh(const int fact=100000000) const;

  //Finding the hanging nodes
  //map vertex to edge or face
  typedef std::map<AOMD::mEntity*, AOMD::mEntity*> hanging_t;
  void findHangingNodes();
  void periodicAssociation();

  static void getPartition(AOMD::mEntity* e, xPartition& partition, xEntityFilter filter=xAcceptAll());
  // Why is this Static ???

  // Ajout K
  void getPartitionMesh(xMesh & partitionMesh);
  void getOctreeMesh(xMesh & octreeMesh);  /// 20-07-09
  /// get the partition attached to entity e
  static void getPartitionEntity(AOMD::mEntity* e, std::set<AOMD::mEntity*>& partition);

  //  static void createSubSimplices(xRegion& support, xMesh* x_interface,
  //				  std::hash_map<AOMD::mEntity*,double,AOMD::EntityHashKey,AOMD::EntityEqualKey>& ls);
  // Where is the implementation ??? not in Xfem thow
  static void createSubSimplices(xRegion& support, xMesh* x_interface,
				             std::unordered_map<AOMD::mEntity*,double,AOMD::EntityHashKey,AOMD::EntityEqualKey>& ls);

  
  // some of those tag are useless !!!!! To be cleaned up !
  static unsigned int was_created_by_tag;
  static unsigned int is_the_creator_of_tag;
  static unsigned int r_on_edge_tag;
  static unsigned int is_duplicated_in_tag;
  static unsigned int was_duplicated_from_tag;
  static unsigned int partition_tag;
  static unsigned int is_in_partition_of_tag;

  static unsigned int was_created_by_tag2;
  static unsigned int is_the_creator_of_tag2;
  static unsigned int is_the_copy_of_tag;
  static unsigned int has_a_copy_tag;
  static unsigned int is_the_element0_of_tag;
  static unsigned int refined_elements_tag;
  static unsigned int levelset_value_tag ;
  static unsigned int cut_edge_tag;
  static unsigned int cut_element_tag;

  //for octree meshes
  static unsigned int is_hanging_on_tag;
  static unsigned int is_hanging_by_tag;
  static unsigned int octree_level_tag;
  static unsigned int down_groupe_tag;
  static unsigned int bnd_groupe_tag;
  
  //for ramp heavisid enrichment
  static unsigned int RH_enrichment_nb_func_tag;
  static unsigned int RH_enrichment_side_tag;
  static unsigned int RH_enrichment_status_tag;

  static void setCurrentPartition(const char* s) 
    {
      partition_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId(s);
    }
 
private:

  xRegularGrid *RegularGrid;

  // access functions to the subsets
  // Subset are implemented as a map between the string that named the subset and an
  // Entity Container. Note :  Could be more efficient too implement it as a tag 
  // associated to the string attached to each Entity that is in the sub
  mutable map<string, xSubMesh *> subsetEntities;
  typedef map<string, xSubMesh *>::const_iterator const_iter_subsets;
  typedef map<string, xSubMesh *>::iterator       iter_subsets;
 
  void lookupSupportBasic(AOMD::mEntity* e, std::set<AOMD::mEntity*> &l);

  /// cutmesh utilities
  void createSubSimplices(const xLevelSet &lst, xMesh* m_interface, xEntityToEntity& class_in,
			                   xEntityToEntity& class_out);
  //  AOMD::mVertex* vertexInBetween(AOMD::mEntity* v0, AOMD::mEntity* v1, AOMD::mEntity* q);
  AOMD::mVertex* vertexInBetween(AOMD::mEntity* v0, AOMD::mEntity* v1, const std::vector<AOMD::mVertex*>& q);
  AOMD::mVertex* vertexInBetween(AOMD::mEntity* v0, AOMD::mEntity* v1);
  void addNewTo(AOMD::mEntity* e, const std::string& side,
		xEntityToEntity& class_in, 
		xEntityToEntity& class_out);
  // For test only :  to create vertex unically .
  //std::map<Trellis_Util::mPoint, AOMD::mVertex *> pointVertexMap;

  public :
  //  AOMD::mVertex *createVertexUnique (const double x, const double y, const double z , pGEntity classif);
  
  void cutMesh(const xLevelSet& lst, xMesh* x_interface,xEntityFilter filter=xAcceptAll(),
	       xEntityToEntity classify_in=xtool::xIdentity<AOMD::mEntity*>(),
	       xEntityToEntity classify_out=xtool::xIdentity<AOMD::mEntity*>(),		     
	       bool create_partition=true,bool keep_old_partition=false,bool recursive=false);
  
  void cleanPartition(void);
  void cleanClassification(void);


  void takeTraceOn(const xLevelSet &ls, xLevelSet& l_interfaced);

  private :
  mutable Octree* octree;

  public:

  //replace by octree when aomd octree available
  void createOctree(void) const;
  void locateElementOctree(const Trellis_Util::mPoint& p, std::set<AOMD::mEntity*>&) const;
  void locateElement(const Trellis_Util::mPoint& p, std::set<AOMD::mEntity*>&);

  static unsigned int get_was_created_by_tag() {return  was_created_by_tag;}
  static unsigned int get_is_the_creator_of_tag(){return is_the_creator_of_tag;}
  static unsigned int get_r_on_edge_tag(){return r_on_edge_tag;}
  static unsigned int get_is_duplicated_in_tag(){return is_duplicated_in_tag;}
  static unsigned int get_was_duplicated_from_tag(){return was_duplicated_from_tag;}
  static unsigned int get_partition_tag(){return partition_tag;}
  static unsigned int get_is_in_partition_of_tag(){return is_in_partition_of_tag;}
  
  static unsigned int get_was_created_by_tag2() {return  was_created_by_tag2;}
  static unsigned int get_is_the_creator_of_tag2(){return is_the_creator_of_tag2;}
  static unsigned int get_is_the_copy_of_tag(){return is_the_copy_of_tag;}
  static unsigned int get_has_a_copy_tag(){return has_a_copy_tag;}
  static unsigned int get_is_the_element0_of_tag(){return is_the_element0_of_tag;}
  static unsigned int get_refined_elements_tag(){return refined_elements_tag ;} 
  static unsigned int get_cut_edge_tag(){return cut_edge_tag ;} 
  static unsigned int get_cut_element_tag(){return cut_element_tag ;} 
    

  //for octree meshes created in MeshAdapt/AMR
  static unsigned int get_is_hanging_on_tag(){return is_hanging_on_tag;}
  static unsigned int get_is_hanging_by_tag(){return is_hanging_by_tag;}
  static unsigned int get_octree_level_tag(){return octree_level_tag;}
  static unsigned int get_down_groupe_tag(){return down_groupe_tag;}
  static unsigned int get_bnd_groupe_tag(){return bnd_groupe_tag;}
  
  int getOctreeLevelMax() const {return level_max;}
  void setOctreeLevelMax(int l) {level_max = l;}
  int level_max;

  void createInterfaceFromLevelSet(const xLevelSet& ls);

  void createInterfaceOnSubmesh(const xLevelSet& ls, xMesh* interface);

  xMesh* cutElement (AOMD::mEntity* e, AOMD::mEntity* e_interface);
  static int nbCutElementCall;
  void createPartitionFromOneLevelSet(const xLevelSet& ls, xMesh* interface, 
				      xEntityToEntity classify_in=xtool::xIdentity<AOMD::mEntity*>(),
				      xEntityToEntity classify_out=xtool::xIdentity<AOMD::mEntity*>(),
				      xEntityToEntity classify_interface=xtool::xIdentity<AOMD::mEntity*>(), bool recursive = false);
    
  /// Used to classify sub-element on one side or the other of a levelset.
  void classifyElement(AOMD::mEntity *e1, AOMD::mEntity* e2, 
		       const xLevelSet& ls, 
		       xEntityToEntity& classify_in,
		       xEntityToEntity& classify_out); 

   
  Trellis_Util::mPoint centroid(AOMD::mEntity* e);
  AOMD::mEntity* getOriginalCreator(AOMD::mEntity *egros);
  void cutAlongInterface(xMesh* interface);
  void cutAlongInterfaceRecursive(xMesh* interface, const xLevelSet& ls);

  void processVertex(AOMD::mVertex *v, const xLevelSet &ls);
  void processEdge(AOMD::mEntity *e, const xLevelSet &ls);
  void attachSimplicies(AOMD::mEntity *e);

public:
  
  ///                          the elements of the mesh are transformed into TRIiangles and TETrahedra elements 
  ///                          with inheritance of tags. Division is made by cutting faces on the diagonal which 
  ///                          contains the node with the lowest Id. In this way the mesh remains conform.
  void Simplexify(xEntityFilter filter=xAcceptAll());
  ///                          transforms a QUAD given by its pointer to TRI added to the xMesh (used by tretriseMesh)
  void QuadToTri(AOMD::mEntity*  );
  ///                          transforms a HEX given by its pointer to TET added to the xMesh (used by tretriseMesh)
  void HexToTet(AOMD::mEntity* );
  ///                          transforms a PRISM given by its pointer to TET added to the xMesh (used by tretriseMesh)
  void PrismToTet(AOMD::mEntity* );
  ///                          transforms a PYRAMID given by its pointer to TET added to the xMesh (used by tretriseMesh)
  void PyramidToTet(AOMD::mEntity*  );
  ///                          returns the index of the minimun value among a set of values : this set if defined by 
  ///                          a subset of the first vector<int> by giving a vector of indices (second vector<int>)
  int IndexOfMinAmong(vector<int>& , vector<int>&);
  /// Refine the mesh n times (only for 2D triangles)
  xMesh* refineInNewMesh(int n, bool debut=true);
  ///Return if an element is on the boundary
  bool isOnBoundary(AOMD::mEntity* e);
  ///Return a subset with the neighbours vertices
  set<AOMD::mVertex*> getNeighbors(AOMD::mVertex* v);
  /// Return the boundary
  xSubMesh* setBoundary();
  xSubMesh* getBoundary();
private:
  xSubMesh* boundary;

}; 
////////////////////////////////// end class xMesh //////////////////////////////////////////////////


/**
   Mesh as attachable data
*/
class xAttachableMesh : public AOMD::mAttachableData
{
 public:
 xAttachableMesh():m(0){};
 ~xAttachableMesh(){};
  xMesh *m;
};

inline std::ostream& operator<<(std::ostream& s, const Trellis_Util::mPoint& p)
{
  s << p(0) << " " << p(1) << " " << p(2);
  return s;
}


/// fill out with a hard copy of in.
/*!
  link between elements  of in and  out are made using attached Entity with  tags xMesh::has_a_copy_tag, and xMesh::is_the_copy_of_tag.
  by defaults, so that the source mesh is not modifyied at all on exit, the attached entity  to get the copied element on the source mesh are deleted. Setting the flag clean_tag_on_source_mesh to false permit to keep then
*/
 void xCopyMesh(xMesh * in  ,  xMesh * out, bool clean_tag_on_source_mesh = true);

void orientQuad(std::vector<AOMD::mVertex*>& p) ;

/// return the list of entity on the boundary of the support the entity in.
/*! Only works for 2d mesh of simplexes. */
void getSupportBoundary(AOMD::mEntity *in, std::set<AOMD::mEntity *> &supportBoundary );

/// return the constitutive entity of entity e 
/*! d is the level of constitutive entity needed
   getSubentities add things to l
  So make sure it is empty before calling it
  It's done on purpose so don't change this behaviour without changing at least xcut::xPhysSurfByTagging::construct_ */
extern void getSubentities(AOMD::mEntity* e, const int d, std::vector<AOMD::mEntity*> &l) ;
 
} // end of namespace

#endif
