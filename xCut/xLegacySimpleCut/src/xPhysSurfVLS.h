#ifndef XPHYSSURFVLS_H
#define XPHYSSURFVLS_H

#include "xIntegrationRule.h"
#include "xMesh.h"
#include "xPartition.h"
#include "xSupportComponent.h"
#include "xVectorLevelSet.h"

namespace xcut
{
class xPhysSurfByTagging;
typedef std::function<void(AOMD::mEntity *, xfem::xPartition &, xfem::xEntityFilter)> xGetPartition;
class xGetSupport;
struct xVLSTriangleCutData;

void edgeCut(const xtensor::xPoint &p0, const xtensor::xPoint &p1, const xfem::xVectorLevelSetData *v0data,
             const xfem::xVectorLevelSetData *v1data, std::vector<xtensor::xPoint> &cutpoints);

class xVLSEdgeCutData
{
  public:
   size_t nbcuts() const { return vertices.size(); }
   void addVertex(AOMD::mVertex *v) { vertices.push_back(v); }
   AOMD::mVertex *operator()(size_t i) const
   {
      if (i > nbcuts())
      {
         throw;
      }
      else
         return vertices[i];
   }

  private:
   std::vector<AOMD::mVertex *> vertices;
};

typedef std::array<AOMD::mVertex *, 3> triangle;
typedef std::vector<triangle> vtriangle;

struct xVLSTriangleCutData
{
   AOMD::mEntity *fromMeshIn;
   AOMD::mEntity *fromMeshOut;
   vtriangle TrianglesIn;
   vtriangle TrianglesOut;
   void clear()
   {
      fromMeshIn = nullptr;
      fromMeshOut = nullptr;
      TrianglesIn.clear();
      TrianglesOut.clear();
   }
   xVLSTriangleCutData() : fromMeshIn(nullptr), fromMeshOut(nullptr)
   {
      TrianglesIn.reserve(4);
      TrianglesOut.reserve(4);
   }
};

AOMD::mEntity *new_mEntityFromTriangle(const triangle &in, pGEntity classification);

class xVLSTriangleCutAttachableData : public AOMD::mAttachableData
{
  public:
#ifdef WITH_XLEGACYSIMPLECUT
   friend class xPhysSurfVLS;
#endif
  private:
   typedef std::set<AOMD::mEntity *> entityContainer;
   entityContainer TrianglesIn;
   entityContainer TrianglesOut;
};

/*!this function cut a triangle (face)
  The results of the cut is attached to face as a xVLSTriangleCutData pointer with tag face_data_tag.
  It is up to the user to destroy this data, using face->deleteData( face_data_tag )
  The resulting  xVLSTriangleCutData mainly consist of a list of triangle pointing to existing AOMD::mVertex *.
  The needed new AOMD::mVertex are stored on the boundary mesh.
  the Function add to the boundary mesh the edge entity that lies on the boundary between in and out of the vector level set.
  !*/

// void cutFacePhysSurfVLS(AOMD::mEntity * face, xMesh &boundarymesh,
//			  const std::vector<const xVectorLevelSetData * > &nodesVLSdata,
//			  const unsigned int edge_data_tag, const unsigned int face_data_tag, const unsigned int filter_tag,
// std::vector<xcut::xPhysSurfByTagging *> &promotors);

void cutFacePhysSurfVLS(AOMD::mEntity *face, xfem::xMesh &boundarymesh,
                        const std::vector<const xfem::xVectorLevelSetData *> &nodesVLSdata,
                        const std::array<xVLSEdgeCutData *, 3> &edges_data, xVLSTriangleCutData &tcutdata);

class xPhysSurfVLS
{
  public:
   // template constructor
   template <class MESHPHYS>
   xPhysSurfVLS(const xfem::xVectorLevelSet &vls, const MESHPHYS &_mesh, const std::string &name,
                const xfem::xGetSupport &_getsupport,
                std::vector<xcut::xPhysSurfByTagging *> _promotors = std::vector<xcut::xPhysSurfByTagging *>());

   ~xPhysSurfVLS();
   // Notes : upon creation, new vertices are stored and owned by the boundarymesh.
   // New faces are stored as array of pointer to vertex, attached to parent entity via xVLSTriangleCutData (attachable data
   // container. ))
   const xfem::xMesh *getIsoZeroMesh() { return &boundarymesh; }
   const xfem::xSupportComponent &getSupportComponent(AOMD::mEntity *e) const;
   std::function<const xfem::xSupportComponent *(AOMD::mEntity *e)> getSupportComponentFunction() const
   {
      return [this](AOMD::mEntity *e) { return &this->getSupportComponent(e); };
   }
   void clearSupportComponent(AOMD::mEntity *e);
   void clearAllSupportComponent();

   void getPartitionOut(AOMD::mEntity *, xfem::xPartition &, xfem::xEntityFilter);
   void getPartition(AOMD::mEntity *, xfem::xPartition &, xfem::xEntityFilter);
   xfem::xGetPartition partitionOutFunction()
   {
      return std::bind(&xPhysSurfVLS::getPartitionOut, std::ref(*this), std::placeholders::_1, std::placeholders::_2,
                       std::placeholders::_3);
   }

   xfem::xGetPartition partitionFunction()
   {
      return std::bind(&xPhysSurfVLS::getPartition, std::ref(*this), std::placeholders::_1, std::placeholders::_2,
                       std::placeholders::_3);
   }

   xfem::xEntityFilter isOutFilter() const;

   xfem::xEntityFilter isInFilter() const;

   xfem::xEntityFilter isCutFilter() const;

   xfem::xEntityFilter isSupportCoversOutFilter() const;

   xfem::xEntityFilter isSupportCoversInFilter() const;

  private:
   void tagAndPromoteEntity(AOMD::mEntity *e, char tag);
   void unTagAndUnPromoteEntity(AOMD::mEntity *e);
   void treatCutData(AOMD::mEntity *face, xVLSTriangleCutData &fcutdata);
   void clearCutData(AOMD::mEntity *face);
   bool isIn(const AOMD::mEntity &e) const;
   bool isCut(const AOMD::mEntity &e) const;
   bool isOut(const AOMD::mEntity &e) const;
   bool isSupportCoversIn(AOMD::mEntity *) const;
   bool isSupportCoversOut(AOMD::mEntity *) const;
   xfem::xMesh boundarymesh;
   const std::string name;
   template <typename T>
   using datamanager = xinterface::aomd::xAttachedDataManagerAOMD<T>;
   datamanager<char> filter_data;
   datamanager<xVLSTriangleCutAttachableData> faces_data;
   mutable datamanager<xfem::xSupportComponent> supportcomponent_data;
   // unsigned int supportcomponent_tag;

   const xfem::xGetSupport &getsupport;
   unsigned int dim;
   std::list<AOMD::mEntity *> facestoclean;
   mutable std::set<AOMD::mEntity *> knowncomponent;
   std::vector<xcut::xPhysSurfByTagging *> promotors;
};

template <class MESHPHYS>
xPhysSurfVLS::xPhysSurfVLS(const xfem::xVectorLevelSet &vls, const MESHPHYS &mesh, const std::string &_name,
                           const xfem::xGetSupport &_getsupport, std::vector<xcut::xPhysSurfByTagging *> _promotors)
    : name(_name), getsupport(_getsupport), promotors(_promotors)
{
   xinterface::aomd::xAttachedDataManagerAOMD<xVLSEdgeCutData> vlsedgecutdata;
   dim = mesh.size(3) != 0 ? 3 : (mesh.size(2) != 0 ? 2 : ((mesh.size(1) != 0) ? 1 : 0));

   AOMD::mMesh &boundarymmesh = boundarymesh.getMesh();

   pGEntity GEn0 = boundarymmesh.getGEntity(0, 0);
   // pGEntity GEn1 = boundarymesh.getGEntity (1, 1);
   auto itedge = mesh.begin(1);
   auto itedgeend = mesh.end(1);
   int vid = 0;
   std::list<AOMD::mEntity *> edgetoclean;

   //  cout << "Start xPhysSurfVLS " << mesh.size(1) << endl;
   for (; itedge != itedgeend; ++itedge)
   {
      AOMD::mEntity &edge = **itedge;
      const AOMD::mVertex &v0 = *static_cast<AOMD::mVertex *>((edge.get(0, 0)));
      const AOMD::mVertex &v1 = *static_cast<AOMD::mVertex *>((edge.get(0, 1)));
      std::vector<xtensor::xPoint> cutpoints;
      edgeCut(v0.point(), v1.point(), vls(v0), vls(v1), cutpoints);
      if (cutpoints.size())
      {
         xVLSEdgeCutData &vlsedgecutdata_e = vlsedgecutdata.setData(edge);
         for (const xtensor::xPoint &P : cutpoints)
         {
            AOMD::mVertex *v = boundarymmesh.createVertex(vid, P(0), P(1), P(2), GEn0);
            vlsedgecutdata_e.addVertex(v);
            vid++;
         }
      }
   }
   // xexport::xExportGmshAscii pExport;
   if (mesh.dim() == 2)
   {
      xVLSTriangleCutData fcutdata;
      auto itfc = mesh.begin(2);
      auto itfcend = mesh.end(2);
      for (; itfc != itfcend; ++itfc)
      {
         AOMD::mEntity *face = (*itfc);
         int nbnodes = face->size(0);
         if (nbnodes != 3) throw;
         std::vector<const xfem::xVectorLevelSetData *> nodesVLSdata(nbnodes);
         for (int i = 0; i < nbnodes; ++i)
         {
            nodesVLSdata[i] = vls(*((AOMD::mVertex *)face->get(0, i)));
         }
         std::array<xVLSEdgeCutData *, 3> edgesVLSdata = {vlsedgecutdata.getData(*face->get(1, 0)),
                                                          vlsedgecutdata.getData(*face->get(1, 1)),
                                                          vlsedgecutdata.getData(*face->get(1, 2))};
         cutFacePhysSurfVLS(face, boundarymesh, nodesVLSdata, edgesVLSdata, fcutdata);
         treatCutData(face, fcutdata);
         facestoclean.push_back(face);
      }  // end cutface loop
   }
   else  // 3d case not treated.
   {
      throw;
      /*
        typename MESHPHYS::iter ittc =  mesh.begin(3);
        typename MESHPHYS::iter ittcend = mesh.end(3);
        for (;ittc!=ittcend; ++ittc){
        AOMD::mEntity * tet = (*ittc);
        int nbnodes = tet->size(0);
        if (nbnodes !=4) throw;
    std::vector<const xfem::xVectorLevelSetData * > nodesVLSdata(nbnodes);
        for (int i=0; i < nbnodes; ++i){
        nodesVLSdata[i] =  vls( *((AOMD::mVertex *) tet->get(0,i)));
        }
      */
      // cutTetPhysSurfVLS(face, boundarymesh, nodesVLSdata, edge_data_tag, face_data_tag, filter_tag, promotors);
      // facestoclean.push_back(face);
   }
}

/*
  template <class MESH>
  void xPhysSurfVLS::exportGmsh(const MESH &mesh) const{
  xexport::xExportGmshAscii pExport;
  Export(&boundarymesh, pExport, name+"_boundarymesh");
  pExport.openFile(name+"partitions");
  pExport.startView ("DeadZone");
  typename MESH::iter itf = mesh.begin(2);
  typename MESH::iter itfend = mesh.end(2);
  for (;itf!=itfend;++itf){
  AOMD::mEntity *face = *itf;
  xVLSTriangleCutAttachableData *tcutdata = (xVLSTriangleCutAttachableData *)face->getData( face_data_tag);
  if (tcutdata){
  const xVLSTriangleCutAttachableData::entityContainer  & TrianglesOut = tcutdata->getTrianglesIn();
  xVLSTriangleCutAttachableData::entityContainer::const_iterator itsubT = TrianglesOut.begin();
  xVLSTriangleCutAttachableData::entityContainer::const_iterator itsubTend = TrianglesOut.end();
  for (;itsubT != itsubTend; ++itsubT){
  const AOMD::mEntity  * subT = *itsubT;
  pExport.exportTriangle ( ((AOMD::mVertex *) (subT->get(0,0)))->point(),
  ((AOMD::mVertex *) (subT->get(0,1)))->point(),
  ((AOMD::mVertex *) (subT->get(0,2)))->point(),
  0., 0., 0.);
  }
  }
  }
  pExport.endView();

  pExport.startView ("DamagedZone");
  itf = mesh.begin(2);
  for (;itf!=itfend;++itf){
  AOMD::mEntity *face = *itf;
  xVLSTriangleCutAttachableData *tcutdata = (xVLSTriangleCutData *)face->getData( face_data_tag);
  if (tcutdata){
  const xVLSTriangleCutAttachableData::entityContainer  & TrianglesOut = tcutdata->getTrianglesOut();
  xVLSTriangleCutAttachableData::entityContainer::const_iterator itsubT = TrianglesOut.begin();
  xVLSTriangleCutAttachableData::entityContainer::const_iterator itsubTend = TrianglesOut.end();
  for (;itsubT != itsubTend; ++itsubT){
  const  AOMD::mEntity *   subT = *itsubT;
  pExport.exportTriangle ( ((AOMD::mVertex *) (subT->get(0,0)))->point(),
  ((AOMD::mVertex *) (subT->get(0,1)))->point(),
  ((AOMD::mVertex *) (subT->get(0,2)))->point(),
  1., 1., 1.);
  }
  }
  }
  pExport.endView();

  pExport.closeFile();
  }
*/

/*
template <class MESH>
void  xPhysSurfVLS::exportSupportComponentGmsh(const MESH &mesh) {
xexport::xExportGmshAscii pExport;
pExport.openFile(name+"_SupportComponent");
typename MESH::iter it = mesh.begin(0);
typename MESH::iter itend = mesh.end(0);
for (;it!=itend;++it){
  std::stringstream vname;
  vname <<(*it)->getId();
  const xSupportComponent & support_component = getSupportComponent(*it);
  if (support_component.getnbComponents()>1)
    {
      pExport.startView ("support_node_" + vname.str() );
      for (int i =0; i< support_component.getnbComponents(); ++i){
        const xSupportComponent::component & current_component = support_component(i);
        xSupportComponent::component::const_iterator itsubT = current_component.begin();
        xSupportComponent::component::const_iterator itsubTend = current_component.end();
        for (;itsubT != itsubTend; ++itsubT){
          const AOMD::mEntity *   subT = *itsubT;
          pExport.exportTriangle ( ((AOMD::mVertex *) (subT->get(0,0)))->point(),
                                   ((AOMD::mVertex *) (subT->get(0,1)))->point(),
                                   ((AOMD::mVertex *) (subT->get(0,2)))->point(),
                                   i, i, i);
        }
      }
      pExport.endView();
    }
}
it = mesh.begin(1);
itend = mesh.end(1);
for (;it!=itend;++it){
  std::stringstream vname;
  vname <<(*it)->get(0,0)->getId()<<"_"<<(*it)->get(0,1)->getId();
  const xSupportComponent & support_component = getSupportComponent(*it);
  if (support_component.getnbComponents()>1)
    {
      pExport.startView ("support_edge_" + vname.str() );
      for (int i =0; i< support_component.getnbComponents(); ++i){
        const xSupportComponent::component & current_component = support_component(i);
        xSupportComponent::component::const_iterator itsubT = current_component.begin();
        xSupportComponent::component::const_iterator itsubTend = current_component.end();
        for (;itsubT != itsubTend; ++itsubT){
          const AOMD::mEntity *   subT = *itsubT;
          pExport.exportTriangle ( ((AOMD::mVertex *) (subT->get(0,0)))->point(),
                                   ((AOMD::mVertex *) (subT->get(0,1)))->point(),
                                   ((AOMD::mVertex *) (subT->get(0,2)))->point(),
                                   i, i, i);
        }
      }
      pExport.endView();
    }
}
}

*/

}  // end namespace xcut

#endif  // XPHYSSURFVLS_H
