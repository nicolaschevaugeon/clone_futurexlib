/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/
#ifndef _FMSK_meshUtil_
#define _FMSK_meshUtil_

#include "xRegion.h"
#include "meshinterfacexRegion.h"
// Utilities Functions :

#include "xIntegrationRule.h"

namespace xfastmarching
{

  namespace skeleton
  {
    typedef meshinterfacexRegion meshinterface;
    typedef AOMD::mMesh   mesh;
    typedef AOMD::mVertex vertex;
    typedef AOMD::mEdge   edge;
    typedef AOMD::mFace   face;
    typedef AOMD::mEntity entity;


    // A collection of simple function to avoid to work with pointers, and simplify some simple mesh queries
    template < class V >

    Trellis_Util::mPoint makePoint(const V & v){
      return Trellis_Util::mPoint(v(0), v(1), v(2));
    } 

    const vertex & getvertex( const edge &e, int i){
      assert(( i==0 ) || (i == 1));
      return static_cast< const vertex &> (*e.get(0,i));
    }

    const vertex & getvertex( const face &f, int i){
      assert(( i==0 ) || (i == 1) || (i == 2));
      return static_cast< const vertex &> (*f.get(0,i));
    }

    const edge & getedge( const face &f, int i){
      assert(( i==0 ) || (i == 1) || (i == 2));
      return static_cast< const edge &> (*f.get(1,i));
    }

    const vertex & othervertex( const edge &e, const vertex &v){
      if (&getvertex(e,0) == &v) return getvertex(e,1);
      if (&getvertex(e,1) == &v) return getvertex(e,0);
      throw;
    }

    const vertex &othervertex(const face &f, const vertex &v0, const vertex &v1){
      if ((&getvertex(f,0) == &v0)  &&  (&getvertex(f,1) == &v1)) return getvertex(f,2);
      if ((&getvertex(f,0) == &v0)  &&  (&getvertex(f,2) == &v1)) return getvertex(f,1);
  
      if ((&getvertex(f,0) == &v1)  &&  (&getvertex(f,1) == &v0)) return getvertex(f,2);
      if ((&getvertex(f,0) == &v1)  &&  (&getvertex(f,2) == &v0)) return getvertex(f,1);

      if ((&getvertex(f,1) == &v0)  &&  (&getvertex(f,2) == &v1)) return getvertex(f,0);
      if ((&getvertex(f,1) == &v1)  &&  (&getvertex(f,2) == &v0)) return getvertex(f,0);
      throw;
    }

    const edge * getpedge( const vertex &v0, const vertex &v1){
      for(int i = 0; i < v0.size(1); ++i){
	const edge & e = static_cast< const edge &>  (*v0.get(1, i));
	if (&othervertex( e, v0) == &v1 ) return &e;
      }
      return nullptr;
    }



    auto range( const mesh & m, int dim) -> xtool::xRange<decltype(m.begin(dim))> {
      assert (dim >= 0 && dim <=3);
      return xtool::xRange< decltype(m.begin(dim)) > (m.begin(dim), m.end(dim)  );
    }

    Trellis_Util::mPoint getPoint( const edge &e, int i){
      assert( (i == 0) || (i ==1));
      return (static_cast<vertex *> (e.get(0,i)))->point();
    }

    Trellis_Util::mPoint getPoint( const vertex &v){
      return v.point();
    }
    vertex & getVertex ( const edge &e, int i){
      assert( (i == 0) || (i ==1));
      return *(  static_cast< vertex *>((e.get(0,i))));
    }

    vertex & getVertex ( const face &f, int i){
      assert( (i == 0) || (i ==1) || (i==2));
      return *(  static_cast< vertex *>((f.get(0,i))));
    }

    Trellis_Util::mPoint getPoint ( const face &f, int i){
      assert( (i == 0) || (i ==1) || (i==2)); 
      return (static_cast<vertex *> (f.get(0,i)))->point();
    }


    edge & getEdge ( const face &f, int i){
      assert( (i == 0) || (i ==1) || (i==2));
      return *(  static_cast< edge *>((f.get(1,i))));
    }


    std::ostream & operator << ( std::ostream & out, const Trellis_Util::mPoint & point ){
      out << "{" << point(0) << ", " <<  point(1)<< ", "<< point(2)<<"}";
      return out;
    }

    /// This class is meant to construct a mesh embeded in a parent mesh.
    /*!  Cloning entitie from teh main mesh, while keeping the rela tion ship. 
      getGetPartition return a function object that permit to get partition of an element of the parent mesh
      !*/
    class embeded_mesh{
      /* Possible policy : MaintainPartition 0 1 2 3
	 Maintain Parent 
	 Maintain son( a unique entity )  atriangle cut in two triangle, as no son but a partition ...
	 son is a copy of an entity in the mesh ...
	 should be called clone son or something ...
	 PArtition of a node is anode.
	 Partition of an edge is one edge or n edge and n-1 vertex
	 partitoion of a face is one face ore 2face and an edge or ... what ever, but only entity whose parent is pe.
      */
    public:
      face *   clone_face( face & f_parent);
      edge *   clone_edge(edge &e_parent);
      vertex * clone_vertex(vertex &v_parent);
      template <class POINT>
	vertex * createVertex(const POINT &point,  entity *parent, AOMD::pGEntity classif);
  
      edge * createEdge_oneLevel( vertex *v0, vertex* v1 , entity * parent, AOMD::pGEntity classif);
      face * createFaceWithEdges_oneLevel( edge *e0, edge * e1 ,  edge *e2, entity * parent, AOMD::pGEntity classif);

      /// the returned function object permits to obtain a parttion of the parent mesh element if it exist, other wise it return the parent element
      /*!  Note only  Triangle are treated properly for now. */
      xfem::xGetPartition getGetPartition() const;
      bool is_son(const entity &e) const;
      bool is_vertex_son( const vertex &v) const;
    public :   
      xfem::xMesh m;
    private :
      xinterface::aomd::xAttachedDataManagerAOMD < entity *> sonManager;
      xinterface::aomd::xAttachedDataManagerAOMD < entity *> parentManager;
      xinterface::aomd::xAttachedDataManagerAOMD < xfem::xPartition > partitionManager;

    };

    /// embeded_mesh implementation
    face * embeded_mesh::clone_face( face & f_parent){
      auto is_pf_son = sonManager.getData( f_parent);
      if (is_pf_son) return dynamic_cast< face *> ( * is_pf_son);
 
 
      std::array< edge *, 3> pe3_parent ={ (edge *)(f_parent.get(1,0)), (edge *)(f_parent.get(1,1)), (edge *) (f_parent.get(1,2)) };
      std::array< edge *, 3> pe3_son = {  clone_edge(*pe3_parent[0]), clone_edge(*pe3_parent[1]), clone_edge(*pe3_parent[2])};
      face * pf_son = m.createFaceWithEdges_oneLevel(pe3_son[0], pe3_son[1], pe3_son[2], f_parent.getClassification());
      sonManager.setData( f_parent) = pf_son;
      parentManager.setData( *pf_son) = &f_parent;
      return pf_son;
    }

    edge *embeded_mesh::clone_edge(edge &e_parent){
      auto is_pe_son = sonManager.getData( e_parent);
      if(is_pe_son) return dynamic_cast< edge *> ( * is_pe_son);
  
      std::array< vertex *, 2> pv2_parent = { (vertex *) ( e_parent.get(0,0)), (vertex *) (e_parent.get(0,1))};
      std::array< vertex *,  2> pv2_son = { clone_vertex(*pv2_parent[0]), clone_vertex(*pv2_parent[1])};
      edge * pe_son = m.createEdge_oneLevel( pv2_son[0], pv2_son[1], e_parent.getClassification());
      sonManager.setData( e_parent) = pe_son;
      parentManager.setData( *pe_son) = &e_parent;
      return pe_son;
    }

    vertex *embeded_mesh::clone_vertex(vertex &v_parent){
      auto is_pv_son = sonManager.getData( v_parent);
      if (is_pv_son) return dynamic_cast< vertex *> ( * is_pv_son);
      Trellis_Util::mPoint p = v_parent.point();
      vertex *pv_son = m.createVertex(p(0), p(1), p(2), v_parent.getClassification());
      sonManager.setData( v_parent) = pv_son;
      parentManager.setData( *pv_son) = &v_parent;
      return pv_son;
    }

    template <class POINT>
      vertex * embeded_mesh::createVertex(const POINT &point,  entity *parent, AOMD::pGEntity classif){
      auto psonv = m.createVertex(point(0), point(1), point(2), classif);
      parentManager.setData(*psonv ) = parent;
      return psonv;
    }


    edge * embeded_mesh::createEdge_oneLevel( vertex *v0, vertex* v1 , entity * parent, AOMD::pGEntity classif){
      auto  e01 = const_cast<edge *>(getpedge( *v0, *v1));
      if(!e01) e01 = m.createEdge_oneLevel(v0, v1, classif);
      parentManager.setData(*e01 ) = parent;
      return e01;
    }

    face * embeded_mesh::createFaceWithEdges_oneLevel( edge *e0, edge * e1 ,  edge *e2, entity * parent, AOMD::pGEntity classif){
      auto psonface = m.createFaceWithEdges_oneLevel(e0, e1, e2,  classif);
      if (parent->getLevel() == 2){
	partitionManager.setData( *parent).insert(psonface);
      }
      parentManager.setData(*psonface) = parent;
      return psonface;
    }

    xfem::xGetPartition embeded_mesh::getGetPartition() const{
      return [this](entity* e, xfem::xPartition& part,   xfem::xEntityFilter filter){
	auto is_partition = partitionManager.getData(*e);
	if(is_partition) {
	  for(auto ppe : (*is_partition))  if (filter (ppe) ) part.insert(ppe);
	  return;
	}
	if(filter( e ) ) part.insert(e);
	return;
      };
    }

    bool embeded_mesh::is_son(const entity &e) const{
      if (parentManager.getData(const_cast<entity &> (e) )) return true;
      return false;
    }

    bool embeded_mesh::is_vertex_son( const vertex &v) const{
      auto p_pv = parentManager.getData(const_cast<vertex &> (v) );
      if (p_pv){
	if (dynamic_cast< vertex * >( *p_pv)) return true;
      }
      return false;
    }


    ////
    /// This class is meant to construct a mesh embeded in a parent mesh.
    /*!  Cloning entitie from teh main mesh, while keeping the rela tion ship. 
      getGetPartition return a function object that permit to get partition of an element of the parent mesh
      !*/
    class embeded_mesh_fit_to_vertex{
      /* Possible policy : MaintainPartition 0 1 2 3
	 Maintain Parent 
	 Maintain son( a unique entity )  atriangle cut in two triangle, as no son but a partition ...
	 son is a copy of an entity in the mesh ...
	 should be called clone son or something ...
	 PArtition of a node is anode.
	 Partition of an edge is one edge or n edge and n-1 vertex
	 partition of a face is one face ore 2face and an edge or ... what ever, but only entity whose parent is pe.
      */
    public:
    embeded_mesh_fit_to_vertex(double _eps):eps(_eps){}
      /// Clone create a face that is a copy of a face supposed to be in the parent mesh and return it
      /*!  if the face already exist in the son mesh, then the input face knows it thanks to the son manager
	and the function just return the pointer to the son face.
	if a new face is really created it is classifyied on the same GENtity than the parent one.
	To create the face, the function use the edges copies of the edge of the parent face if they already exist, other wise it clone       them, using clone edge.
	!*/
      face *   clone_face( face & f_parent);
      edge *   clone_edge(edge &e_parent);
      vertex * clone_vertex(vertex &v_parent);
 
  
      // create or return a vertex at position s on the parent edge in the son mesh 
      /*! 
	the parent of this node is the edge if s is not to clode to one of the node of the parent.
	if s would create a vertex too close to the node of the parent edge (according to fit tol), 
	the vertex returned would be at the position of one of the node, and the parent would be one of the node of the parent edge.
	if classif == nullptr (default) the vertex get the classification of it's parent
	other it get the prescribed classification.
	the parent edge as the created vertex as a son. 
	if a son for the edge already exist, the function throw for now.
	!*/
      vertex * createVertex(const double & ss,  edge *parent, AOMD::pGEntity classif = nullptr);
      vertex * createVertex(const double &u, const double &v, face *parent, AOMD::pGEntity classif= nullptr);
  
      edge * createEdge_oneLevel( vertex *v0, vertex* v1 , entity * parent, AOMD::pGEntity classif = nullptr);
      face * createFaceWithEdges_oneLevel( edge *e0, edge * e1 ,  edge *e2, entity * parent, AOMD::pGEntity classif = nullptr);
 
      /// the returned function object permits to obtain a parttion of the parent mesh element if it exist, other wise it return the parent element
      /*!  Note only  Triangle are treated properly for now. */
      xfem::xGetPartition getGetPartition() const;
      bool is_son(const entity &e) const;
      bool is_vertex_son( const vertex &v) const;
      /// e being supposed to be an entity of the son mesh, the function return a reference to it's parent entity.
      /*! 
	If none are found, it means that the entity was not part of son mesh and the function throw 
	!*/
      const entity & get_parent( const entity &e) const;
  
    public :   
      xfem::xMesh m;
    private :
      xinterface::aomd::xAttachedDataManagerAOMD < entity *> sonManager;
      xinterface::aomd::xAttachedDataManagerAOMD < entity *> parentManager;
      xinterface::aomd::xAttachedDataManagerAOMD < xfem::xPartition > partitionManager;
      const double eps;

    };

    /// embeded_mesh_fit_to_vertex implementation
    face * embeded_mesh_fit_to_vertex::clone_face( face & f_parent){
      auto is_pf_son = sonManager.getData( f_parent);
      if (is_pf_son) return dynamic_cast< face *> ( * is_pf_son);
 
 
      std::array< edge *, 3> pe3_parent ={ (edge *)(f_parent.get(1,0)), (edge *)(f_parent.get(1,1)), (edge *) (f_parent.get(1,2)) };
      std::array< edge *, 3> pe3_son = {  clone_edge(*pe3_parent[0]), clone_edge(*pe3_parent[1]), clone_edge(*pe3_parent[2])};
      face * pf_son = m.createFaceWithEdges_oneLevel(pe3_son[0], pe3_son[1], pe3_son[2], f_parent.getClassification());
      sonManager.setData( f_parent) = pf_son;
      parentManager.setData( *pf_son) = &f_parent;
      return pf_son;
    }

    edge *embeded_mesh_fit_to_vertex::clone_edge(edge &e_parent){
      auto is_pe_son = sonManager.getData( e_parent);
      if(is_pe_son) return dynamic_cast< edge *> ( * is_pe_son);
  
      std::array< vertex *, 2> pv2_parent = { (vertex *) ( e_parent.get(0,0)), (vertex *) (e_parent.get(0,1))};
      std::array< vertex *,  2> pv2_son = { clone_vertex(*pv2_parent[0]), clone_vertex(*pv2_parent[1])};
      edge * pe_son = m.createEdge_oneLevel( pv2_son[0], pv2_son[1], e_parent.getClassification());
      sonManager.setData( e_parent) = pe_son;
      parentManager.setData( *pe_son) = &e_parent;
      return pe_son;
    }

    vertex *embeded_mesh_fit_to_vertex::clone_vertex(vertex &v_parent){
      auto is_pv_son = sonManager.getData( v_parent);
      if (is_pv_son) return dynamic_cast< vertex *> ( * is_pv_son);
      Trellis_Util::mPoint p = v_parent.point();
      vertex *pv_son = m.createVertex(p(0), p(1), p(2), v_parent.getClassification());
      sonManager.setData( v_parent) = pv_son;
      parentManager.setData( *pv_son) = &v_parent;
      return pv_son;
    }

    vertex * embeded_mesh_fit_to_vertex::createVertex(const double &ss,  edge *parent, AOMD::pGEntity classif){
      vertex * psonv = nullptr;
      if (ss <= eps)      psonv = clone_vertex( getVertex(*parent,0) );
      else if (ss >= (1-eps) ) psonv = clone_vertex( getVertex(*parent,1) );
      else{
	auto p = getPoint(*parent,0)*(1-ss) + getPoint(*parent, 1)*(ss);
	psonv = m.createVertex(p(0), p(1), p(2), parent->getClassification());
	if (sonManager.getData( *parent)) throw;
	sonManager.setData( *parent) = psonv;
	parentManager.setData( *psonv) = parent;
      }
      if (classif) psonv->classify(classif);
      return psonv;
    }
 
    vertex * embeded_mesh_fit_to_vertex::createVertex(const double &u, const double &v, face *parent, AOMD::pGEntity classif){
      vertex * psonv = nullptr;
      double L0 = 1.-u-v;
      double L1 = u;
      double L2 = v;
      if ((L1<= eps) && (L2 <= eps))
	psonv = clone_vertex( getVertex( *parent, 0) );
      else if ( (L2 <= eps) && (L0 <=eps) )
	psonv = clone_vertex( getVertex( *parent, 1) );
      else if ( (L0 <= eps) && (L1 <=eps) )
	psonv = clone_vertex( getVertex( *parent, 2) );
      else if ( (L2 <= eps) ){
	auto ppsonv = sonManager.getData( getEdge ( *parent, 0));
	if(!ppsonv) throw;
	psonv = dynamic_cast< vertex *> (*ppsonv);
      }
      else if ( (L0 <= eps) ){
	auto ppsonv = sonManager.getData( getEdge ( *parent, 1));
	if(!ppsonv) throw;
	psonv = dynamic_cast< vertex *> (*ppsonv);
      }
      else if ( (L1 <= eps) ){
	auto ppsonv = sonManager.getData( getEdge ( *parent, 2));
	if(!ppsonv) throw;
	psonv = dynamic_cast< vertex *> (*ppsonv);
    
      }
      else {
	if (sonManager.getData( *parent)) throw;
	auto p = getPoint(*parent, 0)*L0 + getPoint(*parent, 1)*L1 + getPoint(*parent, 2)*L2;
	psonv = m.createVertex(p(0), p(1), p(2), parent->getClassification());
	sonManager.setData( *parent) = psonv;
	parentManager.setData( *psonv) = parent;
      }
      if (classif) psonv->classify(classif);
      return psonv;
    }

    edge * embeded_mesh_fit_to_vertex::createEdge_oneLevel( vertex *v0, vertex* v1 , entity * parent, AOMD::pGEntity classif){
      if (v0 == v1) return nullptr;
      auto  e01 = const_cast<edge *>(getpedge( *v0, *v1));
      if(!e01) {
	e01 = m.createEdge_oneLevel(v0, v1, parent->getClassification());
	parentManager.setData(*e01 ) = parent;
      }
      if (classif)  e01->classify(classif);
      return e01;
    }




    face * embeded_mesh_fit_to_vertex::createFaceWithEdges_oneLevel( edge *e0, edge * e1 ,  edge *e2, entity * parent, AOMD::pGEntity classif){
      if (!e0 || !e1 || !e2) return nullptr;
      if ((e0==e1) || (e1==e2) || (e2==e0)) return nullptr;
      auto psonface = m.createFaceWithEdges_oneLevel(e0, e1, e2, parent->getClassification() );
      if (classif)  psonface->classify(classif);
      if (parent->getLevel() == 2){
	partitionManager.setData( *parent).insert(psonface);
      }
      parentManager.setData(*psonface) = parent;
      return psonface;
    }



    xfem::xGetPartition embeded_mesh_fit_to_vertex::getGetPartition() const{
      return [this](entity* e, xfem::xPartition& part,   xfem::xEntityFilter filter){
	auto is_partition = partitionManager.getData(*e);
	if(is_partition) {
	  for(auto ppe : (*is_partition))  if (filter (ppe) ) part.insert(ppe);
	  return;
	}
	if(filter( e ) ) part.insert(e);
	return;
      };
    }

    bool embeded_mesh_fit_to_vertex::is_son(const entity &e) const{
      if (parentManager.getData(const_cast<entity &> (e) )) return true;
      return false;
    }

    bool embeded_mesh_fit_to_vertex::is_vertex_son( const vertex &v) const{
      auto p_pv = parentManager.getData(const_cast<vertex &> (v) );
      if (p_pv){
	if (dynamic_cast< vertex * >( *p_pv)) return true;
      }
      return false;
    }
 
    const entity & embeded_mesh_fit_to_vertex::get_parent( const entity &e) const{
      auto p_pe = parentManager.getData(const_cast<entity &> (e) );
      if (p_pe){
	return *(*p_pe);
      }
      else throw;
    }
    ///// 
    /*
      vertex * get_parent_vertex(vertex &v) {
      auto p_pv = parentManager.getData(const_cast<vertex &> (v) );
      return dynamic_cast< vertex * >( *p_pv);
      }
      edge * get_parent_edge(vertex &v) {
      auto p_pv = parentManager.getData(const_cast<vertex &> (v) );
      return dynamic_cast< edge * >( *p_pv);
      }
      face* get_parent_face(vertex &v) {
      auto p_pv = parentManager.getData(const_cast<vertex &> (v) );
      return dynamic_cast< face * >( *p_pv);
      }
    */
 
    /// My own version of export ... Multiple view in the same file, related to the same mesh ...
    class exportGMSHMine{
    public:
      exportGMSHMine(const mesh &m,  const meshinterface &mi,   const std::string &filename_base, int shift = 1);
   
      void appendVertexScalarView( const std::string &viewname ,
				   std::function< double ( const vertex & v ) >   vvalue,
				   int timestep =0 , double timevalue =0.);
   

      void appendVertexVectorView(const std::string &viewname ,

					   std::function< std::array< double , 3> ( const vertex & v ) >   vvalue,
					   int timestep = 0 , double timevalue=0. );

      void appendTriangleScalarView(const std::string &viewname,
				    std::function< double ( const face & f ) >   fvalue,
				    int timestep =0 , double timevalue =0.);
      void appendEdgeScalarView(const std::string &viewname,
				std::function< double ( const edge & e ) >   evalue,
				int timestep = 0, double timevalue = 0.);
      void appendEdgeVectorView(const std::string &viewname,
				std::function< std::array< double , 3>  ( const edge & e ) >   evalue,
				int timestep = 0, double timevalue = 0.);
    
    private:
      const mesh &m;
      std::string outfilename;
      // xinterface::aomd::xAttachedDataManagerAOMD< int   >  globalverticesnumbering;
  
      entitystorage< meshinterface, meshinterface::entity, int > entitynumbering;
      int shift;
    };

    exportGMSHMine::exportGMSHMine(const mesh &_m, const meshinterface &mi,
				   // const partitionManager &_part_man,
				   const std::string &filename_base, int _shift):m(_m),  entitynumbering(mi),
      //part_man(_part_man),
      shift(_shift){
	//int mpi_size, mpi_rank;
	//MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	//MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
   
	// globalverticesnumbering = globalVerticesNumbering(const_cast< mesh &>(m), part_man );
	int mpi_rank = 0;
	//  outfilename = filename_base + "_NP"+std::to_string(mpi_size)+"_P"+std::to_string(mpi_rank)+".msh";
	outfilename = filename_base +".msh";
   
	std::ofstream out(outfilename.c_str() );
	out << "$MeshFormat" << std::endl;
	out << "2.2 0 "<< sizeof(double) << std::endl;
	out << "$EndMeshFormat" << std::endl;
	out << "$Nodes" << std::endl;
	out << m.size(0)<< std::endl;
	int nbelem = 0;
	// int i = 0;
	for(auto pv : range(m,0 ) ){
	  const vertex &v = static_cast< const vertex & > (*pv);
	  const Trellis_Util::mPoint p = v.point();
	  // entitynumbering.setVal(v, i++);
	  // out << *globalverticesnumbering.getData(v)+shift  << " " << p(0) << " " << p(1) << " " << p(2) << "\n";
	  out << v.getId()+shift  << " " << p(0) << " " << p(1) << " " << p(2) << "\n";
	  //out << entitynumbering.getVal(v)  +shift  << " " << p(0) << " " << p(1) << " " << p(2) << "\n";
     
	  if ( (!v.getClassification()) || v.getClassification()->dim() == 0 ) {
	    ++nbelem;
	    entitynumbering.set(v, nbelem);
	  }
	  //else if(part_man.hasRemoteCopy(*pv) ) ++nbelem;
	}
	out << "$EndNodes" << std::endl;
   
	for (int dim = 1; dim <= 3; ++dim){
	  for(auto pe : range(m,dim ) ){
	    const entity &e = static_cast< entity & > (*pe);
	    auto classif = e.getClassification();
	    if (dim == 2)
	      entitynumbering.set(e,  ++nbelem);
	    else if (classif)
	      if (classif->dim() == dim)
		//if (e.getClassification()->dim() == dim ) {
		entitynumbering.set(e,  ++nbelem);
	    //}
	    //else if(part_man.hasRemoteCopy(*pe) ) ++nbelem;
	  }
	}
   
	// nbelem +=  m.size(3);
	out << "$Elements" << std::endl;
	out << nbelem << std::endl;
	// int k = shift;
	for(auto pv : range(m,0 ) ){
	  vertex &v = static_cast< vertex& > (*pv);
	  if ( (!v.getClassification()) || (v.getClassification()->dim() == 0) ){// || part_man.hasRemoteCopy(vertex) ) {
	    int tag = 22;
	    if (v.getClassification())
	      tag = v.getClassification()->tag();
	    //const int id0 =   *globalverticesnumbering.getData( vertex) +shift;
	    const int id0 =   v.getId() +shift;
	    out << entitynumbering(v) << " 15 4 " << tag << " " << tag << " 1 " << mpi_rank+1<< " "  << id0 << "\n"; 
	  }
	}
	for(auto pe : range(m,1 ) ){
	  entity &edge = *pe;
	  auto classif = edge.getClassification();
	  if (classif)
	    if (edge.getClassification()->dim() == 1)
	      //if((edge.getClassification()->dim() == 1)){// || part_man.hasRemoteCopy(edge) ){
	      {
		int tag = 22;
		if ( edge.getClassification()) tag = edge.getClassification()->tag();
		// assert(globalverticesnumbering.getData( *edge.get(0,0)));
		//assert(globalverticesnumbering.getData( *edge.get(0,1)));
		//const int id0 =  *globalverticesnumbering.getData( *edge.get(0,0)) +shift;
		//const int id1 =  *globalverticesnumbering.getData( *edge.get(0,1)) +shift;
		const int id0 =  ( *edge.get(0,0)).getId() +shift;
		const int id1 =  ( *edge.get(0,1)).getId() +shift;
		out <<entitynumbering(edge)  << " 1 4 " << tag << " " << tag << " 1 " << mpi_rank+1 << " " << id0 << " " << id1<< "\n"; 
	      }
	}
	for(auto pf : range(m,2 ) ){
	  entity &tri =  *pf; 
	  //if((tri.getClassification()->dim() == 2) ){//|| part_man.hasRemoteCopy(tri) ){
	  {
	    int tag = 22;
	    if (tri.getClassification()) tag = tri.getClassification()->tag();
	    // assert( globalverticesnumbering.getData( *tri.get(0,0)));
	    //assert( globalverticesnumbering.getData( *tri.get(0,1)));
	    //assert( globalverticesnumbering.getData( *tri.get(0,2)));
	 
	    //       const int id0 = *globalverticesnumbering.getData( *tri.get(0,0)) +shift;
	    //const int id1 = *globalverticesnumbering.getData( *tri.get(0,1)) +shift;
	    //const int id2 = *globalverticesnumbering.getData( *tri.get(0,2)) +shift;

	    const int id0 =( *tri.get(0,0)).getId() +shift;
	    const int id1 = ( *tri.get(0,1)).getId() +shift;
	    const int id2 = ( *tri.get(0,2)).getId() +shift;
     
	    out << entitynumbering(tri) << " 2 4 "  << tag << " " << tag << " 1 " << mpi_rank+1<< " " << id0 << " " << id1 << " " << id2 << "\n"; 
	  }
	}
	for(auto pr : range(m,3 ) ){
	  entity &tet =  *pr;
	  int tag = 22;
	  if (tet.getClassification()) tag = tet.getClassification()->tag();
	  //assert( globalverticesnumbering.getData( *tet.get(0,0)));
	  //assert( globalverticesnumbering.getData( *tet.get(0,1)));
	  //assert( globalverticesnumbering.getData( *tet.get(0,2)));
	  //assert( globalverticesnumbering.getData( *tet.get(0,3)));
	  //const int id0 = *globalverticesnumbering.getData( *tet.get(0,0)) +shift;
	  //const int id1 = *globalverticesnumbering.getData( *tet.get(0,1)) +shift;
	  //const int id2 = *globalverticesnumbering.getData( *tet.get(0,2)) +shift;
	  //const int id3 = *globalverticesnumbering.getData( *tet.get(0,3)) +shift;
	  const int id0 = ( *tet.get(0,0)).getId() +shift;
	  const int id1 = ( *tet.get(0,1)).getId() +shift;
	  const int id2 = ( *tet.get(0,2)).getId() +shift;
	  const int id3 = ( *tet.get(0,3)).getId() +shift;
     
	  out <<entitynumbering(tet)  << " 4 4 "  << tag << " " << tag << " 1 " << mpi_rank+1<< " " << id0 << " " << id1 << " " << id2 << " " << id3 << "\n"; 
	}
	out << "$EndElements" << std::endl;
	out.close();
   
      }


    void exportGMSHMine::appendVertexScalarView(const std::string &viewname ,
						std::function< double ( const vertex & v ) >   vvalue,
						int timestep , double timevalue ){
      //int mpi_size, mpi_rank;
      //MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
      //MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      // int mpi_rank = 0;
      std::ifstream in(outfilename.c_str() );
      if(!in.is_open()) {
	std::cout << "file " << outfilename << " Does not exist " << std::endl;
	return;
      }
      in.close();
      std::ofstream out(outfilename.c_str(), std::ofstream::app);
      out << "$NodeData" << std::endl;
      out << "1" << std::endl;
      out << "\""<< viewname << "\"" << std::endl;
      out <<"1" << std::endl;
      out <<timevalue << std::endl;
      out<< "3" << std::endl;
      out << timestep << std::endl;
      out << "1"<< std::endl;
      out << m.size(0) << std::endl;
      for (auto it = m.begin(0); it !=m.end(0); ++it){
	const vertex &v = static_cast< const vertex & > (*(*it));
	const double val = vvalue(v);
    
	//     out << *globalverticesnumbering.getData(v)+shift << " " << val << std::endl;
	out << v.getId()+shift << " " << val << std::endl;
      }
      out << "$EndNodeData" << std::endl;
    }
 

    void exportGMSHMine::appendVertexVectorView(const std::string &viewname ,							 std::function< std::array< double, 3> ( const vertex & v ) >   vvalue,
                 int timestep , double timevalue ){
      // int mpi_size, mpi_rank;
      //MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
      //MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      // int mpi_rank = 0;
      std::ifstream in(outfilename.c_str() );
      if(!in.is_open()) {
	std::cout << "file " << outfilename << " Does not exist " << std::endl;
	return;
      }
      in.close();
      std::ofstream out(outfilename.c_str(), std::ofstream::app);
      out << "$NodeData" << std::endl;
      out << "1" << std::endl; // number of string tag
      out << "\""<< viewname << "\"" << std::endl;
      out <<"1" << std::endl; // number of real tag
      out <<timevalue << std::endl;
      out<< "3" << std::endl; // number of interger tag
      out << timestep << std::endl; //integer tag 1 -> timestep
      out << "3"<< std::endl;      // integer tag 2 -> number of component per node : 1 3 9
      out << m.size(0) << std::endl; // integer tag 3 -> number of nodevalues
      for (auto it = m.begin(0); it !=m.end(0); ++it){
	const vertex &v = static_cast< const vertex & > (*(*it));
	const std::array< double, 3 > val = vvalue(v);
    
	//     out << *globalverticesnumbering.getData(v)+shift << " " << val[0]<< " " << val[1]<<  " "<< val[2] << std::endl;
	out << v.getId()+shift << " " << val[0]<< " " << val[1]<<  " "<< val[2] << std::endl;
      }
      out << "$EndNodeData" << std::endl;
    }


    void exportGMSHMine::appendTriangleScalarView(const std::string &viewname,
						  std::function< double ( const face & f ) >   fvalue,
						  int timestep, double timevalue){
      //int mpi_rank = 0;
      std::ifstream in(outfilename.c_str() );
      if(!in.is_open()) {
	std::cout << "file " << outfilename << " Does not exist " << std::endl;
	return;
      }
      in.close();
      std::ofstream out(outfilename.c_str(), std::ofstream::app);
      out << "$ElementData" << std::endl;
      out << "1" << std::endl;
      out << "\""<< viewname << "\"" << std::endl;
      out <<"1" << std::endl;
      out <<timevalue << std::endl;
      out<< "3" << std::endl;
      out << timestep << std::endl;
      out << "1"<< std::endl;
      out << m.size(2) << std::endl;
      for (auto it = m.begin(2); it !=m.end(2); ++it){
	const face &f = static_cast< const face & > (*(*it));
	const double val = fvalue(f);
	//     out << *globalverticesnumbering.getData(v)+shift << " " << val << std::endl;
	out << entitynumbering(f)  << " " << val << std::endl;
      }
      out << "$EndNodeData" << std::endl;
    }
    void exportGMSHMine::appendEdgeScalarView(const std::string &viewname,
					      std::function< double ( const edge & e ) >   evalue,
					      int timestep, double timevalue){
      //int mpi_rank = 0;
      std::ifstream in(outfilename.c_str() );
      if(!in.is_open()) {
	std::cout << "file " << outfilename << " Does not exist " << std::endl;
	return;
      }
      in.close();
      std::ofstream out(outfilename.c_str(), std::ofstream::app);
      out << "$ElementData" << std::endl;
      out << "1" << std::endl;
      out << "\""<< viewname << "\"" << std::endl;
      out <<"1" << std::endl;
      out <<timevalue << std::endl;
      out<< "3" << std::endl;
      out << timestep << std::endl;
      out << "1"<< std::endl;
      out << m.size(1) << std::endl;
      for (auto it = m.begin(1); it !=m.end(1); ++it){
	const edge &e = static_cast< const edge & > (*(*it));
	const double val = evalue(e);
	//     out << *globalverticesnumbering.getData(v)+shift << " " << val << std::endl;
	out << entitynumbering(e)  << " " << val << std::endl;
      }
      out << "$EndNodeData" << std::endl;
    }

    void exportGMSHMine::appendEdgeVectorView(const std::string &viewname,
					      std::function< std::array< double, 3>  ( const edge & e ) >   evalue,
					      int timestep, double timevalue){
      //int mpi_rank = 0;
      std::ifstream in(outfilename.c_str() );
      if(!in.is_open()) {
	std::cout << "file " << outfilename << " Does not exist " << std::endl;
	return;
      }
      in.close();
      std::ofstream out(outfilename.c_str(), std::ofstream::app);
      out << "$ElementData" << std::endl;
      out << "1" << std::endl;
      out << "\""<< viewname << "\"" << std::endl;
      out <<"1" << std::endl;
      out <<timevalue << std::endl;
      out<< "3" << std::endl;
      out << timestep << std::endl;
      out << "3"<< std::endl;
      out << m.size(1) << std::endl;
      for (auto it = m.begin(1); it !=m.end(1); ++it){
	const edge &e = static_cast< const edge & > (*(*it));
	const std::array< double, 3> val = evalue(e);
	//     out << *globalverticesnumbering.getData(v)+shift << " " << val << std::endl;
	out << entitynumbering(e)  << " " << val[0] << " "<< val[1]<< " "<< val[2] << std::endl;
      }
      out << "$EndNodeData" << std::endl;
    }

    /*double position_on_edge(const Trellis_Util::mPoint &p0, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p){
      point3d p3d0{p0(0), p0(1), p0(2)};
      point3d p3d1{p1(0), p1(1), p1(2)};
      point3d p3d {p(0), p(1), p(2)};
      auto v01 = fromto(p3d0, p3d1);
      double s = scal( v01, fromto(p3d0, p3d)  )/ scal( v01,v01);
      return s;
      }*/
 
  }//end namespace skeleton

}//end namespace  



#endif
