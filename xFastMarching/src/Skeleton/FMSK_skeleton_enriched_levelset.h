/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/
#ifndef _FMSK_skeleton_enriched_levelset_h_
#define _FMSK_skeleton_enriched_levelset_h_
// MESHINTERFACE must have a vertex type.
// VERTEXSCALEVAL must have a return_type representing a scalar real nuber
// VERTEXSCALEVAL must have return_type operator(const MESHINTERFACE::vertex &) const
// VERTEXRANGE must have begin() and end() that return iterator to pointer to vertex 



#include "xEval.h"
#include "FMSK_meshutil.h"
#include "xRegion.h"

namespace xfastmarching
{

  namespace skeleton
  {

    typedef vector3d <double > vec3;
    typedef entitystorage< meshinterface, meshinterface::vertex, double  > lstype;
    typedef entitystorage< meshinterface, meshinterface::vertex, vector3d<double > > vertexvectorstorage;

    std::pair<vector2d<double >  , double > triple_point2 ( vec3 X0, vec3 X1, vec3 X2, double l0, double l1, double l2, vec3 g0, vec3 g1, vec3 g2);


    template <class MESHINTERFACE, class VERTEXSCALEVAL, class VERTEXANYEVAL >
      class skeleton_enriched_ls{
    public:
      using vertex_t          = typename MESHINTERFACE::vertex;
      using scal_t            = typename VERTEXSCALEVAL::result_type;
      using trans_val_t       = typename VERTEXANYEVAL::result_type;
      using vect_t            = vector3d<scal_t>;
      using ls_t              = entitystorage<  MESHINTERFACE, vertex_t,  scal_t >;
      using gls_t             = entitystorage<  MESHINTERFACE, vertex_t,  vect_t >;
      using trans_t           = entitystorage<  MESHINTERFACE, vertex_t,  trans_val_t > ;
      template < class VERTEXRANGEK, class VERTEXRANGET>
	skeleton_enriched_ls(const MESHINTERFACE & _FM_mi, const  VERTEXRANGEK &known_vertex_range, const VERTEXRANGET & trial_vertex_range, const VERTEXSCALEVAL & inputls, const VERTEXANYEVAL & inputtotransport, scal_t shock_treshold, scal_t fit_tol ):enriched_mesh(fit_tol), enriched_mesh_r(&enriched_mesh.m), mi_enr(enriched_mesh_r),  ls0( _FM_mi),
	lsr(mi_enr),
	gls0(_FM_mi), trans0(_FM_mi),
	updater(_FM_mi, ls0, [](const vertex_t &v){return 1.;}, 0., shock_treshold,   gls0){
    
	for(auto pv : known_vertex_range){
	  ls0.set(*pv, inputls(*pv));
	  trans0.set(*pv, inputtotransport(*pv));
	}
	for(auto pv : trial_vertex_range){
	  ls0.set(*pv, inputls(*pv));
	  trans0.set(*pv, inputtotransport(*pv));
	}
	std::list< const vertex_t * > trial_were_known;
	std::function< double ( const vertex_t  &)  > F = [](const vertex_t &v){return 1.;};
   
	fmeik_intern2( known_vertex_range.begin(), known_vertex_range.end(), trial_vertex_range.begin(), trial_vertex_range.end(),    std::back_inserter(trial_were_known), updater );
	build_enriched_mesh_v2( _FM_mi );

	// std::cout << "AFTER FM " << ls0.getSupport().size() << std::endl;
	//std::cout << "AFTER FM " << lsr.getSupport().size() << std::endl;
    
    
      }

  
      std::function<double ( const AOMD::mVertex &v ) >getlsvertex() const{
	return  [this ](const AOMD::mVertex &v ){
	  //std::cout << "BBB " <<std::endl;
	  double val =10000.; this->ls0.get(v, val); return val;
	};
      }
      const xfem::xMesh & get_enriched_mesh() const{ return enriched_mesh.m;}
      const embeded_mesh_fit_to_vertex & get_enriched_embeded_mesh() const{ return enriched_mesh;}
      class xEvalEnrLs : public  xfem::xEval<double > {
      public:
      xEvalEnrLs(const skeleton_enriched_ls &_enrls): enrls(_enrls){};
	void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, result_type& res) const{
	  auto pe_integ = geo_integ->getEntity();
	  face &e_integ = *(static_cast< face *>( pe_integ));
	  auto * ls_integ = &enrls.ls0;
	  if  (enrls.get_enriched_embeded_mesh().is_son(*pe_integ) ) {
	    //res = 0.; return;
	    ls_integ =  &enrls.lsr;
	  }
	  const vertex & v0 = getvertex(e_integ,0); 
	  const vertex & v1 = getvertex(e_integ,1);
	  const vertex & v2 = getvertex(e_integ,2);
          double vls0(0.), vls1(0.), vls2(0.);
	  ls_integ->get(v0, vls0);
	  ls_integ->get(v1, vls1);
	  ls_integ->get(v2, vls2);
	  auto uvw = geo_integ->getUVW();
	  double u = uvw(0);
	  double v = uvw(1);
	  res = vls0*(1.-u-v) + vls1*u+vls2*v;
	  return;
	}
      private:
	const skeleton_enriched_ls & enrls;
      };
  
      xEvalEnrLs getEvalEnrLs() const{
	return  xEvalEnrLs(*this);
      }


    private :
      /*void build_enriched_mesh(const MESHINTERFACE & mi ) {
	std::stringstream tag_name_ske;
	tag_name_ske <<"skeleton_edge_tag"<< this; 
	tag_skeleton_edge  =  getNewTag(mi_enr, tag_name_ske.str()); 
	auto pskeleton_geo =  enriched_mesh.m.getGEntity(tag_skeleton_edge, 1);
	xinterface::aomd::xAttachedDataManagerAOMD <  AOMD::mEntity *> skManager;
	const double eps = 1.e-2;
	std::function < bool (double )> is_close_to_vertex = [&eps]( double s ){ assert(( s >= 0.) && (s <= 1.)); if (s<= eps) return true; if (s >= (1.-eps)) return true; return false; };
	for (auto pe : getEntityRange<1>(mi) ){
	AOMD::mEdge & e = dynamic_cast<AOMD::mEdge &>(*pe);
	if (!skManager.getData(e) ) {
	auto &v0 = dynamic_cast<const AOMD::mVertex & >(*(e.get(0,0)));
	auto &v1 = dynamic_cast<const AOMD::mVertex & >(*(e.get(0,1)));
	vector3d<double> p0,p1,g0,g1;
	getCoord(mi, v0, p0);
	getCoord(mi, v1, p1);
	vector3d<double > x01 = p1-p0;
	bool hasg0 = gls0.get(v0, g0);
	bool hasg1 = gls0.get(v1, g1);
	if (hasg0 && hasg1){
	double l0 = ls0(v0);
	double l1 = ls0(v1);
	double ss, vs;
	bool iss =  updater.compute_shockpoint2(x01, l0, l1, g0, g1, ss, vs);
	if (iss){
	if (is_close_to_vertex(ss)) std::cout << "CLOSE" << std::endl;
	auto p = p0*(1-ss) + p1*(ss);
	auto pv = enriched_mesh.createVertex(p, pe, pskeleton_geo);
	skManager.setData(e) = pv;
	lsr.set(*pv, vs);
	continue;
	}
	if ( (!iss) && (updater.isShockNeeded(x01, g0, g1 )) ){
	bool found = false;
	std::cout << " ARRRRRG " << std::endl;
	for (int i =0 ; i < e.size(2) ; ++i){
	face & fA  =  dynamic_cast< face & > (*e.get(2, i ));
	const vertex &vA = othervertex(fA, v0, v1);
	vector3d<double> pA, gA;
	bool hasgA = gls0.get(vA, gA);
	if (hasgA){
	getCoord(mi, vA, pA);
	double lA = ls0(vA);
	if (updater.isShockNeeded(x01, g0, gA )){
	bool iss =  updater.compute_shockpoint2(x01, l0, lA+dot(gA, p1-pA), g0, gA, ss, vs);
	if (iss) {found =true; break;}
	}
	if (updater.isShockNeeded(x01, gA, g1 )){
	bool iss =  updater.compute_shockpoint2(x01, lA+dot(gA, p0- pA), l1, gA, g1, ss, vs); 
	if (iss) {found =true; break;}
	}
	}
	}
	if(found){
	if (is_close_to_vertex(ss)) std::cout << "CLOSE" << std::endl;
	auto p = p0*(1-ss) + p1*(ss);
	auto pv = enriched_mesh.createVertex(p, pe, pskeleton_geo);
	skManager.setData(e) = pv;
	lsr.set(*pv, vs);
	continue;
	}
	//shock is needed but found nothing !
	std::cout << "REEE ARRRRRG " << std::endl;
	updater.compute_shockpoint2(x01, l0, l1, g0, g1, ss, vs);
	if (is_close_to_vertex(ss)) std::cout << "CLOSE" << std::endl;
	    
	auto p = p0*(1-ss) + p1*(ss);
	auto pv = enriched_mesh.createVertex(p, pe, pskeleton_geo);
	skManager.setData(e) = pv;
	lsr.set(*pv, vs);
	     	    
	}
	}
	}
	}
  
	for (auto pf : getEntityRange<2>(mi) ){
	AOMD::mFace & f = dynamic_cast<AOMD::mFace &>(*pf);
	auto &e01_ = dynamic_cast< AOMD::mEdge & >(*(f.get(1,0)));
	auto &e12_ = dynamic_cast< AOMD::mEdge & >(*(f.get(1,1)));
	auto &e20_ = dynamic_cast< AOMD::mEdge & >(*(f.get(1,2)));
	auto hpv01_son_ = skManager.getData(e01_);
	auto hpv12_son_ = skManager.getData(e12_);
	auto hpv20_son_ = skManager.getData(e20_);
	AOMD::mVertex * pv0_son_ = nullptr;
	AOMD::mVertex * pv1_son_ = nullptr;
	AOMD::mVertex * pv2_son_ = nullptr;
	AOMD::mVertex * pv0_son  = nullptr;
	AOMD::mVertex * pv1_son  = nullptr;
	AOMD::mVertex * pv2_son  = nullptr;
	AOMD::mVertex * pv01_son = nullptr;
	AOMD::mVertex * pv12_son = nullptr;
	AOMD::mEdge   * pe01 = nullptr;
	AOMD::mEdge   * pe12 = nullptr; 
	AOMD::mEdge   * pe20 = nullptr; 
	bool cut2 = false;

	if (hpv01_son_ || hpv12_son_ || hpv20_son_){
	vertex &v0 =  dynamic_cast< AOMD::mVertex &>(*f.get(0,0));
	vertex &v1 =  dynamic_cast< AOMD::mVertex &>(*f.get(0,1));
	vertex &v2 =  dynamic_cast< AOMD::mVertex &>(*f.get(0,2));
	
	pv0_son_ = enriched_mesh.clone_vertex(v0);
	pv1_son_ = enriched_mesh.clone_vertex(v1);
	pv2_son_ = enriched_mesh.clone_vertex(v2);
	double vls0, vls1, vls2;
       	ls0.get( v0, vls0);
       	ls0.get( v1, vls1);
       	ls0.get( v2, vls2);
	
       	lsr.set( *pv0_son_, vls0);
       	lsr.set( *pv1_son_, vls1);
        lsr.set( *pv2_son_, vls2);
	}
	if ((hpv01_son_ && hpv12_son_) && hpv20_son_){
	auto pv01_son_ = dynamic_cast< AOMD::mVertex * > ( *hpv01_son_);
	auto pv12_son_ = dynamic_cast< AOMD::mVertex * > ( *hpv12_son_);
	auto pv20_son_ = dynamic_cast< AOMD::mVertex * > ( *hpv20_son_);
	auto &v0 = dynamic_cast<const AOMD::mVertex & >(*(f.get(0,0)));
	auto &v1 = dynamic_cast<const AOMD::mVertex & >(*(f.get(0,1)));
	auto &v2 = dynamic_cast<const AOMD::mVertex & >(*(f.get(0,2)));
	double vl0 = ls0(v0);
	double vl1 = ls0(v1);
	double vl2 = ls0(v2);
	vector3d<double> p0,p1,p2,g0,g1,g2;
	getCoord(mi, v0, p0);
	getCoord(mi, v1, p1);
	getCoord(mi, v2, p2);
	bool hasg0 = gls0.get(v0, g0);
	bool hasg1 = gls0.get(v1, g1);
	bool hasg2 = gls0.get(v2, g2);
	assert (hasg0 && hasg1 && hasg2);
	
	auto uv_ls = triple_point2 ( p0,  p1, p2, vl0, vl1, vl2, g0, g1, g2);
	auto uv = uv_ls.first;
	double lsT = uv_ls.second;
	double L0 = 1.-uv(0)-uv(1);
	double L1 = uv(0);
	double L2 = uv(1);
	auto p = p0*L0 + p1*L1 + p2*L2;
	auto pe_v0_v01_son  = enriched_mesh.createEdge_oneLevel( pv0_son_, pv01_son_,  &e01_, e01_.getClassification() );
	auto pe_v01_v1_son  = enriched_mesh.createEdge_oneLevel( pv01_son_, pv1_son_,  &e01_, e01_.getClassification() );
	
	auto pe_v1_v12_son  = enriched_mesh.createEdge_oneLevel( pv1_son_, pv12_son_,  &e12_, e12_.getClassification() );
	auto pe_v12_v2_son  = enriched_mesh.createEdge_oneLevel( pv12_son_, pv2_son_,  &e12_, e12_.getClassification() );
	
	auto pe_v2_v20_son  = enriched_mesh.createEdge_oneLevel( pv2_son_, pv20_son_,  &e20_,  e20_.getClassification() );
	auto pe_v20_v0_son  = enriched_mesh.createEdge_oneLevel( pv20_son_, pv0_son_,  &e20_,  e20_.getClassification() );
	if ((L0>=0) && (L1>=0) && (L2>=0)){
	auto pvT_son = enriched_mesh.createVertex(p, pf,  pskeleton_geo);
	lsr.set(*pvT_son, lsT);
	auto pe_v01_vT_son = enriched_mesh.createEdge_oneLevel(pv01_son_, pvT_son , pf, pskeleton_geo);
	auto pe_v12_vT_son = enriched_mesh.createEdge_oneLevel(pv12_son_, pvT_son , pf, pskeleton_geo);
	auto pe_v20_vT_son = enriched_mesh.createEdge_oneLevel(pv20_son_, pvT_son, pf, pskeleton_geo);
	auto pe_v0_vT_son =  enriched_mesh.createEdge_oneLevel(pv0_son_, pvT_son, pf, pf->getClassification() );
	auto pe_v1_vT_son =  enriched_mesh.createEdge_oneLevel(pv1_son_, pvT_son, pf, pf->getClassification() );
	auto pe_v2_vT_son =  enriched_mesh.createEdge_oneLevel(pv2_son_, pvT_son, pf, pf->getClassification() );
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v0_v01_son, pe_v01_vT_son, pe_v0_vT_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son, pe_v1_vT_son, pe_v01_vT_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v1_v12_son, pe_v12_vT_son, pe_v1_vT_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v12_v2_son, pe_v2_vT_son, pe_v12_vT_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v2_v20_son, pe_v20_vT_son, pe_v2_vT_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v20_v0_son, pe_v0_vT_son, pe_v20_vT_son, &f, f.getClassification());
	}
	else if (L0< 0){
	auto pe_v01_v12_son = enriched_mesh.createEdge_oneLevel(pv01_son_, pv12_son_ , pf, pskeleton_geo);
	auto pe_v20_v12_son = enriched_mesh.createEdge_oneLevel(pv20_son_, pv12_son_ , pf, pskeleton_geo);
	auto pe_v0_v12_son = enriched_mesh.createEdge_oneLevel(pv0_son_, pv12_son_, pf, pf->getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v0_v01_son, pe_v01_v12_son, pe_v0_v12_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son, pe_v1_v12_son, pe_v01_v12_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v12_v2_son, pe_v2_v20_son, pe_v20_v12_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v20_v0_son, pe_v0_v12_son, pe_v20_v12_son, &f, f.getClassification());
	}
	else if (L1< 0){
	auto pe_v12_v20_son = enriched_mesh.createEdge_oneLevel(pv12_son_, pv20_son_ , pf, pskeleton_geo);
	auto pe_v01_v20_son = enriched_mesh.createEdge_oneLevel(pv01_son_, pv20_son_ , pf, pskeleton_geo);
	auto pe_v1_v20_son = enriched_mesh.createEdge_oneLevel(pv1_son_, pv20_son_, pf, pf->getClassification() );
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v1_v12_son, pe_v12_v20_son, pe_v1_v20_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v12_v2_son, pe_v2_v20_son, pe_v12_v20_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v20_v0_son, pe_v0_v01_son, pe_v01_v20_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son, pe_v1_v20_son, pe_v01_v20_son, &f, f.getClassification());
	}
	else if (L2< 0){
	auto pe_v12_v01_son = enriched_mesh.createEdge_oneLevel(pv12_son_, pv01_son_ , pf, pskeleton_geo);
	auto pe_v20_v01_son = enriched_mesh.createEdge_oneLevel(pv20_son_, pv01_son_ , pf, pskeleton_geo);
	auto pe_v2_v01_son = enriched_mesh.createEdge_oneLevel(pv2_son_, pv01_son_, pf, pf->getClassification() );

	enriched_mesh.createFaceWithEdges_oneLevel(pe_v2_v20_son, pe_v20_v01_son, pe_v2_v01_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v20_v0_son, pe_v0_v01_son, pe_v20_v01_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son, pe_v1_v12_son, pe_v12_v01_son, &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v12_v2_son, pe_v2_v01_son, pe_v12_v01_son, &f, f.getClassification());
	}
	continue;
	}
	if ((hpv01_son_ && hpv12_son_) && (!hpv20_son_))   {
	pv0_son = pv0_son_;
	pv1_son = pv1_son_;
	pv2_son = pv2_son_;
	pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv01_son_);
	pv12_son = dynamic_cast< AOMD::mVertex * > ( *hpv12_son_);
	pe01 = &e01_;
	pe12 = &e12_;
	pe20 = &e20_;
	cut2 =true;
	}
	if ((hpv01_son_ && (!hpv12_son_)) && (hpv20_son_)) {
	pv0_son = pv2_son_;
	pv1_son = pv0_son_;
	pv2_son = pv1_son_;
	pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv20_son_);
	pv12_son = dynamic_cast< AOMD::mVertex * > ( *hpv01_son_);
	pe20     = &e12_;
	pe01     = &e20_;
	pe12     = &e01_;
	cut2 =true;
	}
	if (((!hpv01_son_)  && (hpv12_son_)) && (hpv20_son_)) {
	pv0_son = pv1_son_;
	pv1_son = pv2_son_;
	pv2_son = pv0_son_;
	 
	// pv01_son = pv12_son_;
	//pv12_son = pv20_son_;
	pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv12_son_);
	pv12_son = dynamic_cast< AOMD::mVertex * > ( *hpv20_son_);
	pe20     = &e01_;
	pe01     = &e12_;
	pe12     = &e20_;
	cut2 =true;
	}
	if (cut2){
	//std::cout << "BLABLABLA" << std::endl;
	auto pe_v01_v12_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv12_son, pf, pskeleton_geo);
	auto pe_v2_v0_son = const_cast< edge * >(getpedge( *pv2_son, *pv0_son));
	if (!pe_v2_v0_son) pe_v2_v0_son =  enriched_mesh.clone_edge(*pe20);
	auto pe_v0_v01_son = enriched_mesh.createEdge_oneLevel(pv0_son, pv01_son, pe01, pe01->getClassification());
	auto pe_v01_v1_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv1_son, pe01, pe01->getClassification());
	auto pe_v1_v12_son = enriched_mesh.createEdge_oneLevel(pv1_son, pv12_son, pe12, pe12->getClassification());
	auto pe_v12_v2_son = enriched_mesh.createEdge_oneLevel(pv12_son, pv2_son, pe12, pe12->getClassification());
	auto pe_v01_v2_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv2_son, &f,   f.getClassification());
	// could choose the best edge here ...
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v0_v01_son,  pe_v01_v2_son, pe_v2_v0_son,   &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v12_son, pe_v12_v2_son, pe_v01_v2_son,  &f, f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son,  pe_v1_v12_son, pe_v01_v12_son, &f, f.getClassification());
	continue;
	}
	bool cut1 = false;

	if (hpv01_son_){
	//	std::cout << "BLA 01" << std::endl;
	//pv01_son = pv01_son_;
        pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv01_son_);
	pv0_son  = pv0_son_;
	pv1_son  = pv1_son_;
	pv2_son  = pv2_son_;
	pe01     = &e01_;
	pe12     = &e12_;
	pe20     = &e20_;
	cut1 = true;
	}
	if (hpv12_son_){
	//	std::cout << "BLA 01" << std::endl;
	//pv01_son = pv12_son_;
        pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv12_son_);
	pv0_son  = pv1_son_;
	pv1_son  = pv2_son_;
	pv2_son  = pv0_son_;
	pe01     = &e12_;
	pe12     = &e20_;
	pe20     = &e01_;
	cut1 = true;
	}
	if (hpv20_son_){
	//	std::cout << "BLA 20" << std::endl;
	//pv01_son = pv20_son_;
	pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv20_son_);
	pv0_son  = pv2_son_;
	pv1_son  = pv0_son_;
	pv2_son  = pv1_son_;
	pe01     = &e20_;
	pe12     = &e01_;
	pe20     = &e12_;
	cut1 = true;
	}
	if (cut1){
	auto pe_v01_v2_son = enriched_mesh.createEdge_oneLevel( pv01_son, pv2_son, &f,   f.getClassification());
	auto pe_v0_v01_son = enriched_mesh.createEdge_oneLevel(pv0_son, pv01_son, pe01,  pe01->getClassification());
	auto pe_v01_v1_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv1_son, pe01,  pe01->getClassification());
	auto pe_v1_v2_son  = const_cast< edge * >(getpedge( *pv1_son, *pv2_son));
	if (!pe_v1_v2_son) pe_v1_v2_son =  enriched_mesh.clone_edge(*pe12);
	auto pe_v2_v0_son  = const_cast< edge * >(getpedge( *pv2_son, *pv0_son));
	if (!pe_v2_v0_son) pe_v2_v0_son =  enriched_mesh.clone_edge(*pe20);
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v0_v01_son, pe_v01_v2_son, pe_v2_v0_son, &f,  f.getClassification());
	enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son, pe_v1_v2_son, pe_v01_v2_son, &f,  f.getClassification());
	continue;
	}
	}

	std::cout << " SKELETON SIZE " << enriched_mesh.m.size(0) << " " <<enriched_mesh.m.size(1) << std::endl;  
	}
      */
      ////
      ///
      void build_enriched_mesh_v2(const MESHINTERFACE & mi ) {
	std::stringstream tag_name_ske;
	tag_name_ske <<"skeleton_edge_tag"<< this; 
	tag_skeleton_edge  =  getNewTag(mi_enr, tag_name_ske.str()); 
	auto pskeleton_geo =  enriched_mesh.m.getGEntity(tag_skeleton_edge, 1);
	xinterface::aomd::xAttachedDataManagerAOMD <  AOMD::mEntity *> skManager;
	const double eps = 1.e-2;
	std::function < bool (double )> is_close_to_vertex = [&eps]( double s ){ assert(( s >= 0.) && (s <= 1.)); if (s<= eps) return true; if (s >= (1.-eps)) return true; return false; };
	for (auto pe : getEntityRange<1>(mi) ){
	  AOMD::mEdge & e = dynamic_cast<AOMD::mEdge &>(*pe);
	  if (!skManager.getData(e) ) {
	    auto &v0 = dynamic_cast<const AOMD::mVertex & >(*(e.get(0,0)));
	    auto &v1 = dynamic_cast<const AOMD::mVertex & >(*(e.get(0,1)));
	    vector3d<double> p0,p1,g0,g1;
	    getCoord(mi, v0, p0);
	    getCoord(mi, v1, p1);
	    vector3d<double > x01 = p1-p0;
	    bool hasg0 = gls0.get(v0, g0);
	    bool hasg1 = gls0.get(v1, g1);
	    if (hasg0 && hasg1){
	      double l0 = ls0(v0);
	      double l1 = ls0(v1);
	      double ss, vs;
	      bool iss =  updater.compute_shockpoint2(x01, l0, l1, g0, g1, ss, vs);
	      if (iss){
		auto pv = enriched_mesh.createVertex(ss, &e, pskeleton_geo);
		skManager.setData(e) = pv;
		lsr.set(*pv, vs);
		continue;
	      }
	      if ( (!iss) && (updater.isShockNeeded(x01, g0, g1 )) ){
		bool found = false;
		std::cout << " ARRRRRG " << std::endl;
		for (int i =0 ; i < e.size(2) ; ++i){
		  face & fA  =  dynamic_cast< face & > (*e.get(2, i ));
		  const vertex &vA = othervertex(fA, v0, v1);
		  vector3d<double> pA, gA;
		  bool hasgA = gls0.get(vA, gA);
		  if (hasgA){
		    getCoord(mi, vA, pA);
		    double lA = ls0(vA);
		    if (updater.isShockNeeded(x01, g0, gA )){
		      bool iss =  updater.compute_shockpoint2(x01, l0, lA+dot(gA, p1-pA), g0, gA, ss, vs);
		      if (iss) {found =true; break;}
		    }
		    if (updater.isShockNeeded(x01, gA, g1 )){
		      bool iss =  updater.compute_shockpoint2(x01, lA+dot(gA, p0- pA), l1, gA, g1, ss, vs); 
		      if (iss) {found =true; break;}
		    }
		  }
		}
		if(found){
		  auto pv = enriched_mesh.createVertex(ss, &e, pskeleton_geo);
		  skManager.setData(e) = pv;
		  lsr.set(*pv, vs);
		  continue;
		}
		//shock is needed but found nothing !
		std::cout << "REEE ARRRRRG " << std::endl;
		updater.compute_shockpoint2(x01, l0, l1, g0, g1, ss, vs);
		auto pv = enriched_mesh.createVertex(ss, &e, pskeleton_geo);
		skManager.setData(e) = pv;
		lsr.set(*pv, vs);
	     	    
	      }
	    }
	  }
	}
  
	for (auto pf : getEntityRange<2>(mi) ){
	  AOMD::mFace & f = dynamic_cast<AOMD::mFace &>(*pf);
	  auto &e01_ = dynamic_cast< AOMD::mEdge & >(*(f.get(1,0)));
	  auto &e12_ = dynamic_cast< AOMD::mEdge & >(*(f.get(1,1)));
	  auto &e20_ = dynamic_cast< AOMD::mEdge & >(*(f.get(1,2)));
	  auto hpv01_son_ = skManager.getData(e01_);
	  auto hpv12_son_ = skManager.getData(e12_);
	  auto hpv20_son_ = skManager.getData(e20_);
	  AOMD::mVertex * pv0_son_ = nullptr;
	  AOMD::mVertex * pv1_son_ = nullptr;
	  AOMD::mVertex * pv2_son_ = nullptr;
	  AOMD::mVertex * pv0_son  = nullptr;
	  AOMD::mVertex * pv1_son  = nullptr;
	  AOMD::mVertex * pv2_son  = nullptr;
	  AOMD::mVertex * pv01_son = nullptr;
	  AOMD::mVertex * pv12_son = nullptr;
	  AOMD::mEdge   * pe01 = nullptr;
	  AOMD::mEdge   * pe12 = nullptr; 
	  AOMD::mEdge   * pe20 = nullptr; 
	  bool cut2 = false;

	  if (hpv01_son_ || hpv12_son_ || hpv20_son_){
	    vertex &v0 =  dynamic_cast< AOMD::mVertex &>(*f.get(0,0));
	    vertex &v1 =  dynamic_cast< AOMD::mVertex &>(*f.get(0,1));
	    vertex &v2 =  dynamic_cast< AOMD::mVertex &>(*f.get(0,2));
	
	    pv0_son_ = enriched_mesh.clone_vertex(v0);
	    pv1_son_ = enriched_mesh.clone_vertex(v1);
	    pv2_son_ = enriched_mesh.clone_vertex(v2);
	    double vls0, vls1, vls2;
	    ls0.get( v0, vls0);
	    ls0.get( v1, vls1);
	    ls0.get( v2, vls2);
	
	    lsr.set( *pv0_son_, vls0);
	    lsr.set( *pv1_son_, vls1);
	    lsr.set( *pv2_son_, vls2);
	  }
	  if ((hpv01_son_ && hpv12_son_) && hpv20_son_){
	    auto pv01_son_ = dynamic_cast< AOMD::mVertex * > ( *hpv01_son_);
	    auto pv12_son_ = dynamic_cast< AOMD::mVertex * > ( *hpv12_son_);
	    auto pv20_son_ = dynamic_cast< AOMD::mVertex * > ( *hpv20_son_);
	    auto &v0 = dynamic_cast<const AOMD::mVertex & >(*(f.get(0,0)));
	    auto &v1 = dynamic_cast<const AOMD::mVertex & >(*(f.get(0,1)));
	    auto &v2 = dynamic_cast<const AOMD::mVertex & >(*(f.get(0,2)));
	    double vl0 = ls0(v0);
	    double vl1 = ls0(v1);
	    double vl2 = ls0(v2);
	    vector3d<double> p0,p1,p2,g0,g1,g2;
	    getCoord(mi, v0, p0);
	    getCoord(mi, v1, p1);
	    getCoord(mi, v2, p2);

            assert (gls0.get(v0, g0) && gls0.get(v1, g1) && gls0.get(v2, g2));
	
	    auto uv_ls = triple_point2 ( p0,  p1, p2, vl0, vl1, vl2, g0, g1, g2);
	    auto uv = uv_ls.first;
	    double lsT = uv_ls.second;
	    auto pe_v0_v01_son  = enriched_mesh.createEdge_oneLevel( pv0_son_, pv01_son_,  &e01_ );
	    auto pe_v01_v1_son  = enriched_mesh.createEdge_oneLevel( pv01_son_, pv1_son_,  &e01_ );
	
	    auto pe_v1_v12_son  = enriched_mesh.createEdge_oneLevel( pv1_son_, pv12_son_,  &e12_ );
	    auto pe_v12_v2_son  = enriched_mesh.createEdge_oneLevel( pv12_son_, pv2_son_,  &e12_ );
	
	    auto pe_v2_v20_son  = enriched_mesh.createEdge_oneLevel( pv2_son_, pv20_son_,  &e20_ );
	    auto pe_v20_v0_son  = enriched_mesh.createEdge_oneLevel( pv20_son_, pv0_son_,  &e20_ );
	
	    auto pvT_son = enriched_mesh.createVertex(uv(0), uv(1) , &f,  pskeleton_geo);
	    lsr.set(*pvT_son, lsT);
	
	    auto pe_v0_vT_son =  enriched_mesh.createEdge_oneLevel(pv0_son_, pvT_son, &f );
	    auto pe_v1_vT_son =  enriched_mesh.createEdge_oneLevel(pv1_son_, pvT_son, &f );
	    auto pe_v2_vT_son =  enriched_mesh.createEdge_oneLevel(pv2_son_, pvT_son, &f );
	    auto pe_v01_vT_son = enriched_mesh.createEdge_oneLevel(pv01_son_, pvT_son , &f, pskeleton_geo);
	    auto pe_v12_vT_son = enriched_mesh.createEdge_oneLevel(pv12_son_, pvT_son , &f, pskeleton_geo);
	    auto pe_v20_vT_son = enriched_mesh.createEdge_oneLevel(pv20_son_, pvT_son, &f, pskeleton_geo);	
	    enriched_mesh.createFaceWithEdges_oneLevel(pe_v0_v01_son, pe_v01_vT_son, pe_v0_vT_son, &f );
	    enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son, pe_v1_vT_son, pe_v01_vT_son, &f);
	    enriched_mesh.createFaceWithEdges_oneLevel(pe_v1_v12_son, pe_v12_vT_son, pe_v1_vT_son, &f);
	    enriched_mesh.createFaceWithEdges_oneLevel(pe_v12_v2_son, pe_v2_vT_son, pe_v12_vT_son, &f);
	    enriched_mesh.createFaceWithEdges_oneLevel(pe_v2_v20_son, pe_v20_vT_son, pe_v2_vT_son, &f);
	    enriched_mesh.createFaceWithEdges_oneLevel(pe_v20_v0_son, pe_v0_vT_son, pe_v20_vT_son, &f);
	
	    continue;
	  }
	  if ((hpv01_son_ && hpv12_son_) && (!hpv20_son_))   {
	    pv0_son = pv0_son_;
	    pv1_son = pv1_son_;
	    pv2_son = pv2_son_;
	    pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv01_son_);
	    pv12_son = dynamic_cast< AOMD::mVertex * > ( *hpv12_son_);
	    pe01 = &e01_;
	    pe12 = &e12_;
	    pe20 = &e20_;
	    cut2 =true;
	  }
	  if ((hpv01_son_ && (!hpv12_son_)) && (hpv20_son_)) {
	    pv0_son = pv2_son_;
	    pv1_son = pv0_son_;
	    pv2_son = pv1_son_;
	    pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv20_son_);
	    pv12_son = dynamic_cast< AOMD::mVertex * > ( *hpv01_son_);
	    pe20     = &e12_;
	    pe01     = &e20_;
	    pe12     = &e01_;
	    cut2 =true;
	  }
	  if (((!hpv01_son_)  && (hpv12_son_)) && (hpv20_son_)) {
	    pv0_son = pv1_son_;
	    pv1_son = pv2_son_;
	    pv2_son = pv0_son_;
	    pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv12_son_);
	    pv12_son = dynamic_cast< AOMD::mVertex * > ( *hpv20_son_);
	    pe20     = &e01_;
	    pe01     = &e12_;
	    pe12     = &e20_;
	    cut2 =true;
	  }
	  if (cut2){

	    auto point01 = getPoint(*pv01_son);	
	    auto point12 =  getPoint(*pv12_son);
	    auto point0 =  getPoint(*pv0_son);
	    auto point2 = getPoint(*pv2_son);
	
	    vector3d<double> vect01_12(point12(0) - point01(0), point12(1) - point01(1), point12(2) - point01(2)  );
	    vector3d<double> vect2_12(point12(0) - point2(0), point12(1) - point2(1), point12(2) - point2(2)  );
	    vector3d<double> vect0_01(point01(0) - point0(0), point01(1) - point0(1), point01(2) - point0(2)  );
	    double ac0 = externprod(vect01_12, vect2_12 ).norm2();
	    double ac1 = externprod(vect0_01, vect01_12 ).norm2();
	
	    if (ac0 > ac1){
	      auto pe_v2_v0_son = const_cast< edge * >(getpedge( *pv2_son, *pv0_son));
	      if (!pe_v2_v0_son) pe_v2_v0_son =  enriched_mesh.clone_edge(*pe20);
	      auto pe_v0_v01_son = enriched_mesh.createEdge_oneLevel(pv0_son, pv01_son, pe01);
	      auto pe_v01_v1_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv1_son, pe01);
	      auto pe_v1_v12_son = enriched_mesh.createEdge_oneLevel(pv1_son, pv12_son, pe12);
	      auto pe_v12_v2_son = enriched_mesh.createEdge_oneLevel(pv12_son, pv2_son, pe12);
	      auto pe_v01_v2_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv2_son, &f);
	      auto pe_v01_v12_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv12_son, &f, pskeleton_geo);
	      enriched_mesh.createFaceWithEdges_oneLevel(pe_v0_v01_son,  pe_v01_v2_son, pe_v2_v0_son,   &f);
	      enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v12_son, pe_v12_v2_son, pe_v01_v2_son,  &f);
	      enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son,  pe_v1_v12_son, pe_v01_v12_son, &f);
	      continue;
	    }
	    else{
	      auto pe_v2_v0_son = const_cast< edge * >(getpedge( *pv2_son, *pv0_son));
	      if (!pe_v2_v0_son) pe_v2_v0_son =  enriched_mesh.clone_edge(*pe20);
	      auto pe_v0_v01_son = enriched_mesh.createEdge_oneLevel(pv0_son, pv01_son, pe01);
	      auto pe_v01_v1_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv1_son, pe01);
	      auto pe_v1_v12_son = enriched_mesh.createEdge_oneLevel(pv1_son, pv12_son, pe12);
	      auto pe_v12_v2_son = enriched_mesh.createEdge_oneLevel(pv12_son, pv2_son, pe12);
	      auto pe_v0_v12_son = enriched_mesh.createEdge_oneLevel(pv0_son, pv12_son, &f);
	      auto pe_v01_v12_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv12_son, &f, pskeleton_geo);
	      enriched_mesh.createFaceWithEdges_oneLevel(pe_v0_v01_son,  pe_v01_v12_son, pe_v0_v12_son,   &f);
	      enriched_mesh.createFaceWithEdges_oneLevel(pe_v0_v12_son, pe_v12_v2_son, pe_v2_v0_son,  &f);
	      enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son,  pe_v1_v12_son, pe_v01_v12_son, &f);
	      continue;
	    }
	    continue;
	  }
	  bool cut1 = false;

	  if (hpv01_son_){
	    pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv01_son_);
	    pv0_son  = pv0_son_;
	    pv1_son  = pv1_son_;
	    pv2_son  = pv2_son_;
	    pe01     = &e01_;
	    pe12     = &e12_;
	    pe20     = &e20_;
	    cut1 = true;
	  }
	  if (hpv12_son_){
	    pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv12_son_);
	    pv0_son  = pv1_son_;
	    pv1_son  = pv2_son_;
	    pv2_son  = pv0_son_;
	    pe01     = &e12_;
	    pe12     = &e20_;
	    pe20     = &e01_;
	    cut1 = true;
	  }
	  if (hpv20_son_){
	    pv01_son = dynamic_cast< AOMD::mVertex * > ( *hpv20_son_);
	    pv0_son  = pv2_son_;
	    pv1_son  = pv0_son_;
	    pv2_son  = pv1_son_;
	    pe01     = &e20_;
	    pe12     = &e01_;
	    pe20     = &e12_;
	    cut1 = true;
	  }
	  if (cut1){
	    auto pe_v01_v2_son = enriched_mesh.createEdge_oneLevel( pv01_son, pv2_son, &f);
	    auto pe_v0_v01_son = enriched_mesh.createEdge_oneLevel(pv0_son, pv01_son, pe01);
	    auto pe_v01_v1_son = enriched_mesh.createEdge_oneLevel(pv01_son, pv1_son, pe01);
	    auto pe_v1_v2_son  = const_cast< edge * >(getpedge( *pv1_son, *pv2_son));
	    if (!pe_v1_v2_son) pe_v1_v2_son =  enriched_mesh.clone_edge(*pe12);
	    auto pe_v2_v0_son  = const_cast< edge * >(getpedge( *pv2_son, *pv0_son));
	    if (!pe_v2_v0_son) pe_v2_v0_son =  enriched_mesh.clone_edge(*pe20);
	    enriched_mesh.createFaceWithEdges_oneLevel(pe_v0_v01_son, pe_v01_v2_son, pe_v2_v0_son, &f);
	    enriched_mesh.createFaceWithEdges_oneLevel(pe_v01_v1_son, pe_v1_v2_son, pe_v01_v2_son, &f);
	    continue;
	  }
	}

	//std::cout << " SKELETON SIZE " << enriched_mesh.m.size(0) << " " <<enriched_mesh.m.size(1) << std::endl;  
      }


  
      /*
	void fit_to_vertex(const double eps){
    
	std::function < bool (double )> is_close_to_vertex = [&eps]( double s ){ assert(( s >= 0.) && (s <= 1.)); if (s<= eps) return true; if (s >= (1.-eps)) return true; return false; };
	auto itv = enriched_mesh.m.begin(0);
	auto itve = enriched_mesh.m.end(0);
	for (;itv; itv!=itve; ++itv){
	vertex * pv = static_cast< vertex *> (*itve);
	auto pe = enriched_mesh.get_parent_edge(pv);
	if (pe){
	double s =  position_on_edge( getPoint(*v, 0),  getPoint(*v, 1), pv->point());
	if is_close_to_vertex(s) {
	    
	}
	}
	}


	}
      */

    public:
      ~skeleton_enriched_ls(){
	releaseTag(mi_enr, tag_skeleton_edge );
      }

      xfem::xClassRegion getSkeletonRegion() const{
	return xfem::xClassRegion(const_cast< xfem::xMesh *>(&enriched_mesh.m), tag_skeleton_edge, 1);
      }
  
 
      embeded_mesh_fit_to_vertex  enriched_mesh;
      xfem::xRegion enriched_mesh_r;
      MESHINTERFACE mi_enr;
    private:
      unsigned int tag_skeleton_edge;
      ls_t    ls0;
      ls_t    lsr; 
      gls_t   gls0;
      trans_t trans0;
      FMupdaterGeneric2<meshinterface, vector3d<double >, double , int> updater;

    };




    std::pair<vector2d<double >  , double > triple_point2 ( vec3 X0, vec3 X1, vec3 X2, double l0, double l1, double l2, vec3 g0, vec3 g1, vec3 g2){
      auto X01 = X1-X0;
      auto X02 = X2-X0;
      double X01g0 = dot(X01, g0);
      double X02g0 = dot(X02, g0);
      double X01g1 = dot(X01, g1);
      double X02g1 = dot(X02, g1);
      double X01g2 = dot(X01, g2);
      double X02g2 = dot(X02, g2);
      double A00 = X01g0-X01g1;
      double A01 = X02g0-X02g1;
      double A10 = X01g0-X01g2;
      double A11 = X02g0-X02g2;
  
      tensor2d<double> A{A00, A10, A01, A11 };
      vector2d<double> Y{-l0+l1 -X01g1, -l0+l2-X02g2};
      vector2d<double> uv;
      solve(A, Y, uv);
      double u = uv(0);
      double v = uv(1);
      vec3   X = X0 + X01*u+X02*v;
      double l0_uv = l0 + dot(X - X0, g0);
      double l1_uv = l1 + dot(X - X1, g1);
      double l2_uv = l2 + dot(X - X2, g2);
      //  std::cout << l0_uv << " "<< l1_uv << " "<< l2_uv << std::endl;
      // std::cout <<" L0 L1 L2 " <<   1-u-v << " " <<u  << " "<< v<< std::endl;
      double LT = std::min(std::min(l0_uv, l1_uv), l2_uv);
      return std::make_pair(uv , LT);
    }
  } //end namespace skeleton
} //end namespace xfastmarching

#endif
