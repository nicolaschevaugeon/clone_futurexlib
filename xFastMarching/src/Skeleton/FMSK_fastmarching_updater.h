/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/
#ifndef _FMSK_fastmarching_updater_
#define _FMSK_fastmarching_updater_
#include "FMSK_fastmarching.h"

namespace xfastmarching
{

  namespace skeleton
  {
    // MESHINTERFACE interface represent the type of the interface to mesh
    // T represent the Type of data that can be transported along the gradient
    // GEOMVECT    is the type used to represent vectors

    // "Shock aware FMUpdater."

    template <class MESHINTERFACE, class VECT, class SCAL, class TRANSPORTED >
      class FMupdaterGeneric2 : public FMupdaterGeneric <MESHINTERFACE, VECT, SCAL, TRANSPORTED >{
    public:
      using scal = SCAL;
      using vect = VECT;
      using transported = TRANSPORTED;
      using vertex = typename MESHINTERFACE::vertex;
      using edge   = typename MESHINTERFACE::edge;
      using face   = typename MESHINTERFACE::face;
      using region = typename MESHINTERFACE::region;
      using r_container_t   = std::vector<const region * >;
      using f_container_t   = std::vector<const face * >;
      using entitystorage_t = entitystorage<  MESHINTERFACE, vertex,  scal >;
      using gls_t           = entitystorage<  MESHINTERFACE, vertex,  vect >;
      using vn_t            = entitystorage<  MESHINTERFACE, vertex,  transported > ;
  
    FMupdaterGeneric2(const MESHINTERFACE &_mi,  entitystorage_t &_ls,
		      std::function< scal (const vertex &) > _Ffunc, scal _epsilon_ratio, scal _shock_treshold_angle ):FMupdaterGeneric<MESHINTERFACE, vect, scal, transported >(_mi, _ls, _Ffunc, _epsilon_ratio), shock_treshold_angle(_shock_treshold_angle) {}
    FMupdaterGeneric2(const MESHINTERFACE &_mi,  entitystorage_t &_ls,
		      std::function< scal (const vertex &) > _Ffunc, scal _epsilon_ratio,scal _shock_treshold_angle,  gls_t &_gls  ):FMupdaterGeneric<MESHINTERFACE, vect, scal, transported >(_mi, _ls, _Ffunc, _epsilon_ratio, _gls), shock_treshold_angle(_shock_treshold_angle){}
    FMupdaterGeneric2(const MESHINTERFACE &_mi,  entitystorage_t &_ls,
		      std::function< scal (const vertex &) > _Ffunc, scal _epsilon_ratio, scal _shock_treshold_angle, gls_t &_gls, vn_t &_vn  ):FMupdaterGeneric<MESHINTERFACE, vect, scal, transported >(_mi, _ls, _Ffunc, _epsilon_ratio, _gls, _vn), shock_treshold_angle(_shock_treshold_angle){}

      bool updatefp(const vertex &vx, const edge & e, const scal &epsilon)  const {
	auto & mi = FMupdaterGeneric<MESHINTERFACE, vect, scal, transported>::mi;
	auto pvn =  FMupdaterGeneric<MESHINTERFACE, vect, scal, transported>::pvn;
	auto pgls = FMupdaterGeneric<MESHINTERFACE, vect, scal, transported>::pgls;
	auto &ls = FMupdaterGeneric<MESHINTERFACE, vect, scal, transported>::ls;

	scal valx = ls( vx);
	vect px = this->getcoord(vx);
	// std::cout << "Treating Point " << px<< std::endl;
	scal Fx = this->Ffunc(vx);
	scal Trial;
	vect GradTrial;
	vect Grad;
	transported Trans;
	bool updated = false;
	const std::function<bool ( const vertex *v) > is_known = [this]( const vertex *v) {
	  //auto s = this->status(v);
	  //bool b = (s==FMStatus::K);
	  //std::cout << s << " " << b << std::endl;
	  return this->status(v) ;//==FMStatus::K  ;
	};
	std::function<bool (double, double)> less=[epsilon](double v, double w){ return v<w+epsilon; };
  
	r_container_t  e2r; 
	f_container_t  e2f;
	e2r.reserve(50);
	e2f.reserve(50);

	getRegions(mi, e, std::back_inserter(e2r));
	for (auto r : e2r){
	  const auto va  = getothernodes(mi, vx, *r);
	  if( std::all_of(va.begin(), va.end(), is_known ) ) {
	    const vect T = { ls( *va[0]), ls( *va[1]), ls( *va[2])};
	    if ( less(T(0), valx) || less(T(1), valx) || less(T(2), valx)){
	      const std::array<vect, 3> a ={this->getcoord(*va[0]) -px, this->getcoord(*va[1]) -px, this->getcoord(*va[2])-px  };
	      const tensor3d<scal> dualG = dual( tensor3d<scal> {a[0], a[1], a[2]} );
	      if (computeTrial(T, dualG, Fx, Trial, GradTrial)){
		if(less(Trial,valx)) {
		  updated =true;
		  valx = Trial;
		  Grad = GradTrial;
		  if (pvn){
		    const std::array<transported, 3> trans_other { (*pvn)(*va[0]), (*pvn)(*va[1]), (*pvn)(*va[2])  }; 
		    auto gc =  Grad * dualG;
		    Trans = (gc(0) *trans_other[0] + gc(1)*trans_other[1]+gc(2)*trans_other[2])/(gc(0) + gc(1) + gc(2));
		  }
		}
	      }
	    }
	  }
	}
    
	getFaces(mi,   e, std::back_inserter(e2f));
	for (auto f : e2f){
	  const auto va = getothernodes(mi, vx, *f);
	  //      std::cout << "   " << this->getcoord(*va[0]) << ", "<< this->getcoord(*va[1]) << std::endl;
	  //std::cout << "   " <<  std::all_of(va.begin(), va.end(), is_known )<< std::endl ;
	  if( std::all_of(va.begin(), va.end(), is_known ) ) {
	    const vector2d<scal > T = { ls( *va[0]), ls( *va[1]) };
	    if ( less(T(0),valx) || less(T(1), valx) ){
	      //std::cout << "   LESS " << std::endl; 
	      const std::array<vect, 2> a  {this->getcoord(*va[0])-px, this->getcoord(*va[1])-px  };
	      // some if gls here ...
	      bool shocked = false;
	      vect g0, g1;
	      vector3d<double > x01  (a[1][0] - a[0][0], a[1][1] - a[0][1], 0.);
	      bool hasg0 = pgls->get(*va[0], g0);
	      bool hasg1 = pgls->get(*va[1], g1);
	      double ss, vs;
	      if (hasg0 && hasg1){
		shocked = compute_shockpoint2(x01, T(0), T(1), g0, g1, ss, vs );
	      }
	      if(!shocked){
		const auto  dualG  =  dual (tensor3d2d<scal>( a[0], a[1]));
                if (computeTrial2(T, dualG, Fx, Trial, GradTrial)){
		  if(less(Trial,valx)) {
		    updated =true;
		    valx = Trial;
		    Grad = GradTrial;
		    if (pvn){
		      const std::array<transported, 2> trans_other { (*pvn)(*va[0]) , (*pvn)(*va[1]) }; 
		      auto gc =  Grad * dualG;
		      Trans = (gc(0) *trans_other[0] + gc(1)*trans_other[1])/(gc(0)+gc(1));
		    }
		  }
		}
	      } // end no shock
	      if (shocked){
		//auto res = compute_shockpoint(x01, T(0), T(1), g0, g1);
		//vect ashock = a[0]*(1-res.s) + a[1]*res.s;
		auto ashock = a[0]*(1-ss) + a[1]*ss;
		{ // start side 0
		  const auto dualG0S = dual( tensor3d2d<scal> (a[0], ashock) );
		  //		const vector2d<scal> T0S = {T(0), res.val};
		  const vector2d<scal> T0S = {T(0), vs};							    
                  if (computeTrial2(T0S, dualG0S, Fx, Trial, GradTrial)){
		    if(less(Trial,valx)) {
		      updated =true;
		      valx = Trial;
		      Grad = GradTrial;
		      if (pvn){
			const std::array<transported, 2> trans_other { (*pvn)(*va[0]) , (*pvn)(*va[1]) };
			transported transS = trans_other[0]*(1-ss) + trans_other[1]*ss;
			auto gc =  Grad * dualG0S;
			Trans = (gc(0) *trans_other[0] + gc(1)*transS)/(gc(0)+gc(1));
		      }
		    } 
		  }
		} //  end side 0
		{ // start side 1
		  const auto dualGS1 = dual( tensor3d2d<scal> ( ashock, a[1]) );
		  //		const vector2d<scal> T0S = {T(0), res.val};
		  const vector2d<scal> TS1 = { vs, T(1) };							    
                  if (computeTrial2(TS1, dualGS1, Fx, Trial, GradTrial)){
		    if(less(Trial,valx)) {
		      updated =true;
		      valx = Trial;
		      Grad = GradTrial;
		      if (pvn){
			const std::array<transported, 2> trans_other { (*pvn)(*va[0]) , (*pvn)(*va[1]) };
			transported transS = trans_other[0]*(1-ss) + trans_other[1]*ss;
			auto gc =  Grad * dualGS1;
			Trans = (gc(0) *transS + gc(1)*trans_other[1])/(gc(0)+gc(1));
		      }
		    } 
		  }
		} //  end side 1 
	      } // end shocked
	    } // end one node is less
	  } // end all known
	} // end face loop
	const vertex * vb = getothernodes(mi, vx, e); 
	if( is_known(vb)){
	  const scal T = ls( *vb);
	  if ( less(T, valx) ){
	    const vect b = this->getcoord( *vb) -px;
	    if (nrm2(b) <= 1.e-8) {
	      Trial =T;
	      if(less(Trial,valx)){
		updated =true;
		Grad = (*pgls)(*vb);
		if(pvn)
		  Trans = (*pvn)(*vb);
	    
	      }
	    }
            else if (computeTrial(T, dual(b), Fx, Trial, GradTrial)){
	      if(less(Trial,valx)) {
		updated =true;
		valx = Trial;
		Grad = GradTrial;
		if (pvn){
		  Trans = (*pvn)(*vb);
		}
	      }
	    }
	  }
	}
    
	if (updated){
	  ls.set(vx, valx);
	  if (pgls) pgls->set(vx, Grad);
	  if (pvn)  pvn->set(vx, Trans);
	  return true;
	}
	return false;
      }

      const double shock_treshold_angle = M_PI/30.;

      bool isShockNeeded(const vector3d<double > &x01, const vector3d<double > &g0, const vector3d<double > &g1 ) const {
	//const double shock_treshold_angle = M_PI/12.;
	/*const double x01g0 = dot(x01,g0);
	  const double x01g1 = dot(x01,g1);
	  if (x01g0 <= 0. && x01g1 >= 0.){ // expansion
	  return false; 
	  }
	  if ( dot(g0,g1)>= cos(shock_treshold_angle) ){
	  return false;
	  }
	  //if (x01g0*x01g1 > 0) 
	  //  return false;
	  return true;
	*/
	//    return (dot( g0 -g1, x01) >=  sqrt(dot(x01,x01)) * cos (shock_treshold_angle ));
	if  (dot( g0 -g1, x01) <= 0.) //expansion
	  return false;
	if ( dot(g0,g1)>= cos(shock_treshold_angle) ){
	  return false;
	}
	return true;
      };

      struct spointonedge{
	double s;
	double val;
      };
  
      spointonedge compute_shockpoint (const vector3d<double > &x01, double l0, double l1, const vector3d<double > &g0, const vector3d<double > &g1 ) const {
	// const double shock_treshold_angle = M_PI/12.;
	spointonedge ret;
	double gx0 = dot(x01, g0);
	double gx1 = dot(x01, g1);
	assert (fabs (gx0-gx1) > 1.e-12);
	ret.s = (l1 - l0 -gx1)/(gx0-gx1);
	if (!( ( ret.s >= 0.) && (ret.s <= 1.))){
	  std::cout << ret.s<< std::endl;
	  std::cout << l0 << " "<< gx0 << std::endl;
	  std::cout << l1 << " "<< gx1 << std::endl;
	  std::cout << dot(g0, g1) << " " <<  cos(shock_treshold_angle) << std::endl;
	  std::cout << gx0-gx1 << std::endl;
	}
	assert(( ret.s >= 0.) && (ret.s <= 1.));
	ret.val = std::min(l0+ret.s*gx0, l1 + gx1*(ret.s-1.));
	return ret;
      }

      bool compute_shockpoint2 (const vector3d<double > &x01, double l0, double l1, const vector3d<double > &g0, const vector3d<double > &g1, double & s, double & val ) const {
	//l0(s) = l0+s*gx0;
	//l1(s) = l1 + (s-1)*gx1;
	//const double shock_treshold_angle = M_PI/12.;
	if  (dot( g0 -g1, x01) <= 0.) //expansion
	  return false;
	if ( dot(g0,g1)>= cos(shock_treshold_angle) ){
	  return false;
	}
	double gx0 = dot(x01, g0);
	double gx1 = dot(x01, g1);
	double l0_0 = l0;
	double l0_1 = l0+gx0;
	double l1_0 = l1 - gx1;
	double l1_1 = l1;
	if ((l1_0 > l0_0) && ( l1_1 < l0_1)){
	  s = (l1 - l0 -gx1)/(gx0-gx1);
	  val = std::min(l0+s*gx0, l1 + gx1*(s-1.));
	  return true;
	}
	//val = (l0_0 + l1_1)*0.5;
	//s = 0.5;
	//return true;

	if (fabs(l0_0-l1_0) <= fabs(l0_1-l1_1))  {
	  val = l0;
	  s = 0.;
	  return false;
	}
	val = l1;
	s= 1.;
  
	return false;
      }
  
  
    };
  } //end namespace FMSK
} //end namespace  
#endif
