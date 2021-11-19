/*
  octree is a subproject of  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include "oKeyManager.h"

namespace xoctree
{

  class oSchemeVector
  {
  public:
  oSchemeVector() : hs(6, 0.) {}
    virtual ~oSchemeVector() {}
    virtual void compute(const oKey* k) = 0;
    const std::vector<double>& getVector() {return vec;}
  protected:
    std::vector<double> vec;
    std::vector<double> hs;
  };

  class oSchemeScalar
  {
  public:
    virtual ~oSchemeScalar() {}
    virtual void compute(const oKey* k) = 0;
    double getScalar() { return sca; }
  protected:
    double sca;
    std::vector<double> hs;
  };

  class oSchemeVectorLaplace : public oSchemeVector
  {
  public:
    void compute(const oKey* key)
    {
      distance(key, key->getStencil(), hs); 
      for (int i = 0; i < key->sizeStencil(); i+=2)
	{
	  double hb = 0.5 * (hs[i]+hs[i+1]);
	  double c  = 1./(hb*hb); 
	  double am = hs[i]/hb;
	  double ap = hs[i+1]/hb;
	  vec[i] = c/am;
	  vec[i+1] = c/ap;
	  vec[hs.size()] += -2.*c/(am*ap);
	}
    }
  };


  class oSchemeScalarUnit : public oSchemeScalar
  {
  public:
    oSchemeScalarUnit {sca = 1.0;}
    void compute(const oKey* key) {}
  };




  //   oMappingCartesian mapping(Dim, box_inf, box_sup);
  //   oOctree octree (mapping, LevelMax, per);

  //   oLevelSetAnalytical lsoct(mapping, LevelMax, ls_circle); 
  //   oRefinementCriteriaOnSign ref_criteria(mapping, lsoct);

  //   octree.optimize(ref_criteria);
  //   oKeyManager key_manager(octree);

  //   oField ls(key_manager, lsoct);

  //   oField ls2(key_manager, lsoct2);


  //   oSchemeVectorLaplace scheme(ls);

  //   key_manager.register("velo_normal");

  //   key_manager.createValues();


 
  //   oField ls(key_manager, lsoct);


  // 	  oField vn(
  //   key_manager.register(

  

  //   ExportGMSHAscii (active_nodes, octree, "active_nodes") ;
  //   ExportGMSHMesh(active_nodes, octree, "active_nodes_mesh");



  // oOctree octree(....);
  // octree.optimi


  // oField 
  // {

  //   double getVal(oKey* k)
  //     {
  //       oKey::setProperVal();
  //       return k->getVal();
  //     }
  // }

  class oSchemeExplicitPropag
  {
  public:
  oSchemeExplicitPropag(const oField& ls_, const oField& velo_) : ls(ls_), velo(velo_) {}
    void compute(const oKey* key)
    {
      const std::vector<oKey*>& stencil = key->getStencil();
      distance(key, stencil, hs);

      double lsc = ls.getVal(key);
      double velo = key->getVelo();
      double grad = 0.;    
      for (int i = 0; i < stencil.size(); i+=2)
	{
	  double lsm = ls.getVal(stencil[i]);
	  double lsp = ls.getVal(stencil[i]);
	  double gm  = (lsc-lsm)/hs[i];
	  double gp  = (lsp-lsc)/hs[i+1];
	  if (velo >= 0.) grad += max(gm,0.)^2 + min(gp,0.)^2;
	  else            grad += max(gp,0.)^2 + min(gm,0.)^2; 
	}
      sca = velo*sqrt(grad);
    }
  private:
    const oField& ls;
    const oField& velo;

  };



  int DeclareDofs(ITER it, ITER end)
  {
    int NumDof = 0;
    for(; it != end; ++it)
      {
	const oKey* key = *it;
	if (key->is_regular() && !key->is_fixed())
	  {
	    key->setDof(NumDof++)
	      }
      }
    return NumDof;
  }



  void Assemble(oSchemeVector& scheme, oAssembler& assembler, ITER it, ITER end)
  {
    for(; it != end; ++it)
      {
	const oKey* key = *it;
	if (key->isDof())
	  {
	    scheme.compute(key);
	    assembler.assemble(key, 
			       key->beginStencil(), 
			       key->endStencil(), 
			       scheme.getVector());
	  }
      }
  }

  template< class M = xCSRMatrix, class V = xCSRVector>
    class oAssembler {
  public:
  typedef M   matrix_type;
  typedef V   vector_type;

  oAssembler(M & A, V & b) :  coeff(1.0), mat(&A), vec(&b) {}
  oAssembler(V & b)        :  coeff(1.0), vec(&b)       {}

  void setCoeff(const double& c) { coeff = c; }
  void assemble(oKey* keyi, key_iterator first, key_iterator last, std::vector<double>& vec)
  {
    for (int i = 0; first != last; ++first, ++i) {
      assembleMatrix(keyi, *firsti, coeff * vec[i]);
    }
  }
  virtual void assembleMatrix(const oKey* keyi, const oKey* keyj, const double& val) 
  {
    if (keyj->isDof() && mat)
      {
	mat->AddMatrix(keyi->getNumDof(), keyj->NumDof(), val);
      }
    else if (!keyj->isRegular() && mat)
      {
	oKey::const_iterator it = keyj->beginTies();
	oKey::const_iterator it = keyj->endTies();
	double coeff = 1./(double) keyj->sizeTies();
	for (it; it !=ite; ++it)
	  {
	    assembleMatrix(keyi, *it, coeff*val);
	  }
      }
    else if (vec)// means keyj is fixed
      {
	vec->AddVector(keyi->getNumDof(), -val*keyj->getVal());
      }
    
  }
  

  private:
  double coeff;
  M * mat;
  V * vec;
  };



  template< class M = xCSRMatrix, class V = xCSRVector>
    class oAssemblerLhsAndRhs : public oAssembler {
  public:
  typedef M   matrix_type;
  typedef V   vector_type;

  oAssemblerLhsAndRhs(M & A, V & b) : oAssembler(A, b) {}
  oAssemblerLhsAndRhs(V & b)        : oAssembler(b) {}      

  virtual void assembleMatrix(const oKey* keyi, const oKey* keyj, const double& val) 
  {
    if (keyj->isDof() && mat && vec)
      {
	mat->AddMatrix(keyi->getNumDof(), keyj->NumDof(), val);
	vec->AddVector(keyi->getNumDof(), val*keyj->getVal());
      }
    else if (!keyj->isRegular() && mat)
      {
	oKey::const_iterator it = keyj->beginTies();
	oKey::const_iterator it = keyj->endTies();
	double coeff = 1./(double) keyj->sizeTies();
	for (it; it !=ite; ++it)
	  {
	    assembleMatrix(keyi, *it, coeff*val);
	  }
      }
    else 
      {
      }
    
  }
  

  private:

  };


} //end namespace xoctree

