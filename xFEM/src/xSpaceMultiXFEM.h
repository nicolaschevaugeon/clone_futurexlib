/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/


#ifndef _XSPACEMULTIXFEM_
#define _XSPACEMULTIXFEM_
#include <map>
#include <vector>
#include <string>

#include "xApproxFctPtr.h"
#include "xValKey.h"
#include "xSpace.h"
#include "xApproxFunction.h"

namespace xfem{
  class xGeneratorApproxFunctionEnrichedXFEMWith
  {
  public:
    typedef std::vector<xValKey> femKeys;
    typedef femFcts_t femFcts;
    
    virtual void generateKeysAndFcts(AOMD::mEntity *e, xValKey &key,ValKeyModifierPtr key_modifier, approxFctPtr_t f, femKeys* keys, femFcts* fcts) = 0 ;
    virtual void generateKeys(xValKey &key,ValKeyModifierPtr key_modifier,femKeys* keys) = 0 ;
    virtual ~xGeneratorApproxFunctionEnrichedXFEMWith()= default;;
  };
  
  class xSpaceMultiXFEM : public xSpace {
  public:
    //  alternate first version of the verion below
    //    template <class T, class E,class M >
    //      xSpaceMultiXFEM(const T& base, const E& enr, const M& mk) :
    //    base_space(new T(base)), generator_approx((xGeneratorApproxFunctionEnrichedXFEMWith* )(new E(enr))), key_modifier(new M(mk) ) {} 
    template <class T, class M >
      xSpaceMultiXFEM(const T& base,   xGeneratorApproxFunctionEnrichedXFEMWith * enr, const M& mk) :
    base_space(new T(base)), generator_approx(enr), key_modifier(new M(mk) ) {} 

    /// Constructor using a spacePtr as argument for the base space
    template < class M >
      xSpaceMultiXFEM(const spacePtr& base,   xGeneratorApproxFunctionEnrichedXFEMWith * enr, const M& mk) :
    base_space(base), generator_approx(enr), key_modifier(new M(mk) ) {} 

    ~xSpaceMultiXFEM() override{
      // memory leak ??? to check
      // delete generator_approx;
      //delete key_modifier;
    }
    
    void getKeys(AOMD::mEntity* e,  femKeys* keys) override; 
    void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;
    
  protected:
    spacePtr            base_space;
    xGeneratorApproxFunctionEnrichedXFEMWith * generator_approx;
    ValKeyModifierPtr    key_modifier;
  };
}
#endif
