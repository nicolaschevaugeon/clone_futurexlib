#ifndef XDATATYPE__H
#define XDATATYPE__H
#include <string>


/*!
Generic data type abstraction
Useful for the creation of template functions dealing with numerics
Type specializations are done in _imp file
xTool MUST NOT DEPEND on any external library
PLEASE, do NOT add specializations other than c++ standard types
Non-standard specializations should be added in your own code
*/

namespace xtool {

template < typename T >
class xDataType{

public:
    static T zero() { throw; }
    static T one() { throw; }
    static T epsilonNumeric() { throw; }
    static std::string  stype() { throw; }

};
}

#include "xDataType_imp.h"

#endif

