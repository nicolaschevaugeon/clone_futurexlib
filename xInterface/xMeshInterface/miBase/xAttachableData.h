#ifndef XATTACHABLEDATA_H
#define XATTACHABLEDATA_H
#include <memory>
//#include "xTypeTraitsStdUniquePtr.h"

namespace  xinterface {

namespace  xmeshinterface {





template<typename T>
class xAttachableData
{
public:
    xAttachableData(const T& _data) : data(_data) {}
    T& getData() { return data;}
    void setData(T _data) { data=_data;}
private:
    T data;
};

}
}


#endif // XATTACHABLEDATA_H
