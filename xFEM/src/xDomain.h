#ifndef XDOMAIN_H
#define XDOMAIN_H

#include "xAttachedDataManagerAOMD.h"
namespace xfem
{
class xDomain
{
  public:
   static int &set(const AOMD::mEntity &e) { return zone().setData(e); }
   static int *get(const AOMD::mEntity &e) { return zone().getData(e); }
   static void del(const AOMD::mEntity &e) { zone().deleteData(e); }
   friend class xMesh;

  private:
   static xinterface::aomd::xAttachedDataManagerAOMD<int> &zone()
   {
      static xinterface::aomd::xAttachedDataManagerAOMD<int> z;
      return z;
   }
};
}  // namespace xfem
#endif  // XDOMAIN_H
