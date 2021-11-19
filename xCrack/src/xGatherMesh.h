/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _BroadCastMesh_
#define _BroadCastMesh_

#include "mExchangeData.h"
#include "mFace.h"
#include "xMesh.h"

class xEntityCopyCallBack
{
  public:
   virtual void operator()(AOMD::mEntity &source, AOMD::mEntity &destination) const = 0;
   virtual ~xEntityCopyCallBack() = default;
};

template <template <class> class DATAMANAGER>
class xEntityCopyCallBackAttachEntity : public xEntityCopyCallBack
{
  public:
   xEntityCopyCallBackAttachEntity(const DATAMANAGER<AOMD::mEntity *> &_sourcedata, DATAMANAGER<AOMD::mEntity *> &_targetdata)
       : sourcedata(_sourcedata), targetdata(_targetdata)
   {
   }

   void operator()(AOMD::mEntity &source, AOMD::mEntity &destination) const override
   {
      AOMD::mEntity *e = sourcedata.at(source);
      targetdata.setData(destination) = e;
   }

   ~xEntityCopyCallBackAttachEntity() override = default;

  private:
   const DATAMANAGER<AOMD::mEntity *> &sourcedata;
   DATAMANAGER<AOMD::mEntity *> &targetdata;
};

/// Copy the mesh partition in the current partition into a global mesh. This Function is collective on al process. The global
/// mesh exist on all the process at the end. xEntityCopyCallBack is called for each mesh entity that is copied. It permit to copy
/// additionnalal data on a mesh entity.
void BroadCastMesh(xfem::xMesh *localMesh, xfem::xMesh *globalMesh, xEntityCopyCallBack *);

//////////////////////////

class xExchangeMaxOnBoundary : public AOMD::AOMD_DataExchanger
{
  public:
   xExchangeMaxOnBoundary(int _data_tag);
   int tag() const override;
   void *AP_alloc_and_fill_buffer(AOMD::mEntity *e, AOMD::AOMD_SharedInfo &si, int tag) override;
   void receiveData(int pid, void *buf) override;

  private:
   int data_tag;
};

#endif
