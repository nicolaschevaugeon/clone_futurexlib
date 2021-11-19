/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef _xAttachedDataManagerAOMD__
#define _xAttachedDataManagerAOMD__
#include <iostream>
#include <sstream>
#include <memory>
#include "mAOMD.h"
#include "mAttachableDataContainer.h"

namespace xfem
{


template < class DATATYPE >
class xAttachedDataManagerAOMD_
{
    private:
        typedef  std::shared_ptr < DATATYPE > data_shr_t;
        class xAttachableSharedPnt : public AOMD::mAttachableData
        {
            public:
               data_shr_t  pnt;
        };
    public:
        inline xAttachedDataManagerAOMD_() : tag( 0)
        {
            std::stringstream tagname;
            tagname << "xAttachedDataManagerAOMD_" << this;
            tag = AOMD::AOMD_Util::Instance()->newMeshDataId(tagname.str().c_str());
        }

        // whitout garbage it doesnt give a real copy
        inline xAttachedDataManagerAOMD_(const xAttachedDataManagerAOMD_ & in) : tag(0)
        {
            std::stringstream tagname;
            tagname << "xAttachedDataManagerAOMD_" << this;
            tag = AOMD::AOMD_Util::Instance()->newMeshDataId(tagname.str().c_str());
        }

        inline const DATATYPE * getData( const AOMD::mAttachableDataContainer &e) const
        {
            const DATATYPE *data = nullptr;
            auto attachable_pnt = static_cast < xAttachableSharedPnt * >( e.getData(tag));
            if (attachable_pnt) data = static_cast < const DATATYPE * >( attachable_pnt->pnt );
            return data;
        }

        inline DATATYPE * getData(AOMD::mAttachableDataContainer &e)
        {
            DATATYPE *data = nullptr;
            auto attachable_pnt = static_cast < xAttachableSharedPnt * >( e.getData(tag));
            if (attachable_pnt) data = attachable_pnt->pnt.get();
            return data;
        }

        inline DATATYPE & setData (AOMD::mAttachableDataContainer &e)
        {
            DATATYPE * data = getData(e);
            if (!data)
            {
                data = new DATATYPE;
                xAttachableSharedPnt *adata = new  xAttachableSharedPnt;
                adata->pnt.reset(data);
                e.attachData(tag, adata);
            }
            return *data;
        }

        inline void   deleteData( AOMD::mAttachableDataContainer &e )
        {
            DATATYPE * data = getData(e);
            //if (data) delete( data );
            e.deleteData(tag);
        }

        /*
        inline void clear()
        {
            for (auto e : garbage)
            {
                deleteData(*e);
            }
            garbage.clear();
        }
        */

        inline ~xAttachedDataManagerAOMD_()
        {
            //clear();
            AOMD::AOMD_Util::Instance()->deleteMeshDataId(tag);
            // tagmanager::instance().releasetag(tag);
        }

    private:
        size_t tag;
};

} // end namespace

#endif
