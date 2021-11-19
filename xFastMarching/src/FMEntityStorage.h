/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _FMENTITYSTORAGE_H
#define _FMENTITYSTORAGE_H


namespace xfastmarching
{



template < class MESHINTERFACE, class ENTITY, class T >
class entitystorage
{
    public:
        entitystorage(const MESHINTERFACE &_mi) : datatag(0), mi(_mi)
        {
            std::stringstream tagname;
            tagname << "entitystorage_" << this;
            datatag = getNewTag(mi, tagname.str());
        };

        entitystorage( const entitystorage < MESHINTERFACE, ENTITY,  T > &in) : datatag(0), mi(in.mi)
        {
            std::stringstream tagname;
            tagname << "entitystorage_" << this;
            datatag = getNewTag(mi, tagname.str());
            for_each(in.vlist.begin(), in.vlist.end(), [&in, this] (const ENTITY *v){
                         T val = xtool::xDataType<T>::zero(); 
                         setData(*v)= *(in.getData(*v));
                     }  );
        }


        const T * getData (const ENTITY &v) const
        {
            void *pp = getAttachedDataPointer(mi,  const_cast < ENTITY & >( v ), datatag);
            assert(pp != nullptr);
            return ( static_cast < const T * > ( pp ));
        }
        T * getData (const ENTITY &v)
        {
            void *pp = getAttachedDataPointer(mi,  const_cast < ENTITY & >( v ), datatag);
            return ( static_cast < T * > ( pp ));
        }
        T& setData (const ENTITY & v)
        {
            void * pp = getAttachedDataPointer(mi, const_cast < ENTITY & >( v ), datatag);
            if (!pp)
            {
                pp = new T;
                attachDataPointer(mi, const_cast < ENTITY & > ( v ), datatag, pp);
                vlist.push_back(&v);
            }
            return  *static_cast < T * > ( pp );
        }
        ~entitystorage()
        {
            std::size_t datatag_ = datatag;
            for_each(vlist.begin(), vlist.end(), [&datatag_, this](const ENTITY *v){
                         void *pp = getAttachedDataPointer( mi, const_cast < ENTITY &  >( *v ), datatag_);
                         delete static_cast < T * >( pp );
                         deleteData( mi,  const_cast < ENTITY &  >( *v ),  datatag_);
                     });
            releaseTag(mi, datatag);
        }

    private:
        std::size_t datatag;
        std::list < const ENTITY * > vlist;
        const MESHINTERFACE &mi;
};

} // enf of namespace xfastmarching
#endif
