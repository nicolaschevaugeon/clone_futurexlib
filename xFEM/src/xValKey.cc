/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include  "mEntity.h"
#include  "xValKey.h"

// uncoment to have a extended version of xValKey output
//#define XFEM_XVALKEY_OPERATOR_FULL 0

namespace xfem
{

xKeyInfo::string_manager_t  xKeyInfo::PhysString;
xKeyInfo::string_manager_t  xKeyInfo::GeomString;


std::ostream & operator << (std::ostream & o, const xValKey& d) {
  const bool debug=false;
  o << xKeyInfo::getPhysName(d.Phys) << " ";
  o << xKeyInfo::getGeomName(d.Geom) << " ";
  o << "Enti: " << d.Enti->getId() << " ";
#ifdef XFEM_XVALKEY_OPERATOR_FULL
  switch ( d.Enti->getType() )
  {
      case AOMD::mEntity::EDGE :
         {
             o << "Enti type: EDGE on nodes : " << d.Enti->get(0,0)->getId() <<" "<<d.Enti->get(0,1)->getId()<< " ";
             break;
         }
      case AOMD::mEntity::TRI :
        {
             o << "Enti type: TRI on nodes : " << d.Enti->get(0,0)->getId() <<" "<<d.Enti->get(0,1)->getId()<<" "<<d.Enti->get(0,2)->getId()<<" ";
             break;
        }
      case AOMD::mEntity::TET :
        {
             o << "Enti type: TET on nodes : " << d.Enti->get(0,0)->getId() <<" "<<d.Enti->get(0,1)->getId()<<" "<<d.Enti->get(0,2)->getId()<<" "<<d.Enti->get(0,3)->getId()<<" ";
             break;
        }
      default:
             break;
  }
#endif

  if (debug) std::cout << "Enti: " << d.Enti->getId() << " ";
  if (debug) d.Enti->print();

  if (d.Refe != -1) o << "Refe: " << d.Refe;
  return o;
}


xValKeyExtend::xValKeyExtend(const std::string& e) : extension(e) {}
void xValKeyExtend::operator()(xValKey& key) 
{ 
  std::string geom_new = xKeyInfo::getGeomName(key.getGeom()) + extension;
  key.setGeom(xKeyInfo::getGeomId(geom_new));
}

} // end of namespace

