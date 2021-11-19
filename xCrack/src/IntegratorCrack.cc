/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.

*/
  
#include <iostream>
#include "IntegratorCrack.h"
#include "xGeomElem.h"
#include "xCommandOnGeomElem.h"
#include "xMesh.h"


using namespace xfem;
using namespace std;
using namespace AOMD;

namespace xcrack
{

void IntegratorCrack_c::accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const
{ 
  const bool debug = false;
  int deg;
  xfem::xPartition partition;
  if (debug) std::cout <<"On ouvre l element ";
  if (debug) e_integ->print();
  xMesh::getPartition(e_integ, partition, filter);
  if (debug) std::cout << "partition size is " << partition.size() << std::endl;
  deg=degree_gen;
  
// dans le cas ou le filtre ne contient que les supports a enrichir
// "supports_touched_by_front"

/*  bool has_support_touched_by_crack=false;
  int i;
    for(i=0;i<e_integ->size(0);i++)
    {
      mEntity* v=e_integ->get(0,i);
      if (filter_front(v)) {has_support_touched_by_crack=true;break;}
    }


  if (has_support_touched_by_crack) deg= degree;
*/


// filtre ad-hoc contenant les entites devant etre integrees proprement
// "supports_integrate_front"

  if (!filter_front(e_integ)) deg= degree_gen; else deg=degree;
  

  for (xPartition::iterator it = partition.begin(); it != partition.end(); ++it) 
    {
      if (debug) std::cout << "partition found inside sub\n";
      mEntity* es_integ = *it;
      xfem::xGeomElem geo_integ_es(es_integ);
      geo_integ_es.SetIntegrationPointNumberForDegree(deg);
      command.execute(&geo_integ_es);
    }
  if (debug) std::cout<<"--------\n";
}

}
