/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/

#include <sstream>
#include <cstdlib>
#include "xExportManager.h"
#include "xParseData.h"

using namespace std;

namespace xexport
{

  xExportManager::xExportManager(const int step_string_length_) : step_string_length(step_string_length_),file_name("") {}

  void xExportManager::init(const xParseData& export_parser_, const string& xParseData_token_EXPORT_MANAGER)
  {
    list<string> string_list = export_parser_.getListString(xParseData_token_EXPORT_MANAGER);
    if (string_list.size()%2!=0)
      {
	cout << "wrong EXPORT_MANAGER format. Must be { name frequency name frequency }" << endl;
	abort();
      }
    list<string>::const_iterator it = string_list.begin();
    int threshold,freq;
    size_t pos;
    while (it!=string_list.end())
      {
	string val = *it;
	pos=(*(++it)).find_first_of('+');
	if (pos!=string::npos)
	  {
	    freq = atoi((*it).substr(0,pos).c_str());
	    threshold = atoi((*it).substr(pos+1,string::npos).c_str());
	  }
	else
	  {
	    freq = atoi((*it).c_str());
	    threshold=0;
	  }
	++it;
	registerExport(val, freq, threshold);
      }
  }

  void xExportManager::registerExport(const string& export_name_, const int frequency_, const int threshold_)
  {
    pair<map<string, pair<int,int> >::iterator,bool> query=
      info_container.insert(make_pair(export_name_, make_pair(threshold_, frequency_)));
    if (!query.second)
      {
        
        int &t=(*(query.first)).second.first;
        int &f=(*(query.first)).second.second;
        cout << "ExportManager Warning : old value ("<<f<<","<<t<<")for "<<export_name_<<" are replaced by "<<frequency_<<" "<<threshold_<< endl;
        t=threshold_;
        f=frequency_;
      }
  }

  void xExportManager::setNumerationSize(const int step_string_length_)
  {
    step_string_length=step_string_length_;
  }

}
