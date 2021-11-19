/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef _EXPORTMANAGER_H
#define _EXPORTMANAGER_H

#include <map>
#include <sstream>
#include <string>
#include <utility>

struct xParseData;

/// Class to manage export in the code.
/// (1) it manages the frequency
/// (2) it manages the file name (adding the step id at the end of the file)
/// To use it, put an export algorithm into an if {} statement
/// that checks the condition toExport()
/// EXAMPLE:
/// if (export_manager.toExport("displacement", 3, ".txt") {
///    export_displacement_algorithm(export_manager.getFileName());
/// }
/// by default (step_string_length), file_name would be "displacement_00003.txt"

namespace xexport
{
class xExportManager
{
  public:
   /// constructor
   xExportManager(const int step_string_length_ = 5);

   /// to read what to export from an xParseData EXPORT_MANAGER label
   /// EXPORT_MANAGER list format must be { name directive name directive .... } where directive is :
   ///       one single integer discribing output frequency
   ///       an integer discribing output frequency followed by '+' character followed by a integer discribing threshold (starting
   ///       potential output)
   void init(const xParseData& export_parser_, const std::string& xParseData_token_EXPORT_MANAGER = "EXPORT_MANAGER");

   /// to manage what to export manually
   void registerExport(const std::string& export_name_, const int frequency_, const int threshold_ = 0);

   /// to know if a specific export is to be exported or not, build the file_name adding the number of the step
   inline bool toExport(const std::string& export_name_, const int step_, const std::string& extension_name_ = "");

   /// to get the current file_name
   inline std::string getFileName() const { return file_name; }

   /// to change default step_string_length
   void setNumerationSize(const int step_string_length_);

  private:
   inline std::string outputId(const int time_id);

   int step_string_length;
   std::string file_name;

   std::map<std::string, std::pair<int, int>> info_container;
};

inline std::string xExportManager::outputId(const int time_id)
{
   std::ostringstream out;
   out.width(step_string_length);
   out.fill('0');
   out << time_id;
   return "_" + out.str();
}

inline bool xExportManager::toExport(const std::string& export_name_, const int step_, const std::string& extension_name_)
{
   std::map<std::string, std::pair<int, int>>::const_iterator it = info_container.find(export_name_);
   if (it != info_container.end())
   {
      const std::pair<int, int>& info = it->second;
      bool test1 = step_ >= info.first;
      bool test2 = step_ % info.second == 0;
      if (test1 && test2)
      {
         file_name = export_name_ + outputId(step_) + extension_name_;
         return true;
      }
   }
   file_name = "";
   return false;
}

}  // namespace xexport

#endif
