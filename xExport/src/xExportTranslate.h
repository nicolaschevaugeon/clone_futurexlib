/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef _translate_ensight_H
#define _translate_ensight_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "xExportEnsight.h"

namespace xexport
{
// export to ensight binary from ascii or binary files.
class xExportTranslate : public xExportEnsightBinary
{
  public:
   xExportTranslate(std::string, std::string, std::string);

  private:
   std::ifstream fpos_;
   std::string project_name_;
   std::string cSrc_;
   bool isBinary_;
   xExportEnsightBinary project_;

   short giveTypeElementGmsh(char a);
   short giveTypeVariable(const char nan);
   bool isFileBinary(const std::string&);
   void readGmshBinary();
   void readGmshAscii();
   void readFileGmsh();
   void writeFileEnsight();

   // second partie : traitement
   std::string giveNomProject(std::string li);
   bool isFinAc(char a);
   bool isSeparateur(char a);
   void traiteVariableScalaire(const std::string& ligne);
   void traiteVariableVector(const std::string& ligne);
   void traiteVariableTensor(const std::string& ligne);
   void traiteVariable(std::string ligne);

   void readPosBinary(int nbr, int cpt);
   void traiteElemBinaire(double* ptr, int cpt);
};
}  // namespace xexport
#endif
