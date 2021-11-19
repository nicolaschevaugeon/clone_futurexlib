/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xcFrontPart_
#define _xcFrontPart_
// std
#include <string>
// xfem
#include "xLevelSet.h"
#include "xMesh.h"
#include "xRegion.h"

struct xParseData;

/// This class define a part of the  front, simply connected or a loop.
class xcFrontPartBase
{
  public:
   enum frontType
   {
      Point,
      LineOpen,
      LineClose,
      None
   };
   xcFrontPartBase(const xfem::xMesh& front_mesh, const std::string& _frontpartname, const xfem::xMesh& mesh,
                   std::function<double(AOMD::mVertex*)> front_distance, const xParseData& parameters);
   virtual ~xcFrontPartBase();
   const std::string& getFrontName() const { return front_part_name; }
   /// The part of the front mesh on which is defined the frontpart
   const xfem::xRegion& getFrontRegion() const { return front_region; }
   double getFrontLenght() const { return totallenght; }
   const xfem::xMesh& getMesh() const { return mesh; }
   const xfem::xMesh& getFrontMesh() const { return front_mesh; }
   const xParseData& getParameters() const { return parameters; }
   std::function<double(AOMD::mVertex*)>& getDistanceFunction() { return front_distance; }
   std::list<AOMD::mVertex*> local_vertices;
   frontType getFrontType() const { return fronttype; }

  private:
   virtual void ParametrizeFront() = 0;
   virtual void restrict_ls1d_support(xfem::xRegion new_support, double init_val = 0.) = 0;

  protected:
   const xfem::xMesh& front_mesh;
   const std::string front_part_name;
   const xfem::xRegion front_region;
   const xfem::xMesh& mesh;
   std::function<double(AOMD::mVertex*)> front_distance;
   frontType fronttype;
   const xParseData& parameters;
   double totallenght;
};

/// Specialisation for front that are just point (2D case)
class xcFrontPartPoint : public xcFrontPartBase
{
  public:
   friend class xcFrontDomainPoint;
   xcFrontPartPoint(const xfem::xMesh& front_mesh, const std::string& _frontname, const xfem::xMesh& mesh,
                    std::function<double(AOMD::mVertex*)> _front_distance, const xParseData& parameters);

  private:
   void restrict_ls1d_support(xfem::xRegion new_support, double init_val = 0.) override;
   void ParametrizeFront() override;
};

/// specialisation for front part that are open (for example a line segment)
class xcFrontPartLineOpen : public xcFrontPartBase
{
  public:
   friend class xcFrontDomainLineOpen;
   xcFrontPartLineOpen(const xfem::xMesh& front, const std::string& _frontname, const xfem::xMesh& mesh,
                       std::function<double(AOMD::mVertex*)> _front_distance, const xParseData& _parameters,
                       AOMD::mVertex* _start, AOMD::mVertex* _end);
   const xfem::xLevelSet& getLss1d() const { return lss1d; }

  private:
   void restrict_ls1d_support(xfem::xRegion new_support, double init_val = 0.) override;
   void ParametrizeFront() override;
   xfem::xLevelSet lss1d;
   AOMD::mVertex *start, *end;
};

/// Specialisation for front that are closed (periodical, for example a circle)
class xcFrontPartLineClose : public xcFrontPartBase
{
  public:
   friend class xcFrontDomainLineClosed;
   xcFrontPartLineClose(const xfem::xMesh& front, const std::string& _frontname, const xfem::xMesh& mesh,
                        std::function<double(AOMD::mVertex*)> _front_distance, const xParseData& _parameters,
                        AOMD::mVertex* _start);
   const xfem::xLevelSet& getLss1dSin() const { return lss1dSin; }
   const xfem::xLevelSet& getLss1dCos() const { return lss1dCos; }

  private:
   void restrict_ls1d_support(xfem::xRegion new_support, double init_val = 0.) override;
   void ParametrizeFront() override;
   xfem::xLevelSet lss1dSin, lss1dCos;
   AOMD::mVertex* start;
};

#endif
