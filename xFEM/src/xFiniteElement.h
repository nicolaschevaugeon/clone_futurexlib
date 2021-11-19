/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _FEM_H
#define _FEM_H

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "mEntity.h"
#include "xApproxFctPtr.h"
#include "xSpacePtr.h"
#include "xValKey.h"

namespace xfem
{
class xGeomElem;
class xApproxFunction;
class xSpace;
class xSpaceBase;

class xFiniteElementKeysOnly
{
  public:
   typedef spaceBasePtr_t spacePtr;
   typedef std::vector<xValKey> femKeys;
   typedef femKeys::iterator iterKey;

   xFiniteElementKeysOnly();
   ~xFiniteElementKeysOnly();
   femKeys* getKeys(const std::string& ityp = "default");
   int sizeKey(const std::string& ityp = "default") const;
   iterKey beginKey(const std::string& ityp = "default");
   iterKey endKey(const std::string& ityp = "default");
   template <class iterFctSpace>
   void setKeys(AOMD::mEntity* e, iterFctSpace it, iterFctSpace last, const std::string& ityp = "default")
   {
      if (e_current != e)
      {
         clear();
         e_current = e;
      }
      femKeys* keys = KeysFor(ityp);
      keys->clear();
      for (; it != last; ++it)
      {
         addKeys(*it, ityp);
      }
   }
   AOMD::mEntity* getEntity() { return e_current; }

  private:
   femKeys* KeysFor(const std::string& ityp = "default");
   void addKeys(spacePtr f, const std::string& appro = "default");

   void clear();
   AOMD::mEntity* e_current;
   typedef std::map<std::string, femKeys*> allKeys_t;
   allKeys_t allKeys;
};

//! The Finite Element
class xFiniteElement
{
  public:
   typedef std::vector<xValKey> femKeys;
   typedef femKeys::iterator iterKey;
   //
   typedef approxFctPtr_t shapeFctPtr;
   typedef femFcts_t femFcts;
   typedef femFcts::iterator iterFct;
   //
   typedef spacePtr_t spacePtr;

   xFiniteElement();
   ~xFiniteElement();

   femFcts* getFcts(const std::string& ityp = "default");
   femKeys* getKeys(const std::string& ityp = "default");

   int sizeKey(const std::string& ityp = "default") const;
   iterKey beginKey(const std::string& ityp = "default");
   iterKey endKey(const std::string& ityp = "default");
   iterFct beginFct(const std::string& ityp = "default");
   iterFct endFct(const std::string& ityp = "default");

   template <class Predicate>
   void filterKey(Predicate pred, const std::string& ityp = "default")
   {
      iterKey new_end = std::remove_if(beginKey(ityp), endKey(ityp), std::not1(pred));
      allKeys[ityp]->erase(new_end, endKey(ityp));
   }

   template <class Predicate>
   void filterKeyAndFct(Predicate pred, const std::string& ityp = "default")
   {
      femKeys fkeys2;
      femFcts ffcts2;
      iterFct itf = beginFct(ityp);
      iterKey itk = beginKey(ityp);
      for (; (itk != endKey(ityp)) && (itf != endFct(ityp)); ++itk, ++itf)
      {
         if (pred(*itk))
         {
            fkeys2.push_back(*itk);
            ffcts2.push_back(*itf);
         }
      }
      allKeys[ityp]->clear();
      allFcts[ityp]->clear();
      *(allKeys[ityp]) = fkeys2;
      *(allFcts[ityp]) = ffcts2;
   }

   void setKeys(AOMD::mEntity* e, spacePtr f, const std::string& appro = "default");

   template <class iterFctSpace>
   void setKeys(AOMD::mEntity* e, iterFctSpace it, iterFctSpace last, const std::string& ityp = "default")
   {
      if (e_current != e)
      {
         clear();
         e_current = e;
      }
      femKeys* keys = KeysFor(ityp);
      keys->clear();
      for (; it != last; ++it)
      {
         addKeys(*it, ityp);
      }
   }

   template <class iterFctSpace>
   void setKeysAndFcts(AOMD::mEntity* e, iterFctSpace it, iterFctSpace last, const std::string& ityp = "default")
   {
      // setEntity(e);
      if (e_current != e)
      {
         clear();
         e_current = e;
      }
      femKeys* keys = KeysFor(ityp);
      keys->clear();
      femFcts* appro = FctsFor(ityp);
      appro->clear();
      for (; it != last; ++it)
      {
         addKeysAndFcts(*it, ityp);
      }
   }

   void setKeysAndFcts(AOMD::mEntity* e, spacePtr f, const std::string& appro = "default");

   AOMD::mEntity* getEntity() { return e_current; }

  private:
   AOMD::mEntity* e_current;

   typedef std::map<std::string, femFcts*> allFcts_t;
   allFcts_t allFcts;
   typedef std::map<std::string, femKeys*> allKeys_t;
   allKeys_t allKeys;
   femFcts* FctsFor(const std::string& ityp = "default");
   femKeys* KeysFor(const std::string& ityp = "default");

   void addKeys(spacePtr f, const std::string& appro = "default");
   void addKeysAndFcts(spacePtr f, const std::string& appro = "default");

   void clear();
   //  void setEntity(AOMD::mEntity* e);
};

}  // namespace xfem

#endif
