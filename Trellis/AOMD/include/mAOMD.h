/****************************************************************************

   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of the Algorithm-Oriented Mesh Database (AOMD) written
   and maintained by the Scientific Computation Research Center (SCOREC) at
   Rensselaer Polytechnic Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.

   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA

*****************************************************************************/
#ifndef _H_AOMD_Util_
#define _H_AOMD_Util_

#include <iosfwd>
#include <iostream>
#include <map>
#include <string>
#include <utility>


namespace AOMD
{
class mMesh;
class AOMD_LoadBalancerCallbacks;

void classifyUnclassifiedVerices(mMesh *);

/**
 AOMD is a singleton i.e. a class having only
 one occurence for a given heap.
 This class contains all global AOMD
 functionalities like import and export files
*/
class AOMD_Util
{
   unsigned int _ATTFMOD, _ATTEMOD, _ATTPARENT, _ATTID, _ATTSIZE, _ATTDN, _ATT1, _ATT2, _ATT3;
   unsigned int _ATTPARAMETRIC;
   unsigned int _ATTWEIGHT;
   int currentProcessor;
   std::map<std::string, unsigned int> attachableDataIds;
   std::map<unsigned int, std::string> attachableDataIds_rev;
#ifdef TSTT_
   std::map<unsigned int, int> attachableDataIds_type;
#endif
  public:
   inline unsigned int getParametric() const { return _ATTPARAMETRIC; }
   inline unsigned int getParent() const { return _ATTPARENT; }
   inline unsigned int getId() const { return _ATTID; }
   inline unsigned int getDn() const { return _ATTDN; }
   inline unsigned int getSize() const { return _ATTSIZE; }
   inline unsigned int getFmod() const { return _ATTFMOD; }
   inline unsigned int getEmod() const { return _ATTEMOD; }
   inline unsigned int getAtt1() const { return _ATT1; }
   inline unsigned int getAtt2() const { return _ATT2; }
   inline unsigned int getAtt3() const { return _ATT3; }
   inline unsigned int getWeight() const { return _ATTWEIGHT; }

   inline void setCurrentProcessor(int p) { currentProcessor = p; }
   /// Get the only instance
   static AOMD_Util *Instance();

   void partition_and_export(const char *baseName, mMesh *theMesh, AOMD_LoadBalancerCallbacks &lb);
   /// export a mesh, check the extension of the files
   void ex_port(const char *, const mMesh *, bool reduce_to_minimum = false);
   /// export a mesh, check the extension of the files
   void print_topology(const char *, mMesh *);
   /// import a mesh, check the extension of the file
   void import(const char *, mMesh *);
   /// specific to sms with <stdio>
   void importSmsFile(const char *, mMesh *);

   void import_oneLevel(const char *, mMesh *);
#ifndef SIM
   void importSmsFile_oneLevel(const char *, mMesh *);
#endif
#ifdef ACIS
   void importAOMDFile_oneLevel(const char *fName, mMesh *theMesh);
   void importAOMDFile(const char *fName, mMesh *theMesh);
#endif
   /// another format
   void importMichelGrid(std::istream &, mMesh *);
   /// yet another one, perhaps the "native" AOMD
   void importDGFile(std::istream &, mMesh *);
   void importDGFile(const char *, mMesh *);
   /// VTU meshes
   void importVTUFile(std::ifstream &, mMesh *);
   void importVTUFile(const char *, mMesh *);
   /// specific to gambit (.neu files)
   void importGambitFile(const char *, mMesh *);
   /// Export functions
   void exportSmsFile(const char *, const mMesh *);
   /// Export functions
   void exportGmshFile(const char *, const mMesh *);
   /// Export functions
   void exportGmshFile(std::ostream &, const mMesh *);
   /// Export functions

   /// Export functions

   void exportDGFile(std::ostream &, const mMesh *, bool reduce_to_minimum = false);
   void exportEnsightFile(const char *, const mMesh *, bool reduce_to_minimum = false);
   void exportVTUFile(const char *, const mMesh *, bool reduce_to_minimum = false);
   void exportICIFile(const char *, const mMesh *, bool reduce_to_minimum = false);
   /// Map from id's to tags
   void setMappingFile(const char *fileName, mMesh *);
   int getTag(int dim, int id);
   inline unsigned int newMeshDataId(const char *tag);
   inline unsigned int lookupMeshDataId(const char *tag);
#ifdef TSTT_
   inline unsigned int newMeshDataId(const char *tag, int type);
   inline int typeMeshDataId(unsigned int);
   inline std::string nameMeshDataId(unsigned int);
   inline int lookupMeshDataId(const char *tag, unsigned int *);
   inline int lookupMeshDataId(unsigned int);
#endif
   inline void deleteMeshDataId(unsigned int id);
   void BuildGTopology(mMesh *m);
   void BuildNullModel(mMesh *m);
   /// Constructor and destructors should be private but I'm fed up of this
   /// c++ compiler always shouting
   ~AOMD_Util();
   static void resolveAssociations(mMesh *theMesh, std::multimap<int, int> &vertexAssociations);

  private:
   /// Constructor and destructors are private
   AOMD_Util();
   /// The only instance
   struct object_creator
   {
      object_creator() { Instance(); }
      inline void do_nothing() const {}
   };
   static object_creator create_object;
   // static AOMD_Util instance;
};

#ifdef TSTT_

inline int AOMD_Util::lookupMeshDataId(const char *tag_name, unsigned int *tag_id)
{
   std::map<std::string, unsigned int>::const_iterator it = attachableDataIds.find(std::string(tag_name));
   if (it == attachableDataIds.end()) return 0;
   *tag_id = (*it).second;
   return 1;
}

inline int AOMD_Util::lookupMeshDataId(unsigned int tag_id)
{
   if (attachableDataIds_type.find(tag_id) == attachableDataIds_type.end())
      return 0;
   else
      return 1;
}

inline int AOMD_Util::typeMeshDataId(unsigned int tag_id)
{
   if (attachableDataIds_type.find(tag_id) == attachableDataIds_type.end())
   {
      std::cerr << "No tag found with tag_id " << tag_id << std::endl;
      throw 1;
   }
   return attachableDataIds_type[tag_id];
}

inline std::string AOMD_Util::nameMeshDataId(unsigned int tag_id)
{
   if (attachableDataIds_rev.find(tag_id) == attachableDataIds_rev.end())
   {
      std::cerr << "No tag found with tag_id " << tag_id << std::endl;
      throw 1;
   }
   return attachableDataIds_rev[tag_id];
}

inline unsigned int AOMD_Util::newMeshDataId(const char *tag, int type)
{
   //    std::cout<<"AOMD_Util::newMeshDataId("<<tag<<", "<<type<<")\n";
   if (attachableDataIds.find(tag) != attachableDataIds.end())
   {
      std::cerr << "AOMD cannot create an existing tag \"" << tag << "\"\n";
      throw 1;
   }
   if (attachableDataIds.empty())
   {
      attachableDataIds_rev[0] = std::string(tag);
      attachableDataIds_type[0] = type;
      attachableDataIds[tag] = 0;
      //        std::cout<<"newMeshDataId creates pMeshDataId "<<0<<" ("<<tag<<")\n";
      return 0;
   }
   else
   {
      unsigned int biggest = (*(--attachableDataIds_rev.end())).first;
      attachableDataIds_rev[biggest + 1] = tag;
      attachableDataIds_type[biggest + 1] = type;
      attachableDataIds[tag] = biggest + 1;
      //        std::cout<<"newMeshDataId creates pMeshDataId "<<biggest+1<<" ("<<tag<<")\n";
      return biggest + 1;
   }
}

#endif

// nico added declaration needed by client code
inline unsigned int AOMD_Util::lookupMeshDataId(const char *tag)
{
   std::map<std::string, unsigned int>::const_iterator it = attachableDataIds.find(std::string(tag));
   if (it == attachableDataIds.end()) return newMeshDataId(tag);
   return (*it).second;
}

inline unsigned int AOMD_Util::newMeshDataId(const char *tag)
{
   if (attachableDataIds.find(tag) != attachableDataIds.end())
   {
      std::cerr << "AOMD cannot create an existing tag \"" << tag << "\"\n";
      throw 1;
   }
   if (attachableDataIds.empty())
   {
      attachableDataIds_rev[0] = tag;
#ifdef TSTT_
      attachableDataIds_type[0] = -1;
#endif
      attachableDataIds[std::string(tag)] = 0;
      //        std::cout<<"newMeshDataId creates pMeshDataId "<<0<<" ("<<tag<<")\n";
      return 0;
   }
   else
   {
      unsigned int biggest = (*(--attachableDataIds_rev.end())).first;
      attachableDataIds_rev[biggest + 1] = tag;

#ifdef TSTT_
      attachableDataIds_type[biggest + 1] = -1;
#endif
      attachableDataIds[tag] = biggest + 1;
      //        std::cout<<"newMeshDataId creates pMeshDataId "<<biggest+1<<" ("<<tag<<")\n";
      return biggest + 1;
   }
}

inline void AOMD_Util::deleteMeshDataId(unsigned int id)
{
   attachableDataIds.erase(attachableDataIds_rev[id]);
   attachableDataIds_rev.erase(id);
#ifdef TSTT_
   attachableDataIds_type.erase(id);
#endif
}

}  // namespace AOMD

#endif
