/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/

#ifndef H_NullModel
#define H_NullModel

#include "SBlock.h"

#include "SGModel.h"

class SString;
class NullRegion;
class NullFace;
class NullEdge;
class NullVertex;
namespace AOMD{
  class SBoundingBox3d;
}
using namespace std;

/** A "fake" model. This so-called model is only a model in the sense
  that it is comprised of a bunch of model entities. The entities have
  no topological adjacency information and no geometric information.
  Despite this fact, it's rather useful in some cases. An example would
  be when you want to load a mesh without loading a "real" model, but
  still want to have classification information available.
  */
class NullModel : public SGModel {
public:
  NullModel(const SString &name); 
  NullModel(); 
  static void setupNull();
  static void registerModeler();
  static SModel * create(const SString &name, int loadOnly, const SSList<SModel*> * const);

  ~NullModel() override = default;;

  /** This will always return a region. If one doesn't exist with the
    given tag, it will be created. */
  GRegion * regionByTag(int n) const override;
  GFace * faceByTag(int n) const override;
  GEdge * edgeByTag(int n) const override;
  GVertex * vertexByTag(int n) const override;

  GRegion * regionByID(int n) const override;
  GFace * faceByID(int n) const override;
  GEdge * edgeByID(int n) const override;
  GVertex * vertexByID(int n) const override;

  SString modeler() const override;

  AOMD::SBoundingBox3d bounds() const override;
  double tolerance() const override;

  void setGeomTolerance(double) override;

protected:
  void buildFromGeom();

  GVertex *createVertex(istream &in) override;
  GEdge *createEdge(istream &in) override;
  GFace *createFace(istream &in) override;
  GRegion *createRegion(istream &in) override;

private:

  static int nullSetup;
};

extern "C" void GM_setupNull();
extern "C" SGModel * GM_createFromNothing(const char *name);

#endif

