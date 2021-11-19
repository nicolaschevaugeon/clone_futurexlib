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
#ifndef H_TopoModel
#define H_TopoModel

#include "SGModel.h"

/** A model consisting only of topology. The TopoModel may be a copy
  of an existing model in which case geometric queries will be
  forwarded on to the copied model. Various operations for modifying
  the copied topology are supported, in particular, augmentation of
  periodic faces to introduce a seam edge. If there is no source
  model for the copy, geometric queries cannot be made.
  */
class TopoModel : public SGModel
{
  public:
   /** Create a TopoModel that is a copy of an existing model. If
     source is given as null, create an empty model. */
   TopoModel(SGModel *source);
   TopoModel(const SString &name, SGModel *source = nullptr);

   // static void registerModeler();
   // static SGModel * create(const SString &name, int loadOnly);

   SString modeler() const override;

   AOMD::SBoundingBox3d bounds() const override;
   void setBounds(const AOMD::SBoundingBox3d &box);
   double tolerance() const override;

   virtual void finalize();

   /** Modify this copy of the model to introduce seam edges on any
     periodic faces. */
   void augmentPeriodicFaces();

  protected:
   GVertex *createVertex(istream &in) override;
   GEdge *createEdge(istream &in) override;
   GFace *createFace(istream &in) override;
   GRegion *createRegion(istream &in) override;
   void readBounds(istream &in) override;

   SGModel *RealModel;
   AOMD::SBoundingBox3d bbox;
};

#endif
