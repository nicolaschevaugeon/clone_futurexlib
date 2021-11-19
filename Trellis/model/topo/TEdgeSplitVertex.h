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
#ifndef H_TEdgeSplitVertex
#define H_TEdgeSplitVertex

#include "GVertex.h"

/** A vertex created during the spliting of an existing edge. The vertex
  lies on the original edges at a given parameter location. */
class TEdgeSplitVertex : public GVertex {
public:
  TEdgeSplitVertex(SGModel *model, int tag, GEdge *e, double param);
  ~TEdgeSplitVertex() override;

  GVPoint point() const override;
  int containsPoint(const AOMD::SPoint3 &pt) const override;

  GeoRep * geometry() override;

  double tolerance() const override;

protected:
  GEdge * SplitEdge;
  double Param;
};

#endif
