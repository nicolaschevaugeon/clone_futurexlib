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

#ifndef H_GInfo
#define H_GInfo

/** Defintion of various topological types. */
class TopoType{
public:
  enum Value{
    Vertex,
    Edge,
    Face,
    Region,
    VertexUse,
    EdgeUse,
    LoopUse,
    FaceUse,
    Shell,
    Model,
    Other
    };
};

/** Definition of various geometric types. */
class GeomType{
public:
  enum Value{
    Unknown,
    Point,
    Line,
    Circle,
    Ellipse,
    ParametricCurve,
    Plane,
    Nurb,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    ParametricSurface
    };
};

/** Definition of underlying representation types for entities. */
class RepType{
public:
  enum Value{
    Unknown,
    Parametric,
    Mesh
  };    
};

#ifdef False
#undef False
#endif

#ifdef True
#undef True
#endif

/** Definition of logical return values. */
class Logical{
public:
  enum Value{
    False = 0,
    True,
    Unknown
    };
};
    
#endif
