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
#ifndef _H_ParUtil
#define _H_ParUtil
#include <string>


namespace AOMD {

/**
   ParUtil is a Singleton. It gives some
   general services for parallel implementation.
*/

class ParUtil {
  ParUtil();
public:
    /// Constructor and destructors should be private but I'm fed up of this
    /// c++ compiler always shouting
  ~ParUtil();
  /// Message severity level
  enum MessageLevel {DEBUG1,DEBUG2,INFO,WARNING,ERROR};
  /// returne the only instance
  static ParUtil* Instance();
  /// initialization, needed for mpi and autopack
  void init(int &argc, char **&argv);
  /// adds a barrier
  void Barrier(int, const char*);
  /// close communications
  void Finalize();
  /// compute wall time
  double wTime () const;
  /// gets the processor name
  std::string processorName() const;
  /// set the verbosity level (0,1,2,3)
  inline void setVertbosityLevel(int i){vl = i;}
  /// prints a message, same format as printf
  void Msg(MessageLevel lev, const char *fmt, ...);
  /// abort a calculation
  void Abort();
  static void setBufferSize (int b) {buffer_size = b;}
  /// gets the processor id
  inline int rank() { return 0; }
  /// gets the number of processors
  inline int size() { return 1; }
  /// tells if it's processor 0
  inline int  master() { return 1; }

  inline void   resetBarrierSensor() { timeSpentOnBarriers = 0; }
  inline double getBarrierSensor()   { return timeSpentOnBarriers; }
private:
  static ParUtil *instance;
  int vl;
  std::string procName;
  double timeSpentOnBarriers;
  static int buffer_size;
};

} // end of namespace

#endif
