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
/* 
   MessageOut.h

   Class definition header file for framework warning/errors

   Parallel Adaptive Framework

   Jim Teresco

   Created: Thu Apr 10 11:42:38 CDT 1997, Jim Teresco

   $Id: MessageOut.h,v 1.4 2005/02/02 00:08:57 acbauer Exp $
   $Log: MessageOut.h,v $
   Revision 1.4  2005/02/02 00:08:57  acbauer
   *** empty log message ***

   Revision 1.3  2003/04/04 10:39:52  remacle
   *** empty log message ***

   Revision 1.2  2002/09/06 19:27:51  xli
   The version deso not depend on any other models
   1. attributes are taken out
   2. inverse classifications are taken out
   3. needed utilities are moved into

   Revision 1.1.2.1  2002/09/04 14:44:06  xli
   *** empty log message ***

   Revision 1.1  1998/03/17 18:18:29  mbeall
   added


*/

#ifndef _H_MessageOut
#define _H_MessageOut

#include <iostream>
#include <fstream>
#include "SString.h"

class MessageOut {
public:
  MessageOut();
  ~MessageOut();
  enum Action {
    Ignore=0,
    Stdout,
    Stderr,
    File 
  };
  enum Level {
    InfoMsg=0,
    DebugMsg,
    WarningMsg,
    ErrorMsg,
    InternalErrorMsg
  };
  inline void SetInfoAction(Action a) { message_action[InfoMsg]=a; }
  inline void SetDebugAction(Action a) { message_action[DebugMsg]=a; }
  inline void SetWarningAction(Action a) { message_action[WarningMsg]=a; }
  inline void SetErrorAction(Action a) { message_action[ErrorMsg]=a; }
  inline void SetInternalErrorAction(Action a) { 
    message_action[InternalErrorMsg]=a;
  }
  void SetFilename(char *filename);

  std::ostream &Message(Level severity);
  void EndMessage(Level severity);

private:
  Action message_action[5];
  SString types[5];
  SString logfilename;
  std::ofstream output;
};

// the one default MessageOut we will instantiate in MessageOut.cc
extern MessageOut FWMessages;

// macros for easy access to the default

// Note: to turn off Informational messages at compile time, define
// the symbol MSGOUT_NOINFO on the compile command line or in source
// before including MessageOut.h

#ifdef MSGOUT_NOINFO
#define Info(x)
#define InfoInOrder(x)
#define InfoMaster(x)
#else
// Info(): Informational messages to be printed on all procs in no particular
// order and without line numbers and filenames
#define Info(x) { FWMessages.Message(MessageOut::InfoMsg) << x << endl; FWMessages.EndMessage(MessageOut::InfoMsg); }

// InfoInOrder(): Same as Info() except when in parallel, messages are
// printed in order by processor rank
// InfoMaster(): Same as Info() except when in parallel, message is printed
// only by the master processor
#endif


// Note: to turn off Informational messages at compile time, define
// the symbol MSGOUT_NODEBUG on the compile command line or in source
// before including MessageOut.h

#ifdef MSGOUT_NODEBUG
#define Debug(x)
#define DebugInOrder(x)
#define DebugMaster(x)
#else
// Debug(): Print a debug message, including line number and filename
// which produced it.
#define Debug(x) { FWMessages.Message(MessageOut::DebugMsg) << "(" << __FILE__ << ")[" << __LINE__ << "] " << x << std::endl; FWMessages.EndMessage(MessageOut::DebugMsg); }

// DebugInOrder(): Same as Debug() except messages are printed in processor
// rank order when in parallel
// DebugMaster(): Same as Debug() except messages are printed only by the
// master processor when in parallel

#endif

// Warning(): Print a warning message, including line number and filename
// which produced it.
#define Warning(x) { FWMessages.Message(MessageOut::WarningMsg) << "(" << __FILE__ << ")[" << __LINE__ << "] " << x << std::endl; FWMessages.EndMessage(MessageOut::WarningMsg); }

// Error(): Print an error message, including line number and filename
// which produced it.
#define Error(x) { FWMessages.Message(MessageOut::ErrorMsg) << "(" << __FILE__ << ")[" << __LINE__ << "] " << x << std::endl; FWMessages.EndMessage(MessageOut::ErrorMsg); }

// InternalError(): Print an error message, including line number and filename
// which produced it.  To be used internal to libraries for errors which
// should not occur in production code
#define InternalError(x) { FWMessages.Message(MessageOut::InternalErrorMsg) << "(" << __FILE__ << ")[" << __LINE__ << "] " << x << std::endl; FWMessages.EndMessage(MessageOut::InternalErrorMsg); }

#endif
