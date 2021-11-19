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
   MessageOut.cc

   Framework warning/error handler

   Parallel Adaptive Framework

   Jim Teresco

   Created: Thu Apr 10 12:21:58 CDT 1997, Jim Teresco

   $Id: MessageOut.cc,v 1.4 2005/02/02 00:08:57 acbauer Exp $
   $Log: MessageOut.cc,v $
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

   Revision 1.1  1998/03/17 18:18:38  mbeall
   added


*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "SString.h"
#include "MessageOut.h"

// instantiate one which we will use as the default
MessageOut FWMessages;

MessageOut::MessageOut() {

  message_action[InfoMsg]=Stdout;
  message_action[DebugMsg]=Stdout;
  message_action[WarningMsg]=Stdout;
  message_action[ErrorMsg]=Stderr;
  message_action[InternalErrorMsg]=Stderr;

  logfilename=SString("messages.log");
  //output.open(logfilename,ios::out);
  //output << "Log file opened...\n";
  //output.close();

  types[InfoMsg]="";
  types[DebugMsg]="Debug:";
  types[WarningMsg]="Warning:";
  types[ErrorMsg]="Error:";
  types[InternalErrorMsg]="Internal Error:";
}

MessageOut::~MessageOut() {

  //if (logfilename) delete logfilename;
}

void MessageOut::SetFilename(char *filename) {

  logfilename=SString(filename);
  output.open(logfilename,std::ios::out);
  output << "Log file " << logfilename << " opened...\n";
  output.close();
}

std::ostream & MessageOut::Message(Level severity) {

  switch (message_action[severity]) {
  case Ignore:
    output.open("/dev/null",std::ios::out);
    return output;
  case Stdout:
    return std::cout << types[severity];
    break;
  case Stderr:
    return std::cerr << types[severity];
    break;
  case File: {
    output.open(logfilename,std::ios::app);
    return output << types[severity];
    break;
  }
  default : throw;
  }
}

void MessageOut::EndMessage(Level severity) {

  switch (message_action[severity]) {
  case File:
    output.flush();
  case Ignore:
    output.close();
    break;
  case Stdout:
    std::cout.flush();
    break;
  case Stderr:
    std::cerr.flush();  // not really needed I suppose
    break;
  default : throw;
  }
}

