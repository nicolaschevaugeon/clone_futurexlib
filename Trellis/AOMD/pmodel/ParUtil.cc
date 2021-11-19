/******************************************************************************** 

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

*********************************************************************************/
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
//#ifdef SUN4
//#include <sys/varargs.h>
//#else
#include <cstdarg>
//#endif
#include "ParUtil.h"
#ifndef MSC
#include <sys/time.h>
#endif

namespace AOMD {

ParUtil* ParUtil::Instance()
{
  if(!instance)
    {
      instance = new ParUtil;
    }
  return instance;
}

ParUtil::~ParUtil() 
{
  Msg(INFO,"Finalizing...\n");
  Finalize();
}

ParUtil::ParUtil() 
{
}

void ParUtil::init(int &argc, char **&argv) {


}

double ParUtil::wTime() const
{


#ifndef MSC
  struct timeval tp;
  struct timezone tzp;
  double timeval;
  
  gettimeofday(&tp,&tzp);
  
  timeval = (double) tp.tv_sec;
  timeval = timeval + (double) ((double) .000001 * (double) tp.tv_usec);
  
  return(timeval);
#else
  return 0;
#endif

}

std::string ParUtil::processorName() const
{

  return std::string("localhost");

}

void ParUtil:: Msg(ParUtil::MessageLevel level, const char *fmt, ...)
{ 
  char buff[1024];
  va_list  args;
  va_start (args, fmt);
  vsprintf(buff, fmt, args);
  va_end (args);

  switch(level)
    {
    case DEBUG1:
      if(vl > 1)
	{
	  /* char logname[256];
	  sprintf(logname,"log-proc%d-%s.dat",myrank,procName);
	  log = fopen (logname,"a");	  
	  fprintf(log,"%s",buff);
	  fclose(log);*/
	}
      break;
    case DEBUG2:
      if(vl > 2)
	{
	  /*	  char logname[256];
	  sprintf(logname,"log-proc%d-%s.dat",myrank,procName);
	  log = fopen (logname,"a");	  
	  fprintf(log,"%s",buff);
	  fclose(log);*/
	}
      break;
    case INFO:
      if(vl >= 0 && master())
	{
	  //	  fprintf(log,"%s",buff);
	  fprintf(stdout,"%s",buff);
	  //	  fflush(log);
	}
      if(vl > 2)
	{
	  //	  fprintf(log,"%s",buff);
	  //	  fflush(log);
	}
      break;
    case WARNING:
      fprintf(stdout,"Processor %d AOMD WARNING : %s",rank(),buff);
      fflush(stdout);
      break;
    case ERROR:
      fprintf(stdout,"AOMD FATAL ERROR : %s",buff);
      fflush(stdout);
      Abort();
      break;
    }
}

void ParUtil::Abort()
{

  abort();

}

void ParUtil::Barrier(int line, const char *fn)
{

}

void ParUtil::Finalize()
{

}

ParUtil* ParUtil::instance = nullptr;
int ParUtil:: buffer_size = 1024;

} // end of namespace
