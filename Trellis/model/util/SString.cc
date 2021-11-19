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
#include "SString.h"

SString::SString()
: len(0)
{ 
  str = new char[1];
  str[0] = '\0';
}

SString::SString(const SString& s)
{
  len = s.len;
  str = new char[len + 1];
  strcpy(str, s.str);
}

SString::SString(const char *s)
{
  len = strlen(s);
  str = new char[len + 1];
  strcpy(str, s);
}

SString::SString(int ln)
{
  len = ln;
  str = new char[len + 1];
  str[0] = '\0';
}

SString::~SString()
{
  delete [] str;
}

std::istream& operator>>(std::istream& os, SString& s)
{
  char st[256];
  os >> st;
  s = st;	
  return os;
}

SString operator+(const SString &a, const SString &b)
{
  SString *s = new SString(a.len+b.len);
  strcpy(s->str,a.str);
  strcat(s->str,b.str);
  return *s;
}

SString& SString::operator=(const SString& s)
{
  len = s.len;
  delete [] str;
  str = new char[len + 1];
  strcpy(str, s.str);
  return *this;
}

SString& SString::operator=(const char *s)
{
  len = strlen(s);
  delete [] str;
  str = new char[len + 1];
  strcpy(str, s);
  return *this;
}

int SString::operator==(const SString &s) const
{
  return strcmp(str, s.str) ? 0 : 1;
}

int SString::operator!=(const SString &s) const
{
  return !operator==(s);
}

int SString::length() const
{ return len; }

