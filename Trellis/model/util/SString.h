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
#ifndef H_SString
#define H_SString

#include <string>
#include <cstring>

#include <iostream>

/** A character string */
class SString {

public:
  SString();
  SString(const SString& s);
  SString(const char *s);
  SString(int ln);
  ~SString();

  int operator==(const SString &s) const;
  int operator!=(const SString &s) const;
  
  friend std::ostream& operator<<(std::ostream& os, const SString& s)
    {return (os << s.str);}
  friend std::istream& operator>>(std::istream& os, SString& s);
  
  friend SString operator+(const SString &a, const SString &b);
  SString& operator=(const SString& s);
  SString& operator=(const char *s);
  
  operator const char*() const { return str; }
  
  operator char*() { return str; }
  
  int length() const;
protected:
  int len;
  char *str;
};

#endif     
