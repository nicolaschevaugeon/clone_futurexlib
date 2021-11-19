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
#ifndef H_AttachableData
#define H_AttachableData

/** A mixin class to be able to attach arbitrary data to an object
ala the original meshdb */

class AttachDataManager;
class AttachDataId;
class AttachDataBase;

class AttachableData {
public:
  void setDataP(const char *, void *);
  void setDataI(char *, long);
  void removeData(char *);
  void *dataP(const char*) const;
  long int dataI(char*) const;
  int modifyDataP(char *tag, void * data);
  int modifyDataI(char *tag, long int data);

  // the section added by Jie starts
  void attachDataInt(AttachDataId *id, int value);
  void attachDataDbl(AttachDataId *id, double value);
  void attachDataPtr(AttachDataId *id, void * value);

  void deleteData(AttachDataId *id, int warn_exist=1);

  int getDataInt(AttachDataId *id, int *value) const;
  int getDataDbl(AttachDataId *id, double *value) const;
  int getDataPtr(AttachDataId *id, void **value) const;

  void modifyDataInt(AttachDataId *id, int value);
  void modifyDataDbl(AttachDataId *id, double value);
  void modifyDataPtr(AttachDataId *id, void * value);
  // the section added by jie ends

protected:
  AttachableData();
  virtual ~AttachableData();
private:
  struct TList *Data; 

  // the section added by jie starts
protected:
  void callCallbacks(int, void *);
  virtual AttachDataManager *dataManager() const=0;

private:
  void attachData(AttachDataBase *);
  AttachDataBase *getData(AttachDataId *) const;
  AttachDataBase * detachData(AttachDataId *, int warn_exist=1);

  // pointer to embedded list of attached data items
  AttachDataBase *Data1; 

  // the section added vy jie ends
};

#endif
