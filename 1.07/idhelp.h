

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __IDHELP_H__
#define __IDHELP_H__

#include <vector>
#include <set>
#include "plink.h"

using namespace std;

class IDFile;

class IDField
{
 public:
  
  string name;
  string jointName;
  int width;

  bool null;
  bool attribute;
  bool equiv;
  bool joint;
  
  
  // Equivalence IDs
  map<string,string > eqid;
  set<string> aliasList;
  set<string> masterList;
  
  
  IDField()
    {
      attribute = false;
      null = false; 
      equiv = false;
      joint = false;
      width = 4;
    }
  
  bool operator< (const IDField & b) const
    {
      if ( name < b.name ) 
	return true;
      return false;
    }

  bool operator== (const IDField & b) const
    {
      return name == b.name; 
    }
};


class IDFile 
{
public:

  string filename;

  int uniqFieldCount;
  
  bool hasHeader;
  
  set<string> missingValues;
  string delimit;
  string alias_delimit;

  vector<IDField*> fields;
  vector< vector<int> > joint;
  vector< vector<int> > equiv;
  map<string,string> injections;

  IDFile() 
    { 
      hasHeader = false;
      missingValues.insert(".");
      uniqFieldCount = 0;
      delimit = " ";
      alias_delimit = ",";
    }
};

class IDValue
{
public:
  
  IDValue()
  {
    value = "";
    field = NULL;
    jointValue = "";
  }
  
  void updateAlias()
  {
    map<string,string>::iterator f = field->eqid.find( value );
    if ( f != field->eqid.end() )
      value = f->second;
  }

  IDField * field;
  
  string jointValue;
  
  // This is main default value
  string value;

  
  bool operator< (const IDValue & b) const
  {
    if ( field->name < b.field->name )
      return true;
    else if ( field->name > b.field->name )
      return false;
    if ( field->joint && b.field->joint )
      return jointValue < b.jointValue;
    else
      return value < b.value; 

  }  


  bool operator!= (const IDValue & b) const
  {
    return ! ( *this == b );
  }


  bool operator== (const IDValue & b) const
  {
        
    // These must refer to the same field to be 
    // comparable

    if ( field->name != b.field->name ) 
      return false;
    
    if ( field->joint && b.field->joint )
      {
	return jointValue == b.jointValue;
      }
    else
      {
	return value == b.value;
      }
    
  }

  bool singleMatch(const IDValue & b) const
    {
      // Do not take joint field values into 
      // account here
      if ( field->name != b.field->name ) 
	return false;
      return value == b.value;
    }


};


class IDGroup
{
public:

  vector<IDValue*> values;
  IDFile * file;

  bool resolved;

  IDGroup()
  {
    resolved = false;
  }

  void display()
  {
    cout << "File = " << file->filename 
	 << ", resolved = " 
	 << resolved << "\n";
    
    for (int k=0; k<values.size(); k++)
      cout << "\t" << values[k]->field->name 
	   << " = " << values[k]->value 
	   << " j= " << values[k]->jointValue 
	   << " joint=" << values[k]->field->joint 
	   << " attrib=" << values[k]->field->attribute 
	   << " null=" << values[k]->field->null 
	   << "\n";
    cout << "\n";
  }

};



class IDHelper
{
 public:

  // Files containing IDs
  vector<IDFile> files;
  
  // The actual IDs we are matching on
  set<IDField> fields;
  map<string,IDField*> fieldMap;
  set<IDField>::iterator iField;
  
  // The basic data we read in, then try to resolve into
  // a smaller set
  vector<IDGroup*> idgroup;
    
  // The lookup table
  map<IDValue,set<IDGroup*> > idmap;
   
  map<string,IDFile> dict;
 
  set< set<string> > uniqueFields;
  map<string,set<string> > jointMap; 
  vector< set<IDField*> > jointField;
  vector< vector<IDField*> > jointOrder;

  set<string> attribFields;

  // Functions

  // Main wrapper
  void idHelp();

  void idReplace();
  void idMatch();
  void idDump();
  void idListAlias();

  // Helper functions

  // Populate joint values for a given set
  void setJointValues( set<IDValue> & ); 
  void setJointValues( IDGroup * ); 

  // Map string ID1+ID2=V1+V2,ID3=V3 
  map<IDField*,set<IDValue> > parseQuery(string );

  // Find unique matching observation
  IDGroup * findUniqueIndividual( set<IDValue> & );

  // Does this person match this template?
  bool matchIndividual(IDGroup * group, map<IDField*,set<IDValue> > & );
  
  // Find all matching observations
  set<IDGroup*> findAllIndividuals( map<IDField*,set<IDValue> > & );

  // Set an alias
  bool setAlias(IDField*, string, int, map<IDField*,string>&);
};



#endif
