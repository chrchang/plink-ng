

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

// Rules:

// ID cannot have any spaces, tabs, commas (,) or plus (+) signs
// Attributes can have commas

#include "idhelp.h"
#include "options.h"
#include "helper.h"
#include "nlist.h"
#include <iomanip>

extern Plink * PP;


map<IDField*,set<IDValue> > IDHelper::parseQuery(string q)
{

  map<IDField*,set<IDValue> > lookup;

  // Convert a query string into a searchable form
  NList tlist(0);
  vector<string> ids = tlist.deparseStringList( q );

  // Consider each term (comma-delimited list)
  map<string,string> qmap;
  map<string,string> qjoint;

  for ( int i = 0 ; i < ids.size() ; i++)
    {
      string s = ids[i];
      if ( s.find("=") == string::npos )
	error("Query should be: ID=value or ATTRIB=value[,ATTRIB=value]");

      string fs = s.substr(0, s.find("="));
      string vs = s.substr(s.find("=")+1);

      // Is this a joint field?
      if ( fs.find("+") != string::npos )
	{
	  // We need to split both fs and vs: they must have the same
	  // number of entries on both sides

	  vector<string> flist =
	    tlist.deparseStringList( searchAndReplace( fs,"+",",") );
	  vector<string> vlist =
	    tlist.deparseStringList( searchAndReplace( vs,"+",",") );

	  if ( flist.size() != vlist.size() )
	    error("Joint query wrong: " + fs + " = " + vs);

	  for (int i=0; i<flist.size(); i++)
	    {
	      qmap.insert(make_pair(flist[i],vlist[i]));
	      qjoint.insert(make_pair(flist[i],vs));
	    }
	}
      else
	{
	  qmap.insert(make_pair(fs,vs));
	}
    }

  // Check that all named fields exist
  map<string,string>::iterator i = qmap.begin();
  while ( i != qmap.end() )
    {

      string fs = i->first;
      map<string,IDField*>::iterator f = fieldMap.find( fs );
      if ( f == fieldMap.end() )
	error("Cannot find field " + fs );

      IDField * thisField = f->second;

      IDValue t;

      t.field = thisField;
      t.value = i->second;
      if ( t.field->equiv )
	t.updateAlias();

      //////////////////////////////////////////////////
      // Is this query being framed as a joint field?

      map<string,string>::iterator k = qjoint.find( fs );
      if ( k != qjoint.end() )
	{
	  t.jointValue = k->second;
	  t.field->joint = true;
	}
      else
	{
	  t.field->joint = false;
	}

      map<IDField*,set<IDValue> >::iterator j = lookup.find( thisField );
      if ( j == lookup.end() )
	{
	  set<IDValue> ts;
	  ts.insert(t);
	  lookup.insert(make_pair(thisField,ts) );
	}
      else
	j->second.insert(t);

      ++i;
    }

  return lookup;
}

bool IDHelper::matchIndividual(IDGroup * group,
			       map<IDField*,set<IDValue> > & matchTemplate )
{

  bool match = true;

  // Compare
  //  map<IDField*, set<IDValue> > lookupValues;
  //  with this group

  int found = 0;

  for (int g=0; g<group->values.size(); g++)
    {

      IDField * f = group->values[g]->field;

      map<IDField*, set<IDValue> >::iterator mf = matchTemplate.find( f );

      if ( mf == matchTemplate.end() )
	continue;

      ++found;

      // The observed value

      IDValue * myValue = group->values[g];

      set<IDValue>::iterator j = mf->second.begin();
      bool matchField = false;
      while ( j != mf->second.end() )
	{

	  if ( *myValue == *j )
	    {
	      matchField = true;
	    }
	  ++j;
	}

      if ( ! matchField )
	return false;
    }

  // Did we get a look at all the match fields?
  if ( found != matchTemplate.size() )
    {
      return false;
    }

  // If still here, we must match
  return true;

}

void IDHelper::setJointValues( set<IDValue> & val )
{
  // Make a dummy ID group, with pointers to the originals
  IDGroup g;
  set<IDValue>::iterator i = val.begin();
  while ( i != val.end() )
    {
      g.values.push_back( (IDValue*)&(*i) );
      ++i;
    }

  // Use main joint-value update function
  setJointValues( &g );

  // We will add "jointValue" attribs to val, but no need
  // to remove any items

}



void IDHelper::setJointValues( IDGroup * group )
{


  // Find the joint fields and edit in values

  // Similar to the above function, except here if the joint values
  // are completely missing, then we need to remove the entire set
  // from the ID group

  for (int j = 0 ; j < jointField.size(); j++ )
    {

      set<IDField*> & jf = jointField[j];
      vector<IDField*> & jo = jointOrder[j];

      // Does this group contain one of these
      // joint fields?

      bool hasJoint = false;
      map<IDField*,int> mapback;
      for (int j = 0 ; j < group->values.size(); j++)
	{
	  if ( jf.find( group->values[j]->field ) != jf.end() )
	    {
	      hasJoint = true;
	      mapback.insert( make_pair( group->values[j]->field , j ) );
	    }
	}

      if ( ! hasJoint )
	continue;


      string jointValue = "";
      bool doneFirst = false;
      bool jointMissing = false;
      bool jointOneSeen = false;

      set<IDField*>::iterator k = jf.begin();
      while ( k != jf.end() )
	{

	  map<IDField*,int>::iterator mi = mapback.find( *k );
	  if ( mi == mapback.end() )
	    {
	      jointMissing = true;
	      if ( doneFirst )
		jointValue += "+.";
	      else
		{
		  jointValue += ".";
		  doneFirst = true;
		}
	    }
	  else
	    {
	      jointOneSeen = true;
	      if ( doneFirst )
		jointValue += "+" + group->values[ mi->second ]->value;
	      else
		{
		  jointValue += group->values[mi->second]->value;
		  doneFirst = true;
		}
	    }
	  ++k;
	}


      // Update values accordingly
      // Except, remove entirely missing values

      vector<bool> mask( group->values.size(), false);

      set<IDField*>::iterator k2 = jf.begin();

      while ( k2 != jf.end() )
	{
	  map<IDField*,int>::iterator mi = mapback.find( *k2 );
	  if ( mi != mapback.end() )
	    {
	      int j = mi->second;
	      if ( jointMissing && ! jointOneSeen )
		mask[j] = true;
	      else
		group->values[j]->jointValue = jointValue;
	    }
	  ++k2;
	}

      // Remove entirely missing values
      if ( jointMissing && ! jointOneSeen )
	{
	  vector<IDValue*> newValues = group->values;
	  group->values.clear();
	  for ( int i = 0 ; i < mask.size() ; i++)
	    {
	      if ( ! mask[i] )
		group->values.push_back( newValues[i] );
	    }
	}
    }

  return;

}



IDGroup * IDHelper::findUniqueIndividual( set<IDValue> & matchTemplate )
{

  set<IDGroup *> thisGroup;
  set<IDValue>::iterator k = matchTemplate.begin();

  while ( k != matchTemplate.end() )
    {

      // Hmm need to decide how to handle: is default joint or single lookup?
      // for IDValues?

      map<IDValue,set<IDGroup*> >::iterator i = idmap.find( *k );
      if ( i != idmap.end() )
	{
	  set<IDGroup*>::iterator j = i->second.begin();
	  while ( j != i->second.end() )
	    {
	      if ( ! (*j)->resolved )
		thisGroup.insert( *j );
	      ++j;
	    }
	}
      ++k;
    }

  if ( thisGroup.size() > 1 )
    error("Internal problem: found more than one match when expecting a unique match");

  return *thisGroup.begin();

}

set<IDGroup*> IDHelper::findAllIndividuals( map<IDField*,set<IDValue> > & matchTemplate )
{
  set<IDGroup*> matched;
  for ( int g = 0 ; g < idgroup.size(); g++ )
    {
      IDGroup * group = idgroup[g];
      if ( matchIndividual( group , matchTemplate ) )
	matched.insert( group );
    }
  return matched;
}


bool IDHelper::setAlias(IDField * myField, string val, int f, map<IDField*,string> & originalEquivalence)
{
  bool needToStore = true;

  map<IDField*,string>::iterator i = originalEquivalence.find( myField );
  if ( i == originalEquivalence.end() )
    {

      set<string>::iterator k = myField->aliasList.find( val );
      if ( k != myField->aliasList.end() )
	error("Alias specified more than once in "
	      + files[f].filename + " : "
	      + myField->name + " " + val );
      myField->masterList.insert( val );
      originalEquivalence.insert(make_pair( myField, val ));
    }
  else
    {

      // We need to add an equivalence value to this field?
      map<string, string>::iterator k = myField->eqid.find( val );

      // We only let each alias be specified once
      if ( k != myField->eqid.end() )
	error("Alias specified more than once in "
	      + files[f].filename + " : "
	      + myField->name + " " + val );

      k = myField->eqid.find( i->second );
      if ( k != myField->eqid.end() )
	error("Alias specified more than once in "
	      + files[f].filename + " : "
	      + myField->name + " " + i->second);

      if ( myField->masterList.find( val ) != myField->masterList.end() )
	error("Alias specified more than once in "
	      + files[f].filename + " : "
	      + myField->name + " " + val );


      // Keep track of this alias
      myField->eqid.insert( make_pair ( val , i->second ) );

      myField->aliasList.insert( val );

      // And now we do not need to store this particular value
      // as a distinct field

      needToStore = false;
    }

  return needToStore;
}


void IDHelper::idHelp()
{

  // Turn off the "-" initiated range-delimiting (i.e. so that
  // IDs can have hyphens)

  par::range_delimiter = " ";


  // Are we only performing a "simple match", performed without
  // a dictionary being specified? If so, jump there now, creating
  // the dictionary on the fly.

  if ( par::idhelp_match && par::idhelp_no_dict )
    {
      idMatch();
      return;
    }


  // 1. Read in dictionary, and make the fields and files

  // Contains files (and can be full path) and description of each field
  // {file name}   {col names } : { rules }
  // e.g.

  // ../files/id1.txt FID IID FID1 FID2 : uniq=FID,IID uniq=FID1,IID2
  // ../names/id.lst ID2 ID23 : missing=NA,---,0
  // ../names/id2.lst ID3 : equiv
  //

  // note: for "equiv" files, assume all IDs on same line are equivalent,
  // only one ID can be specified here


  PP->printLOG("ID helper, with dictionary [ " + par::idhelp_dictionary + " ]\n");
  checkFileExists( par::idhelp_dictionary );


  ifstream DICT( par::idhelp_dictionary.c_str() , ios::in );

  while ( ! DICT.eof() )
    {
      vector<string> tokens = tokenizeLine( DICT );
      // Needs atleast three fields:
      // filename, and two IDs

      if ( tokens.size() == 0 )
	continue;

      if ( tokens.size() < 2 || tokens[1] == ":" )
	error("Expecting at least 2 fields w/ 1 ID in every dictionary row\n");

      IDFile d;
      checkFileExists( tokens[0] );
      d.filename = tokens[0];

      int p = 1;

      while (1)
	{

	  if ( tokens[p] == ":" )
	    {
	      ++p;
	      break;
	    }

	  // Is this a new field name?
	  IDField f;
	  f.name = tokens[p];

	  // Set to skip this field?
	  if ( f.name == "." )
	    f.null = f.attribute = true;

	  iField = fields.find( f );

	  if ( iField == fields.end() )
	    {
	      fields.insert( f );
	      iField = fields.find( f );
	    }

	  // Track that this field was found in this file

	  IDField * ip = (IDField*)&( *iField );
	  d.fields.push_back( ip );
	  ++d.uniqFieldCount;

	  // Consider the next field in this file
	  if ( ++p == tokens.size() )
	    break;
	}

      bool seenJoint = false;
      bool seenAttrib = false;

      // Now read any rules
      if ( p < tokens.size() )
	{
	  while (1)
	    {
	      string cmd = tokens[p];

	      // Commands
	      // attrib=X,Y,Z
	      // joint=X,Y
	      // set:X=Y
	      // header
	      // missing=NA,-9
	      // alias-delimit=|

	      bool parsed = false;

	      if ( cmd.size() > 6 && cmd.substr(0,6) == "joint=" )
		{
		  parsed = true;
		  bool seenJoint = true;
		  string u = cmd.substr(6);
		  string jName = searchAndReplace(u,","," ");

		  // This should contain at least two ID fields
		  // These fields must always then appear together in
		  // any file that features at least one

		  NList tlist(0);
		  vector<string> ids = tlist.deparseStringList( u );
		  if ( ids.size() < 2 )
		    error("Problem with specification of : " + cmd );
		  set<string> t;
		  for (int i = 0 ; i < ids.size() ; i++)
		    {
		      t.insert(ids[i]);
		    }

		  for (int i = 0 ; i < ids.size() ; i++)
		    jointMap.insert( make_pair( ids[i] , t ) );

		  set<IDField*> jointf;
		  vector<IDField*> jointo;

		  for (int i = 0 ; i < ids.size() ; i++)
		    {
		      IDField f;
		      f.name = ids[i];
		      iField = fields.find( f );
		      if ( iField == fields.end() )
			{
			  error("Could not find field " + ids[i]
				+ " which is specified as joint");
			}

		      IDField * ip = (IDField*)&( *iField );
		      ip->joint = true;
		      ip->jointName = jName;

		      jointf.insert( ip );
		      jointo.push_back( ip );
		    }
		  jointField.push_back( jointf );
		  jointOrder.push_back( jointo );

		}


	      if ( cmd.size() > 7 && cmd.substr(0,7) == "attrib=" )
		{
		  parsed = true;
		  seenAttrib = true;
		  string u = cmd.substr(7);
		  NList tlist(0);
		  vector<string> ids = tlist.deparseStringList( u );
		  for (int i = 0 ; i < ids.size() ; i++)
		    attribFields.insert( ids[i] );
		}

	      if ( cmd.size() > 8 && cmd.substr(0,8) == "missing=" )
		{
		  parsed = true;
		  string u = cmd.substr(8);
		  NList tlist(0);
		  vector<string> ids = tlist.deparseStringList( u );
		  for (int i = 0 ; i < ids.size() ; i++)
		    d.missingValues.insert( ids[i] );
		}

	      if ( cmd.size() > 4 && cmd.substr(0,4) == "set:" )
		{
		  parsed = true;
		  if ( seenJoint )
		    error("Must specify set:X=Y before joint=X,Z");
		  if ( seenAttrib )
		    error("Must specify set:X=Y before attrib=X");

		  string u = cmd.substr(4);
		  bool okay = true;

		  if ( u.find("=") == string::npos )
		    okay = false;
		  else
		    {
		      string u1 = u.substr(0,u.find("="));
		      string u2 = u.substr(u.find("=")+1);
		      if ( u1.size() < 1 )
			okay = false;
		      if ( u2.size() < 1 )
			okay = false;

		      if ( okay )
			{
			  d.injections.insert(make_pair(u1,u2));
			  IDField f;
			  f.name = u1;
			  if ( f.name == "." )
			    error("Cannot set field value name to .");

			  iField = fields.find( f );

			  if ( iField == fields.end() )
			    {
			      fields.insert( f );
			      iField = fields.find( f );
			    }

			  // Track that this field was found in this file
			  IDField * ip = (IDField*)&( *iField );
			  d.fields.push_back( ip );
			  ++d.uniqFieldCount;

			}

		    }

		  if ( ! okay )
		    error("Badly formed set:X=V command");

		}

	      if( cmd == "header" || cmd == "hasHeader" )
		{
		  parsed = true;
		  d.hasHeader = true;
		}

// 	      if ( cmd.size() > 10 && cmd.substr(0,10) == "delimiter=" )
// 		{
// 		  string u = cmd.substr(10);
// 		  if ( u.size() != 1 )
// 		    error("Delimiters can only be a single character length");
// 		  d.delimit = u;
// 		}

	      if ( cmd.size() > 16 && cmd.substr(0,16) == "alias-delimiter=" )
		{
		  parsed = true;
		  string u = cmd.substr(16);
		  if ( u.size() != 1 )
		    error("Delimiters can only be a single character length");
		  if ( u == d.delimit )
		    error("Delimiter and alias delimiter cannot be the same value");
		  d.alias_delimit = u;
		}

	      if ( ! parsed )
		error("Could not parse the following rule in the ID dictionary: " + cmd );

	      if ( ++p == tokens.size() )
		break;
	    }
	}

      files.push_back(d);
      // Next line
    }

  DICT.close();

  PP->printLOG("Read " + int2str( fields.size() ) + " unique fields\n");

  set<IDField>::iterator i = fields.begin();
  while ( i != fields.end() )
    {
      fieldMap.insert(make_pair( i->name , (IDField*)&(*i) ));
      ++i;
    }


  /////////////////////////////////////////////////
  // Set flags for any attribute fields

  if ( attribFields.size() > 0 )
    {
      PP->printLOG("   Attribute fields: ");
      set<string>::iterator i = attribFields.begin();
      while ( i != attribFields.end() )
	{
	  if ( fieldMap.find( *i ) == fieldMap.end() )
	    error("Cannot find specified attribute "
		  + *i
		  + " -- please check your dictionary file");

	  // Set as an attribute field
	  fieldMap.find(*i)->second->attribute = true;

	  PP->printLOG( *i + " " );
	  ++i;
	}
      PP->printLOG("\n");
    }



  /////////////////////////////////////////////////
  // If set fields, check either joint or attrib

  for ( int f = 0 ; f < files.size() ; f++ )
    {
      map<string,string>::iterator i = files[f].injections.begin();
      while ( i != files[f].injections.end() )
	{
	  IDField * f = fieldMap.find( i->first )->second;
	  if ( ! ( f->joint || f->attribute ) )
	    error("Any set:field=value should be an attribute or a joint field");
	  ++i;
	}
    }

  ////////////////////////////////////////////////
  // If joint fields, check always all specified

  for (int j = 0 ; j < jointField.size(); j++ )
    {
      set<IDField*> & jf = jointField[j];
      for ( int f = 0 ; f < files.size() ; f++ )
	{
	  for ( int j = 0; j < files[f].fields.size(); j++)
	    {
	      if ( jointMap.find( files[f].fields[j]->name ) != jointMap.end() )
		{
		  map<string,set<string> >::iterator i = jointMap.find( files[f].fields[j]->name );
		  set<string> & ss = i->second;
		  set<string>::iterator is = ss.begin();
		  while ( is != ss.end() )
		    {
		      bool okay = false;
		      for ( int j = 0; j < files[f].fields.size(); j++)
			if ( *is == files[f].fields[j]->name )
			  okay =true;
		      if ( ! okay )
			error("Need to specify all joint fields in dictionary, [" + files[f].filename + " ]");
		      ++is;
		    }
		}
	    }
	}
    }

  if ( jointField.size() > 0 )
    {
      PP->printLOG("   Joint fields:");
      for (int j = 0 ; j < jointField.size(); j++ )
	{
	  set<IDField*> & jf = jointField[j];
	  set<IDField*>::iterator jiter = jf.begin();
	  PP->printLOG(" { ");
	  while ( jiter != jf.end() )
	    {
	      PP->printLOG( (*jiter)->name + " " );
	      ++jiter;
	    }
	  PP->printLOG(" }");
	}
      PP->printLOG("\n");
    }

  /////////////////////////////////
  // 2. Read in and index all IDs

  for ( int f = 0 ; f < files.size() ; f++ )
    {
      IDFile * file = &files[f];
      PP->printLOG("Reading [ " + files[f].filename + " ] with fields : ");
      for ( int j = 0 ; j < files[f].fields.size(); j++ )
	{
	  if (j>0 )
	    PP->printLOG(", ");
	  PP->printLOG( files[f].fields[j]->name );
	}


      ////////////////////////////
      // Find any equivalence sets

      bool foundEquiv = false;
      map<string, vector<int> > equiv;
      map<string, map<int,int> > equivMap;
      for ( int j = 0 ; j < files[f].fields.size(); j++)
	{
	  string n = files[f].fields[j]->name;
	  map<string, vector<int> >::iterator i = equiv.find( n );
	  if ( i == equiv.end() )
	    {
	      vector<int> t;
	      t.push_back(j);
	      equiv.insert(make_pair( n , t ));
	    }
	  else
	    {
	      i->second.push_back(j);
	      files[f].fields[j]->equiv = true;
	      foundEquiv = true;
	    }
	}
      if ( foundEquiv )
	{
	  PP->printLOG(" : ");
	  map<string, vector<int> >::iterator i = equiv.begin();
	  while ( i != equiv.end() )
	    {
	      if ( i->second.size() > 1 )
		{
		  map<int,int> tmap;

		  PP->printLOG(" "+files[f].fields[i->second[0]]->name + "(");
		  for (int k = 0 ; k < i->second.size(); k++)
		    {
		      if ( k>0 )
			{
			  if ( k==1 )
			    PP->printLOG("<-");
			  else
			    PP->printLOG(",");
			  tmap.insert(make_pair( i->second[k] , i->second[0] ) );
			}
		      PP->printLOG(int2str( i->second[k]+1 ));
		    }
		  PP->printLOG(")");

		  equivMap.insert( make_pair( files[f].fields[i->second[0]]->name , tmap ));
		}

	      ++i;
	    }
	}

      PP->printLOG("\n");


      ////////////////////////////////////
      // Read the raw data

      ifstream ID1( files[f].filename.c_str() , ios::in );

      if ( files[f].hasHeader )
	{
	  vector<string> header = tokenizeLine( ID1 );
	}

      while ( !ID1.eof() )
	{

	  vector<string> tokens = tokenizeLine( ID1 );

	  if ( tokens.size() == 0 )
	    continue;

	  // Insert SET:X=Y values here
	  map<string,string>::iterator i = files[f].injections.begin();
	  while ( i != files[f].injections.end() )
	    {
	      tokens.push_back( i->second );
	      // Note -- this won't be same order -- need to check/fix this???
	      ++i;

	    }

	  if ( tokens.size() != file->uniqFieldCount )
	    {
	      PP->printLOG("\n\nIn [ " + file->filename
		       + " ] encountered a row with the wrong number of fields\n");
	      PP->printLOG("Found " + int2str( tokens.size() )
		       + " fields but expecting " + int2str( file->fields.size() ) + "\n");
	      int mx = tokens.size() > file->uniqFieldCount ? tokens.size() : file->uniqFieldCount ;

	      for (int j = 0 ; j < mx ; j++)
		{
		  if ( j < file->uniqFieldCount )
		    PP->printLOG( "   " + file->fields[j]->name + ":\t" );
		  else
		    PP->printLOG( "   {?}:\t" );

		  if ( j < tokens.size() )
		    PP->printLOG( "   " + tokens[j] + "\n" );
		  else
		    PP->printLOG( "   {?}\n" );

		}


	      error("Problem with [ " + file->filename + " ]\n" );
	    }

	  IDGroup * g = new IDGroup;

	  // Track which file this group of IDs came form
	  g->file = file;

	  // Track what the original eq-value is for this line
	  map<IDField*,string> originalEquivalence;

	  for ( int j = 0 ; j < files[f].fields.size(); j++ )
	    {

	      // Is this a missing value?

	      if ( files[f].missingValues.find( tokens[j] ) != files[f].missingValues.end() )
		continue;

	      IDField * myField = files[f].fields[j];

	      if ( myField->null )
		continue;

	      string val = tokens[j];

	      // Handle equivalence specifications (aliases)

	      bool needToStore = true;

	      if ( val.find( files[f].alias_delimit ) != string::npos )
		{
		  // This is now an equiv. field
		  myField->equiv = true;
		  NList nl(0);
		  nl.setDelimiter( files[f].alias_delimit );
		  nl.setRangeChar(" "); // allow hyphens
		  vector<string> atoken = nl.deparseStringList( val );
		  string first = "";
		  for ( int i=0; i<atoken.size(); i++)
		    {
		      if ( files[f].missingValues.find( atoken[i] ) == files[f].missingValues.end() )
			{
			  setAlias( myField, atoken[i] , f, originalEquivalence );
			  if ( first=="" ) first = atoken[i];
			}
		    }
		  // Store the first value as the main value
		  if ( first != "" )
		    tokens[j] = first;
		  else
		    continue;  // if all missing
		}
	      else if ( myField->equiv )
		needToStore = setAlias( myField, tokens[j] , f , originalEquivalence );

	      if ( needToStore )
		{
		  IDValue * v = new IDValue;
		  v->field = myField;
		  v->value = tokens[j];
		  g->values.push_back(v);

		  // for pretty-printing
		  if ( v->value.size() + 3 > myField->width )
		    myField->width = v->value.size() + 3 ;
		}
	    }

	  idgroup.push_back(g);
	}

      ID1.close();
    }


  /////////////////////////////////
  // Done reading in the raw data


//   cout << "DISPLAY LEVEL 0\n";
//   for ( int g = 0 ; g < idgroup.size(); g++ )
//     idgroup[g]->display();
//   cout << "-------------------------------------------\n";


  //////////////////////////////////////////////////////////////
  // 1. Swap in preferred values for any equiv fields

  for ( int g = 0 ; g < idgroup.size(); g++ )
    {
      IDGroup * group = idgroup[g];
      for (int j = 0 ; j < group->values.size(); j++)
	if ( group->values[j]->field->equiv )
	  group->values[j]->updateAlias();
    }



  //////////////////////////////////////////////////////////////
  // 1b. Compile joint fields into joint values

  for ( int g = 0 ; g < idgroup.size(); g++ )
    {
      IDGroup * group = idgroup[g];
      setJointValues( idgroup[g] );
    }



  ////////////////////////////////
  // 2. Create the idmap


  for ( int g = 0 ; g < idgroup.size(); g++ )
    {
      IDGroup * group = idgroup[g];
      for (int j = 0 ; j < group->values.size(); j++)
	{

	  IDValue & v = *(group->values[j]);

	  map<IDValue,set<IDGroup*> >::iterator i = idmap.find( *(group->values[j]) );
	  if ( i == idmap.end() )
	    {
	      set<IDGroup*> t;
	      t.insert(group);
	      idmap.insert(make_pair( *(group->values[j]) , t) );
	    }
	  else
	    {
	      i->second.insert( group );
	    }
	}
    }


  ///////////////////////////////////////////////
  // 2.5 Simple match in any file?

  if ( par::idhelp_dump_from_dict )
    {
      idDump();
      return;
    }


  ///////////////////////////////////////////////
  // 3. Attempt to resolve into a single table


  bool okay = true;
  map<string,string> problem ;

  while (1)
    {

      bool allDone = true;

      for ( int g = 0 ; g < idgroup.size(); g++ )
	{

	  IDGroup * group = idgroup[g];

	  // Has this group already been assigned to a person?

	  if ( group->resolved )
	    continue;

	  // Find all other groups (resolved or otherwise) that this group
	  // matches with, but ignoring attributes

	  set<IDGroup*> matches;

	  for (int j = 0 ; j < group->values.size(); j++)
	    {
	      // Skip matching on attributes
	      if ( group->values[j]->field->attribute )
		continue;

	      map<IDValue, set<IDGroup*> >::iterator i = idmap.find( *(group->values[j]) );
	      if ( i != idmap.end() )
		{
		  set<IDGroup*>::iterator i2 = i->second.begin();
		  while ( i2 != i->second.end() )
		    {
		      matches.insert( *i2 );
		      ++i2;
		    }
		}
	    }


	  //////////////////////////////
	  // Merge into the key group

	  // Make a set of the key groups IDValues

	  map<IDField*,IDValue*> keyValues;
	  for ( int j = 0; j < group->values.size(); j++ )
	    {
	      keyValues.insert( make_pair( group->values[j]->field, group->values[j] ));
	    }

	  set<IDGroup*>::iterator i0 = matches.begin();
	  while ( i0 != matches.end() )
	    {

	      if ( *i0 == group )
		{
		  ++i0;
		  continue;
		}


	      // Step through all the values in this matching group

	      for ( int k = 0; k < (*i0)->values.size(); k++)
		{
		  // Does the key have this field?

		  IDField * f = (*i0)->values[k]->field;

		  map<IDField*,IDValue*>::iterator i = keyValues.find( f );

		  if ( i == keyValues.end() )
		    {
		      // Insert this key value into place
		      IDValue * t = new IDValue;
		      t->field = f;
		      t->value = (*i0)->values[k]->value;
		      t->jointValue = (*i0)->values[k]->jointValue;
		      group->values.push_back(t);
		      allDone = false;

		      // Keep track of what has been added, so we don't add twice
		      keyValues.insert( make_pair( f , (*i0)->values[k] ));
		    }
		  else
		    {

		      // Something already exists -- check it is not inconsistent
		      // as if is an ID, i.e. attributes are allowed to not be
		      // unique by definition

		      if (  *((*i0)->values[k]) != *(i->second) )
			{

			  okay = false;

			  string title = "Two unique entries [ " + (*i0)->values[k]->field->name + " = ";

			  if ( (*i0)->values[k]->field->joint )
			    {
			      if ( (*i0)->values[k]->jointValue <=  i->second->jointValue )
				title += (*i0)->values[k]->jointValue + " and " + i->second->jointValue + " (joint)";
			      else
				title += i->second->jointValue + " and " + (*i0)->values[k]->jointValue + " (joint)";
			    }
			  else
			    {
			      if ( (*i0)->values[k]->value <=  i->second->value )
				title += (*i0)->values[k]->value + " and " + i->second->value;
			      else
				title += i->second->value + " and " + (*i0)->values[k]->value;
			    }

			  title += " ] that match elsewhere";

			  if ( problem.find( title ) == problem.end() )
			    {
			      string p = "\n a) ";
			      for (int z=0; z<group->values.size(); z++)
				p += group->values[z]->field->name + "="
				  + group->values[z]->value + " ";
			      p += "\n";

			      p += " b) ";
			      for (int z=0; z<(*i0)->values.size(); z++)
				p += (*i0)->values[z]->field->name + "="
				  + (*i0)->values[z]->value + " ";
			      p += "\n\n";

			      problem.insert(make_pair(title,p));
			    }

			}
		    }
		}


	      // We are now done with this IDGroup

	      (*i0)->resolved = true;

	      ++i0;
	    }
	}

      if ( allDone )
	break;

    }

  if ( ! okay )
    {
      PP->printLOG("\n\n*** Problems were detected in the ID lists:\n\n");
      map<string,string>::iterator p = problem.begin();
      while ( p != problem.end() )
	{
	  PP->printLOG( p->first + p->second );
	  ++p;
	}
      error("You need to fix the above problems");
    }




  //////////////////////////////////////////////////////////////
  // Update the IDMAP, now only for resolved values
  //  value -> idgroup

  idmap.clear();

  for ( int g = 0 ; g < idgroup.size(); g++ )
    {
      if ( idgroup[g]->resolved )
	continue;

      IDGroup * group = idgroup[g];
      for (int j = 0 ; j < group->values.size(); j++)
	{

	  IDValue & v = *(group->values[j]);
	  map<IDValue,set<IDGroup*> >::iterator i = idmap.find( *(group->values[j]) );
	  if ( i == idmap.end() )
	    {
	      set<IDGroup*> t;
	      t.insert(group);
	      idmap.insert(make_pair( *(group->values[j]) , t) );
	    }
	  else
	    {
	      i->second.insert( group );
	    }
	}
    }


//   cout << "DISPLAY LEVEL 1\n";
//   for ( int g = 0 ; g < idgroup.size(); g++ )
//     idgroup[g]->display();
//   cout << "-------------------------------------------\n";


  //////////////////////////////////////////////////////////////
  // 4. Figure out actual transformation required, and perform

  // Rules:  ID cannot contain whitespace or commas
  // Can contain "-", "_", "=", ".", etc.
  // IDs are case-sensitive
  // Cannot contain "+" or ","
  // Values if "." are taken to mean not known, n/a

  // Functions: dump all fields on this person, or group
  if ( par::idhelp_list_aliases )

  // Lookup a single person       --id-lookup ID=27364883-1
  //  or match on an attribute    --id-lookup SITE=Boston

  // Dump the entire table        (DEFAULT)
  //     or subset of cols        --id-table ID,CLIN_ID,BSP_ID

  // Take an existing file and replace a field

  //                              --id-replace [header|noheader,field=N,skip|miss|warn|list] mydata.txt ID1 ID2

  // Default behavior is to put 'missing' as the field
  // default = to autodetect a header field
  //           skip  (do not print these lines)
  //           miss  (print a missing ID code)
  //           warn  (do not allow if >1 missing)
  //           list  (print only these lines, ID in file not in DB)


  // Take an existing file and replace a field  {file} {field} {old ID} {new ID}
  //                              --id-replace mydata.txt ID1 CLIN_ID

  // par::idhelp_command =
  // dump_table,  dump_subtable, replace, lookup


    if ( par::idhelp_list_aliases )
      {
	idListAlias();
	return;
      }



  ////////////////////////////////////////////////////////////////////////
  // Replace mode

  if ( par::idhelp_replace )
    {
      idReplace();
      return;
    }


  ///////////////////////////////////////////////////////////////
  // Line up 1 or more files based on the first file

  if ( par::idhelp_match )
    {
      idMatch();
      return;
    }




  /////////////////////////////////////////////////////////////////////////
  // Lookup functions


  set<string>  subsetFields;
  map<IDField*, set<IDValue> > lookupValues;


  if ( par::idhelp_subset )
    {
      NList tlist(0);
      vector<string> ids = tlist.deparseStringList( par::idhelp_subset_string );
      for (int i=0; i<ids.size(); i++)
	{
	  if ( fieldMap.find(ids[i]) == fieldMap.end() )
	    error("Cannot find field " + ids[i] );
	  subsetFields.insert(ids[i]);
	}
    }


  if ( par::idhelp_lookup )
    {

      // Input should be in form of a comma delimited list ID=value,
      // or, attrib=value.   More than on attrib is allowed, as a comma-
      // delimited list. Only one ID can be specified (although it may be
      // a joint one)

      // If same attrib C=2,C=1,D=1
      // e.g. ( C==2 OR C==1 ) AND (D==1)

      // Use "+" to identify a joint ID set

      // FID+IID=6363723+8383293

      lookupValues = parseQuery( par::idhelp_lookup_string );

      PP->printLOG("Looking up items matching: ");

      // These will be sorted in field name order, so we
      // can easily figure out OR versus AND conditions

      map<IDField*,set<IDValue> >::iterator i = lookupValues.begin();
      while ( i != lookupValues.end() )
	{
	  set<IDValue>::iterator j = i->second.begin();
	  PP->printLOG( "\n  " + i->first->name + " = " );
	  while ( j != i->second.end() )
	    {
	      PP->printLOG( j->value + " " );
	      ++j;
	    }
	  if ( i->first->attribute ) PP->printLOG(" (attribute)");
	  else PP->printLOG(" (id)");
	  ++i;
	}
      PP->printLOG("\n");

    }



  //////////////////////////////////////////////////////////
  // Main output routine

  PP->printLOG("Writing output to [ "
	       + par::output_file_name + ".id ]\n");
  ofstream OFILE( (par::output_file_name+".id").c_str() , ios::out );

  // Header row
  set<IDField>::iterator f = fields.begin();
  while ( f != fields.end() )
    {
      if ( f->null )
	{
	  ++f;
	  continue;
	}

      if ( (! par::idhelp_subset )
	   || subsetFields.find( f->name )!=subsetFields.end() )
	OFILE << setw(f->width) << f->name << par::idhelp_output_delimit << " ";
      ++f;
    }
  OFILE << "\n";

  // Keep track of how many records we retrieve
  int numMatched = 0;

  for ( int g = 0 ; g < idgroup.size(); g++ )
    {

      IDGroup * group = idgroup[g];

      // Has this group already been assigned to a person?
      if ( group->resolved )
	continue;


      // Make a set of the key groups IDValues
      // If doing a lookup, we might also need to include
      // all fields here

      map<IDField*,IDValue*> keyValues;

      for ( int j = 0; j < group->values.size(); j++ )
	{
	  if ( par::idhelp_lookup ||
	       ( ! par::idhelp_subset ) ||
	       subsetFields.find( group->values[j]->field->name )!=subsetFields.end() )

	    keyValues.insert( make_pair( group->values[j]->field, group->values[j] ));
	}



      ///////////////////////////////////////////////////
      // If we are filtering, does this person match?

      if ( par::idhelp_lookup )
	{
	  if ( ! matchIndividual( group, lookupValues ) )
	    continue;
	}

      ++numMatched;

      ///////////////////////////////////////////
      // Print row, in same order for all fields

      set<IDField>::iterator f = fields.begin();
      while ( f != fields.end() )
	{

	  if ( f->null )
	    {
	      ++f;
	      continue;
	    }

	  if ( (! par::idhelp_subset)
	       || subsetFields.find( f->name ) != subsetFields.end() )
	    {
	      map<IDField*,IDValue*>::iterator k = keyValues.find( (IDField*)&(*f) );
	      if ( k == keyValues.end() )
		OFILE << setw( f->width ) << "." << par::idhelp_output_delimit << " ";
	      else
		{
		  OFILE << setw( f->width ) << k->second->value << par::idhelp_output_delimit << " ";
		}
	    }
	  ++f;
	}
      OFILE << "\n";
    }
  OFILE.close();

  PP->printLOG( int2str( numMatched ) + " unique records retrieved\n");
}





void IDHelper::idListAlias()
{

  PP->printLOG("Listing ID equivalents/aliases to [ " + par::output_file_name + ".id.eq ]\n");
  ofstream O1( (par::output_file_name + ".id.eq").c_str() , ios::out );

  O1 << setw(20) << "FIELD" << " "
     << setw(20) << "PREF" << " "
     << setw(20) << "EQUIV" << "\n";

  map<string,IDField*>::iterator i1 = fieldMap.begin();
  while ( i1 != fieldMap.end() )
    {
      if ( i1->second->equiv )
	{
	  IDField * f = i1->second;
	  map<string,string>::iterator j = f->eqid.begin();
	  while ( j != f->eqid.end() )
	    {
	      O1 << setw(20) << i1->first << " "
		 << setw(20) << j->second << " "
		 << setw(20) << j->first << "\n";
	      ++j;
	    }
	}
      ++i1;
    }
  O1.close();

  return;
}


void IDHelper::idReplace()
{

  NList nl(0);

  vector<string> tok = nl.deparseStringList( par::idhelp_replace_string );

  if ( tok.size() != 3 )
    error("Problem with --id-replace string format\n");

  checkFileExists( tok[0] );


  // This/these are the fields are
  vector<int> rep_field;

  NList tlist(0);
  string t = searchAndReplace( tok[1] , "+" , "," );
  vector<string> targetFields = tlist.deparseStringList( t );
  for ( int i=0; i<targetFields.size(); i++)
    if ( fieldMap.find( targetFields[i] ) == fieldMap.end() )
      error("Cannot find target field [ " + targetFields[i] + " ]");

  if ( fieldMap.find( tok[2] ) == fieldMap.end() )
    error("Cannot find replacement field [ " + tok[2] + " ]" );

  string fname = tok[2];
  IDField * f = fieldMap.find( tok[2] )->second;
  if ( f->joint )
    fname = f->jointName;

  PP->printLOG("Replacing " + tok[1] + " with "
	       + fname + " from [ " + tok[0] + " ]\n");
  PP->printLOG("Writing new file to [ " + par::output_file_name + ".rep ]\n");

  OptionSet * id_opt = par::opt.getOptions("IDHELP");

  bool skipMode = id_opt->isSet("skip");
  bool missMode = id_opt->isSet("miss");
  bool warnMode = id_opt->isSet("warn");
  bool listMode = id_opt->isSet("list");

  int c = 0;

  if ( skipMode )
    {
      PP->printLOG("Set to skip unmatched observations\n");
      ++c;
    }

  if ( missMode )
    {
      PP->printLOG("Set to set unmatched observations to NA\n");
      ++c;
    }

  if ( warnMode )
    {
      PP->printLOG("Set to give error for first unmatched observation\n");
      ++c;
    }

  if ( listMode )
    {
      PP->printLOG("Set to list only IDs in file but not in database\n");
      ++c;
    }

  if ( c == 0 )
    PP->printLOG("Set to keep original value for unmatched observations\n");

  if ( c>1 )
    error("Can only specify one of [miss|warn|skip|list] options in --id-replace");


  // Do we have a header row?
  bool header = false;
  if ( id_opt->isSet("header") )
    header = true;


  // If not, we need a number specified
  if ( ! header )
    {
      string field_str;

      if ( id_opt->isSet("field") )
	field_str = id_opt->getValue("field");
      else
	error("Need to specify field={N,N} if no header");

      // Convert field str to vector of fields
      NList tlist(0);
      string t = searchAndReplace( field_str , "+" , "," );
      vector<string> ids = tlist.deparseStringList( t );

      if ( ids.size() != targetFields.size() )
	error("Must specify the same number of fields/cols");
      rep_field.resize( ids.size() , -1 );
      for (int i=0; i<ids.size(); i++)
	{
	  if ( ! from_string<int>( rep_field[i], ids[i] , std::dec ) )
	    error("Problem with field specified in --id-replace options");
	  // Make 0-based
	  --rep_field[i];
	}
    }
  else
    {

      // Lookup in header row

      ifstream IN1( tok[0].c_str() , ios::in );
      vector<string> tokens = tokenizeLine( IN1 );
      rep_field.resize( targetFields.size() , -1 );

      for (int f = 0 ; f < targetFields.size(); f++)
	for (int i = 0 ; i < tokens.size(); i++)
	  {
	    if ( tokens[i] == targetFields[f] )
	      rep_field[f] = i;
	  }
      for (int i=0; i<rep_field.size(); i++)
	if ( rep_field[i] == -1 )
	  error("Could not find field in header for --id-replace");
      IN1.close();
      IN1.clear();
    }

  int insertField = rep_field[0];


  /////////////////////////////////////////////
  // Load input file

  ifstream IN1( tok[0].c_str() , ios::in );
  ofstream OUT1( (par::output_file_name+".rep").c_str() , ios::out );
  int notFound = 0;
  bool readHeader = false;
  if ( ! header )
    readHeader = true;

  int maxfield = 0;
  for (int i=0; i<rep_field.size(); i++)
    if ( rep_field[i] > maxfield )
      maxfield = rep_field[i];

  while( ! IN1.eof() )
    {
      vector<string> tokens = tokenizeLine( IN1 );
      if ( tokens.size() == 0 )
	continue;
      if ( tokens.size() <= maxfield )
	error("Not enough columns here");

      bool changed = false;

      // Deal with header row?
      if ( ! readHeader )
	{
	  for (int i=0; i<rep_field.size(); i++)
	    {
	      if ( i == 0 )
		tokens[ rep_field[i] ] = fname;
	      else
		tokens[ rep_field[i] ] = "";
	    }

	  readHeader = true;
	  changed = true;
	}
      else
	{

	  // For actual data, read all fields, and find the
	  // individual that matches

	  set<IDValue> myTemplate;
	  for (int i=0; i<rep_field.size(); i++)
	    {
	      IDValue findField;
	      findField.field = fieldMap.find( targetFields[i] )->second;
	      findField.value = tokens[ rep_field[i] ];
	      if ( findField.field->equiv )
		findField.updateAlias();
	      myTemplate.insert( findField );
	    }

	  // Connect up any joint fields
	  setJointValues( myTemplate );

	  // Find the match
	  IDGroup * thisGroup = findUniqueIndividual( myTemplate );

	  if ( thisGroup == NULL )
	    {
	      ++notFound;
	    }
	  else
	    {
	      // Find each old field
	      for (int f=0; f<targetFields.size(); f++)
		{
		  for (int i = 0 ; i < thisGroup->values.size(); i++)
		    if ( thisGroup->values[i]->field->name == tok[2] )
		      {
			if ( f>0 )
			  tokens[ rep_field[f] ] = ".";
			else
			  {
			    // Replace this item: handle if the replacing
			    // field is itself a joint one

			    if ( thisGroup->values[i]->field->joint )
			      {
				tokens[ rep_field[f] ] =
				  searchAndReplace( thisGroup->values[i]->jointValue,"+"," " );
			      }
			    else
			      tokens[ rep_field[f] ] = thisGroup->values[i]->value;
			  }

			changed = true;
		      }
		}
	    }

	  // Done reading, processng this line
	}


      ///////////////////////////////
      // Output this line

      if ( ! changed )
	{
	  if ( listMode )
	    {
	      for (int i=0; i<rep_field.size(); i++)
		OUT1 << tokens[ rep_field[i] ] << "\t";
	      OUT1 << "\n";
	      continue;
	    }

	  if ( skipMode )
	    continue;

	  if ( warnMode )
	    error("Could not find replacement for " + tokens[insertField] + "\n");

	  if ( missMode )
	    {
	      for (int f=0; f<targetFields.size(); f++)
		{
		  if ( f == 0 )
		    tokens[ rep_field[0] ] = "NA";
		  else
		    tokens[ rep_field[f] ] = ".";
		}
	    }
	}

      if ( listMode )
	continue;


      ////////////////////////////////////////////////////////
      // Print out line, which may or may not be modified

      for (int i = 0 ; i < tokens.size() ; i++ )
	OUT1 << tokens[i] << par::idhelp_output_delimit << " ";
      OUT1 << "\n";

    }

  if ( notFound > 0 )
    PP->printLOG("Could not find matches for " + int2str( notFound ) + " lines\n");

  OUT1.close();
  IN1.close();

  return;
}



void IDHelper::idMatch()
{

  // e.g. --id-match myfile.fam FID+IID 1+2 file1.txt CLIN_ID,1 file2.txt ID

  // in form {file} {id,{col}}
  // where joint IDs are ID1+ID2, or with fields: ID1+ID2,5+7

  // If joint IDs, then all must be specified.
  // Cannot specify more than 1 non-joint ID though
  // Can be different IDs in different files

  if ( par::idhelp_match_string.size() < 4 )
    error("Must specify more than 1 file to match");

  // Assemble all data here:
  vector<map<IDGroup*, vector<string> > > table;
  vector<int> tableSize;

  // And keep track of which files which individuals are in
  map<IDGroup*,set<int> > foundIn;

  // Keep track of order from first file
  map<int,IDGroup*> fileOrder;
  set<IDGroup*> seenBefore;

  // do we see at least 1 header
  bool atleastOneHeader = false;
  vector< vector<string> > headers;


  // Create a single field on the fly, for use in simple match mode
  IDField * nullField = new IDField;
  nullField->name = "tmp1";

  //////////////////
  // Read each file

  for (int s=0; s< par::idhelp_match_string.size(); s+=2)
    {

      PP->printLOG("Matching [ " + par::idhelp_match_string[s]
		   + " ] on " + par::idhelp_match_string[s+1] + "\n");

      int t = (int)s/2;

      if ( par::idhelp_no_dict )
	{
	  IDFile f;
	  f.filename = "F" + int2str(t);
	  files.push_back(f);
	}

      map<IDGroup*, vector<string> > inserts;

      // Each element should be in form: filename ID   filename2 ID+ID   filename3 ID ...

      checkFileExists( par::idhelp_match_string[s] );
      ifstream I1( (par::idhelp_match_string[s]).c_str() , ios::in );

      string id = par::idhelp_match_string[s+1];

      // A         Implies header row
      // A,2       Implies no header row
      // A+B       Implies header row
      // A+B,2+3   Implies no header row

      // If in "quick-match" mode, then we assume that the ID column is always the
      // same, whether or not it is explicitly named differently here


      bool jointQuery = id.find("+") != string::npos;

      // Fields to match on
      vector<int> fieldCodes;
      vector<string> fieldNames;
      int maxF = -1;

      if ( id.find(",") != string::npos )
	{

	  NList nl(0);
	  nl.setDelimiter("+");
	  nl.setRangeChar(" ");
	  vector<string> fstr = nl.deparseStringList(  id.substr( id.find(",")+1 ) );
	  for (int i=0; i<fstr.size(); i++)
	    {
	      int myf;
	      if ( ! from_string<int>( myf, id.substr( id.find(",")+1 ) , std::dec ) )
		error("Trouble converting to a field number");
	      // Make zero-based
	      --myf;
	      if ( myf < 0 )
		error("Invalid value for field # specified");
	      if ( myf > maxF ) maxF = myf;
	      fieldCodes.push_back( myf );
	    }

	  string tmp = id.substr( 0, id.find(",") );

	  NList nl2(0);
	  nl2.setDelimiter("+");
	  nl2.setRangeChar(" ");
	  fieldNames = nl2.deparseStringList( tmp );

	  if ( fieldNames.size() != fieldCodes.size() )
	    error("Problem with joint ID specification in: " + id );

	  if ( ! par::idhelp_no_dict )
	    for (int f = 0; fieldNames.size(); f++)
	      if ( fieldMap.find( fieldNames[f] ) == fieldMap.end() )
		error("Field " + fieldNames[f] + " does not exist in the database");

	  // Read rest of line; insert dummy headers
	  vector<string> h = tokenizeLine(I1);
	  I1.close();
	  I1.clear();
	  I1.open( (par::idhelp_match_string[s]).c_str() , ios::in );
	  vector<string> header;
	  for (int k=0; k<h.size(); k++)
	    {
	      header.push_back( "F"+int2str(t)+"_"+int2str(k) );
	    }
	  headers.push_back(header);

	}
      else
	{

	  // Read first line, which should be a header
	  atleastOneHeader = true;
	  vector<string> header = tokenizeLine(I1);

	  NList nl2(0);
	  nl2.setDelimiter("+");
	  nl2.setRangeChar(" ");
	  fieldNames = nl2.deparseStringList( id );

	  // Find each field
	  for (int f = 0; f<fieldNames.size(); f++)
	    {
	      if ( ! par::idhelp_no_dict )
		if ( fieldMap.find( fieldNames[f] ) == fieldMap.end() )
		  error("Field " + fieldNames[f] + " does not exist in the database");

	      bool foundField = false;
	      for (int i=0; i<header.size(); i++)
		{
		  if ( header[i] == fieldNames[f] )
		    {
		      fieldCodes.push_back( i );
		      if ( i > maxF ) maxF = i;
		      foundField = true;
		      break;
		    }
		}

	      if ( ! foundField )
		error("Could not find field " + fieldNames[f]
		      + " in [ " + par::idhelp_match_string[s] + " ]");
	    }

	  headers.push_back(header);
	}

      vector<IDField *> thisField;
      if ( ! par::idhelp_no_dict )
	{
	  for (int f = 0 ; f < fieldNames.size(); f++)
	    thisField.push_back( fieldMap.find( fieldNames[f] )->second );
	}
      else
	{
	  thisField.push_back( nullField );
	}


      // Keep track of the # of columns per file, as a check
      // for rectangular files

      tableSize.push_back(-1);



      ////////////////////////////////////////
      // Read each row of the data files:

      int notFound = 0;
      string missingList = "";

      while ( ! I1.eof() )
	{

	  vector<string> tok = tokenizeLine(I1);
	  if ( tok.size() == 0 )
	    continue;

	  if ( tableSize[t] == -1 )
	    tableSize[t] = tok.size();
	  else
	    if ( tok.size() != tableSize[t] )
	      {
		string msg = "Problem with non-rectangular file [ "
		  + par::idhelp_match_string[s] + " ]\n";
		msg += "Execting " + int2str( tableSize[t] )
		  + " fields but found " + int2str(tok.size()) + "\n";
		for (int k=0; k<tok.size(); k++)
		  msg += tok[k] + " ";
		error( msg );
	      }

	  // Warn when line does not contain enough columns
	  if ( maxF >= tok.size() )
	    error("Line does not contain enough columns");

	  set<IDValue> myTemplate;
	  for (int f=0; f<fieldCodes.size(); f++)
	    {
	      IDValue findField;
	      findField.field = thisField[f];
	      findField.value = tok[ fieldCodes[f] ];
	      if ( findField.field->equiv )
		findField.updateAlias();
	      myTemplate.insert( findField );
	    }

	  // Connect up any joint fields
	  setJointValues( myTemplate );


	  // If running without a dictionary, now add this person in
	  if ( par::idhelp_no_dict )
	    {
	      IDGroup * thisGroup = findUniqueIndividual( myTemplate );
	      if ( thisGroup == NULL )
		{
		  IDGroup * g = new IDGroup;
		  g->resolved = false;
		  g->file = &files[t];
		  set<IDValue>::iterator k = myTemplate.begin();
		  while ( k != myTemplate.end() )
		    {
		      IDValue * nv = new IDValue;
		      *nv = *k;
		      g->values.push_back( nv );
		      idgroup.push_back(g);
		      ++k;
		    }

		  // Add to the ID Map
		  set<IDGroup*> t;
		  t.insert( idgroup[ idgroup.size()-1] );
		  for (int v=0; v<g->values.size(); v++)
		    idmap.insert(make_pair( *(g->values[v]) ,t ));

		}
	    }


	  // Find the match
	  IDGroup * thisGroup = findUniqueIndividual( myTemplate );

	  if ( thisGroup != NULL )
	    {
	      // Add this row to the collection
	      inserts.insert(make_pair( thisGroup, tok ) );

	      // Keep track of order in which we came across this person for the
	      // first time

	      if ( seenBefore.find( thisGroup ) == seenBefore.end() )
		{
		  int sz = fileOrder.size();
		  fileOrder.insert(make_pair(sz,thisGroup));
		  seenBefore.insert( thisGroup );
		}

	      // Keep track that this person was in this file
	      map<IDGroup*,set<int> >::iterator i1 = foundIn.find( thisGroup );
	      if ( i1 == foundIn.end() )
		{
		  set<int> t;
		  t.insert(int(s/2));
		  foundIn.insert(make_pair( thisGroup, t ) );
		}
	      else
		i1->second.insert(int(s/2));
	    }
	  else
	    {
	      set<IDValue>::iterator i = myTemplate.begin();
	      while ( i != myTemplate.end() )
		{
		  missingList += i->field->name + " = " + i->value + "\t";
		  ++i;
		}
	      missingList += "\n";
	      ++notFound;
	    }

	  // Get next line from this file
	}

      // Add this file to big collection
      table.push_back( inserts );


      if ( notFound > 0 )
	{
	  PP->printLOG("Could not find "
		       + int2str( notFound )
		       + " individuals from [ "
		       + par::idhelp_match_string[s]
		       + " ] in database\n");
	  PP->printLOG("Writing this list to [ " + par::output_file_name + ".noid ]\n");
	  ofstream O1( ( par::output_file_name+".noid").c_str() , ios::out );
	  O1 << "FIELD = VALUE\n";
	  O1 << missingList ;
	  O1.close();
	}

      I1.close();

      // Get next file
    }


  // Now output all files tied together
  // Output, in order of the order in which we encountered
  // each unique individual

  // Only output complete rows?
  OptionSet * id_opt = par::opt.getOptions("IDHELP");

  bool complete = false;
  if ( id_opt->isSet("complete") )
    complete = true;
  if ( id_opt->isSet("noheader") )
    atleastOneHeader = false;

  PP->printLOG("Writing output file to [ " + par::output_file_name + ".matched ]\n");

  ofstream O1( (par::output_file_name+".matched").c_str() , ios::out );

  // What about header row?
  // skip for now...

  if ( atleastOneHeader )
    {
      for (int k=0; k<headers.size(); k++)
	for (int j=0; j<headers[k].size(); j++)
	  O1 << setw(12) << headers[k][j] << " ";
      O1 << "\n";
    }


  int tfiles = table.size();
  map<int,IDGroup*>::iterator j = fileOrder.begin();
  while ( j != fileOrder.end() )
    {

      map<IDGroup*,set<int> >::iterator i = foundIn.find( (IDGroup* const)j->second );

      int f = i->second.size();
      if ( complete && f != tfiles )
	{
	  ++j;
	  continue;
	}

      for (int t = 0 ; t < tfiles; t++)
	{
	  // Does this individual exist for this file?
	  map<IDGroup*,vector<string> > & thisTable = table[t];

	  if ( thisTable.find( i->first ) == thisTable.end() )
	    {
	      for (int k=0; k<tableSize[t]; k++)
		O1 << setw(12) << "NA" << " ";
	    }
	  else
	    {
	      vector<string> & thisLine = thisTable.find( i->first )->second;
	      for (int k = 0 ; k < thisLine.size() ; k++ )
		O1 << setw(12) << thisLine[k] << " ";
	    }
	}
      O1 << "\n";

      // Next individual
      ++j;
    }

  O1.close();

  return;
} // end of id-match routine



void IDHelper::idDump()
{

  // Note: use of parseQuery modifies the data structure, setting
  // fields to non-joint potentially, and so subsequent attempts to
  // line things up will not work; therefore, for now we stop here.

  PP->printLOG("\nReporting rows that match [ " + par::idhelp_dump_from_dict_cmd + " ] \n\n");

  map<IDField*,set<IDValue> > myTemplate = parseQuery( par::idhelp_dump_from_dict_cmd );

  for ( int g = 0 ; g < idgroup.size(); g++ )
    {
      IDGroup * group = idgroup[g];
      if ( matchIndividual( group , myTemplate ) )
	{
	  for (int j = 0 ; j < group->values.size(); j++ )
	    {
	      PP->printLOG( group->file->filename + " : " );
	      PP->printLOG( group->values[j]->field->name + " = " );
	      PP->printLOG( group->values[j]->value + "\n" );
	    }
	  PP->printLOG("\n");
	}
    }

  PP->printLOG("---------------------------------\n");

  // Important: joint status of fields will have changed; in any
  // case, let's stop here
  return;
}
