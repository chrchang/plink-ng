

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <vector>

#include "zed.h"
#include "helper.h"
#include "nlist.h"

extern Plink * PP;

ZInput::ZInput(string f, bool cmode)
{
  open(f,cmode);
}

ZInput::ZInput()
{
  //
}


void ZInput::open(string f, bool cmode)
{
  filename = f;
  compressed = cmode;

#ifndef WITH_ZLIB
  if ( compressed )
    {
      error("ZLIB support is not currently compiled in");
    }
#endif  

  if ( compressed )
    {
      zinf.open( filename.c_str() );
      if ( ! zinf.is_open() )
	error("Problem opening " + filename + "\n");
    }
  else
    {
      inf.open( filename.c_str() );
      if ( ! inf.is_open() )
	error("Problem opening " + filename + "\n");
    }
}

string ZInput::readLine()
{
  
  if ( compressed ) 
    {
      zinf.getline(buf,MAX_LINE_LENGTH,'\n');
    }
  else
    {
      inf.getline(buf,MAX_LINE_LENGTH,'\n');
    }
  
  return buf;

  // std::cerr << buf 
  //           << "\t(" << inf.rdbuf()->in_avail() 
  //           << " chars left in buffer)  ";
}

char ZInput::readChar()
{
  char c;
  if ( compressed ) 
    {
      zinf.get(c);
    }
  else
    {
      inf.get(c);
    }  
  return c;
}

vector<string> ZInput::tokenizeLine()
{  
  string s = readLine();
  vector<string> tok;
  string buf; 
  stringstream ss(s); 
  while (ss >> buf)
    tok.push_back(buf);
  return tok;
}

void ZInput::close()
{
  if ( compressed ) 
    inf.close();
  else
    zinf.close();
}

bool ZInput::endOfFile()
{
  // Check -- eof() doesn't work here -- look up the differences between 
  // these different file states

  if ( compressed )
    return zinf.fail() || ( ! zinf.good() ) ;
  else
    return inf.fail() || ( ! inf.good() ) ;
}

void ZInput::unbuffered()
{
  if ( compressed )
    zinf.rdbuf()->pubsetbuf(0,0);
}

void ZOutput::open(string f, bool cmode)
{
  filename = f;
  compressed = cmode;

#ifndef WITH_ZLIB
  if ( compressed )
    {
      PP->printLOG("Warning: ZLIB support not enabled, so writing uncompressed file\n");
      compressed = false;
    }
#endif  

  if ( compressed )
    {
      zoutf.open( filename.c_str() );
      if ( ! zoutf.is_open() )
	{
	  error("Problem opening " + filename );
	}
    }
  else
    {
      outf.open( filename.c_str() );
      if ( ! outf.is_open() )
	{
	  error("Problem opening " + filename );
	}
    }
}

ZOutput::ZOutput(string f, bool cmode)
{
  open(f,cmode);
}

ZOutput::ZOutput()
{
  //
}

void ZOutput::write(string s)
{
  if ( compressed )
    zoutf << s; 
  else
    outf << s;
}

void ZOutput::writeLine(string s)
{
  if ( compressed )
    zoutf << s << endl;
  else
    outf << s << endl;
}

void ZOutput::close()
{
  if ( compressed )
    zoutf.close();
  else
    outf.close();
}

void ZOutput::unbuffered()
{
  if ( compressed )
    zoutf.rdbuf()->pubsetbuf(0,0);
}

void fileCompress()
{

#ifndef WITH_ZLIB
  error("ZLIB support is not compiled in");
#endif

  PP->printLOG("Compressing [ " + par::compress_filename + " ]...\n");
  ZInput zin( par::compress_filename , false );
  ZOutput zout( par::compress_filename+".gz", true );
  while ( ! zin.endOfFile() )
    {
      string line = zin.readLine();
      if ( zin.endOfFile() )
	break;
      zout.writeLine(line);
    }
  zin.close();
  zout.close();
  PP->printLOG("Wrote compressed file to [ " + par::compress_filename + ".gz ]\n");
}


void fileUncompress()
{

#ifndef WITH_ZLIB
  error("ZLIB support is not compiled in");
#endif
      
  PP->printLOG("Uncompressing [ " + par::compress_filename + " ]...\n");

  int s = par::compress_filename.size();

  if ( s < 3 || par::compress_filename.substr(s-3,3) != ".gz" )
    error("Filename must end if .gz");

  ZInput zin( par::compress_filename , true );
  ZOutput zout( par::compress_filename.substr(0,s-3), false );
  while ( ! zin.endOfFile() )
    {
      string line = zin.readLine();
      if ( zin.endOfFile() )
	break;
      zout.writeLine(line);
    }
  zin.close();
  zout.close();
  PP->printLOG("Wrote uncompressed file to [ " 
	       + par::compress_filename.substr(0,s-3) + " ]\n");
}
