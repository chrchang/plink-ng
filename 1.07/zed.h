#ifndef __ZED_H__
#define __ZED_H__

#include <cstring>
#include <vector>
#include <fstream>


#ifdef WITH_ZLIB
#include "zfstream.h"
#endif

using namespace std;

// A lightweight wrapper around zfstream wrapper to zlib

const int MAX_LINE_LENGTH = 1000000;

class ZInput
{
  string filename;
  bool compressed;  
  char buf[MAX_LINE_LENGTH];

#ifdef WITH_ZLIB
  gzifstream zinf;
#else
  ifstream zinf;
#endif

  ifstream inf;

 public:
  ZInput(string, bool);
  ZInput();
  void open(string, bool);
  char readChar();
  string readLine();
  vector<string> tokenizeLine();
  void close();
  bool endOfFile();
  void unbuffered();
};

class ZOutput
{

#ifdef WITH_ZLIB
  gzofstream zoutf;
#else
  ofstream zoutf;
#endif

  ofstream outf;

  string filename;
  bool compressed;
  char buf[MAX_LINE_LENGTH];
 public:
  ZOutput(string, bool);
  ZOutput();
  void open(string,bool);
  void write(string);
  ZOutput & operator<< (const string & s) { write(s); return *this; }
  void writeLine(string);
  void close();
  void unbuffered();
};

void fileCompress();
void fileUncompress();


#endif
