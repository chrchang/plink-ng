/*
 *  C++ Interface to Rserve
 *  Copyright (C) 2004-8 Simon Urbanek, All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; version 2.1 of the License
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Leser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  Although this code is licensed under LGPL v2.1, we strongly encourage
 *  everyone modifying this software to contribute back any improvements and
 *  bugfixes to the project for the benefit all other users. Thank you.
 *
 *  $Id$
 */

/* external defines:
   SWAPEND  - needs to be defined for platforms with inverse endianess related to Intel
   MAIN     - should be defined in just one file that will contain the fn definitions and variables
              (this is inherited from Rsrv.h and sisocks.h)
*/

#ifndef __RCONNECTION_H__
#define __RCONNECTION_H__

#if defined __GNUC__ && !defined unix && !defined Win32 && !defined WIN32
#define unix
#endif

#include <iostream>

#include <stdio.h>
#include "sisocks.h"
#include "Rsrv.h"

typedef unsigned int Rsize_t;

//=== Rconnection error codes

#define CERR_connect_failed    -1
#define CERR_handshake_failed  -2
#define CERR_invalid_id        -3
#define CERR_protocol_not_supp -4
#define CERR_not_connected     -5
#define CERR_peer_closed       -7
#define CERR_malformed_packet  -8
#define CERR_send_error        -9
#define CERR_out_of_mem       -10
#define CERR_not_supported    -11
#define CERR_io_error         -12

// this one is custom - authentication method required by
// the server is not supported in this client
#define CERR_auth_unsupported -20


#define A_required 0x001
#define A_crypt    0x002
#define A_plain    0x004

//===================================== Rmessage ---- QAP1 storage

class Rmessage {
 public:
    struct phdr head;
    char *data;
    Rsize_t len;
    int complete;
    
    // the following is avaliable only for parsed messages (max 16 pars)
    int pars;
    unsigned int *par[16];

    Rmessage();
    Rmessage(int cmd); // 0 data
    Rmessage(int cmd, const char *txt); // DT_STRING data
    Rmessage(int cmd, int i); // DT_INT data (1 entry)
    Rmessage(int cmd, const void *buf, int len, int raw_data=0); // raw data or DT_BYTESTREAM
    virtual ~Rmessage();
        
    int command() { return complete?head.cmd:-1; }
    Rsize_t length() { return complete?head.len:-1; }
    int is_complete() { return complete; }
    
    int read(int s);
    void parse();
    int send(int s);    
};

//===================================== Rexp --- basis for all SEXPs

class Rexp {
public:
    Rmessage *msg;
    unsigned int *pos;
    Rsize_t len;
    Rexp *attr;
    int type;
    /* memory manegement for data/len:
        - content is in a message and this Rexp is the master of that message:
          master=0; msg=<source message>;
        - content is in a message, but this Rexp is not the master
          master=<master Rexp>; msg=0
        - content is all self-allocated with no message associated
          master=this; msg=0 */
    char *data, *next;

protected:
    // the next two are only cached if requested, no direct access allowed
    int attribs; 
    const char **attrnames;
    
    Rexp *master; // if this is set then this Rexp allocated the memory for us, so we are not supposed to free anything; if this is set to "this" then the content is self-allocated, including any data
    int rcount;  // reference count - only for a master - it counts how many children still exist
    
public:
    Rexp(Rmessage *msg);
    Rexp(unsigned int *pos, Rmessage *msg=0);
    Rexp(int type, const char *data=0, int len=0, Rexp *attr=0);
    
    virtual ~Rexp();
    
    void set_master(Rexp *m);
    char *parse(unsigned int *pos);

    virtual Rsize_t storageSize() { return len+((len>0x7fffff)?8:4); }
    
    virtual void store(char *buf);
    Rexp *attribute(const char *name);
    const char **attributeNames();
    
    virtual Rsize_t length() { return len; }
    
    friend std::ostream& operator<< (std::ostream& os, const Rexp& exp) {
        return ((Rexp&)exp).os_print(os);
    }

    friend std::ostream& operator<< (std::ostream& os, const Rexp* exp) {
        return ((Rexp*)exp)->os_print(os);
    }
    
    virtual std::ostream& os_print(std::ostream& os) {
        return os << "Rexp[type=" << type << ",len=" << len <<"]";
    }
};

//===================================== Rint --- XT_INT/XT_ARRAY_INT

class Rinteger : public Rexp {
public:
    Rinteger(Rmessage *msg) : Rexp(msg) { fix_content(); }
    Rinteger(unsigned int *ipos, Rmessage *imsg) : Rexp(ipos, imsg) { fix_content(); }
    Rinteger(int *array, int count) : Rexp(XT_ARRAY_INT, (char*)array, count*sizeof(int)) { fix_content(); }
    
    int *intArray() { return (int*) data; }
    int intAt(int pos) { return (pos>=0 && (unsigned)pos<len/4)?((int*)data)[pos]:0; }
    virtual Rsize_t length() { return len/4; }

    virtual std::ostream& os_print (std::ostream& os) {
        return os << "Rinteger[" << (len/4) <<"]";
    }
    
private:
    void fix_content();
};

//===================================== Rdouble --- XT_DOUBLE/XT_ARRAY_DOUBLE

class Rdouble : public Rexp {
public:
    Rdouble(Rmessage *msg) : Rexp(msg) { fix_content(); }
    Rdouble(unsigned int *ipos, Rmessage *imsg) : Rexp(ipos, imsg) { fix_content(); }
    Rdouble(double *array, int count) : Rexp(XT_ARRAY_DOUBLE, (char*)array, count*sizeof(double)) { fix_content(); }
    
    double *doubleArray() { return (double*) data; }
    double doubleAt(int pos) { return (pos>=0 && (unsigned)pos<len/8)?((double*)data)[pos]:0; }
    virtual Rsize_t length() { return len/8; }

    virtual std::ostream& os_print (std::ostream& os) {
        return os << "Rdouble[" << (len/8) <<"]";
    }
    
private:
    void fix_content();
};

//===================================== Rsymbol --- XT_SYM

class Rsymbol : public Rexp {
protected:
    const char *name;
    
public:
    Rsymbol(Rmessage *msg) : Rexp(msg)
    { name=""; fix_content(); }
    
    Rsymbol(unsigned int *ipos, Rmessage *imsg) : Rexp(ipos, imsg)
    { name=""; fix_content(); }
    
    const char *symbolName() { return name; }
    
    virtual std::ostream& os_print (std::ostream& os) {
        return os << "Rsymbol[" << symbolName() <<"]";
    }
    
private:
    void fix_content();
};

//===================================== Rstrings --- XT_ARRAY_STR
// NOTE: XT_ARRAY_STR is new in 0103 and ths class is just a
//       very crude implementation. It replaces Rstring because
//       XT_STR has been deprecated.
// FIXME: it should be a subclass of Rvector!
class Rstrings : public Rexp {
    char **cont;
    unsigned int nel;
public:
   Rstrings(Rmessage *msg) : Rexp(msg) { decode(); }
   Rstrings(unsigned int *ipos, Rmessage *imsg) : Rexp(ipos, imsg) { decode(); }
    /*Rstring(const char *str) : Rexp(XT_STR, str, strlen(str)+1) {}*/
    
    char **strings() { return cont; }
    char *stringAt(unsigned int i) { return (i < 0 || i >= nel) ? 0 : cont[i]; }
    char *string() { return stringAt(0); }
    virtual Rsize_t length() { return nel; }

    unsigned int count() { return nel; }
    int indexOfString(const char *str);

    virtual std::ostream& os_print (std::ostream& os) {
        return os << "char*[" << nel <<"]\"" << string() <<"\"..";
    }
 private:
    void decode() {
      char *c = (char*) data;
      unsigned int i = 0;
      nel = 0;
      while (i < len) { if (!c[i]) nel++; i++; }
      if (nel) {
	i = 0;
	cont = (char**) malloc(sizeof(char*)*nel);
	while (i < nel) {
	  cont[i] = strdup(c);
	  while (*c) c++;
	  c++; i++;
	}	
      } else
	cont = 0;
    }
};

//===================================== Rstring --- XT_STR

class Rstring : public Rexp {
public:
    Rstring(Rmessage *msg) : Rexp(msg) {}
    Rstring(unsigned int *ipos, Rmessage *imsg) : Rexp(ipos, imsg) {}
    Rstring(const char *str) : Rexp(XT_STR, str, strlen(str)+1) {}
    
    char *string() { return (char*) data; }

    virtual std::ostream& os_print (std::ostream& os) {
        return os << "\"" << string() <<"\"";
    }
};



//===================================== Rlist --- XT_LIST (CONS lists)

class Rlist : public Rexp {
public:
    Rexp *head, *tag;
    Rlist *tail;

    Rlist(Rmessage *msg) : Rexp(msg)
    { head=tag=0; tail=0; fix_content(); }
    
    Rlist(unsigned int *ipos, Rmessage *imsg) : Rexp(ipos, imsg)
    { head=tag=0; tail=0; fix_content(); }

    /* this is a sort of special constructor that allows to create a Rlist
       based solely on its content. This is necessary since 0.5 because
       each LISTSXP is no longer represented by its own encoded SEXP
       but they are packed in one content list instead */
    Rlist(int type, Rexp *head, Rexp *tag, char *next, Rmessage *imsg) : Rexp(type, 0, 0, 0) { this->head = head; this->tag = tag; tail = 0; this->next = next; this->msg = imsg; master = 0; }

    virtual ~Rlist();
    
    Rexp *entryByTagName(const char *tagName)  {
      if (tag && (tag->type==XT_SYM || tag->type==XT_SYMNAME) && !strcmp(((Rsymbol*)tag)->symbolName(),tagName)) return head;
        if (tail) return tail->entryByTagName(tagName);
        return 0;
    }

    virtual std::ostream& os_print (std::ostream& os) {
        os << "Rlist[tag=";
        if (tag) os << *tag; else os << "<none>";
        os << ",head=";
        if (head) os << *head; else os << "<none>";
        if (tail) os << ",tail=" << *tail;
        return os << "]";
    }
    
private:
    void fix_content();
};

//===================================== Rvector --- XT_VECTOR (general lists)

class Rvector : public Rexp {
protected:
    Rexp **cont;
    int count;

    // cached
    char **strs;
public:
    Rvector(Rmessage *msg) : Rexp(msg)
    { cont=0; count=0; strs=0; fix_content(); }
    
    Rvector(unsigned int *ipos, Rmessage *imsg) : Rexp(ipos, imsg)
    { cont=0; count=0; strs=0; fix_content(); }
  
    virtual ~Rvector();
    
    char **strings();
    int indexOf(Rexp *exp);
    int indexOfString(const char *str);
    virtual Rsize_t length() { return (Rsize_t) count; }

    char *stringAt(int i) {
        if (i < 0 || i >= count || !cont[i] || cont[i]->type != XT_STR) return 0;
        return ((Rstring*)cont[i])->string();
    }

    Rexp *elementAt(int i) {
	return (i < 0 || i >= count || !cont[i]) ? 0 : cont[i];
    }
    
    Rexp* byName(const char *name);

    virtual std::ostream& os_print (std::ostream& os) {
        os << "Rvector[count=" << count << ":";
        int i=0;
        while (i<count) {
            if (cont[i]) os << *cont[i]; else os << "NULL";
            i++;
            if (i<count) os << ",";
        }
        return os << "]";
    }
private:
    int capacity;
    void fix_content();
};

//===================================== Rconnection ---- Rserve interface class

class Rconnection;

class Rsession {
protected:
    char *host_;
    int port_;
    char key_[32];
    
public:
    Rsession(const char *host, int port, const char key[32]) {
	host_ = host ? strdup(host) : 0;
	port_ = port;
	memcpy(key_, key, 32);
    }
    
    ~Rsession() {
	if (host_) free(host_);
    }
    
    const char *host() { return host_; }
    int port() { return port_; }
    const char *key() { return key_; }
};

class Rconnection {
protected:
    char *host;
    int  port;
    SOCKET s;
    int  family;
    int auth;
    char salt[2];
    char *session_key;

public:
    /** host - either host name or unix socket path
        port - either TCP port or -1 if unix sockets should be used */
    Rconnection(const char *host="127.0.0.1", int port=default_Rsrv_port);
    Rconnection(Rsession *session);
    
    virtual ~Rconnection();
    
    int connect();
    int disconnect();
    
    /**--- low-level functions (should not be used directly) --- */
    
    int request(Rmessage *msg, int cmd, int len=0, void *par=0);
    int request(Rmessage *targetMsg, Rmessage *contents);
    
    /** --- high-level functions --- */
    
    int assign(const char *symbol, Rexp *exp);
    int voidEval(const char *cmd);
    Rexp *eval(const char *cmd, int *status=0, int opt=0);
    int login(const char *user, const char *pwd);
    int shutdown(const char *key);

    /*      ( I/O functions )     */
    int openFile(const char *fn);
    int createFile(const char *fn);
    int readFile(char *buf, unsigned int len);
    int writeFile(const char *buf, unsigned int len);
    int closeFile();
    int removeFile(const char *fn);

    /* session methods - results of detach [if not NULL] must be deleted by the caller when no longer needed! */
    Rsession *detachedEval(const char *cmd, int *status = 0);
    Rsession *detach(int *status = 0);
    // sessions are resumed using resume() method of the Rsession object
    
#ifdef CMD_ctrl
    /* server control functions (need Rserve 0.6-0 or higher) */
    int serverEval(const char *cmd);
    int serverSource(const char *fn);
    int serverShutdown();
#endif
};

#endif
