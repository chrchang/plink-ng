/*  system independent sockets (basically for unix and Win)
 *  Copyright (C) 2000,1 Simon Urbanek
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation; version 2.1 of the License
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
   
   conditional defines: 

   MAIN
     should be defined in just one file that will contain the fn definitions and variables

   USE_SNPRINTF 
     emulate snprintf on Win platforms (you will
     lose the security which is provided under unix of course)

   SOCK_ERRORS
     include error code handling and checking functions
*/

#ifndef __SISOCKS_H__
#define __SISOCKS_H__

#if defined __GNUC__ && !defined unix && !defined Win32 /* MacOS X hack (gcc on any platform should behave as unix - except for Win32, where we need to keep using winsock) */
#define unix
#endif

#if defined STANDALONE_RSERVE && defined RSERVE_PKG
#undef RSERVE_PKG
#endif

#if defined SOCK_ERRORS || defined USE_SNPRINTF
#include <stdio.h>
#endif
#include <string.h>

#ifdef unix
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <errno.h>
#include <stdlib.h>

#define sockerrno errno

#define SOCKET int
#define INVALID_SOCKET (-1)
#define closesocket(A) close(A)

#else

#define windows
#include <winsock2.h>
#include <windows.h>
#include <string.h>
#include <stdlib.h>
#define inet_aton(A,B) (0, B.s_addr=inet_addr(A))

#define sockerrno WSAGetLastError()

#ifndef WIN64
#define ECONNREFUSED WSAECONNREFUSED
#define EADDRINUSE WSAEADDRINUSE
#define ENOTSOCK WSAENOTSOCK
#define EISCONN WSAEISCONN
#define ETIMEDOUT WSAETIMEDOUT
#define ENETUNREACH WSAENETUNREACH
#define EINPROGRESS WSAEINPROGRESS
#define EALREADY WSAEALREADY
#define EAFNOSUPPORT WSAEAFNOSUPPORT
#define EBADF WSAEBADF
#define EINVAL WSAEINVAL
#define EOPNOTSUPP WSAEOPNOTSUPP
#define EFAULT WSAEFAULT
#define EWOULDBLOCK WSAEWOULDBLOCK
#define EACCES WSAEACCES
#else
#define ECONNREFUSED WSAECONNREFUSED
#define EADDRINUSE WSAEADDRINUSE
#define ENOTSOCK WSAENOTSOCK
#define EISCONN WSAEISCONN
#define ETIMEDOUT WSAETIMEDOUT
#define ENETUNREACH WSAENETUNREACH
#define EINPROGRESS WSAEINPROGRESS
#define EALREADY WSAEALREADY
#define EOPNOTSUPP WSAEOPNOTSUPP
#define EWOULDBLOCK WSAEWOULDBLOCK
#endif

#ifdef USE_SNPRINTF
#ifdef MAIN
int snprintf(char *buf, int len, char *fmt, ...)
{
   va_list argptr;
   int cnt;

   va_start(argptr, fmt);
   cnt = vsprintf(buf, fmt, argptr);
   va_end(argptr);

   return(cnt);
}
#else
extern int snprintf(char *buf, int len, char *fmt, ...);
#endif
#endif

#endif

#define SA struct sockaddr
#define SAIN struct sockaddr_in

#ifdef windows

#ifdef MAIN
int initsocks(void)
{
  WSADATA dt;
  /* initialize WinSock 1.1 */
  return (WSAStartup(0x0101,&dt))?-1:0;
}
#else
extern int initsocks(void);
#endif

#define donesocks() WSACleanup()
#else
 
/* no stupid stuff necessary for unix */
#define initsocks()
#define donesocks()

#endif

#ifdef SOCK_ERRORS

#ifdef MAIN
int suppmode=0;
int socklasterr;
FILE *sockerrlog=0;

/* copy error description to buf or set *buf=0 if none */
int sockerrorchecks(char *buf, int blen, int res) {
  *buf=0;
  if (res==-1) {
    switch(sockerrno) {
    case EBADF: strncpy(buf,"bad descriptor",blen); break;
    case EINVAL: strncpy(buf,"already in use",blen); break;
    case EACCES: strncpy(buf,"access denied",blen); break;
    case ENOTSOCK: strncpy(buf,"descriptor is not a socket",blen); break;
    case EOPNOTSUPP: strncpy(buf,"operation not supported",blen); break;
    case EFAULT: strncpy(buf,"fault",blen); break;
    case EWOULDBLOCK: strncpy(buf,"operation would block",blen); break;
    case EISCONN: strncpy(buf,"is already connected",blen); break;
    case ECONNREFUSED: strncpy(buf,"connection refused",blen); break;
    case ETIMEDOUT: strncpy(buf,"operation timed out",blen); break;
    case ENETUNREACH: strncpy(buf,"network is unreachable",blen); break;
    case EADDRINUSE: strncpy(buf,"address already in use",blen); break;
    case EINPROGRESS: strncpy(buf,"in progress",blen); break;
    case EALREADY: strncpy(buf,"previous connect request not completed yet",blen); break;
#ifdef unix
    default: snprintf(buf,blen,"unknown socket error %d",sockerrno);
#else
    default: sprintf(buf,"unknown socket error %d",sockerrno);
#endif
    }
  }
  return res;
}

#ifdef RSERVE_PKG /* use error instead of exit */

#include <R.h>

int sockerrorcheck(char *sn, int rtb, int res) {
    if ((signed int)res == -1) {
	char sock_err_buf[72];
	sockerrorchecks(sock_err_buf, sizeof(sock_err_buf), res);
	if (rtb)
	    Rf_error("%s socket error #%d (%s)", sn, (int) sockerrno, sock_err_buf);
	else
	    Rf_warning("%s socket error #%d (%s)", sn, (int) sockerrno, sock_err_buf);
    }
    return res;
}

#else

/* check socket error and add to log file if necessary */
int sockerrorcheck(char *sn, int rtb, int res) {
  if (!sockerrlog) sockerrlog=stderr;
  if ((signed int)res==-1) {
    if (socklasterr==sockerrno) {
      suppmode++;
    } else {
      if (suppmode>0) {
        fprintf(sockerrlog,"##> REP: (last error has been repeated %d times.)\n",suppmode);
        suppmode=0;
      }
      fprintf(sockerrlog,"##> SOCK_ERROR: %s error #%d",sn,sockerrno);
      switch(sockerrno) {
      case EBADF: fprintf(sockerrlog,"(bad descriptor)"); break;
      case EINVAL: fprintf(sockerrlog,"(already in use)"); break;
      case EACCES: fprintf(sockerrlog,"(access denied)"); break;
      case ENOTSOCK: fprintf(sockerrlog,"(descriptor is not a socket)"); break;
      case EOPNOTSUPP: fprintf(sockerrlog,"(operation not supported)"); break;
      case EFAULT: fprintf(sockerrlog,"(fault)"); break;
      case EWOULDBLOCK: fprintf(sockerrlog,"(operation would block)"); break;
      case EISCONN: fprintf(sockerrlog,"(is already connected)"); break;
      case ECONNREFUSED: fprintf(sockerrlog,"(connection refused)"); break;
      case ETIMEDOUT: fprintf(sockerrlog,"(operation timed out)"); break;
      case ENETUNREACH: fprintf(sockerrlog,"(network is unreachable)"); break;
      case EADDRINUSE: fprintf(sockerrlog,"(address already in use)"); break;
      case EINPROGRESS: fprintf(sockerrlog,"(in progress)"); break;
      case EALREADY: fprintf(sockerrlog,"(previous connect request not completed yet)"); break;
      default: fprintf(sockerrlog,"(?)");
      }
      fprintf(sockerrlog,"\n"); fflush(sockerrlog);
      socklasterr=sockerrno;
    }
    if (rtb) exit(1);
  }
  return res;
}
#endif

#else
extern int suppmode;
extern int socklasterr;
extern FILE *sockerrlog;

int sockerrorchecks(char *buf, int blen, int res);
int sockerrorcheck(char *sn, int rtb, int res);
#endif
    
#define FCF(X,F) sockerrorcheck(X,1,F)
#define CF(X,F) sockerrorcheck(X,0,F)

#endif

#ifdef MAIN
struct sockaddr *build_sin(struct sockaddr_in *sa,char *ip,int port) {
  memset(sa,0,sizeof(struct sockaddr_in));
  sa->sin_family=AF_INET;
  sa->sin_port=htons(port);
  sa->sin_addr.s_addr=(ip)?inet_addr(ip):htonl(INADDR_ANY);
  return (struct sockaddr*)sa;
}
#else
struct sockaddr *build_sin(struct sockaddr_in *sa,char *ip,int port);
#endif
          
#endif /* __SISOCKS_H__ */
