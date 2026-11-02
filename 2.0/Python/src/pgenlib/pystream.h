#ifndef __PGENLIB_PYSTREAM_H__
#define __PGENLIB_PYSTREAM_H__

// Wraps a Python binary file-like object (io.BytesIO, an open file, an
// fsspec/s3fs file, ...) in a read-only FILE*, so that pgenlib_read can use it
// without modification.  Uses funopen() on macOS/BSD and fopencookie() on
// Linux; other platforms (Windows) are not supported.
//
// The callbacks may run while the GIL is released, so each one reacquires it.
// A Python exception raised inside a callback is appended to the caller's
// `errors` list and reported to stdio as an I/O error.
//
// fclose() drops our references; it does not close the Python object, which
// remains owned by the caller.
//
// A reader that is never closed keeps its stream open past interpreter
// shutdown, and glibc's exit() then syncs every open stream, which calls the
// seek callback after Py_Finalize().  So once Python is finalizing, the
// callbacks fail without touching it.

#include <Python.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct PgenlibPyStreamStruct {
  PyObject* obj;
  PyObject* errors;
  int has_readinto;
} PgenlibPyStream;

static int PgenlibPyStreamPythonUnavailable(void) {
#if PY_VERSION_HEX >= 0x030D0000
  return (!Py_IsInitialized()) || Py_IsFinalizing();
#else
  return (!Py_IsInitialized()) || _Py_IsFinalizing();
#endif
}

// Assumes the GIL is held and a Python exception is set.
static void PgenlibPyStreamStashError(PgenlibPyStream* psp) {
  PyObject* exc_type;
  PyObject* exc_value;
  PyObject* exc_tb;
  PyErr_Fetch(&exc_type, &exc_value, &exc_tb);
  PyErr_NormalizeException(&exc_type, &exc_value, &exc_tb);
  if (exc_value) {
    if (exc_tb) {
      PyException_SetTraceback(exc_value, exc_tb);
    }
    if (PyList_Append(psp->errors, exc_value)) {
      PyErr_Clear();
    }
  }
  Py_XDECREF(exc_type);
  Py_XDECREF(exc_value);
  Py_XDECREF(exc_tb);
  errno = EIO;
}

// Fills buf completely unless EOF is reached first.  Returns the number of
// bytes read, or -1 on error.
static int64_t PgenlibPyStreamRead(void* cookie, char* buf, size_t size) {
  PgenlibPyStream* psp = (PgenlibPyStream*)cookie;
  if (PgenlibPyStreamPythonUnavailable()) {
    errno = EIO;
    return -1;
  }
  PyGILState_STATE gstate = PyGILState_Ensure();
  int64_t total = 0;
  while ((size_t)total < size) {
    const size_t remaining = size - (size_t)total;
    Py_ssize_t cur_len;
    if (psp->has_readinto) {
      PyObject* view = PyMemoryView_FromMemory(&(buf[total]), (Py_ssize_t)remaining, PyBUF_WRITE);
      if (!view) {
        goto PgenlibPyStreamRead_fail;
      }
      PyObject* result = PyObject_CallMethod(psp->obj, "readinto", "O", view);
      if (!result) {
        Py_DECREF(view);
        goto PgenlibPyStreamRead_fail;
      }
      // Don't let the callee keep a view into stdio's buffer.
      PyObject* release_result = PyObject_CallMethod(view, "release", NULL);
      Py_DECREF(view);
      if (!release_result) {
        Py_DECREF(result);
        goto PgenlibPyStreamRead_fail;
      }
      Py_DECREF(release_result);
      if (result == Py_None) {
        Py_DECREF(result);
        PyErr_SetString(PyExc_BlockingIOError, "pgenlib: readinto() returned None (non-blocking stream?)");
        goto PgenlibPyStreamRead_fail;
      }
      cur_len = PyLong_AsSsize_t(result);
      Py_DECREF(result);
      if ((cur_len == -1) && PyErr_Occurred()) {
        goto PgenlibPyStreamRead_fail;
      }
    } else {
      PyObject* result = PyObject_CallMethod(psp->obj, "read", "n", (Py_ssize_t)remaining);
      if (!result) {
        goto PgenlibPyStreamRead_fail;
      }
      char* data;
      if (PyBytes_AsStringAndSize(result, &data, &cur_len)) {
        Py_DECREF(result);
        goto PgenlibPyStreamRead_fail;
      }
      if ((size_t)cur_len <= remaining) {
        memcpy(&(buf[total]), data, (size_t)cur_len);
      }
      Py_DECREF(result);
    }
    if ((cur_len < 0) || ((size_t)cur_len > remaining)) {
      PyErr_SetString(PyExc_OSError, "pgenlib: file-like object returned an invalid read length");
      goto PgenlibPyStreamRead_fail;
    }
    if (!cur_len) {
      break;
    }
    total += cur_len;
  }
  PyGILState_Release(gstate);
  return total;
 PgenlibPyStreamRead_fail:
  PgenlibPyStreamStashError(psp);
  PyGILState_Release(gstate);
  return -1;
}

// Returns the new absolute position, or -1 on error.
static int64_t PgenlibPyStreamSeek(void* cookie, int64_t offset, int whence) {
  PgenlibPyStream* psp = (PgenlibPyStream*)cookie;
  if (PgenlibPyStreamPythonUnavailable()) {
    errno = EIO;
    return -1;
  }
  PyGILState_STATE gstate = PyGILState_Ensure();
  int64_t new_pos = -1;
  PyObject* result = PyObject_CallMethod(psp->obj, "seek", "Li", (long long)offset, whence);
  if (result) {
    new_pos = PyLong_AsLongLong(result);
    Py_DECREF(result);
  }
  if (new_pos < 0) {
    if (!PyErr_Occurred()) {
      PyErr_SetString(PyExc_OSError, "pgenlib: file-like object returned an invalid seek position");
    }
    PgenlibPyStreamStashError(psp);
    new_pos = -1;
  }
  PyGILState_Release(gstate);
  return new_pos;
}

static int PgenlibPyStreamClose(void* cookie) {
  PgenlibPyStream* psp = (PgenlibPyStream*)cookie;
  if (!PgenlibPyStreamPythonUnavailable()) {
    PyGILState_STATE gstate = PyGILState_Ensure();
    Py_DECREF(psp->obj);
    Py_DECREF(psp->errors);
    PyGILState_Release(gstate);
  }
  // else the references are deliberately leaked.
  free(psp);
  return 0;
}

#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__DragonFly__)
static int PgenlibPyStreamReadBsd(void* cookie, char* buf, int size) {
  return (int)PgenlibPyStreamRead(cookie, buf, (size_t)size);
}

static fpos_t PgenlibPyStreamSeekBsd(void* cookie, fpos_t offset, int whence) {
  return (fpos_t)PgenlibPyStreamSeek(cookie, (int64_t)offset, whence);
}
#  define PGENLIB_PYSTREAM_SUPPORTED
#elif defined(__linux__)
// Python.h defines _GNU_SOURCE on Linux, so fopencookie() is declared.
static ssize_t PgenlibPyStreamReadGnu(void* cookie, char* buf, size_t size) {
  return (ssize_t)PgenlibPyStreamRead(cookie, buf, size);
}

#  ifdef __GLIBC__
typedef off64_t PgenlibPyStreamOff;
#  else
// musl's off_t is always 64-bit, and its fopencookie() uses it directly.
typedef off_t PgenlibPyStreamOff;
#  endif

static int PgenlibPyStreamSeekGnu(void* cookie, PgenlibPyStreamOff* offset_ptr, int whence) {
  const int64_t new_pos = PgenlibPyStreamSeek(cookie, (int64_t)(*offset_ptr), whence);
  if (new_pos < 0) {
    return -1;
  }
  *offset_ptr = (PgenlibPyStreamOff)new_pos;
  return 0;
}
#  define PGENLIB_PYSTREAM_SUPPORTED
#endif

// Assumes the GIL is held.  Returns nullptr with a Python exception set on
// failure.
static FILE* PgenlibPyStreamOpen(PyObject* obj, PyObject* errors) {
#ifdef PGENLIB_PYSTREAM_SUPPORTED
  if (!PyList_Check(errors)) {
    PyErr_SetString(PyExc_TypeError, "pgenlib: errors must be a list");
    return NULL;
  }
  static const char* kRequiredMethods[] = {"read", "seek"};
  for (size_t uii = 0; uii != sizeof(kRequiredMethods) / sizeof(kRequiredMethods[0]); ++uii) {
    if (!PyObject_HasAttrString(obj, kRequiredMethods[uii])) {
      PyErr_Format(PyExc_TypeError, "file-like object must have a %s() method", kRequiredMethods[uii]);
      return NULL;
    }
  }
  PgenlibPyStream* psp = (PgenlibPyStream*)malloc(sizeof(PgenlibPyStream));
  if (!psp) {
    PyErr_NoMemory();
    return NULL;
  }
  psp->obj = obj;
  psp->errors = errors;
  psp->has_readinto = PyObject_HasAttrString(obj, "readinto");
#  ifdef __linux__
  cookie_io_functions_t io_funcs;
  io_funcs.read = PgenlibPyStreamReadGnu;
  io_funcs.write = NULL;
  io_funcs.seek = PgenlibPyStreamSeekGnu;
  io_funcs.close = PgenlibPyStreamClose;
  FILE* ff = fopencookie(psp, "rb", io_funcs);
#  else
  FILE* ff = funopen(psp, PgenlibPyStreamReadBsd, NULL, PgenlibPyStreamSeekBsd, PgenlibPyStreamClose);
#  endif
  if (!ff) {
    free(psp);
    PyErr_SetFromErrno(PyExc_OSError);
    return NULL;
  }
  Py_INCREF(obj);
  Py_INCREF(errors);
  return ff;
#else
  (void)obj;
  (void)errors;
  PyErr_SetString(PyExc_NotImplementedError, "pgenlib: file-like objects are not supported on this platform; pass a filename instead");
  return NULL;
#endif
}

#endif  // __PGENLIB_PYSTREAM_H__
