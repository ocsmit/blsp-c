#include <Python.h>
#include <stdio.h>
#include "libblsp.h"

#if PY_MAJOR_VERSION < 3
#error "Requires Python 3"
#include "stopcompilation" // Non existant file ensures stop
#endif                     // PY_MAJOR_VERSION

static PyObject *hello(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args)) {
  printf("Hello world\n");
  Py_RETURN_NONE;
}

static PyMethodDef methods[] = {
    {"hello", &hello, METH_VARARGS, "Hello world function"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT, /* always required */
    "blsp",
    /* module name */
    "testing module",      /* description */
    -1,                    /* module size, -1 means we don't use this feature */
    methods,               /* method table */
};

PyMODINIT_FUNC PyInit_blsp() {
  printf("initializing module\n");
  return PyModule_Create(&module_def);
}
