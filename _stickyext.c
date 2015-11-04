#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

static PyObject* fastsub(PyObject* self, PyObject* args) {

    PyArrayObject *py_x;
    PyArrayObject      *py_r;
    PyArrayIterObject *itr;
    double *p1,*res;
    double g,d;
    int axis = 1;
    int go;
    int i;
    
    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &py_x, &PyArray_Type, &py_r))
        return NULL;
    
    g = 0;
    d = 0;
    itr = (PyArrayIterObject *) PyArray_IterAllButAxis((PyObject*)py_x,&axis);
    while(PyArray_ITER_NOTDONE(itr)) {
        go = py_x->strides[axis]/sizeof(double);
        p1 = (double *) PyArray_ITER_DATA(itr);
        res = (double *) PyArray_GETPTR1(py_r,itr->index);
        g = 0;
        d = 0;
        for (i = 0; i < py_x->dimensions[axis]; i++) {
            d+=*p1;
            if (d>g) g=d;
            if ((*p1)==0) d=0;
            p1+=go;
        }
        *res = g;
        PyArray_ITER_NEXT(itr);
    }
    Py_RETURN_NONE;
    
}

/*  define functions in module */
static PyMethodDef FastSub[] =
{
     {"fastsub", fastsub, METH_VARARGS,
         "fastsub!"},
     {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC

init_stickyext(void)
{
     (void) Py_InitModule("_stickyext", FastSub);
     /* IMPORTANT: this must be called */
     import_array();
}

