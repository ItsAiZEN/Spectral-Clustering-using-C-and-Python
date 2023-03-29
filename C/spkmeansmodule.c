#include "spkmeans.h"

//static PyObject* fit(PyObject *self, PyObject *args) {
//    int num_of_clusters;
//    int num_of_iterations;
//    int vector_dimension;
//    int count;
//    PyObject *pyObjVector_list;
//    double eps;
//    PyObject *pyObjInit_centroids;
//    PyObject *item_row;
//    PyObject *item_col;
//    double num;
//    if (!PyArg_ParseTuple(args, "iiiiOdO", &num_of_clusters, &num_of_iterations, &vector_dimension, &count, &pyObjVector_list, &eps, &pyObjInit_centroids)) {
//        return NULL;
//    }
//
//    if (vector_dimension < 0 || count < 0){
//        return NULL;
//    }
//
//    double vector_list[count][vector_dimension];
//    int i;
//    int j;
//    for (i = 0; i < count; i++) {
//        item_row = PyList_GetItem(pyObjVector_list, i);
//        for (j = 0; j < vector_dimension; j++) {
//            item_col = PyList_GetItem(item_row, j);
//            num = PyFloat_AsDouble(item_col);
//            vector_list[i][j] = num;
//        }
//    }
//
//    double init_centroids[count][vector_dimension];
//    for (i = 0; i < num_of_clusters; i++) {
//        item_row = PyList_GetItem(pyObjInit_centroids, i);
//        for (j = 0; j < vector_dimension; j++) {
//            item_col = PyList_GetItem(item_row, j);
//            num = PyFloat_AsDouble(item_col);
//            init_centroids[i][j] = num;
//        }
//    }
//
//    kmeans(num_of_clusters, num_of_iterations, vector_dimension, count, vector_list, eps, init_centroids);
//
//    Py_RETURN_NONE;
//}
//
//static PyMethodDef kmeansMethods[] = {
//        {
//                "fit", // name exposed to Python
//                fit, // C wrapper function
//                     METH_VARARGS, // received variable args (but really just 1)
//                "Performs the kmeans algorithm according to the given centroid initialization" // documentation
//        },
//        {NULL, NULL, 0, NULL}
//};
//
//static struct PyModuleDef mykmeanssp = {
//        PyModuleDef_HEAD_INIT,
//        "mykmeanssp",     // name of module exposed to Python
//        "A module that performs the kmeans algorithm", // module documentation
//        -1,
//        kmeansMethods
//};
//
//PyMODINIT_FUNC PyInit_mykmeanssp(void)
//{
//    PyObject *m;
//    m = PyModule_Create(&mykmeanssp);
//    if (!m) {
//        return NULL;
//    }
//    return m;
//}