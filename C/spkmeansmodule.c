#include "spkmeans.h"

static PyObject* wam(PyObject *self, PyObject *args) {
    double **pyObj_vectors;
    int num_of_vectors;
    int vector_dimension;
    double **Cwam_matrix;
    double **vector_list;
    int i;
    int j;
    double num;
    double wam_num;
    PyObject *item_row;
    PyObject *item_col;
    PyObject *py_wam_matrix;
    PyObject *py_wam_matrix_row;
    PyObject *return_value;

    // check args are as expected
    if (!PyArg_ParseTuple(args, "Oii", &pyObj_vectors, &num_of_vectors, &vector_dimension)) {
        return NULL;
    }

    // transform input vector list (given in Python) to C so that we can use it as an argument in the C wam function
    vector_list = (double **)malloc(num_of_vectors * sizeof(double *));
    if (vector_list == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        vector_list[i] = (double *)malloc(num_of_vectors * sizeof(double));
        if (vector_list[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        item_row = PyList_GetItem(pyObj_vectors, i);
        for (j = 0; j < vector_dimension; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            vector_list[i][j] = num;
        }
    }

    // allocate memory for returned wam_matrix in C (for returned val of wam function)
    Cwam_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (Cwam_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        Cwam_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (Cwam_matrix[i] == NULL) {
            print_error();
        }
    }

    Cwam_matrix = wam(vector_list, num_of_vectors, vector_dimension);

    // create new python list (to return to python) and insert Cwam_matrix into the list
    py_wam_matrix = PyList_New(num_of_vectors); // PyList_New returns a list of size num_of_vectors on success, and NULL on failure
    if (py_wam_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        py_wam_matrix_row = PyList_New(vector_dimension);
        PyList_SetItem(py_wam_matrix, i, py_wam_index_row);
        for (j = 0; j < vector_dimension; j++) {
            wam_num = PyFloat_FromDouble(Cwam_matrix[i][j]);
            PyList_SetItem(py_wam_matrix_row, j, wam_num);
        }
    }

    // free memory (Cwam_matrix + vector_list in C)
    for (i = 0; i < num_of_vectors; i++) {
        free(Cwam_matrix[i]);
    }
    free(Cwam_matrix);
    for (i = 0; i < num_of_vectors; i++) {
        free(vector_list[i]);
    }
    free(vector_list);

    // we are returning the python list that contains the wam matrix
    return_val = Py_BuildValue("O", py_wam_matrix); // not sure if this should be "O" or maybe "O&";
    return return_val;

}

static PyObject* ddg(PyObject *self, PyObject *args) {
    double **Cddg_matrix;
    int i;
    int j;
    int num_of_vectors;
    int vector_dimension;
    double num;
    double ddg_num;
    double **C_wam_matrix;
    PyObject *return_val;
    PyObject *py_wam_matrix;
    PyObject *item_row;
    PyObject *item_col;
    PyObject *py_ddg_matrix;
    PyObject *py_ddg_matrix_row;

    // check args are as expected
    if (!PyArg_ParseTuple(args, "Oi", &py_wam_matrix, &num_of_vectors)) {
        return NULL;
    }

    // transform input wam matrix (given as list in Python) to C so that we can use it as an argument in the C ddg function
    C_wam_matrix = (double **)malloc(num_of_vectors * sizeof(double *));
    if (C_wam_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        C_wam_matrix[i] = (double *)malloc(num_of_vectors * sizeof(double));
        if (C_wam_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        item_row = PyList_GetItem(py_wam_matrix, i);
        for (j = 0; j < vector_dimension; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            C_wam_matrix[i][j] = num;
        }
    }

    // allocate memory for returned ddg_matrix in C (for returned val of ddg function)
    Cddg_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (Cddg_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        Cddg_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (Cddg_matrix[i] == NULL) {
            print_error();
        }
    }

    Cddg_matrix = ddg(C_wam_matrix, num_of_vectors);

    // create new python list (to return to python) and insert Cddg_matrix into the list
    py_ddg_matrix = PyList_New(num_of_vectors); // PyList_New returns a list of size num_of_vectors on success, and NULL on failure
    if (py_ddg_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        py_ddg_matrix_row = PyList_New(vector_dimension);
        PyList_SetItem(py_ddg_matrix, i, py_ddg_index_row);
        for (j = 0; j < vector_dimension; j++) {
            ddg_num = PyFloat_FromDouble(Cddg_matrix[i][j]);
            PyList_SetItem(py_ddg_matrix_row, j, ddg_num);
        }
    }

    // free memory (Cddg_matrix + C_wam_matrix)
    for (i = 0; i < num_of_vectors; i++) {
        free(Cddg_matrix[i]);
    }
    free(Cddg_matrix);
    for (i = 0; i < num_of_vectors; i++) {
        free(C_wam_matrix[i]);
    }
    free(C_wam_matrix);

    // we are returning the python list that contains the ddg matrix
    return_val = Py_BuildValue("O", py_ddg_matrix); // not sure if this should be "O" or maybe "O&";
    return return_val;

}

static PyObject* gl(PyObject *self, PyObject *args) {
    int num_of_vectors;
    PyObject *py_ddg_matrix;
    PyObject *py_wam_matrix;

    // check args are as expected
    if (!PyArg_ParseTuple(args, "OOi", &py_ddg_matrix, &py_wam_matrix, &num_of_vectors)) {
        return NULL;
    }

    // not finished. maybe we can use the pyObj functions for wam and ddg above instead of copying everything ?
}

static PyObject* jacobi(PyObject *self, PyObject *args) {

}

static PyObject* spk(PyObject *self, PyObject *args) {

}



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