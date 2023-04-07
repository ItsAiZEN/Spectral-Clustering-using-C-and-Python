
#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include "spkmeans.h"

static PyObject *wam1(PyObject *self, PyObject *args) {
    int num_of_vectors;
    int vector_dimension;
    double **Cwam_matrix;
    double **vector_list;
    int i;
    int j;
    double num;
    PyObject *wam_num;
    PyObject *pyObj_vectors;
    PyObject *item_row;
    PyObject *item_col;
    PyObject *py_wam_matrix;
    PyObject *py_wam_matrix_row;
    PyObject *return_val;

    // check args are as expected
    if (!PyArg_ParseTuple(args, "Oii", &pyObj_vectors, &num_of_vectors, &vector_dimension)) {
        return NULL;
    }

    // transform input vector list (given in Python) to C so that we can use it as an argument in the C wam function
    vector_list = (double **) malloc(num_of_vectors * sizeof(double *));
    if (vector_list == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        vector_list[i] = (double *) malloc(num_of_vectors * sizeof(double));
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
        PyList_SetItem(py_wam_matrix, i, py_wam_matrix_row);
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

static PyObject *ddg1(PyObject *self, PyObject *args) {
    double **Cddg_matrix;
    int i;
    int j;
    int num_of_vectors;
    double num;
    double **C_wam_matrix;
    PyObject *ddg_num;
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
    C_wam_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (C_wam_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        C_wam_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (C_wam_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        item_row = PyList_GetItem(py_wam_matrix, i);
        for (j = 0; j < num_of_vectors; j++) {
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
    py_ddg_matrix = PyList_New(
            num_of_vectors); // PyList_New returns a list of size num_of_vectors on success, and NULL on failure
    if (py_ddg_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        py_ddg_matrix_row = PyList_New(num_of_vectors);
        PyList_SetItem(py_ddg_matrix, i, py_ddg_matrix_row);
        for (j = 0; j < num_of_vectors; j++) {
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

static PyObject *gl1(PyObject *self, PyObject *args) {
    int num_of_vectors;
    double **C_wam_matrix;
    double **C_ddg_matrix;
    double **C_gl_matrix;
    int i;
    int j;
    double num;
    PyObject *gl_num;
    PyObject *item_row;
    PyObject *item_col;
    PyObject *return_val;
    PyObject *py_ddg_matrix;
    PyObject *py_wam_matrix;
    PyObject *py_gl_matrix;
    PyObject *py_gl_matrix_row;


    // check args are as expected
    if (!PyArg_ParseTuple(args, "OOi", &py_ddg_matrix, &py_wam_matrix, &num_of_vectors)) {
        return NULL;
    }

    // not finished. maybe we can use the pyObj functions for wam and ddg above instead of copying everything ?

    // transform input wam matrix (given as list in Python) to C so that we can use it as an argument in the C ddg function
    C_wam_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (C_wam_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        C_wam_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (C_wam_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        item_row = PyList_GetItem(py_wam_matrix, i);
        for (j = 0; j < num_of_vectors; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            C_wam_matrix[i][j] = num;
        }
    }

    // transform input ddg matrix (given as list in Python) to C so that we can use it as an argument in the C ddg function
    C_ddg_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (C_ddg_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        C_ddg_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (C_ddg_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        item_row = PyList_GetItem(py_ddg_matrix, i);
        for (j = 0; j < num_of_vectors; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            C_ddg_matrix[i][j] = num;
        }
    }

    // allocate memory for returned gl_matrix in C (for returned val of gl function)
    C_gl_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (C_gl_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        C_gl_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (C_gl_matrix[i] == NULL) {
            print_error();
        }
    }

    C_gl_matrix = gl(C_wam_matrix, C_ddg_matrix, num_of_vectors);

    // create new python list (to return to python) and insert C_gl_matrix into the list
    py_gl_matrix = PyList_New(
            num_of_vectors); // PyList_New returns a list of size num_of_vectors on success, and NULL on failure
    if (py_gl_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        py_gl_matrix_row = PyList_New(num_of_vectors);
        PyList_SetItem(py_gl_matrix, i, py_gl_matrix_row);
        for (j = 0; j < num_of_vectors; j++) {
            gl_num = PyFloat_FromDouble(C_gl_matrix[i][j]);
            PyList_SetItem(py_gl_matrix_row, j, gl_num);
        }
    }

    // free memory (C_ddg_matrix + C_wam_matrix)
    for (i = 0; i < num_of_vectors; i++) {
        free(C_ddg_matrix[i]);
    }
    free(C_ddg_matrix);
    for (i = 0; i < num_of_vectors; i++) {
        free(C_wam_matrix[i]);
    }
    free(C_wam_matrix);

    // we are returning the python list that contains the gl matrix
    return_val = Py_BuildValue("O", py_gl_matrix); // not sure if this should be "O" or maybe "O&";
    return return_val;
}

static PyObject *jacobi1(PyObject *self, PyObject *args) {
    int num_of_vectors;
    double **C_matrix;
    double num;
    PyObject *jacobi_num;
    int i;
    int j;
    double **C_jacobi_matrix;
    PyObject *return_val;
    PyObject *py_matrix;
    PyObject *item_row;
    PyObject *item_col;
    PyObject *py_jacobi_matrix;
    PyObject *py_jacobi_matrix_row;


    // check args are as expected
    if (!PyArg_ParseTuple(args, "Oi", &py_matrix, &num_of_vectors)) {
        return NULL;
    }

    // transform input matrix (given as list in Python) to C so that we can use it as an argument in the C ddg function
    C_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (C_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        C_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (C_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        item_row = PyList_GetItem(py_matrix, i);
        for (j = 0; j < num_of_vectors; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            C_matrix[i][j] = num;
        }
    }

    // allocate memory for returned gl_matrix in C (for returned val of gl function)
    C_jacobi_matrix = (double **) malloc((num_of_vectors + 1) * sizeof(double *));
    if (C_jacobi_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < (num_of_vectors + 1); i++) {
        C_jacobi_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (C_jacobi_matrix[i] == NULL) {
            print_error();
        }
    }

    C_jacobi_matrix = jacobi(C_matrix, num_of_vectors);

    // create new python list (to return to python) and insert C_jacobi_matrix into the list

    py_jacobi_matrix = PyList_New(
            num_of_vectors); // PyList_New returns a list of size num_of_vectors on success, and NULL on failure
    if (py_jacobi_matrix == NULL) {
        print_error();
    }

    for (i = 0; i < num_of_vectors; i++) {
        py_jacobi_matrix_row = PyList_New(num_of_vectors);
        PyList_SetItem(py_jacobi_matrix, i, py_jacobi_matrix_row);
        for (j = 0; j < num_of_vectors; j++) {
            jacobi_num = PyFloat_FromDouble(C_jacobi_matrix[i][j]);
            PyList_SetItem(py_jacobi_matrix_row, j, jacobi_num);
        }
    }

    // free memory (C_matrix + C_jacobi_matrix)
    for (i = 0; i < num_of_vectors; i++) {
        free(C_matrix[i]);
    }
    free(C_matrix);
    for (i = 0; i < (num_of_vectors + 1); i++) {
        free(C_jacobi_matrix[i]);
    }
    free(C_jacobi_matrix);


    // we are returning the python list that contains the gl matrix
    return_val = Py_BuildValue("O", py_jacobi_matrix); // not sure if this should be "O" or maybe "O&";
    return return_val;
}

static PyObject *spk(PyObject *self, PyObject *args) {
    int n;
    int k;
    PyObject *pyObjVector_list;
    PyObject *pyObjInit_centroids;
    PyObject *item_row;
    PyObject *item_col;
    double num;
    if (!PyArg_ParseTuple(args, "iiOO", &k, &n, &pyObjVector_list, &pyObjInit_centroids)) {
        return NULL;
    }

    if (n < 0 || k < 0) {
        return NULL;
    }

    // transform py list with vectors (from U) to C matrix
    double vector_list[n][k];
    int i;
    int j;
    for (i = 0; i < n; i++) {
        item_row = PyList_GetItem(pyObjVector_list, i);
        for (j = 0; j < k; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            vector_list[i][j] = num;
        }
    }

    // transform py list with init centroids to C matrix
    double init_centroids[k][k];
    for (i = 0; i < k; i++) {
        item_row = PyList_GetItem(pyObjInit_centroids, i);
        for (j = 0; j < k; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            init_centroids[i][j] = num;
        }
    }

    kmeanspp(k, 300, k, n, vector_list, 0, init_centroids);

    Py_RETURN_NONE;
}

static PyObject *eigengap_heuristic1(PyObject *self, PyObject *args) {
    int num_of_vectors;
    int i;
    int j;
    int k;
    double num;
    double **C_jacobi_matrix;
    PyObject *item_row;
    PyObject *item_col;
    PyObject *return_val;
    PyObject *py_jacobi_matrix;


    if (!PyArg_ParseTuple(args, "Oi", &py_jacobi_matrix, &num_of_vectors)) {
        return NULL;
    }

    // transform input jacobi matrix (given as list in Python) to C so that we can use it as an argument in the C ddg function
    C_jacobi_matrix = (double **) malloc((num_of_vectors + 1) * sizeof(double *));
    if (C_jacobi_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        C_jacobi_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (C_jacobi_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors + 1; i++) {
        item_row = PyList_GetItem(py_jacobi_matrix, i);
        for (j = 0; j < num_of_vectors; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            C_jacobi_matrix[i][j] = num;
        }
    }

    // call eigengap heuristic function
    k = eigengap_heuristic(C_jacobi_matrix, num_of_vectors);

    // free memory (C_jacobi_matrix)
    for (i = 0; i < (num_of_vectors + 1); i++) {
        free(C_jacobi_matrix[i]);
    }
    free(C_jacobi_matrix);

    // we are returning the python list that contains the gl matrix
    return_val = Py_BuildValue("i", k); // not sure if this should be "O" or maybe "O&";
    return return_val;
}

static PyObject *calculateUmatrix1(PyObject *self, PyObject *args) {
    int k;
    int i;
    int j;
    int num_of_vectors;
    double num;
    double **C_jacobi_matrix;
    double **C_U_matrix;
    PyObject *u_num;
    PyObject *py_jacobi_matrix;
    PyObject *item_row;
    PyObject *item_col;
    PyObject *py_u_matrix;
    PyObject *py_u_matrix_row;
    PyObject *return_val;

    if (!PyArg_ParseTuple(args, "Oii", &py_jacobi_matrix, &num_of_vectors, &k)) {
        return NULL;
    }

    // transform input jacobi matrix (given as list in Python) to C so that we can use it as an argument in the C ddg function
    C_jacobi_matrix = (double **) malloc((num_of_vectors + 1) * sizeof(double *));
    if (C_jacobi_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        C_jacobi_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (C_jacobi_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors + 1; i++) {
        item_row = PyList_GetItem(py_jacobi_matrix, i);
        for (j = 0; j < num_of_vectors; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            C_jacobi_matrix[i][j] = num;
        }
    }

    // allocate memory for returned C_U_matrix in C (for returned val of calcU function)
    C_U_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (C_U_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        C_U_matrix[i] = (double *) malloc(k * sizeof(double));
        if (C_U_matrix[i] == NULL) {
            print_error();
        }
    }

    C_U_matrix = calculateUmatrix(C_jacobi_matrix, num_of_vectors, k);

    py_u_matrix = PyList_New(
            num_of_vectors); // PyList_New returns a list of size num_of_vectors on success, and NULL on failure
    if (py_u_matrix == NULL) {
        print_error();
    }

    for (i = 0; i < num_of_vectors; i++) {
        py_u_matrix_row = PyList_New(num_of_vectors);
        PyList_SetItem(py_u_matrix, i, py_u_matrix_row);
        for (j = 0; j < k; j++) {
            u_num = PyFloat_FromDouble(C_U_matrix[i][j]);
            PyList_SetItem(py_u_matrix_row, j, u_num);
        }
    }

    // free memory (C_jacobi_matrix)
    for (i = 0; i < (num_of_vectors + 1); i++) {
        free(C_jacobi_matrix[i]);
    }
    free(C_jacobi_matrix);

    // free memory (C_U_matrix)
    for (i = 0; i < num_of_vectors; i++) {
        free(C_U_matrix[i]);
    }
    free(C_U_matrix);

    // we are returning the python list that contains the U matrix
    return_val = Py_BuildValue("O", py_u_matrix); // not sure if this should be "O" or maybe "O&";
    return return_val;
}


static PyMethodDef kmeanssp_methods[] = {
        {
                "wam", // name exposed to Python
                wam1, // C wrapper function
                      METH_VARARGS, // received variable args (but really just 1)
                "Returns the wam matrix" // documentation
        },
        {
                "ddg", // name exposed to Python
                ddg1, // C wrapper function
                      METH_VARARGS, // received variable args (but really just 1)
                "Returns the ddg matrix" // documentation
        },
        {
                "gl", // name exposed to Python
                gl1, // C wrapper function
                      METH_VARARGS, // received variable args (but really just 1)
                "Returns the gl matrix" // documentation
        },
        {
                "jacobi", // name exposed to Python
                jacobi1, // C wrapper function
                      METH_VARARGS, // received variable args (but really just 1)
                "Returns the jacobi matrix" // documentation
        },
        {
                "spk", // name exposed to Python
                spk, // C wrapper function
                      METH_VARARGS, // received variable args (but really just 1)
                "performs spkmeans (kmeanspp with specific values)" // documentation
        },
        {
                "eigengap_heuristic", // name exposed to Python
                eigengap_heuristic1, // C wrapper function
                      METH_VARARGS, // received variable args (but really just 1)
                "Returns relevant k" // documentation
        },
        {
                "calculateUmatrix", // name exposed to Python
                calculateUmatrix1, // C wrapper function
                      METH_VARARGS, // received variable args (but really just 1)
                "Returns the U matrix" // documentation
        },
        { NULL, NULL, 0, NULL }
};

static struct PyModuleDef mykmeanssp = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",     // name of module exposed to Python
        "A module that performs spectral clustering and returns wam, ddg, gl, jacobi matrices", // module documentation
        -1,
        kmeanssp_methods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject *m;
    m = PyModule_Create(&mykmeanssp);
    if (!m) {
        return NULL;
    }
    return m;
}
