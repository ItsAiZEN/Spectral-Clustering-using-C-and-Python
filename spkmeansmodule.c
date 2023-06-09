#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include "spkmeans.h"

static PyObject *create_pyList_from_C(double **C_matrix, int num_of_rows, int num_of_cols);

double **create_Cmatrix_from_Py(PyObject *pyList, int num_of_rows, int num_of_cols);


static PyObject *create_pyList_from_C(double **C_matrix, int num_of_rows, int num_of_cols) {
    int i, j;
    PyObject *num;
    PyObject *py_matrix;
    PyObject *py_matrix_row;

    py_matrix = PyList_New(num_of_rows);
    if (py_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_rows; i++) {
        py_matrix_row = PyList_New(num_of_cols);
        if (py_matrix_row == NULL) {
            print_error();
        }
        PyList_SetItem(py_matrix, i, py_matrix_row);
        for (j = 0; j < num_of_cols; j++) {
            num = PyFloat_FromDouble(C_matrix[i][j]);
            PyList_SetItem(py_matrix_row, j, num);
        }
    }

    return py_matrix;
}

double **create_Cmatrix_from_Py(PyObject *pyList, int num_of_rows, int num_of_cols) {
    double **C_matrix;
    double num;
    int i, j;
    PyObject *item_row;
    PyObject *item_col;

    C_matrix = (double **) malloc(num_of_rows * sizeof(double *));
    if (C_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_rows; i++) {
        C_matrix[i] = (double *) malloc(num_of_cols * sizeof(double));
        if (C_matrix[i] == NULL) {
            print_error();
        }
    }

    for (i = 0; i < num_of_rows; i++) {
        item_row = PyList_GetItem(pyList, i);
        for (j = 0; j < num_of_cols; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            C_matrix[i][j] = num;
        }
    }
    return C_matrix;
}

static PyObject *wam_module(PyObject *self, PyObject *args) {
    double **Cwam_matrix;
    double **vector_list;
    int num_of_vectors, vector_dimension, i;
    PyObject *pyObj_vectors;
    PyObject *py_wam_matrix;
    PyObject *return_val;

    /* check args are as expected */
    if (!PyArg_ParseTuple(args, "Oii", &pyObj_vectors, &num_of_vectors, &vector_dimension)) {
        return NULL;
    }

    /* transform input vector_list matrix (given as list in Python) to C so that we can use it as an argument in the C ddg function*/
    vector_list = create_Cmatrix_from_Py(pyObj_vectors, num_of_vectors, vector_dimension);

    /* call wam function in C */
    Cwam_matrix = wam(vector_list, num_of_vectors, vector_dimension);

    /* create new python list (to return to python) and insert Cwam_matrix into the list */
    py_wam_matrix = create_pyList_from_C(Cwam_matrix, num_of_vectors, num_of_vectors);

    /* free memory (Cwam_matrix + vector_list in C) */
    for (i = 0; i < num_of_vectors; i++) {
        free(Cwam_matrix[i]);
    }
    free(Cwam_matrix);

    for (i = 0; i < num_of_vectors; i++) {
        free(vector_list[i]);
    }
    free(vector_list);

    /* we are returning the python list that contains the wam matrix */
    return_val = Py_BuildValue("O", py_wam_matrix);
    return return_val;

}

static PyObject *ddg_module(PyObject *self, PyObject *args) {
    double **Cddg_matrix;
    double **C_wam_matrix;
    int num_of_vectors, i;
    PyObject *return_val;
    PyObject *py_wam_matrix;
    PyObject *py_ddg_matrix;

    /* check args are as expected */
    if (!PyArg_ParseTuple(args, "Oi", &py_wam_matrix, &num_of_vectors)) {
        return NULL;
    }

    /* transform input wam matrix (given as list in Python) to C so that we can use it as an argument in the C ddg function*/
    C_wam_matrix = create_Cmatrix_from_Py(py_wam_matrix, num_of_vectors, num_of_vectors);

    /* call ddg function in C */
    Cddg_matrix = ddg(C_wam_matrix, num_of_vectors);

    /* create new python list (to return to python) and insert Cddg_matrix into the list */
    py_ddg_matrix = create_pyList_from_C(Cddg_matrix, num_of_vectors, num_of_vectors);

    /* free memory (Cddg_matrix + C_wam_matrix) */
    for (i = 0; i < num_of_vectors; i++) {
        free(Cddg_matrix[i]);
    }
    free(Cddg_matrix);

    for (i = 0; i < num_of_vectors; i++) {
        free(C_wam_matrix[i]);
    }
    free(C_wam_matrix);

    /* we are returning the python list that contains the ddg matrix */
    return_val = Py_BuildValue("O", py_ddg_matrix);
    return return_val;

}

static PyObject *gl_module(PyObject *self, PyObject *args) {
    double **C_wam_matrix;
    double **C_ddg_matrix;
    double **C_gl_matrix;
    int num_of_vectors, i;
    PyObject *return_val;
    PyObject *py_ddg_matrix;
    PyObject *py_wam_matrix;
    PyObject *py_gl_matrix;

    /* check args are as expected */
    if (!PyArg_ParseTuple(args, "OOi", &py_ddg_matrix, &py_wam_matrix, &num_of_vectors)) {
        return NULL;
    }

    /* transform input wam matrix (given as list in Python) to C so that we can use it as an argument in the C gl function */
    C_wam_matrix = create_Cmatrix_from_Py(py_wam_matrix, num_of_vectors, num_of_vectors);

    /* transform input ddg matrix (given as list in Python) to C so that we can use it as an argument in the C gl function */
    C_ddg_matrix = create_Cmatrix_from_Py(py_ddg_matrix, num_of_vectors, num_of_vectors);

    /* call gl function in C */
    C_gl_matrix = gl(C_wam_matrix, C_ddg_matrix, num_of_vectors);

    /* create new python list (to return to python) and insert C_gl_matrix into the list */
    py_gl_matrix = create_pyList_from_C(C_gl_matrix, num_of_vectors, num_of_vectors);

    /* free memory (C_ddg_matrix + C_wam_matrix) */
    for (i = 0; i < num_of_vectors; i++) {
        free(C_ddg_matrix[i]);
    }
    free(C_ddg_matrix);

    for (i = 0; i < num_of_vectors; i++) {
        free(C_wam_matrix[i]);
    }
    free(C_wam_matrix);

    /* we are returning the python list that contains the gl matrix */
    return_val = Py_BuildValue("O", py_gl_matrix);
    return return_val;
}

static PyObject *jacobi_module(PyObject *self, PyObject *args) {
    double **C_matrix;
    double **C_jacobi_matrix;
    int num_of_vectors, i;
    PyObject *return_val;
    PyObject *py_matrix;
    PyObject *py_jacobi_matrix;

    /* check args are as expected */
    if (!PyArg_ParseTuple(args, "Oi", &py_matrix, &num_of_vectors)) {
        return NULL;
    }

    /* transform input matrix (given as list in Python) to C so that we can use it as an argument in the C jacobi function */
    C_matrix = create_Cmatrix_from_Py(py_matrix, num_of_vectors, num_of_vectors);

    /* call jacobi function in C */
    C_jacobi_matrix = jacobi(C_matrix, num_of_vectors);

    /* create new python list (to return to python) and insert C_jacobi_matrix into the list */
    py_jacobi_matrix = create_pyList_from_C(C_jacobi_matrix, num_of_vectors + 1, num_of_vectors);

    /* free memory (C_matrix + C_jacobi_matrix) */
    for (i = 0; i < num_of_vectors; i++) {
        free(C_matrix[i]);
    }
    free(C_matrix);

    for (i = 0; i < (num_of_vectors + 1); i++) {
        free(C_jacobi_matrix[i]);
    }
    free(C_jacobi_matrix);

    /* we are returning the python list that contains the jacobi matrix */
    return_val = Py_BuildValue("O", py_jacobi_matrix);
    return return_val;
}

static PyObject *spk_module(PyObject *self, PyObject *args) {
    double **vector_list;
    double **init_centroids;
    int n, k, i;
    PyObject *pyObjVector_list;
    PyObject *pyObjInit_centroids;


    /* check args are as expected */
    if (!PyArg_ParseTuple(args, "iiOO", &k, &n, &pyObjVector_list, &pyObjInit_centroids)) {
        return NULL;
    }

    /* check args are valid */
    if (n < 0 || k < 0) {
        return NULL;
    }

    /* transform input vector list (given in Python) to C so that we can use it as an argument in the C kmeanspp function*/
    vector_list = create_Cmatrix_from_Py(pyObjVector_list, n, k);

    /* transform input init_centroids list (given in Python) to C so that we can use it as an argument in the C kmeanspp function*/
    init_centroids = create_Cmatrix_from_Py(pyObjInit_centroids, k, k);

    /* call kmeanspp function in C */
    kmeanspp(k, 300, k, n, vector_list, 0, init_centroids);

    /* free memory (vector_list + init_centroids) */
    for (i = 0; i < n; i++) {
        free(vector_list[i]);
    }
    free(vector_list);

    for (i = 0; i < k; i++) {
        free(init_centroids[i]);
    }
    free(init_centroids);

    /* kmeanspp prints the result and is of void type so we do not need to return a value to Python */
    Py_RETURN_NONE;
}

static PyObject *eigengap_heuristic_module(PyObject *self, PyObject *args) {
    double **C_jacobi_matrix;
    int num_of_vectors, i, k;
    PyObject *return_val;
    PyObject *py_jacobi_matrix;

    /* check args are as expected */
    if (!PyArg_ParseTuple(args, "Oi", &py_jacobi_matrix, &num_of_vectors)) {
        return NULL;
    }

    /* transform input jacobi matrix (given as list in Python) to C so that we can use it as an argument in the C eigengap_heuristic function */
    C_jacobi_matrix = create_Cmatrix_from_Py(py_jacobi_matrix, num_of_vectors + 1, num_of_vectors);

    /* call eigengap heuristic function in C */
    k = eigengap_heuristic(C_jacobi_matrix, num_of_vectors);

    /* free memory (C_jacobi_matrix) */
    for (i = 0; i < (num_of_vectors + 1); i++) {
        free(C_jacobi_matrix[i]);
    }
    free(C_jacobi_matrix);

    return_val = Py_BuildValue("i", k);
    return return_val;
}

static PyObject *calculateUmatrix_module(PyObject *self, PyObject *args) {
    double **C_jacobi_matrix;
    double **C_U_matrix;
    int num_of_vectors, i, k;
    PyObject *py_jacobi_matrix;
    PyObject *py_u_matrix;
    PyObject *return_val;

    /* check args are as expected */
    if (!PyArg_ParseTuple(args, "Oii", &py_jacobi_matrix, &num_of_vectors, &k)) {
        return NULL;
    }

    /* transform input jacobi matrix (given as list in Python) to C so that we can use it as an argument in the C calculateUmatrix function */
    C_jacobi_matrix = create_Cmatrix_from_Py(py_jacobi_matrix, num_of_vectors + 1, num_of_vectors);

    /* call calculateUmatrix function in C */
    C_U_matrix = calculateUmatrix(C_jacobi_matrix, num_of_vectors, k);

    /* create new python list (to return to python) and insert C_U_matrix into the list */
    py_u_matrix = create_pyList_from_C(C_U_matrix, num_of_vectors, k);

    /* free memory (C_jacobi_matrix + C_U_matrix) */
    for (i = 0; i < (num_of_vectors + 1); i++) {
        free(C_jacobi_matrix[i]);
    }
    free(C_jacobi_matrix);

    for (i = 0; i < num_of_vectors; i++) {
        free(C_U_matrix[i]);
    }
    free(C_U_matrix);

    /* we are returning the python list that contains the U matrix */
    return_val = Py_BuildValue("O", py_u_matrix);
    return return_val;
}


static PyMethodDef kmeanssp_methods[] = {
        {
                "wam", /* name exposed to Python */
                      wam_module, /* C wrapper function */
                            METH_VARARGS,
                               "Returns the wam matrix"
        },
        {
                "ddg", /* name exposed to Python */
                      ddg_module, /* C wrapper function */
                            METH_VARARGS,
                               "Returns the ddg matrix"
        },
        {
                "gl", /* name exposed to Python */
                      gl_module, /* C wrapper function */
                            METH_VARARGS,
                               "Returns the gl matrix"
        },
        {
                "jacobi", /* name exposed to Python */
                      jacobi_module, /* C wrapper function */
                            METH_VARARGS,
                               "Returns the jacobi matrix"
        },
        {
                "spk", /* name exposed to Python */
                      spk_module, /* C wrapper function */
                            METH_VARARGS,
                               "performs spkmeans (kmeanspp with specific values)"
        },
        {
                "eigengap_heuristic", /* name exposed to Python */
                      eigengap_heuristic_module, /* C wrapper function */
                            METH_VARARGS,
                               "Returns relevant k"
        },
        {
                "calculateUmatrix", /* name exposed to Python */
                      calculateUmatrix_module, /* C wrapper function */
                            METH_VARARGS,
                               "Returns the U matrix"
        },
        {       NULL, NULL, 0, NULL}
};

static struct PyModuleDef mykmeanssp = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",     /* name of module exposed to Python */
        "A module that performs spectral clustering and can return wam, ddg, gl, jacobi matrices",
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
